#include "cartesian3d_assembly.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace cartesian3d {
namespace {

std::vector<std::int64_t> resolve_block_ids(const SpatialDefinition& definition, const UnstructuredHex8Mesh& source_mesh) {
    if (definition.regions.empty()) throw std::invalid_argument("Cartesian three-dimensional assembly requires at least one region");
    if (!definition.contacts.empty()) throw std::invalid_argument("Cartesian three-dimensional stage B does not support contact");
    std::vector<std::int64_t> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        if (region.name.empty() || region.block.empty()) throw std::invalid_argument("Cartesian three-dimensional region names and blocks must be nonempty");
        if (region.strain_formulation != StrainFormulation::small) throw std::invalid_argument("Cartesian three-dimensional stage B supports only small strain");
        const std::int64_t block_id = source_mesh.element_block(region.block).id;
        if (std::find(result.begin(), result.end(), block_id) != result.end()) throw std::invalid_argument("Cartesian three-dimensional regions must use distinct element blocks");
        result.push_back(block_id);
    }
    return result;
}

std::vector<Hex8RegionMesh> build_meshes(const SpatialDefinition& definition, const UnstructuredHex8Mesh& source_mesh) {
    std::vector<Hex8RegionMesh> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) result.push_back(Hex8RegionMesh::from_unstructured_block(source_mesh, region.block));
    return result;
}

std::size_t total_nodes(const std::vector<Hex8RegionMesh>& meshes) {
    std::size_t result = 0;
    for (const Hex8RegionMesh& mesh : meshes) {
        if (mesh.nodes().size() > std::numeric_limits<std::size_t>::max() - result) throw std::length_error("Cartesian three-dimensional node count overflows");
        result += mesh.nodes().size();
    }
    return result;
}

Hex8Coordinates element_coordinates(const Hex8RegionMesh& mesh, const Hex8Element& element) {
    Hex8Coordinates result{};
    for (std::size_t node = 0; node < 8; ++node) result[node] = mesh.nodes().at(element.nodes[node]);
    return result;
}

Quad4FaceCoordinates face_coordinates(const Hex8RegionMesh& mesh, const Quad4FaceElement& face) {
    Quad4FaceCoordinates result{};
    for (std::size_t node = 0; node < 4; ++node) result[node] = mesh.nodes().at(face.nodes[node]);
    return result;
}

void validate_dirichlet(std::vector<DirichletCondition>& conditions) {
    std::sort(conditions.begin(), conditions.end(), [](const DirichletCondition& lhs, const DirichletCondition& rhs) { return lhs.dof < rhs.dof; });
    for (std::size_t index = 1; index < conditions.size(); ++index) {
        if (conditions[index - 1].dof != conditions[index].dof) continue;
        if (conditions[index - 1].value != conditions[index].value) throw std::invalid_argument("Cartesian three-dimensional boundary has conflicting Dirichlet values");
        throw std::invalid_argument("Cartesian three-dimensional boundary has duplicate Dirichlet values");
    }
}

} // namespace

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _definition(std::move(definition)), _block_ids(resolve_block_ids(_definition, source_mesh)), _meshes(build_meshes(_definition, source_mesh)), _dof_map(total_nodes(_meshes)), _load_factor(0.0), _time(0.0) {
    _node_offsets.push_back(0);
    _element_offsets.push_back(0);
    for (const Hex8RegionMesh& mesh : _meshes) {
        _node_offsets.push_back(_node_offsets.back() + mesh.nodes().size());
        _element_offsets.push_back(_element_offsets.back() + mesh.elements().size());
    }
    _geometries.resize(_meshes.size());
    for (std::size_t region_value = 0; region_value < _meshes.size(); ++region_value) {
        for (const Hex8Element& element : _meshes[region_value].elements()) _geometries[region_value].push_back(make_hex8_geometry(element_coordinates(_meshes[region_value], element)));
    }

    for (const BoundaryConditionDefinition& boundary : _definition.boundary_conditions) {
        if (boundary.name.empty() || boundary.boundary.empty()) throw std::invalid_argument("Cartesian three-dimensional boundary names must be nonempty");
        if (boundary.scale_with_load && !boundary.function.empty()) throw std::invalid_argument("Boundary condition cannot combine scale_with_load and a time function: " + boundary.name);
        const std::int64_t block_id = source_mesh.side_set_block_id(boundary.boundary);
        const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
        if (found == _block_ids.end()) throw std::invalid_argument("Boundary belongs to an undeclared three-dimensional block: " + boundary.boundary);
        const std::size_t region_value = static_cast<std::size_t>(found - _block_ids.begin());
        const Hex8RegionBoundary mapped = _meshes[region_value].map_side_set(source_mesh, boundary.boundary);
        if (boundary.type == BoundaryConditionType::dirichlet) {
            for (std::size_t local_node : mapped.nodes) {
                const std::size_t dof = _dof_map.dof(boundary.field, global_node(region_value, local_node));
                _dirichlet_conditions.push_back({dof, boundary.value});
                if (boundary.scale_with_load || !boundary.function.empty()) _controlled_dirichlet.push_back({dof, {boundary.value, boundary.scale_with_load, boundary.function}});
            }
            continue;
        }
        if (boundary.use_displaced_geometry) throw std::invalid_argument("Cartesian three-dimensional stage B loads use the reference configuration");
        std::size_t kernel = 0;
        if (boundary.type == BoundaryConditionType::pressure) {
            kernel = _pressure_kernels.size();
            _pressure_controls.push_back({boundary.value, boundary.scale_with_load, boundary.function});
            _pressure_kernels.emplace_back(boundary.value);
        } else if (boundary.type == BoundaryConditionType::traction) {
            CartesianTractionComponent component = CartesianTractionComponent::x;
            if (boundary.field == Field::displacement_y)
                component = CartesianTractionComponent::y;
            else if (boundary.field == Field::displacement_z)
                component = CartesianTractionComponent::z;
            else if (boundary.field != Field::displacement_x)
                throw std::invalid_argument("Three-dimensional traction requires displacement_x, displacement_y, or "
                                            "displacement_z");
            kernel = _traction_kernels.size();
            _traction_controls.push_back({boundary.value, boundary.scale_with_load, boundary.function});
            _traction_kernels.emplace_back(component, boundary.value);
        } else {
            kernel = _convection_kernels.size();
            _convection_controls.push_back({boundary.heat_transfer_coefficient, boundary.ambient_temperature, boundary.coefficient_function, boundary.ambient_temperature_function});
            _convection_kernels.emplace_back(boundary.heat_transfer_coefficient, boundary.ambient_temperature);
        }
        for (const Quad4FaceElement& face : mapped.faces) {
            std::array<std::size_t, 4> nodes{};
            for (std::size_t node = 0; node < 4; ++node) nodes[node] = global_node(region_value, face.nodes[node]);
            _boundary_contributions.push_back({boundary.type == BoundaryConditionType::pressure ? SpatialContributionType::pressure : (boundary.type == BoundaryConditionType::traction ? SpatialContributionType::traction : SpatialContributionType::convection),
                                               kernel, nodes, make_quad4_face_geometry(face_coordinates(_meshes[region_value], face))});
        }
    }
    validate_dirichlet(_dirichlet_conditions);
    refresh_controls();
}

const SpatialDefinition& SpatialAssembly::definition() const noexcept { return _definition; }
const Hex8DofMap& SpatialAssembly::dof_map() const noexcept { return _dof_map; }
std::size_t SpatialAssembly::region_count() const noexcept { return _meshes.size(); }
std::size_t SpatialAssembly::region_index(const std::string& name) const {
    for (std::size_t region_value = 0; region_value < _definition.regions.size(); ++region_value) {
        if (_definition.regions[region_value].name == name) return region_value;
    }
    throw std::invalid_argument("Unknown three-dimensional region: " + name);
}
const RegionDefinition& SpatialAssembly::region(std::size_t index) const { return _definition.regions.at(index); }
const Hex8RegionMesh& SpatialAssembly::region_mesh(std::size_t index) const { return _meshes.at(index); }
std::size_t SpatialAssembly::region_node_offset(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Three-dimensional region index is out of range");
    return _node_offsets[index];
}
std::size_t SpatialAssembly::region_element_count(std::size_t index) const { return _meshes.at(index).elements().size(); }
std::size_t SpatialAssembly::region_element_offset(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Three-dimensional region index is out of range");
    return _element_offsets[index];
}
std::size_t SpatialAssembly::volume_contribution_count() const noexcept { return _element_offsets.back(); }
SpatialContributionType SpatialAssembly::contribution_type(std::size_t contribution_index) const {
    if (contribution_index < volume_contribution_count()) return SpatialContributionType::volume;
    return _boundary_contributions.at(contribution_index - volume_contribution_count()).type;
}
std::pair<std::size_t, std::size_t> SpatialAssembly::element_location(std::size_t contribution_index) const {
    if (contribution_index >= volume_contribution_count()) throw std::out_of_range("Three-dimensional contribution is not a volume element");
    const auto upper = std::upper_bound(_element_offsets.begin(), _element_offsets.end(), contribution_index);
    const std::size_t region_value = static_cast<std::size_t>(upper - _element_offsets.begin() - 1);
    return {region_value, contribution_index - _element_offsets[region_value]};
}
const Hex8Geometry& SpatialAssembly::region_element_geometry(std::size_t region_value, std::size_t element_index) const { return _geometries.at(region_value).at(element_index); }
double SpatialAssembly::region_heat_source(std::size_t region_value) const {
    const RegionDefinition& region_definition = region(region_value);
    const double multiplier = region_definition.heat_source_function.empty() ? _load_factor : function_value(region_definition.heat_source_function);
    return multiplier * region_definition.volumetric_heat_source;
}

void SpatialAssembly::set_load_factor(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Three-dimensional load factor must be finite");
    _load_factor = value;
    refresh_controls();
}
double SpatialAssembly::load_factor() const noexcept { return _load_factor; }
void SpatialAssembly::set_time(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Three-dimensional time must be finite");
    _time = value;
    refresh_controls();
}

std::vector<double> SpatialAssembly::initial_state() const {
    std::vector<double> result(dof_count(), 0.0);
    for (std::size_t region_value = 0; region_value < region_count(); ++region_value) {
        for (std::size_t local_node = 0; local_node < _meshes[region_value].nodes().size(); ++local_node) result[_dof_map.temperature(global_node(region_value, local_node))] = region(region_value).initial_temperature;
    }
    return result;
}

std::size_t SpatialAssembly::dof_count() const noexcept { return _dof_map.dof_count(); }
std::size_t SpatialAssembly::contribution_count() const noexcept { return volume_contribution_count() + _boundary_contributions.size(); }
const std::vector<DirichletCondition>& SpatialAssembly::dirichlet_conditions() const noexcept { return _dirichlet_conditions; }
void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count()) throw std::invalid_argument("Three-dimensional state size does not match the problem");
    if (!std::all_of(state.begin(), state.end(), [](double value) { return std::isfinite(value); })) throw std::domain_error("Three-dimensional state must contain only finite values");
}
std::size_t SpatialAssembly::contribution_dof_count(std::size_t contribution_index) const {
    if (contribution_index >= contribution_count()) throw std::out_of_range("Three-dimensional contribution index is out of range");
    return contribution_index < volume_contribution_count() ? hex8_local_dof_count : quad4_face_local_dof_count;
}
void SpatialAssembly::fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const {
    if (contribution_index < volume_contribution_count()) {
        const auto location = element_location(contribution_index);
        std::array<std::size_t, 8> nodes{};
        const Hex8Element& element = _meshes[location.first].elements()[location.second];
        for (std::size_t node = 0; node < 8; ++node) nodes[node] = global_node(location.first, element.nodes[node]);
        const Hex8LocalDofs fixed = _dof_map.local_dofs(nodes);
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const BoundaryContribution& contribution = _boundary_contributions.at(contribution_index - volume_contribution_count());
    const Quad4FaceLocalDofs fixed = _dof_map.face_local_dofs(contribution.nodes);
    dofs.assign(fixed.begin(), fixed.end());
}

Quad4FaceLocalResidual SpatialAssembly::boundary_residual(std::size_t contribution_index, const Quad4FaceLocalValues& state) const {
    const BoundaryContribution& contribution = _boundary_contributions.at(contribution_index - volume_contribution_count());
    if (contribution.type == SpatialContributionType::pressure) return _pressure_kernels[contribution.kernel].residual(contribution.geometry, state);
    if (contribution.type == SpatialContributionType::traction) return _traction_kernels[contribution.kernel].residual(contribution.geometry, state);
    return _convection_kernels[contribution.kernel].residual(contribution.geometry, state);
}

Quad4FaceLocalSystem SpatialAssembly::boundary_system(std::size_t contribution_index, const Quad4FaceLocalValues& state) const {
    const BoundaryContribution& contribution = _boundary_contributions.at(contribution_index - volume_contribution_count());
    if (contribution.type == SpatialContributionType::pressure) return _pressure_kernels[contribution.kernel].linearize(contribution.geometry, state);
    if (contribution.type == SpatialContributionType::traction) return _traction_kernels[contribution.kernel].linearize(contribution.geometry, state);
    return _convection_kernels[contribution.kernel].linearize(contribution.geometry, state);
}

double SpatialAssembly::function_value(const std::string& name) const {
    const auto found = std::find_if(_definition.time_tables.begin(), _definition.time_tables.end(), [&name](const PiecewiseLinearTimeTable& table) { return table.name() == name; });
    if (found == _definition.time_tables.end()) throw std::invalid_argument("Unknown time-table function: " + name);
    return found->value(_time);
}
double SpatialAssembly::controlled_value(const ControlledScalar& control) const {
    if (!control.function.empty()) return control.base_value * function_value(control.function);
    return control.base_value * (control.scale_with_load ? _load_factor : 1.0);
}
void SpatialAssembly::refresh_controls() {
    for (const ControlledDirichlet& controlled : _controlled_dirichlet) {
        const auto condition = std::lower_bound(_dirichlet_conditions.begin(), _dirichlet_conditions.end(), controlled.dof, [](const DirichletCondition& candidate, std::size_t dof) { return candidate.dof < dof; });
        if (condition == _dirichlet_conditions.end() || condition->dof != controlled.dof) throw std::logic_error("Three-dimensional controlled Dirichlet mapping is invalid");
        condition->value = controlled_value(controlled.control);
    }
    for (std::size_t kernel = 0; kernel < _pressure_kernels.size(); ++kernel) _pressure_kernels[kernel].set_pressure(controlled_value(_pressure_controls[kernel]));
    for (std::size_t kernel = 0; kernel < _traction_kernels.size(); ++kernel) _traction_kernels[kernel].set_traction(controlled_value(_traction_controls[kernel]));
    for (std::size_t kernel = 0; kernel < _convection_kernels.size(); ++kernel) {
        const ConvectionControl& control = _convection_controls[kernel];
        const double coefficient = control.coefficient_function.empty() ? control.coefficient : control.coefficient * function_value(control.coefficient_function);
        const double ambient = control.ambient_function.empty() ? control.ambient : control.ambient * function_value(control.ambient_function);
        _convection_kernels[kernel].set_properties(coefficient, ambient);
    }
}
std::size_t SpatialAssembly::global_node(std::size_t region_value, std::size_t local_node) const {
    if (local_node >= _meshes.at(region_value).nodes().size()) throw std::out_of_range("Three-dimensional local node is out of range");
    return _node_offsets.at(region_value) + local_node;
}

} // namespace cartesian3d
} // namespace fuelsim
