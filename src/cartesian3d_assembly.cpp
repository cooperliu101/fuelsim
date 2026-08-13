#include "cartesian3d_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
namespace fuelsim::cartesian {
namespace {
std::vector<std::int64_t> resolve_block_ids(
    const SpatialDefinition& definition, const UnstructuredHex8Mesh& source_mesh) {
    if (definition.regions.empty())
        throw std::invalid_argument("Cartesian three-dimensional assembly requires at least one region");
    if (!definition.contacts.empty())
        throw std::invalid_argument("Cartesian three-dimensional stage B does not support contact");
    std::vector<std::int64_t> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        if (region.name.empty() || region.block.empty())
            throw std::invalid_argument("Cartesian three-dimensional region names and blocks must be nonempty");
        if (region.strain_formulation != StrainFormulation::small)
            throw std::invalid_argument("Cartesian three-dimensional stage B supports only small strain");
        const std::int64_t block_id = source_mesh.element_block(region.block).id;
        if (std::find(result.begin(), result.end(), block_id) != result.end())
            throw std::invalid_argument("Cartesian three-dimensional regions must use distinct element blocks");
        result.push_back(block_id);
    }
    return result;
}
std::vector<Hex8RegionMesh> build_meshes(const SpatialDefinition& definition, const UnstructuredHex8Mesh& source_mesh) {
    std::vector<Hex8RegionMesh> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions)
        result.push_back(Hex8RegionMesh::from_unstructured_block(source_mesh, region.block));
    return result;
}
std::size_t total_nodes(const std::vector<Hex8RegionMesh>& meshes) {
    std::size_t result = 0;
    for (const Hex8RegionMesh& mesh : meshes) {
        if (mesh.nodes().size() > std::numeric_limits<std::size_t>::max() - result)
            throw std::length_error("Cartesian three-dimensional node count overflows");
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
    std::sort(conditions.begin(), conditions.end(),
        [](const DirichletCondition& lhs, const DirichletCondition& rhs) { return lhs.dof < rhs.dof; });
    for (std::size_t index = 1; index < conditions.size(); ++index) {
        if (conditions[index - 1].dof != conditions[index].dof) continue;
        if (conditions[index - 1].value != conditions[index].value)
            throw std::invalid_argument("Cartesian three-dimensional boundary has conflicting Dirichlet values");
        throw std::invalid_argument("Cartesian three-dimensional boundary has duplicate Dirichlet values");
    }
}
Hex8LocalValues hex8_values(const std::vector<double>& values) {
    if (values.size() != hex8_local_dof_count)
        throw std::invalid_argument("HEX8 contribution state must contain 32 DOFs");
    Hex8LocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}
Quad4FaceLocalValues face_values(const std::vector<double>& values) {
    if (values.size() != quad4_face_local_dof_count)
        throw std::invalid_argument("Three-dimensional face state must contain 16 DOFs");
    Quad4FaceLocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}
} // namespace
SpatialAssembly::SpatialAssembly(
    SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh, std::vector<double> heat_capacities)
    : _definition(std::move(definition)), _block_ids(resolve_block_ids(_definition, source_mesh)),
      _meshes(build_meshes(_definition, source_mesh)), _dof_map(total_nodes(_meshes), DofLayout::cartesian_3d),
      _load_factor(0.0), _time(0.0) {
    _node_offsets.push_back(0);
    _element_offsets.push_back(0);
    for (const Hex8RegionMesh& mesh : _meshes) {
        _node_offsets.push_back(_node_offsets.back() + mesh.nodes().size());
        _element_offsets.push_back(_element_offsets.back() + mesh.elements().size());
    }
    _geometries.resize(_meshes.size());
    for (std::size_t region = 0; region < _meshes.size(); ++region)
        for (const Hex8Element& element : _meshes[region].elements())
            _geometries[region].push_back(make_hex8_geometry(element_coordinates(_meshes[region], element)));
    for (const BoundaryConditionDefinition& boundary : _definition.boundary_conditions) {
        if (boundary.name.empty() || boundary.boundary.empty())
            throw std::invalid_argument("Cartesian three-dimensional boundary names must be nonempty");
        if (boundary.scale_with_load && !boundary.function.empty())
            throw std::invalid_argument(
                "Boundary condition cannot combine scale_with_load and a time function: " + boundary.name);
        const std::int64_t block_id = source_mesh.side_set_block_id(boundary.boundary);
        const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
        if (found == _block_ids.end())
            throw std::invalid_argument(
                "Boundary belongs to an undeclared three-dimensional block: " + boundary.boundary);
        const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
        const Hex8RegionBoundary mapped = _meshes[region].map_side_set(source_mesh, boundary.boundary);
        if (boundary.type == BoundaryConditionType::dirichlet) {
            for (std::size_t local_node : mapped.nodes) {
                const std::size_t dof = _dof_map.dof(boundary.field, global_node(region, local_node));
                _dirichlet_conditions.push_back({dof, boundary.value});
                if (boundary.scale_with_load || !boundary.function.empty())
                    _controlled_dirichlet.push_back(
                        {dof, {boundary.value, boundary.scale_with_load, boundary.function}});
            }
            continue;
        }
        if (boundary.use_displaced_geometry)
            throw std::invalid_argument("Cartesian three-dimensional stage B loads use the reference configuration");
        const std::size_t kernel = _boundary_kernels.size();
        if (boundary.type == BoundaryConditionType::pressure) {
            _boundary_controls.emplace_back(
                ControlledScalar{boundary.value, boundary.scale_with_load, boundary.function});
            _boundary_kernels.emplace_back(boundary.value);
        } else if (boundary.type == BoundaryConditionType::traction) {
            CartesianTractionComponent component = CartesianTractionComponent::x;
            if (boundary.field == Field::displacement_y)
                component = CartesianTractionComponent::y;
            else if (boundary.field == Field::displacement_z)
                component = CartesianTractionComponent::z;
            else if (boundary.field != Field::displacement_x)
                throw std::invalid_argument("Three-dimensional traction requires a displacement field");
            _boundary_controls.emplace_back(
                ControlledScalar{boundary.value, boundary.scale_with_load, boundary.function});
            _boundary_kernels.emplace_back(component, boundary.value);
        } else {
            _boundary_controls.emplace_back(ConvectionControl{boundary.heat_transfer_coefficient,
                boundary.ambient_temperature, boundary.coefficient_function, boundary.ambient_temperature_function});
            _boundary_kernels.emplace_back(boundary.heat_transfer_coefficient, boundary.ambient_temperature);
        }
        for (const Quad4FaceElement& face : mapped.faces) {
            std::array<std::size_t, 4> nodes{};
            for (std::size_t node = 0; node < 4; ++node) nodes[node] = global_node(region, face.nodes[node]);
            _boundary_contributions.push_back(
                {boundary.type == BoundaryConditionType::pressure
                        ? SpatialContributionType::pressure
                        : (boundary.type == BoundaryConditionType::traction ? SpatialContributionType::traction
                                                                            : SpatialContributionType::convection),
                    kernel, nodes, make_quad4_face_geometry(face_coordinates(_meshes[region], face))});
        }
    }
    validate_dirichlet(_dirichlet_conditions);
    if (!heat_capacities.empty() && heat_capacities.size() != region_count())
        throw std::invalid_argument("Three-dimensional heat-capacity count must match the regions");
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index) {
        IsotropicThermoelasticMaterial material(region(region_index).material);
        if (heat_capacities.empty())
            _kernels.emplace_back(std::move(material), region_heat_source(region_index));
        else
            _kernels.emplace_back(std::move(material), heat_capacities[region_index], region_heat_source(region_index));
    }
    refresh_controls();
}
std::size_t SpatialAssembly::region_index(const std::string& name) const {
    for (std::size_t region = 0; region < _definition.regions.size(); ++region)
        if (_definition.regions[region].name == name) return region;
    throw std::invalid_argument("Unknown three-dimensional region: " + name);
}
std::size_t SpatialAssembly::region_node_offset(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Three-dimensional region index is out of range");
    return _node_offsets[index];
}
std::size_t SpatialAssembly::region_element_count(std::size_t index) const {
    return _meshes.at(index).elements().size();
}
std::size_t SpatialAssembly::region_element_offset(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Three-dimensional region index is out of range");
    return _element_offsets[index];
}
SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    if (index < volume_contribution_count()) return SpatialContributionType::volume;
    return _boundary_contributions.at(index - volume_contribution_count()).type;
}
std::pair<std::size_t, std::size_t> SpatialAssembly::element_location(std::size_t index) const {
    if (index >= volume_contribution_count())
        throw std::out_of_range("Three-dimensional contribution is not a volume element");
    const auto upper = std::upper_bound(_element_offsets.begin(), _element_offsets.end(), index);
    const std::size_t region = static_cast<std::size_t>(upper - _element_offsets.begin() - 1);
    return {region, index - _element_offsets[region]};
}
const Hex8Geometry& SpatialAssembly::region_element_geometry(std::size_t region, std::size_t element_index) const {
    return _geometries.at(region).at(element_index);
}
double SpatialAssembly::region_heat_source(std::size_t region_index) const {
    const RegionDefinition& region_definition = region(region_index);
    const double multiplier = region_definition.heat_source_function.empty()
                                  ? _load_factor
                                  : function_value(region_definition.heat_source_function);
    return multiplier * region_definition.volumetric_heat_source;
}
void SpatialAssembly::set_load_factor(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Three-dimensional load factor must be finite");
    _load_factor = value;
    refresh_controls();
}
void SpatialAssembly::set_time(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Three-dimensional time must be finite");
    _time = value;
    for (Hex8ThermoelasticKernel& kernel : _kernels) kernel.set_time(value);
    refresh_controls();
}
std::vector<double> SpatialAssembly::initial_state() const {
    std::vector<double> result(dof_count(), 0.0);
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index)
        for (std::size_t local_node = 0; local_node < _meshes[region_index].nodes().size(); ++local_node)
            result[_dof_map.temperature(global_node(region_index, local_node))] =
                region(region_index).initial_temperature;
    return result;
}
void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional state size does not match the problem");
    if (!std::all_of(state.begin(), state.end(), [](double value) { return std::isfinite(value); }))
        throw std::domain_error("Three-dimensional state must contain only finite values");
}
void SpatialAssembly::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    if (index < volume_contribution_count()) {
        const auto location = element_location(index);
        std::array<std::size_t, 8> nodes{};
        const Hex8Element& element = _meshes[location.first].elements()[location.second];
        for (std::size_t node = 0; node < 8; ++node) nodes[node] = global_node(location.first, element.nodes[node]);
        Hex8LocalDofs fixed{};
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            fixed[node] = _dof_map.temperature(nodes[node]);
            fixed[8 + node] = _dof_map.displacement_x(nodes[node]);
            fixed[16 + node] = _dof_map.displacement_y(nodes[node]);
            fixed[24 + node] = _dof_map.displacement_z(nodes[node]);
        }
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - volume_contribution_count());
    Quad4FaceLocalDofs fixed{};
    for (std::size_t node = 0; node < entry.nodes.size(); ++node) {
        fixed[node] = _dof_map.temperature(entry.nodes[node]);
        fixed[4 + node] = _dof_map.displacement_x(entry.nodes[node]);
        fixed[8 + node] = _dof_map.displacement_y(entry.nodes[node]);
        fixed[12 + node] = _dof_map.displacement_z(entry.nodes[node]);
    }
    dofs.assign(fixed.begin(), fixed.end());
}
Hex8LocalValues SpatialAssembly::volume_state(std::size_t index, const std::vector<double>& global_state) const {
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    if (dofs.size() != hex8_local_dof_count)
        throw std::logic_error("HEX8 volume contribution has an invalid DOF layout");
    Hex8LocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}
void SpatialAssembly::contribution_residual(std::size_t index, const std::vector<double>& state,
    const std::vector<double>* committed_solution, double time_step, std::vector<double>& residual) const {
    if (index < volume_contribution_count()) {
        const auto location = element_location(index);
        const Hex8LocalValues current = hex8_values(state);
        const Hex8LocalResidual result =
            committed_solution == nullptr
                ? _kernels[location.first].residual(region_element_geometry(location.first, location.second), current)
                : _kernels[location.first].residual(region_element_geometry(location.first, location.second), current,
                      volume_state(index, *committed_solution), time_step);
        residual.assign(result.begin(), result.end());
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - volume_contribution_count());
    const Quad4FaceLocalValues current = face_values(state);
    const Quad4FaceLocalResidual result = _boundary_kernels[entry.kernel].residual(entry.geometry, current);
    residual.assign(result.begin(), result.end());
}
void SpatialAssembly::contribution_system(std::size_t index, const std::vector<double>& state,
    const std::vector<double>* committed_solution, double time_step, std::vector<double>& residual,
    std::vector<double>& jacobian) const {
    if (index < volume_contribution_count()) {
        const auto location = element_location(index);
        const Hex8LocalValues current = hex8_values(state);
        const Hex8LocalSystem result =
            committed_solution == nullptr
                ? _kernels[location.first].linearize(region_element_geometry(location.first, location.second), current)
                : _kernels[location.first].linearize(region_element_geometry(location.first, location.second), current,
                      volume_state(index, *committed_solution), time_step);
        residual.assign(result.residual.begin(), result.residual.end());
        jacobian.assign(result.jacobian.begin(), result.jacobian.end());
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - volume_contribution_count());
    const Quad4FaceLocalValues current = face_values(state);
    const Quad4FaceLocalSystem result = _boundary_kernels[entry.kernel].linearize(entry.geometry, current);
    residual.assign(result.residual.begin(), result.residual.end());
    jacobian.assign(result.jacobian.begin(), result.jacobian.end());
}
std::array<SymmetricTensor3Values, 8> SpatialAssembly::stress(
    std::size_t region, std::size_t element, const std::vector<double>& state) const {
    return _kernels.at(region).stress_values(
        region_element_geometry(region, element), volume_state(region_element_offset(region) + element, state));
}
double SpatialAssembly::heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const {
    return _kernels.at(region).heat_capacity(temperature, position.x, position.y, position.z);
}
double SpatialAssembly::function_value(const std::string& name) const {
    const auto found = std::find_if(_definition.time_tables.begin(), _definition.time_tables.end(),
        [&name](const PiecewiseLinearTimeTable& table) { return table.name() == name; });
    if (found == _definition.time_tables.end()) throw std::invalid_argument("Unknown time-table function: " + name);
    return found->value(_time);
}
double SpatialAssembly::controlled_value(const ControlledScalar& control) const {
    if (!control.function.empty()) return control.base_value * function_value(control.function);
    return control.base_value * (control.scale_with_load ? _load_factor : 1.0);
}
void SpatialAssembly::refresh_controls() {
    for (std::size_t region = 0; region < region_count(); ++region)
        _kernels[region].set_volumetric_heat_source(region_heat_source(region));
    for (const ControlledDirichlet& controlled : _controlled_dirichlet) {
        const auto condition = std::lower_bound(_dirichlet_conditions.begin(), _dirichlet_conditions.end(),
            controlled.dof, [](const DirichletCondition& candidate, std::size_t dof) { return candidate.dof < dof; });
        if (condition == _dirichlet_conditions.end() || condition->dof != controlled.dof)
            throw std::logic_error("Three-dimensional controlled Dirichlet mapping is invalid");
        condition->value = controlled_value(controlled.control);
    }
    for (std::size_t kernel = 0; kernel < _boundary_kernels.size(); ++kernel) {
        if (std::holds_alternative<ControlledScalar>(_boundary_controls[kernel])) {
            _boundary_kernels[kernel].set_load(
                controlled_value(std::get<ControlledScalar>(_boundary_controls[kernel])));
            continue;
        }
        const ConvectionControl& control = std::get<ConvectionControl>(_boundary_controls[kernel]);
        const double coefficient =
            control.coefficient *
            (control.coefficient_function.empty() ? 1.0 : function_value(control.coefficient_function));
        const double ambient =
            control.ambient * (control.ambient_function.empty() ? 1.0 : function_value(control.ambient_function));
        _boundary_kernels[kernel].set_convection(coefficient, ambient);
    }
}
std::size_t SpatialAssembly::global_node(std::size_t region, std::size_t local_node) const {
    if (local_node >= _meshes.at(region).nodes().size())
        throw std::out_of_range("Three-dimensional local node is out of range");
    return _node_offsets.at(region) + local_node;
}
} // namespace fuelsim::cartesian
