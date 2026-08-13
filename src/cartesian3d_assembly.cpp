#include "cartesian3d_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace fuelsim::cartesian {
SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : SpatialLayout(definition, spatial_detail::resolve_block_ids(definition, source_mesh, false, false),
          spatial_detail::DofLayout::cartesian_3d) {
    _meshes.reserve(_block_ids.size());
    for (const std::int64_t block_id : _block_ids)
        _meshes.push_back(Hex8RegionMesh::from_unstructured_block(source_mesh, block_id));
    std::vector<std::size_t> node_counts, element_counts;
    node_counts.reserve(_meshes.size());
    element_counts.reserve(_meshes.size());
    for (const Hex8RegionMesh& mesh : _meshes) {
        node_counts.push_back(mesh.nodes().size());
        element_counts.push_back(mesh.elements().size());
    }
    initialize_counts(node_counts, element_counts);
    _geometries.resize(_meshes.size());
    for (std::size_t region = 0; region < _meshes.size(); ++region) {
        for (const Hex8Element& element : _meshes[region].elements()) {
            Hex8Coordinates coordinates{};
            for (std::size_t node = 0; node < coordinates.size(); ++node)
                coordinates[node] = _meshes[region].nodes().at(element.nodes[node]);
            _geometries[region].push_back(make_hex8_geometry(coordinates));
        }
    }
    for (std::size_t boundary_index = 0; boundary_index < _definition.boundary_conditions.size(); ++boundary_index) {
        const BoundaryConditionDefinition& boundary = _definition.boundary_conditions[boundary_index];
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
                const std::size_t dof = this->dof(boundary.field, global_node(region, local_node));
                add_dirichlet(dof, boundary_index);
            }
            continue;
        }
        if (boundary.use_displaced_geometry)
            throw std::invalid_argument("Cartesian three-dimensional stage B loads use the reference configuration");
        const std::size_t kernel = _boundary_data.size();
        if (boundary.type == BoundaryConditionType::pressure) {
            _boundary_data.push_back(
                {Quad4FaceBoundaryKind::pressure, CartesianTractionComponent::x, boundary.value, 0.0});
        } else if (boundary.type == BoundaryConditionType::traction) {
            CartesianTractionComponent component = CartesianTractionComponent::x;
            if (boundary.field == Field::displacement_y)
                component = CartesianTractionComponent::y;
            else if (boundary.field == Field::displacement_z)
                component = CartesianTractionComponent::z;
            else if (boundary.field != Field::displacement_x)
                throw std::invalid_argument("Three-dimensional traction requires a displacement field");
            _boundary_data.push_back({Quad4FaceBoundaryKind::traction, component, boundary.value, 0.0});
        } else {
            _boundary_data.push_back({Quad4FaceBoundaryKind::convection, CartesianTractionComponent::x,
                boundary.heat_transfer_coefficient, boundary.ambient_temperature});
        }
        _boundary_definition_indices.push_back(boundary_index);
        for (const Quad4FaceElement& face : mapped.faces) {
            std::array<std::size_t, 4> nodes{};
            Quad4FaceCoordinates coordinates{};
            for (std::size_t node = 0; node < 4; ++node) nodes[node] = global_node(region, face.nodes[node]);
            for (std::size_t node = 0; node < 4; ++node)
                coordinates[node] = _meshes[region].nodes().at(face.nodes[node]);
            _boundary_contributions.push_back(
                {boundary.type == BoundaryConditionType::pressure
                        ? SpatialContributionType::pressure
                        : (boundary.type == BoundaryConditionType::traction ? SpatialContributionType::traction
                                                                            : SpatialContributionType::convection),
                    kernel, nodes, make_quad4_face_geometry(coordinates)});
        }
    }
    spatial_detail::validate_dirichlet_conditions(_dirichlet_conditions,
        "Cartesian three-dimensional boundary has conflicting Dirichlet values",
        "Cartesian three-dimensional boundary has duplicate Dirichlet values");
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index) {
        const MaterialFunctionSet& functions = *region(region_index).material.functions;
        if (functions.has_creep() || functions.has_plasticity())
            throw std::invalid_argument("Cartesian three-dimensional stage B supports only elastic materials");
        _kernel_data.push_back(
            {IsotropicThermoelasticMaterial(region(region_index).material), region_heat_source(region_index), 0.0});
    }
    refresh_controls();
}
SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    if (index < volume_contribution_count()) return SpatialContributionType::volume;
    return _boundary_contributions.at(index - volume_contribution_count()).type;
}
const Hex8Geometry& SpatialAssembly::region_element_geometry(std::size_t region, std::size_t element_index) const {
    return _geometries.at(region).at(element_index);
}
void SpatialAssembly::set_load_factor(double value) {
    set_load_factor_value(value);
    refresh_controls();
}
void SpatialAssembly::set_time(double value) {
    set_time_value(value);
    for (Hex8ThermoelasticData& kernel_data : _kernel_data) kernel_data.time = value;
    refresh_controls();
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
            fixed[node] = dof(Field::temperature, nodes[node]);
            fixed[8 + node] = dof(Field::displacement_x, nodes[node]);
            fixed[16 + node] = dof(Field::displacement_y, nodes[node]);
            fixed[24 + node] = dof(Field::displacement_z, nodes[node]);
        }
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - volume_contribution_count());
    Quad4FaceLocalDofs fixed{};
    for (std::size_t node = 0; node < entry.nodes.size(); ++node) {
        fixed[node] = dof(Field::temperature, entry.nodes[node]);
        fixed[4 + node] = dof(Field::displacement_x, entry.nodes[node]);
        fixed[8 + node] = dof(Field::displacement_y, entry.nodes[node]);
        fixed[12 + node] = dof(Field::displacement_z, entry.nodes[node]);
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
void SpatialAssembly::compute_contribution(std::size_t index, const std::vector<double>& state,
    const std::vector<double>* committed_solution, double time_step, std::vector<double>& residual,
    std::vector<double>* jacobian) const {
    if (index < volume_contribution_count()) {
        if (state.size() != hex8_local_dof_count)
            throw std::invalid_argument("HEX8 contribution state must contain 32 DOFs");
        const auto location = element_location(index);
        Hex8LocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        const Hex8LocalValues committed =
            committed_solution == nullptr ? Hex8LocalValues{} : volume_state(index, *committed_solution);
        Hex8LocalJacobian local_jacobian{};
        const Hex8LocalResidual result = compute_hex8_thermoelastic(_kernel_data[location.first],
            region_element_geometry(location.first, location.second), current,
            committed_solution == nullptr ? nullptr : &committed, time_step,
            jacobian == nullptr ? nullptr : &local_jacobian);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    if (state.size() != quad4_face_local_dof_count)
        throw std::invalid_argument("Three-dimensional face state must contain 16 DOFs");
    const BoundaryContribution& entry = _boundary_contributions.at(index - volume_contribution_count());
    Quad4FaceLocalValues current{};
    std::copy(state.begin(), state.end(), current.begin());
    Quad4FaceLocalJacobian local_jacobian{};
    const Quad4FaceLocalResidual result = compute_quad4_face_boundary(
        _boundary_data[entry.kernel], entry.geometry, current, jacobian == nullptr ? nullptr : &local_jacobian);
    residual.assign(result.begin(), result.end());
    if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
}
std::array<SymmetricTensor3Values, 8> SpatialAssembly::stress(
    std::size_t region, std::size_t element, const std::vector<double>& state) const {
    return compute_hex8_stress(_kernel_data.at(region), region_element_geometry(region, element),
        volume_state(region_element_offset(region) + element, state));
}
double SpatialAssembly::heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const {
    const Hex8ThermoelasticData& data = _kernel_data.at(region);
    return data.material.heat_capacity(temperature, {data.time, position.x, position.y, position.z}).value();
}
void SpatialAssembly::refresh_controls() {
    for (std::size_t region = 0; region < region_count(); ++region)
        _kernel_data[region].volumetric_heat_source = region_heat_source(region);
    refresh_dirichlet_values();
    for (std::size_t kernel = 0; kernel < _boundary_data.size(); ++kernel) {
        const BoundaryConditionDefinition& boundary =
            _definition.boundary_conditions[_boundary_definition_indices[kernel]];
        if (boundary.type != BoundaryConditionType::convection) {
            _boundary_data[kernel].load = spatial_detail::controlled_value(
                _definition, _time, _load_factor, boundary.value, boundary.scale_with_load, boundary.function);
            continue;
        }
        const spatial_detail::ConvectionValues values = convection_values(boundary);
        _boundary_data[kernel].load = values.coefficient;
        _boundary_data[kernel].ambient_temperature = values.ambient;
    }
}
} // namespace fuelsim::cartesian
