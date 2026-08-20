#include "cartesian3d_assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace fuelsim::cartesian {
namespace {
CartesianPoint3 face_centroid(const Quad4FaceCoordinates& coordinates) {
    CartesianPoint3 result{0.0, 0.0, 0.0};
    for (const CartesianPoint3& point : coordinates) {
        result.x += 0.25 * point.x;
        result.y += 0.25 * point.y;
        result.z += 0.25 * point.z;
    }
    return result;
}

CartesianPoint3 element_centroid(const Hex8RegionMesh& mesh, std::size_t element_index) {
    CartesianPoint3 result{0.0, 0.0, 0.0};
    for (std::size_t node : mesh.elements().at(element_index).nodes) {
        result.x += 0.125 * mesh.nodes().at(node).x;
        result.y += 0.125 * mesh.nodes().at(node).y;
        result.z += 0.125 * mesh.nodes().at(node).z;
    }
    return result;
}

Quad4FaceCoordinates face_coordinates(const Hex8RegionMesh& mesh, const Quad4FaceElement& face) {
    Quad4FaceCoordinates result{};
    for (std::size_t node = 0; node < 4; ++node) result[node] = mesh.nodes().at(face.nodes[node]);
    return result;
}

CartesianPoint3 subtract(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

CartesianPoint3 cross(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y, first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

double dot(const CartesianPoint3& first, const CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

std::array<double, 3> thermal_search_point(
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    std::array<double, 3> result{};
    for (std::size_t node = 0; node < 4; ++node) {
        result[0] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].x + state[8 + node]);
        result[1] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].y + state[16 + node]);
        result[2] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].z + state[24 + node]);
    }
    return result;
}

std::array<double, 3> mechanical_search_point(
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    const std::size_t node = geometry.secondary_local_node;
    return {geometry.secondary_coordinates[node].x + state[8 + node],
        geometry.secondary_coordinates[node].y + state[16 + node],
        geometry.secondary_coordinates[node].z + state[24 + node]};
}

double normal_orientation(const Quad4FaceCoordinates& primary_coordinates,
    const CartesianPoint3& secondary_parent_centroid, const CartesianPoint3& primary_parent_centroid) {
    const CartesianPoint3 tangent_xi = subtract(primary_coordinates[1], primary_coordinates[0]),
                          tangent_eta = subtract(primary_coordinates[3], primary_coordinates[0]),
                          area = cross(tangent_xi, tangent_eta),
                          material_direction = subtract(primary_parent_centroid, secondary_parent_centroid);
    const double measure = std::sqrt(dot(area, area)), orientation = dot(area, material_direction);
    if (!std::isfinite(measure) || !(measure > 0.0) || !std::isfinite(orientation) || orientation == 0.0)
        throw std::invalid_argument(
            "Three-dimensional contact faces require nondegenerate opposing material-side centroids");
    return orientation > 0.0 ? 1.0 : -1.0;
}

double minimum_normal_length(const Hex8RegionMesh& mesh, const Hex8RegionBoundary& boundary, const std::string& name) {
    double result = std::numeric_limits<double>::infinity();
    for (const Quad4FaceElement& face : boundary.faces) {
        const Quad4FaceCoordinates coordinates = face_coordinates(mesh, face);
        const CartesianPoint3 face_center = face_centroid(coordinates),
                              parent = element_centroid(mesh, face.parent_element),
                              first = subtract(coordinates[1], coordinates[0]),
                              second = subtract(coordinates[3], coordinates[0]), area = cross(first, second);
        const double measure = std::sqrt(dot(area, area));
        if (!std::isfinite(measure) || !(measure > 0.0))
            throw std::invalid_argument("Three-dimensional contact boundary has a degenerate face: " + name);
        const double height = 2.0 * std::abs(dot(subtract(parent, face_center), area)) / measure;
        if (!std::isfinite(height) || !(height > 0.0))
            throw std::invalid_argument("Three-dimensional contact boundary has a nonpositive normal length: " + name);
        result = std::min(result, height);
    }
    if (!std::isfinite(result)) throw std::invalid_argument("Three-dimensional contact boundary is empty: " + name);
    return result;
}

void set_pattern_block(std::vector<unsigned char>& pattern, std::size_t size, std::size_t row_begin,
    std::size_t row_end, std::size_t column_begin, std::size_t column_end) {
    for (std::size_t row = row_begin; row < row_end; ++row)
        std::fill(pattern.begin() + static_cast<std::ptrdiff_t>(row * size + column_begin),
            pattern.begin() + static_cast<std::ptrdiff_t>(row * size + column_end), 1U);
}

std::pair<std::size_t, std::size_t> offset_location(
    const std::vector<std::size_t>& offsets, std::size_t index, const char* message) {
    if (offsets.empty() || index >= offsets.back()) throw std::out_of_range(message);
    const auto upper = std::upper_bound(offsets.begin(), offsets.end(), index);
    const std::size_t group = static_cast<std::size_t>(upper - offsets.begin() - 1);
    return {group, index - offsets[group]};
}
} // namespace

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : SpatialLayout(definition, spatial_detail::resolve_block_ids(definition, source_mesh, true, true),
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
    std::vector<std::vector<std::size_t>> region_source_node_ids;
    region_source_node_ids.reserve(_meshes.size());
    for (const Hex8RegionMesh& mesh : _meshes) region_source_node_ids.push_back(mesh.source_node_ids());
    initialize_shared_nodes(region_source_node_ids);
    _geometries.resize(_meshes.size());
    for (std::size_t region = 0; region < _meshes.size(); ++region) {
        for (const Hex8Element& element : _meshes[region].elements()) {
            Hex8Coordinates coordinates{};
            for (std::size_t node = 0; node < coordinates.size(); ++node)
                coordinates[node] = _meshes[region].nodes().at(element.nodes[node]);
            _geometries[region].push_back(make_hex8_geometry(coordinates));
        }
    }
    build_contacts(source_mesh);
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
        const bool finite_strain = this->region(region).strain_formulation == StrainFormulation::finite;
        if (boundary.use_displaced_geometry && !finite_strain)
            throw std::invalid_argument(
                "Cartesian current-configuration loads require a finite-strain region: " + boundary.name);
        if (finite_strain && boundary.type == BoundaryConditionType::pressure && !boundary.use_displaced_geometry)
            throw std::invalid_argument(
                "Cartesian finite-strain pressure requires configuration = current: " + boundary.name);
        const std::size_t kernel = _boundary_data.size();
        if (boundary.type == BoundaryConditionType::pressure) {
            _boundary_data.push_back({Quad4FaceBoundaryKind::pressure, CartesianTractionComponent::x, boundary.value,
                0.0, boundary.use_displaced_geometry});
        } else if (boundary.type == BoundaryConditionType::traction) {
            CartesianTractionComponent component = CartesianTractionComponent::x;
            if (boundary.field == Field::displacement_y)
                component = CartesianTractionComponent::y;
            else if (boundary.field == Field::displacement_z)
                component = CartesianTractionComponent::z;
            else if (boundary.field != Field::displacement_x)
                throw std::invalid_argument("Three-dimensional traction requires a displacement field");
            _boundary_data.push_back(
                {Quad4FaceBoundaryKind::traction, component, boundary.value, 0.0, boundary.use_displaced_geometry});
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
        _kernel_data.push_back({IsotropicThermoelasticMaterial(region(region_index).material),
            region_heat_source(region_index), 0.0, region(region_index).strain_formulation});
    }
    _committed_contact_solution = initial_state();
    validate_local_state(0, contribution_count(), _committed_contact_solution);
    refresh_controls();
}

SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    const ContributionRanges ranges = contribution_ranges();
    if (index >= ranges.end) throw std::out_of_range("Three-dimensional contribution index is out of range");
    if (index < ranges.thermal_begin) return SpatialContributionType::volume;
    if (index < ranges.mechanical_begin) return SpatialContributionType::thermal_contact;
    if (index < ranges.boundary_begin) return SpatialContributionType::mechanical_contact;
    return _boundary_contributions.at(index - ranges.boundary_begin).type;
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
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index) {
        if (region(region_index).strain_formulation != StrainFormulation::finite) continue;
        for (std::size_t element = 0; element < region_element_count(region_index); ++element) {
            const Hex8LocalValues local = volume_state(region_element_offset(region_index) + element, state);
            for (const Hex8QuadraturePoint& point : region_element_geometry(region_index, element).points)
                validate_cartesian_deformation(point, local);
        }
    }
    validate_local_state(0, contribution_count(), state);
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
    const ContributionRanges ranges = contribution_ranges();
    if (index < ranges.mechanical_begin) {
        const std::size_t point = index - ranges.thermal_begin;
        const Quad4SurfaceContactLocalDofs fixed =
            contact_dofs(thermal_candidate(point, _thermal_active_primary.at(point)).nodes);
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    if (index < ranges.boundary_begin) {
        const std::size_t point = index - ranges.mechanical_begin;
        const Quad4SurfaceContactLocalDofs fixed =
            contact_dofs(mechanical_candidate(point, _mechanical_active_primary.at(point)).nodes);
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - ranges.boundary_begin);
    Quad4FaceLocalDofs fixed{};
    for (std::size_t node = 0; node < entry.nodes.size(); ++node) {
        fixed[node] = dof(Field::temperature, entry.nodes[node]);
        fixed[4 + node] = dof(Field::displacement_x, entry.nodes[node]);
        fixed[8 + node] = dof(Field::displacement_y, entry.nodes[node]);
        fixed[12 + node] = dof(Field::displacement_z, entry.nodes[node]);
    }
    dofs.assign(fixed.begin(), fixed.end());
}

void SpatialAssembly::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    const ContributionRanges ranges = contribution_ranges();
    if (index >= ranges.end) throw std::out_of_range("Three-dimensional contribution index is out of range");
    if (index < ranges.thermal_begin) {
        pattern.assign(hex8_local_dof_count * hex8_local_dof_count, 0U);
        set_pattern_block(pattern, hex8_local_dof_count, 0, 8, 0, 8);
        set_pattern_block(pattern, hex8_local_dof_count, 8, hex8_local_dof_count, 0, hex8_local_dof_count);
        return;
    }
    if (index < ranges.mechanical_begin) {
        pattern.assign(quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count, 0U);
        set_pattern_block(
            pattern, quad4_surface_contact_local_dof_count, 0, 8, 0, quad4_surface_contact_local_dof_count);
        return;
    }
    if (index < ranges.boundary_begin) {
        pattern.assign(quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count, 0U);
        set_pattern_block(pattern, quad4_surface_contact_local_dof_count, 8, quad4_surface_contact_local_dof_count, 8,
            quad4_surface_contact_local_dof_count);
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - ranges.boundary_begin);
    const Quad4FaceBoundaryData& data = _boundary_data.at(entry.kernel);
    pattern.assign(quad4_face_local_dof_count * quad4_face_local_dof_count, 0U);
    if (data.kind == Quad4FaceBoundaryKind::convection) {
        set_pattern_block(pattern, quad4_face_local_dof_count, 0, 4, 0, 4);
    } else if (data.use_displaced_geometry) {
        const std::size_t row_begin = data.kind == Quad4FaceBoundaryKind::pressure
                                          ? 4
                                          : (data.component == CartesianTractionComponent::x
                                                    ? 4
                                                    : (data.component == CartesianTractionComponent::y ? 8 : 12));
        const std::size_t row_end = data.kind == Quad4FaceBoundaryKind::pressure ? 16 : row_begin + 4;
        set_pattern_block(pattern, quad4_face_local_dof_count, row_begin, row_end, 4, 16);
    }
}

void SpatialAssembly::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    if (index < volume_contribution_count()) {
        contribution_dofs(index, dofs);
        return;
    }
    index -= volume_contribution_count();
    const Quad4SurfaceContactLocalDofs fixed = contact_dofs(sparsity_contact(index).nodes);
    dofs.assign(fixed.begin(), fixed.end());
}

void SpatialAssembly::sparsity_contribution_jacobian_pattern(
    std::size_t index, std::vector<unsigned char>& pattern) const {
    if (index < volume_contribution_count()) {
        contribution_jacobian_pattern(index, pattern);
        return;
    }
    index -= volume_contribution_count();
    const SparsityContact contact = sparsity_contact(index);
    pattern.assign(quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count, 0U);
    if (contact.thermal)
        set_pattern_block(
            pattern, quad4_surface_contact_local_dof_count, 0, 8, 0, quad4_surface_contact_local_dof_count);
    if (contact.mechanical)
        set_pattern_block(pattern, quad4_surface_contact_local_dof_count, 8, quad4_surface_contact_local_dof_count, 8,
            quad4_surface_contact_local_dof_count);
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
    const std::vector<double>* committed_solution, const Hex8MaterialHistory* committed_material, double time_step,
    std::vector<double>& residual, std::vector<double>* jacobian, bool include_thermal_time_term) const {
    if (index < volume_contribution_count()) {
        if (state.size() != hex8_local_dof_count)
            throw std::invalid_argument("HEX8 contribution state must contain 32 DOFs");
        const auto location = element_location(index);
        Hex8LocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        const Hex8LocalValues committed =
            committed_solution == nullptr ? Hex8LocalValues{} : volume_state(index, *committed_solution);
        Hex8LocalJacobian local_jacobian{};
        const Hex8LocalResidual result =
            committed_material == nullptr
                ? compute_hex8_thermoelastic(_kernel_data[location.first],
                      region_element_geometry(location.first, location.second), current,
                      committed_solution == nullptr ? nullptr : &committed, time_step,
                      jacobian == nullptr ? nullptr : &local_jacobian)
                : compute_hex8_transient(_kernel_data[location.first],
                      region_element_geometry(location.first, location.second), current, committed, *committed_material,
                      time_step, jacobian == nullptr ? nullptr : &local_jacobian, include_thermal_time_term);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    const ContributionRanges ranges = contribution_ranges();
    if (index < ranges.mechanical_begin) {
        if (state.size() != quad4_surface_contact_local_dof_count)
            throw std::invalid_argument("Three-dimensional thermal-contact state must contain 32 DOFs");
        const std::size_t point = index - ranges.thermal_begin;
        const ThermalCandidate entry = thermal_candidate(point, _thermal_active_primary.at(point));
        Quad4SurfaceContactLocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        Quad4SurfaceContactLocalJacobian local_jacobian{};
        const Quad4SurfaceContactLocalResidual result =
            compute_quad4_to_quad4_gap_heat(_thermal_properties[entry.contact], entry.geometry, current,
                jacobian == nullptr ? nullptr : &local_jacobian);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    if (index < ranges.boundary_begin) {
        if (state.size() != quad4_surface_contact_local_dof_count)
            throw std::invalid_argument("Three-dimensional mechanical-contact state must contain 32 DOFs");
        const std::size_t point = index - ranges.mechanical_begin;
        const MechanicalCandidate entry = mechanical_candidate(point, _mechanical_active_primary.at(point));
        Quad4SurfaceContactLocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        const Quad4SurfaceContactLocalValues committed = contribution_state(index, _committed_contact_solution);
        Quad4SurfaceContactLocalJacobian local_jacobian{};
        const Quad4SurfaceContactLocalResidual result =
            compute_node_to_quad4_contact(_mechanical_properties[entry.contact], entry.geometry, current, committed,
                _contact_histories[entry.contact][entry.secondary], jacobian == nullptr ? nullptr : &local_jacobian);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    if (state.size() != quad4_face_local_dof_count)
        throw std::invalid_argument("Three-dimensional face state must contain 16 DOFs");
    const BoundaryContribution& entry = _boundary_contributions.at(index - ranges.boundary_begin);
    Quad4FaceLocalValues current{};
    std::copy(state.begin(), state.end(), current.begin());
    Quad4FaceLocalJacobian local_jacobian{};
    const Quad4FaceLocalResidual result = compute_quad4_face_boundary(
        _boundary_data[entry.kernel], entry.geometry, current, jacobian == nullptr ? nullptr : &local_jacobian);
    residual.assign(result.begin(), result.end());
    if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
}

Hex8MaterialHistory SpatialAssembly::transient_update(std::size_t region, std::size_t element,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material,
    double time_step) const {
    return compute_hex8_transient_update(_kernel_data.at(region), region_element_geometry(region, element), state,
        committed_state, committed_material, time_step);
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

SpatialAssembly::ContributionRanges SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count(),
                      mechanical_begin = thermal_begin + _thermal_contact_offsets.back(),
                      boundary_begin = mechanical_begin + _mechanical_points.size();
    return {thermal_begin, mechanical_begin, boundary_begin, boundary_begin + _boundary_contributions.size()};
}

std::size_t SpatialAssembly::sparsity_contribution_count() const noexcept {
    // Every face boundary block is a subset of its adjacent volume block. Contact quadrature points and face nodes
    // that share one secondary-face/primary-face pair also have the same 32-DOF graph, so one representative preserves
    // the complete graph without repeating it for each runtime contribution.
    return volume_contribution_count() + _sparsity_contact_offsets.back();
}

std::size_t SpatialAssembly::contribution_work(std::size_t index, std::size_t partition_count) const {
    // These relative units track the measured AD-local work of eight-point finite-strain volume integration and the
    // three surface kernels. They affect only the contiguous MPI ownership boundary, never the residual or Jacobian.
    const ContributionRanges ranges = contribution_ranges();
    if (index < ranges.thermal_begin) return 64;
    if (partition_count < 4) {
        if (index < ranges.mechanical_begin) return 10;
        if (index < ranges.boundary_begin) return 11;
        if (index < ranges.end) return 1;
    } else {
        // With four or more partitions the last contiguous interval owns every surface contribution. Account for
        // residual work as well as AD Jacobian work so that this interval does not become the synchronization tail.
        if (index < ranges.mechanical_begin) return 13;
        if (index < ranges.boundary_begin) return 14;
        if (index < ranges.end) return 2;
    }
    throw std::out_of_range("Three-dimensional contribution work index is out of range");
}

std::pair<std::size_t, std::size_t> SpatialAssembly::contribution_partition(
    std::size_t partition, std::size_t partition_count) const {
    if (partition_count == 0 || partition >= partition_count)
        throw std::out_of_range("Three-dimensional contribution partition is out of range");
    const std::size_t count = contribution_count();
    std::size_t total_work = 0;
    for (std::size_t entry = 0; entry < count; ++entry) {
        const std::size_t work = contribution_work(entry, partition_count);
        if (total_work > std::numeric_limits<std::size_t>::max() - work)
            throw std::overflow_error("Three-dimensional contribution work exceeds size_t range");
        total_work += work;
    }
    const auto boundary = [&](std::size_t boundary_partition) {
        if (boundary_partition == 0) return std::size_t{0};
        if (boundary_partition == partition_count) return count;
        const std::size_t target = (total_work / partition_count) * boundary_partition +
                                   ((total_work % partition_count) * boundary_partition) / partition_count;
        std::size_t accumulated = 0;
        for (std::size_t entry = 0; entry < count; ++entry) {
            const std::size_t next = accumulated + contribution_work(entry, partition_count);
            if (next >= target) return target - accumulated < next - target ? entry : entry + 1U;
            accumulated = next;
        }
        return count;
    };
    return {boundary(partition), boundary(partition + 1U)};
}

ResolvedBoundary SpatialAssembly::resolve_boundary(
    const UnstructuredHex8Mesh& source_mesh, const std::string& name) const {
    const std::int64_t block_id = source_mesh.side_set_block_id(name);
    const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
    if (found == _block_ids.end())
        throw std::invalid_argument("Contact boundary belongs to an undeclared three-dimensional block: " + name);
    const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
    Hex8RegionBoundary boundary = _meshes[region].map_side_set(source_mesh, name);
    if (boundary.faces.empty()) throw std::invalid_argument("Three-dimensional contact side set is empty: " + name);
    return {region, std::move(boundary)};
}

Quad4SurfaceContactLocalDofs SpatialAssembly::contact_dofs(const std::array<std::size_t, 8>& nodes) const {
    Quad4SurfaceContactLocalDofs result{};
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        result[node] = dof(Field::temperature, nodes[node]);
        result[8 + node] = dof(Field::displacement_x, nodes[node]);
        result[16 + node] = dof(Field::displacement_y, nodes[node]);
        result[24 + node] = dof(Field::displacement_z, nodes[node]);
    }
    return result;
}

Quad4SurfaceContactLocalValues SpatialAssembly::contribution_state(
    std::size_t index, const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional contact global state has the wrong size");
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    if (dofs.size() != quad4_surface_contact_local_dof_count)
        throw std::logic_error("Three-dimensional contact contribution has an invalid DOF layout");
    Quad4SurfaceContactLocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

Quad4SurfaceContactLocalValues SpatialAssembly::contact_state(
    const std::array<std::size_t, 8>& nodes, const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional contact state has the wrong size");
    const Quad4SurfaceContactLocalDofs dofs = contact_dofs(nodes);
    Quad4SurfaceContactLocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

void SpatialAssembly::build_contacts(const UnstructuredHex8Mesh& source_mesh) {
    _primary_boundaries.reserve(_definition.contacts.size());
    _secondary_boundaries.reserve(_definition.contacts.size());
    _primary_contact_faces.reserve(_definition.contacts.size());
    _secondary_contact_faces.reserve(_definition.contacts.size());
    _thermal_properties.reserve(_definition.contacts.size());
    _mechanical_properties.reserve(_definition.contacts.size());
    _thermal_point_counts.reserve(_definition.contacts.size());
    _contact_histories.resize(_definition.contacts.size());
    for (std::size_t contact_value = 0; contact_value < _definition.contacts.size(); ++contact_value) {
        ContactDefinition& definition = _definition.contacts[contact_value];
        ResolvedBoundary primary = resolve_boundary(source_mesh, definition.primary),
                         secondary = resolve_boundary(source_mesh, definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Three-dimensional self-contact is not supported: " + definition.name);
        const Hex8RegionMesh& primary_mesh = _meshes[primary.region];
        const Hex8RegionMesh& secondary_mesh = _meshes[secondary.region];
        for (const std::size_t secondary_local : secondary.boundary.nodes) {
            const std::size_t source_node = secondary_mesh.source_node_ids().at(secondary_local);
            if (std::any_of(
                    primary.boundary.nodes.begin(), primary.boundary.nodes.end(), [&](std::size_t primary_local) {
                        return primary_mesh.source_node_ids().at(primary_local) == source_node;
                    }))
                throw std::invalid_argument(
                    "Three-dimensional contact boundaries must not share source nodes: " + definition.name);
        }
        if (definition.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian)
            throw std::invalid_argument(
                "Three-dimensional contact currently supports the penalty formulation only: " + definition.name);
        if (definition.mechanical && definition.automatic_penalty) {
            const double primary_length = minimum_normal_length(primary_mesh, primary.boundary, definition.primary),
                         secondary_length =
                             minimum_normal_length(secondary_mesh, secondary.boundary, definition.secondary),
                         primary_modulus = _definition.regions[primary.region].material.reference_young_modulus,
                         secondary_modulus = _definition.regions[secondary.region].material.reference_young_modulus;
            definition.penalty =
                definition.penalty_factor / (primary_length / primary_modulus + secondary_length / secondary_modulus);
            if (!std::isfinite(definition.penalty) || !(definition.penalty > 0.0))
                throw std::overflow_error(
                    "Automatic three-dimensional contact penalty is not finite and positive: " + definition.name);
        }
        _thermal_properties.push_back({definition.thermal ? definition.gap_conductivity : 1.0,
            definition.thermal ? definition.minimum_gap : 1.0});
        _mechanical_properties.push_back({definition.mechanical ? definition.penalty : 1.0,
            definition.mechanical ? definition.friction_coefficient : 0.0, false});
        _contact_histories[contact_value].resize(secondary.boundary.nodes.size());
        std::vector<PrimaryContactFace> primary_faces;
        primary_faces.reserve(primary.boundary.faces.size());
        for (const Quad4FaceElement& primary_face : primary.boundary.faces) {
            const Quad4FaceCoordinates coordinates = face_coordinates(primary_mesh, primary_face);
            std::array<std::size_t, 4> nodes{};
            for (std::size_t node = 0; node < nodes.size(); ++node)
                nodes[node] = global_node(primary.region, primary_face.nodes[node]);
            primary_faces.push_back({nodes, coordinates, element_centroid(primary_mesh, primary_face.parent_element)});
        }
        std::vector<SecondaryContactFace> secondary_faces;
        secondary_faces.reserve(secondary.boundary.faces.size());
        std::size_t thermal_point_count = 0;
        for (const Quad4FaceElement& secondary_face : secondary.boundary.faces) {
            const Quad4FaceCoordinates secondary_coordinates = face_coordinates(secondary_mesh, secondary_face);
            const Quad4FaceGeometry secondary_geometry = make_quad4_face_geometry(secondary_coordinates);
            const CartesianPoint3 secondary_parent = element_centroid(secondary_mesh, secondary_face.parent_element);
            std::array<std::size_t, 4> secondary_nodes{};
            for (std::size_t node = 0; node < 4; ++node)
                secondary_nodes[node] = global_node(secondary.region, secondary_face.nodes[node]);
            const std::size_t secondary_face_index = secondary_faces.size();
            secondary_faces.push_back({secondary_nodes, secondary_coordinates, secondary_geometry, secondary_parent});
            if (definition.thermal) thermal_point_count += secondary_geometry.points.size();
            if (definition.mechanical)
                for (std::size_t secondary_local_node = 0; secondary_local_node < 4; ++secondary_local_node) {
                    const auto found = std::find(secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(),
                        secondary_face.nodes[secondary_local_node]);
                    if (found == secondary.boundary.nodes.end())
                        throw std::logic_error("Three-dimensional secondary contact node mapping failed");
                    const std::size_t secondary_node_index =
                        static_cast<std::size_t>(found - secondary.boundary.nodes.begin());
                    _mechanical_points.push_back(
                        {contact_value, secondary_node_index, secondary_face_index, secondary_local_node});
                }
        }
        _thermal_point_counts.push_back(thermal_point_count);
        _primary_contact_faces.push_back(std::move(primary_faces));
        _secondary_contact_faces.push_back(std::move(secondary_faces));
        _primary_boundaries.push_back(std::move(primary));
        _secondary_boundaries.push_back(std::move(secondary));
    }
    _thermal_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    _mechanical_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    _sparsity_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        _thermal_contact_offsets[contact + 1] = _thermal_contact_offsets[contact] + _thermal_point_counts[contact];
        _mechanical_contact_offsets[contact + 1] =
            _mechanical_contact_offsets[contact] + _contact_histories[contact].size();
        const std::size_t primary_count = _primary_contact_faces[contact].size(),
                          secondary_count = _secondary_contact_faces[contact].size();
        if (secondary_count > std::numeric_limits<std::size_t>::max() / primary_count)
            throw std::overflow_error("Three-dimensional contact sparsity candidate count exceeds size_t range");
        const std::size_t pairs = (_definition.contacts[contact].thermal || _definition.contacts[contact].mechanical)
                                      ? secondary_count * primary_count
                                      : 0;
        if (_sparsity_contact_offsets[contact] > std::numeric_limits<std::size_t>::max() - pairs)
            throw std::overflow_error("Three-dimensional contact sparsity candidate offset exceeds size_t range");
        _sparsity_contact_offsets[contact + 1] = _sparsity_contact_offsets[contact] + pairs;
    }
    _touched_thermal_points.resize(_thermal_contact_offsets.back());
    _thermal_minimum_distance.resize(_thermal_contact_offsets.back());
    _thermal_active_primary.resize(_thermal_contact_offsets.back());
    _thermal_cached_primary.assign(_thermal_contact_offsets.back(), std::numeric_limits<std::size_t>::max());
    _touched_mechanical_nodes.resize(_mechanical_contact_offsets.back());
    _mechanical_minimum_distance.resize(_mechanical_contact_offsets.back());
    _mechanical_selected_primary.resize(_mechanical_contact_offsets.back());
    _mechanical_cached_primary.assign(_mechanical_contact_offsets.back(), std::numeric_limits<std::size_t>::max());
    _mechanical_active_primary.resize(_mechanical_points.size());
    _contact_search_trees.resize(_definition.contacts.size());
}

SpatialAssembly::ThermalCandidate SpatialAssembly::thermal_candidate(std::size_t point, std::size_t primary) const {
    const auto location =
        offset_location(_thermal_contact_offsets, point, "Three-dimensional thermal-contact point is out of range");
    const std::size_t contact = location.first, local_point = location.second, secondary_face = local_point / 4,
                      quadrature_point = local_point % 4;
    const SecondaryContactFace& secondary = _secondary_contact_faces.at(contact).at(secondary_face);
    const PrimaryContactFace& primary_face = _primary_contact_faces.at(contact).at(primary);
    std::array<std::size_t, 8> nodes{};
    std::copy(secondary.nodes.begin(), secondary.nodes.end(), nodes.begin());
    std::copy(primary_face.nodes.begin(), primary_face.nodes.end(), nodes.begin() + 4);
    const Quad4FaceQuadraturePoint& quadrature = secondary.geometry.points[quadrature_point];
    return {contact, nodes,
        {secondary.coordinates, primary_face.coordinates, quadrature.shape, quadrature.derivative_xi,
            quadrature.derivative_eta,
            normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid)},
        primary};
}

SpatialAssembly::MechanicalCandidate SpatialAssembly::mechanical_candidate(
    std::size_t point, std::size_t primary) const {
    const MechanicalPoint& metadata = _mechanical_points.at(point);
    const SecondaryContactFace& secondary = _secondary_contact_faces.at(metadata.contact).at(metadata.secondary_face);
    const PrimaryContactFace& primary_face = _primary_contact_faces.at(metadata.contact).at(primary);
    std::array<std::size_t, 8> nodes{};
    std::copy(secondary.nodes.begin(), secondary.nodes.end(), nodes.begin());
    std::copy(primary_face.nodes.begin(), primary_face.nodes.end(), nodes.begin() + 4);
    std::array<std::array<double, 4>, 4> shapes{}, derivatives_xi{}, derivatives_eta{};
    for (std::size_t q = 0; q < secondary.geometry.points.size(); ++q) {
        shapes[q] = secondary.geometry.points[q].shape;
        derivatives_xi[q] = secondary.geometry.points[q].derivative_xi;
        derivatives_eta[q] = secondary.geometry.points[q].derivative_eta;
    }
    return {metadata.contact, nodes,
        {secondary.coordinates, primary_face.coordinates, shapes, derivatives_xi, derivatives_eta,
            metadata.secondary_local_node,
            normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid)},
        metadata.secondary, primary};
}

SpatialAssembly::SparsityContact SpatialAssembly::sparsity_contact(std::size_t index) const {
    const auto location = offset_location(
        _sparsity_contact_offsets, index, "Three-dimensional sparsity contribution index is out of range");
    const std::size_t contact = location.first, primary_count = _primary_contact_faces.at(contact).size(),
                      secondary_face = location.second / primary_count, primary = location.second % primary_count;
    std::array<std::size_t, 8> nodes{};
    const std::array<std::size_t, 4>& secondary_nodes = _secondary_contact_faces.at(contact).at(secondary_face).nodes;
    const std::array<std::size_t, 4>& primary_nodes = _primary_contact_faces.at(contact).at(primary).nodes;
    std::copy(secondary_nodes.begin(), secondary_nodes.end(), nodes.begin());
    std::copy(primary_nodes.begin(), primary_nodes.end(), nodes.begin() + 4);
    return {nodes, _definition.contacts[contact].thermal, _definition.contacts[contact].mechanical};
}

void SpatialAssembly::update_contact_search_trees(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional contact-search tree state size mismatch");
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        const ResolvedBoundary& primary = _primary_boundaries[contact];
        if (primary.boundary.faces.size() <= spatial_detail::contact_search_tree_minimum_items) continue;
        const Hex8RegionMesh& mesh = _meshes[primary.region];
        _contact_search_boxes.clear();
        _contact_search_boxes.reserve(primary.boundary.faces.size());
        for (std::size_t face = 0; face < primary.boundary.faces.size(); ++face) {
            spatial_detail::ContactSearchBox box;
            box.minimum.fill(std::numeric_limits<double>::infinity());
            box.maximum.fill(-std::numeric_limits<double>::infinity());
            box.item = face;
            for (std::size_t local_node : primary.boundary.faces[face].nodes) {
                const std::size_t node = global_node(primary.region, local_node);
                const CartesianPoint3& reference = mesh.nodes().at(local_node);
                const std::array<double, 3> current = {reference.x + state[dof(Field::displacement_x, node)],
                    reference.y + state[dof(Field::displacement_y, node)],
                    reference.z + state[dof(Field::displacement_z, node)]};
                for (std::size_t component = 0; component < 3; ++component) {
                    box.minimum[component] = std::min(box.minimum[component], current[component]);
                    box.maximum[component] = std::max(box.maximum[component], current[component]);
                }
            }
            _contact_search_boxes.push_back(box);
        }
        if (_contact_search_trees[contact].can_refit(_contact_search_boxes.size()))
            _contact_search_trees[contact].refit(_contact_search_boxes);
        else
            _contact_search_trees[contact].build(_contact_search_boxes);
    }
}

std::size_t SpatialAssembly::mechanical_node_index(std::size_t contact, std::size_t node) const noexcept {
    return _mechanical_contact_offsets[contact] + node;
}

bool SpatialAssembly::mark_touched_thermal_points(std::size_t first, std::size_t last) const {
    std::fill(_touched_thermal_points.begin(), _touched_thermal_points.end(), 0U);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin = std::max(first, ranges.thermal_begin), end = std::min(last, ranges.mechanical_begin);
    if (begin >= end) return false;
    for (std::size_t entry = begin; entry < end; ++entry) _touched_thermal_points[entry - ranges.thermal_begin] = 1U;
    return true;
}

bool SpatialAssembly::mark_touched_mechanical_nodes(std::size_t first, std::size_t last) const {
    std::fill(_touched_mechanical_nodes.begin(), _touched_mechanical_nodes.end(), 0U);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin = std::max(first, ranges.mechanical_begin), end = std::min(last, ranges.boundary_begin);
    if (begin >= end) return false;
    for (std::size_t entry = begin; entry < end; ++entry) {
        const MechanicalPoint& point = _mechanical_points[entry - ranges.mechanical_begin];
        _touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] = 1U;
    }
    return true;
}

void SpatialAssembly::update_thermal_candidates(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional thermal-contact shadow state size mismatch");
    if (!mark_touched_thermal_points(first, last)) return;
    for (std::size_t point = 0; point < _thermal_active_primary.size(); ++point) {
        if (_touched_thermal_points[point] == 0U) continue;
        const auto location =
            offset_location(_thermal_contact_offsets, point, "Three-dimensional thermal-contact point is out of range");
        const std::size_t contact = location.first, primary_count = _primary_contact_faces[contact].size();
        _thermal_minimum_distance[point] = std::numeric_limits<double>::infinity();
        _thermal_active_primary[point] = std::numeric_limits<std::size_t>::max();
        const auto consider = [this, point, &state](std::size_t primary) {
            const ThermalCandidate candidate = thermal_candidate(point, primary);
            const ContactProjectionValue value =
                compute_quad4_to_quad4_heat_projection(candidate.geometry, contact_state(candidate.nodes, state));
            if (!value.projected) return;
            const double distance = std::abs(value.gap);
            const std::size_t selected = _thermal_active_primary[point];
            if (distance < _thermal_minimum_distance[point] ||
                (distance == _thermal_minimum_distance[point] &&
                    (selected == std::numeric_limits<std::size_t>::max() || primary < selected))) {
                _thermal_minimum_distance[point] = distance;
                _thermal_active_primary[point] = primary;
            }
        };
        const std::size_t cached_primary = _thermal_cached_primary[point];
        if (cached_primary < primary_count) consider(cached_primary);
        if (primary_count > spatial_detail::contact_search_tree_minimum_items) {
            const ThermalCandidate representative = thermal_candidate(point, 0);
            const Quad4SurfaceContactLocalValues representative_state = contact_state(representative.nodes, state);
            _contact_search_trees[contact].begin_query(
                thermal_search_point(representative.geometry, representative_state), _contact_search_query);
            std::size_t primary = 0;
            while (_contact_search_trees[contact].next_candidate(
                _contact_search_query, _thermal_minimum_distance[point], primary)) {
                if (primary != cached_primary) consider(primary);
            }
        } else
            for (std::size_t primary = 0; primary < primary_count; ++primary)
                if (primary != cached_primary) consider(primary);
        if (_thermal_active_primary[point] != std::numeric_limits<std::size_t>::max())
            _thermal_cached_primary[point] = _thermal_active_primary[point];
    }
}

void SpatialAssembly::update_mechanical_candidates(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional mechanical-contact shadow state size mismatch");
    if (!mark_touched_mechanical_nodes(first, last)) return;
    for (std::size_t node = 0; node < _mechanical_selected_primary.size(); ++node)
        if (_touched_mechanical_nodes[node] != 0U) {
            _mechanical_minimum_distance[node] = std::numeric_limits<double>::infinity();
            _mechanical_selected_primary[node] = std::numeric_limits<std::size_t>::max();
        }
    for (std::size_t point = 0; point < _mechanical_points.size(); ++point) {
        const MechanicalPoint& metadata = _mechanical_points[point];
        const std::size_t node = mechanical_node_index(metadata.contact, metadata.secondary);
        _mechanical_active_primary[point] = std::numeric_limits<std::size_t>::max();
        if (_touched_mechanical_nodes[node] == 0U) continue;
        const std::size_t primary_count = _primary_contact_faces[metadata.contact].size();
        const auto consider = [this, point, node, &state](std::size_t primary) {
            const MechanicalCandidate candidate = mechanical_candidate(point, primary);
            const ContactProjectionValue value =
                compute_node_to_quad4_contact_projection(candidate.geometry, contact_state(candidate.nodes, state));
            if (!value.projected) return;
            const double distance = std::abs(value.gap);
            if (distance < _mechanical_minimum_distance[node] ||
                (distance == _mechanical_minimum_distance[node] && primary < _mechanical_selected_primary[node])) {
                _mechanical_minimum_distance[node] = distance;
                _mechanical_selected_primary[node] = primary;
            }
        };
        const std::size_t cached_primary = _mechanical_cached_primary[node];
        if (cached_primary < primary_count) consider(cached_primary);
        if (primary_count > spatial_detail::contact_search_tree_minimum_items) {
            const MechanicalCandidate representative = mechanical_candidate(point, 0);
            const Quad4SurfaceContactLocalValues representative_state = contact_state(representative.nodes, state);
            _contact_search_trees[metadata.contact].begin_query(
                mechanical_search_point(representative.geometry, representative_state), _contact_search_query);
            std::size_t primary = 0;
            while (_contact_search_trees[metadata.contact].next_candidate(
                _contact_search_query, _mechanical_minimum_distance[node], primary)) {
                if (primary != cached_primary) consider(primary);
            }
        } else
            for (std::size_t primary = 0; primary < primary_count; ++primary)
                if (primary != cached_primary) consider(primary);
    }
    for (std::size_t point = 0; point < _mechanical_points.size(); ++point) {
        const MechanicalPoint& metadata = _mechanical_points[point];
        const std::size_t node = mechanical_node_index(metadata.contact, metadata.secondary);
        const std::size_t primary = _mechanical_selected_primary[node];
        if (_touched_mechanical_nodes[node] == 0U || primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        const ContactProjectionValue value =
            compute_node_to_quad4_contact_projection(candidate.geometry, contact_state(candidate.nodes, state));
        if (value.projected) _mechanical_active_primary[point] = primary;
    }
    for (std::size_t node = 0; node < _mechanical_selected_primary.size(); ++node)
        if (_touched_mechanical_nodes[node] != 0U &&
            _mechanical_selected_primary[node] != std::numeric_limits<std::size_t>::max())
            _mechanical_cached_primary[node] = _mechanical_selected_primary[node];
}

std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("Three-dimensional contribution range is out of bounds");
    std::vector<std::size_t> result;
    result.reserve(dof_count());
    std::vector<std::size_t> dofs;
    const auto append_face = [this, &result](const std::array<std::size_t, 4>& nodes) {
        for (const Field field :
            {Field::temperature, Field::displacement_x, Field::displacement_y, Field::displacement_z})
            for (const std::size_t node : nodes) result.push_back(dof(field, node));
    };
    for (std::size_t entry = first; entry < last; ++entry) {
        contribution_dofs(entry, dofs);
        result.insert(result.end(), dofs.begin(), dofs.end());
    }
    std::vector<unsigned char> touched_contacts(_definition.contacts.size(), 0U);
    if (mark_touched_thermal_points(first, last)) {
        for (std::size_t point = 0; point < _touched_thermal_points.size(); ++point) {
            if (_touched_thermal_points[point] == 0U) continue;
            const auto location = offset_location(
                _thermal_contact_offsets, point, "Three-dimensional thermal-contact point is out of range");
            touched_contacts[location.first] = 1U;
            append_face(_secondary_contact_faces[location.first][location.second / 4].nodes);
        }
    }
    if (mark_touched_mechanical_nodes(first, last)) {
        for (const MechanicalPoint& point : _mechanical_points) {
            if (_touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] == 0U) continue;
            touched_contacts[point.contact] = 1U;
            append_face(_secondary_contact_faces[point.contact][point.secondary_face].nodes);
        }
    }
    for (std::size_t contact = 0; contact < touched_contacts.size(); ++contact)
        if (touched_contacts[contact] != 0U)
            for (const PrimaryContactFace& face : _primary_contact_faces[contact]) append_face(face.nodes);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void SpatialAssembly::validate_local_state(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("Three-dimensional contribution range is out of bounds");
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional local validation state has the wrong size");
    update_contact_search_trees(state);
    update_thermal_candidates(first, last, state);
    update_mechanical_candidates(first, last, state);
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        if (_definition.contacts[contact].thermal) {
            std::size_t unprojected = 0;
            for (std::size_t point = _thermal_contact_offsets[contact]; point < _thermal_contact_offsets[contact + 1];
                ++point)
                if (_touched_thermal_points[point] != 0U &&
                    _thermal_active_primary[point] == std::numeric_limits<std::size_t>::max())
                    ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Three-dimensional thermal contact '" + _definition.contacts[contact].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        " secondary integration points after searching every primary face");
        }
        if (_definition.contacts[contact].mechanical) {
            std::size_t unprojected = 0;
            for (std::size_t node = _mechanical_contact_offsets[contact];
                node < _mechanical_contact_offsets[contact + 1]; ++node)
                if (_touched_mechanical_nodes[node] != 0U &&
                    _mechanical_selected_primary[node] == std::numeric_limits<std::size_t>::max())
                    ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Three-dimensional mechanical contact '" + _definition.contacts[contact].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        " secondary nodes after searching every primary face");
        }
    }
}

void SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional committed contact state size mismatch");
    update_contact_search_trees(state);
    update_mechanical_candidates(0, contribution_count(), state);
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    std::vector<std::vector<bool>> updated(_definition.contacts.size());
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact)
        updated[contact].resize(_contact_histories[contact].size(), false);
    for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
        const std::size_t primary = _mechanical_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        const CartesianContactPointValue value =
            compute_node_to_quad4_contact_value(_mechanical_properties[candidate.contact], candidate.geometry,
                contact_state(candidate.nodes, state), contact_state(candidate.nodes, _committed_contact_solution),
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        ContactPointHistory trial = _contact_histories[candidate.contact][candidate.secondary];
        trial.sliding = value.sliding;
        trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
        if (updated[candidate.contact][candidate.secondary]) {
            const ContactPointHistory& prior = staged[candidate.contact][candidate.secondary];
            bool equal = prior.sliding == trial.sliding;
            for (std::size_t component = 0; component < 3; ++component) {
                const double scale = std::max({1.0, std::abs(prior.cartesian_elastic_tangential_slip[component]),
                    std::abs(trial.cartesian_elastic_tangential_slip[component])});
                equal = equal && std::abs(prior.cartesian_elastic_tangential_slip[component] -
                                          trial.cartesian_elastic_tangential_slip[component]) <=
                                     64.0 * std::numeric_limits<double>::epsilon() * scale;
            }
            if (!equal)
                throw std::logic_error(
                    "Three-dimensional secondary half-faces disagree on contact-node friction history");
            continue;
        }
        staged[candidate.contact][candidate.secondary] = trial;
        updated[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        if (!_definition.contacts[contact].mechanical) continue;
        if (std::find(updated[contact].begin(), updated[contact].end(), false) != updated[contact].end())
            throw std::domain_error("Cannot commit three-dimensional friction history for an unprojected node");
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
}

void SpatialAssembly::restore_contact_state(
    const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories) {
    if (state.size() != dof_count() || histories.size() != _definition.contacts.size())
        throw std::invalid_argument("Three-dimensional restored contact state layout mismatch");
    for (std::size_t contact = 0; contact < histories.size(); ++contact) {
        if (histories[contact].size() != _contact_histories[contact].size())
            throw std::invalid_argument("Three-dimensional restored contact history layout mismatch");
        for (const ContactPointHistory& history : histories[contact]) {
            if (!std::isfinite(history.normal_multiplier) || history.normal_multiplier < 0.0)
                throw std::invalid_argument("Three-dimensional restored normal contact history is invalid");
            for (double component : history.cartesian_elastic_tangential_slip)
                if (!std::isfinite(component))
                    throw std::invalid_argument("Three-dimensional restored friction history is invalid");
        }
    }
    _committed_contact_solution = state;
    _contact_histories = std::move(histories);
}

std::vector<CartesianContactNodeSummary> SpatialAssembly::summarize_contact_nodes(
    std::size_t contact_value, const std::vector<double>& state) const {
    if (contact_value >= _definition.contacts.size())
        throw std::out_of_range("Three-dimensional contact index is out of range");
    update_contact_search_trees(state);
    update_mechanical_candidates(0, contribution_count(), state);
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const Hex8RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<CartesianContactNodeSummary> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) {
        const CartesianPoint3& point = mesh.nodes().at(node);
        result.push_back({point.x, point.y, point.z, false, std::numeric_limits<std::size_t>::max(),
            std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0, {}, false});
    }
    for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
        const std::size_t primary = _mechanical_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        if (candidate.contact != contact_value) continue;
        const std::size_t node = mechanical_node_index(candidate.contact, candidate.secondary);
        if (candidate.primary != _mechanical_selected_primary[node]) continue;
        const CartesianContactPointValue value =
            compute_node_to_quad4_contact_value(_mechanical_properties[candidate.contact], candidate.geometry,
                contact_state(candidate.nodes, state), contact_state(candidate.nodes, _committed_contact_solution),
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        CartesianContactNodeSummary& summary = result.at(candidate.secondary);
        if (summary.projected && summary.primary_face != candidate.primary)
            throw std::logic_error("Three-dimensional contact node has more than one active primary face");
        summary.projected = true;
        summary.primary_face = candidate.primary;
        summary.gap = std::min(summary.gap, value.gap);
        summary.tributary_area += value.tributary_area;
        summary.contact_force += value.contact_force;
        summary.tangential_force += value.tangential_force;
        summary.elastic_tangential_slip = value.elastic_tangential_slip;
        summary.sliding = value.sliding;
    }
    for (CartesianContactNodeSummary& summary : result)
        if (summary.tributary_area > 0.0) {
            summary.pressure = summary.contact_force / summary.tributary_area;
            summary.tangential_traction = summary.tangential_force / summary.tributary_area;
        }
    return result;
}

std::vector<std::size_t> SpatialAssembly::contact_secondary_source_nodes(std::size_t contact_value) const {
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const Hex8RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<std::size_t> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) result.push_back(mesh.source_node_ids().at(node));
    return result;
}

InterfaceSummary SpatialAssembly::summarize_interface(
    std::size_t contact_value, const std::vector<double>& state) const {
    if (contact_value >= _definition.contacts.size())
        throw std::out_of_range("Three-dimensional contact index is out of range");
    validate_state(state);
    InterfaceSummary result;
    for (std::size_t point = _thermal_contact_offsets[contact_value];
        point < _thermal_contact_offsets[contact_value + 1]; ++point) {
        const std::size_t primary = _thermal_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const ThermalCandidate candidate = thermal_candidate(point, primary);
        const CartesianHeatQuadratureValue value = compute_quad4_to_quad4_gap_heat_value(
            _thermal_properties[contact_value], candidate.geometry, contact_state(candidate.nodes, state));
        if (!value.projected)
            throw std::logic_error("Active three-dimensional thermal-contact candidate is not projected");
        result.minimum_gap = std::min(result.minimum_gap, value.gap);
        result.total_heat_rate += value.weighted_measure * value.heat_flux;
    }
    if (!_definition.contacts[contact_value].thermal) result.minimum_gap = 0.0;
    if (_definition.contacts[contact_value].mechanical) {
        for (const CartesianContactNodeSummary& node : summarize_contact_nodes(contact_value, state)) {
            if (!node.projected) {
                ++result.unprojected_contact_nodes;
                continue;
            }
            ++result.projected_contact_nodes;
            result.minimum_contact_gap = std::min(result.minimum_contact_gap, node.gap);
            result.maximum_contact_pressure = std::max(result.maximum_contact_pressure, node.pressure);
            result.total_contact_force += node.contact_force;
            result.total_tangential_force += node.tangential_force;
            if (node.pressure > 0.0) ++result.active_contact_nodes;
        }
    } else {
        result.minimum_contact_gap = 0.0;
    }
    return result;
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
