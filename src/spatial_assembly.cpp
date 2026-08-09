#include "spatial_assembly.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

Line2InterfaceSideCoordinates
edge_coordinates(const RegionMesh& mesh, const Line2BoundaryElement& edge) {
    return {{mesh.nodes().at(edge.nodes[0]), mesh.nodes().at(edge.nodes[1])}};
}

double planar_quad_area(const RegionMesh& mesh, const Quad4Element& element) {
    double twice_area = 0.0;
    for (std::size_t node = 0; node < element.nodes.size(); ++node) {
        const RzPoint& current = mesh.nodes().at(element.nodes[node]);
        const RzPoint& next =
            mesh.nodes().at(element.nodes[(node + 1U) % element.nodes.size()]);
        twice_area += current.r * next.z - next.r * current.z;
    }
    return 0.5 * std::abs(twice_area);
}

double minimum_boundary_normal_length(const SpatialLayout& layout,
                                      std::size_t region,
                                      const RegionBoundary& boundary) {
    const RegionMesh& mesh = layout.region_mesh(region);
    double result = std::numeric_limits<double>::infinity();
    for (const Line2BoundaryElement& edge : boundary.elements) {
        const auto parent = layout.edge_parent(region, edge);
        const Line2InterfaceSideCoordinates coordinates =
            edge_coordinates(mesh, edge);
        const double edge_length =
            std::hypot(coordinates[1].r - coordinates[0].r,
                       coordinates[1].z - coordinates[0].z);
        const double area =
            planar_quad_area(mesh, mesh.elements().at(parent.first));
        const double normal_length = area / edge_length;
        if (!std::isfinite(normal_length) || !(normal_length > 0.0))
            throw std::invalid_argument(
                "Contact boundary has a nonpositive characteristic element "
                "length");
        result = std::min(result, normal_length);
    }
    return result;
}

RzPoint element_centroid(const RegionMesh& mesh, const Quad4Element& element) {
    RzPoint centroid{0.0, 0.0};
    for (const std::size_t node : element.nodes) {
        centroid.r += mesh.nodes().at(node).r;
        centroid.z += mesh.nodes().at(node).z;
    }
    centroid.r /= static_cast<double>(element.nodes.size());
    centroid.z /= static_cast<double>(element.nodes.size());
    return centroid;
}

// Signed distance of `point` from the primary line along the base normal
// (tangent_z, -tangent_r)/length. This is the same normal convention as the
// raw gap in reference_normal_orientation in interface.cpp, so the sign of
// this value for a point that rides on the line is exactly the raw-gap sign
// the point would have if it were moved to that side of the line.
double
primary_line_side(const RzPoint& point,
                  const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double tangent_r =
        primary_coordinates[1].r - primary_coordinates[0].r;
    const double tangent_z =
        primary_coordinates[1].z - primary_coordinates[0].z;
    const double length = std::hypot(tangent_r, tangent_z);
    return ((point.r - primary_coordinates[0].r) * tangent_z -
            (point.z - primary_coordinates[0].z) * tangent_r) /
           length;
}

// Zero-gap orientation hint from the material-side topology: the signed side
// of the secondary parent-element centroid relative to the primary line. If
// both parent centroids fall on the same strict side of the line, the two
// material bodies overlap next to the interface and no meaningful hint
// exists; returning zero keeps the explicit construction error for any
// secondary point that actually rides on the segment. If either centroid
// lies on the line, the corresponding element is degenerate and the same
// explicit error remains in force.
double zero_gap_orientation_hint(
    const RzPoint& secondary_centroid, const RzPoint& primary_centroid,
    const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double secondary_side =
        primary_line_side(secondary_centroid, primary_coordinates);
    const double primary_side =
        primary_line_side(primary_centroid, primary_coordinates);
    if (secondary_side == 0.0 || primary_side == 0.0)
        return 0.0;
    const bool same_side = (secondary_side > 0.0 && primary_side > 0.0) ||
                           (secondary_side < 0.0 && primary_side < 0.0);
    return same_side ? 0.0 : secondary_side;
}

void validate_definitions(const SpatialDefinition& definition) {
    if (definition.regions.empty())
        throw std::invalid_argument(
            "SpatialAssembly requires at least one region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        const RegionDefinition& value = definition.regions[region];
        if (value.name.empty() || (value.block.empty() && value.block_id < 0) ||
            (!value.block.empty() && value.block_id >= 0))
            throw std::invalid_argument(
                "SpatialAssembly regions require a name and exactly one block "
                "selector");
        if (!std::isfinite(value.volumetric_heat_source) ||
            value.volumetric_heat_source < 0.0)
            throw std::invalid_argument("SpatialAssembly region heat sources "
                                        "must be finite and nonnegative");
        if (!std::isfinite(value.initial_temperature) ||
            !(value.initial_temperature > 0.0))
            throw std::invalid_argument("SpatialAssembly region temperatures "
                                        "must be finite and positive");
        for (std::size_t previous = 0; previous < region; ++previous) {
            if (definition.regions[previous].name == value.name)
                throw std::invalid_argument("Duplicate region name: " +
                                            value.name);
            if ((!value.block.empty() &&
                 definition.regions[previous].block == value.block) ||
                (value.block_id >= 0 &&
                 definition.regions[previous].block_id == value.block_id))
                throw std::invalid_argument(
                    "Duplicate region block: " +
                    (value.block.empty() ? std::to_string(value.block_id)
                                         : value.block));
        }
    }
    for (std::size_t table = 0; table < definition.time_tables.size();
         ++table) {
        for (std::size_t previous = 0; previous < table; ++previous) {
            if (definition.time_tables[previous].name() ==
                definition.time_tables[table].name())
                throw std::invalid_argument(
                    "Duplicate time-table name: " +
                    definition.time_tables[table].name());
        }
    }
    for (std::size_t contact = 0; contact < definition.contacts.size();
         ++contact) {
        const ContactDefinition& value = definition.contacts[contact];
        if (value.name.empty() || value.primary.empty() ||
            value.secondary.empty())
            throw std::invalid_argument(
                "SpatialAssembly contact names and boundaries must not be empty");
        if (!value.thermal && !value.mechanical)
            throw std::invalid_argument(
                "Contact must enable thermal or mechanical coupling: " +
                value.name);
        if (value.primary == value.secondary)
            throw std::invalid_argument(
                "Contact primary and secondary must differ: " + value.name);
        if (value.thermal &&
            (!(value.gap_conductivity > 0.0) || !(value.minimum_gap > 0.0)))
            throw std::invalid_argument(
                "Thermal contact parameters must be positive: " + value.name);
        if (value.mechanical && !value.automatic_penalty &&
            (!std::isfinite(value.penalty) || !(value.penalty > 0.0)))
            throw std::invalid_argument(
                "Mechanical contact penalty must be positive: " + value.name);
        if (value.mechanical && (!std::isfinite(value.penalty_factor) ||
                                 !(value.penalty_factor > 0.0)))
            throw std::invalid_argument(
                "Mechanical contact penalty factor must be finite and "
                "positive: " +
                value.name);
        if (value.mechanical &&
            value.mechanical_formulation ==
                MechanicalContactFormulation::augmented_lagrangian &&
            (!std::isfinite(value.penetration_tolerance) ||
             !(value.penetration_tolerance > 0.0) ||
             value.maximum_augmented_iterations == 0))
            throw std::invalid_argument(
                "Augmented contact requires a positive penetration tolerance "
                "and iteration limit: " +
                value.name);
        if (!std::isfinite(value.friction_coefficient) ||
            value.friction_coefficient < 0.0)
            throw std::invalid_argument(
                "Mechanical contact friction coefficient must be finite and "
                "nonnegative: " +
                value.name);
        for (std::size_t previous = 0; previous < contact; ++previous) {
            if (definition.contacts[previous].name == value.name)
                throw std::invalid_argument("Duplicate contact name: " +
                                            value.name);
            if (value.mechanical && definition.contacts[previous].mechanical &&
                definition.contacts[previous].secondary == value.secondary)
                throw std::invalid_argument(
                    "Mechanical secondary boundary is reused: " +
                    value.secondary);
        }
    }
}

RegionBoundary ordered_connected_boundary(const RegionMesh& mesh,
                                          RegionBoundary boundary,
                                          const std::string& name) {
    if (boundary.elements.empty())
        throw std::invalid_argument("Contact side set is empty: " + name);

    std::vector<std::size_t> degree(mesh.nodes().size(), 0);
    for (const Line2BoundaryElement& edge : boundary.elements) {
        if (edge.nodes[0] == edge.nodes[1])
            throw std::invalid_argument("Contact side set has a zero edge: " +
                                        name);
        for (std::size_t node : edge.nodes) {
            if (++degree.at(node) > 2)
                throw std::invalid_argument(
                    "Contact side set must not branch: " + name);
        }
    }
    std::vector<std::size_t> endpoints;
    for (std::size_t node = 0; node < degree.size(); ++node) {
        if (degree[node] == 1)
            endpoints.push_back(node);
    }
    if (endpoints.size() != 2)
        throw std::invalid_argument(
            "Contact side set must be one open connected chain: " + name);

    const auto coordinate_less = [&](std::size_t lhs, std::size_t rhs) {
        const RzPoint& left = mesh.nodes().at(lhs);
        const RzPoint& right = mesh.nodes().at(rhs);
        if (left.r != right.r)
            return left.r < right.r;
        if (left.z != right.z)
            return left.z < right.z;
        return lhs < rhs;
    };
    std::size_t current = coordinate_less(endpoints[0], endpoints[1])
                              ? endpoints[0]
                              : endpoints[1];
    std::vector<bool> used(boundary.elements.size(), false);
    std::vector<Line2BoundaryElement> ordered;
    std::vector<std::size_t> nodes = {current};
    ordered.reserve(boundary.elements.size());
    nodes.reserve(boundary.elements.size() + 1);
    for (std::size_t position = 0; position < boundary.elements.size();
         ++position) {
        std::size_t selected = boundary.elements.size();
        Line2BoundaryElement oriented{};
        for (std::size_t edge = 0; edge < boundary.elements.size(); ++edge) {
            if (used[edge])
                continue;
            const Line2BoundaryElement& candidate = boundary.elements[edge];
            if (candidate.nodes[0] == current) {
                selected = edge;
                oriented = candidate;
                break;
            }
            if (candidate.nodes[1] == current) {
                selected = edge;
                oriented = {{{candidate.nodes[1], candidate.nodes[0]}}};
                break;
            }
        }
        if (selected == boundary.elements.size())
            throw std::invalid_argument(
                "Contact side set must be one connected chain: " + name);
        used[selected] = true;
        ordered.push_back(oriented);
        current = oriented.nodes[1];
        nodes.push_back(current);
    }
    boundary.nodes = std::move(nodes);
    boundary.elements = std::move(ordered);
    return boundary;
}

double projection_fraction(const RzPoint& point,
                           const Line2InterfaceSideCoordinates& segment) {
    const double dr = segment[1].r - segment[0].r;
    const double dz = segment[1].z - segment[0].z;
    return ((point.r - segment[0].r) * dr + (point.z - segment[0].z) * dz) /
           (dr * dr + dz * dz);
}

bool projection_interval(const Line2InterfaceSideCoordinates& secondary,
                         const Line2InterfaceSideCoordinates& primary,
                         double& lower, double& upper) {
    const RzPoint midpoint = {
        0.5 * (secondary[0].r + secondary[1].r),
        0.5 * (secondary[0].z + secondary[1].z),
    };
    const RzPoint half = {
        0.5 * (secondary[1].r - secondary[0].r),
        0.5 * (secondary[1].z - secondary[0].z),
    };
    const double center = projection_fraction(midpoint, primary);
    const double slope =
        projection_fraction({midpoint.r + half.r, midpoint.z + half.z},
                            primary) -
        center;
    lower = -1.0;
    upper = 1.0;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(slope) <= tolerance) {
        return center >= -tolerance && center <= 1.0 + tolerance;
    }
    double first = (0.0 - center) / slope;
    double second = (1.0 - center) / slope;
    if (first > second)
        std::swap(first, second);
    lower = std::max(lower, first);
    upper = std::min(upper, second);
    lower = std::max(-1.0, std::min(1.0, lower));
    upper = std::max(-1.0, std::min(1.0, upper));
    return upper - lower > tolerance;
}

} // namespace

std::vector<std::int64_t>
SpatialAssembly::resolve_block_ids(const SpatialDefinition& definition,
                                 const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<std::int64_t> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        const std::int64_t id =
            region.block_id >= 0 ? region.block_id
                                 : source_mesh.element_block(region.block).id;
        const bool known = std::any_of(
            source_mesh.element_blocks().begin(),
            source_mesh.element_blocks().end(),
            [id](const ElementBlockInfo& block) { return block.id == id; });
        if (!known)
            throw std::invalid_argument("Unknown element block ID: " +
                                        std::to_string(id));
        if (std::find(result.begin(), result.end(), id) != result.end())
            throw std::invalid_argument("Duplicate resolved region block ID: " +
                                        std::to_string(id));
        result.push_back(id);
    }
    return result;
}

std::vector<RegionMesh>
SpatialAssembly::build_meshes(const SpatialDefinition& definition,
                            const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<RegionMesh> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        if (region.block_id >= 0)
            result.push_back(RegionMesh::from_unstructured_block(
                source_mesh, region.block_id));
        else
            result.push_back(
                RegionMesh::from_unstructured_block(source_mesh, region.block));
    }
    return result;
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh)
    : SpatialAssembly(definition, source_mesh,
                    resolve_block_ids(definition, source_mesh),
                    build_meshes(definition, source_mesh)) {}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh,
                             std::vector<std::int64_t> block_ids,
                             std::vector<RegionMesh> meshes)
    : _layout(std::move(definition), std::move(block_ids), std::move(meshes)) {
    validate_definitions(_layout._definition);

    _layout.build_volume_geometries();
    build_contacts(source_mesh);
    _boundary.build(source_mesh, _layout);
    _contact._contact_histories.resize(contact_count());
    _contact._projected_mechanical_nodes.resize(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        _contact._contact_histories[contact_value].resize(
            _contact._secondary_boundaries[contact_value].boundary.nodes.size());
        _contact._projected_mechanical_nodes[contact_value].resize(
            _contact._secondary_boundaries[contact_value].boundary.nodes.size(), false);
    }
    _contact._committed_contact_solution = initial_state();
    update_mechanical_candidates(_contact._committed_contact_solution);
}

std::vector<double> SpatialAssembly::initial_state() const {
    std::vector<double> result(dof_count(), 0.0);
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const std::size_t offset = _layout._node_offsets[region_value];
        for (std::size_t node = 0; node < _layout._meshes[region_value].nodes().size();
             ++node) {
            result[_layout._dof_map.temperature(offset + node)] =
                _layout._definition.regions[region_value].initial_temperature;
        }
    }
    for (const DirichletCondition& condition :
         _boundary.dirichlet_conditions())
        result[condition.dof] = condition.value;
    return result;
}

std::size_t SpatialAssembly::dof_count() const noexcept {
    return _layout._dof_map.dof_count();
}

std::size_t SpatialAssembly::contribution_count() const noexcept {
    return contribution_ranges().end;
}

const std::vector<DirichletCondition>&
SpatialAssembly::dirichlet_conditions() const noexcept {
    return _boundary.dirichlet_conditions();
}

std::vector<std::size_t>
SpatialAssembly::required_state_dofs(std::size_t contribution_begin,
                                     std::size_t contribution_end) const {
    if (contribution_begin > contribution_end ||
        contribution_end > contribution_count())
        throw std::out_of_range(
            "SpatialAssembly contribution range is out of bounds");
    std::vector<std::size_t> result;
    result.reserve((contribution_end - contribution_begin) * local_dof_count);
    for (std::size_t contribution = contribution_begin;
         contribution < contribution_end; ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        for (const std::size_t dof : dofs) {
            if (dof >= dof_count())
                throw std::out_of_range(
                    "SpatialAssembly contribution has an invalid DOF");
            result.push_back(dof);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t local_begin =
        std::max(contribution_begin, ranges.mechanical_begin);
    const std::size_t local_end =
        std::min(contribution_end, ranges.pressure_begin);
    if (local_begin >= local_end)
        return result;

    std::vector<std::vector<bool>> touched(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value)
        touched[contact_value].resize(_contact._contact_histories[contact_value].size(),
                                      false);
    for (std::size_t full = local_begin; full < local_end; ++full) {
        const std::size_t contribution = full - ranges.mechanical_begin;
        const MechanicalContribution& candidate =
            _contact._mechanical_contributions[contribution];
        touched[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contribution = 0;
         contribution < _contact._mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _contact._mechanical_contributions[contribution];
        if (!touched[candidate.contact][candidate.secondary])
            continue;
        const LocalDofs dofs =
            contribution_dofs(ranges.mechanical_begin + contribution);
        result.insert(result.end(), dofs.begin(), dofs.end());
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::pair<std::size_t, std::size_t>
SpatialAssembly::element_location(std::size_t contribution_index) const {
    if (contribution_index >= _layout._element_offsets.back())
        throw std::out_of_range(
            "SpatialAssembly volume contribution is out of range");
    const auto upper = std::upper_bound(
        _layout._element_offsets.begin(), _layout._element_offsets.end(), contribution_index);
    const std::size_t region_value =
        static_cast<std::size_t>(upper - _layout._element_offsets.begin() - 1);
    return {region_value, contribution_index - _layout._element_offsets[region_value]};
}

LocalDofs
SpatialAssembly::contribution_dofs(std::size_t contribution_index) const {
    const ContributionLocation contribution =
        locate_contribution(contribution_index);
    switch (contribution.type) {
    case SpatialContributionType::volume: {
        const auto location = element_location(contribution.local_index);
        const Quad4Element& element =
            _layout._meshes[location.first].elements().at(location.second);
        std::array<std::size_t, 4> nodes{};
        for (std::size_t node = 0; node < nodes.size(); ++node)
            nodes[node] =
                _layout.global_node(location.first, element.nodes[node]);
        return _layout._dof_map.local_dofs(nodes);
    }
    case SpatialContributionType::thermal_contact:
        return _layout._dof_map.local_dofs(
            _contact._thermal_contributions.at(contribution.local_index).nodes);
    case SpatialContributionType::mechanical_contact:
        return _layout._dof_map.local_dofs(
            _contact._mechanical_contributions.at(contribution.local_index).nodes);
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection:
        return _boundary.contribution_dofs(
            contribution.type, contribution.local_index, _layout);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}

LocalValues SpatialAssembly::contribution_state(
    std::size_t contribution_index,
    const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly global state has the wrong size");
    return contribution_state(contribution_index,
                              GlobalStateView(global_state));
}

LocalValues SpatialAssembly::contribution_state(
    std::size_t contribution_index,
    const GlobalStateView& global_state) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly shadow state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(contribution_index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local)
        result[local] = global_state.value(dofs[local]);
    return result;
}

LocalResidual
SpatialAssembly::contribution_residual(std::size_t contribution_index,
                                       const LocalValues& state) const {
    const ContributionLocation location =
        locate_contribution(contribution_index);
    switch (location.type) {
    case SpatialContributionType::volume:
        throw std::logic_error(
            "SpatialAssembly does not own volume residual physics");
    case SpatialContributionType::thermal_contact: {
        const ThermalContribution& contribution =
            _contact._thermal_contributions[location.local_index];
        return _contact._thermal_kernels[contribution.contact].residual(
            contribution.geometry, state);
    }
    case SpatialContributionType::mechanical_contact: {
        const MechanicalContribution& contribution =
            _contact._mechanical_contributions.at(location.local_index);
        if (!contribution.active)
            return {};
        return _contact._mechanical_kernels[contribution.contact].residual(
            contribution.geometry, state,
            contribution_state(contribution_index,
                               _contact._committed_contact_solution),
            _contact._contact_histories[contribution.contact][contribution.secondary]);
    }
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection:
        return _boundary.contribution_residual(
            location.type, location.local_index, state);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}

LocalSystem
SpatialAssembly::linearize_contribution(std::size_t contribution_index,
                                        const LocalValues& state) const {
    const ContributionLocation location =
        locate_contribution(contribution_index);
    switch (location.type) {
    case SpatialContributionType::volume:
        throw std::logic_error(
            "SpatialAssembly does not own volume Jacobian physics");
    case SpatialContributionType::thermal_contact: {
        const ThermalContribution& contribution =
            _contact._thermal_contributions[location.local_index];
        return _contact._thermal_kernels[contribution.contact].linearize(
            contribution.geometry, state);
    }
    case SpatialContributionType::mechanical_contact: {
        const MechanicalContribution& contribution =
            _contact._mechanical_contributions.at(location.local_index);
        if (!contribution.active)
            return {};
        return _contact._mechanical_kernels[contribution.contact].linearize(
            contribution.geometry, state,
            contribution_state(contribution_index,
                               _contact._committed_contact_solution),
            _contact._contact_histories[contribution.contact][contribution.secondary]);
    }
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection:
        return _boundary.linearize_contribution(
            location.type, location.local_index, state);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}

void SpatialAssembly::build_contacts(const UnstructuredQuad4Mesh& source_mesh) {
    _contact._primary_boundaries.reserve(contact_count());
    _contact._secondary_boundaries.reserve(contact_count());
    _contact._thermal_kernels.reserve(contact_count());
    _contact._mechanical_kernels.reserve(contact_count());

    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        ContactDefinition& contact_definition =
            _layout._definition.contacts[contact_value];
        ResolvedBoundary primary =
            _layout.resolve_boundary(source_mesh, contact_definition.primary);
        ResolvedBoundary secondary =
            _layout.resolve_boundary(source_mesh,
                                     contact_definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Self-contact is not supported: " +
                                        contact_definition.name);
        const RegionMesh& primary_mesh = _layout._meshes[primary.region];
        const RegionMesh& secondary_mesh = _layout._meshes[secondary.region];
        primary.boundary = ordered_connected_boundary(
            primary_mesh, std::move(primary.boundary),
            contact_definition.primary);
        secondary.boundary = ordered_connected_boundary(
            secondary_mesh, std::move(secondary.boundary),
            contact_definition.secondary);

        // Parent-element centroids give the material side of each boundary
        // edge; they form the zero-gap orientation hint passed to every
        // thermal and mechanical contact geometry below.
        std::vector<RzPoint> primary_parent_centroids;
        primary_parent_centroids.reserve(primary.boundary.elements.size());
        for (const Line2BoundaryElement& edge : primary.boundary.elements)
            primary_parent_centroids.push_back(element_centroid(
                primary_mesh, primary_mesh.elements().at(
                                  _layout.edge_parent(primary.region, edge)
                                      .first)));
        std::vector<RzPoint> secondary_parent_centroids;
        secondary_parent_centroids.reserve(secondary.boundary.elements.size());
        for (const Line2BoundaryElement& edge : secondary.boundary.elements)
            secondary_parent_centroids.push_back(element_centroid(
                secondary_mesh, secondary_mesh.elements().at(
                                    _layout.edge_parent(secondary.region, edge)
                                        .first)));

        if (contact_definition.mechanical &&
            contact_definition.automatic_penalty) {
            const double primary_length =
                minimum_boundary_normal_length(_layout, primary.region,
                                               primary.boundary);
            const double secondary_length = minimum_boundary_normal_length(
                _layout, secondary.region, secondary.boundary);
            const double primary_modulus =
                _layout._definition.regions[primary.region].material.young_modulus;
            const double secondary_modulus =
                _layout._definition.regions[secondary.region].material.young_modulus;
            contact_definition.penalty = contact_definition.penalty_factor /
                                         (primary_length / primary_modulus +
                                          secondary_length / secondary_modulus);
            if (!std::isfinite(contact_definition.penalty) ||
                !(contact_definition.penalty > 0.0))
                throw std::overflow_error(
                    "Automatic contact penalty is not finite and positive: " +
                    contact_definition.name);
        }

        _contact._thermal_kernels.emplace_back(GapHeatProperties{
            contact_definition.thermal ? contact_definition.gap_conductivity
                                       : 1.0,
            contact_definition.thermal ? contact_definition.minimum_gap : 1.0});
        _contact._mechanical_kernels.emplace_back(NormalContactProperties{
            contact_definition.mechanical ? contact_definition.penalty : 1.0,
            contact_definition.mechanical
                ? contact_definition.friction_coefficient
                : 0.0,
            contact_definition.mechanical &&
                contact_definition.mechanical_formulation ==
                    MechanicalContactFormulation::augmented_lagrangian});

        if (contact_definition.thermal) {
            for (std::size_t edge_index = 0;
                 edge_index < secondary.boundary.elements.size();
                 ++edge_index) {
                const Line2BoundaryElement& secondary_edge =
                    secondary.boundary.elements[edge_index];
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                struct ThermalProjection final {
                    std::size_t primary_edge;
                    double lower;
                    double upper;
                };
                std::vector<ThermalProjection> projections;
                for (std::size_t primary_edge = 0;
                     primary_edge < primary.boundary.elements.size();
                     ++primary_edge) {
                    double lower = 0.0;
                    double upper = 0.0;
                    if (projection_interval(
                            secondary_coordinates,
                            edge_coordinates(
                                primary_mesh,
                                primary.boundary.elements[primary_edge]),
                            lower, upper))
                        projections.push_back({primary_edge, lower, upper});
                }
                std::sort(projections.begin(), projections.end(),
                          [](const ThermalProjection& lhs,
                             const ThermalProjection& rhs) {
                              return lhs.lower < rhs.lower;
                          });
                double covered = -1.0;
                constexpr double coverage_tolerance = 1.0e-10;
                for (const ThermalProjection& projection : projections) {
                    if (projection.lower > covered + coverage_tolerance ||
                        projection.lower < covered - coverage_tolerance)
                        throw std::invalid_argument(
                            "Secondary thermal edge projection has a gap or "
                            "overlap on the primary surface: " +
                            contact_definition.name);
                    covered = projection.upper;
                }
                if (covered < 1.0 - coverage_tolerance)
                    throw std::invalid_argument(
                        "Secondary thermal edge is outside the primary "
                        "surface projection: " +
                        contact_definition.name);
                for (const ThermalProjection& projection : projections) {
                    const Line2BoundaryElement& primary_edge =
                        primary.boundary.elements[projection.primary_edge];
                    const Line2InterfaceSideCoordinates primary_coordinates =
                        edge_coordinates(primary_mesh, primary_edge);
                    _contact._thermal_contributions.push_back(
                        {contact_value,
                         {_layout.global_node(secondary.region,
                                              secondary_edge.nodes[0]),
                          _layout.global_node(secondary.region,
                                              secondary_edge.nodes[1]),
                          _layout.global_node(primary.region,
                                              primary_edge.nodes[0]),
                          _layout.global_node(primary.region,
                                              primary_edge.nodes[1])},
                         make_line2_rz_heat_geometry(
                             secondary_coordinates, primary_coordinates,
                             projection.lower, projection.upper,
                             zero_gap_orientation_hint(
                                 secondary_parent_centroids[edge_index],
                                 primary_parent_centroids[projection
                                                              .primary_edge],
                                 primary_coordinates))});
                }
            }
        }

        if (contact_definition.mechanical) {
            for (std::size_t edge_index = 0;
                 edge_index < secondary.boundary.elements.size();
                 ++edge_index) {
                const Line2BoundaryElement& secondary_edge =
                    secondary.boundary.elements[edge_index];
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                for (std::size_t secondary_node = 0;
                     secondary_node < line2_interface_side_node_count;
                     ++secondary_node) {
                    for (std::size_t candidate = 0;
                         candidate < primary.boundary.elements.size();
                         ++candidate) {
                        const Line2BoundaryElement& primary_edge =
                            primary.boundary.elements[candidate];
                        const Line2InterfaceSideCoordinates
                            primary_coordinates =
                                edge_coordinates(primary_mesh, primary_edge);
                        _contact._mechanical_contributions.push_back(
                            {contact_value,
                             {_layout.global_node(secondary.region,
                                                  secondary_edge.nodes[0]),
                              _layout.global_node(secondary.region,
                                                  secondary_edge.nodes[1]),
                              _layout.global_node(primary.region,
                                                  primary_edge.nodes[0]),
                              _layout.global_node(primary.region,
                                                  primary_edge.nodes[1])},
                             make_node_to_line_rz_contact_geometry(
                                 secondary_coordinates, primary_coordinates,
                                 secondary_node, candidate == 0,
                                 candidate + 1 ==
                                     primary.boundary.elements.size(),
                                 zero_gap_orientation_hint(
                                     secondary_parent_centroids[edge_index],
                                     primary_parent_centroids[candidate],
                                     primary_coordinates)),
                             edge_index + secondary_node,
                             candidate});
                    }
                }
            }
        }
        _contact._primary_boundaries.push_back(std::move(primary));
        _contact._secondary_boundaries.push_back(std::move(secondary));
    }
}

} // namespace fuelsim
