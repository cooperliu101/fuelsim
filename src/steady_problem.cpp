#include "fuelsim/steady_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss = 0.577350269189625764509148780501957456;

std::size_t checked_total_nodes(const std::vector<RegionMesh>& meshes) {
    std::size_t total = 0;
    for (const RegionMesh& mesh : meshes) {
        if (mesh.nodes().size() >
            std::numeric_limits<std::size_t>::max() - total)
            throw std::length_error("SteadyProblem node count overflows");
        total += mesh.nodes().size();
    }
    return total;
}

Quad4Coordinates element_coordinates(const RegionMesh& mesh,
                                     const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}

Line2InterfaceSideCoordinates
edge_coordinates(const RegionMesh& mesh, const Line2BoundaryElement& edge) {
    return {{mesh.nodes().at(edge.nodes[0]), mesh.nodes().at(edge.nodes[1])}};
}

void validate_definitions(const SteadyProblemDefinition& definition) {
    if (definition.regions.empty())
        throw std::invalid_argument(
            "SteadyProblem requires at least one region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        const RegionDefinition& value = definition.regions[region];
        if (value.name.empty() || value.block.empty())
            throw std::invalid_argument(
                "SteadyProblem region names and blocks must not be empty");
        if (!std::isfinite(value.volumetric_heat_source) ||
            value.volumetric_heat_source < 0.0)
            throw std::invalid_argument("SteadyProblem region heat sources "
                                        "must be finite and nonnegative");
        if (!std::isfinite(value.initial_temperature) ||
            !(value.initial_temperature > 0.0))
            throw std::invalid_argument("SteadyProblem region temperatures "
                                        "must be finite and positive");
        for (std::size_t previous = 0; previous < region; ++previous) {
            if (definition.regions[previous].name == value.name)
                throw std::invalid_argument("Duplicate region name: " +
                                            value.name);
            if (definition.regions[previous].block == value.block)
                throw std::invalid_argument("Duplicate region block: " +
                                            value.block);
        }
    }
    for (std::size_t contact = 0; contact < definition.contacts.size();
         ++contact) {
        const ContactDefinition& value = definition.contacts[contact];
        if (value.name.empty() || value.primary.empty() ||
            value.secondary.empty())
            throw std::invalid_argument(
                "SteadyProblem contact names and boundaries must not be empty");
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
        if (value.mechanical && !(value.penalty > 0.0))
            throw std::invalid_argument(
                "Mechanical contact penalty must be positive: " + value.name);
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

void validate_dirichlet_conditions(
    std::vector<DirichletCondition>& conditions) {
    std::sort(conditions.begin(), conditions.end(),
              [](const DirichletCondition& lhs, const DirichletCondition& rhs) {
                  return lhs.dof < rhs.dof;
              });
    for (std::size_t index = 1; index < conditions.size(); ++index) {
        const DirichletCondition& previous = conditions[index - 1];
        const DirichletCondition& current = conditions[index];
        if (previous.dof != current.dof)
            continue;
        if (previous.value != current.value)
            throw std::invalid_argument(
                "SteadyProblem has conflicting Dirichlet conditions");
        throw std::invalid_argument(
            "SteadyProblem has duplicate Dirichlet conditions");
    }
}

std::size_t
find_containing_segment(double z, const RegionMesh& mesh,
                        const std::vector<Line2BoundaryElement>& edges) {
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const double lower = mesh.nodes().at(edges[edge].nodes[0]).z;
        const double upper = mesh.nodes().at(edges[edge].nodes[1]).z;
        const bool includes_upper = edge + 1 == edges.size();
        if (z >= lower && (z < upper || (includes_upper && z <= upper)))
            return edge;
    }
    throw std::invalid_argument("Contact point is outside the primary surface");
}

void validate_connected_radial_boundary(const RegionMesh& mesh,
                                        const RegionBoundary& boundary,
                                        const std::string& name) {
    if (boundary.kind != RegionBoundaryKind::radial_inner &&
        boundary.kind != RegionBoundaryKind::radial_outer)
        throw std::invalid_argument("Contact requires radial RZ side sets: " +
                                    name);
    if (boundary.elements.empty())
        throw std::invalid_argument("Contact side set is empty: " + name);
    for (std::size_t edge = 1; edge < boundary.elements.size(); ++edge) {
        const double previous_upper =
            mesh.nodes().at(boundary.elements[edge - 1].nodes[1]).z;
        const double current_lower =
            mesh.nodes().at(boundary.elements[edge].nodes[0]).z;
        const double scale =
            std::max({1.0, std::abs(previous_upper), std::abs(current_lower)});
        if (std::abs(previous_upper - current_lower) > 1.0e-12 * scale)
            throw std::invalid_argument(
                "Contact side set must be one connected axial chain: " + name);
    }
}

} // namespace

std::vector<std::int64_t>
SteadyProblem::resolve_block_ids(const SteadyProblemDefinition& definition,
                                 const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<std::int64_t> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions)
        result.push_back(source_mesh.element_block(region.block).id);
    return result;
}

std::vector<RegionMesh>
SteadyProblem::build_meshes(const SteadyProblemDefinition& definition,
                            const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<RegionMesh> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions)
        result.push_back(
            RegionMesh::from_unstructured_block(source_mesh, region.block));
    return result;
}

SteadyProblem::SteadyProblem(SteadyProblemDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh)
    : SteadyProblem(definition, source_mesh,
                    resolve_block_ids(definition, source_mesh),
                    build_meshes(definition, source_mesh)) {}

SteadyProblem::SteadyProblem(SteadyProblemDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh,
                             std::vector<std::int64_t> block_ids,
                             std::vector<RegionMesh> meshes)
    : _definition(std::move(definition)), _block_ids(std::move(block_ids)),
      _meshes(std::move(meshes)), _dof_map(checked_total_nodes(_meshes)),
      _load_factor(1.0) {
    validate_definitions(_definition);

    _node_offsets.reserve(_meshes.size() + 1);
    _element_offsets.reserve(_meshes.size() + 1);
    _node_offsets.push_back(0);
    _element_offsets.push_back(0);
    for (const RegionMesh& mesh : _meshes) {
        _node_offsets.push_back(_node_offsets.back() + mesh.nodes().size());
        _element_offsets.push_back(_element_offsets.back() +
                                   mesh.elements().size());
    }

    _region_kernels.reserve(_definition.regions.size());
    for (const RegionDefinition& region : _definition.regions) {
        _region_kernels.emplace_back(
            IsotropicThermoelasticMaterial(region.material),
            region.volumetric_heat_source);
    }
    build_volume_geometries();
    build_contacts(source_mesh);
    build_boundary_conditions(source_mesh);
}

const SteadyProblemDefinition& SteadyProblem::definition() const noexcept {
    return _definition;
}

const DofMap& SteadyProblem::dof_map() const noexcept {
    return _dof_map;
}

std::size_t SteadyProblem::region_count() const noexcept {
    return _definition.regions.size();
}

std::size_t SteadyProblem::region_index(const std::string& name) const {
    for (std::size_t region = 0; region < region_count(); ++region) {
        if (_definition.regions[region].name == name)
            return region;
    }
    throw std::invalid_argument("Unknown region: " + name);
}

const RegionDefinition& SteadyProblem::region(std::size_t index) const {
    return _definition.regions.at(index);
}

const RegionMesh& SteadyProblem::region_mesh(std::size_t index) const {
    return _meshes.at(index);
}

const Quad4RzThermoelasticKernel&
SteadyProblem::region_kernel(std::size_t index) const {
    return _region_kernels.at(index);
}

std::size_t SteadyProblem::region_node_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("SteadyProblem region index is out of range");
    return _node_offsets[index];
}

std::size_t SteadyProblem::region_element_count(std::size_t index) const {
    return _meshes.at(index).elements().size();
}

std::size_t SteadyProblem::region_element_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("SteadyProblem region index is out of range");
    return _element_offsets[index];
}

std::size_t SteadyProblem::volume_contribution_count() const noexcept {
    return _element_offsets.back();
}

const Quad4RzGeometry&
SteadyProblem::region_element_geometry(std::size_t region_value,
                                       std::size_t element_index) const {
    return _region_geometries.at(region_value).at(element_index);
}

std::size_t SteadyProblem::contact_count() const noexcept {
    return _definition.contacts.size();
}

const ContactDefinition& SteadyProblem::contact(std::size_t index) const {
    return _definition.contacts.at(index);
}

void SteadyProblem::set_load_factor(double value) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument(
            "SteadyProblem load factor must be finite and nonnegative");
    _load_factor = value;
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        _region_kernels[region_value].set_volumetric_heat_source(
            value * _definition.regions[region_value].volumetric_heat_source);
    }
}

double SteadyProblem::load_factor() const noexcept {
    return _load_factor;
}

std::vector<double> SteadyProblem::initial_state() const {
    std::vector<double> result(dof_count(), 0.0);
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const std::size_t offset = _node_offsets[region_value];
        for (std::size_t node = 0; node < _meshes[region_value].nodes().size();
             ++node) {
            result[_dof_map.temperature(offset + node)] =
                _definition.regions[region_value].initial_temperature;
        }
    }
    for (const DirichletCondition& condition : _dirichlet_conditions)
        result[condition.dof] = condition.value;
    return result;
}

std::size_t SteadyProblem::dof_count() const noexcept {
    return _dof_map.dof_count();
}

std::size_t SteadyProblem::contribution_count() const noexcept {
    return _element_offsets.back() + _thermal_geometries.size() +
           _mechanical_geometries.size();
}

const std::vector<DirichletCondition>&
SteadyProblem::dirichlet_conditions() const noexcept {
    return _dirichlet_conditions;
}

std::pair<std::size_t, std::size_t>
SteadyProblem::element_location(std::size_t contribution_index) const {
    if (contribution_index >= _element_offsets.back())
        throw std::out_of_range(
            "SteadyProblem volume contribution is out of range");
    const auto upper = std::upper_bound(
        _element_offsets.begin(), _element_offsets.end(), contribution_index);
    const std::size_t region_value =
        static_cast<std::size_t>(upper - _element_offsets.begin() - 1);
    return {region_value, contribution_index - _element_offsets[region_value]};
}

LocalDofs
SteadyProblem::contribution_dofs(std::size_t contribution_index) const {
    if (contribution_index < _element_offsets.back()) {
        const auto location = element_location(contribution_index);
        const Quad4Element& element =
            _meshes[location.first].elements().at(location.second);
        std::array<std::size_t, 4> nodes{};
        for (std::size_t node = 0; node < nodes.size(); ++node)
            nodes[node] = global_node(location.first, element.nodes[node]);
        return _dof_map.local_dofs(nodes);
    }
    contribution_index -= _element_offsets.back();
    if (contribution_index < _thermal_nodes.size())
        return _dof_map.local_dofs(_thermal_nodes.at(contribution_index));
    contribution_index -= _thermal_nodes.size();
    return _dof_map.local_dofs(_mechanical_nodes.at(contribution_index));
}

LocalResidual
SteadyProblem::contribution_residual(std::size_t contribution_index,
                                     const LocalValues& state) const {
    if (contribution_index < _element_offsets.back()) {
        const auto location = element_location(contribution_index);
        return _region_kernels[location.first].residual(
            _region_geometries[location.first][location.second], state);
    }
    contribution_index -= _element_offsets.back();
    if (contribution_index < _thermal_geometries.size()) {
        const std::size_t contact_value =
            _thermal_contact_indices[contribution_index];
        return _thermal_kernels[contact_value].residual(
            _thermal_geometries[contribution_index], state);
    }
    contribution_index -= _thermal_geometries.size();
    const std::size_t contact_value =
        _mechanical_contact_indices.at(contribution_index);
    return _mechanical_kernels[contact_value].residual(
        _mechanical_geometries.at(contribution_index), state);
}

LocalSystem
SteadyProblem::linearize_contribution(std::size_t contribution_index,
                                      const LocalValues& state) const {
    if (contribution_index < _element_offsets.back()) {
        const auto location = element_location(contribution_index);
        return _region_kernels[location.first].linearize(
            _region_geometries[location.first][location.second], state);
    }
    contribution_index -= _element_offsets.back();
    if (contribution_index < _thermal_geometries.size()) {
        const std::size_t contact_value =
            _thermal_contact_indices[contribution_index];
        return _thermal_kernels[contact_value].linearize(
            _thermal_geometries[contribution_index], state);
    }
    contribution_index -= _thermal_geometries.size();
    const std::size_t contact_value =
        _mechanical_contact_indices.at(contribution_index);
    return _mechanical_kernels[contact_value].linearize(
        _mechanical_geometries.at(contribution_index), state);
}

std::size_t SteadyProblem::global_node(std::size_t region_value,
                                       std::size_t local_node) const {
    if (local_node >= _meshes.at(region_value).nodes().size())
        throw std::out_of_range("SteadyProblem local node is out of range");
    return _node_offsets.at(region_value) + local_node;
}

SteadyProblem::ResolvedBoundary
SteadyProblem::resolve_boundary(const UnstructuredQuad4Mesh& source_mesh,
                                const std::string& name) const {
    const std::int64_t block_id = source_mesh.side_set_block_id(name);
    const auto found =
        std::find(_block_ids.begin(), _block_ids.end(), block_id);
    if (found == _block_ids.end())
        throw std::invalid_argument(
            "Boundary belongs to an undeclared block: " + name);
    const std::size_t region_value =
        static_cast<std::size_t>(found - _block_ids.begin());
    return {region_value,
            _meshes[region_value].map_side_set(source_mesh, name)};
}

void SteadyProblem::build_volume_geometries() {
    _region_geometries.resize(region_count());
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const RegionMesh& mesh = _meshes[region_value];
        std::vector<Quad4RzGeometry>& geometries =
            _region_geometries[region_value];
        geometries.reserve(mesh.elements().size());
        for (const Quad4Element& element : mesh.elements())
            geometries.push_back(
                make_quad4_rz_geometry(element_coordinates(mesh, element)));
    }
}

void SteadyProblem::build_contacts(const UnstructuredQuad4Mesh& source_mesh) {
    _primary_boundaries.reserve(contact_count());
    _secondary_boundaries.reserve(contact_count());
    _thermal_kernels.reserve(contact_count());
    _mechanical_kernels.reserve(contact_count());

    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        const ContactDefinition& contact_definition =
            _definition.contacts[contact_value];
        ResolvedBoundary primary =
            resolve_boundary(source_mesh, contact_definition.primary);
        ResolvedBoundary secondary =
            resolve_boundary(source_mesh, contact_definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Self-contact is not supported: " +
                                        contact_definition.name);
        validate_connected_radial_boundary(_meshes[primary.region],
                                           primary.boundary,
                                           contact_definition.primary);
        validate_connected_radial_boundary(_meshes[secondary.region],
                                           secondary.boundary,
                                           contact_definition.secondary);

        const RegionMesh& primary_mesh = _meshes[primary.region];
        const RegionMesh& secondary_mesh = _meshes[secondary.region];
        const double primary_radius =
            primary_mesh.nodes().at(primary.boundary.nodes.front()).r;
        const double secondary_radius =
            secondary_mesh.nodes().at(secondary.boundary.nodes.front()).r;
        if (!(primary_radius > secondary_radius))
            throw std::invalid_argument(
                "RZ contact currently requires primary outside secondary: " +
                contact_definition.name);

        _thermal_kernels.emplace_back(GapHeatProperties{
            contact_definition.thermal ? contact_definition.gap_conductivity
                                       : 1.0,
            contact_definition.thermal ? contact_definition.minimum_gap : 1.0});
        _mechanical_kernels.emplace_back(NormalContactProperties{
            contact_definition.mechanical ? contact_definition.penalty : 1.0});

        if (contact_definition.thermal) {
            for (const Line2BoundaryElement& secondary_edge :
                 secondary.boundary.elements) {
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                const double lower_gauss_z =
                    0.5 * (1.0 + gauss) * secondary_coordinates[0].z +
                    0.5 * (1.0 - gauss) * secondary_coordinates[1].z;
                const double upper_gauss_z =
                    0.5 * (1.0 - gauss) * secondary_coordinates[0].z +
                    0.5 * (1.0 + gauss) * secondary_coordinates[1].z;
                const std::size_t lower_segment = find_containing_segment(
                    lower_gauss_z, primary_mesh, primary.boundary.elements);
                const std::size_t upper_segment = find_containing_segment(
                    upper_gauss_z, primary_mesh, primary.boundary.elements);
                if (lower_segment != upper_segment)
                    throw std::invalid_argument(
                        "Each secondary thermal edge must project to one "
                        "primary segment: " +
                        contact_definition.name);
                const Line2BoundaryElement& primary_edge =
                    primary.boundary.elements[lower_segment];
                _thermal_contact_indices.push_back(contact_value);
                _thermal_nodes.push_back({
                    global_node(secondary.region, secondary_edge.nodes[0]),
                    global_node(secondary.region, secondary_edge.nodes[1]),
                    global_node(primary.region, primary_edge.nodes[0]),
                    global_node(primary.region, primary_edge.nodes[1]),
                });
                _thermal_geometries.push_back(make_line2_rz_heat_geometry(
                    secondary_coordinates,
                    edge_coordinates(primary_mesh, primary_edge)));
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
                    const double z = secondary_coordinates[secondary_node].z;
                    const std::size_t containing = find_containing_segment(
                        z, primary_mesh, primary.boundary.elements);
                    const std::size_t first =
                        containing == 0 ? 0 : containing - 1;
                    const std::size_t last = std::min(
                        containing + 1, primary.boundary.elements.size() - 1);
                    for (std::size_t candidate = first; candidate <= last;
                         ++candidate) {
                        const Line2BoundaryElement& primary_edge =
                            primary.boundary.elements[candidate];
                        _mechanical_contact_indices.push_back(contact_value);
                        _mechanical_nodes.push_back({
                            global_node(secondary.region,
                                        secondary_edge.nodes[0]),
                            global_node(secondary.region,
                                        secondary_edge.nodes[1]),
                            global_node(primary.region, primary_edge.nodes[0]),
                            global_node(primary.region, primary_edge.nodes[1]),
                        });
                        _mechanical_geometries.push_back(
                            make_node_to_line_rz_contact_geometry(
                                secondary_coordinates,
                                edge_coordinates(primary_mesh, primary_edge),
                                secondary_node,
                                candidate + 1 ==
                                    primary.boundary.elements.size()));
                        _mechanical_secondary_indices.push_back(edge_index +
                                                                secondary_node);
                    }
                }
            }
        }
        _primary_boundaries.push_back(std::move(primary));
        _secondary_boundaries.push_back(std::move(secondary));
    }
}

void SteadyProblem::build_boundary_conditions(
    const UnstructuredQuad4Mesh& source_mesh) {
    for (const BoundaryConditionDefinition& definition :
         _definition.boundary_conditions) {
        if (definition.name.empty() || definition.boundary.empty() ||
            !std::isfinite(definition.value))
            throw std::invalid_argument("Boundary-condition names, boundaries, "
                                        "and values must be valid");
        ResolvedBoundary resolved =
            resolve_boundary(source_mesh, definition.boundary);
        if (definition.type == BoundaryConditionType::dirichlet) {
            for (std::size_t local_node : resolved.boundary.nodes) {
                _dirichlet_conditions.push_back(
                    {_dof_map.dof(definition.field,
                                  global_node(resolved.region, local_node)),
                     definition.value});
            }
        } else {
            if (definition.value < 0.0)
                throw std::invalid_argument(
                    "Pressure boundary conditions must be nonnegative");
            if (resolved.boundary.kind != RegionBoundaryKind::radial_inner &&
                resolved.boundary.kind != RegionBoundaryKind::radial_outer)
                throw std::invalid_argument(
                    "Pressure currently requires a radial boundary: " +
                    definition.boundary);
            _pressure_loads.push_back({resolved.region,
                                       std::move(resolved.boundary),
                                       definition.value});
        }
    }
    validate_dirichlet_conditions(_dirichlet_conditions);
}

void SteadyProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    add_pressure_residual(residual);
}

void SteadyProblem::add_external_residual(std::vector<double>& residual) const {
    add_pressure_residual(residual);
}

void SteadyProblem::add_pressure_residual(std::vector<double>& residual) const {
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const PressureLoad& load : _pressure_loads) {
        if (load.pressure == 0.0)
            continue;
        const double normal_r =
            load.boundary.kind == RegionBoundaryKind::radial_inner ? -1.0 : 1.0;
        const RegionMesh& mesh = _meshes[load.region];
        for (const Line2BoundaryElement& edge : load.boundary.elements) {
            const RzPoint& first = mesh.nodes().at(edge.nodes[0]);
            const RzPoint& second = mesh.nodes().at(edge.nodes[1]);
            const double dr = second.r - first.r;
            const double dz = second.z - first.z;
            const double line_jacobian = 0.5 * std::sqrt(dr * dr + dz * dz);
            for (double xi : locations) {
                const std::array<double, 2> shape = {0.5 * (1.0 - xi),
                                                     0.5 * (1.0 + xi)};
                const double radius = shape[0] * first.r + shape[1] * second.r;
                const double measure = 2.0 * pi * radius * line_jacobian;
                for (std::size_t node = 0; node < 2; ++node) {
                    const std::size_t dof = _dof_map.radial_displacement(
                        global_node(load.region, edge.nodes[node]));
                    residual[dof] +=
                        measure * load.pressure * normal_r * shape[node];
                }
            }
        }
    }
}

std::vector<ContactNodeSummary>
SteadyProblem::summarize_contact_nodes(std::size_t contact_value,
                                       const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SteadyProblem contact summary state size mismatch");
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _meshes[secondary.region];
    std::vector<ContactNodeSummary> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) {
        result.push_back({mesh.nodes().at(node).z, false,
                          std::numeric_limits<double>::infinity(), 0.0, 0.0,
                          0.0, 0.0});
    }
    const std::size_t first_mechanical =
        _element_offsets.back() + _thermal_geometries.size();
    for (std::size_t contribution = 0;
         contribution < _mechanical_geometries.size(); ++contribution) {
        if (_mechanical_contact_indices[contribution] != contact_value)
            continue;
        const LocalValues local_state =
            contribution_state(first_mechanical + contribution, state);
        const ContactPointValue value =
            _mechanical_kernels[contact_value].value(
                _mechanical_geometries[contribution], local_state);
        if (!value.projected)
            continue;
        ContactNodeSummary& node =
            result.at(_mechanical_secondary_indices[contribution]);
        node.projected = true;
        node.gap = std::min(node.gap, value.gap);
        node.tributary_area += value.tributary_area;
        node.tributary_length += value.tributary_length;
        node.contact_force += value.contact_force;
    }
    for (ContactNodeSummary& node : result) {
        if (node.tributary_area > 0.0)
            node.pressure = node.contact_force / node.tributary_area;
    }
    return result;
}

InterfaceSummary
SteadyProblem::summarize_interface(std::size_t contact_value,
                                   const std::vector<double>& state) const {
    if (contact_value >= contact_count())
        throw std::out_of_range("SteadyProblem contact index is out of range");
    InterfaceSummary summary = {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        0,
        0,
        0.0,
    };
    const std::size_t first_thermal = _element_offsets.back();
    bool has_thermal = false;
    for (std::size_t contribution = 0;
         contribution < _thermal_geometries.size(); ++contribution) {
        if (_thermal_contact_indices[contribution] != contact_value)
            continue;
        has_thermal = true;
        const LocalValues local_state =
            contribution_state(first_thermal + contribution, state);
        const HeatQuadratureValues values =
            _thermal_kernels[contact_value].quadrature_values(
                _thermal_geometries[contribution], local_state);
        for (const HeatQuadratureValue& value : values) {
            summary.minimum_gap = std::min(summary.minimum_gap, value.gap);
            summary.maximum_gap = std::max(summary.maximum_gap, value.gap);
            summary.total_heat_rate += value.weighted_measure * value.heat_flux;
        }
    }
    if (!has_thermal) {
        summary.minimum_gap = 0.0;
        summary.maximum_gap = 0.0;
    }

    bool has_mechanical = false;
    for (std::size_t index : _mechanical_contact_indices) {
        if (index == contact_value) {
            has_mechanical = true;
            break;
        }
    }
    if (has_mechanical) {
        for (const ContactNodeSummary& node :
             summarize_contact_nodes(contact_value, state)) {
            if (!node.projected)
                continue;
            ++summary.projected_contact_nodes;
            summary.minimum_contact_gap =
                std::min(summary.minimum_contact_gap, node.gap);
            summary.maximum_contact_pressure =
                std::max(summary.maximum_contact_pressure, node.pressure);
            summary.total_contact_force += node.contact_force;
            if (node.pressure > 0.0) {
                ++summary.active_contact_nodes;
                summary.active_contact_length += node.tributary_length;
            }
        }
    } else {
        summary.minimum_contact_gap = 0.0;
    }
    return summary;
}

} // namespace fuelsim
