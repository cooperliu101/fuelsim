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

std::pair<std::size_t, std::array<std::size_t, 2>>
edge_parent(const RegionMesh& mesh, const Line2BoundaryElement& edge) {
    std::size_t parent = mesh.elements().size();
    std::array<std::size_t, 2> local_nodes{};
    for (std::size_t element_index = 0; element_index < mesh.elements().size();
         ++element_index) {
        const Quad4Element& element = mesh.elements()[element_index];
        std::array<std::size_t, 2> candidate{};
        bool contains = true;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
            const auto found =
                std::find(element.nodes.begin(), element.nodes.end(),
                          edge.nodes[edge_node]);
            if (found == element.nodes.end()) {
                contains = false;
                break;
            }
            candidate[edge_node] =
                static_cast<std::size_t>(found - element.nodes.begin());
        }
        if (!contains)
            continue;
        if (parent != mesh.elements().size())
            throw std::invalid_argument(
                "Boundary edge has more than one adjacent region element");
        parent = element_index;
        local_nodes = candidate;
    }
    if (parent == mesh.elements().size())
        throw std::invalid_argument(
            "Boundary edge has no adjacent region element");
    return {parent, local_nodes};
}

void validate_definitions(const SteadyProblemDefinition& definition) {
    if (definition.regions.empty())
        throw std::invalid_argument(
            "SteadyProblem requires at least one region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        const RegionDefinition& value = definition.regions[region];
        if (value.name.empty() || (value.block.empty() && value.block_id < 0) ||
            (!value.block.empty() && value.block_id >= 0))
            throw std::invalid_argument(
                "SteadyProblem regions require a name and exactly one block "
                "selector");
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
    return ((point.r - segment[0].r) * dr +
            (point.z - segment[0].z) * dz) /
           (dr * dr + dz * dz);
}

bool projection_interval(
    const Line2InterfaceSideCoordinates& secondary,
    const Line2InterfaceSideCoordinates& primary, double& lower,
    double& upper) {
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
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon();
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

std::size_t closest_primary_segment(
    const RzPoint& point, const RegionMesh& primary_mesh,
    const std::vector<Line2BoundaryElement>& primary_edges) {
    std::size_t selected = primary_edges.size();
    double minimum_distance = std::numeric_limits<double>::infinity();
    constexpr double tolerance = 1.0e-12;
    for (std::size_t edge = 0; edge < primary_edges.size(); ++edge) {
        const Line2InterfaceSideCoordinates coordinates =
            edge_coordinates(primary_mesh, primary_edges[edge]);
        const double fraction = projection_fraction(point, coordinates);
        if (fraction < -tolerance || fraction > 1.0 + tolerance)
            continue;
        const double clipped = std::max(0.0, std::min(1.0, fraction));
        const double dr = coordinates[0].r +
                              clipped *
                                  (coordinates[1].r - coordinates[0].r) -
                          point.r;
        const double dz = coordinates[0].z +
                              clipped *
                                  (coordinates[1].z - coordinates[0].z) -
                          point.z;
        const double distance = std::hypot(dr, dz);
        if (distance < minimum_distance) {
            selected = edge;
            minimum_distance = distance;
        }
    }
    if (selected == primary_edges.size())
        throw std::invalid_argument(
            "Contact node is outside the primary surface projection");
    return selected;
}

} // namespace

std::vector<std::int64_t>
SteadyProblem::resolve_block_ids(const SteadyProblemDefinition& definition,
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
SteadyProblem::build_meshes(const SteadyProblemDefinition& definition,
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
      _load_factor(1.0), _time(0.0) {
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
    refresh_controlled_values();
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
    refresh_controlled_values();
}

double SteadyProblem::load_factor() const noexcept {
    return _load_factor;
}

void SteadyProblem::set_time(double value) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument(
            "SteadyProblem time must be finite and nonnegative");
    _time = value;
    refresh_controlled_values();
}

double SteadyProblem::time() const noexcept {
    return _time;
}

double SteadyProblem::function_value(const std::string& name) const {
    const auto found = std::find_if(
        _definition.time_tables.begin(), _definition.time_tables.end(),
        [&name](const PiecewiseLinearTimeTable& table) {
            return table.name() == name;
        });
    if (found == _definition.time_tables.end())
        throw std::invalid_argument("Unknown time-table function: " + name);
    return found->value(_time);
}

double SteadyProblem::load_multiplier(bool scale_with_load,
                                      const std::string& function) const {
    if (!function.empty())
        return function_value(function);
    return scale_with_load ? _load_factor : 1.0;
}

void SteadyProblem::refresh_controlled_values() {
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const RegionDefinition& region = _definition.regions[region_value];
        const double multiplier =
            region.heat_source_function.empty()
                ? _load_factor
                : function_value(region.heat_source_function);
        _region_kernels[region_value].set_volumetric_heat_source(
            multiplier * region.volumetric_heat_source);
    }
    for (const ControlledDirichlet& controlled :
         _controlled_dirichlet_conditions) {
        const auto condition = std::lower_bound(
            _dirichlet_conditions.begin(), _dirichlet_conditions.end(),
            controlled.dof,
            [](const DirichletCondition& candidate, std::size_t dof) {
                return candidate.dof < dof;
            });
        if (condition == _dirichlet_conditions.end() ||
            condition->dof != controlled.dof)
            throw std::logic_error(
                "SteadyProblem controlled Dirichlet mapping is invalid");
        condition->value =
            load_multiplier(controlled.scale_with_load, controlled.function) *
            controlled.value;
    }
    for (std::size_t load = 0; load < _convection_loads.size(); ++load) {
        const ConvectionLoad& convection = _convection_loads[load];
        const double coefficient_multiplier =
            convection.coefficient_function.empty()
                ? 1.0
                : function_value(convection.coefficient_function);
        const double ambient_multiplier =
            convection.ambient_temperature_function.empty()
                ? 1.0
                : function_value(convection.ambient_temperature_function);
        _convection_kernels[load].set_properties(
            {coefficient_multiplier * convection.heat_transfer_coefficient,
             ambient_multiplier * convection.ambient_temperature});
    }
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
           _mechanical_geometries.size() + _convection_geometries.size();
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
    if (contribution_index < _mechanical_nodes.size())
        return _dof_map.local_dofs(_mechanical_nodes.at(contribution_index));
    contribution_index -= _mechanical_nodes.size();
    return _dof_map.local_dofs(_convection_nodes.at(contribution_index));
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
    if (contribution_index >= _mechanical_geometries.size()) {
        contribution_index -= _mechanical_geometries.size();
        const std::size_t load =
            _convection_load_indices.at(contribution_index);
        return _convection_kernels[load].residual(
            _convection_geometries.at(contribution_index), state);
    }
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
    if (contribution_index >= _mechanical_geometries.size()) {
        contribution_index -= _mechanical_geometries.size();
        const std::size_t load =
            _convection_load_indices.at(contribution_index);
        return _convection_kernels[load].linearize(
            _convection_geometries.at(contribution_index), state);
    }
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
        const RegionMesh& primary_mesh = _meshes[primary.region];
        const RegionMesh& secondary_mesh = _meshes[secondary.region];
        primary.boundary = ordered_connected_boundary(
            primary_mesh, std::move(primary.boundary),
            contact_definition.primary);
        secondary.boundary = ordered_connected_boundary(
            secondary_mesh, std::move(secondary.boundary),
            contact_definition.secondary);

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
                        projections.push_back(
                            {primary_edge, lower, upper});
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
                    _thermal_contact_indices.push_back(contact_value);
                    _thermal_nodes.push_back({
                        global_node(secondary.region, secondary_edge.nodes[0]),
                        global_node(secondary.region, secondary_edge.nodes[1]),
                        global_node(primary.region, primary_edge.nodes[0]),
                        global_node(primary.region, primary_edge.nodes[1]),
                    });
                    _thermal_geometries.push_back(make_line2_rz_heat_geometry(
                        secondary_coordinates,
                        edge_coordinates(primary_mesh, primary_edge),
                        projection.lower, projection.upper));
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
                    const std::size_t containing = closest_primary_segment(
                        secondary_coordinates[secondary_node], primary_mesh,
                        primary.boundary.elements);
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
        if (definition.scale_with_load && !definition.function.empty())
            throw std::invalid_argument(
                "Boundary condition cannot combine scale_with_load and a "
                "time function: " +
                definition.name);
        if (definition.type == BoundaryConditionType::dirichlet) {
            for (std::size_t local_node : resolved.boundary.nodes) {
                const std::size_t dof = _dof_map.dof(
                    definition.field, global_node(resolved.region, local_node));
                _dirichlet_conditions.push_back(
                    {dof, load_multiplier(definition.scale_with_load,
                                          definition.function) *
                              definition.value});
                if (definition.scale_with_load || !definition.function.empty())
                    _controlled_dirichlet_conditions.push_back(
                        {dof, definition.value, definition.scale_with_load,
                         definition.function});
            }
        } else if (definition.type == BoundaryConditionType::pressure) {
            if (definition.value < 0.0)
                throw std::invalid_argument(
                    "Pressure boundary conditions must be nonnegative");
            if (resolved.boundary.kind != RegionBoundaryKind::radial_inner &&
                resolved.boundary.kind != RegionBoundaryKind::radial_outer)
                throw std::invalid_argument(
                    "Pressure currently requires a radial boundary: " +
                    definition.boundary);
            _pressure_loads.push_back(
                {resolved.region, std::move(resolved.boundary),
                 definition.value, definition.scale_with_load,
                 definition.function});
        } else if (definition.type == BoundaryConditionType::traction) {
            if (definition.field == Field::temperature)
                throw std::invalid_argument(
                    "Traction requires a displacement field: " +
                    definition.boundary);
            _traction_loads.push_back(
                {resolved.region, std::move(resolved.boundary),
                 definition.field, definition.value, definition.scale_with_load,
                 definition.function});
        } else {
            if (!(definition.heat_transfer_coefficient > 0.0) ||
                !(definition.ambient_temperature > 0.0))
                throw std::invalid_argument(
                    "Convection coefficient and ambient temperature must be "
                    "positive: " +
                    definition.name);
            const std::size_t load = _convection_loads.size();
            _convection_loads.push_back(
                {definition.heat_transfer_coefficient,
                 definition.ambient_temperature,
                 definition.coefficient_function,
                 definition.ambient_temperature_function});
            _convection_kernels.emplace_back(
                ConvectionProperties{definition.heat_transfer_coefficient,
                                     definition.ambient_temperature});
            const RegionMesh& mesh = _meshes[resolved.region];
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const auto parent = edge_parent(mesh, edge);
                const Quad4Element& element = mesh.elements().at(parent.first);
                std::array<std::size_t, 4> nodes{};
                for (std::size_t node = 0; node < nodes.size(); ++node)
                    nodes[node] =
                        global_node(resolved.region, element.nodes[node]);
                _convection_nodes.push_back(nodes);
                _convection_geometries.push_back(
                    make_line2_rz_convection_geometry(
                        {{mesh.nodes().at(edge.nodes[0]),
                          mesh.nodes().at(edge.nodes[1])}},
                        parent.second));
                _convection_load_indices.push_back(load);
            }
        }
    }
    validate_dirichlet_conditions(_dirichlet_conditions);
}

void SteadyProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    add_pressure_residual(residual);
    add_traction_residual(residual);
}

void SteadyProblem::add_external_residual(std::vector<double>& residual) const {
    add_pressure_residual(residual);
    add_traction_residual(residual);
}

void SteadyProblem::add_pressure_residual(std::vector<double>& residual) const {
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const PressureLoad& load : _pressure_loads) {
        const double pressure =
            load_multiplier(load.scale_with_load, load.function) *
            load.pressure;
        if (pressure < 0.0)
            throw std::domain_error(
                "Pressure time function produced a negative load");
        if (pressure == 0.0)
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
                        measure * pressure * normal_r * shape[node];
                }
            }
        }
    }
}

void SteadyProblem::add_traction_residual(std::vector<double>& residual) const {
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const TractionLoad& load : _traction_loads) {
        const double traction =
            load_multiplier(load.scale_with_load, load.function) *
            load.traction;
        if (traction == 0.0)
            continue;
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
                    const std::size_t dof = _dof_map.dof(
                        load.field, global_node(load.region, edge.nodes[node]));
                    residual[dof] -= measure * traction * shape[node];
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
        result.push_back({mesh.nodes().at(node).r, mesh.nodes().at(node).z,
                          false,
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

std::vector<std::size_t>
SteadyProblem::contact_secondary_source_nodes(std::size_t contact_value) const {
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<std::size_t> result;
    result.reserve(secondary.boundary.nodes.size());
    for (const std::size_t node : secondary.boundary.nodes)
        result.push_back(mesh.source_node_ids().at(node));
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
