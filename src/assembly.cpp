#include "assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
namespace fuelsim::spatial_detail {
namespace {
void validate_definition(const SpatialDefinition& definition, bool allow_contacts, bool allow_finite_strain) {
    if (definition.regions.empty()) throw std::invalid_argument("SpatialAssembly requires at least one region");
    if (!allow_contacts && !definition.contacts.empty())
        throw std::invalid_argument("Cartesian three-dimensional stage B does not support contact");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        const RegionDefinition& value = definition.regions[region];
        if (value.name.empty() || (value.block.empty() && value.block_id < 0) ||
            (!value.block.empty() && value.block_id >= 0))
            throw std::invalid_argument("SpatialAssembly regions require a name and exactly one block selector");
        if (!allow_finite_strain && value.strain_formulation != StrainFormulation::small)
            throw std::invalid_argument("Cartesian three-dimensional stage B supports only small strain");
        for (std::size_t previous = 0; previous < region; ++previous) {
            if (definition.regions[previous].name == value.name)
                throw std::invalid_argument("Duplicate region name: " + value.name);
            if ((!value.block.empty() && definition.regions[previous].block == value.block) ||
                (value.block_id >= 0 && definition.regions[previous].block_id == value.block_id))
                throw std::invalid_argument(
                    "Duplicate region block: " + (value.block.empty() ? std::to_string(value.block_id) : value.block));
        }
    }
    for (std::size_t table = 0; table < definition.time_tables.size(); ++table)
        for (std::size_t previous = 0; previous < table; ++previous)
            if (definition.time_tables[previous].name() == definition.time_tables[table].name())
                throw std::invalid_argument("Duplicate time-table name: " + definition.time_tables[table].name());
    for (std::size_t contact = 0; contact < definition.contacts.size(); ++contact) {
        const ContactDefinition& value = definition.contacts[contact];
        if (value.name.empty() || value.primary.empty() || value.secondary.empty())
            throw std::invalid_argument("SpatialAssembly contact names and boundaries must not be empty");
        if (!value.thermal && !value.mechanical)
            throw std::invalid_argument("Contact must enable thermal or mechanical coupling: " + value.name);
        if (value.primary == value.secondary)
            throw std::invalid_argument("Contact primary and secondary must differ: " + value.name);
        for (std::size_t previous = 0; previous < contact; ++previous) {
            if (definition.contacts[previous].name == value.name)
                throw std::invalid_argument("Duplicate contact name: " + value.name);
            if (value.mechanical && definition.contacts[previous].mechanical &&
                definition.contacts[previous].secondary == value.secondary)
                throw std::invalid_argument("Mechanical secondary boundary is reused: " + value.secondary);
        }
    }
}
} // namespace
std::vector<std::int64_t> resolve_block_ids(const SpatialDefinition& definition,
    const UnstructuredMeshMetadata& source_mesh, bool allow_contacts, bool allow_finite_strain) {
    validate_definition(definition, allow_contacts, allow_finite_strain);
    std::vector<std::int64_t> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        const std::int64_t id = region.block_id >= 0 ? region.block_id : source_mesh.element_block(region.block).id;
        const bool known = std::any_of(source_mesh.element_blocks().begin(), source_mesh.element_blocks().end(),
            [id](const ElementBlockInfo& block) { return block.id == id; });
        if (!known) throw std::invalid_argument("Unknown element block ID: " + std::to_string(id));
        if (std::find(result.begin(), result.end(), id) != result.end())
            throw std::invalid_argument("Duplicate resolved region block ID: " + std::to_string(id));
        result.push_back(id);
    }
    return result;
}
void validate_dirichlet_conditions(
    std::vector<DirichletCondition>& conditions, const char* conflict_message, const char* duplicate_message) {
    std::sort(conditions.begin(), conditions.end(),
        [](const DirichletCondition& lhs, const DirichletCondition& rhs) { return lhs.dof < rhs.dof; });
    for (std::size_t index = 1; index < conditions.size(); ++index) {
        if (conditions[index - 1].dof != conditions[index].dof) continue;
        if (conditions[index - 1].value != conditions[index].value) throw std::invalid_argument(conflict_message);
        throw std::invalid_argument(duplicate_message);
    }
}
double function_value(const SpatialDefinition& definition, double time, const std::string& name) {
    const auto found = std::find_if(definition.time_tables.begin(), definition.time_tables.end(),
        [&name](const PiecewiseLinearTimeTable& table) { return table.name() == name; });
    if (found == definition.time_tables.end()) throw std::invalid_argument("Unknown time-table function: " + name);
    return found->value(time);
}
double controlled_value(const SpatialDefinition& definition, double time, double load_factor, double value,
    bool scale_with_load, const std::string& function) {
    if (!function.empty()) return value * function_value(definition, time, function);
    return value * (scale_with_load ? load_factor : 1.0);
}
SpatialLayout::SpatialLayout(SpatialDefinition definition, std::vector<std::int64_t> block_ids, DofLayout layout)
    : _definition(std::move(definition)), _block_ids(std::move(block_ids)), _layout(layout) {}
void SpatialLayout::initialize_counts(
    const std::vector<std::size_t>& node_counts, const std::vector<std::size_t>& element_counts) {
    if (node_counts.size() != region_count() || element_counts.size() != region_count())
        throw std::invalid_argument("Spatial layout counts must match the declared regions");
    _node_offsets = {0};
    _element_offsets = {0};
    for (std::size_t region = 0; region < region_count(); ++region) {
        if (node_counts[region] > std::numeric_limits<std::size_t>::max() - _node_offsets.back() ||
            element_counts[region] > std::numeric_limits<std::size_t>::max() - _element_offsets.back())
            throw std::length_error("Spatial layout count overflows");
        _node_offsets.push_back(_node_offsets.back() + node_counts[region]);
        _element_offsets.push_back(_element_offsets.back() + element_counts[region]);
    }
    _dof_map.emplace(_node_offsets.back(), _layout);
}
std::size_t SpatialLayout::region_node_offset(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Spatial layout region index is out of range");
    return _node_offsets[index];
}
std::size_t SpatialLayout::region_element_count(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Spatial layout region index is out of range");
    return _element_offsets[index + 1] - _element_offsets[index];
}
std::size_t SpatialLayout::region_element_offset(std::size_t index) const {
    if (index >= region_count()) throw std::out_of_range("Spatial layout region index is out of range");
    return _element_offsets[index];
}
std::pair<std::size_t, std::size_t> SpatialLayout::element_location(std::size_t index) const {
    if (index >= volume_contribution_count()) throw std::out_of_range("Volume contribution is out of range");
    const auto upper = std::upper_bound(_element_offsets.begin(), _element_offsets.end(), index);
    const std::size_t region = static_cast<std::size_t>(upper - _element_offsets.begin() - 1);
    return {region, index - _element_offsets[region]};
}
double SpatialLayout::region_heat_source(std::size_t index) const {
    const RegionDefinition& value = region(index);
    return controlled_value(
        _definition, _time, _load_factor, value.volumetric_heat_source, true, value.heat_source_function);
}
std::vector<double> SpatialLayout::initial_state() const {
    std::vector<double> result(dof_count(), 0.0);
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index)
        for (std::size_t node = _node_offsets[region_index]; node < _node_offsets[region_index + 1]; ++node)
            result[dof_map().temperature(node)] = region(region_index).initial_temperature;
    for (const DirichletCondition& condition : _dirichlet_conditions) result[condition.dof] = condition.value;
    return result;
}
std::size_t SpatialLayout::global_node(std::size_t region_index, std::size_t local_node) const {
    if (region_index >= region_count()) throw std::out_of_range("Spatial layout region index is out of range");
    if (local_node >= _node_offsets[region_index + 1] - _node_offsets[region_index])
        throw std::out_of_range("Spatial layout local node is out of range");
    return _node_offsets[region_index] + local_node;
}
void SpatialLayout::add_dirichlet(std::size_t dof, std::size_t boundary_index) {
    const BoundaryConditionDefinition& boundary = _definition.boundary_conditions.at(boundary_index);
    _dirichlet_conditions.push_back({dof, controlled_value(_definition, _time, _load_factor, boundary.value,
                                              boundary.scale_with_load, boundary.function)});
    if (boundary.scale_with_load || !boundary.function.empty())
        _controlled_dirichlet_conditions.push_back({dof, boundary_index});
}
void SpatialLayout::set_load_factor_value(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Spatial layout load factor must be finite");
    _load_factor = value;
}
void SpatialLayout::set_time_value(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Spatial layout time must be finite");
    _time = value;
}
void SpatialLayout::refresh_dirichlet_values() {
    for (const ControlledDirichlet& controlled : _controlled_dirichlet_conditions) {
        const auto condition = std::lower_bound(_dirichlet_conditions.begin(), _dirichlet_conditions.end(),
            controlled.dof, [](const DirichletCondition& candidate, std::size_t dof) { return candidate.dof < dof; });
        if (condition == _dirichlet_conditions.end() || condition->dof != controlled.dof)
            throw std::logic_error("Spatial layout controlled Dirichlet mapping is invalid");
        const BoundaryConditionDefinition& boundary = _definition.boundary_conditions.at(controlled.boundary_index);
        condition->value = controlled_value(
            _definition, _time, _load_factor, boundary.value, boundary.scale_with_load, boundary.function);
    }
}
ConvectionValues SpatialLayout::convection_values(const BoundaryConditionDefinition& boundary) const {
    const double coefficient =
        boundary.heat_transfer_coefficient *
        (boundary.coefficient_function.empty() ? 1.0
                                               : function_value(_definition, _time, boundary.coefficient_function));
    const double ambient = boundary.ambient_temperature *
                           (boundary.ambient_temperature_function.empty()
                                   ? 1.0
                                   : function_value(_definition, _time, boundary.ambient_temperature_function));
    return {coefficient, ambient};
}
} // namespace fuelsim::spatial_detail
namespace fuelsim::rz {
namespace {
Quad4Coordinates element_coordinates(const RegionMesh& mesh, const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}
void check_state_size(std::size_t state_size, std::size_t dof_count, const char* message) {
    if (state_size != dof_count) throw std::invalid_argument(message);
}
} // namespace
const Quad4RzGeometry& SpatialAssembly::region_element_geometry(std::size_t region, std::size_t element) const {
    return _region_geometries.at(region).at(element);
}
ResolvedBoundary SpatialAssembly::resolve_boundary(
    const UnstructuredQuad4Mesh& source_mesh, const std::string& name) const {
    const std::int64_t block_id = source_mesh.side_set_block_id(name);
    const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
    if (found == _block_ids.end()) throw std::invalid_argument("Boundary belongs to an undeclared block: " + name);
    const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
    return {region, _meshes[region].map_side_set(source_mesh, name)};
}
std::pair<std::size_t, std::array<std::size_t, 2>> SpatialAssembly::edge_parent(
    std::size_t region, const Line2BoundaryElement& edge) const {
    const RegionMesh& mesh = region_mesh(region);
    std::size_t parent = mesh.elements().size();
    std::array<std::size_t, 2> local_nodes{};
    for (std::size_t element_index = 0; element_index < mesh.elements().size(); ++element_index) {
        const Quad4Element& element = mesh.elements()[element_index];
        std::array<std::size_t, 2> candidate{};
        bool contains = false;
        for (std::size_t side = 0; side < element.nodes.size(); ++side) {
            const std::size_t next = (side + 1U) % element.nodes.size();
            const bool forward = element.nodes[side] == edge.nodes[0] && element.nodes[next] == edge.nodes[1],
                       reverse = element.nodes[side] == edge.nodes[1] && element.nodes[next] == edge.nodes[0];
            if (forward || reverse) {
                candidate = {{side, next}};
                contains = true;
                break;
            }
        }
        if (!contains) continue;
        if (parent != mesh.elements().size())
            throw std::invalid_argument("Boundary edge has more than one adjacent region element");
        parent = element_index;
        local_nodes = candidate;
    }
    if (parent == mesh.elements().size()) throw std::invalid_argument("Boundary edge has no adjacent region element");
    return {parent, local_nodes};
}
void SpatialAssembly::build_volume_geometries() {
    _region_geometries.resize(region_count());
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index) {
        const RegionMesh& mesh = _meshes[region_index];
        std::vector<Quad4RzGeometry>& geometries = _region_geometries[region_index];
        geometries.reserve(mesh.elements().size());
        for (const Quad4Element& element : mesh.elements())
            geometries.push_back(make_quad4_rz_geometry(element_coordinates(mesh, element)));
    }
}
SpatialAssembly::ContributionRanges SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count(),
                      mechanical_begin = thermal_begin + _thermal_contributions.size(),
                      boundary_begin = mechanical_begin + _mechanical_contributions.size();
    return {thermal_begin, mechanical_begin, boundary_begin, boundary_begin + _boundary_contributions.size()};
}
SpatialAssembly::ContributionLocation SpatialAssembly::locate_contribution(std::size_t index) const {
    const ContributionRanges ranges = contribution_ranges();
    if (index >= ranges.end) throw std::out_of_range("SpatialAssembly contribution index is out of range");
    if (index < ranges.thermal_begin) return {SpatialContributionType::volume, index};
    if (index < ranges.mechanical_begin)
        return {SpatialContributionType::thermal_contact, index - ranges.thermal_begin};
    if (index < ranges.boundary_begin)
        return {SpatialContributionType::mechanical_contact, index - ranges.mechanical_begin};
    const std::size_t local = index - ranges.boundary_begin;
    return {_boundary_contributions[local].type, local};
}
SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    return locate_contribution(index).type;
}
namespace {
struct BoundaryEdgeData final {
    std::array<std::size_t, 4> nodes;
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};
BoundaryEdgeData boundary_edge_data(
    const SpatialAssembly& layout, const ResolvedBoundary& boundary, const Line2BoundaryElement& edge) {
    const RegionMesh& mesh = layout.region_mesh(boundary.region);
    const auto parent = layout.edge_parent(boundary.region, edge);
    const Quad4Element& element = mesh.elements().at(parent.first);
    BoundaryEdgeData result{{}, {}, parent.second};
    for (std::size_t node = 0; node < result.nodes.size(); ++node)
        result.nodes[node] = layout.global_node(boundary.region, element.nodes[node]);
    for (std::size_t node = 0; node < result.coordinates.size(); ++node)
        result.coordinates[node] = mesh.nodes().at(element.nodes[parent.second[node]]);
    return result;
}
} // namespace
void SpatialAssembly::build_boundaries(const UnstructuredQuad4Mesh& source_mesh) {
    for (std::size_t boundary_index = 0; boundary_index < _definition.boundary_conditions.size(); ++boundary_index) {
        const BoundaryConditionDefinition& definition = _definition.boundary_conditions[boundary_index];
        if (definition.name.empty() || definition.boundary.empty())
            throw std::invalid_argument("Boundary-condition names and boundaries must be valid");
        const ResolvedBoundary resolved = resolve_boundary(source_mesh, definition.boundary);
        if (definition.scale_with_load && !definition.function.empty())
            throw std::invalid_argument(
                "Boundary condition cannot combine scale_with_load and a time function: " + definition.name);
        if (definition.type == BoundaryConditionType::dirichlet) {
            for (const std::size_t local_node : resolved.boundary.nodes) {
                const std::size_t dof = dof_map().dof(definition.field, global_node(resolved.region, local_node));
                add_dirichlet(dof, boundary_index);
            }
            continue;
        }
        const std::size_t kernel = _boundary_data.size();
        SpatialContributionType type;
        if (definition.type == BoundaryConditionType::pressure) {
            type = SpatialContributionType::pressure;
            const bool displaced = region(resolved.region).strain_formulation == StrainFormulation::finite;
            _boundary_data.push_back(
                make_line2_rz_pressure_data(spatial_detail::controlled_value(_definition, _time, _load_factor,
                                                definition.value, definition.scale_with_load, definition.function),
                    displaced));
        } else if (definition.type == BoundaryConditionType::traction) {
            type = SpatialContributionType::traction;
            if (definition.field == Field::temperature)
                throw std::invalid_argument("Traction requires a displacement field: " + definition.boundary);
            if (definition.use_displaced_geometry &&
                region(resolved.region).strain_formulation != StrainFormulation::finite)
                throw std::invalid_argument(
                    "Current-configuration traction requires finite strain: " + definition.name);
            _boundary_data.push_back(make_line2_rz_traction_data(
                definition.field == Field::radial_displacement ? TractionComponent::radial : TractionComponent::axial,
                spatial_detail::controlled_value(_definition, _time, _load_factor, definition.value,
                    definition.scale_with_load, definition.function),
                definition.use_displaced_geometry));
        } else {
            type = SpatialContributionType::convection;
            _boundary_data.push_back(
                make_line2_rz_convection_data(definition.heat_transfer_coefficient, definition.ambient_temperature));
        }
        _boundary_definition_indices.push_back(boundary_index);
        for (const Line2BoundaryElement& edge : resolved.boundary.elements) {
            const BoundaryEdgeData data = boundary_edge_data(*this, resolved, edge);
            _boundary_contributions.push_back(
                {type, kernel, data.nodes, make_line2_rz_boundary_geometry(data.coordinates, data.local_nodes)});
        }
    }
    spatial_detail::validate_dirichlet_conditions(_dirichlet_conditions,
        "BoundaryAssembly has conflicting Dirichlet conditions", "BoundaryAssembly has duplicate Dirichlet conditions");
    refresh_controlled_values();
}
void SpatialAssembly::set_load_factor(double value) {
    set_load_factor_value(value);
    refresh_controlled_values();
}
void SpatialAssembly::set_time(double value) {
    set_time_value(value);
    refresh_controlled_values();
}
void SpatialAssembly::refresh_controlled_values() {
    refresh_dirichlet_values();
    for (std::size_t load = 0; load < _boundary_definition_indices.size(); ++load) {
        const BoundaryConditionDefinition& boundary =
            _definition.boundary_conditions[_boundary_definition_indices[load]];
        if (boundary.type == BoundaryConditionType::convection) {
            const spatial_detail::ConvectionValues values = convection_values(boundary);
            _boundary_data[load].load = values.coefficient;
            _boundary_data[load].ambient = values.ambient;
        } else {
            _boundary_data[load].load = spatial_detail::controlled_value(
                _definition, _time, _load_factor, boundary.value, boundary.scale_with_load, boundary.function);
        }
    }
}
LocalDofs SpatialAssembly::local_dofs(const std::array<std::size_t, 4>& nodes) const {
    LocalDofs result{};
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        result[node] = dof_map().temperature(nodes[node]);
        result[4 + node] = dof_map().radial_displacement(nodes[node]);
        result[8 + node] = dof_map().axial_displacement(nodes[node]);
    }
    return result;
}
bool SpatialAssembly::uses_augmented_contact() const noexcept {
    return std::any_of(_definition.contacts.begin(), _definition.contacts.end(), [](const ContactDefinition& value) {
        return value.mechanical && value.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian;
    });
}
AugmentedContactUpdate SpatialAssembly::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed_updates) {
    check_state_size(state.size(), dof_count(), "SpatialAssembly augmented-contact state size mismatch");
    AugmentedContactUpdate result;
    result.penetration_tolerance = std::numeric_limits<double>::infinity();
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        const ContactDefinition& definition = _definition.contacts[contact_value];
        if (!definition.mechanical ||
            definition.mechanical_formulation != MechanicalContactFormulation::augmented_lagrangian)
            continue;
        result.penetration_tolerance = std::min(result.penetration_tolerance, definition.penetration_tolerance);
        const std::vector<ContactNodeSummary> nodes = summarize_contact_nodes(contact_value, state);
        double contact_penetration = 0.0, contact_constraint_violation = 0.0;
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            contact_penetration = std::max(contact_penetration, std::max(-nodes[node].gap, 0.0));
            const bool captured = staged[contact_value][node].normal_multiplier > 0.0 || nodes[node].gap <= 0.0;
            if (captured)
                contact_constraint_violation = std::max(contact_constraint_violation, std::abs(nodes[node].gap));
        }
        result.maximum_penetration = std::max(result.maximum_penetration, contact_penetration);
        result.maximum_constraint_violation =
            std::max(result.maximum_constraint_violation, contact_constraint_violation);
        if (contact_constraint_violation <= definition.penetration_tolerance) continue;
        result.converged = false;
        if (completed_updates >= definition.maximum_augmented_iterations) {
            result.update_allowed = false;
            continue;
        }
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            ContactPointHistory& history = staged[contact_value][node];
            history.normal_multiplier = std::max(0.0, history.normal_multiplier - definition.penalty * nodes[node].gap);
        }
    }
    if (!std::isfinite(result.penetration_tolerance)) result.penetration_tolerance = 0.0;
    if (!result.converged && result.update_allowed) _contact_histories.swap(staged);
    return result;
}
void SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    check_state_size(state.size(), dof_count(), "SpatialAssembly committed contact state size mismatch");
    update_mechanical_candidates(state);
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    std::vector<std::vector<bool>> updated(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value)
        updated[contact_value].resize(_contact_histories[contact_value].size(), false);
    const std::size_t first_mechanical = contribution_ranges().mechanical_begin;
    for (std::size_t entry = 0; entry < _mechanical_contributions.size(); ++entry) {
        const MechanicalContribution& candidate = _mechanical_contributions[entry];
        if (!candidate.active) continue;
        const LocalValues local_state = contribution_state(first_mechanical + entry, state);
        const LocalValues committed_state = contribution_state(first_mechanical + entry, _committed_contact_solution);
        const ContactPointValue value =
            compute_node_to_line_rz_contact_value(_mechanical_properties[candidate.contact], candidate.geometry,
                local_state, committed_state, _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        const ContactPointHistory trial =
            compute_node_to_line_rz_contact_trial_history(_mechanical_properties[candidate.contact], candidate.geometry,
                local_state, committed_state, _contact_histories[candidate.contact][candidate.secondary]);
        if (updated[candidate.contact][candidate.secondary]) {
            const ContactPointHistory& prior = staged[candidate.contact][candidate.secondary];
            const double scale =
                std::max({1.0, std::abs(prior.elastic_tangential_slip), std::abs(trial.elastic_tangential_slip)});
            if (std::abs(prior.elastic_tangential_slip - trial.elastic_tangential_slip) >
                    32.0 * std::numeric_limits<double>::epsilon() * scale ||
                prior.sliding != trial.sliding)
                throw std::logic_error("Mechanical half-edges disagree on contact-node friction history");
            continue;
        }
        staged[candidate.contact][candidate.secondary] = trial;
        updated[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        if (!_definition.contacts[contact_value].mechanical) continue;
        if (std::find(updated[contact_value].begin(), updated[contact_value].end(), false) !=
            updated[contact_value].end())
            throw std::domain_error("Cannot commit friction history for an unprojected contact node");
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
}
void SpatialAssembly::restore_contact_state(
    const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories) {
    check_state_size(state.size(), dof_count(), "SpatialAssembly restored contact state layout mismatch");
    if (histories.size() != contact_count())
        throw std::invalid_argument("SpatialAssembly restored contact state layout mismatch");
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        if (histories[contact_value].size() != _contact_histories[contact_value].size())
            throw std::invalid_argument("SpatialAssembly restored contact history layout mismatch");
        for (const ContactPointHistory& history : histories[contact_value])
            if (!std::isfinite(history.elastic_tangential_slip) || !std::isfinite(history.normal_multiplier) ||
                history.normal_multiplier < 0.0)
                throw std::invalid_argument("SpatialAssembly restored contact history is invalid");
    }
    _committed_contact_solution = state;
    _contact_histories = std::move(histories);
}
void SpatialAssembly::update_thermal_candidates(const std::vector<double>& state) const {
    check_state_size(state.size(), dof_count(), "SpatialAssembly thermal-contact search state size mismatch");
    const ContributionRanges ranges = contribution_ranges();
    const GlobalStateView state_view(state);
    update_thermal_candidates(ranges.thermal_begin, ranges.mechanical_begin, state_view);
}
void SpatialAssembly::update_thermal_candidates(
    std::size_t first, std::size_t last, const GlobalStateView& state) const {
    check_state_size(state.global_size(), dof_count(), "SpatialAssembly thermal-contact shadow state size mismatch");
    for (const ThermalContribution& entry : _thermal_contributions) entry.active = false;
    std::fill(_projected_thermal_points.begin(), _projected_thermal_points.end(), 0U);
    if (!mark_touched_thermal_points(first, last)) return;
    const ContributionRanges ranges = contribution_ranges();
    std::fill(
        _thermal_minimum_distance.begin(), _thermal_minimum_distance.end(), std::numeric_limits<double>::infinity());
    std::fill(
        _thermal_selected_primary.begin(), _thermal_selected_primary.end(), std::numeric_limits<std::size_t>::max());
    std::fill(_projected_thermal_candidates.begin(), _projected_thermal_candidates.end(), 0U);
    for (std::size_t entry = 0; entry < _thermal_contributions.size(); ++entry) {
        const ThermalContribution& candidate = _thermal_contributions[entry];
        const std::size_t point = thermal_point_index(candidate.contact, candidate.integration_point);
        if (_touched_thermal_points[point] == 0U) continue;
        const LocalValues local_state = contribution_state(ranges.thermal_begin + entry, state);
        const HeatQuadratureValue value =
            compute_line2_rz_gap_heat_value(_thermal_properties[candidate.contact], candidate.geometry, local_state);
        if (!value.projected) continue;
        _projected_thermal_candidates[entry] = 1U;
        const double distance = std::abs(value.gap);
        if (distance < _thermal_minimum_distance[point] ||
            (distance == _thermal_minimum_distance[point] && candidate.primary < _thermal_selected_primary[point])) {
            _thermal_minimum_distance[point] = distance;
            _thermal_selected_primary[point] = candidate.primary;
        }
    }
    for (std::size_t entry = 0; entry < _thermal_contributions.size(); ++entry) {
        if (_projected_thermal_candidates[entry] == 0U) continue;
        const ThermalContribution& candidate = _thermal_contributions[entry];
        const std::size_t point = thermal_point_index(candidate.contact, candidate.integration_point);
        if (candidate.primary != _thermal_selected_primary[point]) continue;
        candidate.active = true;
        _projected_thermal_points[point] = 1U;
    }
}
void SpatialAssembly::update_mechanical_candidates(const std::vector<double>& state) const {
    check_state_size(state.size(), dof_count(), "SpatialAssembly contact-search state size mismatch");
    const ContributionRanges ranges = contribution_ranges();
    const GlobalStateView state_view(state);
    update_mechanical_candidates(ranges.mechanical_begin, ranges.boundary_begin, state_view);
}
void SpatialAssembly::update_mechanical_candidates(
    std::size_t first, std::size_t last, const GlobalStateView& state) const {
    check_state_size(state.global_size(), dof_count(), "SpatialAssembly contact-search shadow state size mismatch");
    for (const MechanicalContribution& entry : _mechanical_contributions) entry.active = false;
    std::fill(_projected_mechanical_nodes.begin(), _projected_mechanical_nodes.end(), 0U);
    if (!mark_touched_mechanical_nodes(first, last)) return;
    const ContributionRanges ranges = contribution_ranges();
    std::fill(_mechanical_minimum_distance.begin(), _mechanical_minimum_distance.end(),
        std::numeric_limits<double>::infinity());
    std::fill(_mechanical_selected_primary.begin(), _mechanical_selected_primary.end(),
        std::numeric_limits<std::size_t>::max());
    std::fill(_projected_mechanical_candidates.begin(), _projected_mechanical_candidates.end(), 0U);
    for (std::size_t entry = 0; entry < _mechanical_contributions.size(); ++entry) {
        const MechanicalContribution& candidate = _mechanical_contributions[entry];
        const std::size_t node = mechanical_node_index(candidate.contact, candidate.secondary);
        if (_touched_mechanical_nodes[node] == 0U) continue;
        const LocalValues local_state = contribution_state(ranges.mechanical_begin + entry, state);
        const LocalValues committed_state =
            contribution_state(ranges.mechanical_begin + entry, _committed_contact_solution);
        const ContactPointValue value =
            compute_node_to_line_rz_contact_value(_mechanical_properties[candidate.contact], candidate.geometry,
                local_state, committed_state, _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        _projected_mechanical_candidates[entry] = 1U;
        const double distance = std::abs(value.gap);
        if (distance < _mechanical_minimum_distance[node] ||
            (distance == _mechanical_minimum_distance[node] &&
                candidate.primary < _mechanical_selected_primary[node])) {
            _mechanical_minimum_distance[node] = distance;
            _mechanical_selected_primary[node] = candidate.primary;
        }
    }
    for (std::size_t entry = 0; entry < _mechanical_contributions.size(); ++entry) {
        if (_projected_mechanical_candidates[entry] == 0U) continue;
        const MechanicalContribution& candidate = _mechanical_contributions[entry];
        const std::size_t node = mechanical_node_index(candidate.contact, candidate.secondary);
        if (candidate.primary != _mechanical_selected_primary[node]) continue;
        candidate.active = true;
        _projected_mechanical_nodes[node] = 1U;
    }
}
bool SpatialAssembly::mark_touched_thermal_points(std::size_t first, std::size_t last) const {
    std::fill(_touched_thermal_points.begin(), _touched_thermal_points.end(), 0U);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin = std::max(first, ranges.thermal_begin), end = std::min(last, ranges.mechanical_begin);
    if (begin >= end) return false;
    for (std::size_t full = begin; full < end; ++full) {
        const ThermalContribution& candidate = _thermal_contributions[full - ranges.thermal_begin];
        _touched_thermal_points[thermal_point_index(candidate.contact, candidate.integration_point)] = 1U;
    }
    return true;
}
bool SpatialAssembly::mark_touched_mechanical_nodes(std::size_t first, std::size_t last) const {
    std::fill(_touched_mechanical_nodes.begin(), _touched_mechanical_nodes.end(), 0U);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin = std::max(first, ranges.mechanical_begin), end = std::min(last, ranges.boundary_begin);
    if (begin >= end) return false;
    for (std::size_t full = begin; full < end; ++full) {
        const MechanicalContribution& candidate = _mechanical_contributions[full - ranges.mechanical_begin];
        _touched_mechanical_nodes[mechanical_node_index(candidate.contact, candidate.secondary)] = 1U;
    }
    return true;
}
std::size_t SpatialAssembly::thermal_point_index(std::size_t contact, std::size_t point) const noexcept {
    return _thermal_contact_offsets[contact] + point;
}
std::size_t SpatialAssembly::mechanical_node_index(std::size_t contact, std::size_t node) const noexcept {
    return _mechanical_contact_offsets[contact] + node;
}
std::vector<ContactNodeSummary> SpatialAssembly::summarize_contact_nodes(
    std::size_t contact_value, const std::vector<double>& state) const {
    check_state_size(state.size(), dof_count(), "SpatialAssembly contact summary state size mismatch");
    update_mechanical_candidates(state);
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _meshes[secondary.region];
    std::vector<ContactNodeSummary> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) {
        result.push_back(
            {mesh.nodes().at(node).r, mesh.nodes().at(node).z, false, std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false});
    }
    const std::size_t first_mechanical = contribution_ranges().mechanical_begin;
    for (std::size_t entry = 0; entry < _mechanical_contributions.size(); ++entry) {
        const MechanicalContribution& candidate = _mechanical_contributions[entry];
        if (candidate.contact != contact_value) continue;
        if (!candidate.active) continue;
        const LocalValues local_state = contribution_state(first_mechanical + entry, state);
        const LocalValues committed_state = contribution_state(first_mechanical + entry, _committed_contact_solution);
        const std::size_t secondary_index = candidate.secondary;
        const ContactPointValue value = compute_node_to_line_rz_contact_value(_mechanical_properties[contact_value],
            candidate.geometry, local_state, committed_state, _contact_histories[contact_value][secondary_index]);
        if (!value.projected) continue;
        ContactNodeSummary& node = result.at(secondary_index);
        const std::size_t primary_segment = candidate.primary;
        if (node.projected && node.primary_segment != primary_segment)
            throw std::logic_error("Mechanical contact node has more than one active primary segment");
        node.projected = true;
        node.primary_segment = primary_segment;
        node.gap = std::min(node.gap, value.gap);
        node.tributary_area += value.tributary_area;
        node.tributary_length += value.tributary_length;
        node.contact_force += value.contact_force;
        node.tangential_force += value.tangential_force;
        node.elastic_tangential_slip = value.elastic_tangential_slip;
        node.sliding = value.sliding;
    }
    for (ContactNodeSummary& node : result) {
        if (node.tributary_area > 0.0) {
            node.pressure = node.contact_force / node.tributary_area;
            node.tangential_traction = node.tangential_force / node.tributary_area;
        }
    }
    return result;
}
void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    check_state_size(state.size(), dof_count(), "SpatialAssembly global state has the wrong size");
    update_thermal_candidates(state);
    update_mechanical_candidates(state);
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        if (_definition.contacts[contact_value].thermal) {
            std::size_t unprojected = 0;
            for (std::size_t point = _thermal_contact_offsets[contact_value];
                point < _thermal_contact_offsets[contact_value + 1]; ++point)
                if (_projected_thermal_points[point] == 0U) ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Thermal contact '" + _definition.contacts[contact_value].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        " secondary integration points after searching the complete primary chain");
        }
        if (!_definition.contacts[contact_value].mechanical) continue;
        std::size_t unprojected = 0;
        for (std::size_t node = _mechanical_contact_offsets[contact_value];
            node < _mechanical_contact_offsets[contact_value + 1]; ++node)
            if (_projected_mechanical_nodes[node] == 0U) ++unprojected;
        if (unprojected != 0)
            throw std::domain_error("Mechanical contact '" + _definition.contacts[contact_value].name +
                                    "' lost projection for " + std::to_string(unprojected) +
                                    " secondary nodes after searching the complete primary chain");
    }
}
void SpatialAssembly::validate_local_state(std::size_t first, std::size_t last, const GlobalStateView& state) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("SpatialAssembly contribution range is out of bounds");
    check_state_size(state.global_size(), dof_count(), "SpatialAssembly local state has the wrong global size");
    update_thermal_candidates(first, last, state);
    update_mechanical_candidates(first, last, state);
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        if (_definition.contacts[contact_value].thermal) {
            std::size_t unprojected = 0;
            for (std::size_t point = _thermal_contact_offsets[contact_value];
                point < _thermal_contact_offsets[contact_value + 1]; ++point)
                if (_touched_thermal_points[point] != 0U && _projected_thermal_points[point] == 0U) ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Thermal contact '" + _definition.contacts[contact_value].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        " locally owned secondary integration points after searching the complete "
                                        "primary chain");
        }
        if (_definition.contacts[contact_value].mechanical) {
            std::size_t unprojected = 0;
            for (std::size_t node = _mechanical_contact_offsets[contact_value];
                node < _mechanical_contact_offsets[contact_value + 1]; ++node)
                if (_touched_mechanical_nodes[node] != 0U && _projected_mechanical_nodes[node] == 0U) ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Mechanical contact '" + _definition.contacts[contact_value].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        " locally owned secondary nodes after searching the complete primary chain");
        }
    }
}
std::vector<std::size_t> SpatialAssembly::contact_secondary_source_nodes(std::size_t contact_value) const {
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<std::size_t> result;
    result.reserve(secondary.boundary.nodes.size());
    for (const std::size_t node : secondary.boundary.nodes) result.push_back(mesh.source_node_ids().at(node));
    return result;
}
InterfaceSummary SpatialAssembly::summarize_interface(
    std::size_t contact_value, const std::vector<double>& state) const {
    if (contact_value >= contact_count()) throw std::out_of_range("SpatialAssembly contact index is out of range");
    check_state_size(state.size(), dof_count(), "SpatialAssembly interface summary state size mismatch");
    update_thermal_candidates(state);
    if (_definition.contacts[contact_value].thermal) {
        std::size_t unprojected = 0;
        for (std::size_t point = _thermal_contact_offsets[contact_value];
            point < _thermal_contact_offsets[contact_value + 1]; ++point)
            if (_projected_thermal_points[point] == 0U) ++unprojected;
        if (unprojected != 0)
            throw std::domain_error("Thermal contact '" + _definition.contacts[contact_value].name +
                                    "' lost projection for " + std::to_string(unprojected) +
                                    " secondary integration points while summarizing the complete primary chain");
    }
    InterfaceSummary summary = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0,
        0,
    };
    const std::size_t first_thermal = contribution_ranges().thermal_begin;
    bool has_thermal = false;
    for (std::size_t entry = 0; entry < _thermal_contributions.size(); ++entry) {
        const ThermalContribution& candidate = _thermal_contributions[entry];
        if (candidate.contact != contact_value) continue;
        has_thermal = true;
        if (!candidate.active) continue;
        const LocalValues local_state = contribution_state(first_thermal + entry, state);
        const HeatQuadratureValue value =
            compute_line2_rz_gap_heat_value(_thermal_properties[contact_value], candidate.geometry, local_state);
        if (!value.projected) throw std::logic_error("Active thermal-contact candidate is not projected");
        summary.minimum_gap = std::min(summary.minimum_gap, value.gap);
        summary.total_heat_rate += value.weighted_measure * value.heat_flux;
    }
    if (!has_thermal) summary.minimum_gap = 0.0;
    bool has_mechanical = false;
    for (const MechanicalContribution& entry : _mechanical_contributions) {
        if (entry.contact == contact_value) {
            has_mechanical = true;
            break;
        }
    }
    if (has_mechanical) {
        for (const ContactNodeSummary& node : summarize_contact_nodes(contact_value, state)) {
            if (!node.projected) {
                ++summary.unprojected_contact_nodes;
                continue;
            }
            ++summary.projected_contact_nodes;
            summary.minimum_contact_gap = std::min(summary.minimum_contact_gap, node.gap);
            summary.maximum_contact_pressure = std::max(summary.maximum_contact_pressure, node.pressure);
            summary.total_contact_force += node.contact_force;
            summary.total_tangential_force += node.tangential_force;
            if (node.pressure > 0.0) ++summary.active_contact_nodes;
        }
    } else {
        summary.minimum_contact_gap = 0.0;
    }
    return summary;
}
namespace {
Line2InterfaceSideCoordinates edge_coordinates(const RegionMesh& mesh, const Line2BoundaryElement& edge) {
    return {{mesh.nodes().at(edge.nodes[0]), mesh.nodes().at(edge.nodes[1])}};
}
double planar_quad_area(const RegionMesh& mesh, const Quad4Element& element) {
    double twice_area = 0.0;
    for (std::size_t node = 0; node < element.nodes.size(); ++node) {
        const RzPoint &current = mesh.nodes().at(element.nodes[node]),
                      &next = mesh.nodes().at(element.nodes[(node + 1U) % element.nodes.size()]);
        twice_area += current.r * next.z - next.r * current.z;
    }
    return 0.5 * std::abs(twice_area);
}
double minimum_boundary_normal_length(
    const SpatialAssembly& layout, std::size_t region, const RegionBoundary& boundary) {
    const RegionMesh& mesh = layout.region_mesh(region);
    double result = std::numeric_limits<double>::infinity();
    for (const Line2BoundaryElement& edge : boundary.elements) {
        const auto parent = layout.edge_parent(region, edge);
        const Line2InterfaceSideCoordinates coordinates = edge_coordinates(mesh, edge);
        const double edge_length = std::hypot(coordinates[1].r - coordinates[0].r, coordinates[1].z - coordinates[0].z),
                     area = planar_quad_area(mesh, mesh.elements().at(parent.first)),
                     normal_length = area / edge_length;
        if (!std::isfinite(normal_length) || !(normal_length > 0.0))
            throw std::invalid_argument("Contact boundary has a nonpositive characteristic element length");
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
std::vector<RzPoint> boundary_parent_centroids(
    const SpatialAssembly& layout, std::size_t region, const RegionBoundary& boundary) {
    const RegionMesh& mesh = layout.region_mesh(region);
    std::vector<RzPoint> result;
    result.reserve(boundary.elements.size());
    for (const Line2BoundaryElement& edge : boundary.elements)
        result.push_back(element_centroid(mesh, mesh.elements().at(layout.edge_parent(region, edge).first)));
    return result;
}
std::array<std::size_t, 4> interface_nodes(const SpatialAssembly& layout, std::size_t secondary_region,
    const Line2BoundaryElement& secondary_edge, std::size_t primary_region, const Line2BoundaryElement& primary_edge) {
    return {layout.global_node(secondary_region, secondary_edge.nodes[0]),
        layout.global_node(secondary_region, secondary_edge.nodes[1]),
        layout.global_node(primary_region, primary_edge.nodes[0]),
        layout.global_node(primary_region, primary_edge.nodes[1])};
}
double primary_line_side(const RzPoint& point, const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double tangent_r = primary_coordinates[1].r - primary_coordinates[0].r,
                 tangent_z = primary_coordinates[1].z - primary_coordinates[0].z,
                 length = std::hypot(tangent_r, tangent_z);
    return ((point.r - primary_coordinates[0].r) * tangent_z - (point.z - primary_coordinates[0].z) * tangent_r) /
           length;
}
double zero_gap_orientation_hint(const RzPoint& secondary_centroid, const RzPoint& primary_centroid,
    const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double secondary_side = primary_line_side(secondary_centroid, primary_coordinates),
                 primary_side = primary_line_side(primary_centroid, primary_coordinates);
    if (secondary_side == 0.0 || primary_side == 0.0) return 0.0;
    const bool same_side = (secondary_side > 0.0 && primary_side > 0.0) || (secondary_side < 0.0 && primary_side < 0.0);
    return same_side ? 0.0 : secondary_side;
}
RegionBoundary ordered_connected_boundary(const RegionMesh& mesh, RegionBoundary boundary, const std::string& name) {
    if (boundary.elements.empty()) throw std::invalid_argument("Contact side set is empty: " + name);
    std::vector<std::size_t> degree(mesh.nodes().size(), 0);
    for (const Line2BoundaryElement& edge : boundary.elements) {
        if (edge.nodes[0] == edge.nodes[1]) throw std::invalid_argument("Contact side set has a zero edge: " + name);
        for (std::size_t node : edge.nodes)
            if (++degree.at(node) > 2) throw std::invalid_argument("Contact side set must not branch: " + name);
    }
    std::vector<std::size_t> endpoints;
    for (std::size_t node = 0; node < degree.size(); ++node)
        if (degree[node] == 1) endpoints.push_back(node);
    if (endpoints.size() != 2)
        throw std::invalid_argument("Contact side set must be one open connected chain: " + name);
    const auto coordinate_less = [&](std::size_t lhs, std::size_t rhs) {
        const RzPoint &left = mesh.nodes().at(lhs), &right = mesh.nodes().at(rhs);
        if (left.r != right.r) return left.r < right.r;
        if (left.z != right.z) return left.z < right.z;
        return lhs < rhs;
    };
    std::size_t current = coordinate_less(endpoints[0], endpoints[1]) ? endpoints[0] : endpoints[1];
    std::vector<bool> used(boundary.elements.size(), false);
    std::vector<Line2BoundaryElement> ordered;
    std::vector<std::size_t> nodes = {current};
    ordered.reserve(boundary.elements.size());
    nodes.reserve(boundary.elements.size() + 1);
    for (std::size_t position = 0; position < boundary.elements.size(); ++position) {
        std::size_t selected = boundary.elements.size();
        Line2BoundaryElement oriented{};
        for (std::size_t edge = 0; edge < boundary.elements.size(); ++edge) {
            if (used[edge]) continue;
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
            throw std::invalid_argument("Contact side set must be one connected chain: " + name);
        used[selected] = true;
        ordered.push_back(oriented);
        current = oriented.nodes[1];
        nodes.push_back(current);
    }
    boundary.nodes = std::move(nodes);
    boundary.elements = std::move(ordered);
    return boundary;
}
double projection_fraction(const RzPoint& point, const Line2InterfaceSideCoordinates& segment) {
    const double dr = segment[1].r - segment[0].r, dz = segment[1].z - segment[0].z;
    return ((point.r - segment[0].r) * dr + (point.z - segment[0].z) * dz) / (dr * dr + dz * dz);
}
bool projection_interval(const Line2InterfaceSideCoordinates& secondary, const Line2InterfaceSideCoordinates& primary,
    double& lower, double& upper) {
    const RzPoint midpoint = {
        0.5 * (secondary[0].r + secondary[1].r),
        0.5 * (secondary[0].z + secondary[1].z),
    };
    const RzPoint half = {
        0.5 * (secondary[1].r - secondary[0].r),
        0.5 * (secondary[1].z - secondary[0].z),
    };
    const double center = projection_fraction(midpoint, primary),
                 slope = projection_fraction({midpoint.r + half.r, midpoint.z + half.z}, primary) - center;
    lower = -1.0;
    upper = 1.0;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(slope) <= tolerance) return center >= -tolerance && center <= 1.0 + tolerance;
    double first = (0.0 - center) / slope, second = (1.0 - center) / slope;
    if (first > second) std::swap(first, second);
    lower = std::max(lower, first);
    upper = std::min(upper, second);
    lower = std::max(-1.0, std::min(1.0, lower));
    upper = std::max(-1.0, std::min(1.0, upper));
    return upper - lower > tolerance;
}
} // namespace
std::vector<RegionMesh> SpatialAssembly::build_meshes(
    const SpatialDefinition& definition, const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<RegionMesh> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions)
        if (region.block_id >= 0)
            result.push_back(RegionMesh::from_unstructured_block(source_mesh, region.block_id));
        else
            result.push_back(RegionMesh::from_unstructured_block(source_mesh, region.block));
    return result;
}
void SpatialAssembly::initialize_contact_search_workspace() {
    _thermal_contact_offsets.assign(contact_count() + 1, 0);
    _mechanical_contact_offsets.assign(contact_count() + 1, 0);
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        _thermal_contact_offsets[contact_value + 1] =
            _thermal_contact_offsets[contact_value] + _thermal_point_counts[contact_value];
        _mechanical_contact_offsets[contact_value + 1] =
            _mechanical_contact_offsets[contact_value] + _contact_histories[contact_value].size();
    }
    const std::size_t thermal_point_count = _thermal_contact_offsets.back();
    _touched_thermal_points.resize(thermal_point_count);
    _projected_thermal_points.resize(thermal_point_count);
    _thermal_minimum_distance.resize(thermal_point_count);
    _thermal_selected_primary.resize(thermal_point_count);
    _projected_thermal_candidates.resize(_thermal_contributions.size());
    const std::size_t mechanical_node_count = _mechanical_contact_offsets.back();
    _touched_mechanical_nodes.resize(mechanical_node_count);
    _projected_mechanical_nodes.resize(mechanical_node_count);
    _mechanical_minimum_distance.resize(mechanical_node_count);
    _mechanical_selected_primary.resize(mechanical_node_count);
    _projected_mechanical_candidates.resize(_mechanical_contributions.size());
}
SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : SpatialAssembly(definition, source_mesh, spatial_detail::resolve_block_ids(definition, source_mesh, true, true),
          build_meshes(definition, source_mesh)) {}
SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh,
    std::vector<std::int64_t> block_ids, std::vector<RegionMesh> meshes)
    : SpatialLayout(std::move(definition), std::move(block_ids), DofLayout::axisymmetric_rz),
      _meshes(std::move(meshes)) {
    std::vector<std::size_t> node_counts, element_counts;
    node_counts.reserve(_meshes.size());
    element_counts.reserve(_meshes.size());
    for (const RegionMesh& mesh : _meshes) {
        node_counts.push_back(mesh.nodes().size());
        element_counts.push_back(mesh.elements().size());
    }
    initialize_counts(node_counts, element_counts);
    _load_factor = 1.0;
    build_volume_geometries();
    build_contacts(source_mesh);
    build_boundaries(source_mesh);
    _contact_histories.resize(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value)
        _contact_histories[contact_value].resize(_secondary_boundaries[contact_value].boundary.nodes.size());
    initialize_contact_search_workspace();
    _committed_contact_solution = initial_state();
    update_thermal_candidates(_committed_contact_solution);
    update_mechanical_candidates(_committed_contact_solution);
}
std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("SpatialAssembly contribution range is out of bounds");
    std::vector<std::size_t> result;
    result.reserve((last - first) * local_dof_count);
    for (std::size_t entry = first; entry < last; ++entry) {
        const LocalDofs dofs = contribution_dofs(entry);
        for (const std::size_t dof : dofs) {
            if (dof >= dof_count()) throw std::out_of_range("SpatialAssembly contribution has an invalid DOF");
            result.push_back(dof);
        }
    }
    if (mark_touched_thermal_points(first, last)) {
        const ContributionRanges ranges = contribution_ranges();
        for (std::size_t entry = 0; entry < _thermal_contributions.size(); ++entry) {
            const ThermalContribution& candidate = _thermal_contributions[entry];
            if (_touched_thermal_points[thermal_point_index(candidate.contact, candidate.integration_point)] == 0U)
                continue;
            const LocalDofs dofs = contribution_dofs(ranges.thermal_begin + entry);
            result.insert(result.end(), dofs.begin(), dofs.end());
        }
    }
    if (mark_touched_mechanical_nodes(first, last)) {
        const ContributionRanges ranges = contribution_ranges();
        for (std::size_t entry = 0; entry < _mechanical_contributions.size(); ++entry) {
            const MechanicalContribution& candidate = _mechanical_contributions[entry];
            if (_touched_mechanical_nodes[mechanical_node_index(candidate.contact, candidate.secondary)] == 0U)
                continue;
            const LocalDofs dofs = contribution_dofs(ranges.mechanical_begin + entry);
            result.insert(result.end(), dofs.begin(), dofs.end());
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}
LocalDofs SpatialAssembly::contribution_dofs(std::size_t index) const {
    const ContributionLocation entry = locate_contribution(index);
    switch (entry.type) {
    case SpatialContributionType::volume: {
        const auto location = element_location(entry.local_index);
        const Quad4Element& element = _meshes[location.first].elements().at(location.second);
        std::array<std::size_t, 4> nodes{};
        for (std::size_t node = 0; node < nodes.size(); ++node)
            nodes[node] = global_node(location.first, element.nodes[node]);
        return local_dofs(nodes);
    }
    case SpatialContributionType::thermal_contact:
        return local_dofs(_thermal_contributions.at(entry.local_index).nodes);
    case SpatialContributionType::mechanical_contact:
        return local_dofs(_mechanical_contributions.at(entry.local_index).nodes);
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection: return local_dofs(_boundary_contributions.at(entry.local_index).nodes);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}
LocalValues SpatialAssembly::contribution_state(std::size_t index, const std::vector<double>& global_state) const {
    check_state_size(global_state.size(), dof_count(), "SpatialAssembly global state has the wrong size");
    return contribution_state(index, GlobalStateView(global_state));
}
LocalValues SpatialAssembly::contribution_state(std::size_t index, const GlobalStateView& global_state) const {
    check_state_size(global_state.global_size(), dof_count(), "SpatialAssembly shadow state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = global_state.value(dofs[local]);
    return result;
}
LocalResidual SpatialAssembly::contribution_residual(std::size_t index, const LocalValues& state) const {
    const ContributionLocation location = locate_contribution(index);
    switch (location.type) {
    case SpatialContributionType::volume:
        throw std::logic_error("SpatialAssembly does not own volume residual physics");
    case SpatialContributionType::thermal_contact: {
        const ThermalContribution& entry = _thermal_contributions[location.local_index];
        if (!entry.active) return {};
        return compute_line2_rz_gap_heat_residual(_thermal_properties[entry.contact], entry.geometry, state);
    }
    case SpatialContributionType::mechanical_contact: {
        const MechanicalContribution& entry = _mechanical_contributions.at(location.local_index);
        if (!entry.active) return {};
        return compute_node_to_line_rz_contact_residual(_mechanical_properties[entry.contact], entry.geometry, state,
            contribution_state(index, _committed_contact_solution), _contact_histories[entry.contact][entry.secondary]);
    }
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection: {
        const BoundaryContribution& entry = _boundary_contributions.at(location.local_index);
        return compute_line2_rz_boundary_residual(_boundary_data[entry.kernel], entry.geometry, state);
    }
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}
LocalSystem SpatialAssembly::linearize_contribution(std::size_t index, const LocalValues& state) const {
    const ContributionLocation location = locate_contribution(index);
    switch (location.type) {
    case SpatialContributionType::volume:
        throw std::logic_error("SpatialAssembly does not own volume Jacobian physics");
    case SpatialContributionType::thermal_contact: {
        const ThermalContribution& entry = _thermal_contributions[location.local_index];
        if (!entry.active) return {};
        return compute_line2_rz_gap_heat_system(_thermal_properties[entry.contact], entry.geometry, state);
    }
    case SpatialContributionType::mechanical_contact: {
        const MechanicalContribution& entry = _mechanical_contributions.at(location.local_index);
        if (!entry.active) return {};
        return compute_node_to_line_rz_contact_system(_mechanical_properties[entry.contact], entry.geometry, state,
            contribution_state(index, _committed_contact_solution), _contact_histories[entry.contact][entry.secondary]);
    }
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection: {
        const BoundaryContribution& entry = _boundary_contributions.at(location.local_index);
        return compute_line2_rz_boundary_system(_boundary_data[entry.kernel], entry.geometry, state);
    }
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}
void SpatialAssembly::build_contacts(const UnstructuredQuad4Mesh& source_mesh) {
    _primary_boundaries.reserve(contact_count());
    _secondary_boundaries.reserve(contact_count());
    _thermal_properties.reserve(contact_count());
    _mechanical_properties.reserve(contact_count());
    _thermal_point_counts.reserve(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count(); ++contact_value) {
        std::size_t thermal_point_count = 0;
        ContactDefinition& contact_definition = _definition.contacts[contact_value];
        ResolvedBoundary primary = resolve_boundary(source_mesh, contact_definition.primary);
        ResolvedBoundary secondary = resolve_boundary(source_mesh, contact_definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Self-contact is not supported: " + contact_definition.name);
        const RegionMesh &primary_mesh = _meshes[primary.region], &secondary_mesh = _meshes[secondary.region];
        primary.boundary =
            ordered_connected_boundary(primary_mesh, std::move(primary.boundary), contact_definition.primary);
        secondary.boundary =
            ordered_connected_boundary(secondary_mesh, std::move(secondary.boundary), contact_definition.secondary);
        const std::vector<RzPoint> primary_parent_centroids =
            boundary_parent_centroids(*this, primary.region, primary.boundary);
        const std::vector<RzPoint> secondary_parent_centroids =
            boundary_parent_centroids(*this, secondary.region, secondary.boundary);
        if (contact_definition.mechanical && contact_definition.automatic_penalty) {
            const double primary_length = minimum_boundary_normal_length(*this, primary.region, primary.boundary),
                         secondary_length = minimum_boundary_normal_length(*this, secondary.region, secondary.boundary),
                         primary_modulus = _definition.regions[primary.region].material.reference_young_modulus,
                         secondary_modulus = _definition.regions[secondary.region].material.reference_young_modulus;
            contact_definition.penalty = contact_definition.penalty_factor /
                                         (primary_length / primary_modulus + secondary_length / secondary_modulus);
            if (!std::isfinite(contact_definition.penalty) || !(contact_definition.penalty > 0.0))
                throw std::overflow_error(
                    "Automatic contact penalty is not finite and positive: " + contact_definition.name);
        }
        _thermal_properties.push_back({contact_definition.thermal ? contact_definition.gap_conductivity : 1.0,
            contact_definition.thermal ? contact_definition.minimum_gap : 1.0});
        _mechanical_properties.push_back({contact_definition.mechanical ? contact_definition.penalty : 1.0,
            contact_definition.mechanical ? contact_definition.friction_coefficient : 0.0,
            contact_definition.mechanical &&
                contact_definition.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian});
        if (contact_definition.thermal) {
            for (std::size_t edge_index = 0; edge_index < secondary.boundary.elements.size(); ++edge_index) {
                const Line2BoundaryElement& secondary_edge = secondary.boundary.elements[edge_index];
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                struct ThermalProjection final {
                    std::size_t primary_edge;
                    double lower, upper;
                };
                std::vector<ThermalProjection> projections;
                for (std::size_t primary_edge = 0; primary_edge < primary.boundary.elements.size(); ++primary_edge) {
                    double lower = 0.0, upper = 0.0;
                    if (projection_interval(secondary_coordinates,
                            edge_coordinates(primary_mesh, primary.boundary.elements[primary_edge]), lower, upper))
                        projections.push_back({primary_edge, lower, upper});
                }
                std::sort(projections.begin(), projections.end(),
                    [](const ThermalProjection& lhs, const ThermalProjection& rhs) { return lhs.lower < rhs.lower; });
                double covered = -1.0;
                constexpr double coverage_tolerance = 1.0e-10;
                for (const ThermalProjection& projection : projections) {
                    if (projection.lower > covered + coverage_tolerance ||
                        projection.lower < covered - coverage_tolerance)
                        throw std::invalid_argument("Secondary thermal edge projection has a gap or "
                                                    "overlap on the primary surface: " +
                                                    contact_definition.name);
                    covered = projection.upper;
                }
                if (covered < 1.0 - coverage_tolerance)
                    throw std::invalid_argument(
                        "Secondary thermal edge is outside the primary surface projection: " + contact_definition.name);
                for (const ThermalProjection& projection : projections) {
                    const Line2BoundaryElement& primary_edge = primary.boundary.elements[projection.primary_edge];
                    const Line2InterfaceSideCoordinates primary_coordinates =
                        edge_coordinates(primary_mesh, primary_edge);
                    const Line2RzHeatGeometry integration_geometry = make_line2_rz_heat_geometry(secondary_coordinates,
                        primary_coordinates, projection.lower, projection.upper,
                        zero_gap_orientation_hint(secondary_parent_centroids[edge_index],
                            primary_parent_centroids[projection.primary_edge], primary_coordinates));
                    for (const Line2RzHeatQuadraturePoint& point : integration_geometry.points) {
                        for (std::size_t candidate = 0; candidate < primary.boundary.elements.size(); ++candidate) {
                            const Line2BoundaryElement& candidate_edge = primary.boundary.elements[candidate];
                            const Line2InterfaceSideCoordinates candidate_coordinates =
                                edge_coordinates(primary_mesh, candidate_edge);
                            _thermal_contributions.push_back({contact_value,
                                interface_nodes(
                                    *this, secondary.region, secondary_edge, primary.region, candidate_edge),
                                make_line2_rz_heat_point_geometry(secondary_coordinates, candidate_coordinates,
                                    point.secondary_shape, point.integration_weight,
                                    candidate + 1 == primary.boundary.elements.size(),
                                    zero_gap_orientation_hint(secondary_parent_centroids[edge_index],
                                        primary_parent_centroids[candidate], candidate_coordinates)),
                                thermal_point_count, candidate});
                        }
                        ++thermal_point_count;
                    }
                }
            }
        }
        _thermal_point_counts.push_back(thermal_point_count);
        if (contact_definition.mechanical) {
            for (std::size_t edge_index = 0; edge_index < secondary.boundary.elements.size(); ++edge_index) {
                const Line2BoundaryElement& secondary_edge = secondary.boundary.elements[edge_index];
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                for (std::size_t secondary_node = 0; secondary_node < line2_interface_side_node_count;
                    ++secondary_node) {
                    for (std::size_t candidate = 0; candidate < primary.boundary.elements.size(); ++candidate) {
                        const Line2BoundaryElement& primary_edge = primary.boundary.elements[candidate];
                        const Line2InterfaceSideCoordinates primary_coordinates =
                            edge_coordinates(primary_mesh, primary_edge);
                        _mechanical_contributions.push_back({contact_value,
                            interface_nodes(*this, secondary.region, secondary_edge, primary.region, primary_edge),
                            make_node_to_line_rz_contact_geometry(secondary_coordinates, primary_coordinates,
                                secondary_node, candidate == 0, candidate + 1 == primary.boundary.elements.size(),
                                zero_gap_orientation_hint(secondary_parent_centroids[edge_index],
                                    primary_parent_centroids[candidate], primary_coordinates)),
                            edge_index + secondary_node, candidate});
                    }
                }
            }
        }
        _primary_boundaries.push_back(std::move(primary));
        _secondary_boundaries.push_back(std::move(secondary));
    }
}
} // namespace fuelsim::rz
