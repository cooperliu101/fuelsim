#include "spatial_assembly.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {

std::size_t SpatialAssembly::contact_count() const noexcept {
    return _definition.contacts.size();
}

const ContactDefinition& SpatialAssembly::contact(std::size_t index) const {
    return _definition.contacts.at(index);
}

const std::vector<std::vector<ContactPointHistory>>&
SpatialAssembly::committed_contact_histories() const noexcept {
    return _contact_histories;
}

bool SpatialAssembly::uses_augmented_contact() const noexcept {
    return std::any_of(
        _definition.contacts.begin(), _definition.contacts.end(),
        [](const ContactDefinition& contact) {
            return contact.mechanical &&
                   contact.mechanical_formulation ==
                       MechanicalContactFormulation::augmented_lagrangian;
        });
}

AugmentedContactUpdate SpatialAssembly::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed_updates) {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly augmented-contact state size mismatch");
    AugmentedContactUpdate result;
    result.penetration_tolerance = std::numeric_limits<double>::infinity();
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        const ContactDefinition& definition =
            _definition.contacts[contact_value];
        if (!definition.mechanical ||
            definition.mechanical_formulation !=
                MechanicalContactFormulation::augmented_lagrangian)
            continue;
        result.penetration_tolerance = std::min(
            result.penetration_tolerance, definition.penetration_tolerance);
        const std::vector<ContactNodeSummary> nodes =
            summarize_contact_nodes(contact_value, state);
        double contact_penetration = 0.0;
        double contact_constraint_violation = 0.0;
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            contact_penetration =
                std::max(contact_penetration, std::max(-nodes[node].gap, 0.0));
            const bool captured =
                staged[contact_value][node].normal_multiplier > 0.0 ||
                nodes[node].gap <= 0.0;
            if (captured)
                contact_constraint_violation = std::max(
                    contact_constraint_violation, std::abs(nodes[node].gap));
        }
        result.maximum_penetration =
            std::max(result.maximum_penetration, contact_penetration);
        result.maximum_constraint_violation = std::max(
            result.maximum_constraint_violation, contact_constraint_violation);
        if (contact_constraint_violation <= definition.penetration_tolerance)
            continue;
        result.converged = false;
        if (completed_updates >= definition.maximum_augmented_iterations) {
            result.update_allowed = false;
            continue;
        }
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            ContactPointHistory& history = staged[contact_value][node];
            history.normal_multiplier =
                std::max(0.0, history.normal_multiplier -
                                  definition.penalty * nodes[node].gap);
        }
    }
    if (!std::isfinite(result.penetration_tolerance))
        result.penetration_tolerance = 0.0;
    if (!result.converged && result.update_allowed)
        _contact_histories.swap(staged);
    return result;
}

void SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly committed contact state size mismatch");
    update_mechanical_candidates(state);
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    std::vector<std::vector<bool>> updated(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value)
        updated[contact_value].resize(_contact_histories[contact_value].size(),
                                      false);

    const std::size_t first_mechanical =
        contribution_ranges().mechanical_begin;
    for (std::size_t contribution = 0;
         contribution < _mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        if (!candidate.active)
            continue;
        const LocalValues local_state =
            contribution_state(first_mechanical + contribution, state);
        const LocalValues committed_state = contribution_state(
            first_mechanical + contribution, _committed_contact_solution);
        const ContactPointValue value =
            _mechanical_kernels[candidate.contact].value(
                candidate.geometry, local_state, committed_state,
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected)
            continue;
        const ContactPointHistory trial =
            _mechanical_kernels[candidate.contact].trial_history(
                candidate.geometry, local_state, committed_state,
                _contact_histories[candidate.contact][candidate.secondary]);
        if (updated[candidate.contact][candidate.secondary]) {
            const ContactPointHistory& prior =
                staged[candidate.contact][candidate.secondary];
            const double scale =
                std::max({1.0, std::abs(prior.elastic_tangential_slip),
                          std::abs(trial.elastic_tangential_slip)});
            if (std::abs(prior.elastic_tangential_slip -
                         trial.elastic_tangential_slip) >
                    32.0 * std::numeric_limits<double>::epsilon() * scale ||
                prior.sliding != trial.sliding)
                throw std::logic_error(
                    "Mechanical half-edge contributions disagree on the "
                    "contact-node friction history");
            continue;
        }
        staged[candidate.contact][candidate.secondary] = trial;
        updated[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        if (!_definition.contacts[contact_value].mechanical)
            continue;
        if (std::find(updated[contact_value].begin(),
                      updated[contact_value].end(),
                      false) != updated[contact_value].end())
            throw std::domain_error(
                "Cannot commit friction history for an unprojected contact "
                "node");
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
}

void SpatialAssembly::restore_contact_state(
    const std::vector<double>& state,
    std::vector<std::vector<ContactPointHistory>> histories) {
    if (state.size() != dof_count() || histories.size() != contact_count())
        throw std::invalid_argument(
            "SpatialAssembly restored contact state layout mismatch");
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        if (histories[contact_value].size() !=
            _contact_histories[contact_value].size())
            throw std::invalid_argument(
                "SpatialAssembly restored contact history layout mismatch");
        for (const ContactPointHistory& history : histories[contact_value]) {
            if (!std::isfinite(history.elastic_tangential_slip) ||
                !std::isfinite(history.normal_multiplier) ||
                history.normal_multiplier < 0.0)
                throw std::invalid_argument(
                    "SpatialAssembly restored contact history is invalid");
        }
    }
    _committed_contact_solution = state;
    _contact_histories = std::move(histories);
}


void SpatialAssembly::update_mechanical_candidates(
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly contact-search state size mismatch");
    const ContributionRanges ranges = contribution_ranges();
    const GlobalStateView state_view(state);
    update_mechanical_candidates(ranges.mechanical_begin,
                                 ranges.pressure_begin, state_view);
}

void SpatialAssembly::update_mechanical_candidates(
    std::size_t contribution_begin, std::size_t contribution_end,
    const GlobalStateView& state) const {
    if (state.global_size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly contact-search shadow state size mismatch");
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t local_begin =
        std::max(contribution_begin, ranges.mechanical_begin);
    const std::size_t local_end =
        std::min(contribution_end, ranges.pressure_begin);

    for (const MechanicalContribution& contribution : _mechanical_contributions)
        contribution.active = false;
    for (std::vector<bool>& nodes : _projected_mechanical_nodes)
        std::fill(nodes.begin(), nodes.end(), false);
    if (local_begin >= local_end)
        return;

    std::vector<std::vector<bool>> touched(contact_count());
    std::vector<std::vector<double>> minimum_distance(contact_count());
    std::vector<std::vector<std::size_t>> selected_primary(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        const std::size_t node_count = _contact_histories[contact_value].size();
        touched[contact_value].resize(node_count, false);
        minimum_distance[contact_value].assign(
            node_count, std::numeric_limits<double>::infinity());
        selected_primary[contact_value].assign(
            node_count, std::numeric_limits<std::size_t>::max());
    }
    for (std::size_t full = local_begin; full < local_end; ++full) {
        const std::size_t contribution = full - ranges.mechanical_begin;
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        touched[candidate.contact][candidate.secondary] = true;
    }

    std::vector<bool> projected(_mechanical_contributions.size(), false);
    for (std::size_t contribution = 0;
         contribution < _mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        if (!touched[candidate.contact][candidate.secondary])
            continue;
        const LocalValues local_state =
            contribution_state(ranges.mechanical_begin + contribution, state);
        const LocalValues committed_state = contribution_state(
            ranges.mechanical_begin + contribution,
            _committed_contact_solution);
        const ContactPointValue value =
            _mechanical_kernels[candidate.contact].value(
                candidate.geometry, local_state, committed_state,
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected)
            continue;
        projected[contribution] = true;
        const double distance = std::abs(value.gap);
        if (distance <
                minimum_distance[candidate.contact][candidate.secondary] ||
            (distance ==
                 minimum_distance[candidate.contact][candidate.secondary] &&
             candidate.primary <
                 selected_primary[candidate.contact][candidate.secondary])) {
            minimum_distance[candidate.contact][candidate.secondary] = distance;
            selected_primary[candidate.contact][candidate.secondary] =
                candidate.primary;
        }
    }

    for (std::size_t contribution = 0;
         contribution < _mechanical_contributions.size(); ++contribution) {
        if (!projected[contribution])
            continue;
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        if (candidate.primary !=
            selected_primary[candidate.contact][candidate.secondary])
            continue;
        candidate.active = true;
        _projected_mechanical_nodes[candidate.contact][candidate.secondary] =
            true;
    }
}

std::vector<ContactNodeSummary>
SpatialAssembly::summarize_contact_nodes(std::size_t contact_value,
                                       const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly contact summary state size mismatch");
    update_mechanical_candidates(state);
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _meshes[secondary.region];
    std::vector<ContactNodeSummary> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) {
        result.push_back({mesh.nodes().at(node).r, mesh.nodes().at(node).z,
                          false, std::numeric_limits<std::size_t>::max(),
                          std::numeric_limits<double>::infinity(), 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, false});
    }
    const std::size_t first_mechanical =
        contribution_ranges().mechanical_begin;
    for (std::size_t contribution = 0;
         contribution < _mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        if (candidate.contact != contact_value)
            continue;
        if (!candidate.active)
            continue;
        const LocalValues local_state =
            contribution_state(first_mechanical + contribution, state);
        const LocalValues committed_state = contribution_state(
            first_mechanical + contribution, _committed_contact_solution);
        const std::size_t secondary_index = candidate.secondary;
        const ContactPointValue value =
            _mechanical_kernels[contact_value].value(
                candidate.geometry, local_state, committed_state,
                _contact_histories[contact_value][secondary_index]);
        if (!value.projected)
            continue;
        ContactNodeSummary& node = result.at(secondary_index);
        const std::size_t primary_segment = candidate.primary;
        if (node.projected && node.primary_segment != primary_segment)
            throw std::logic_error(
                "Mechanical contact node has more than one active primary "
                "segment");
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
        if (node.tributary_area > 0.0)
            node.pressure = node.contact_force / node.tributary_area;
        if (node.tributary_area > 0.0)
            node.tangential_traction =
                node.tangential_force / node.tributary_area;
    }
    return result;
}

void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly global state has the wrong size");
    update_mechanical_candidates(state);
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        if (!_definition.contacts[contact_value].mechanical)
            continue;
        const std::size_t unprojected = static_cast<std::size_t>(std::count(
            _projected_mechanical_nodes[contact_value].begin(),
            _projected_mechanical_nodes[contact_value].end(), false));
        if (unprojected != 0)
            throw std::domain_error(
                "Mechanical contact '" +
                _definition.contacts[contact_value].name +
                "' lost projection " + "for " + std::to_string(unprojected) +
                " secondary nodes after searching the complete primary chain");
    }
}

void SpatialAssembly::validate_local_state(std::size_t contribution_begin,
                                           std::size_t contribution_end,
                                           const GlobalStateView& state) const {
    if (contribution_begin > contribution_end ||
        contribution_end > contribution_count())
        throw std::out_of_range(
            "SpatialAssembly contribution range is out of bounds");
    if (state.global_size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly local state has the wrong global size");
    update_mechanical_candidates(contribution_begin, contribution_end, state);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t local_begin =
        std::max(contribution_begin, ranges.mechanical_begin);
    const std::size_t local_end =
        std::min(contribution_end, ranges.pressure_begin);
    std::vector<std::vector<bool>> touched(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value)
        touched[contact_value].resize(_contact_histories[contact_value].size(),
                                      false);
    for (std::size_t full = local_begin; full < local_end; ++full) {
        const std::size_t contribution = full - ranges.mechanical_begin;
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        touched[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        if (!_definition.contacts[contact_value].mechanical)
            continue;
        std::size_t unprojected = 0;
        for (std::size_t node = 0; node < touched[contact_value].size();
             ++node) {
            if (touched[contact_value][node] &&
                !_projected_mechanical_nodes[contact_value][node])
                ++unprojected;
        }
        if (unprojected != 0)
            throw std::domain_error(
                "Mechanical contact '" +
                _definition.contacts[contact_value].name +
                "' lost projection for " + std::to_string(unprojected) +
                " locally owned secondary nodes after searching the complete "
                "primary chain");
    }
}

std::vector<std::size_t>
SpatialAssembly::contact_secondary_source_nodes(std::size_t contact_value) const {
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<std::size_t> result;
    result.reserve(secondary.boundary.nodes.size());
    for (const std::size_t node : secondary.boundary.nodes)
        result.push_back(mesh.source_node_ids().at(node));
    return result;
}

InterfaceSummary
SpatialAssembly::summarize_interface(std::size_t contact_value,
                                   const std::vector<double>& state) const {
    if (contact_value >= contact_count())
        throw std::out_of_range("SpatialAssembly contact index is out of range");
    InterfaceSummary summary = {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0,
        0,
        0.0,
    };
    const std::size_t first_thermal = contribution_ranges().thermal_begin;
    bool has_thermal = false;
    for (std::size_t contribution = 0;
         contribution < _thermal_contributions.size(); ++contribution) {
        const ThermalContribution& candidate =
            _thermal_contributions[contribution];
        if (candidate.contact != contact_value)
            continue;
        has_thermal = true;
        const LocalValues local_state =
            contribution_state(first_thermal + contribution, state);
        const HeatQuadratureValues values =
            _thermal_kernels[contact_value].quadrature_values(
                candidate.geometry, local_state);
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
    for (const MechanicalContribution& contribution :
         _mechanical_contributions) {
        if (contribution.contact == contact_value) {
            has_mechanical = true;
            break;
        }
    }
    if (has_mechanical) {
        for (const ContactNodeSummary& node :
             summarize_contact_nodes(contact_value, state)) {
            if (!node.projected) {
                ++summary.unprojected_contact_nodes;
                continue;
            }
            ++summary.projected_contact_nodes;
            summary.minimum_contact_gap =
                std::min(summary.minimum_contact_gap, node.gap);
            summary.maximum_contact_pressure =
                std::max(summary.maximum_contact_pressure, node.pressure);
            summary.total_contact_force += node.contact_force;
            summary.total_tangential_force += node.tangential_force;
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
