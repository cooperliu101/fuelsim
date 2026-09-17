#include "plane_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim::plane {
void SpatialAssembly::build_contacts(const UnstructuredPlaneQuad8Mesh& mesh) {
    const auto sides = [&](const std::string& name) {
        std::vector<ContactSide> result;
        for (const auto& side : mesh.side_set(name).sides) {
            bool found = false;
            for (std::size_t r = 0; r < region_count(); ++r) {
                const auto it = std::find(_source_elements[r].begin(), _source_elements[r].end(), side.element);
                if (it != _source_elements[r].end()) {
                    result.push_back({r, static_cast<std::size_t>(it - _source_elements[r].begin()), side.local_side});
                    found = true;
                }
            }
            if (!found)
                throw std::invalid_argument("CPEG8T contact side lies outside selected regions");
        }
        if (result.empty())
            throw std::invalid_argument("CPEG8T contact side set must not be empty");
        return result;
    };
    const auto compliance = [&](const std::vector<ContactSide>& surface) {
        double total = 0.0;
        for (const auto& side : surface) {
            const auto& g = geometry(side.region, side.element);
            const auto& a = g.coordinates[side.side];
            const auto& b = g.coordinates[(side.side + 1) % 4];
            const double length = std::hypot(b[0] - a[0], b[1] - a[1]);
            double volume = 0.0;
            for (const auto& point : g.points)
                volume += point.measure;
            total += volume / (g.thickness * length * region(side.region).material.reference_young_modulus);
        }
        return total / static_cast<double>(surface.size());
    };
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> locations{-g, 0.0, g}, weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (std::size_t c = 0; c < _definition.contacts.size(); ++c) {
        const auto& contact = _definition.contacts[c];
        if (contact.mechanical_formulation != MechanicalContactFormulation::penalty
            || contact.friction_coefficient != 0.0
            || (contact.mechanical && contact.mechanical_sliding != MechanicalContactSliding::finite)
            || contact.mechanical_discretization == MechanicalContactDiscretization::node_to_surface
            || contact.thermal_discretization != ThermalContactDiscretization::surface_to_surface)
            throw std::invalid_argument(
                "CPEG8T contact requires frictionless finite-sliding surface-to-surface penalty and heat transfer");
        if (contact.primary == contact.secondary)
            throw std::invalid_argument("CPEG8T primary and secondary surfaces must differ");
        const auto primary = sides(contact.primary), secondary = sides(contact.secondary);
        const double penalty = contact.automatic_penalty
                                   ? contact.penalty_factor / (compliance(primary) + compliance(secondary))
                                   : contact.penalty;
        if (contact.mechanical && (!std::isfinite(penalty) || !(penalty > 0)))
            throw std::invalid_argument("CPEG8T mechanical contact requires positive penalty");
        for (const auto& secondary_side : secondary) {
            const auto first = _candidates.size();
            for (const auto& primary_side : primary)
                for (std::size_t segment = 0; segment < 5; ++segment)
                    for (std::size_t q = 0; q < 3; ++q) {
                        if (secondary_side.region == primary_side.region
                            && secondary_side.element == primary_side.element)
                            throw std::invalid_argument("CPEG8T contact sides cannot belong to the same element");
                        Candidate
                            candidate{secondary_side, primary_side, {}, c, locations[q], weights[q], penalty, segment};
                        for (std::size_t side_index = 0; side_index < 2; ++side_index) {
                            const auto& side = side_index == 0 ? secondary_side : primary_side;
                            const auto& volume_dofs = _dofs.at(region_element_offset(side.region) + side.element);
                            const std::array<std::size_t, 3> nodes{side.side, (side.side + 1) % 4, 4 + side.side};
                            for (std::size_t n = 0; n < 2; ++n)
                                candidate.dofs[2 * side_index + n] = volume_dofs[nodes[n]];
                            for (std::size_t n = 0; n < 3; ++n) {
                                candidate.dofs[4 + 6 * side_index + n] = volume_dofs[4 + nodes[n]];
                                candidate.dofs[7 + 6 * side_index + n] = volume_dofs[12 + nodes[n]];
                                candidate.dofs[16 + 3 * side_index + n] = volume_dofs[20 + n];
                            }
                        }
                        _candidates.push_back(candidate);
                    }
            _constraints.emplace_back(first, _candidates.size());
        }
    }
}

elements::Line3PlaneContactResult
SpatialAssembly::evaluate_candidate(std::size_t index, const elements::Line3PlaneValues& state, bool jacobian) const {
    const auto& candidate = _candidates.at(index);
    const auto& secondary = geometry(candidate.secondary.region, candidate.secondary.element);
    const auto& primary = geometry(candidate.primary.region, candidate.primary.element);
    const auto& contact = _definition.contacts[candidate.contact];
    elements::Line3PlaneContactInput input{{},
        {},
        secondary.reference_point,
        primary.reference_point,
        secondary.thickness,
        primary.thickness,
        region(candidate.secondary.region).strain_formulation == StrainFormulation::finite,
        region(candidate.primary.region).strain_formulation == StrainFormulation::finite,
        state,
        candidate.coordinate,
        candidate.weight,
        {contact.gap_conductivity,
            contact.minimum_gap,
            contact.gap_heat_conductance_law,
            contact.gap_conductance,
            contact.gap_conductance_clearance_derivative,
            contact.gap_conductance_pressure_derivative,
            contact.gap_conductance_temperature_derivative,
            contact.gap_conductance_reference_temperature,
            candidate.penalty},
        candidate.penalty,
        contact.thermal,
        contact.mechanical,
        candidate.segment};
    for (std::size_t n = 0; n < 3; ++n) {
        const auto secondary_node = n < 2 ? (candidate.secondary.side + n) % 4 : 4 + candidate.secondary.side;
        const auto primary_node = n < 2 ? (candidate.primary.side + n) % 4 : 4 + candidate.primary.side;
        input.secondary[n] = secondary.coordinates[secondary_node];
        input.primary[n] = primary.coordinates[primary_node];
    }
    return elements::evaluate_line3_plane_contact(input, jacobian);
}

void SpatialAssembly::validate_contacts(std::size_t first, std::size_t last, const std::vector<double>& state) const {
    for (std::size_t constraint = 0; constraint < _constraints.size(); ++constraint) {
        const auto [begin, end] = _constraints[constraint];
        if (last <= contact_offset() + begin || first >= contact_offset() + end)
            continue;
        std::vector<std::pair<double, double>> intervals;
        for (std::size_t candidate = begin; candidate < end; candidate += 3) {
            elements::Line3PlaneValues local{};
            for (std::size_t i = 0; i < local.size(); ++i)
                local[i] = state.at(_candidates[candidate].dofs[i]);
            const auto result = evaluate_candidate(candidate, local, false);
            for (std::size_t q = 1; q < 3; ++q)
                if (evaluate_candidate(candidate + q, local, false).projected != result.projected)
                    throw std::domain_error("CPEG8T contact projection changes inside a segment");
            if (result.projected)
                intervals.emplace_back(result.interval_begin, result.interval_end);
        }
        std::sort(intervals.begin(), intervals.end());
        if (intervals.empty())
            throw std::domain_error("CPEG8T secondary contact edge has no valid primary projection");
        // Clip the ends to the actual overlap. Curvature can move an endpoint
        // projection just outside the opposing edge even without gross sliding.
        double covered = intervals.front().first;
        for (const auto& interval : intervals) {
            if (interval.first > covered + 1e-10)
                throw std::domain_error("CPEG8T secondary contact edge has an uncovered projection interval");
            if (interval.first < covered - 1e-10)
                throw std::domain_error("CPEG8T primary projection intervals overlap ambiguously");
            covered = interval.second;
        }
    }
}

elements::Line3PlaneContactResult
SpatialAssembly::compute_contact(std::size_t index, const std::vector<double>& local, bool jacobian) const {
    const auto candidate = index - contact_offset();
    if (local.size() != 22)
        throw std::invalid_argument("CPEG8T contact requires 22 interface-local values");
    elements::Line3PlaneValues values{};
    std::copy(local.begin(), local.end(), values.begin());
    return evaluate_candidate(candidate, values, jacobian);
}

std::vector<elements::Line3PlaneContactResult> SpatialAssembly::contact_points(std::size_t contact,
    const std::vector<double>& state) const {
    if (contact >= _definition.contacts.size())
        throw std::out_of_range("CPEG8T contact index");
    validate_contacts(contact_offset(), contribution_count(), state);
    std::vector<elements::Line3PlaneContactResult> result;
    for (std::size_t selected = 0; selected < _candidates.size(); ++selected) {
        const auto& candidate = _candidates.at(selected);
        if (candidate.contact != contact)
            continue;
        elements::Line3PlaneValues local{};
        for (std::size_t i = 0; i < local.size(); ++i)
            local[i] = state.at(candidate.dofs[i]);
        result.push_back(evaluate_candidate(selected, local, false));
    }
    return result;
}

InterfaceSummary SpatialAssembly::summarize_interface(std::size_t contact, const std::vector<double>& state) const {
    InterfaceSummary result;
    for (const auto& point : contact_points(contact, state)) {
        if (!point.projected)
            continue;
        result.minimum_gap = std::min(result.minimum_gap, point.gap);
        result.minimum_contact_gap = std::min(result.minimum_contact_gap, point.gap);
        result.maximum_contact_pressure = std::max(result.maximum_contact_pressure, point.pressure);
        result.total_heat_rate += point.heat_rate;
        result.total_contact_force += point.force;
        ++result.projected_contact_nodes;
        if (point.pressure > 0)
            ++result.active_contact_nodes;
    }
    return result;
}
} // namespace fuelsim::plane
