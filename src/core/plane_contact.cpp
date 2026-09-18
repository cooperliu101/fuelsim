#include "plane_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <map>
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
        for (const auto& secondary_side : secondary)
            for (const auto& primary_side : primary)
                if (secondary_side.region == primary_side.region && secondary_side.element == primary_side.element)
                    throw std::invalid_argument("CPEG8T contact sides cannot belong to the same element");
        for (const bool thermal : {false, true}) {
            if (thermal ? !contact.thermal : !contact.mechanical)
                continue;
            std::map<std::size_t, std::vector<std::pair<ContactSide, std::size_t>>> supports;
            for (const auto& side : secondary) {
                const auto& dofs = _dofs.at(region_element_offset(side.region) + side.element);
                for (std::size_t n = 0; n < (thermal ? 2U : 3U); ++n) {
                    const auto local = n < 2 ? (side.side + n) % 4 : 4 + side.side;
                    supports[dofs[4 + local]].push_back({side, n});
                }
            }
            for (const auto& support : supports) {
                ContactConstraint constraint{};
                constraint.contact = c;
                constraint.geometry.penalty = penalty;
                constraint.geometry.heat = {contact.gap_conductivity,
                    contact.minimum_gap,
                    contact.gap_heat_conductance_law,
                    contact.gap_conductance,
                    contact.gap_conductance_clearance_derivative,
                    contact.gap_conductance_pressure_derivative,
                    contact.gap_conductance_temperature_derivative,
                    contact.gap_conductance_reference_temperature,
                    penalty};
                std::map<std::size_t, std::size_t> local_nodes, temperature_dofs;
                const auto append_edge = [&](const ContactSide& side) {
                    std::array<std::size_t, 3> nodes{};
                    const auto& dofs = _dofs.at(region_element_offset(side.region) + side.element);
                    const auto& g = geometry(side.region, side.element);
                    for (std::size_t n = 0; n < 3; ++n) {
                        const auto local = n < 2 ? (side.side + n) % 4 : 4 + side.side;
                        const auto inserted = local_nodes.emplace(dofs[4 + local], local_nodes.size());
                        nodes[n] = inserted.first->second;
                        if (thermal && n < 2)
                            temperature_dofs.emplace(nodes[n], dofs[local]);
                        if (inserted.second) {
                            constraint.geometry.coordinates.push_back(g.coordinates[local]);
                            constraint.dofs.push_back(dofs[4 + local]);
                            constraint.dofs.push_back(dofs[12 + local]);
                        }
                    }
                    return nodes;
                };
                for (const auto& entry : support.second)
                    constraint.geometry.secondary.push_back({append_edge(entry.first),
                        entry.second,
                        geometry(entry.first.region, entry.first.element).thickness});
                for (const auto& side : primary)
                    constraint.geometry.primary.push_back(append_edge(side));
                for (const auto& entry : temperature_dofs) {
                    constraint.geometry.temperature_nodes.push_back(entry.first);
                    constraint.dofs.push_back(entry.second);
                }
                _contact_constraints.push_back(std::move(constraint));
            }
        }
    }
}

void SpatialAssembly::validate_contacts(std::size_t first, std::size_t last, const std::vector<double>& state) const {
    for (auto index = std::max(first, contact_offset()); index < last; ++index) {
        const auto& constraint = _contact_constraints.at(index - contact_offset());
        std::vector<double> local;
        for (auto dof : constraint.dofs)
            local.push_back(state.at(dof));
        (void)elements::evaluate_plane_averaged_contact(constraint.geometry, local, false);
    }
}

elements::PlaneSurfaceContactResult
SpatialAssembly::compute_contact(std::size_t index, const std::vector<double>& local, bool jacobian) const {
    return elements::evaluate_plane_averaged_contact(_contact_constraints.at(index - contact_offset()).geometry,
        local,
        jacobian);
}

std::vector<elements::Line3PlaneContactResult> SpatialAssembly::contact_points(std::size_t contact,
    const std::vector<double>& state) const {
    if (contact >= _definition.contacts.size())
        throw std::out_of_range("CPEG8T contact index");
    validate_contacts(contact_offset(), contribution_count(), state);
    std::vector<elements::Line3PlaneContactResult> result;
    for (const auto& constraint : _contact_constraints) {
        if (constraint.contact != contact)
            continue;
        std::vector<double> local;
        for (auto dof : constraint.dofs)
            local.push_back(state.at(dof));
        result.push_back(elements::evaluate_plane_averaged_contact(constraint.geometry, local, false).point);
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
