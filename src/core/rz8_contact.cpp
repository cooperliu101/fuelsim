#include "rz8_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::rz8 {
namespace {
constexpr auto invalid = std::numeric_limits<std::size_t>::max();

void order_chain(Quad8RegionBoundary& boundary) {
    auto remaining = boundary.elements;
    if (remaining.empty()) throw std::invalid_argument("CAX8T contact boundary is empty");
    std::size_t first = invalid;
    for (std::size_t i = 0; i < remaining.size(); ++i) {
        bool predecessor = false;
        for (const auto& e : remaining)
            if (e.nodes[1] == remaining[i].nodes[0]) predecessor = true;
        if (!predecessor) {
            if (first != invalid) throw std::invalid_argument("CAX8T contact must form one connected open chain");
            first = i;
        }
    }
    if (first == invalid) throw std::invalid_argument("CAX8T contact must form an open chain");
    boundary.elements.clear();
    while (!remaining.empty()) {
        const auto next = remaining[first];
        boundary.elements.push_back(next);
        remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(first));
        if (remaining.empty()) break;
        first = invalid;
        for (std::size_t i = 0; i < remaining.size(); ++i)
            if (remaining[i].nodes[0] == next.nodes[1]) {
                if (first != invalid) throw std::invalid_argument("CAX8T contact chain branches");
                first = i;
            }
        if (first == invalid)
            throw std::invalid_argument("CAX8T contact chain is disconnected or inconsistently oriented");
    }
}
} // namespace

void SpatialAssembly::build_contacts(const UnstructuredQuad8Mesh& source) {
    _contact_histories.resize(_definition.contacts.size());
    for (std::size_t c = 0; c < _definition.contacts.size(); ++c) {
        const auto& contact = _definition.contacts[c];
        if (contact.mechanical_discretization == MechanicalContactDiscretization::surface_to_surface)
            throw std::invalid_argument("Axisymmetric CAX8T mechanical contact uses node-to-surface discretization");
        for (const auto* name : {&contact.primary, &contact.secondary}) {
            const auto block = source.side_set_block_id(*name);
            const auto found = std::find(_block_ids.begin(), _block_ids.end(), block);
            if (found == _block_ids.end())
                throw std::invalid_argument("CAX8T contact boundary is outside selected regions");
            const auto r = static_cast<std::size_t>(found - _block_ids.begin());
            auto boundary = _meshes[r].map_side_set(source, *name);
            order_chain(boundary);
            (name == &contact.primary ? _primary : _secondary).push_back({r, std::move(boundary)});
        }
        const auto &primary = _primary.back(), &secondary = _secondary.back();
        _contact_histories[c].resize(secondary.boundary.displacement_nodes.size());
        double penalty = contact.penalty;
        if (contact.automatic_penalty) {
            double compliance = 0;
            for (const auto* boundary : {&primary, &secondary}) {
                const auto& material = region(boundary->region).material;
                const auto properties = IsotropicThermoelasticMaterial(material).active_properties(
                    region(boundary->region).initial_temperature);
                const double modulus = (properties.lame_lambda + 2 * properties.shear_modulus).value();
                double size = 0;
                for (const auto& edge : boundary->boundary.elements) {
                    const auto& geometry = _geometries[boundary->region][edge.parent_element];
                    double volume = 0, radius = 0;
                    for (std::size_t q = 0; q < geometry.point_count; ++q) {
                        const auto& p = geometry.points[q];
                        volume += p.weighted_measure;
                        radius += p.weighted_measure * p.radius;
                    }
                    radius /= volume;
                    const auto& a = _meshes[boundary->region].nodes()[edge.nodes[0]];
                    const auto& b = _meshes[boundary->region].nodes()[edge.nodes[1]];
                    size += volume / (2 * std::acos(-1.0) * radius * std::hypot(b.r - a.r, b.z - a.z));
                }
                size /= static_cast<double>(boundary->boundary.elements.size());
                compliance += size / modulus;
            }
            penalty = contact.penalty_factor / compliance;
        }
        if (contact.mechanical && (!(penalty > 0) || !std::isfinite(penalty)))
            throw std::invalid_argument("CAX8T contact penalty must be positive");
        double maximum_elastic_slip = 0.0;
        if (contact.friction_slip_tolerance > 0.0) {
            double reference_length = 0.0;
            for (const auto& edge : secondary.boundary.elements) {
                std::array<RzPoint, 3> coordinates;
                for (std::size_t n = 0; n < coordinates.size(); ++n)
                    coordinates[n] = _meshes[secondary.region].nodes()[edge.nodes[n]];
                const double g = std::sqrt(3.0 / 5.0);
                const std::array<double, 3> locations = {-g, 0.0, g};
                const std::array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
                for (std::size_t q = 0; q < locations.size(); ++q) {
                    const double x = locations[q];
                    const std::array<double, 3> derivative = {x - 0.5, x + 0.5, -2.0 * x};
                    double dr = 0.0, dz = 0.0;
                    for (std::size_t n = 0; n < coordinates.size(); ++n) {
                        dr += derivative[n] * coordinates[n].r;
                        dz += derivative[n] * coordinates[n].z;
                    }
                    reference_length += weights[q] * std::hypot(dr, dz);
                }
            }
            maximum_elastic_slip = contact.friction_slip_tolerance * reference_length /
                                   static_cast<double>(secondary.boundary.elements.size());
            if (!std::isfinite(maximum_elastic_slip) || !(maximum_elastic_slip > 0.0))
                throw std::invalid_argument("CAX8T slip_tolerance gives an invalid elastic slip: " + contact.name);
        }
        _mechanical.push_back({penalty, contact.friction_coefficient,
            contact.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian,
            maximum_elastic_slip});
        _heat.push_back({contact.gap_conductivity, contact.minimum_gap, contact.gap_heat_conductance_law,
            contact.gap_conductance, contact.gap_conductance_clearance_derivative,
            contact.gap_conductance_pressure_derivative, contact.gap_conductance_temperature_derivative,
            contact.gap_conductance_reference_temperature, penalty});
        for (const auto& edge : secondary.boundary.elements) {
            std::array<RzPoint, 3> coordinates;
            for (std::size_t n = 0; n < 3; ++n) coordinates[n] = _meshes[secondary.region].nodes()[edge.nodes[n]];
            std::vector<double> cuts = {-1, 1};
            for (const auto& pe : primary.boundary.elements)
                for (std::size_t n = 0; n < 2; ++n) {
                    const auto projection = project_line3(coordinates, _meshes[primary.region].nodes()[pe.nodes[n]]);
                    if (projection.first && projection.second > -1 + 1e-12 && projection.second < 1 - 1e-12)
                        cuts.push_back(projection.second);
                }
            std::sort(cuts.begin(), cuts.end());
            cuts.erase(
                std::unique(cuts.begin(), cuts.end(), [](double a, double b) { return std::abs(a - b) < 1e-12; }),
                cuts.end());
            const auto add_point = [&](bool mechanical, double coordinate, double weight, std::size_t node) {
                const auto point = _contact_points.size(), first = _candidates.size();
                const auto sn =
                    static_cast<std::size_t>(std::find(secondary.boundary.displacement_nodes.begin(),
                                                 secondary.boundary.displacement_nodes.end(), edge.nodes[node]) -
                                             secondary.boundary.displacement_nodes.begin());
                for (std::size_t p = 0; p < primary.boundary.elements.size(); ++p) {
                    const auto& pe = primary.boundary.elements[p];
                    Candidate candidate{{coordinates, {}, p == 0, p + 1 == primary.boundary.elements.size(), coordinate,
                                            weight, node, mechanical},
                        {}, c, point, sn, p};
                    for (std::size_t n = 0; n < 3; ++n)
                        candidate.geometry.primary[n] = _meshes[primary.region].nodes()[pe.nodes[n]];
                    for (std::size_t n = 0; n < 2; ++n) {
                        candidate.dofs[n] =
                            dof(Field::temperature, global_temperature_node(secondary.region, edge.nodes[n]));
                        candidate.dofs[2 + n] =
                            dof(Field::temperature, global_temperature_node(primary.region, pe.nodes[n]));
                    }
                    for (std::size_t n = 0; n < 3; ++n) {
                        candidate.dofs[4 + n] =
                            dof(Field::radial_displacement, global_node(secondary.region, edge.nodes[n]));
                        candidate.dofs[7 + n] =
                            dof(Field::axial_displacement, global_node(secondary.region, edge.nodes[n]));
                        candidate.dofs[10 + n] =
                            dof(Field::radial_displacement, global_node(primary.region, pe.nodes[n]));
                        candidate.dofs[13 + n] =
                            dof(Field::axial_displacement, global_node(primary.region, pe.nodes[n]));
                    }
                    _candidates.push_back(candidate);
                }
                _contact_points.push_back({first, _candidates.size()});
            };
            if (contact.mechanical)
                for (std::size_t n = 0; n < 3; ++n) add_point(true, 0, 0, n);
            if (contact.thermal) {
                const double g = std::sqrt(3.0 / 5.0);
                const std::array<double, 3> qs = {-g, 0, g}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9};
                for (std::size_t i = 1; i < cuts.size(); ++i)
                    for (std::size_t q = 0; q < 3; ++q)
                        add_point(false, (cuts[i] + cuts[i - 1]) / 2 + qs[q] * (cuts[i] - cuts[i - 1]) / 2,
                            weights[q] * (cuts[i] - cuts[i - 1]) / 2, 0);
            }
        }
    }
    _active_candidates.assign(_contact_points.size(), invalid);
}

Line3ContactResult SpatialAssembly::candidate_value(
    std::size_t index, const std::vector<double>& state, bool jacobian) const {
    const auto& candidate = _candidates.at(index);
    std::vector<double> local(16), old(16);
    for (std::size_t i = 0; i < 16; ++i) {
        local[i] = state[candidate.dofs[i]];
        old[i] = _committed_contact_solution[candidate.dofs[i]];
    }
    return compute_line3_contact(candidate.geometry, _heat[candidate.contact], _mechanical[candidate.contact], local,
        old, _contact_histories[candidate.contact][candidate.secondary], jacobian);
}

void SpatialAssembly::validate_contact_state(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    const auto offset = volume_contribution_count() + _boundaries.size();
    for (std::size_t p = 0; p < _contact_points.size(); ++p) {
        const auto& point = _contact_points[p];
        if (last <= offset + point.first || first >= offset + point.last) continue;
        auto& selected = _active_candidates[p];
        selected = invalid;
        double distance = std::numeric_limits<double>::infinity();
        for (auto c = point.first; c < point.last; ++c) {
            const auto value = candidate_value(c, state);
            const bool mechanical = _candidates[c].geometry.mechanical;
            const bool projected = mechanical ? value.mechanical.projected : value.thermal.projected;
            const double gap = mechanical ? value.mechanical.gap : value.thermal.gap;
            if (projected && std::abs(gap) < distance) {
                selected = c;
                distance = std::abs(gap);
            }
        }
        if (selected == invalid) throw std::domain_error("CAX8T contact point lost all current primary projections");
    }
}

std::vector<std::size_t> SpatialAssembly::contact_secondary_source_nodes(std::size_t contact) const {
    const auto& boundary = _secondary.at(contact);
    std::vector<std::size_t> result;
    for (auto n : boundary.boundary.displacement_nodes) result.push_back(_meshes[boundary.region].source_node_ids()[n]);
    return result;
}

std::vector<ContactNodeSummary> SpatialAssembly::summarize_contact_nodes(
    std::size_t contact, const std::vector<double>& state) const {
    validate_contact_state(0, contribution_count(), state);
    const auto& boundary = _secondary.at(contact);
    std::vector<ContactNodeSummary> result(boundary.boundary.displacement_nodes.size());
    for (std::size_t n = 0; n < result.size(); ++n) {
        const auto& p = _meshes[boundary.region].nodes()[boundary.boundary.displacement_nodes[n]];
        result[n].r = p.r;
        result[n].z = p.z;
    }
    for (auto c : _active_candidates) {
        const auto& candidate = _candidates[c];
        if (candidate.contact != contact || !candidate.geometry.mechanical) continue;
        const auto value = candidate_value(c, state).mechanical;
        auto& row = result[candidate.secondary];
        row.projected = value.projected;
        row.primary_segment = candidate.primary;
        row.gap = value.gap;
        row.pressure = value.pressure;
        row.tributary_area += value.tributary_area;
        row.tributary_length += value.tributary_length;
        row.contact_force += value.contact_force;
        row.tangential_force += value.tangential_force;
        row.tangential_traction = value.tangential_traction;
        row.elastic_tangential_slip = value.elastic_tangential_slip;
        row.sliding = value.sliding;
        row.total_tangential_slip = value.total_tangential_slip;
    }
    return result;
}

std::vector<std::array<double, 2>> SpatialAssembly::recover_contact_tractions(
    std::size_t contact, const std::vector<ContactNodeSummary>& nodes, const std::vector<double>& state) const {
    const auto& boundary = _secondary.at(contact).boundary;
    if (nodes.size() != boundary.displacement_nodes.size())
        throw std::invalid_argument("CAX8T contact recovery node count mismatch");
    std::vector<std::array<double, 2>> result(nodes.size());
    if (!_definition.contacts.at(contact).mechanical) return result;
    std::vector<std::size_t> counts(nodes.size());
    std::vector<std::size_t> indices(_meshes[_secondary[contact].region].nodes().size(), invalid);
    std::array<double, 2> minimum = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    std::array<double, 2> maximum = {
        -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (std::size_t n = 0; n < nodes.size(); ++n) {
        indices[boundary.displacement_nodes[n]] = n;
        if (!std::isfinite(nodes[n].pressure) || nodes[n].pressure < 0 || !std::isfinite(nodes[n].tangential_traction))
            throw std::domain_error("CAX8T contact recovery requires finite compressive tractions");
        maximum[0] = std::max(maximum[0], nodes[n].pressure);
        maximum[1] = std::max(maximum[1], nodes[n].tangential_traction);
        minimum[0] = std::min(minimum[0], nodes[n].pressure);
        minimum[1] = std::min(minimum[1], nodes[n].tangential_traction);
    }
    // The pair-wide range includes zero at primary nodes that participate in no
    // active constraint. B10.14 distinguishes this from secondary-only bounds.
    const auto& primary = _primary.at(contact);
    std::vector<bool> participating(_meshes[primary.region].nodes().size(), false);
    for (const auto c : _active_candidates) {
        const auto& candidate = _candidates[c];
        if (candidate.contact != contact || !candidate.geometry.mechanical || nodes[candidate.secondary].pressure <= 0)
            continue;
        const double x = candidate_value(c, state).primary_coordinate;
        const std::array<double, 3> shape = {x * (x - 1) / 2, x * (x + 1) / 2, 1 - x * x};
        const auto& edge = primary.boundary.elements[candidate.primary];
        for (std::size_t n = 0; n < 3; ++n)
            if (shape[n] != 0.0) participating[edge.nodes[n]] = true;
    }
    for (const auto node : primary.boundary.displacement_nodes)
        if (!participating[node])
            for (std::size_t field = 0; field < 2; ++field) {
                minimum[field] = std::min(minimum[field], 0.0);
                maximum[field] = std::max(maximum[field], 0.0);
            }
    // Equal-weight least-squares projection from quadratic nodal values onto
    // {1, xi}, with xi = {-1, 1, 0}. Shared edge endpoints are averaged.
    constexpr std::array<std::array<double, 3>, 3> projection = {
        {{5.0 / 6, -1.0 / 6, 1.0 / 3}, {-1.0 / 6, 5.0 / 6, 1.0 / 3}, {1.0 / 3, 1.0 / 3, 1.0 / 3}}};
    for (const auto& edge : boundary.elements) {
        std::array<std::size_t, 3> output;
        for (std::size_t i = 0; i < 3; ++i) {
            output[i] = indices.at(edge.nodes[i]);
            if (output[i] == invalid) throw std::logic_error("CAX8T recovery edge node missing from boundary");
        }
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                result[output[i]][0] += projection[i][j] * nodes[output[j]].pressure;
                result[output[i]][1] += projection[i][j] * nodes[output[j]].tangential_traction;
            }
            ++counts[output[i]];
        }
    }
    for (std::size_t n = 0; n < nodes.size(); ++n) {
        if (counts[n] == 0) throw std::logic_error("CAX8T recovery node has no adjacent edge");
        // This limiter is part of output recovery, not material/geometry clipping.
        // It prevents new extrema in the pair-wide traction range.
        for (std::size_t field = 0; field < 2; ++field)
            result[n][field] =
                std::clamp(result[n][field] / static_cast<double>(counts[n]), minimum[field], maximum[field]);
        result[n][1] = -result[n][1];
    }
    return result;
}

InterfaceSummary SpatialAssembly::summarize_interface(std::size_t contact, const std::vector<double>& state) const {
    InterfaceSummary summary;
    const auto nodes = summarize_contact_nodes(contact, state);
    for (const auto& n : nodes) {
        if (n.projected) {
            ++summary.projected_contact_nodes;
            summary.minimum_contact_gap = std::min(summary.minimum_contact_gap, n.gap);
        } else
            ++summary.unprojected_contact_nodes;
        if (n.pressure > 0) ++summary.active_contact_nodes;
        if (n.projected) summary.minimum_gap = std::min(summary.minimum_gap, n.gap);
        summary.maximum_contact_pressure = std::max(summary.maximum_contact_pressure, n.pressure);
        summary.total_contact_force += n.contact_force;
        summary.total_tangential_force += n.tangential_force;
    }
    for (auto c : _active_candidates) {
        const auto& candidate = _candidates[c];
        if (candidate.contact != contact || candidate.geometry.mechanical) continue;
        const auto value = candidate_value(c, state).thermal;
        summary.minimum_gap = std::min(summary.minimum_gap, value.gap);
        summary.total_heat_rate += value.heat_flux * value.weighted_measure;
    }
    return summary;
}

void SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    validate_state(state);
    auto staged = _contact_histories;
    for (std::size_t c = 0; c < _definition.contacts.size(); ++c) {
        if (!_definition.contacts[c].mechanical) continue;
        const auto rows = summarize_contact_nodes(c, state);
        for (std::size_t n = 0; n < rows.size(); ++n) {
            staged[c][n].elastic_tangential_slip = rows[n].elastic_tangential_slip;
            staged[c][n].total_tangential_slip = rows[n].total_tangential_slip;
            staged[c][n].sliding = rows[n].sliding;
        }
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
}

void SpatialAssembly::restore_contact_state(
    const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories) {
    if (state.size() != dof_count() || histories.size() != _contact_histories.size())
        throw std::invalid_argument("CAX8T contact history layout mismatch");
    for (std::size_t c = 0; c < histories.size(); ++c) {
        if (histories[c].size() != _contact_histories[c].size())
            throw std::invalid_argument("CAX8T contact node history layout mismatch");
        for (const auto& h : histories[c])
            if (!std::isfinite(h.elastic_tangential_slip) || !std::isfinite(h.total_tangential_slip) ||
                !std::isfinite(h.normal_multiplier) || h.normal_multiplier < 0)
                throw std::invalid_argument("CAX8T contact history is invalid");
    }
    _contact_histories = std::move(histories);
    _committed_contact_solution = state;
}

bool SpatialAssembly::uses_augmented_contact() const noexcept {
    for (const auto& c : _definition.contacts)
        if (c.mechanical && c.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian) return true;
    return false;
}

AugmentedContactUpdate SpatialAssembly::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed) {
    AugmentedContactUpdate result;
    for (std::size_t c = 0; c < _definition.contacts.size(); ++c) {
        const auto& contact = _definition.contacts[c];
        if (!contact.mechanical || !_mechanical[c].augmented_lagrangian) continue;
        const auto rows = summarize_contact_nodes(c, state);
        for (std::size_t n = 0; n < rows.size(); ++n) {
            const double penetration = std::max(0., -rows[n].gap),
                         violation =
                             _contact_histories[c][n].normal_multiplier > 0 ? std::abs(rows[n].gap) : penetration;
            result.maximum_penetration = std::max(result.maximum_penetration, penetration);
            result.maximum_constraint_violation = std::max(result.maximum_constraint_violation, violation);
            result.penetration_tolerance = std::max(result.penetration_tolerance, contact.penetration_tolerance);
            if (violation > contact.penetration_tolerance) result.converged = false;
            if (completed >= contact.maximum_augmented_iterations)
                result.update_allowed = false;
            else
                _contact_histories[c][n].normal_multiplier =
                    std::max(0., _contact_histories[c][n].normal_multiplier - _mechanical[c].penalty * rows[n].gap);
        }
    }
    return result;
}
} // namespace fuelsim::rz8
