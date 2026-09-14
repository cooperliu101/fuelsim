#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Real = long double;
using fuelsim::test::ExodusResults;
constexpr Real pi = 3.141592653589793238462643383279502884L;
constexpr Real young = 1.0e9L, poisson = 0.3L, mu = 0.3L;
constexpr Real lame = young * poisson / ((1.0L + poisson) * (1.0L - 2.0L * poisson));
constexpr Real shear = young / (2.0L * (1.0L + poisson));
constexpr Real axial_modulus = lame + 2.0L * shear;
constexpr std::array<Real, 2> heights{0.01L, 0.02L};
constexpr Real total_height = heights[0] + heights[1];
constexpr std::array<std::array<Real, 2>, 2> radii{{{{0.004L, 0.005L}}, {{0.00501L, 0.006L}}}};
constexpr std::array<std::array<std::size_t, 3>, 2> controls{{{{8, 9, 10}}, {{11, 12, 13}}}};
constexpr std::array<std::array<std::size_t, 2>, 4> radial_nodes{{{{0, 1}}, {{4, 5}}, {{2, 3}}, {{6, 7}}}};
constexpr Real elastic_slip_limit = 1.0e-6L * total_height / 2.0L;

void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

class Audit final {
  public:
    void add(const std::string& name, double actual, Real reference, double zero_tolerance) {
        auto& metric = _metrics[name];
        metric.values.add(actual, static_cast<double>(reference));
        metric.zero_tolerance = std::max(metric.zero_tolerance, zero_tolerance);
    }

    bool passed() const {
        bool result = true;
        for (const auto& entry : _metrics) {
            const auto& values = entry.second.values;
            if (values.has_relative_norm())
                fuelsim::test::print_relative_metrics(entry.first, values);
            else
                fuelsim::test::print_absolute_metrics(entry.first, values);
            result = (!values.has_relative_norm() || fuelsim::test::relative_metrics_below(values, 0.001))
                     && values.maximum_zero_reference_difference <= entry.second.zero_tolerance && result;
        }
        return result;
    }

  private:
    struct Metric final {
        fuelsim::test::FieldErrorMetrics values;
        double zero_tolerance = 0.0;
    };

    std::map<std::string, Metric> _metrics;
};

void check_summary(const std::string& path) {
    std::ifstream stream(path);
    require(static_cast<bool>(stream), "Cannot read chain production summary");
    std::string line;
    std::getline(stream, line);
    require(line == "metric,value", "Unexpected chain production summary header");
    std::map<std::string, std::string> values;
    while (std::getline(stream, line)) {
        const auto comma = line.find(',');
        require(comma != std::string::npos, "Incomplete chain production summary row");
        require(values.emplace(line.substr(0, comma), line.substr(comma + 1)).second,
            "Duplicate chain production summary metric");
    }
    require(values.at("completed") == "true" || values.at("completed") == "1", "Chain solve did not complete");
    require(std::stod(values.at("accepted_steps")) == 10.0 && std::stod(values.at("rejected_steps")) == 0.0
                && std::abs(std::stod(values.at("committed_time")) - 1.0) < 1.0e-12,
        "Chain validation requires all ten fixed increments without rejected steps");
}

struct BodyReference final {
    std::array<Real, 3> displacement{};
    std::array<Real, 2> height = heights;
    std::array<Real, 2> axial_elastic{};
    std::array<Real, 2> radial_elastic{};
    Real radial_stretch = 1.0L, area = 0.0L;
    std::array<std::array<Real, 4>, 2> stress{};
};

struct Reference final {
    std::array<BodyReference, 2> bodies{};
    std::array<std::array<Real, 2>, 2> total_slip{};
    std::array<std::array<Real, 2>, 2> elastic_slip{};
    std::array<Real, 2> pressure{}, gap{};
    std::array<bool, 2> active{};
};

// Given old layer lengths a,b, new total length L, and d equal to the
// prescribed increment of (sigma_z1-sigma_z2)/(2*Czz), diagonal Hughes-Winget
// integration gives
//   d*(h+a)*(L-h+b) = 2*((a+b)*h-a*L).
// This quadratic has exactly one root in (0,L). No equilibrium iteration or
// production material/element/solver code is used by this reference.
Real first_layer_height(Real a, Real b, Real length, Real d) {
    if (d == 0.0L)
        return a * length / (a + b);
    const Real coefficient_b = 2.0L * (a + b) - d * (length + b - a);
    const Real coefficient_c = -2.0L * a * length - d * a * (length + b);
    const Real discriminant = coefficient_b * coefficient_b - 4.0L * d * coefficient_c;
    require(discriminant > 0.0L, "Independent chain quadratic has no distinct real roots");
    const Real q = -0.5L * (coefficient_b + std::copysign(std::sqrt(discriminant), coefficient_b));
    const std::array<Real, 2> roots{q / d, coefficient_c / q};
    Real selected = 0.0L;
    std::size_t admissible = 0;
    for (const Real root : roots)
        if (root > 0.0L && root < length) {
            selected = root;
            ++admissible;
        }
    require(admissible == 1, "Independent chain quadratic must have one positive physical layer-length root");
    return selected;
}

Real prescribed_radial_strain(bool finite, std::size_t layer) {
    return finite || layer == 0 ? 0.004L : 0.0038L;
}

Reference next_reference(const Reference& old, std::size_t step, bool finite, bool touching_active) {
    Reference result;
    const Real time = static_cast<Real>(step) / 10.0L;
    for (std::size_t layer = 0; layer < 2; ++layer) {
        result.gap[layer] = 1.0e-5L - radii[0][1] * prescribed_radial_strain(finite, layer) * time;
        if (step == 5 && (finite || layer == 0))
            result.gap[layer] = 0.0L;
        result.pressure[layer] = std::max(0.0L, -1.0e10L * result.gap[layer]);
        result.active[layer] = result.gap[layer] < 0.0L || (result.gap[layer] == 0.0L && touching_active);
    }
    const Real top = (finite ? 0.003L : 0.0001L) * time;
    const Real current_inner_radius = radii[0][1] * (1.0L + 0.004L * time);
    Real middle_friction = 0.0L;
    for (std::size_t layer = 0; layer < 2; ++layer)
        middle_friction += 2.0L * pi * (finite ? current_inner_radius : radii[0][1]) * mu * result.pressure[layer]
                           * heights[layer] / 2.0L;
    if (finite)
        middle_friction *= (total_height + top) / total_height;
    for (std::size_t body = 0; body < result.bodies.size(); ++body) {
        auto& current = result.bodies[body];
        const auto& previous = old.bodies[body];
        const Real sign = body == 0 ? -1.0L : 1.0L;
        current.displacement[2] = body == 0 ? top : 0.0L;
        current.radial_stretch = body == 0 ? 1.0L + 0.004L * time : 1.0L;
        const Real reference_area = pi * (radii[body][1] * radii[body][1] - radii[body][0] * radii[body][0]);
        current.area = reference_area * (finite ? current.radial_stretch * current.radial_stretch : 1.0L);
        if (!finite) {
            for (std::size_t layer = 0; layer < 2; ++layer)
                current.radial_elastic[layer] = body == 0 ? prescribed_radial_strain(finite, layer) * time : 0.0L;
            const Real radial_stress_difference = 2.0L * lame * (current.radial_elastic[0] - current.radial_elastic[1]);
            current.displacement[1] = current.displacement[2] * heights[0] / total_height
                                      + (sign * middle_friction / reference_area - radial_stress_difference)
                                            / (axial_modulus * (1.0L / heights[0] + 1.0L / heights[1]));
        } else {
            for (std::size_t layer = 0; layer < 2; ++layer)
                current.radial_elastic[layer] = previous.radial_elastic[layer]
                                                + 2.0L * (current.radial_stretch - previous.radial_stretch)
                                                      / (current.radial_stretch + previous.radial_stretch);
            if (step <= 5) {
                current.displacement[1] = current.displacement[2] * heights[0] / total_height;
            } else {
                const Real target_difference = sign * middle_friction / current.area;
                const Real old_difference = previous.stress[0][1] - previous.stress[1][1];
                const Real d = (target_difference - old_difference) / (2.0L * axial_modulus);
                current.displacement[1] = first_layer_height(previous.height[0],
                                              previous.height[1],
                                              total_height + current.displacement[2],
                                              d)
                                          - heights[0];
            }
        }
        for (std::size_t layer = 0; layer < 2; ++layer) {
            current.height[layer] = heights[layer] + current.displacement[layer + 1] - current.displacement[layer];
            if (!finite)
                current.axial_elastic[layer] =
                    (current.displacement[layer + 1] - current.displacement[layer]) / heights[layer];
            else if (step <= 5) {
                const Real stretch = 1.0L + current.displacement[2] / total_height;
                const Real old_stretch = 1.0L + previous.displacement[2] / total_height;
                current.axial_elastic[layer] =
                    previous.axial_elastic[layer] + 2.0L * (stretch - old_stretch) / (stretch + old_stretch);
            } else
                current.axial_elastic[layer] = previous.axial_elastic[layer]
                                               + 2.0L * (current.height[layer] - previous.height[layer])
                                                     / (current.height[layer] + previous.height[layer]);
            const Real trace = 2.0L * current.radial_elastic[layer] + current.axial_elastic[layer];
            current.stress[layer] = {{lame * trace + 2.0L * shear * current.radial_elastic[layer],
                lame * trace + 2.0L * shear * current.axial_elastic[layer],
                lame * trace + 2.0L * shear * current.radial_elastic[layer],
                0.0L}};
        }
    }
    constexpr Real gauss = 0.577350269189625764509148780501957456L;
    for (std::size_t layer = 0; layer < 2; ++layer)
        for (std::size_t q = 0; q < 2; ++q) {
            const Real upper = 0.5L * (1.0L + (q == 0 ? -gauss : gauss));
            Real increment = 0.0L;
            for (std::size_t end = 0; end < 2; ++end) {
                const Real shape = end == 0 ? 1.0L - upper : upper;
                increment +=
                    shape
                    * ((result.bodies[0].displacement[layer + end] - old.bodies[0].displacement[layer + end])
                        - (result.bodies[1].displacement[layer + end] - old.bodies[1].displacement[layer + end]));
            }
            result.total_slip[layer][q] = old.total_slip[layer][q];
            if (result.active[layer]) {
                require(increment > elastic_slip_limit, "Independent reference must be on the forward sliding branch");
                result.total_slip[layer][q] += increment;
                result.elastic_slip[layer][q] = elastic_slip_limit;
            }
        }
    return result;
}

void check_geometry(const ExodusResults& frame) {
    require(frame.nodes.size() == 14 && frame.block_names == std::vector<std::string>{"inner", "outer"}
                && frame.block_element_counts == std::vector<std::size_t>{2, 2},
        "Chain acceptance requires the complete two-body, two-layer BAR2 mesh");
    for (std::size_t body = 0; body < 2; ++body) {
        for (std::size_t station = 0; station < 3; ++station) {
            const auto& node = frame.nodes[controls[body][station]];
            const Real z = station == 0 ? 0.0L : (station == 1 ? heights[0] : total_height);
            require(node[0] == 0.0 && std::abs(static_cast<Real>(node[1]) - z) < 1.0e-16L,
                "Axial control coordinates changed");
        }
        for (std::size_t layer = 0; layer < 2; ++layer)
            for (std::size_t end = 0; end < 2; ++end) {
                const auto& node = frame.nodes[radial_nodes[2 * body + layer][end]];
                const Real z = layer == 0 ? heights[0] / 2.0L : heights[0] + heights[1] / 2.0L;
                require(std::abs(static_cast<Real>(node[0]) - radii[body][end]) < 1.0e-16L
                            && std::abs(static_cast<Real>(node[1]) - z) < 1.0e-16L,
                    "Radial mesh radii or layer centers changed");
            }
    }
}

bool check_chain(bool finite, const std::string& output, const std::string& summary) {
    check_summary(summary);
    const auto frames = fuelsim::test::read_exodus_history(output);
    require(frames.size() == 11 && frames.front().time == 0.0,
        "Chain acceptance requires the initial frame and all ten increments");
    check_geometry(frames.front());
    // The exactly touching frame has a zero physical traction. Its cumulative
    // slip activation can nevertheless depend on the sign of the independent
    // double-precision geometry subtraction; record that boundary explicitly.
    const double touching_gap = frames.front().nodes[2][0] - frames.front().nodes[1][0] - 20.0e-6 * frames[5].time;
    const bool touching_active = touching_gap < 0.0;
    std::cout << "touching_frame_roundoff_contact_active=" << (touching_active ? 1 : 0) << '\n';
    Audit audit;
    Reference previous;
    Real previous_elastic_energy = 0.0L;
    std::size_t coupled_frames = 0;
    constexpr std::array<const char*, 4> components{"rr", "zz", "hoop", "rz"};
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        const Real time = static_cast<Real>(step) / 10.0L;
        require(std::abs(static_cast<Real>(frame.time) - time) < 1.0e-12L, "Chain increment time changed");
        const auto reference = next_reference(previous, step, finite, touching_active);
        const auto& temperature = frame.nodal("temperature");
        const auto& ur = frame.nodal("displacement_r");
        const auto& uz = frame.nodal("displacement_z");
        const auto& role = frame.nodal("node_role");
        const auto& reaction_r = frame.nodal("reaction_force_r");
        const auto& reaction_z = frame.nodal("reaction_force_z");
        const auto& reaction_heat = frame.nodal("reaction_heat_flux");
        if (!finite) {
            audit.add("cross_slice_inner_radial_difference", ur.at(1) - ur.at(5), 1.0e-6L * time, 1.0e-14);
            require(ur.at(1) > ur.at(5), "Inner slices must develop different radial displacements");
        }
        std::array<Real, 2> contact_area{}, friction_force{};
        Real total_normal = 0.0L, total_friction = 0.0L, body_work = 0.0L, contact_work = 0.0L, dirichlet_work = 0.0L,
             friction_dissipation = 0.0L, elastic_energy = 0.0L;
        for (std::size_t layer = 0; layer < 2; ++layer) {
            const Real radius = radii[0][1] * (finite ? reference.bodies[0].radial_stretch : 1.0L);
            contact_area[layer] = 2.0L * pi * radius * (finite ? reference.bodies[0].height[layer] : heights[layer]);
            friction_force[layer] = mu * reference.pressure[layer] * contact_area[layer];
            total_normal += reference.pressure[layer] * contact_area[layer];
            total_friction += friction_force[layer];
            const Real radial_increment = radii[0][1] * prescribed_radial_strain(finite, layer) / 10.0L;
            contact_work += reference.pressure[layer] * contact_area[layer] * radial_increment;
            for (std::size_t q = 0; q < 2; ++q) {
                const auto suffix = "_interface_q" + std::to_string(q);
                audit.add("contact_gap",
                    frame.element("contact_gap" + suffix).at(layer),
                    reference.gap[layer],
                    1.0e-14);
                audit.add("contact_pressure",
                    frame.element("contact_pressure" + suffix).at(layer),
                    reference.pressure[layer],
                    1.0e-5);
                audit.add("contact_tangential_traction",
                    frame.element("contact_tangential_traction" + suffix).at(layer),
                    mu * reference.pressure[layer],
                    1.0e-5);
                audit.add("contact_point_area",
                    frame.element("contact_area" + suffix).at(layer),
                    contact_area[layer] / 2.0L,
                    1.0e-15);
                audit.add("contact_heat_flux", frame.element("contact_heat_flux" + suffix).at(layer), 0.0L, 1.0e-9);
                audit.add("contact_elastic_slip",
                    frame.element("contact_elastic_tangential_slip" + suffix).at(layer),
                    reference.elastic_slip[layer][q],
                    1.0e-13);
                audit.add("contact_accumulated_slip",
                    frame.element("contact_total_tangential_slip" + suffix).at(layer),
                    reference.total_slip[layer][q],
                    1.0e-13);
                const double sliding = frame.element("contact_sliding" + suffix).at(layer);
                require(std::isfinite(sliding) && (sliding == 0.0 || sliding == 1.0),
                    "Contact sliding flags must be finite Boolean values, including the exactly touching frame");
                if (step > 5) {
                    if (finite)
                        require(sliding == 1.0, "Finite chain closed-contact points must report forward sliding");
                    // Output re-evaluates the accepted state with zero slip increment.
                    // At the friction limit its Boolean branch is roundoff-sensitive;
                    // establish sliding during this step from the actual saved history.
                    const Real saved_increment =
                        frame.element("contact_total_tangential_slip" + suffix).at(layer)
                        - frames[step - 1].element("contact_total_tangential_slip" + suffix).at(layer);
                    require(saved_increment > elastic_slip_limit,
                        "Every closed-contact point must accumulate forward slip beyond the elastic limit");
                } else if (reference.gap[layer] > 0.0L)
                    require(sliding == 0.0, "Open contact must not report sliding");
                const Real increment = reference.total_slip[layer][q] - previous.total_slip[layer][q];
                contact_work += mu * reference.pressure[layer] * contact_area[layer] / 2.0L * increment;
                friction_dissipation +=
                    mu * reference.pressure[layer] * contact_area[layer] / 2.0L
                    * (increment - reference.elastic_slip[layer][q] + previous.elastic_slip[layer][q]);
                for (std::size_t outer = 2; outer < 4; ++outer)
                    for (const char* field : {"gap",
                             "pressure",
                             "heat_flux",
                             "area",
                             "tangential_traction",
                             "elastic_tangential_slip",
                             "total_tangential_slip",
                             "sliding"})
                        require(std::isnan(frame.element("contact_" + std::string(field) + suffix).at(outer)),
                            "Non-secondary elements must not invent contact quadrature values");
            }
        }
        for (std::size_t body = 0; body < 2; ++body) {
            const auto& current = reference.bodies[body];
            const auto& old = previous.bodies[body];
            const Real friction_sign = body == 0 ? 1.0L : -1.0L;
            std::array<Real, 3> expected_axial_reaction{};
            for (std::size_t station = 0; station < 3; ++station) {
                const auto node = controls[body][station];
                require(role.at(node) == 2.0 && std::isnan(temperature.at(node)) && std::isnan(ur.at(node))
                            && std::isnan(reaction_r.at(node)) && std::isnan(reaction_heat.at(node)),
                    "Axial controls must expose only their active axial field");
                audit.add("axial_control_displacement", uz.at(node), current.displacement[station], 1.0e-12);
                if (station == 1)
                    audit.add("free_middle_displacement", uz.at(node), current.displacement[station], 1.0e-12);
            }
            for (std::size_t layer = 0; layer < 2; ++layer) {
                const auto element = 2 * body + layer;
                require(frame.element("material_point_count").at(element) == 2.0,
                    "Each radial body must retain both material points");
                const Real force = current.area * current.stress[layer][1];
                // The independent diagonal reference is uniform within each
                // layer and has zero shear. Integrate 0.5 sigma:elastic over
                // its reference volume for small strain or current volume for
                // finite strain; retain the previous energy at its own volume.
                const Real volume = current.area * (finite ? current.height[layer] : heights[layer]);
                elastic_energy +=
                    0.5L * volume
                    * ((current.stress[layer][0] + current.stress[layer][2]) * current.radial_elastic[layer]
                        + current.stress[layer][1] * current.axial_elastic[layer]);
                audit.add("axial_section_force", frame.element("axial_force").at(element), force, 1.0e-6);
                audit.add("layer_axial_strain",
                    frame.element("axial_strain").at(element),
                    (current.displacement[layer + 1] - current.displacement[layer]) / heights[layer],
                    1.0e-12);
                expected_axial_reaction[layer] -= force;
                expected_axial_reaction[layer + 1] += force;
                expected_axial_reaction[layer] += friction_sign * friction_force[layer] / 2.0L;
                expected_axial_reaction[layer + 1] += friction_sign * friction_force[layer] / 2.0L;
                body_work += force
                             * ((current.displacement[layer + 1] - old.displacement[layer + 1])
                                 - (current.displacement[layer] - old.displacement[layer]));
                for (std::size_t end = 0; end < 2; ++end) {
                    const auto node = radial_nodes[2 * body + layer][end];
                    const Real prescribed_ur =
                        body == 0 ? prescribed_radial_strain(finite, layer) * time * radii[body][end] : 0.0L;
                    require(role.at(node) == 1.0 && std::isnan(reaction_z.at(node)),
                        "Radial nodes must expose radial fields and interpolated axial displacement");
                    audit.add("radial_temperature", temperature.at(node), 600.0L, 1.0e-10);
                    audit.add("radial_displacement", ur.at(node), prescribed_ur, 1.0e-13);
                    audit.add("radial_node_axial_displacement",
                        uz.at(node),
                        (current.displacement[layer] + current.displacement[layer + 1]) / 2.0L,
                        1.0e-12);
                    audit.add("thermal_reaction", reaction_heat.at(node), 0.0L, 1.0e-8);
                    const Real radius = radii[body][end] * (finite ? current.radial_stretch : 1.0L);
                    const Real height = finite ? current.height[layer] : heights[layer];
                    const Real body_radial =
                        (end == 0 ? -1.0L : 1.0L) * 2.0L * pi * radius * height * current.stress[layer][0];
                    Real reaction = body_radial;
                    if ((body == 0 && end == 1) || (body == 1 && end == 0))
                        reaction += friction_sign * reference.pressure[layer] * contact_area[layer];
                    audit.add("radial_reaction", reaction_r.at(node), reaction, 1.0e-6);
                    const Real increment =
                        body == 0 ? prescribed_radial_strain(finite, layer) / 10.0L * radii[body][end] : 0.0L;
                    body_work += body_radial * increment;
                    dirichlet_work += reaction * increment;
                }
                for (std::size_t q = 0; q < 2; ++q) {
                    const auto suffix = "_q" + std::to_string(q);
                    audit.add("material_temperature",
                        frame.element("temperature" + suffix).at(element),
                        600.0L,
                        1.0e-10);
                    for (std::size_t component = 0; component < components.size(); ++component) {
                        audit.add("stress_" + std::string(components[component]),
                            frame.element("stress_" + std::string(components[component]) + suffix).at(element),
                            current.stress[layer][component],
                            1.0e-3);
                        const Real elastic = component == 3 ? 0.0L
                                                            : (component == 1 ? current.axial_elastic[layer]
                                                                              : current.radial_elastic[layer]);
                        audit.add("elastic_" + std::string(components[component]),
                            frame.element("elastic_" + std::string(components[component]) + suffix).at(element),
                            elastic,
                            1.0e-12);
                        for (const char* mechanism : {"plastic_", "creep_"})
                            audit.add("zero_inelastic_history",
                                frame.element(std::string(mechanism) + components[component] + suffix).at(element),
                                0.0L,
                                1.0e-13);
                    }
                    for (const char* mechanism : {"equiv_plastic", "equiv_creep"})
                        audit.add("zero_inelastic_history",
                            frame.element(std::string(mechanism) + suffix).at(element),
                            0.0L,
                            1.0e-13);
                }
            }
            // The independently predicted intermediate station is unconstrained.
            audit.add("middle_force_balance", reaction_z.at(controls[body][1]), 0.0L, 1.0e-6);
            const Real transmitted = (current.stress[0][1] - current.stress[1][1]) * current.area;
            audit.add("cross_layer_friction_transfer",
                frame.element("axial_force").at(2 * body) - frame.element("axial_force").at(2 * body + 1),
                -friction_sign * (friction_force[0] + friction_force[1]) / 2.0L,
                1.0e-6);
            require(std::abs(transmitted + friction_sign * (friction_force[0] + friction_force[1]) / 2.0L) < 1.0e-8L,
                "Independent closed-form reference violates middle force balance");
            for (const auto station : {0U, 2U}) {
                audit.add("axial_end_reaction",
                    reaction_z.at(controls[body][station]),
                    expected_axial_reaction[station],
                    1.0e-6);
                dirichlet_work +=
                    expected_axial_reaction[station] * (current.displacement[station] - old.displacement[station]);
            }
        }
        const Real observed_end_force =
            static_cast<Real>(reaction_z.at(8)) + reaction_z.at(10) + reaction_z.at(11) + reaction_z.at(13);
        audit.add("global_axial_force_balance", static_cast<double>(observed_end_force), 0.0L, 1.0e-6);
        audit.add("contact_normal_force", frame.global("contact_force_interface"), total_normal, 1.0e-6);
        audit.add("contact_tangential_force",
            frame.global("contact_tangential_force_interface"),
            total_friction,
            1.0e-6);
        audit.add("contact_heat_rate", frame.global("contact_heat_rate_interface"), 0.0L, 1.0e-9);
        audit.add("body_work", frame.global("conservation_internal_mechanical_work_increment"), body_work, 1.0e-10);
        audit.add("elastic_energy_change",
            frame.global("conservation_elastic_energy_change"),
            elastic_energy - previous_elastic_energy,
            1.0e-10);
        audit.add("contact_work", frame.global("conservation_contact_work_increment"), contact_work, 1.0e-10);
        audit.add("dirichlet_work",
            frame.global("conservation_dirichlet_reaction_work_increment"),
            dirichlet_work,
            1.0e-10);
        audit.add("friction_dissipation",
            frame.global("conservation_friction_dissipation_increment"),
            friction_dissipation,
            1.0e-10);
        for (const char* name : {"global_thermal_balance",
                 "interface_heat_imbalance",
                 "mechanical_work_balance",
                 "generated_heat_rate",
                 "stored_heat_rate",
                 "dirichlet_heat_input_rate",
                 "pressure_traction_work_increment",
                 "plastic_dissipation_increment",
                 "creep_dissipation_increment"})
            audit.add(std::string("conservation_") + name,
                frame.global(std::string("conservation_") + name),
                0.0L,
                1.0e-8);
        require(std::abs(frame.global("conservation_relative_mechanical_work_balance")) < 0.001
                    && std::abs(frame.global("conservation_relative_thermal_balance")) < 0.001,
            "Global thermal and mechanical conservation must satisfy the independent balance gate");
        if (step > 5) {
            const Real interpolation = reference.bodies[0].displacement[2] * heights[0] / total_height;
            require(static_cast<Real>(uz.at(9)) < interpolation - 1.0e-9L && uz.at(12) > 1.0e-9,
                "Friction must move both unconstrained intermediate stations away from uncoupled interpolation");
            require(std::abs(frame.element("axial_strain").at(0) - frame.element("axial_strain").at(1)) > 1.0e-7
                        && std::abs(frame.element("axial_strain").at(2) - frame.element("axial_strain").at(3)) > 1.0e-7,
                "Coupled neighboring layers must develop different axial strains");
            ++coupled_frames;
        }
        previous = reference;
        previous_elastic_energy = elastic_energy;
    }
    require(coupled_frames == 5, "Every post-closure frame must demonstrate axial force transfer across both layers");
    const bool passed = audit.passed();
    std::cout << "radial_chain_coupled_frames=" << coupled_frames << '\n'
              << "radial_chain_reference="
              << (finite ? "closed_form_incremental_hughes_winget" : "closed_form_small_strain") << '\n'
              << "radial_chain_validation=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: production_radial_chain <small|finite> <results.e> <summary.csv>\n";
        return 2;
    }
    try {
        const std::string mode = argv[1];
        require(mode == "small" || mode == "finite", "Unknown chain validation mode");
        return check_chain(mode == "finite", argv[2], argv[3]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
