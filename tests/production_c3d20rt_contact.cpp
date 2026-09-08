#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
std::vector<double> numbers(const std::string& line) {
    std::vector<double> values;
    std::istringstream stream(line);
    std::string token;
    while (std::getline(stream, token, ',')) {
        const double value = std::stod(token);
        if (!std::isfinite(value))
            throw std::runtime_error("Nonfinite Abaqus reference");
        values.push_back(value);
    }
    return values;
}

void skip_header(std::ifstream& input) {
    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("Missing reference header");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 6 && argc != 7 && argc != 8)
            throw std::invalid_argument(
                "Expected results.e nodes.csv contact.csv points.csv increments [material_points "
                "[--qualified-finite-frictionless|--qualified-finite-friction]]");
        const std::string policy = argc == 8 ? argv[7] : "";
        const bool qualified_frictionless = policy == "--qualified-finite-frictionless";
        const bool qualified_friction = policy == "--qualified-finite-friction";
        if (!policy.empty() && !qualified_frictionless && !qualified_friction)
            throw std::invalid_argument("Unknown contact qualification option");
        const std::string point_text = argc >= 7 ? argv[6] : "8";
        if (point_text != "8" && point_text != "27")
            throw std::invalid_argument("Contact material point count must be 8 or 27");
        const std::size_t material_points = std::stoul(point_text), points_per_frame = 6 * material_points;
        const std::string increment_text(argv[5]);
        if (increment_text.empty() || increment_text.find_first_not_of("0123456789") != std::string::npos)
            throw std::invalid_argument("Increment count must be a positive integer");
        const std::size_t increments = std::stoul(increment_text);
        if (increments == 0 || increments > 100000)
            throw std::invalid_argument("Increment count is outside the supported range");
        if ((qualified_frictionless || qualified_friction) && material_points != 8)
            throw std::invalid_argument("Finite contact qualification requires C3D20RT");
        if (qualified_friction && increments != 280)
            throw std::invalid_argument("Finite friction qualification requires 280 increments");
        if (qualified_frictionless && increments != 80 && increments != 280)
            throw std::invalid_argument("Finite frictionless qualification requires 80 or 280 increments");
        const auto history = fuelsim::test::read_exodus_history(argv[1]);
        if (history.size() != increments + 1)
            throw std::runtime_error("Expected initial state and every requested increment");
        fuelsim::test::FieldErrorMetrics temperature, thermal_reaction, pressure, gap, heat;
        fuelsim::test::GroupedFieldErrorMetrics displacement, reaction, normal, tangent, stress, slip;
        const double tangent_pointwise_tolerance = qualified_friction ? 0.0015 : 0.001;
        const double slip_pointwise_tolerance = qualified_friction ? 0.006 : (qualified_frictionless ? 0.01 : 0.001);
        const double reaction_relative_tolerance = 0.0015;
        const double absolute_reaction_tolerance = 10.0; // N, Euclidean vector difference
        std::size_t qualified_friction_samples = 0, reaction_pointwise_failures = 0;
        double qualified_friction_maximum_difference = 0.0;
        std::vector<double> reference_heat(increments, 0.0);
        std::string line;
        std::ifstream nodes(argv[2]), contacts(argv[3]), points(argv[4]);
        skip_header(nodes);
        skip_header(contacts);
        skip_header(points);
        std::size_t count = 0;
        const std::array<std::size_t, 32> temperature_nodes = {1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            21,
            22,
            23,
            24,
            33,
            34,
            35,
            36,
            45,
            46,
            47,
            48,
            57,
            58,
            59,
            60,
            61,
            62,
            63,
            64,
            77,
            78,
            79,
            80};
        while (std::getline(nodes, line)) {
            const auto values = numbers(line);
            const std::size_t frame_index = count / 88 + 1, node = count % 88;
            if (values.size() != 10 || frame_index > increments || values[1] != static_cast<double>(node + 1)
                || std::abs(values[0] - history[frame_index].time) > 1.0e-6)
                throw std::runtime_error("Node time/label sequence mismatch");
            const auto& frame = history[frame_index];
            temperature.add(frame.nodal("temperature")[node], values[2]);
            std::array<double, 3> actual_u{}, actual_r{}, reference_u{}, reference_r{};
            for (std::size_t component = 0; component < 3; ++component) {
                const std::string axis(1, "xyz"[component]);
                actual_u[component] = frame.nodal("displacement_" + axis)[node];
                actual_r[component] = frame.nodal("reaction_force_" + axis)[node];
                reference_u[component] = values[3 + component];
                reference_r[component] = values[6 + component];
                if (frame.nodes[node][0] == 0.0) {
                    if (std::abs(reference_u[component]) > 1.0e-12)
                        throw std::runtime_error("Fixed primary boundary moved in Abaqus reference");
                    reference_u[component] = 0.0;
                }
            }
            displacement.add(actual_u.data(), reference_u.data(), 3);
            if (std::find(temperature_nodes.begin(), temperature_nodes.end(), node + 1) != temperature_nodes.end())
                thermal_reaction.add(frame.nodal("reaction_heat_flux")[node], values[9]);
            else if (!std::isnan(frame.nodal("reaction_heat_flux")[node]) || values[9] != 0.0)
                throw std::runtime_error("Inactive midpoint thermal reaction is not marked as unavailable");
            reaction.add(actual_r.data(), reference_r.data(), 3);
            if (qualified_friction) {
                const double difference = std::hypot(actual_r[0] - reference_r[0],
                    actual_r[1] - reference_r[1],
                    actual_r[2] - reference_r[2]);
                const double reference_norm = std::hypot(reference_r[0], reference_r[1], reference_r[2]);
                // Apply the OR rule at each node and time, not to global maxima.
                if (!(reference_norm > 0.0 && difference / reference_norm < reaction_relative_tolerance)) {
                    ++qualified_friction_samples;
                    qualified_friction_maximum_difference = std::max(qualified_friction_maximum_difference, difference);
                    if (!(difference < absolute_reaction_tolerance))
                        ++reaction_pointwise_failures;
                }
            }
            ++count;
        }
        if (count != 88 * increments)
            throw std::runtime_error("Expected every node of every increment");
        count = 0;
        std::size_t sticking = 0, sliding = 0;
        bool compressed_sliding_observed = false;
        const std::array<std::size_t, 13> secondary = {57, 60, 61, 64, 68, 69, 72, 76, 78, 80, 83, 85, 88};
        while (std::getline(contacts, line)) {
            const auto values = numbers(line);
            const std::size_t frame_index = count / 13 + 1, node = secondary[count % 13] - 1;
            if (values.size() != 14 || frame_index > increments || values[1] != static_cast<double>(node + 1)
                || std::abs(values[0] - history[frame_index].time) > 1.0e-6)
                throw std::runtime_error("Contact time/label sequence mismatch");
            const auto& frame = history[frame_index];
            gap.add(frame.nodal("contact_gap_interface")[node], values[2]);
            pressure.add(frame.nodal("contact_pressure_interface")[node], values[3]);
            // Fuelsim stores the contact residual; Abaqus reports the applied contact force.
            std::array<double, 3> actual_normal{}, actual_tangent{}, reference_normal{}, reference_tangent{};
            for (std::size_t component = 0; component < 3; ++component) {
                const std::string axis(1, "xyz"[component]);
                actual_normal[component] = -frame.nodal("contact_normal_force_" + axis + "_interface")[node];
                actual_tangent[component] = -frame.nodal("contact_tangential_force_" + axis + "_interface")[node];
                reference_normal[component] = values[4 + component];
                reference_tangent[component] = values[7 + component];
            }
            normal.add(actual_normal.data(), reference_normal.data(), 3);
            tangent.add(actual_tangent.data(), reference_tangent.data(), 3);
            std::array<double, 3> actual_slip{}, reference_slip{};
            for (std::size_t component = 0; component < 3; ++component) {
                const std::string axis(1, "xyz"[component]);
                actual_slip[component] = frame.nodal("contact_total_slip_" + axis + "_interface")[node];
                reference_slip[component] = values[10 + component];
            }
            slip.add(actual_slip.data(), reference_slip.data(), 3);
            if (frame.time > 1.0 && frame.nodal("contact_pressure_interface")[node] > 0.0
                && std::hypot(actual_slip[0], actual_slip[1], actual_slip[2]) > 1e-6)
                compressed_sliding_observed = true;
            reference_heat[frame_index - 1] += values[13];
            if (frame.nodal("contact_projected_interface")[node] != 1)
                throw std::runtime_error("Contact projection was lost");
            if (frame.nodal("contact_sliding_interface")[node] == 1)
                ++sliding;
            else
                ++sticking;
            ++count;
        }
        if (count != 13 * increments)
            throw std::runtime_error("Expected every secondary contact node at every increment");
        for (std::size_t step = 1; step <= increments; ++step)
            heat.add(history[step].global("contact_heat_rate_interface"), reference_heat[step - 1]);
        count = 0;
        const std::array<std::string, 6> components = {"xx", "yy", "zz", "xy", "yz", "xz"};
        while (std::getline(points, line)) {
            const auto values = numbers(line);
            const std::size_t frame_index = count / points_per_frame + 1,
                              element = count % points_per_frame / material_points, q = count % material_points;
            if (values.size() != 9 || frame_index > increments || values[1] != static_cast<double>(element + 1)
                || values[2] != static_cast<double>(q + 1) || std::abs(values[0] - history[frame_index].time) > 1.0e-6)
                throw std::runtime_error("Stress time/element/point sequence mismatch");
            const auto& frame = history[frame_index];
            std::array<double, 6> actual{}, reference{};
            for (std::size_t c = 0; c < 6; ++c) {
                actual[c] = frame.element("stress_" + components[c] + "_q" + std::to_string(q))[element];
                reference[c] = values[3 + c];
            }
            stress.add(actual.data(), reference.data(), 6);
            if (frame.element("material_point_count")[element] != static_cast<double>(material_points))
                throw std::runtime_error("Contact material output has wrong active point count");
            ++count;
        }
        if (count != points_per_frame * increments)
            throw std::runtime_error("Expected all material points at all increments");
        std::cout << std::scientific << std::setprecision(12);
        // User-approved exception for small reaction/slip values in this finite
        // frictionless case. Keep raw metrics and aggregate tolerances unchanged.
        const double pointwise_tolerance = qualified_frictionless ? 0.01 : 0.001;
        // User-approved per-sample relative-first, absolute-fallback reaction policy.
        const double zero_reaction_tolerance =
            qualified_friction ? absolute_reaction_tolerance : (qualified_frictionless ? 0.1 : 1e-5);
        if (qualified_friction && !tangent.has_relative_norm())
            throw std::runtime_error("Absolute reaction policy requires a frictional reference");
        if (qualified_frictionless && (tangent.has_relative_norm() || tangent.maximum_zero_reference_difference != 0.0))
            throw std::runtime_error("Finite frictionless qualification requires strictly zero tangential forces");
        std::cout << "acceptance_policy="
                  << (qualified_friction ? "qualified_finite_friction"
                                         : (qualified_frictionless ? "qualified_finite_frictionless" : "strict"))
                  << '\n'
                  << "aggregate_relative_tolerance=" << 0.001 << '\n'
                  << "reaction_pointwise_relative_tolerance="
                  << (qualified_friction ? reaction_relative_tolerance : pointwise_tolerance) << '\n'
                  << "tangential_force_pointwise_relative_tolerance=" << tangent_pointwise_tolerance << '\n'
                  << "slip_pointwise_relative_tolerance=" << slip_pointwise_tolerance << '\n'
                  << "reaction_zero_reference_absolute_tolerance=" << zero_reaction_tolerance << '\n';
        if (qualified_friction)
            std::cout << "reaction_pointwise_absolute_tolerance_N=" << absolute_reaction_tolerance << '\n'
                      << "reaction_absolute_fallback_samples=" << qualified_friction_samples << '\n'
                      << "reaction_absolute_fallback_maximum_difference_N=" << qualified_friction_maximum_difference
                      << '\n'
                      << "reaction_pointwise_failures=" << reaction_pointwise_failures << '\n';
        if (qualified_frictionless && increments == 80) {
            if (std::abs(history.back().time - 2.0) > 1e-10 || !compressed_sliding_observed || !heat.has_relative_norm()
                || !normal.has_relative_norm())
                throw std::runtime_error(
                    "Short frictionless case must cover compression, geometric sliding and heat transfer");
            std::cout << "short_frictionless_coverage_passed=1\n";
        }
        // The sliding flag denotes Coulomb friction activation, not geometric
        // motion. A frictionless history has zero tangential force and can keep
        // this flag unset while accumulating nonzero geometric slip.
        const bool frictionless_reference = !tangent.has_relative_norm();
        bool passed = frictionless_reference
                          ? tangent.maximum_zero_reference_difference == 0.0 && slip.has_relative_norm()
                          : sticking > 0 && sliding > 0;
        const std::array<const fuelsim::test::FieldErrorMetrics*, 5> scalars = {&temperature,
            &thermal_reaction,
            &pressure,
            &gap,
            &heat};
        const std::array<std::string, 5> names = {"temperature", "thermal_reaction", "pressure", "gap", "contact_heat"};
        const std::array<double, 5> zero_tolerances = {1e-10, 1e-7, 1e-4, 1e-10, 1e-7};
        for (std::size_t i = 0; i < 5; ++i) {
            passed = scalars[i]->maximum_zero_reference_difference < zero_tolerances[i] && passed;
            if (scalars[i]->has_relative_norm()) {
                fuelsim::test::print_relative_metrics(names[i], *scalars[i]);
                passed = fuelsim::test::relative_metrics_below(*scalars[i], 0.001) && passed;
            } else
                fuelsim::test::print_absolute_metrics(names[i], *scalars[i]);
        }
        for (const auto& item : std::array<std::pair<std::string, const fuelsim::test::GroupedFieldErrorMetrics*>, 6>{
                 {{"displacement", &displacement},
                     {"reaction", &reaction},
                     {"normal_force", &normal},
                     {"tangential_force", &tangent},
                     {"stress", &stress},
                     {"slip", &slip}}}) {
            if (item.second->has_relative_norm()) {
                fuelsim::test::print_grouped_relative_metrics(item.first, *item.second);
                if (qualified_friction && item.first == "reaction") {
                    passed = reaction.relative_l2() < 0.001 && reaction.relative_absolute_peak() < 0.001
                             && reaction_pointwise_failures == 0 && passed;
                } else if (qualified_friction && (item.first == "tangential_force" || item.first == "slip")) {
                    const double tolerance =
                        item.first == "slip" ? slip_pointwise_tolerance : tangent_pointwise_tolerance;
                    passed = item.second->relative_l2() < 0.001 && item.second->relative_absolute_peak() < 0.001
                             && item.second->maximum_pointwise_relative < tolerance && passed;
                } else if (qualified_frictionless && (item.first == "reaction" || item.first == "slip")) {
                    passed = item.second->relative_l2() < 0.001 && item.second->relative_absolute_peak() < 0.001
                             && item.second->maximum_pointwise_relative < pointwise_tolerance && passed;
                } else {
                    passed = fuelsim::test::grouped_relative_metrics_below(*item.second, 0.001) && passed;
                }
            } else {
                std::cout << item.first << "_zero_reference_count=" << item.second->zero_reference_count << '\n'
                          << item.first << "_maximum_zero_reference_absolute_difference="
                          << item.second->maximum_zero_reference_difference << '\n';
            }
        }
        passed = displacement.maximum_zero_reference_difference < 1e-10
                 && reaction.maximum_zero_reference_difference < zero_reaction_tolerance
                 && normal.maximum_zero_reference_difference < 1e-5 && tangent.maximum_zero_reference_difference < 1e-5
                 && stress.maximum_zero_reference_difference < 1e-3 && slip.maximum_zero_reference_difference < 1e-10
                 && passed;
        std::cout << "sticking_samples=" << sticking << "\nsliding_samples=" << sliding << '\n';
        std::cout << "acceptance_passed=" << (passed ? 1 : 0) << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
