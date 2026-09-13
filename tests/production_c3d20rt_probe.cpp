#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
constexpr std::array<std::array<double, 3>, 20> probe_nodes = {{{0, 0, 0},
    {1, 0, 0},
    {1, 1, 0},
    {0, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {1, 1, 1},
    {0, 1, 1},
    {0.5, 0, 0},
    {1, 0.5, 0},
    {0.5, 1, 0},
    {0, 0.5, 0},
    {0.5, 0, 1},
    {1, 0.5, 1},
    {0.5, 1, 1},
    {0, 0.5, 1},
    {0, 0, 0.5},
    {1, 0, 0.5},
    {1, 1, 0.5},
    {0, 1, 0.5}}};

struct ProbeThermalReference final {
    std::array<double, 8> initial_mass_capacity{}, current_volume_capacity{}, conduction_increment{}, conduction_hold{};
};

ProbeThermalReference thermal_reference(const std::array<double, 8>& temperature, bool finite, bool nonaffine) {
    // Independent integration of the tracked unit-cube probe, without a Fuelsim element or material call.
    // Corner order is shared by the native CSV and Exodus. At t=0 all temperatures are 300 K.
    constexpr std::array<std::array<int, 3>, 8> corners = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    const double gauss = 1.0 / std::sqrt(3.0);
    ProbeThermalReference result;
    for (double x : {0.5 * (1.0 - gauss), 0.5 * (1.0 + gauss)}) {
        for (double y : {0.5 * (1.0 - gauss), 0.5 * (1.0 + gauss)}) {
            for (double z : {0.5 * (1.0 - gauss), 0.5 * (1.0 + gauss)}) {
                std::array<double, 8> shape{};
                std::array<std::array<double, 3>, 8> gradient{};
                std::array<double, 3> temperature_gradient{};
                double temperature_rate = 0.0;
                for (std::size_t node = 0; node < 8; ++node) {
                    const auto& corner = corners[node];
                    const std::array<double, 3> factor = {corner[0] ? x : 1.0 - x,
                        corner[1] ? y : 1.0 - y,
                        corner[2] ? z : 1.0 - z};
                    shape[node] = factor[0] * factor[1] * factor[2];
                    for (std::size_t d = 0; d < 3; ++d) {
                        gradient[node][d] = (corner[d] ? 1.0 : -1.0) * factor[(d + 1) % 3] * factor[(d + 2) % 3];
                        temperature_gradient[d] += gradient[node][d] * temperature[node];
                    }
                    temperature_rate += shape[node] * (temperature[node] - 300.0); // dt=1 s.
                }
                // u_x=.2X for the affine case; u_x=.4X-.2X^2 for the prescribed midside variant.
                // y=Y and z=Z. The first increment uses the midpoint gradient; the hold uses current geometry.
                const double current_stretch = finite ? (nonaffine ? 1.4 - 0.4 * x : 1.2) : 1.0;
                const double midpoint_stretch = 0.5 * (1.0 + current_stretch);
                for (std::size_t node = 0; node < 8; ++node) {
                    const double capacity = 6.0 * shape[node] * temperature_rate / 8.0; // rho0=2, cp=3.
                    result.initial_mass_capacity[node] += capacity;
                    result.current_volume_capacity[node] += current_stretch * capacity;
                    const double transverse =
                        gradient[node][1] * temperature_gradient[1] + gradient[node][2] * temperature_gradient[2];
                    result.conduction_increment[node] +=
                        10.0 * current_stretch / 8.0
                        * (gradient[node][0] * temperature_gradient[0] / (midpoint_stretch * midpoint_stretch)
                            + transverse);
                    result.conduction_hold[node] +=
                        10.0 * current_stretch / 8.0
                        * (gradient[node][0] * temperature_gradient[0] / (current_stretch * current_stretch)
                            + transverse);
                }
            }
        }
    }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Expected results.e and Abaqus nodes.csv");
        const std::string reference_name = std::filesystem::path(argv[2]).filename().string();
        const bool nonaffine = reference_name == "c3d20rt_nonaffine_probe_nodes.csv";
        const bool finite = nonaffine || reference_name == "c3d20rt_thermal_probe_nodes.csv";
        if (!finite && reference_name != "c3d20rt_small_probe_nodes.csv")
            throw std::invalid_argument("Unknown C3D20RT prescribed thermal probe");
        const auto frames = fuelsim::test::read_exodus_history(argv[1]);
        if (frames.size() != 4)
            throw std::runtime_error("Expected initial state and three completed steps");
        for (std::size_t step = 0; step < frames.size(); ++step) {
            const auto& frame = frames[step];
            if (frame.nodes.size() != probe_nodes.size() || frame.block_element_counts.size() != 1
                || frame.block_element_counts[0] != 1 || frame.time != static_cast<double>(step))
                throw std::runtime_error("Prescribed probe node/element/time coverage changed");
            for (std::size_t node = 0; node < probe_nodes.size(); ++node) {
                if (frame.nodes[node] != probe_nodes[node])
                    throw std::runtime_error("Prescribed probe reference geometry changed");
                const double x = probe_nodes[node][0];
                const double prescribed_x = step == 0 ? 0.0 : (nonaffine || !finite ? 0.4 * x - 0.2 * x * x : 0.2 * x);
                if (!std::isfinite(frame.nodal("displacement_x").at(node))
                    || !std::isfinite(frame.nodal("displacement_y").at(node))
                    || !std::isfinite(frame.nodal("displacement_z").at(node))
                    || std::abs(frame.nodal("displacement_x").at(node) - prescribed_x) >= 1e-10
                    || std::abs(frame.nodal("displacement_y").at(node)) >= 1e-10
                    || std::abs(frame.nodal("displacement_z").at(node)) >= 1e-10)
                    throw std::runtime_error(
                        "Prescribed corner/midside deformation no longer matches the thermal reference");
                if (step == 0 && node < 8 && frame.nodal("temperature").at(node) != 300.0)
                    throw std::runtime_error("Prescribed probe initial temperature changed");
            }
        }
        std::ifstream input(argv[2]);
        std::string line;
        if (!std::getline(input, line))
            throw std::runtime_error("Missing Abaqus node reference");
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line != "time,node,temperature,ux,uy,uz,reaction_heat")
            throw std::runtime_error("Unexpected Abaqus reference schema");
        std::array<std::array<fuelsim::test::FieldErrorMetrics, 5>, 3> metrics{};
        std::array<std::array<double, 8>, 3> reference_heat{};
        std::array<double, 8> reference_temperature{};
        std::size_t count = 0;
        while (std::getline(input, line)) {
            for (char& value : line)
                if (value == ',')
                    value = ' ';
            std::istringstream row(line);
            double time;
            std::size_t node;
            std::array<double, 5> reference{};
            if (!(row >> time >> node))
                throw std::runtime_error("Incomplete reference row");
            for (double& value : reference)
                if (!(row >> value) || !std::isfinite(value))
                    throw std::runtime_error("Invalid reference field");
            const std::size_t step = count / 8 + 1;
            if (step > 3 || node != count % 8 + 1 || time != static_cast<double>(step))
                throw std::runtime_error("Reference time/node sequence mismatch");
            const auto& frame = frames[step];
            if (std::abs(frame.time - time) > 1.0e-12)
                throw std::runtime_error("Production time mismatch");
            const std::array<std::string, 5> names = {"temperature",
                "displacement_x",
                "displacement_y",
                "displacement_z",
                "reaction_heat_flux"};
            // These exact zero constraints are explicit in both tracked input cards.
            // Abaqus stores tiny constraint-solver roundoff (about 1e-33), not a nonzero physical reference.
            for (std::size_t field = 1; field < 4; ++field) {
                const bool zero_constraint = field > 1 || node == 1 || node == 4 || node == 5 || node == 8;
                if (zero_constraint) {
                    if (std::abs(reference[field]) > 1.0e-25)
                        throw std::runtime_error("Abaqus reference violates the prescribed zero displacement");
                    reference[field] = 0.0;
                }
            }
            for (std::size_t field = 0; field < 4; ++field)
                metrics[step - 1][field].add(frame.nodal(names[field]).at(node - 1), reference[field]);
            if (reference[0] != 300.0 + static_cast<double>(node))
                throw std::runtime_error("Native probe temperature no longer matches the prescribed hold history");
            reference_temperature[node - 1] = reference[0];
            reference_heat[step - 1][node - 1] = reference[4];
            ++count;
        }
        if (count != 24)
            throw std::runtime_error("Expected 24 reference node samples");
        const auto reference = thermal_reference(reference_temperature, finite, nonaffine);
        fuelsim::test::FieldErrorMetrics native_capacity, native_conduction;
        for (std::size_t node = 0; node < 8; ++node) {
            native_capacity.add(reference_heat[0][node] - reference.conduction_increment[node],
                reference.current_volume_capacity[node]);
            native_conduction.add(reference_heat[1][node], reference.conduction_hold[node]);
        }
        fuelsim::test::print_relative_metrics("native_current_capacity_identity", native_capacity);
        fuelsim::test::print_relative_metrics("native_hold_conduction_identity", native_conduction);
        if (native_capacity.maximum_absolute_difference >= 1.0e-10
            || native_conduction.maximum_absolute_difference >= 1.0e-10)
            throw std::runtime_error(
                "Native thermal probe does not satisfy the independent capacity/conduction identity");
        std::cout << "step_1_heat_reference="
                  << (finite ? "derived_initial_mass: native_RFL-current_capacity+reference_capacity"
                             : "native_small_strain")
                  << '\n';
        for (std::size_t step = 0; step < 3; ++step) {
            for (std::size_t node = 0; node < 8; ++node) {
                double expected = reference_heat[step][node];
                if (step == 0 && finite)
                    expected += reference.initial_mass_capacity[node] - reference.current_volume_capacity[node];
                metrics[step][4].add(frames[step + 1].nodal("reaction_heat_flux").at(node), expected);
            }
        }
        bool passed = true;
        for (std::size_t step = 0; step < 3; ++step) {
            for (std::size_t field = 0; field < 5; ++field) {
                const auto& metric = metrics[step][field];
                const auto name = "step_" + std::to_string(step + 1) + "_field_" + std::to_string(field)
                                  + (step == 0 && field == 4 && finite ? "_derived_initial_mass" : "");
                if (metric.nonzero_reference_count > 0) {
                    fuelsim::test::print_relative_metrics(name, metric);
                    passed = fuelsim::test::relative_metrics_below(metric, 0.001) && passed;
                } else {
                    fuelsim::test::print_absolute_metrics(name, metric);
                }
                passed = metric.maximum_zero_reference_difference < 1.0e-10 && passed;
            }
            const auto& frame = frames[step + 1];
            double expected_stored = 0.0;
            if (step == 0)
                for (double value : reference.initial_mass_capacity)
                    expected_stored += value;
            const double expected_generated = step == 2 ? 10.0 * (finite ? 1.2 : 1.0) : 0.0;
            std::cout << "step_" << step + 1 << "_stored_heat_rate=" << frame.global("conservation_stored_heat_rate")
                      << " expected=" << expected_stored << '\n'
                      << "step_" << step + 1
                      << "_generated_heat_rate=" << frame.global("conservation_generated_heat_rate")
                      << " expected=" << expected_generated << '\n';
            passed = std::abs(frame.global("conservation_stored_heat_rate") - expected_stored) < 1.0e-10
                     && std::abs(frame.global("conservation_generated_heat_rate") - expected_generated) < 1.0e-10
                     && std::abs(frame.global("conservation_global_thermal_balance")) < 1.0e-10 && passed;
            passed = frame.element("material_point_count").at(0) == 8 && passed;
            for (std::size_t q = 8; q < 27; ++q)
                for (const std::string prefix : {"stress_xx_q", "equiv_plastic_q", "equiv_creep_q"})
                    passed = std::isnan(frame.element(prefix + std::to_string(q)).at(0)) && passed;
        }
        fuelsim::test::FieldErrorMetrics source;
        for (std::size_t node = 0; node < 8; ++node)
            source.add(frames[3].nodal("reaction_heat_flux")[node] - frames[2].nodal("reaction_heat_flux")[node],
                reference_heat[2][node] - reference_heat[1][node]);
        fuelsim::test::print_relative_metrics("isolated_source", source);
        passed = fuelsim::test::relative_metrics_below(source, 0.001) && passed;
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
