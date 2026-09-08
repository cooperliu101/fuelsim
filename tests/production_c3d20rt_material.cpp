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
#include <vector>

namespace {
std::vector<double> row_values(const std::string& line) {
    std::istringstream row(line);
    std::vector<double> result;
    std::string value;
    while (std::getline(row, value, ',')) {
        const double number = std::stod(value);
        if (!std::isfinite(number))
            throw std::runtime_error("Nonfinite reference value");
        result.push_back(number);
    }
    return result;
}

double exact_zero(double value, double tolerance) {
    if (std::abs(value) > tolerance)
        throw std::runtime_error("Reference violates uniaxial symmetry/zero constraint");
    return 0.0;
}

void header(std::ifstream& file) {
    std::string line;
    if (!std::getline(file, line))
        throw std::runtime_error("Missing reference header");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5)
            throw std::invalid_argument("Expected results.e nodes.csv points.csv coupled|plastic|creep|creep_hold");
        const std::string mode = argv[4];
        const bool hold = mode == "creep_hold";
        const bool plastic = mode == "coupled" || mode == "plastic",
                   creep = mode == "coupled" || mode == "creep" || hold;
        if (!plastic && !creep)
            throw std::invalid_argument("Unknown material mode");
        const auto frames = fuelsim::test::read_exodus_history(argv[1]);
        if (frames.size() != 11)
            throw std::runtime_error("Expected initial state plus ten increments");
        std::array<fuelsim::test::FieldErrorMetrics, 7> nodal{};
        std::array<fuelsim::test::FieldErrorMetrics, 26> material{};
        std::array<fuelsim::test::FieldErrorMetrics, 4> hold_analytic{};
        double reference_zero_force = 0.0, reference_zero_stress = 0.0;
        const std::array<std::string, 7> node_names = {"temperature",
            "displacement_x",
            "displacement_y",
            "displacement_z",
            "reaction_force_x",
            "reaction_force_y",
            "reaction_force_z"};
        const std::array<std::string, 6> components = {"xx", "yy", "zz", "xy", "yz", "xz"};
        const std::array<std::string, 4> prefixes = {"stress_", "elastic_", "plastic_", "creep_"};
        std::ifstream nodes(argv[2]), points(argv[3]);
        header(nodes);
        header(points);
        std::string line;
        std::size_t count = 0;
        while (std::getline(nodes, line)) {
            const auto values = row_values(line);
            const std::size_t step = count / 20 + 1, node = count % 20;
            if (values.size() != 9 || step > 10 || values[1] != static_cast<double>(node + 1)
                || std::abs(values[0] - frames[step].time) > 1.0e-7)
                throw std::runtime_error("Node reference sequence mismatch");
            for (std::size_t field = 0; field < 7; ++field) {
                double reference = values[field + 2];
                if (field >= 1 && field <= 3 && frames[step].nodes[node][field - 1] == 0.0)
                    reference = exact_zero(reference, 1.0e-12);
                if (field >= 5 || (field == 4 && frames[step].nodes[node][0] == 0.5)) {
                    reference_zero_force = std::max(reference_zero_force, std::abs(reference));
                    reference = exact_zero(reference, 10.0);
                }
                nodal[field].add(frames[step].nodal(node_names[field])[node], reference);
            }
            ++count;
        }
        if (count != 200)
            throw std::runtime_error("Expected all 200 node samples");
        count = 0;
        std::array<double, 8> previous_p{}, previous_c{}, reference_previous_p{}, reference_previous_c{};
        while (std::getline(points, line)) {
            const auto values = row_values(line);
            const std::size_t step = count / 8 + 1, q = count % 8;
            if (values.size() != 29 || step > 10 || values[1] != 1 || values[2] != static_cast<double>(q + 1)
                || std::abs(values[0] - frames[step].time) > 1.0e-7)
                throw std::runtime_error("Material reference sequence mismatch");
            const auto& frame = frames[step];
            if (hold && std::abs(frame.time - (1e-9 + 0.1 * static_cast<double>(step - 1))) > 1e-12)
                throw std::runtime_error("Hold case has an unexpected physical time");
            if (frame.element("material_point_count")[0] != 8)
                throw std::runtime_error("Wrong active point count");
            for (std::size_t field = 0; field < 26; ++field) {
                const std::string name = field < 24 ? prefixes[field / 6] + components[field % 6] + "_q"
                                                    : (field == 24 ? "equiv_plastic_q" : "equiv_creep_q");
                const double actual = frame.element(name + std::to_string(q))[0];
                double reference = values[field + 3];
                // This single-axis loading has exactly zero transverse/shear stress and shear histories.
                if (field > 0 && field < 6) {
                    reference_zero_stress = std::max(reference_zero_stress, std::abs(reference));
                    reference = exact_zero(reference, 10.0);
                }
                if (field >= 6 && field < 24 && field % 6 >= 3)
                    reference = exact_zero(reference, 1.0e-12);
                material[field].add(actual, reference);
                if (hold && field == 0) {
                    hold_analytic[0].add(actual, 1e8);
                    hold_analytic[1].add(reference, 1e8);
                }
                if (field == 24) {
                    if (plastic && (!(actual > previous_p[q]) || !(reference > reference_previous_p[q])))
                        throw std::runtime_error("Plasticity did not advance in every prescribed loading increment");
                    previous_p[q] = actual;
                    reference_previous_p[q] = reference;
                }
                if (field == 25) {
                    if (creep && (!hold || step > 1)
                        && (!(actual > previous_c[q]) || !(reference > reference_previous_c[q])))
                        throw std::runtime_error(
                            "Creep did not advance in every prescribed loading increment: step=" + std::to_string(step)
                            + " actual=" + std::to_string(actual) + " reference=" + std::to_string(reference));
                    previous_c[q] = actual;
                    reference_previous_c[q] = reference;
                    if (hold) {
                        const double exact_creep = 1e-4 * 0.1 * static_cast<double>(step - 1);
                        hold_analytic[2].add(actual, exact_creep);
                        hold_analytic[3].add(reference, exact_creep);
                    }
                }
            }
            ++count;
        }
        if (count != 80)
            throw std::runtime_error("Expected all 80 material point samples");
        bool passed = true;
        std::cout << std::scientific << std::setprecision(12);
        for (std::size_t field = 0; field < 33; ++field) {
            const auto& metric = field < 7 ? nodal[field] : material[field - 7];
            const std::string name = field < 7 ? node_names[field] : "material_" + std::to_string(field - 7);
            if (metric.nonzero_reference_count > 0) {
                fuelsim::test::print_relative_metrics(name, metric);
                passed = fuelsim::test::relative_metrics_below(metric, 0.001) && passed;
            } else
                fuelsim::test::print_absolute_metrics(name, metric);
            const double zero_tolerance = (field >= 4 && field <= 12) ? 0.1 : 1.0e-12;
            passed = metric.maximum_zero_reference_difference < zero_tolerance && passed;
        }
        if (hold) {
            const std::array<std::string, 4> names = {"hold_fuelsim_stress",
                "hold_abaqus_stress",
                "hold_fuelsim_creep",
                "hold_abaqus_creep"};
            for (std::size_t i = 0; i < hold_analytic.size(); ++i) {
                fuelsim::test::print_relative_metrics(names[i], hold_analytic[i]);
                passed = fuelsim::test::relative_metrics_below(hold_analytic[i], 0.001)
                         && hold_analytic[i].maximum_zero_reference_difference < 1e-12 && passed;
            }
            std::cout << "creep_active_every_hold_increment=true\n";
        }
        std::cout << "plastic_active_every_increment=" << plastic
                  << "\ncreep_active_every_increment=" << (creep && !hold) << '\n';
        std::cout << "abaqus_analytic_zero_force_maximum_absolute=" << reference_zero_force << '\n'
                  << "abaqus_analytic_zero_stress_maximum_absolute=" << reference_zero_stress << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
