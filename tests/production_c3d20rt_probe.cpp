#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Expected results.e and Abaqus nodes.csv");
        const auto frames = fuelsim::test::read_exodus_history(argv[1]);
        if (frames.size() != 4)
            throw std::runtime_error("Expected initial state and three completed steps");
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
            for (std::size_t field = 0; field < names.size(); ++field)
                metrics[step - 1][field].add(frame.nodal(names[field]).at(node - 1), reference[field]);
            reference_heat[step - 1][node - 1] = reference[4];
            ++count;
        }
        if (count != 24)
            throw std::runtime_error("Expected 24 reference node samples");
        bool passed = true;
        for (std::size_t step = 0; step < 3; ++step) {
            for (std::size_t field = 0; field < 5; ++field) {
                const auto& metric = metrics[step][field];
                const auto name = "step_" + std::to_string(step + 1) + "_field_" + std::to_string(field);
                if (metric.nonzero_reference_count > 0) {
                    fuelsim::test::print_relative_metrics(name, metric);
                    passed = fuelsim::test::relative_metrics_below(metric, 0.001) && passed;
                } else {
                    fuelsim::test::print_absolute_metrics(name, metric);
                }
                passed = metric.maximum_zero_reference_difference < 1.0e-10 && passed;
            }
            const auto& frame = frames[step + 1];
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
