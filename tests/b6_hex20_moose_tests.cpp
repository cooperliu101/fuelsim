#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
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
struct ReferenceRow final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> values;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    return values;
}

std::vector<ReferenceRow> read_reference(const std::string& path, bool temperature) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read HEX20 MOOSE reference: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected = temperature ? "T,id,x,y,z" : "disp_x,disp_y,disp_z,id,x,y,z";
    if (line != expected) throw std::invalid_argument("Unexpected HEX20 MOOSE reference header: " + path);
    std::vector<ReferenceRow> result;
    while (std::getline(input, line)) {
        const auto values = split(line);
        const std::size_t offset = temperature ? 1U : 3U;
        ReferenceRow row{};
        row.id = static_cast<std::size_t>(std::stoull(values.at(offset)));
        row.point = {
            std::stod(values.at(offset + 1)), std::stod(values.at(offset + 2)), std::stod(values.at(offset + 3))};
        if (temperature)
            row.values[0] = std::stod(values.at(0));
        else
            for (std::size_t component = 0; component < 3; ++component)
                row.values[component] = std::stod(values.at(component));
        result.push_back(row);
    }
    return result;
}

bool run(const std::string& case_path, const std::string& temperature_path, const std::string& displacement_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.preconditioner = definition.solver.preconditioner;
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks_per_step, definition.steady_execution.minimum_load_increment},
        options);
    bool passed = check(solve.completed && solve.solve.converged, "HEX20 fuelsim-to-MOOSE comparison solve converges");
    const auto temperature = read_reference(temperature_path, true);
    const auto displacement = read_reference(displacement_path, false);
    if (temperature.size() != 8 || displacement.size() != 20)
        throw std::invalid_argument("HEX20 MOOSE reference has the wrong mixed-order node counts");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    fuelsim::test::FieldErrorMetrics temperature_error;
    std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_error;
    double coordinate_error = 0.0;
    for (const ReferenceRow& row : temperature) {
        const auto& point = mesh.nodes().at(row.id);
        coordinate_error = std::max(coordinate_error,
            std::max(
                {std::abs(point.x - row.point.x), std::abs(point.y - row.point.y), std::abs(point.z - row.point.z)}));
        temperature_error.add(
            solve.solve.state.at(fields[0].begin + spatial.global_temperature_node(0, row.id)), row.values[0]);
    }
    for (const ReferenceRow& row : displacement) {
        const auto& point = mesh.nodes().at(row.id);
        coordinate_error = std::max(coordinate_error,
            std::max(
                {std::abs(point.x - row.point.x), std::abs(point.y - row.point.y), std::abs(point.z - row.point.z)}));
        const std::size_t global = spatial.global_node(0, row.id);
        for (std::size_t component = 0; component < 3; ++component)
            displacement_error[component].add(
                solve.solve.state.at(fields[component + 1].begin + global), row.values[component]);
    }
    fuelsim::test::print_relative_metrics("b6_hex20_temperature", temperature_error);
    passed = check(fuelsim::test::relative_metrics_below(temperature_error, 1.0e-10),
                 "HEX20 temperature relative L2, relative absolute peak, and maximum pointwise errors are below 1e-8 "
                 "percent") &&
             passed;
    const std::array<std::string, 3> names = {"displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t component = 0; component < 3; ++component) {
        fuelsim::test::print_relative_metrics("b6_hex20_" + names[component], displacement_error[component]);
        passed = check(fuelsim::test::relative_metrics_below(displacement_error[component], 1.0e-8) &&
                           displacement_error[component].maximum_zero_reference_difference < 1.0e-12,
                     "HEX20 " + names[component] +
                         " relative L2, relative absolute peak, and maximum pointwise errors are below 1e-6 percent") &&
                 passed;
    }
    return check(coordinate_error < 1.0e-14, "HEX20 comparison uses identical tracked MOOSE mesh coordinates") &&
           passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_b6_hex20_moose_tests <case.fsi> <temperature.csv> <displacement.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim HEX20-U2/T1 MOOSE comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3]);
        if (passed && session.rank() == 0) std::cout << "[PASS] HEX20-U2/T1 MOOSE comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] HEX20 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
