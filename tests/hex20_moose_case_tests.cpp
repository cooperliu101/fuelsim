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
    double temperature = 0.0;
    std::array<double, 3> displacement{};
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
        if (line.empty()) continue;
        const auto values = split(line);
        const std::size_t offset = temperature ? 1U : 3U;
        ReferenceRow row{};
        row.id = static_cast<std::size_t>(std::stoull(values.at(offset)));
        row.point = {
            std::stod(values.at(offset + 1)), std::stod(values.at(offset + 2)), std::stod(values.at(offset + 3))};
        if (temperature)
            row.temperature = std::stod(values.at(0));
        else
            for (std::size_t component = 0; component < 3; ++component)
                row.displacement[component] = std::stod(values.at(component));
        result.push_back(row);
    }
    return result;
}

bool compare_case(const std::string& case_path, const std::string& temperature_path,
    const std::string& displacement_path, double tolerance, bool compare_displacement) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.geometry != fuelsim::CaseGeometry::cartesian_3d ||
        definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument("HEX20 MOOSE case comparison requires a steady Cartesian three-dimensional case");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.preconditioner = definition.solver.preconditioner;
    options.field_residual_scaling = definition.solver.field_residual_scaling;
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem, definition.steady_execution, options);
    bool passed = check(solve.completed && solve.solve.converged, "HEX20 MOOSE comparison solve converges");
    const auto temperature = read_reference(temperature_path, true);
    const auto displacement = read_reference(displacement_path, false);
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
            solve.solve.state.at(fields[0].begin + spatial.global_temperature_node(0, row.id)), row.temperature);
    }
    for (const ReferenceRow& row : displacement) {
        const auto& point = mesh.nodes().at(row.id);
        coordinate_error = std::max(coordinate_error,
            std::max(
                {std::abs(point.x - row.point.x), std::abs(point.y - row.point.y), std::abs(point.z - row.point.z)}));
        const std::size_t global = spatial.global_node(0, row.id);
        for (std::size_t component = 0; component < 3; ++component)
            displacement_error[component].add(
                solve.solve.state.at(fields[component + 1].begin + global), row.displacement[component]);
    }
    fuelsim::test::print_relative_metrics("hex20_case_temperature", temperature_error);
    passed = check(fuelsim::test::relative_metrics_below(temperature_error, tolerance),
                 "HEX20 temperature matches the MOOSE reference") &&
             passed;
    if (compare_displacement) {
        const std::array<std::string, 3> names = {"displacement_x", "displacement_y", "displacement_z"};
        for (std::size_t component = 0; component < 3; ++component) {
            fuelsim::test::print_relative_metrics("hex20_case_" + names[component], displacement_error[component]);
            passed = check(fuelsim::test::relative_metrics_below(displacement_error[component], tolerance),
                         "HEX20 " + names[component] + " matches the MOOSE reference") &&
                     passed;
        }
    }
    return check(coordinate_error < 1.0e-14, "HEX20 comparison uses identical MOOSE mesh coordinates") && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5 && argc != 6) {
        std::cerr << "Usage: fuelsim_hex20_moose_case_tests <case.fsi> <temperature.csv> <displacement.csv> "
                     "<tolerance> [temperature_only]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim HEX20 MOOSE case comparison\n");
        const bool passed = compare_case(argv[1], argv[2], argv[3], std::stod(argv[4]), argc == 5);
        if (passed && session.rank() == 0) std::cout << "[PASS] HEX20 MOOSE case comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] HEX20 MOOSE comparison raised: " << error.what() << '\n';
        return 1;
    }
}
