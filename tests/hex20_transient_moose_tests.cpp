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
struct NodeReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> fields;
};

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

std::size_t column(const std::vector<std::string>& header, const std::string& name, const std::string& path) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("Missing column '" + name + "' in " + path);
    return static_cast<std::size_t>(found - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete MOOSE row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_temperature(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read HEX20 transient MOOSE nodes: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split(line);
    const std::array<std::size_t, 5> columns = {column(header, "id", path), column(header, "x", path),
        column(header, "y", path), column(header, "z", path), column(header, "T", path)};
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({static_cast<std::size_t>(number(values, columns[0], path)),
            {number(values, columns[1], path), number(values, columns[2], path), number(values, columns[3], path)},
            {number(values, columns[4], path), 0.0, 0.0, 0.0}});
    }
    return result;
}

std::vector<NodeReference> read_displacement(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read HEX20 transient MOOSE displacement nodes: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split(line);
    const std::array<std::size_t, 7> columns = {column(header, "id", path), column(header, "x", path),
        column(header, "y", path), column(header, "z", path), column(header, "disp_x", path),
        column(header, "disp_y", path), column(header, "disp_z", path)};
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({static_cast<std::size_t>(number(values, columns[0], path)),
            {number(values, columns[1], path), number(values, columns[2], path), number(values, columns[3], path)},
            {0.0, number(values, columns[4], path), number(values, columns[5], path),
                number(values, columns[6], path)}});
    }
    return result;
}

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool run(const std::string& case_path, const std::string& temperature_path, const std::string& displacement_path,
    double tolerance) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::transient ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d)
        throw std::invalid_argument("HEX20 transient comparison requires a transient Cartesian three-dimensional case");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.preconditioner = definition.solver.preconditioner;
    const auto& execution = definition.transient_execution;
    const fuelsim::TransientResult solve = fuelsim::solve_transient(problem,
        {execution.end_time, execution.initial_time_step, execution.minimum_time_step, execution.maximum_time_step,
            execution.growth_factor, execution.cutback_factor, execution.maximum_cutbacks_per_step,
            execution.load_ramp_time},
        options);
    bool passed = check(solve.completed, "HEX20 transient comparison solve converges") &&
                  check(solve.accepted_steps.size() == 10, "HEX20 transient comparison accepts ten fixed steps");
    const auto temperatures = read_temperature(temperature_path);
    const auto displacements = read_displacement(displacement_path);
    if (temperatures.size() != 8 || displacements.size() != mesh.nodes().size())
        throw std::invalid_argument("MOOSE and fuelsim node counts differ");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    std::array<fuelsim::test::FieldErrorMetrics, 4> metrics;
    double coordinate_error = 0.0;
    for (const NodeReference& reference : temperatures) {
        const auto& point = mesh.nodes().at(reference.id);
        coordinate_error = std::max(
            coordinate_error, std::max({std::abs(point.x - reference.point.x), std::abs(point.y - reference.point.y),
                                  std::abs(point.z - reference.point.z)}));
        metrics[0].add(
            problem.committed_solution().at(fields[0].begin + spatial.global_temperature_node(0, reference.id)),
            reference.fields[0]);
    }
    for (const NodeReference& reference : displacements) {
        const auto& point = mesh.nodes().at(reference.id);
        coordinate_error = std::max(
            coordinate_error, std::max({std::abs(point.x - reference.point.x), std::abs(point.y - reference.point.y),
                                  std::abs(point.z - reference.point.z)}));
        const std::size_t global = spatial.global_node(0, reference.id);
        for (std::size_t component = 0; component < 3; ++component)
            metrics[component + 1].add(
                problem.committed_solution().at(fields[component + 1].begin + global), reference.fields[component + 1]);
    }
    const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t field = 0; field < metrics.size(); ++field) {
        fuelsim::test::print_relative_metrics("hex20_transient_" + names[field], metrics[field]);
        passed = check(fuelsim::test::relative_metrics_below(metrics[field], tolerance),
                     "HEX20 transient " + names[field] + " matches the MOOSE reference") &&
                 passed;
    }
    return check(coordinate_error < 1.0e-14, "HEX20 transient comparison uses identical mesh coordinates") && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_hex20_transient_moose_tests <case.fsi> <temperature.csv> "
                     "<displacement.csv> <tolerance>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim HEX20 transient MOOSE comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3], std::stod(argv[4]));
        if (passed && session.rank() == 0) std::cout << "[PASS] HEX20 transient MOOSE comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] HEX20 transient comparison raised: " << error.what() << '\n';
        return 1;
    }
}
