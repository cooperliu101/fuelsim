#include "fuelsim/case_input.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"
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
struct StressReference final {
    std::size_t id;
    fuelsim::SymmetricTensor3Values stress;
};
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}
double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete three-dimensional MOOSE row in " + path);
    return std::stod(values[index]);
}
std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "T,disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional MOOSE nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 4, path)),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path), number(values, 3, path)}});
    }
    return result;
}
std::vector<StressReference> read_stresses(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE stresses: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "id,stress_xx,stress_xy,stress_xz,stress_yy,stress_yz,stress_zz,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional MOOSE stress header in " + path);
    std::vector<StressReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 0, path)),
            {number(values, 1, path), number(values, 4, path), number(values, 6, path), number(values, 2, path),
                number(values, 5, path), number(values, 3, path)}});
    }
    return result;
}
std::array<fuelsim::test::FieldErrorMetrics, 4> compare_nodes(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::SteadyProblem& problem, const std::vector<double>& state,
    const std::vector<NodeReference>& reference, double& maximum_coordinate_difference) {
    if (reference.size() != mesh.nodes().size())
        throw std::invalid_argument("Three-dimensional MOOSE and fuelsim node counts differ");
    std::array<fuelsim::test::FieldErrorMetrics, 4> result;
    std::vector<bool> present(mesh.nodes().size(), false);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (source >= reference.size() || reference[source].id != source || present[source])
                throw std::invalid_argument("Three-dimensional source-node mapping is not unique");
            present[source] = true;
            const auto& actual_point = mesh.nodes()[source];
            const auto& expected = reference[source];
            maximum_coordinate_difference = std::max(maximum_coordinate_difference,
                std::max({std::abs(actual_point.x - expected.point.x), std::abs(actual_point.y - expected.point.y),
                    std::abs(actual_point.z - expected.point.z)}));
            result[0].add(state[dofs.temperature(offset + local)], expected.fields[0]);
            result[1].add(state[dofs.displacement_x(offset + local)], expected.fields[1]);
            result[2].add(state[dofs.displacement_y(offset + local)], expected.fields[2]);
            result[3].add(state[dofs.displacement_z(offset + local)], expected.fields[3]);
        }
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Three-dimensional comparison did not visit every source node");
    return result;
}
std::array<fuelsim::test::FieldErrorMetrics, 6> compare_stresses(const fuelsim::SteadyProblem& problem,
    const std::vector<double>& state, const std::vector<StressReference>& reference) {
    std::array<fuelsim::test::FieldErrorMetrics, 6> result;
    std::vector<bool> present(reference.size(), false);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            const std::size_t source = mesh.source_element_ids()[element];
            if (source >= reference.size() || reference[source].id != source || present[source])
                throw std::invalid_argument("Three-dimensional source-element mapping is not unique");
            present[source] = true;
            const auto stresses = fuelsim::cartesian::ProblemAccess::stress(problem, state, region, element);
            const auto expected = reference[source].stress;
            for (const auto& actual : stresses) {
                result[0].add(actual.xx, expected.xx);
                result[1].add(actual.yy, expected.yy);
                result[2].add(actual.zz, expected.zz);
                result[3].add(actual.xy, expected.xy);
                result[4].add(actual.yz, expected.yz);
                result[5].add(actual.xz, expected.xz);
            }
        }
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Three-dimensional comparison did not visit every source element");
    return result;
}
bool run(const std::string& input_path, const std::string& nodal_path, const std::string& stress_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d)
        throw std::invalid_argument("Stage B comparison requires a steady Cartesian three-dimensional input card");
    const auto mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    const auto solve = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks_per_step, definition.steady_execution.minimum_load_increment},
        options);
    bool passed = check(solve.completed && solve.solve.converged, "stage B input-card solve converges");
    double coordinate_error = 0.0;
    const auto fields = compare_nodes(mesh, problem, solve.solve.state, read_nodes(nodal_path), coordinate_error);
    const auto stresses = compare_stresses(problem, solve.solve.state, read_stresses(stress_path));
    const std::array<std::string, 4> field_names = {
        "temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t field = 0; field < fields.size(); ++field) {
        fuelsim::test::print_relative_metrics("b3_" + field_names[field], fields[field]);
        passed = check(fuelsim::test::relative_metrics_below(fields[field], 1.0e-3) &&
                           fields[field].maximum_zero_reference_difference < 1.0e-10,
                     "stage B " + field_names[field] + " three metrics are below 0.1 percent") &&
                 passed;
    }
    fuelsim::test::print_relative_metrics("b3_stress_xx", stresses[0]);
    passed = check(fuelsim::test::relative_metrics_below(stresses[0], 1.0e-3),
                 "stage B nonzero stress three metrics are below 0.1 percent") &&
             passed;
    for (std::size_t component = 1; component < stresses.size(); ++component)
        passed = check(stresses[component].maximum_absolute_difference < 1.0e-6,
                     "stage B near-zero stress component satisfies its absolute tolerance") &&
                 passed;
    return check(coordinate_error < 1.0e-12, "stage B compares all nodes at matching coordinates") && passed;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_b3_hex8_moose_tests <case.fsi> <nodes.csv> <stresses.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim stage B three-dimensional MOOSE comparison\n");
        if (!run(argv[1], argv[2], argv[3])) return 1;
        std::cout << "[PASS] stage B three-dimensional MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] stage B comparison raised: " << error.what() << '\n';
        return 1;
    }
}
