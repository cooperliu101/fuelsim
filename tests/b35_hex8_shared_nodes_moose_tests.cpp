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
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> fields;
};

struct NodalComparison final {
    std::array<fuelsim::test::FieldErrorMetrics, 4> fields;
    std::vector<std::size_t> source_global_nodes;
    double maximum_coordinate_difference = 0.0;
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

std::vector<NodeReference> read_nodes(const std::string& path, std::size_t node_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "T,disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional MOOSE nodal header in " + path);
    std::vector<NodeReference> result(node_count);
    std::vector<bool> present(node_count, false);
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        const std::size_t id = static_cast<std::size_t>(number(values, 4, path));
        if (id >= node_count) throw std::invalid_argument("MOOSE node identifier is outside the Exodus mesh");
        const NodeReference candidate = {{number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path), number(values, 3, path)}};
        if (present[id]) {
            const NodeReference& existing = result[id];
            const bool identical = std::abs(existing.point.x - candidate.point.x) < 1.0e-12 &&
                                   std::abs(existing.point.y - candidate.point.y) < 1.0e-12 &&
                                   std::abs(existing.point.z - candidate.point.z) < 1.0e-12 &&
                                   std::equal(
                                       existing.fields.begin(), existing.fields.end(), candidate.fields.begin(),
                                       [](double left, double right) { return std::abs(left - right) < 1.0e-12; });
            if (!identical)
                throw std::invalid_argument("MOOSE emitted inconsistent duplicate values for a shared node");
            continue;
        }
        present[id] = true;
        result[id] = candidate;
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("MOOSE nodal reference does not contain every Exodus node");
    return result;
}

NodalComparison compare_nodes(const fuelsim::UnstructuredHex8Mesh& mesh, const fuelsim::SteadyProblem& problem,
    const std::vector<double>& state, const std::vector<NodeReference>& reference) {
    if (reference.size() != mesh.nodes().size())
        throw std::invalid_argument("Three-dimensional MOOSE and fuelsim node counts differ");
    NodalComparison result;
    result.source_global_nodes.assign(mesh.nodes().size(), std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> source_occurrences(mesh.nodes().size(), 0U);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (source >= mesh.nodes().size()) throw std::invalid_argument("Three-dimensional source node is invalid");
            const std::size_t global = dofs.global_node(region, local);
            if (source_occurrences[source]++ != 0U) {
                if (result.source_global_nodes[source] != global)
                    throw std::invalid_argument("Shared source node does not map to one global four-field node");
                continue;
            }
            result.source_global_nodes[source] = global;
            const auto& actual_point = mesh.nodes()[source];
            const auto& expected = reference[source];
            result.maximum_coordinate_difference = std::max(result.maximum_coordinate_difference,
                std::max({std::abs(actual_point.x - expected.point.x), std::abs(actual_point.y - expected.point.y),
                    std::abs(actual_point.z - expected.point.z)}));
            result.fields[0].add(state[dofs.dof(fuelsim::Field::temperature, global)], expected.fields[0]);
            result.fields[1].add(state[dofs.dof(fuelsim::Field::displacement_x, global)], expected.fields[1]);
            result.fields[2].add(state[dofs.dof(fuelsim::Field::displacement_y, global)], expected.fields[2]);
            result.fields[3].add(state[dofs.dof(fuelsim::Field::displacement_z, global)], expected.fields[3]);
        }
    }
    if (std::find(source_occurrences.begin(), source_occurrences.end(), 0U) != source_occurrences.end())
        throw std::invalid_argument("Three-dimensional comparison did not visit every Exodus source node");
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        const bool is_interface = std::abs(mesh.nodes()[source].x - 1.0) < 1.0e-12;
        const std::size_t expected_occurrences = is_interface ? 2U : 1U;
        if (source_occurrences[source] != expected_occurrences)
            throw std::invalid_argument("The two-block shared-node interface has an unexpected node ownership count");
    }
    if (dofs.node_count() != mesh.nodes().size())
        throw std::invalid_argument("Shared-node blocks did not produce exactly one global node per Exodus node");
    return result;
}

double average_temperature_at_x(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::spatial_detail::SpatialLayout& dofs, const std::vector<std::size_t>& source_global_nodes,
    const std::vector<double>& state, double x) {
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        if (std::abs(mesh.nodes()[source].x - x) >= 1.0e-12) continue;
        sum += state[dofs.dof(fuelsim::Field::temperature, source_global_nodes[source])];
        ++count;
    }
    if (count == 0U) throw std::invalid_argument("Expected temperature plane does not exist in shared-node mesh");
    return sum / static_cast<double>(count);
}

bool check_interface_balance(const fuelsim::UnstructuredHex8Mesh& mesh, const fuelsim::SteadyProblem& problem,
    const std::vector<double>& state, const std::vector<std::size_t>& source_global_nodes) {
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    std::vector<std::size_t> interface_dofs;
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        if (std::abs(mesh.nodes()[source].x - 1.0) >= 1.0e-12) continue;
        for (const fuelsim::Field field : {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
                 fuelsim::Field::displacement_y, fuelsim::Field::displacement_z})
            interface_dofs.push_back(dofs.dof(field, source_global_nodes[source]));
    }
    std::vector<double> balance(interface_dofs.size(), 0.0), scale(interface_dofs.size(), 0.0);
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            for (std::size_t item = 0; item < interface_dofs.size(); ++item)
                if (workspace.dofs[local] == interface_dofs[item]) {
                    balance[item] += workspace.residual[local];
                    scale[item] += std::abs(workspace.residual[local]);
                }
    }
    bool passed = true;
    for (std::size_t item = 0; item < interface_dofs.size(); ++item) {
        const double relative = std::abs(balance[item]) / std::max(scale[item], 1.0);
        std::cout << "b35_interface_residual_" << item << '=' << relative << '\n';
        passed =
            check(relative < 1.0e-8, "the shared interface has balanced assembled thermal or mechanical residuals") &&
            passed;
    }
    const double left = average_temperature_at_x(mesh, dofs, source_global_nodes, state, 0.0);
    const double interface = average_temperature_at_x(mesh, dofs, source_global_nodes, state, 1.0);
    const double right = average_temperature_at_x(mesh, dofs, source_global_nodes, state, 2.0);
    const double meat_heat_flux = 10.0 * (interface - left);
    const double clad_heat_flux = 20.0 * (right - interface);
    const double heat_flux_difference =
        std::abs(meat_heat_flux - clad_heat_flux) / std::max({std::abs(meat_heat_flux), std::abs(clad_heat_flux), 1.0});
    std::cout << "b35_interface_temperature=" << interface << " b35_heat_flux=" << meat_heat_flux
              << " b35_heat_flux_relative_difference=" << heat_flux_difference << '\n';
    return check(std::abs(interface - 500.0) < 1.0e-9,
               "the shared material interface has the analytic 500 K temperature") &&
           check(heat_flux_difference < 1.0e-10, "the two material blocks have continuous interface heat flux") &&
           passed;
}

bool run(const std::string& input_path, const std::string& nodal_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d)
        throw std::invalid_argument("B3.5 comparison requires a steady Cartesian three-dimensional input card");
    if (definition.spatial.contacts.size() != 0U || definition.spatial.regions.size() != 2U)
        throw std::invalid_argument("B3.5 comparison requires two conforming blocks without contact");
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
    bool passed = check(solve.completed && solve.solve.converged, "B3.5 shared-node input-card solve converges");
    const auto comparison =
        compare_nodes(mesh, problem, solve.solve.state, read_nodes(nodal_path, mesh.nodes().size()));
    const std::array<std::string, 4> field_names = {
        "temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t field = 0; field < comparison.fields.size(); ++field) {
        fuelsim::test::print_relative_metrics("b35_" + field_names[field], comparison.fields[field]);
        if (field == 0U)
            passed = check(fuelsim::test::relative_metrics_below(comparison.fields[field], 1.0e-3),
                         "B3.5 temperature three metrics are below 0.1 percent") &&
                     passed;
        passed = check(comparison.fields[field].maximum_zero_reference_difference < 1.0e-10,
                     "B3.5 " + field_names[field] + " zero-reference values remain absolutely bounded") &&
                 passed;
    }
    std::cout << "b35_moose_displacement_formulation_difference=diagnostic_only\n";
    return check(comparison.maximum_coordinate_difference < 1.0e-12,
               "B3.5 compares every Exodus node at matching coordinates") &&
           check_interface_balance(mesh, problem, solve.solve.state, comparison.source_global_nodes) && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_b35_hex8_shared_nodes_moose_tests <case.fsi> <nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B3.5 shared-node MOOSE comparison\n");
        if (!run(argv[1], argv[2])) return 1;
        std::cout << "[PASS] B3.5 shared-node MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B3.5 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
