#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool run(const std::string& input_path, const std::string& reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady ||
        definition.geometry != fuelsim::CaseGeometry::axisymmetric_rz)
        throw std::invalid_argument("RZ shared-node comparison requires a steady axisymmetric input card");
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SolverOptions options = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks_per_step, definition.steady_execution.minimum_load_increment},
        options);
    bool passed = check(solve.completed && solve.solve.converged, "RZ shared-node input-card solve converges");
    const std::filesystem::path results = std::filesystem::temp_directory_path() / "fuelsim_rz_shared_nodes_results.e";
    try {
        fuelsim::write_steady_results(results.string(), mesh, problem, solve.solve.state);
        std::remove(results.c_str());
    } catch (const std::exception& error) {
        (void)std::remove(results.c_str());
        passed = check(false, std::string("RZ shared-node Exodus output succeeds: ") + error.what()) && passed;
    }
    const auto& dofs = fuelsim::rz::ProblemAccess::dof_map(problem);
    passed = check(dofs.node_count() == mesh.nodes().size() && problem.dof_count() == 3 * mesh.nodes().size(),
                 "six Exodus nodes form six global RZ temperature-displacement tuples") &&
             passed;
    for (const std::size_t source : {1U, 2U}) {
        std::size_t occurrences = 0;
        std::size_t global = 0;
        for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
            const auto& region_mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, region);
            const auto found =
                std::find(region_mesh.source_node_ids().begin(), region_mesh.source_node_ids().end(), source);
            if (found == region_mesh.source_node_ids().end()) continue;
            const std::size_t local = static_cast<std::size_t>(found - region_mesh.source_node_ids().begin());
            const std::size_t candidate = dofs.global_node(region, local);
            if (occurrences++ == 0U)
                global = candidate;
            else
                passed = check(global == candidate, "each shared Exodus interface node maps to one RZ node") && passed;
        }
        passed = check(occurrences == 2U, "each conforming interface node belongs to both material blocks") && passed;
    }
    const auto comparison = fuelsim::test::compare_moose_nodal_fields(
        problem, solve.solve.state, fuelsim::test::read_moose_nodal_reference(reference_path));
    constexpr double tolerance = 1.0e-3;
    for (const auto& field :
        {std::pair<std::string, const fuelsim::test::FieldErrorMetrics>{"temperature", comparison.temperature},
            {"radial_displacement", comparison.radial_displacement},
            {"axial_displacement", comparison.axial_displacement}}) {
        fuelsim::test::print_relative_metrics("rz_shared_" + field.first, field.second);
        passed = check(fuelsim::test::relative_metrics_below(field.second, tolerance) &&
                           field.second.maximum_zero_reference_difference < 1.0e-10,
                     "RZ shared-node MOOSE field metrics are below 0.1 percent") &&
                 passed;
    }
    return check(comparison.node_count == mesh.nodes().size() && comparison.maximum_coordinate_difference < 1.0e-12,
               "RZ shared-node comparison covers every Exodus node at matching coordinates") &&
           passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_rz_shared_nodes_moose_tests <case.fsi> <all-nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim RZ shared-node MOOSE comparison\n");
        if (!run(argv[1], argv[2])) return 1;
        std::cout << "[PASS] RZ shared-node MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] RZ shared-node comparison raised: " << error.what() << '\n';
        return 1;
    }
}
