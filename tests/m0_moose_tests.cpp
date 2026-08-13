#include "fuelsim/case_input.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"
#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
bool run_comparison(const std::string& input_path, const std::string& nodal_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument("M0 comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::read_exodus_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
    const fuelsim::SolverOptions options = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        options);
    bool passed = check(result.completed && result.solve.converged, "M0 input-card solve converged");
    passed = check(fuelsim::rz::ProblemAccess::region_count(problem) == 1 &&
                       fuelsim::rz::ProblemAccess::region_mesh(problem, 0).elements().size() == 400,
                 "M0 input card selects the complete MOOSE block") &&
             passed;
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, result.solve.state, reference);
    constexpr double tolerance = 1.0e-10;
    passed = check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                 "M0 compares every MOOSE node at matching coordinates") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.temperature, tolerance),
                 "M0 full-field temperature three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.radial_displacement, tolerance),
                 "M0 full-field radial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.axial_displacement, tolerance),
                 "M0 full-field axial displacement three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m0_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m0_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m0_axial_displacement", fields.axial_displacement);
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_m0_moose_tests <m0.fsi> "
                     "<all-nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim input-card M0 MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2])) return 1;
        std::cout << "[PASS] input-card M0 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M0 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
