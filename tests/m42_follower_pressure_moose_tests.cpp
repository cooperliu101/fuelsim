#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"
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
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument("M4.2 follower-pressure comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        {definition.solver.absolute_tolerance, definition.solver.relative_tolerance, definition.solver.step_tolerance,
            definition.solver.maximum_iterations});
    bool passed = check(result.completed && result.solve.converged, "M4.2 follower-pressure solve converged");
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, result.solve.state, reference);
    constexpr double tolerance = 5.0e-3;
    passed = check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                 "M4.2 compares every MOOSE node at matching coordinates") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.temperature, tolerance),
                 "M4.2 temperature three full-field errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.radial_displacement, tolerance),
                 "M4.2 radial-displacement three full-field errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.axial_displacement, tolerance),
                 "M4.2 axial-displacement three full-field errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m42_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m42_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m42_axial_displacement", fields.axial_displacement);
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_m42_follower_pressure_moose_tests "
                     "<outer-case.fsi> <outer-all-nodes.csv> "
                     "<left-top-case.fsi> <left-top-all-nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M4.2 follower-pressure MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2]) || !run_comparison(argv[3], argv[4])) return 1;
        std::cout << "[PASS] M4.2 follower-pressure MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M4.2 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
