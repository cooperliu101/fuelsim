#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool test_time_event_alignment(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblemDefinition definition =
        input.transient_definition();
    definition.spatial.time_tables.emplace_back(
        "events", std::vector<double>{0.0, 0.75, 2.0},
        std::vector<double>{1.0, 1.0, 1.0});
    fuelsim::TransientProblem problem(std::move(definition), mesh);
    const fuelsim::TransientTimeOptions time_options = {2.0, 2.0, 0.125, 2.0,
                                                        1.0, 0.5, 2,     20.0};
    const fuelsim::SolverOptions solver_options = {
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options, solver_options);
    return check(result.completed && result.accepted_steps.size() == 2 &&
                     result.accepted_steps[0].time == 0.75 &&
                     result.accepted_steps[1].time == 2.0,
                 "time stepping lands exactly on table events without losing "
                 "the controller step");
}

bool test_moose_time_table_convection(const std::string& input_path,
                                      const std::string& nodal_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.transient_definition(), mesh);
    const fuelsim::TransientTimeOptions time_options = {
        definition.transient_execution.end_time,
        definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step,
        definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor,
        definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks,
        definition.transient_execution.load_ramp_time};
    const fuelsim::SolverOptions solver_options = {
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance, definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options, solver_options);
    bool passed =
        check(result.completed && result.accepted_steps.size() == 4 &&
                  result.accepted_steps[0].time == 2.5 &&
                  result.accepted_steps[1].time == 5.0 &&
                  result.accepted_steps[2].time == 8.0 &&
                  result.accepted_steps[3].time == 10.0,
              "M3.1 lands on power-table events and preserves nominal dt") &&
        check(problem.contribution_count() == 10,
              "M3.1 adds two local convection contributions to eight Quad4s");
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(
            problem, result.committed_state, reference);
    passed = check(fields.node_count == mesh.nodes().size() &&
                       fields.maximum_coordinate_difference < 1.0e-12,
                   "M3.1 compares every MOOSE convection node") &&
             passed;
    passed =
        check(fuelsim::test::relative_metrics_below(fields.temperature, 1.0e-3),
              "M3.1 convection temperature full-field errors pass") &&
        passed;
    passed = check(fuelsim::test::absolute_metrics_below(
                       fields.radial_displacement, 1.0e-12) &&
                       fuelsim::test::absolute_metrics_below(
                           fields.axial_displacement, 1.0e-12),
                   "M3.1 zero displacement fields pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m31_convection_temperature",
                                          fields.temperature);
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m3_load_boundary_tests <m21.fsi> "
                     "<m31.fsi> <m31-all-nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M3.1 time loads and boundary test\n");
        if (!test_time_event_alignment(argv[1]) ||
            !test_moose_time_table_convection(argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] fuelsim M3.1 time loads and boundary test\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.1 test raised: " << error.what() << '\n';
        return 1;
    }
}
