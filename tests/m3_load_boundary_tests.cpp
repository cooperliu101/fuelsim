#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <exception>
#include <cmath>
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
    const fuelsim::TransientTimeOptions time_options = {
        2.0, 0.5, 0.125, 2.0, 2.0, 0.5, 2, 20.0, 4, 1};
    const fuelsim::SolverOptions solver_options = {
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options, solver_options);
    return check(result.completed && result.accepted_steps.size() == 3 &&
                     result.accepted_steps[0].time == 0.5 &&
                     result.accepted_steps[1].time == 0.75 &&
                     result.accepted_steps[2].time == 2.0 &&
                     result.accepted_steps[2].time_step == 1.25,
                 "iteration-adaptive stepping grows, lands on an event, and "
                 "preserves the controller step");
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

bool test_failure_diagnostics(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    const fuelsim::TransientTimeOptions time_options = {1.0, 1.0, 0.125, 1.0,
                                                        1.0, 0.5, 0,     20.0};
    const fuelsim::SolverOptions solver_options = {
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, 1};
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options, solver_options);
    bool passed =
        check(!result.completed && result.rejected_steps.size() == 1 &&
                  result.termination_reason ==
                      fuelsim::TransientTerminationReason::maximum_cutbacks &&
                  result.rejected_steps[0].time_step == 1.0 &&
                  result.committed_time == 0.0,
              "failed nonlinear step is categorized and leaves committed time "
              "unchanged");
    fuelsim::TransientProblem minimum_problem(input.transient_definition(),
                                              mesh);
    const fuelsim::TransientTimeOptions minimum_options = {
        1.0, 1.0, 0.75, 1.0, 1.0, 0.5, 3, 20.0};
    const fuelsim::TransientResult minimum_result = fuelsim::solve_transient(
        minimum_problem, minimum_options, solver_options);
    passed =
        check(!minimum_result.completed &&
                  minimum_result.termination_reason ==
                      fuelsim::TransientTerminationReason::minimum_time_step &&
                  minimum_result.rejected_steps.size() == 2 &&
                  minimum_result.rejected_steps.back().time_step == 0.75 &&
                  minimum_result.committed_time == 0.0,
              "failed step attempts dt_min once before reporting the limit") &&
        passed;
    return passed;
}

bool test_steady_load_cutback(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::SteadyProblem problem(input.steady_definition(), mesh);
    fuelsim::SolverOptions solver_options{
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, 30};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem, {1, 0.5, 12, 1.0e-4}, solver_options);
    std::cout << "steady_cutback_rejected_steps="
              << result.rejected_steps.size() << '\n';
    std::cout << "steady_cutback_total_cutbacks=" << result.total_cutbacks
              << '\n';
    if (!result.rejected_steps.empty()) {
        std::cout << "steady_cutback_first_attempted_load="
                  << result.rejected_steps.front().attempted_load_factor
                  << '\n';
        std::cout << "steady_cutback_last_attempted_load="
                  << result.rejected_steps.back().attempted_load_factor
                  << '\n';
        std::cout << "steady_cutback_last_failure="
                  << fuelsim::solve_failure_category_name(
                         result.rejected_steps.back().failure_category)
                  << '\n';
        std::cout << "steady_cutback_last_message="
                  << result.rejected_steps.back().failure_message << '\n';
    }
    return check(result.completed && result.solve.converged &&
                     result.total_cutbacks > 0 &&
                     !result.rejected_steps.empty(),
                 "steady loading bisects a failed nominal increment and "
                 "continues from the accepted state");
}

bool test_pressure_production_path(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::SteadyProblem problem(input.steady_definition(), mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem,
        {input.steady_execution.load_steps,
         input.steady_execution.cutback_factor,
         input.steady_execution.maximum_cutbacks,
         input.steady_execution.minimum_load_increment},
        {input.solver.absolute_tolerance, input.solver.relative_tolerance,
         input.solver.step_tolerance, input.solver.maximum_iterations});
    if (!check(result.completed && result.solve.converged,
               "pressure input-card production solve converges"))
        return false;

    double radial_stress_sum = 0.0;
    double hoop_stress_sum = 0.0;
    std::size_t stress_points = 0;
    for (std::size_t element = 0;
         element < problem.region_element_count(0); ++element) {
        const fuelsim::LocalValues local = problem.contribution_state(
            problem.region_element_offset(0) + element, result.solve.state);
        for (const fuelsim::AxisymmetricStressValues& stress :
             problem.region_kernel(0).stress_values(
                 problem.region_element_geometry(0, element), local)) {
            radial_stress_sum += stress.rr;
            hoop_stress_sum += stress.hoop;
            ++stress_points;
        }
    }
    const double average_radial_stress =
        radial_stress_sum / static_cast<double>(stress_points);
    const double average_hoop_stress =
        hoop_stress_sum / static_cast<double>(stress_points);
    const double relative_error = std::max(
        std::abs(average_radial_stress + 1.0e6) / 1.0e6,
        std::abs(average_hoop_stress + 1.0e6) / 1.0e6);
    std::cout << "pressure_average_radial_stress=" << average_radial_stress
              << '\n';
    std::cout << "pressure_average_hoop_stress=" << average_hoop_stress
              << '\n';
    std::cout << "pressure_cylinder_stress_relative_error="
              << relative_error << '\n';
    return check(relative_error < 1.0e-10,
                 "radial pressure input produces the solid-cylinder stress");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: fuelsim_m3_load_boundary_tests <m21.fsi> "
                     "<m31.fsi> <m31-all-nodes.csv> <pcmi.fsi> "
                     "<pressure.fsi>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M3.1 time loads and boundary test\n");
        if (!test_time_event_alignment(argv[1]) ||
            !test_moose_time_table_convection(argv[2], argv[3]) ||
            !test_failure_diagnostics(argv[4]) ||
            !test_steady_load_cutback(argv[4]) ||
            !test_pressure_production_path(argv[5]))
            return 1;
        std::cout << "[PASS] fuelsim M3.1 time loads and boundary test\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.1 test raised: " << error.what() << '\n';
        return 1;
    }
}
