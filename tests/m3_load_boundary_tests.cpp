#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <array>
#include <exception>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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

double temperature_relative_l2(const fuelsim::TransientProblem& problem,
                               const std::vector<double>& actual,
                               const std::vector<double>& reference) {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    const std::size_t node_count = problem.dof_map().node_count();
    for (std::size_t node = 0; node < node_count; ++node) {
        const std::size_t dof = problem.dof_map().temperature(node);
        const double difference = actual[dof] - reference[dof];
        difference_squared += difference * difference;
        reference_squared += reference[dof] * reference[dof];
    }
    return std::sqrt(difference_squared / reference_squared);
}

struct ConvergenceAccumulator final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_absolute_difference = 0.0;
};

struct ConvergenceMetric final {
    double absolute_l2 = 0.0;
    double relative_l2 = 0.0;
    double maximum_absolute_difference = 0.0;
    bool zero_reference = false;
};

void accumulate_convergence(ConvergenceAccumulator& accumulator,
                            double actual, double reference) {
    const double difference = actual - reference;
    accumulator.difference_squared += difference * difference;
    accumulator.reference_squared += reference * reference;
    accumulator.maximum_absolute_difference =
        std::max(accumulator.maximum_absolute_difference,
                 std::abs(difference));
}

ConvergenceMetric finish_convergence(
    const ConvergenceAccumulator& accumulator) {
    ConvergenceMetric result;
    result.absolute_l2 = std::sqrt(accumulator.difference_squared);
    result.maximum_absolute_difference =
        accumulator.maximum_absolute_difference;
    result.zero_reference = accumulator.reference_squared == 0.0;
    if (!result.zero_reference)
        result.relative_l2 = std::sqrt(accumulator.difference_squared /
                                       accumulator.reference_squared);
    return result;
}

std::array<ConvergenceMetric, 9> compare_committed_states(
    const fuelsim::TransientCommittedState& actual,
    const fuelsim::TransientCommittedState& reference) {
    if (actual.solution.size() != reference.solution.size() ||
        actual.solution.size() % 3 != 0 ||
        actual.material_histories.size() !=
            reference.material_histories.size() ||
        actual.material_stresses.size() !=
            reference.material_stresses.size())
        throw std::logic_error(
            "time-convergence committed-state layouts differ");

    std::array<ConvergenceAccumulator, 9> accumulators{};
    const std::size_t node_count = actual.solution.size() / 3;
    for (std::size_t field = 0; field < 3; ++field) {
        for (std::size_t node = 0; node < node_count; ++node) {
            const std::size_t dof = field * node_count + node;
            accumulate_convergence(accumulators[field], actual.solution[dof],
                                   reference.solution[dof]);
        }
    }

    for (std::size_t region = 0;
         region < actual.material_histories.size(); ++region) {
        if (actual.material_histories[region].size() !=
                reference.material_histories[region].size() ||
            actual.material_stresses[region].size() !=
                reference.material_stresses[region].size())
            throw std::logic_error(
                "time-convergence material-state region layouts differ");
        for (std::size_t element = 0;
             element < actual.material_histories[region].size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const fuelsim::MaterialPointState& actual_history =
                    actual.material_histories[region][element][q];
                const fuelsim::MaterialPointState& reference_history =
                    reference.material_histories[region][element][q];
                const fuelsim::AxisymmetricStressValues& actual_stress =
                    actual.material_stresses[region][element][q];
                const fuelsim::AxisymmetricStressValues& reference_stress =
                    reference.material_stresses[region][element][q];
                const std::array<double, 4> actual_stress_values = {
                    actual_stress.rr, actual_stress.zz, actual_stress.hoop,
                    actual_stress.rz};
                const std::array<double, 4> reference_stress_values = {
                    reference_stress.rr, reference_stress.zz,
                    reference_stress.hoop, reference_stress.rz};
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_convergence(
                        accumulators[3], actual_stress_values[component],
                        reference_stress_values[component]);
                    accumulate_convergence(
                        accumulators[4],
                        actual_history.elastic_strain[component],
                        reference_history.elastic_strain[component]);
                    accumulate_convergence(
                        accumulators[5],
                        actual_history.plastic_strain[component],
                        reference_history.plastic_strain[component]);
                    accumulate_convergence(
                        accumulators[6],
                        actual_history.creep_strain[component],
                        reference_history.creep_strain[component]);
                }
                accumulate_convergence(
                    accumulators[7],
                    actual_history.equivalent_plastic_strain,
                    reference_history.equivalent_plastic_strain);
                accumulate_convergence(
                    accumulators[8],
                    actual_history.equivalent_creep_strain,
                    reference_history.equivalent_creep_strain);
            }
        }
    }

    std::array<ConvergenceMetric, 9> result{};
    for (std::size_t field = 0; field < result.size(); ++field)
        result[field] = finish_convergence(accumulators[field]);
    return result;
}

fuelsim::TransientCommittedState solve_fixed_pcmi(
    const fuelsim::FuelSimCaseDefinition& input,
    const fuelsim::UnstructuredQuad4Mesh& mesh, double time_step) {
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    fuelsim::SolverOptions solver_options = {
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    solver_options.temperature_residual_scale = 1.0e4;
    solver_options.mechanical_residual_scale = 1.0e3;
    const fuelsim::TransientTimeOptions time_options = {
        100.0, time_step, time_step, time_step, 1.0, 0.5, 0, 20.0};
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options, solver_options);
    if (!result.completed || result.aggregate_timing.workspace_setups != 1)
        throw std::runtime_error(
            "fixed-step PCMI time-convergence solve did not complete with "
            "one PETSc workspace");
    return problem.committed_state();
}

bool test_long_transient_time_convergence(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    const fuelsim::TransientCommittedState coarse =
        solve_fixed_pcmi(input, mesh, 1.0);
    const fuelsim::TransientCommittedState medium =
        solve_fixed_pcmi(input, mesh, 0.5);
    const fuelsim::TransientCommittedState fine =
        solve_fixed_pcmi(input, mesh, 0.25);
    const fuelsim::TransientCommittedState reference =
        solve_fixed_pcmi(input, mesh, 0.125);
    const std::array<ConvergenceMetric, 9> coarse_error =
        compare_committed_states(coarse, reference);
    const std::array<ConvergenceMetric, 9> medium_error =
        compare_committed_states(medium, reference);
    const std::array<ConvergenceMetric, 9> fine_error =
        compare_committed_states(fine, reference);
    const std::array<ConvergenceMetric, 9> coarse_to_medium_difference =
        compare_committed_states(coarse, medium);
    const std::array<ConvergenceMetric, 9> medium_to_fine_difference =
        compare_committed_states(medium, fine);
    const std::array<ConvergenceMetric, 9> fine_to_reference_difference =
        compare_committed_states(fine, reference);
    const std::array<const char*, 9> names = {
        "temperature", "radial_displacement", "axial_displacement",
        "stress", "elastic_strain", "plastic_strain", "creep_strain",
        "equivalent_plastic_strain", "equivalent_creep_strain"};

    bool passed = true;
    double minimum_coarse_to_medium_order =
        std::numeric_limits<double>::infinity();
    double minimum_medium_to_fine_order =
        std::numeric_limits<double>::infinity();
    std::size_t rate_evidence_fields = 0;
    std::size_t first_order_trend_fields = 0;
    for (std::size_t field = 0; field < names.size(); ++field) {
        std::cout << "long_time_convergence_" << names[field]
                  << "_relative_l2=" << coarse_error[field].relative_l2 << ','
                  << medium_error[field].relative_l2 << ','
                  << fine_error[field].relative_l2 << '\n';
        std::cout << "long_time_convergence_" << names[field]
                  << "_maximum_absolute_difference="
                  << coarse_error[field].maximum_absolute_difference << ','
                  << medium_error[field].maximum_absolute_difference << ','
                  << fine_error[field].maximum_absolute_difference << '\n';
        std::cout << "long_time_convergence_" << names[field]
                  << "_zero_reference="
                  << coarse_error[field].zero_reference << '\n';
        if (coarse_error[field].zero_reference) {
            passed = check(
                         coarse_error[field].maximum_absolute_difference ==
                                 0.0 &&
                             medium_error[field]
                                     .maximum_absolute_difference == 0.0 &&
                             fine_error[field].maximum_absolute_difference ==
                                 0.0,
                         std::string("zero-reference ") + names[field] +
                             " field remains exactly zero under time "
                             "refinement") &&
                     passed;
            continue;
        }
        const bool roundoff_limited_temperature =
            field == 0 &&
            coarse_error[field].maximum_absolute_difference < 1.0e-8;
        std::cout << "long_time_convergence_" << names[field]
                  << "_roundoff_limited=" << roundoff_limited_temperature
                  << '\n';
        if (roundoff_limited_temperature)
            continue;
        ++rate_evidence_fields;
        const double coarse_to_medium_order = std::log2(
            coarse_to_medium_difference[field].absolute_l2 /
            medium_to_fine_difference[field].absolute_l2);
        const double medium_to_fine_order = std::log2(
            medium_to_fine_difference[field].absolute_l2 /
            fine_to_reference_difference[field].absolute_l2);
        minimum_coarse_to_medium_order =
            std::min(minimum_coarse_to_medium_order,
                     coarse_to_medium_order);
        minimum_medium_to_fine_order =
            std::min(minimum_medium_to_fine_order, medium_to_fine_order);
        if (coarse_to_medium_order > 0.75 &&
            medium_to_fine_order > 0.75)
            ++first_order_trend_fields;
        std::cout << "long_time_convergence_" << names[field]
                  << "_observed_orders=" << coarse_to_medium_order << ','
                  << medium_to_fine_order << '\n';
        passed = check(coarse_error[field].relative_l2 >
                               medium_error[field].relative_l2 &&
                           medium_error[field].relative_l2 >
                               fine_error[field].relative_l2,
                       std::string("100-second PCMI ") + names[field] +
                           " error decreases under time-step refinement") &&
                 passed;
    }
    std::cout << "long_time_convergence_rate_evidence_fields="
              << rate_evidence_fields << '\n';
    std::cout << "long_time_convergence_first_order_trend_fields="
              << first_order_trend_fields << '\n';
    std::cout << "long_time_convergence_minimum_observed_orders="
              << minimum_coarse_to_medium_order << ','
              << minimum_medium_to_fine_order << '\n';
    return check(rate_evidence_fields == 8 &&
                     first_order_trend_fields >= 7 &&
                     minimum_coarse_to_medium_order > 0.4 &&
                     minimum_medium_to_fine_order > 0.4,
                 "100-second PCMI nodal, stress, and complete inelastic "
                 "history fields monotonically approach the fine-step "
                 "reference, with the nonsmooth radial contact response "
                 "reported separately") &&
           passed;
}

bool test_time_error_control(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    const fuelsim::SolverOptions solver_options = {
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};

    fuelsim::TransientProblem reference_problem(
        input.transient_definition(), mesh);
    const fuelsim::TransientTimeOptions reference_options = {
        10.0, 0.009765625, 0.009765625, 0.009765625,
        1.0, 0.5, 0, 20.0};
    const fuelsim::TransientResult reference = fuelsim::solve_transient(
        reference_problem, reference_options, solver_options);

    fuelsim::TransientProblem coarse_problem(input.transient_definition(),
                                             mesh);
    const fuelsim::TransientTimeOptions coarse_options = {
        10.0, 2.5, 2.5, 2.5, 1.0, 0.5, 0, 20.0};
    const fuelsim::TransientResult coarse = fuelsim::solve_transient(
        coarse_problem, coarse_options, solver_options);

    fuelsim::TransientProblem adaptive_problem(input.transient_definition(),
                                               mesh);
    fuelsim::TransientTimeOptions adaptive_options = {
        10.0, 2.5, 0.01953125, 2.5, 2.0, 0.5, 20, 20.0};
    adaptive_options.time_error_relative_tolerance = 2.0e-4;
    adaptive_options.temperature_time_absolute_tolerance = 1.0e-3;
    adaptive_options.displacement_time_absolute_tolerance = 1.0e-8;
    const fuelsim::TransientResult adaptive = fuelsim::solve_transient(
        adaptive_problem, adaptive_options, solver_options);

    const double coarse_error = temperature_relative_l2(
        adaptive_problem, coarse.committed_state, reference.committed_state);
    const double adaptive_error = temperature_relative_l2(
        adaptive_problem, adaptive.committed_state,
        reference.committed_state);
    fuelsim::TransientProblem two_half_problem(input.transient_definition(),
                                               mesh);
    const double first_accepted_step =
        adaptive.accepted_steps.empty()
            ? 2.0 * adaptive_options.minimum_time_step
            : adaptive.accepted_steps.front().time_step;
    const fuelsim::TransientTimeOptions two_half_options = {
        first_accepted_step,
        0.5 * first_accepted_step,
        0.5 * first_accepted_step,
        0.5 * first_accepted_step,
        1.0,
        0.5,
        0,
        20.0};
    const fuelsim::TransientResult two_half = fuelsim::solve_transient(
        two_half_problem, two_half_options, solver_options);
    bool full_interval_conservation = false;
    if (!adaptive.accepted_steps.empty() &&
        two_half.accepted_steps.size() == 2) {
        const fuelsim::TransientConservationSummary& actual =
            adaptive.accepted_steps.front().conservation;
        const fuelsim::TransientConservationSummary& first =
            two_half.accepted_steps[0].conservation;
        const fuelsim::TransientConservationSummary& second =
            two_half.accepted_steps[1].conservation;
        const double expected_stored =
            0.5 * (first.stored_heat_rate + second.stored_heat_rate);
        const double expected_convection =
            0.5 * (first.convection_heat_rate +
                   second.convection_heat_rate);
        const double scale = std::max(
            {1.0, std::abs(expected_stored),
             std::abs(expected_convection)});
        full_interval_conservation =
            std::abs(actual.stored_heat_rate - expected_stored) <=
                1.0e-12 * scale &&
            std::abs(actual.convection_heat_rate - expected_convection) <=
                1.0e-12 * scale;
    }
    double maximum_accepted_estimate = 0.0;
    bool accepted_cutback_observed = false;
    bool accepted_cutback_regrew_immediately = false;
    for (const fuelsim::TransientAcceptedStep& step :
         adaptive.accepted_steps) {
        maximum_accepted_estimate =
            std::max(maximum_accepted_estimate, step.time_error_estimate);
        if (step.cutbacks > 0) {
            accepted_cutback_observed = true;
            accepted_cutback_regrew_immediately =
                accepted_cutback_regrew_immediately ||
                step.next_time_step >
                    step.time_step *
                        (1.0 + 16.0 *
                                   std::numeric_limits<double>::epsilon());
        }
    }
    double maximum_rejected_estimate = 0.0;
    double minimum_rejected_estimate =
        std::numeric_limits<double>::infinity();
    for (const fuelsim::TransientRejectedStep& step :
         adaptive.rejected_steps) {
        if (step.failure_category !=
            fuelsim::SolveFailureCategory::time_discretization)
            continue;
        maximum_rejected_estimate =
            std::max(maximum_rejected_estimate, step.time_error_estimate);
        minimum_rejected_estimate =
            std::min(minimum_rejected_estimate, step.time_error_estimate);
    }
    std::cout << "time_control_coarse_temperature_relative_l2="
              << coarse_error << '\n';
    std::cout << "time_control_adaptive_temperature_relative_l2="
              << adaptive_error << '\n';
    std::cout << "time_control_rejections="
              << adaptive.time_error_rejections << '\n';
    std::cout << "time_control_maximum_accepted_estimate="
              << maximum_accepted_estimate << '\n';
    std::cout << "time_control_minimum_rejected_estimate="
              << minimum_rejected_estimate << '\n';
    std::cout << "time_control_maximum_rejected_estimate="
              << maximum_rejected_estimate << '\n';
    std::cout << "time_control_completed=" << adaptive.completed << '\n';
    std::cout << "time_control_accepted_steps="
              << adaptive.accepted_steps.size() << '\n';
    return check(reference.completed && coarse.completed &&
                     adaptive.completed &&
                     adaptive.time_error_rejections > 0 &&
                     accepted_cutback_observed &&
                     !accepted_cutback_regrew_immediately &&
                     full_interval_conservation &&
                     maximum_accepted_estimate <= 1.0 &&
                     adaptive_error < 0.5 * coarse_error &&
                     adaptive.aggregate_timing.workspace_setups == 1,
                 "BE step-doubling rejects inaccurate steps, reuses the "
                 "PETSc workspace, and reduces temporal error");
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

bool test_history_time_error_control(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    fuelsim::TransientTimeOptions time_options = {
        20.0, 4.0, 0.125, 4.0, 1.0, 0.5, 12, 20.0};
    time_options.time_error_relative_tolerance = 5.0e-1;
    time_options.temperature_time_absolute_tolerance = 1.0e-3;
    time_options.displacement_time_absolute_tolerance = 1.0e-9;
    time_options.strain_history_time_absolute_tolerance = 1.0e-10;
    time_options.stress_history_time_absolute_tolerance = 10.0;
    const fuelsim::SolverOptions solver_options = {
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options, solver_options);
    double maximum_nodal = 0.0;
    double maximum_history = 0.0;
    for (const fuelsim::TransientAcceptedStep& step : result.accepted_steps) {
        const fuelsim::TransientTimeErrorEstimate& error =
            step.time_error_components;
        maximum_nodal = std::max(
            maximum_nodal,
            std::max({error.temperature, error.radial_displacement,
                      error.axial_displacement}));
        maximum_history = std::max(
            maximum_history,
            std::max({error.elastic_strain, error.plastic_strain,
                      error.creep_strain, error.equivalent_plastic_strain,
                      error.equivalent_creep_strain, error.stress}));
    }
    double maximum_rejected_history = 0.0;
    bool history_only_rejection = false;
    for (const fuelsim::TransientRejectedStep& step : result.rejected_steps) {
        const fuelsim::TransientTimeErrorEstimate& error =
            step.time_error_components;
        const double rejected_nodal =
            std::max({error.temperature, error.radial_displacement,
                      error.axial_displacement});
        const double rejected_history =
            std::max({error.elastic_strain, error.plastic_strain,
                      error.creep_strain, error.equivalent_plastic_strain,
                      error.equivalent_creep_strain, error.stress});
        maximum_rejected_history =
            std::max(maximum_rejected_history, rejected_history);
        history_only_rejection = history_only_rejection ||
                                 (rejected_nodal <= 1.0 &&
                                  rejected_history > 1.0);
    }
    std::cout << "history_time_error_maximum_nodal=" << maximum_nodal << '\n';
    std::cout << "history_time_error_maximum_material=" << maximum_history
              << '\n';
    std::cout << "history_time_error_rejections="
              << result.time_error_rejections << '\n';
    std::cout << "history_time_error_maximum_rejected_material="
              << maximum_rejected_history << '\n';
    std::cout << "history_time_error_completed=" << result.completed << '\n';
    std::cout << "history_time_error_accepted_steps="
              << result.accepted_steps.size() << '\n';
    return check(result.completed && result.time_error_rejections > 0 &&
                     history_only_rejection &&
                     maximum_history > maximum_nodal &&
                     maximum_history <= 1.0,
                 "step-doubling controls committed inelastic histories in "
                 "addition to nodal fields");
}

bool test_long_transient_diagnostics(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    const fuelsim::TransientTimeOptions time_options = {
        100.0, 0.5, 0.5, 0.5, 1.0, 0.5, 0, 20.0};
    fuelsim::SolverOptions solver_options;
    solver_options.absolute_tolerance = input.solver.absolute_tolerance;
    solver_options.relative_tolerance = input.solver.relative_tolerance;
    solver_options.step_tolerance = input.solver.step_tolerance;
    solver_options.maximum_iterations = input.solver.maximum_iterations;
    solver_options.temperature_residual_scale = 1.0e4;
    solver_options.mechanical_residual_scale = 1.0e3;
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options, solver_options);
    double maximum_thermal_balance = 0.0;
    double maximum_mechanical_balance = 0.0;
    double maximum_absolute_mechanical_balance = 0.0;
    double accumulated_absolute_mechanical_balance = 0.0;
    double accumulated_mechanical_scale = 0.0;
    double maximum_interface_imbalance = 0.0;
    double minimum_plastic_dissipation =
        std::numeric_limits<double>::infinity();
    double minimum_creep_dissipation =
        std::numeric_limits<double>::infinity();
    double accumulated_plastic_dissipation = 0.0;
    double accumulated_creep_dissipation = 0.0;
    for (const fuelsim::TransientAcceptedStep& step : result.accepted_steps) {
        const fuelsim::TransientConservationSummary& summary =
            step.conservation;
        maximum_thermal_balance =
            std::max(maximum_thermal_balance,
                     summary.relative_thermal_balance);
        maximum_mechanical_balance =
            std::max(maximum_mechanical_balance,
                     summary.relative_mechanical_work_balance);
        maximum_absolute_mechanical_balance =
            std::max(maximum_absolute_mechanical_balance,
                     std::abs(summary.mechanical_work_balance));
        accumulated_absolute_mechanical_balance +=
            std::abs(summary.mechanical_work_balance);
        accumulated_mechanical_scale +=
            std::abs(summary.internal_mechanical_work_increment) +
            std::abs(summary.contact_work_increment) +
            std::abs(summary.pressure_traction_work_increment) +
            std::abs(summary.dirichlet_reaction_work_increment);
        maximum_interface_imbalance =
            std::max(maximum_interface_imbalance,
                     std::abs(summary.interface_heat_imbalance));
        minimum_plastic_dissipation =
            std::min(minimum_plastic_dissipation,
                     summary.plastic_dissipation_increment);
        minimum_creep_dissipation =
            std::min(minimum_creep_dissipation,
                     summary.creep_dissipation_increment);
        accumulated_plastic_dissipation +=
            summary.plastic_dissipation_increment;
        accumulated_creep_dissipation +=
            summary.creep_dissipation_increment;
    }
    std::cout << "long_transient_steps=" << result.accepted_steps.size()
              << '\n';
    std::cout << "long_transient_maximum_relative_thermal_balance="
              << maximum_thermal_balance << '\n';
    std::cout << "long_transient_maximum_relative_mechanical_balance="
              << maximum_mechanical_balance << '\n';
    std::cout << "long_transient_maximum_absolute_mechanical_balance="
              << maximum_absolute_mechanical_balance << '\n';
    std::cout << "long_transient_accumulated_relative_mechanical_balance="
              << accumulated_absolute_mechanical_balance /
                     accumulated_mechanical_scale
              << '\n';
    std::cout << "long_transient_maximum_interface_heat_imbalance="
              << maximum_interface_imbalance << '\n';
    std::cout << "long_transient_minimum_plastic_dissipation="
              << minimum_plastic_dissipation << '\n';
    std::cout << "long_transient_minimum_creep_dissipation="
              << minimum_creep_dissipation << '\n';
    std::cout << "long_transient_accumulated_plastic_dissipation="
              << accumulated_plastic_dissipation << '\n';
    std::cout << "long_transient_accumulated_creep_dissipation="
              << accumulated_creep_dissipation << '\n';
    return check(
        result.completed && result.accepted_steps.size() == 200 &&
            result.aggregate_timing.workspace_setups == 1 &&
            result.last_attempt.field_residual_scalings[0] == 1.0e-4 &&
            result.last_attempt.field_residual_scalings[1] == 1.0e-3 &&
            result.last_attempt.field_residual_scalings[2] == 1.0e-3 &&
            maximum_thermal_balance < 1.0e-8 &&
            accumulated_absolute_mechanical_balance /
                    accumulated_mechanical_scale <
                1.0e-6 &&
            maximum_absolute_mechanical_balance < 1.0e-12 &&
            maximum_interface_imbalance < 1.0e-8 &&
            minimum_plastic_dissipation >= -1.0e-12 &&
            minimum_creep_dissipation >= -1.0e-12 &&
            accumulated_plastic_dissipation > 0.0 &&
            accumulated_creep_dissipation > 0.0,
        "200-step PCMI preserves the PETSc workspace, fixed physical "
        "residual scales, global balances, and nonnegative dissipation");
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
    solver_options.line_search =
        fuelsim::SolverOptions::LineSearch::backtracking;
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
    bool passed = check(result.completed && result.solve.converged &&
                            result.total_cutbacks > 0 &&
                            !result.rejected_steps.empty(),
                        "steady loading bisects a failed nominal increment "
                        "and continues from the accepted state");

    fuelsim::SteadyProblem minimum_problem(input.steady_definition(), mesh);
    fuelsim::SolverOptions minimum_solver = solver_options;
    minimum_solver.maximum_iterations = 1;
    minimum_solver.line_search =
        fuelsim::SolverOptions::LineSearch::basic;
    minimum_solver.backtracking_fallback = false;
    const fuelsim::SteadyResult minimum = fuelsim::solve_steady(
        minimum_problem, {1, 0.5, 3, 0.75}, minimum_solver);
    passed = check(!minimum.completed &&
                       minimum.rejected_steps.size() == 2 &&
                       minimum.total_cutbacks == 1 &&
                       minimum.rejected_steps.back().load_increment == 0.75,
                   "steady loading attempts the minimum load increment once "
                   "before terminating") &&
             passed;
    return passed;
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
            !test_time_error_control(argv[2]) ||
            !test_history_time_error_control(argv[4]) ||
            !test_long_transient_time_convergence(argv[4]) ||
            !test_long_transient_diagnostics(argv[4]) ||
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
