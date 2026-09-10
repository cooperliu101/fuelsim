#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
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
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::SpatialDefinition definition = input.spatial;
    definition.time_tables.emplace_back("events",
        std::vector<double>{0.0, 0.75, 2.0},
        std::vector<double>{1.0, 1.0, 1.0});
    fuelsim::TransientProblem problem(std::move(definition), mesh);
    const fuelsim::TransientTimeOptions time_options = {2.0, 0.5, 0.125, 2.0, 2.0, 0.5, 2, 20.0, 4, 1};
    const fuelsim::SolverOptions solver_options = {input.solver.absolute_tolerance,
        input.solver.relative_tolerance,
        input.solver.step_tolerance,
        input.solver.maximum_iterations};
    const fuelsim::TransientResult result = fuelsim::solve_transient(problem, time_options, solver_options);
    return check(result.completed && result.accepted_steps.size() == 3 && result.accepted_steps[0].time == 0.5
                     && result.accepted_steps[1].time == 0.75 && result.accepted_steps[2].time == 2.0
                     && result.accepted_steps[2].time_step == 1.25,
        "iteration-adaptive stepping grows, lands on an event, and "
        "preserves the controller step");
}

double temperature_relative_l2(const fuelsim::TransientProblem& problem,
    const std::vector<double>& actual,
    const std::vector<double>& reference) {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    const std::size_t node_count = fuelsim::rz::ProblemAccess::dof_map(problem).node_count();
    for (std::size_t node = 0; node < node_count; ++node) {
        const std::size_t dof = fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, node);
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

void accumulate_convergence(ConvergenceAccumulator& accumulator, double actual, double reference) {
    const double difference = actual - reference;
    accumulator.difference_squared += difference * difference;
    accumulator.reference_squared += reference * reference;
    accumulator.maximum_absolute_difference = std::max(accumulator.maximum_absolute_difference, std::abs(difference));
}

ConvergenceMetric finish_convergence(const ConvergenceAccumulator& accumulator) {
    ConvergenceMetric result;
    result.absolute_l2 = std::sqrt(accumulator.difference_squared);
    result.maximum_absolute_difference = accumulator.maximum_absolute_difference;
    result.zero_reference = accumulator.reference_squared == 0.0;
    if (!result.zero_reference)
        result.relative_l2 = std::sqrt(accumulator.difference_squared / accumulator.reference_squared);
    return result;
}

std::array<ConvergenceMetric, 9> compare_committed_states(const fuelsim::TransientCommittedState& actual,
    const fuelsim::TransientCommittedState& reference) {
    if (actual.solution.size() != reference.solution.size() || actual.solution.size() % 3 != 0
        || actual.material_histories.size() != reference.material_histories.size())
        throw std::logic_error("time-convergence committed-state layouts differ");
    std::array<ConvergenceAccumulator, 9> accumulators{};
    const std::size_t node_count = actual.solution.size() / 3;
    for (std::size_t field = 0; field < 3; ++field) {
        for (std::size_t node = 0; node < node_count; ++node) {
            const std::size_t dof = field * node_count + node;
            accumulate_convergence(accumulators[field], actual.solution[dof], reference.solution[dof]);
        }
    }
    for (std::size_t region = 0; region < actual.material_histories.size(); ++region) {
        if (actual.material_histories[region].size() != reference.material_histories[region].size())
            throw std::logic_error("time-convergence material-state region layouts differ");
        for (std::size_t element = 0; element < actual.material_histories[region].size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const fuelsim::MaterialPointState& actual_history = actual.material_histories[region][element][q];
                const fuelsim::MaterialPointState& reference_history = reference.material_histories[region][element][q];
                const fuelsim::AxisymmetricStressValues& actual_stress = actual_history.stress;
                const fuelsim::AxisymmetricStressValues& reference_stress = reference_history.stress;
                const std::array<double, 4> actual_stress_values = {actual_stress.rr,
                    actual_stress.zz,
                    actual_stress.hoop,
                    actual_stress.rz};
                const std::array<double, 4> reference_stress_values = {reference_stress.rr,
                    reference_stress.zz,
                    reference_stress.hoop,
                    reference_stress.rz};
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_convergence(accumulators[3],
                        actual_stress_values[component],
                        reference_stress_values[component]);
                    accumulate_convergence(accumulators[4],
                        actual_history.elastic_strain[component],
                        reference_history.elastic_strain[component]);
                    accumulate_convergence(accumulators[5],
                        actual_history.plastic_strain[component],
                        reference_history.plastic_strain[component]);
                    accumulate_convergence(accumulators[6],
                        actual_history.creep_strain[component],
                        reference_history.creep_strain[component]);
                }
                accumulate_convergence(accumulators[7],
                    actual_history.equivalent_plastic_strain,
                    reference_history.equivalent_plastic_strain);
                accumulate_convergence(accumulators[8],
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

bool test_opaque_state_snapshot(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.spatial, mesh);
    const fuelsim::TransientCommittedState reference = fuelsim::rz::ProblemAccess::committed_state(problem);
    const fuelsim::ProblemStateSnapshot snapshot = problem.capture_state();
    bool empty_rejected = false;
    try {
        problem.restore_state(fuelsim::ProblemStateSnapshot{});
    } catch (const std::invalid_argument&) {
        empty_rejected = true;
    }
    fuelsim::TransientProblem other_problem(input.spatial, mesh);
    bool foreign_snapshot_rejected = false;
    try {
        other_problem.restore_state(snapshot);
    } catch (const std::invalid_argument&) {
        foreign_snapshot_rejected = true;
    }
    fuelsim::TransientCommittedState changed = reference;
    changed.time += 0.5;
    changed.load_factor = 0.5;
    fuelsim::rz::ProblemAccess::restore_committed_state(problem, std::move(changed));
    problem.restore_state(snapshot);
    const fuelsim::TransientCommittedState restored = fuelsim::rz::ProblemAccess::committed_state(problem);
    const std::array<ConvergenceMetric, 9> state_difference = compare_committed_states(restored, reference);
    bool identical = restored.time == reference.time && restored.load_factor == reference.load_factor
                     && restored.contact_histories.size() == reference.contact_histories.size();
    for (const ConvergenceMetric& metric : state_difference)
        identical = identical && metric.maximum_absolute_difference == 0.0;
    for (const fuelsim::TransientConservationField& field : fuelsim::transient_conservation_fields)
        identical = identical && restored.conservation.*field.member == reference.conservation.*field.member;
    if (restored.contact_histories.size() == reference.contact_histories.size()) {
        for (std::size_t contact = 0; contact < restored.contact_histories.size(); ++contact) {
            identical =
                identical && restored.contact_histories[contact].size() == reference.contact_histories[contact].size();
            if (restored.contact_histories[contact].size() != reference.contact_histories[contact].size())
                continue;
            for (std::size_t point = 0; point < restored.contact_histories[contact].size(); ++point) {
                const fuelsim::ContactPointHistory& actual = restored.contact_histories[contact][point];
                const fuelsim::ContactPointHistory& expected = reference.contact_histories[contact][point];
                identical = identical && actual.elastic_tangential_slip == expected.elastic_tangential_slip
                            && actual.total_tangential_slip == expected.total_tangential_slip
                            && actual.normal_multiplier == expected.normal_multiplier
                            && actual.sliding == expected.sliding;
            }
        }
    }
    return check(!snapshot.empty() && empty_rejected && foreign_snapshot_rejected && identical,
        "opaque snapshots reject invalid use and restore the complete RZ committed state exactly");
}

bool test_time_error_control(const std::string& input_path) {
    fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    // Exercise nonzero convection diagnostics on the same cladding boundary.
    auto& coolant = input.spatial.boundary_conditions.back();
    coolant.type = fuelsim::BoundaryConditionType::convection;
    coolant.heat_transfer_coefficient = 1000.0;
    coolant.ambient_temperature = 500.0;
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    const fuelsim::SolverOptions solver_options = {input.solver.absolute_tolerance,
        input.solver.relative_tolerance,
        input.solver.step_tolerance,
        input.solver.maximum_iterations};
    fuelsim::TransientProblem reference_problem(input.spatial, mesh);
    const fuelsim::TransientTimeOptions reference_options = {2.0, 0.03125, 0.03125, 0.03125, 1.0, 0.5, 0, 20.0};
    const fuelsim::TransientResult reference =
        fuelsim::solve_transient(reference_problem, reference_options, solver_options);
    fuelsim::TransientProblem coarse_problem(input.spatial, mesh);
    const fuelsim::TransientTimeOptions coarse_options = {2.0, 2.0, 2.0, 2.0, 1.0, 0.5, 0, 20.0};
    const fuelsim::TransientResult coarse = fuelsim::solve_transient(coarse_problem, coarse_options, solver_options);
    fuelsim::TransientProblem adaptive_problem(input.spatial, mesh);
    fuelsim::TransientTimeOptions adaptive_options = {2.0, 2.0, 0.03125, 2.0, 2.0, 0.5, 20, 20.0};
    adaptive_options.time_error_relative_tolerance = 2.0e-4;
    adaptive_options.temperature_time_absolute_tolerance = 1.0e-3;
    // Isolate thermal error control; inelastic-history control is tested separately.
    adaptive_options.displacement_time_absolute_tolerance = 1.0;
    adaptive_options.strain_history_time_absolute_tolerance = 1.0;
    adaptive_options.stress_history_time_absolute_tolerance = 1.0e12;
    const fuelsim::TransientResult adaptive =
        fuelsim::solve_transient(adaptive_problem, adaptive_options, solver_options);
    const double coarse_error =
        temperature_relative_l2(adaptive_problem, coarse.committed_state, reference.committed_state);
    const double adaptive_error =
        temperature_relative_l2(adaptive_problem, adaptive.committed_state, reference.committed_state);
    fuelsim::TransientProblem two_half_problem(input.spatial, mesh);
    const double first_accepted_step = adaptive.accepted_steps.empty() ? 2.0 * adaptive_options.minimum_time_step
                                                                       : adaptive.accepted_steps.front().time_step;
    const fuelsim::TransientTimeOptions two_half_options = {first_accepted_step,
        0.5 * first_accepted_step,
        0.5 * first_accepted_step,
        0.5 * first_accepted_step,
        1.0,
        0.5,
        0,
        20.0};
    const fuelsim::TransientResult two_half =
        fuelsim::solve_transient(two_half_problem, two_half_options, solver_options);
    bool full_interval_conservation = false;
    if (!adaptive.accepted_steps.empty() && two_half.accepted_steps.size() == 2) {
        const fuelsim::TransientConservationSummary& actual = adaptive.accepted_steps.front().conservation;
        const fuelsim::TransientConservationSummary& first = two_half.accepted_steps[0].conservation;
        const fuelsim::TransientConservationSummary& second = two_half.accepted_steps[1].conservation;
        const double expected_stored = 0.5 * (first.stored_heat_rate + second.stored_heat_rate);
        const double expected_convection = 0.5 * (first.convection_heat_rate + second.convection_heat_rate);
        const double scale = std::max({1.0, std::abs(expected_stored), std::abs(expected_convection)});
        full_interval_conservation = std::abs(expected_convection) > 0.0
                                     && std::abs(actual.stored_heat_rate - expected_stored) <= 1.0e-12 * scale
                                     && std::abs(actual.convection_heat_rate - expected_convection) <= 1.0e-12 * scale;
    }
    double maximum_accepted_estimate = 0.0;
    bool accepted_cutback_observed = false;
    bool accepted_cutback_regrew_immediately = false;
    for (const fuelsim::TransientAcceptedStep& step : adaptive.accepted_steps) {
        maximum_accepted_estimate = std::max(maximum_accepted_estimate, step.time_error_estimate);
        if (step.cutbacks > 0) {
            accepted_cutback_observed = true;
            accepted_cutback_regrew_immediately =
                accepted_cutback_regrew_immediately
                || step.next_time_step > step.time_step * (1.0 + 16.0 * std::numeric_limits<double>::epsilon());
        }
    }
    double maximum_rejected_estimate = 0.0;
    double minimum_rejected_estimate = std::numeric_limits<double>::infinity();
    for (const fuelsim::TransientRejectedStep& step : adaptive.rejected_steps) {
        if (step.failure_category != fuelsim::SolveFailureCategory::time_discretization)
            continue;
        maximum_rejected_estimate = std::max(maximum_rejected_estimate, step.time_error_estimate);
        minimum_rejected_estimate = std::min(minimum_rejected_estimate, step.time_error_estimate);
    }
    std::cout << "time_control_coarse_temperature_relative_l2=" << coarse_error << '\n';
    std::cout << "time_control_adaptive_temperature_relative_l2=" << adaptive_error << '\n';
    std::cout << "time_control_rejections=" << adaptive.time_error_rejections << '\n';
    std::cout << "time_control_maximum_accepted_estimate=" << maximum_accepted_estimate << '\n';
    std::cout << "time_control_minimum_rejected_estimate=" << minimum_rejected_estimate << '\n';
    std::cout << "time_control_maximum_rejected_estimate=" << maximum_rejected_estimate << '\n';
    std::cout << "time_control_completed=" << adaptive.completed << '\n';
    std::cout << "time_control_accepted_steps=" << adaptive.accepted_steps.size() << '\n';
    return check(reference.completed && coarse.completed && adaptive.completed && adaptive.time_error_rejections > 0
                     && accepted_cutback_observed && !accepted_cutback_regrew_immediately && full_interval_conservation
                     && maximum_accepted_estimate <= 1.0 && adaptive_error < 0.5 * coarse_error
                     && adaptive.aggregate_timing.workspace_setups == 1,
        "BE step-doubling rejects inaccurate steps, reuses the "
        "PETSc workspace, and reduces temporal error");
}

bool test_failure_diagnostics(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.spatial, mesh);
    const fuelsim::TransientTimeOptions time_options = {1.0, 1.0, 0.125, 1.0, 1.0, 0.5, 0, 20.0};
    const fuelsim::SolverOptions solver_options = {input.solver.absolute_tolerance,
        input.solver.relative_tolerance,
        input.solver.step_tolerance,
        1};
    const fuelsim::TransientResult result = fuelsim::solve_transient(problem, time_options, solver_options);
    bool passed = check(!result.completed && result.rejected_steps.size() == 1
                            && result.termination_reason == fuelsim::TransientTerminationReason::maximum_cutbacks
                            && result.rejected_steps[0].time_step == 1.0 && result.committed_time == 0.0,
        "failed nonlinear step is categorized and leaves committed time "
        "unchanged");
    fuelsim::TransientProblem minimum_problem(input.spatial, mesh);
    const fuelsim::TransientTimeOptions minimum_options = {1.0, 1.0, 0.75, 1.0, 1.0, 0.5, 3, 20.0};
    const fuelsim::TransientResult minimum_result =
        fuelsim::solve_transient(minimum_problem, minimum_options, solver_options);
    passed =
        check(!minimum_result.completed
                  && minimum_result.termination_reason == fuelsim::TransientTerminationReason::minimum_time_step
                  && minimum_result.rejected_steps.size() == 2 && minimum_result.rejected_steps.back().time_step == 0.75
                  && minimum_result.committed_time == 0.0,
            "failed step attempts dt_min once before reporting the limit")
        && passed;
    return passed;
}

bool test_history_time_error_control(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.spatial, mesh);
    fuelsim::TransientTimeOptions time_options = {20.0, 4.0, 0.125, 4.0, 1.0, 0.5, 12, 20.0};
    time_options.time_error_relative_tolerance = 5.0e-1;
    time_options.temperature_time_absolute_tolerance = 1.0e-3;
    time_options.displacement_time_absolute_tolerance = 1.0e-9;
    time_options.strain_history_time_absolute_tolerance = 1.0e-10;
    time_options.stress_history_time_absolute_tolerance = 10.0;
    const fuelsim::SolverOptions solver_options = {input.solver.absolute_tolerance,
        input.solver.relative_tolerance,
        input.solver.step_tolerance,
        input.solver.maximum_iterations};
    const fuelsim::TransientResult result = fuelsim::solve_transient(problem, time_options, solver_options);
    double maximum_nodal = 0.0;
    double maximum_history = 0.0;
    for (const fuelsim::TransientAcceptedStep& step : result.accepted_steps) {
        const fuelsim::TransientTimeErrorEstimate& error = step.time_error_components;
        for (const fuelsim::TransientFieldTimeError& field : error.nodal_fields)
            maximum_nodal = std::max(maximum_nodal, field.value);
        maximum_history = std::max(maximum_history,
            std::max({error.elastic_strain,
                error.plastic_strain,
                error.creep_strain,
                error.equivalent_plastic_strain,
                error.equivalent_creep_strain,
                error.stress}));
    }
    double maximum_rejected_history = 0.0;
    bool history_only_rejection = false;
    for (const fuelsim::TransientRejectedStep& step : result.rejected_steps) {
        const fuelsim::TransientTimeErrorEstimate& error = step.time_error_components;
        double rejected_nodal = 0.0;
        for (const fuelsim::TransientFieldTimeError& field : error.nodal_fields)
            rejected_nodal = std::max(rejected_nodal, field.value);
        const double rejected_history = std::max({error.elastic_strain,
            error.plastic_strain,
            error.creep_strain,
            error.equivalent_plastic_strain,
            error.equivalent_creep_strain,
            error.stress});
        maximum_rejected_history = std::max(maximum_rejected_history, rejected_history);
        history_only_rejection = history_only_rejection || (rejected_nodal <= 1.0 && rejected_history > 1.0);
    }
    std::cout << "history_time_error_maximum_nodal=" << maximum_nodal << '\n';
    std::cout << "history_time_error_maximum_material=" << maximum_history << '\n';
    std::cout << "history_time_error_rejections=" << result.time_error_rejections << '\n';
    std::cout << "history_time_error_maximum_rejected_material=" << maximum_rejected_history << '\n';
    std::cout << "history_time_error_completed=" << result.completed << '\n';
    std::cout << "history_time_error_accepted_steps=" << result.accepted_steps.size() << '\n';
    return check(result.completed && result.time_error_rejections > 0 && history_only_rejection
                     && maximum_history > maximum_nodal && maximum_history <= 1.0,
        "step-doubling controls committed inelastic histories in "
        "addition to nodal fields");
}

bool test_long_transient_diagnostics(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.spatial, mesh);
    const fuelsim::TransientTimeOptions time_options = {20.0, 0.5, 0.5, 0.5, 1.0, 0.5, 0, 20.0};
    fuelsim::SolverOptions solver_options;
    solver_options.absolute_tolerance = input.solver.absolute_tolerance;
    solver_options.relative_tolerance = input.solver.relative_tolerance;
    solver_options.step_tolerance = input.solver.step_tolerance;
    solver_options.maximum_iterations = input.solver.maximum_iterations;
    solver_options.temperature_residual_scale = 1.0e4;
    solver_options.mechanical_residual_scale = 1.0e3;
    const fuelsim::TransientResult result = fuelsim::solve_transient(problem, time_options, solver_options);
    double maximum_thermal_balance = 0.0;
    double maximum_mechanical_balance = 0.0;
    double maximum_absolute_mechanical_balance = 0.0;
    double accumulated_absolute_mechanical_balance = 0.0;
    double accumulated_mechanical_scale = 0.0;
    double maximum_interface_imbalance = 0.0;
    double minimum_plastic_dissipation = std::numeric_limits<double>::infinity();
    double minimum_creep_dissipation = std::numeric_limits<double>::infinity();
    double accumulated_plastic_dissipation = 0.0;
    double accumulated_creep_dissipation = 0.0;
    for (const fuelsim::TransientAcceptedStep& step : result.accepted_steps) {
        const fuelsim::TransientConservationSummary& summary = step.conservation;
        maximum_thermal_balance = std::max(maximum_thermal_balance, summary.relative_thermal_balance);
        maximum_mechanical_balance = std::max(maximum_mechanical_balance, summary.relative_mechanical_work_balance);
        maximum_absolute_mechanical_balance =
            std::max(maximum_absolute_mechanical_balance, std::abs(summary.mechanical_work_balance));
        accumulated_absolute_mechanical_balance += std::abs(summary.mechanical_work_balance);
        accumulated_mechanical_scale +=
            std::abs(summary.internal_mechanical_work_increment) + std::abs(summary.contact_work_increment)
            + std::abs(summary.pressure_traction_work_increment) + std::abs(summary.dirichlet_reaction_work_increment);
        maximum_interface_imbalance = std::max(maximum_interface_imbalance, std::abs(summary.interface_heat_imbalance));
        minimum_plastic_dissipation = std::min(minimum_plastic_dissipation, summary.plastic_dissipation_increment);
        minimum_creep_dissipation = std::min(minimum_creep_dissipation, summary.creep_dissipation_increment);
        accumulated_plastic_dissipation += summary.plastic_dissipation_increment;
        accumulated_creep_dissipation += summary.creep_dissipation_increment;
    }
    std::cout << "long_transient_steps=" << result.accepted_steps.size() << '\n';
    std::cout << "long_transient_maximum_relative_thermal_balance=" << maximum_thermal_balance << '\n';
    std::cout << "long_transient_maximum_relative_mechanical_balance=" << maximum_mechanical_balance << '\n';
    std::cout << "long_transient_maximum_absolute_mechanical_balance=" << maximum_absolute_mechanical_balance << '\n';
    std::cout << "long_transient_accumulated_relative_mechanical_balance="
              << accumulated_absolute_mechanical_balance / accumulated_mechanical_scale << '\n';
    std::cout << "long_transient_maximum_interface_heat_imbalance=" << maximum_interface_imbalance << '\n';
    std::cout << "long_transient_minimum_plastic_dissipation=" << minimum_plastic_dissipation << '\n';
    std::cout << "long_transient_minimum_creep_dissipation=" << minimum_creep_dissipation << '\n';
    std::cout << "long_transient_accumulated_plastic_dissipation=" << accumulated_plastic_dissipation << '\n';
    std::cout << "long_transient_accumulated_creep_dissipation=" << accumulated_creep_dissipation << '\n';
    return check(result.completed && result.accepted_steps.size() == 40 && result.aggregate_timing.workspace_setups == 1
                     && result.last_attempt.field_residual_scalings[0] == 1.0e-4
                     && result.last_attempt.field_residual_scalings[1] == 1.0e-3
                     && result.last_attempt.field_residual_scalings[2] == 1.0e-3 && maximum_thermal_balance < 1.0e-8
                     && accumulated_absolute_mechanical_balance / accumulated_mechanical_scale < 1.0e-6
                     && maximum_absolute_mechanical_balance < 1.0e-12 && maximum_interface_imbalance < 1.0e-8
                     && minimum_plastic_dissipation >= -1.0e-12 && minimum_creep_dissipation >= -1.0e-12
                     && accumulated_plastic_dissipation > 0.0 && accumulated_creep_dissipation > 0.0,
        "40-step PCMI preserves the PETSc workspace, fixed physical "
        "residual scales, global balances, and nonnegative dissipation");
}

bool test_steady_load_cutback(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::SteadyProblem problem(input.spatial, mesh);
    const fuelsim::ProblemStateSnapshot snapshot = problem.capture_internal_state();
    bool empty_snapshot_rejected = false;
    try {
        problem.restore_internal_state(fuelsim::ProblemStateSnapshot{}, problem.initial_state());
    } catch (const std::invalid_argument&) {
        empty_snapshot_rejected = true;
    }
    fuelsim::SteadyProblem other_problem(input.spatial, mesh);
    bool foreign_snapshot_rejected = false;
    try {
        other_problem.restore_internal_state(snapshot, other_problem.initial_state());
    } catch (const std::invalid_argument&) {
        foreign_snapshot_rejected = true;
    }
    fuelsim::SolverOptions solver_options{input.solver.absolute_tolerance,
        input.solver.relative_tolerance,
        input.solver.step_tolerance,
        30};
    solver_options.line_search = fuelsim::SolverOptions::LineSearch::backtracking;
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 12, 1.0e-4}, solver_options);
    std::cout << "steady_cutback_rejected_steps=" << result.rejected_steps.size() << '\n';
    std::cout << "steady_cutback_total_cutbacks=" << result.total_cutbacks << '\n';
    if (!result.rejected_steps.empty()) {
        std::cout << "steady_cutback_first_attempted_load=" << result.rejected_steps.front().attempted_load_factor
                  << '\n';
        std::cout << "steady_cutback_last_attempted_load=" << result.rejected_steps.back().attempted_load_factor
                  << '\n';
        std::cout << "steady_cutback_last_failure="
                  << fuelsim::solve_failure_category_name(result.rejected_steps.back().failure_category) << '\n';
        std::cout << "steady_cutback_last_message=" << result.rejected_steps.back().failure_message << '\n';
    }
    bool passed = check(!snapshot.empty() && empty_snapshot_rejected && foreign_snapshot_rejected,
                      "steady snapshots reject empty and foreign problem state")
                  && check(result.completed && result.solve.converged && result.total_cutbacks > 0
                               && !result.rejected_steps.empty(),
                      "steady loading bisects a failed nominal increment "
                      "and continues from the accepted state");
    fuelsim::SteadyProblem predicted_problem(input.spatial, mesh);
    fuelsim::SteadyLoadOptions predicted_loading{1, 0.5, 18, 1.0e-6};
    predicted_loading.use_linear_load_predictor = true;
    fuelsim::SolverOptions predicted_solver = solver_options;
    predicted_solver.maximum_iterations = 1;
    const auto predicted = fuelsim::solve_steady(predicted_problem, predicted_loading, predicted_solver);
    std::cout << "steady_load_predictor_attempts=" << predicted.load_predictor_attempts << '\n';
    std::cout << "steady_load_predictor_fallbacks=" << predicted.load_predictor_fallbacks << '\n';
    bool predicted_fields_equal = true;
    for (const auto& field : problem.field_layout()) {
        double difference = 0.0, scale = 0.0;
        for (std::size_t i = field.begin; i < field.end; ++i) {
            difference = std::hypot(difference, predicted.solve.state[i] - result.solve.state[i]);
            scale = std::hypot(scale, result.solve.state[i]);
        }
        const double absolute = field.category == fuelsim::FieldCategory::thermal ? 1.0e-8 : 1.0e-12;
        predicted_fields_equal = predicted_fields_equal && difference < absolute + 1.0e-7 * scale;
    }
    passed = check(predicted.completed && predicted.load_predictor_attempts > 0
                       && predicted.load_predictor_fallbacks > 0 && predicted_fields_equal,
                 "load prediction across unequal accepted increments preserves the converged equilibrium")
             && passed;
    fuelsim::SteadyProblem minimum_problem(input.spatial, mesh);
    fuelsim::SolverOptions minimum_solver = solver_options;
    minimum_solver.maximum_iterations = 1;
    minimum_solver.line_search = fuelsim::SolverOptions::LineSearch::basic;
    minimum_solver.backtracking_fallback = false;
    const fuelsim::SteadyResult minimum = fuelsim::solve_steady(minimum_problem, {1, 0.5, 3, 0.75}, minimum_solver);
    passed = check(!minimum.completed && minimum.rejected_steps.size() == 2 && minimum.total_cutbacks == 1
                       && minimum.rejected_steps.back().load_increment == 0.75,
                 "steady loading attempts the minimum load increment once "
                 "before terminating")
             && passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m3_load_boundary_tests <internal.fsi>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M3.1 time loads and boundary test\n");
        if (!test_time_event_alignment(argv[1]) || !test_opaque_state_snapshot(argv[1])
            || !test_time_error_control(argv[1]) || !test_history_time_error_control(argv[1])
            || !test_long_transient_diagnostics(argv[1]) || !test_failure_diagnostics(argv[1])
            || !test_steady_load_cutback(argv[1]))
            return 1;
        std::cout << "[PASS] fuelsim M3.1 time loads and boundary test\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.1 test raised: " << error.what() << '\n';
        return 1;
    }
}
