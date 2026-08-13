#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
fuelsim::FuelSimCaseDefinition augmented_pcmi_input(
    const std::string& input_path, double penetration_tolerance, std::size_t maximum_augmented_iterations) {
    fuelsim::FuelSimCaseDefinition input = fuelsim::CaseInputReader::read(input_path);
    fuelsim::ContactDefinition& contact = input.contacts.at(0);
    contact.mechanical_formulation = fuelsim::MechanicalContactFormulation::augmented_lagrangian;
    contact.automatic_penalty = false;
    contact.penalty = 1.0e14;
    contact.penetration_tolerance = penetration_tolerance;
    contact.maximum_augmented_iterations = maximum_augmented_iterations;
    return input;
}
fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& input) {
    const fuelsim::TransientExecutionInput& execution = input.transient_execution;
    return {execution.end_time, execution.initial_time_step, execution.minimum_time_step, execution.maximum_time_step,
        execution.growth_factor, execution.cutback_factor, execution.maximum_cutbacks, execution.load_ramp_time,
        execution.target_nonlinear_iterations, execution.iteration_window, execution.time_error_relative_tolerance,
        execution.temperature_time_absolute_tolerance, execution.displacement_time_absolute_tolerance,
        execution.time_error_safety_factor, execution.strain_history_time_absolute_tolerance,
        execution.stress_history_time_absolute_tolerance};
}
fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& input) {
    fuelsim::SolverOptions options{input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    options.backtracking_fallback = input.solver.backtracking_fallback;
    options.field_residual_scaling = input.solver.field_residual_scaling;
    options.residual_reduction_tolerance = input.solver.residual_reduction_tolerance;
    options.temperature_residual_absolute_tolerance = input.solver.temperature_residual_absolute_tolerance;
    options.mechanical_residual_absolute_tolerance = input.solver.mechanical_residual_absolute_tolerance;
    options.temperature_residual_scale = input.solver.temperature_residual_scale;
    options.mechanical_residual_scale = input.solver.mechanical_residual_scale;
    return options;
}
double maximum_committed_multiplier(const fuelsim::TransientProblem& problem) {
    double maximum = 0.0;
    for (const std::vector<fuelsim::ContactPointHistory>& histories :
        fuelsim::rz::ProblemAccess::committed_state(problem).contact_histories)
        for (const fuelsim::ContactPointHistory& history : histories)
            maximum = std::max(maximum, history.normal_multiplier);
    return maximum;
}
class AugmentedContactObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep&) override {
        const fuelsim::InterfaceSummary interface =
            fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, problem.committed_solution());
        _maximum_penetration = std::max(_maximum_penetration, std::max(-interface.minimum_contact_gap, 0.0));
        _maximum_multiplier = std::max(_maximum_multiplier, maximum_committed_multiplier(problem));
        _final_active_contact_nodes = interface.active_contact_nodes;
        _final_total_contact_force = interface.total_contact_force;
    }
    double maximum_penetration() const noexcept { return _maximum_penetration; }
    double maximum_multiplier() const noexcept { return _maximum_multiplier; }
    std::size_t final_active_contact_nodes() const noexcept { return _final_active_contact_nodes; }
    double final_total_contact_force() const noexcept { return _final_total_contact_force; }

  private:
    double _maximum_penetration = 0.0;
    double _maximum_multiplier = 0.0;
    std::size_t _final_active_contact_nodes = 0;
    double _final_total_contact_force = 0.0;
};
bool test_transient_augmented_contact(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = augmented_pcmi_input(input_path, 1.0e-9, 50);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    AugmentedContactObserver observer;
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options(input), solver_options(input), &observer);
    std::cout << "transient_augmented_completed=" << result.completed << '\n'
              << "transient_augmented_accepted_steps=" << result.accepted_steps.size() << '\n'
              << "transient_augmented_rejected_steps=" << result.rejected_steps.size() << '\n'
              << "transient_augmented_solve_calls=" << result.aggregate_timing.solve_calls << '\n'
              << "transient_augmented_last_step_multiplier_updates="
              << result.last_attempt.augmented_lagrangian_iterations << '\n'
              << "transient_augmented_maximum_committed_multiplier=" << observer.maximum_multiplier() << '\n'
              << "transient_augmented_maximum_penetration=" << observer.maximum_penetration() << '\n'
              << "transient_augmented_penetration_tolerance=" << input.contacts.at(0).penetration_tolerance << '\n'
              << "transient_augmented_final_active_contact_nodes=" << observer.final_active_contact_nodes() << '\n'
              << "transient_augmented_final_total_contact_force=" << observer.final_total_contact_force() << '\n';
    return check(problem.uses_augmented_contact(), "transient fixture enables the augmented contact "
                                                   "formulation") &&
           check(result.completed && result.termination_reason == fuelsim::TransientTerminationReason::completed &&
                     result.rejected_steps.empty(),
               "transient augmented-contact solve completes without a "
               "rejected step") &&
           check(observer.maximum_multiplier() > 0.0, "at least one outer multiplier update is committed during "
                                                      "the transient solve") &&
           check(observer.maximum_penetration() <= input.contacts.at(0).penetration_tolerance,
               "every accepted transient step satisfies the 1 nm "
               "penetration tolerance") &&
           check(observer.final_active_contact_nodes() > 0 && observer.final_total_contact_force() > 0.0,
               "the transient augmented solve activates contact nodes") &&
           check(result.aggregate_timing.workspace_setups == 1, "transient augmented outer iterations reuse one PETSc "
                                                                "workspace");
}
bool test_transient_augmented_failure_rollback(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = augmented_pcmi_input(input_path, 1.0e-20, 1);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options(input), solver_options(input));
    std::cout << "transient_augmented_failure_committed_time=" << result.committed_time << '\n'
              << "transient_augmented_failure_termination="
              << fuelsim::transient_termination_reason_name(result.termination_reason) << '\n'
              << "transient_augmented_failure_rejected_steps=" << result.rejected_steps.size() << '\n'
              << "transient_augmented_failure_committed_multiplier=" << maximum_committed_multiplier(problem) << '\n';
    const bool cutback_terminated =
        result.termination_reason == fuelsim::TransientTerminationReason::maximum_cutbacks ||
        result.termination_reason == fuelsim::TransientTerminationReason::minimum_time_step;
    return check(!result.completed && cutback_terminated && !result.accepted_steps.empty() &&
                     result.committed_time > 0.0 && result.committed_time < input.transient_execution.end_time,
               "unreachable penetration tolerance stops the transient "
               "solve at the cutback or minimum-step limit") &&
           check(result.last_attempt.failure_category == fuelsim::SolveFailureCategory::contact_constraint,
               "the failed transient step is categorized as an augmented "
               "contact constraint failure") &&
           check(maximum_committed_multiplier(problem) == 0.0, "failed transient attempts roll back every committed "
                                                               "normal multiplier") &&
           check(result.aggregate_timing.workspace_setups == 1, "failed transient attempts reuse one PETSc workspace");
}
bool test_transient_augmented_time_error(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition input = augmented_pcmi_input(input_path, 1.0e-9, 50);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    fuelsim::TransientProblem problem(input.transient_definition(), mesh);
    fuelsim::TransientTimeOptions options = {
        input.transient_execution.end_time, 4.0, 0.125, 4.0, 1.0, 0.5, 12, input.transient_execution.load_ramp_time};
    options.time_error_relative_tolerance = 5.0e-1;
    options.temperature_time_absolute_tolerance = 1.0e-3;
    options.displacement_time_absolute_tolerance = 1.0e-9;
    options.strain_history_time_absolute_tolerance = 1.0e-10;
    options.stress_history_time_absolute_tolerance = 10.0;
    AugmentedContactObserver observer;
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, options, solver_options(input), &observer);
    double maximum_accepted_multiplier_error = 0.0;
    for (const fuelsim::TransientAcceptedStep& step : result.accepted_steps)
        maximum_accepted_multiplier_error =
            std::max(maximum_accepted_multiplier_error, step.time_error_components.contact_normal_multiplier);
    double maximum_rejected_multiplier_error = 0.0;
    bool multiplier_dominated_rejection = false;
    std::size_t time_discretization_rejections = 0;
    for (const fuelsim::TransientRejectedStep& step : result.rejected_steps) {
        if (step.failure_category != fuelsim::SolveFailureCategory::time_discretization) continue;
        ++time_discretization_rejections;
        const fuelsim::TransientTimeErrorEstimate& components = step.time_error_components;
        maximum_rejected_multiplier_error =
            std::max(maximum_rejected_multiplier_error, components.contact_normal_multiplier);
        multiplier_dominated_rejection =
            multiplier_dominated_rejection ||
            (components.contact_normal_multiplier > 1.0 && components.contact_normal_multiplier >= components.maximum);
    }
    std::cout << "transient_augmented_time_error_completed=" << result.completed << '\n'
              << "transient_augmented_time_error_accepted_steps=" << result.accepted_steps.size() << '\n'
              << "transient_augmented_time_error_rejections=" << result.time_error_rejections << '\n'
              << "transient_augmented_time_error_maximum_accepted_multiplier=" << maximum_accepted_multiplier_error
              << '\n'
              << "transient_augmented_time_error_maximum_rejected_multiplier=" << maximum_rejected_multiplier_error
              << '\n'
              << "transient_augmented_time_error_multiplier_dominated=" << multiplier_dominated_rejection << '\n';
    return check(result.completed && result.termination_reason == fuelsim::TransientTerminationReason::completed,
               "error-controlled transient augmented solve completes") &&
           check(time_discretization_rejections > 0, "step-doubling rejects at least one augmented-contact "
                                                     "step") &&
           check(maximum_rejected_multiplier_error > 0.0, "the normal-multiplier error component contributes to a "
                                                          "rejected step estimate") &&
           check(maximum_accepted_multiplier_error <= 1.0, "accepted steps keep the multiplier error component inside "
                                                           "the normalized bound") &&
           check(observer.maximum_penetration() <= input.contacts.at(0).penetration_tolerance,
               "error-controlled accepted steps satisfy the penetration "
               "tolerance") &&
           check(observer.maximum_multiplier() > 0.0 && observer.final_active_contact_nodes() > 0,
               "error-controlled solve commits active augmented "
               "multipliers") &&
           check(result.aggregate_timing.workspace_setups == 1, "step-doubling retries reuse one PETSc workspace");
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m54_transient_augmented_contact_tests "
                     "<transient.fsi>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M5.4 transient augmented-contact tests\n");
        bool passed = test_transient_augmented_contact(argv[1]);
        passed = test_transient_augmented_failure_rollback(argv[1]) && passed;
        passed = test_transient_augmented_time_error(argv[1]) && passed;
        if (!passed) return 1;
        std::cout << "[PASS] M5.4 transient augmented-contact end-to-end\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.4 transient augmented tests raised: " << error.what() << '\n';
        return 1;
    }
}
