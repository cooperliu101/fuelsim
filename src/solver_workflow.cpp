#include "solver_workflow.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace fuelsim::solver_workflow {
namespace {

void combine_attempt(SolveResult& aggregate, const SolveResult& addition) {
    const int nonlinear_iterations =
        aggregate.nonlinear_iterations + addition.nonlinear_iterations;
    const int linear_iterations =
        aggregate.linear_iterations + addition.linear_iterations;
    const std::size_t nonlinear_attempts =
        aggregate.nonlinear_attempts + addition.nonlinear_attempts;
    const std::size_t augmented_iterations =
        aggregate.augmented_lagrangian_iterations +
        addition.augmented_lagrangian_iterations;
    const double maximum_penetration =
        std::max(aggregate.maximum_contact_penetration,
                 addition.maximum_contact_penetration);
    SolveTiming timing = aggregate.timing;
    accumulate_timing(timing, addition.timing);
    const bool used_backtracking = aggregate.used_backtracking_fallback ||
                                   addition.used_backtracking_fallback;
    const SolveFailureCategory basic_failure =
        aggregate.basic_failure_category != SolveFailureCategory::none
            ? aggregate.basic_failure_category
            : addition.basic_failure_category;
    const std::string basic_message =
        !aggregate.basic_failure_message.empty()
            ? aggregate.basic_failure_message
            : addition.basic_failure_message;
    aggregate = addition;
    aggregate.nonlinear_iterations = nonlinear_iterations;
    aggregate.linear_iterations = linear_iterations;
    aggregate.nonlinear_attempts = nonlinear_attempts;
    aggregate.timing = timing;
    aggregate.used_backtracking_fallback = used_backtracking;
    aggregate.basic_failure_category = basic_failure;
    aggregate.basic_failure_message = basic_message;
    aggregate.augmented_lagrangian_iterations = augmented_iterations;
    aggregate.maximum_contact_penetration = maximum_penetration;
}

void mark_augmented_failure(SolveResult& result,
                            const AugmentedContactUpdate& status,
                            std::size_t completed_updates) {
    result.converged = false;
    result.failure_category = SolveFailureCategory::contact_constraint;
    result.failure_message =
        "Augmented contact did not reach penetration tolerance after " +
        std::to_string(completed_updates) +
        " multiplier updates; maximum constraint violation=" +
        std::to_string(status.maximum_constraint_violation) +
        ", maximum penetration=" +
        std::to_string(status.maximum_penetration) +
        ", tolerance=" + std::to_string(status.penetration_tolerance);
}

} // namespace

double elapsed_seconds(const SteadyClock::time_point& start) {
    return std::chrono::duration<double>(SteadyClock::now() - start).count();
}

void accumulate_timing(SolveTiming& total, const SolveTiming& step) {
    total.setup_seconds += step.setup_seconds;
    total.nonlinear_solve_seconds += step.nonlinear_solve_seconds;
    total.residual_callback_seconds += step.residual_callback_seconds;
    total.jacobian_callback_seconds += step.jacobian_callback_seconds;
    total.total_seconds += step.total_seconds;
    total.residual_evaluations += step.residual_evaluations;
    total.jacobian_evaluations += step.jacobian_evaluations;
    total.workspace_setups += step.workspace_setups;
    total.solve_calls += step.solve_calls;
}

void merge_attempt(SolveResult& aggregate, const SolveResult& addition) {
    combine_attempt(aggregate, addition);
}

std::vector<double>
initial_guess_with_dirichlet_values(const NonlinearProblem& problem,
                                    const std::vector<double>& state) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument(
            "Dirichlet initial-guess state size does not match problem");
    std::vector<double> result = state;
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = condition.value;
    return result;
}

SolveResult solve_contact_equilibrium(
    PetscSolver& solver, SteadyProblem& problem,
    const std::vector<double>& initial_guess, const SolverOptions& options) {
    SolveResult result = solver.solve(problem, initial_guess, options);
    if (!problem.uses_augmented_contact())
        return result;
    std::size_t updates = 0;
    while (result.converged) {
        const AugmentedContactUpdate status =
            problem.update_augmented_contact_multipliers(result.state, updates);
        result.maximum_contact_penetration = status.maximum_penetration;
        result.augmented_lagrangian_iterations = updates;
        if (status.converged)
            return result;
        if (!status.update_allowed) {
            mark_augmented_failure(result, status, updates);
            return result;
        }
        ++updates;
        SolveResult next = solver.solve(
            problem,
            initial_guess_with_dirichlet_values(problem, result.state),
            options);
        combine_attempt(result, next);
    }
    result.augmented_lagrangian_iterations = updates;
    return result;
}

SolveResult solve_contact_equilibrium(
    PetscSolver& solver, TransientProblem& problem,
    const std::vector<double>& initial_guess, const SolverOptions& options) {
    SolveResult result = solver.solve(problem, initial_guess, options);
    if (!problem.uses_augmented_contact())
        return result;
    std::size_t updates = 0;
    while (result.converged) {
        const AugmentedContactUpdate status =
            problem.update_augmented_contact_multipliers(result.state, updates);
        result.maximum_contact_penetration = status.maximum_penetration;
        result.augmented_lagrangian_iterations = updates;
        if (status.converged)
            return result;
        if (!status.update_allowed) {
            mark_augmented_failure(result, status, updates);
            return result;
        }
        ++updates;
        SolveResult next = solver.solve(
            problem,
            initial_guess_with_dirichlet_values(problem, result.state),
            options);
        combine_attempt(result, next);
    }
    result.augmented_lagrangian_iterations = updates;
    return result;
}

} // namespace fuelsim::solver_workflow
