#include "fuelsim/solver/solve_workflows.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include "solver_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
using SteadyClock = solver_detail::Clock;
using solver_detail::accumulate_timing;
using solver_detail::seconds_since;

namespace solver_workflow {
namespace {
void merge_attempt(SolveResult& aggregate, const SolveResult& addition) {
    const int nonlinear_iterations = aggregate.nonlinear_iterations + addition.nonlinear_iterations,
              linear_iterations = aggregate.linear_iterations + addition.linear_iterations;
    const std::size_t nonlinear_attempts = aggregate.nonlinear_attempts + addition.nonlinear_attempts;
    const std::size_t augmented_iterations =
        aggregate.augmented_lagrangian_iterations + addition.augmented_lagrangian_iterations;
    const double maximum_penetration =
        std::max(aggregate.maximum_contact_penetration, addition.maximum_contact_penetration);
    SolveTiming timing = aggregate.timing;
    accumulate_timing(timing, addition.timing);
    const bool used_backtracking = aggregate.used_backtracking_fallback || addition.used_backtracking_fallback;
    const SolveFailureCategory basic_failure = aggregate.basic_failure_category != SolveFailureCategory::none
                                                   ? aggregate.basic_failure_category
                                                   : addition.basic_failure_category;
    const std::string basic_message =
        !aggregate.basic_failure_message.empty() ? aggregate.basic_failure_message : addition.basic_failure_message;
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

void mark_augmented_failure(SolveResult& result, const AugmentedContactUpdate& status, std::size_t completed_updates) {
    result.converged = false;
    result.failure_category = SolveFailureCategory::contact_constraint;
    result.failure_message =
        "Augmented contact did not reach penetration tolerance after " + std::to_string(completed_updates) +
        " multiplier updates; maximum constraint violation=" + std::to_string(status.maximum_constraint_violation) +
        ", maximum penetration=" + std::to_string(status.maximum_penetration) +
        ", tolerance=" + std::to_string(status.penetration_tolerance);
}
} // namespace

std::vector<double> initial_guess_with_dirichlet_values(
    const NonlinearProblem& problem, const std::vector<double>& state) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument("Dirichlet initial-guess state size does not match problem");
    std::vector<double> result = state;
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = condition.value;
    return result;
}

SolveResult solve_contact_equilibrium(PetscSolver& solver, NonlinearProblem& problem,
    const std::vector<double>& initial_guess, const SolverOptions& options) {
    SolveResult result = solver.solve(problem, initial_guess, options);
    if (!problem.uses_augmented_contact()) return result;
    std::size_t updates = 0;
    while (result.converged) {
        const AugmentedContactUpdate status = problem.update_augmented_contact_multipliers(result.state, updates);
        result.maximum_contact_penetration = status.maximum_penetration;
        result.augmented_lagrangian_iterations = updates;
        if (status.converged) return result;
        if (!status.update_allowed) {
            mark_augmented_failure(result, status, updates);
            return result;
        }
        ++updates;
        SolveResult next = solver.solve(problem, initial_guess_with_dirichlet_values(problem, result.state), options);
        merge_attempt(result, next);
    }
    result.augmented_lagrangian_iterations = updates;
    return result;
}
} // namespace solver_workflow

using solver_workflow::initial_guess_with_dirichlet_values;
using solver_workflow::solve_contact_equilibrium;

SteadyResult solve_steady(SteadyProblem& problem, const SteadyLoadOptions& load_options, const SolverOptions& options) {
    const SteadyClock::time_point start = SteadyClock::now();
    SteadyResult result;
    PetscSolver solver;
    std::vector<double> state = problem.initial_state();
    double accepted_load_factor = 0.0;
    for (std::size_t step = 1; step <= load_options.load_steps; ++step) {
        const double target_load_factor = static_cast<double>(step) / static_cast<double>(load_options.load_steps);
        while (accepted_load_factor < target_load_factor) {
            std::size_t cutbacks = 0;
            double attempted_load_factor = target_load_factor,
                   load_increment = attempted_load_factor - accepted_load_factor;
            SolveResult attempt;
            for (;;) {
                const ProblemStateSnapshot internal_state = problem.capture_internal_state();
                attempt = SolveResult{};
                try {
                    problem.set_load_factor(attempted_load_factor);
                    attempt = solve_contact_equilibrium(
                        solver, problem, initial_guess_with_dirichlet_values(problem, state), options);
                    if (attempt.converged) problem.commit_internal_state(attempt.state);
                } catch (const std::domain_error& error) {
                    attempt.converged = false;
                    attempt.failure_category = SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                } catch (const std::overflow_error& error) {
                    attempt.converged = false;
                    attempt.failure_category = SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                }
                if (!attempt.converged) {
                    problem.restore_internal_state(internal_state, state);
                    problem.set_load_factor(accepted_load_factor);
                }
                accumulate_timing(result.aggregate_timing, attempt.timing);
                result.total_nonlinear_iterations += attempt.nonlinear_iterations;
                result.total_linear_iterations += attempt.linear_iterations;
                if (attempt.converged) break;
                result.rejected_steps.push_back({attempted_load_factor, load_increment, cutbacks,
                    attempt.failure_category, attempt.failure_message});
                result.solve = attempt;
                if (cutbacks >= load_options.maximum_cutbacks_per_step) {
                    result.total_seconds = seconds_since(start);
                    return result;
                }
                const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
                if (load_increment <= load_options.minimum_load_increment + tolerance) {
                    result.total_seconds = seconds_since(start);
                    return result;
                }
                load_increment =
                    std::max(load_increment * load_options.cutback_factor, load_options.minimum_load_increment);
                attempted_load_factor = accepted_load_factor + load_increment;
                ++cutbacks;
                ++result.total_cutbacks;
            }
            accepted_load_factor = attempted_load_factor;
            state = attempt.state;
            result.solve = std::move(attempt);
        }
        result.completed_steps = step;
    }
    result.completed = true;
    result.total_seconds = seconds_since(start);
    return result;
}

double step_factor(const TransientTimeOptions& options, double error) {
    if (!(error > 0.0)) return options.growth_factor;
    return std::clamp(options.time_error_safety_factor / std::sqrt(error), 0.1, options.growth_factor);
}

namespace {
using solver_workflow::merge_attempt;

bool reaches_end(double time, double end_time) {
    const double scale = std::max({1.0, std::abs(time), std::abs(end_time)});
    return end_time - time <= 16.0 * std::numeric_limits<double>::epsilon() * scale;
}

bool times_equal(double first, double second) {
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <= 16.0 * std::numeric_limits<double>::epsilon() * scale;
}

std::vector<double> linear_transient_predictor(const std::vector<double>& committed,
    const std::vector<double>* previous, double committed_time, double previous_time, double target_time) {
    if (previous == nullptr || previous->size() != committed.size() || !(committed_time > previous_time))
        return committed;
    const double factor = (target_time - committed_time) / (committed_time - previous_time);
    if (!std::isfinite(factor)) return committed;
    std::vector<double> result(committed.size());
    for (std::size_t dof = 0; dof < committed.size(); ++dof)
        result[dof] = committed[dof] + factor * (committed[dof] - (*previous)[dof]);
    return result;
}

void validate_time_options(const TransientProblem& problem, const TransientTimeOptions& options) {
    if (problem.time_step_active()) throw std::logic_error("solve_transient cannot start with an active time step");
    const double time_scale = std::max({1.0, std::abs(problem.committed_time()), std::abs(options.end_time)}),
                 time_tolerance = 16.0 * std::numeric_limits<double>::epsilon() * time_scale;
    if (!std::isfinite(options.end_time) || options.end_time < problem.committed_time() - time_tolerance)
        throw std::invalid_argument("solve_transient end time must not precede committed time");
}

double load_factor_at_time(const TransientTimeOptions& options, double time) {
    if (options.load_ramp_time == 0.0) return 1.0;
    return std::min(time / options.load_ramp_time, 1.0);
}

double accepted_next_time_step(const TransientTimeOptions& options, double actual_time_step,
    double controller_time_step, bool event_truncated, std::size_t cutbacks, int nonlinear_iterations) {
    const double base = event_truncated && cutbacks == 0 ? controller_time_step : actual_time_step;
    if (options.target_nonlinear_iterations == 0) {
        if (cutbacks > 0) return std::clamp(base, options.minimum_time_step, options.maximum_time_step);
        return std::min(options.maximum_time_step, base * options.growth_factor);
    }
    const std::size_t iterations = nonlinear_iterations < 0 ? 0 : static_cast<std::size_t>(nonlinear_iterations),
                      lower = options.target_nonlinear_iterations - options.iteration_window;
    const std::size_t upper =
        options.target_nonlinear_iterations > std::numeric_limits<std::size_t>::max() - options.iteration_window
            ? std::numeric_limits<std::size_t>::max()
            : options.target_nonlinear_iterations + options.iteration_window;
    if (iterations < lower && cutbacks == 0) return std::min(options.maximum_time_step, base * options.growth_factor);
    if (iterations > upper) return std::max(options.minimum_time_step, base * options.cutback_factor);
    return std::clamp(base, options.minimum_time_step, options.maximum_time_step);
}
} // namespace

TransientResult solve_transient(TransientProblem& problem, const TransientTimeOptions& options,
    const SolverOptions& solver_options, TransientStepObserver* observer) {
    validate_time_options(problem, options);
    const SteadyClock::time_point start = SteadyClock::now();
    TransientResult result;
    PetscSolver solver;
    const std::vector<double> events = problem.time_events();
    double next_time_step = options.initial_time_step;
    const std::vector<double> predictor_reference_state = problem.initial_solution();
    const auto run_step = [&](double target_time, const std::vector<double>& initial_guess) {
        problem.begin_time_step(
            {target_time, load_factor_at_time(options, target_time), options.include_thermal_time_term});
        try {
            SolveResult step_result = solve_contact_equilibrium(
                solver, problem, initial_guess_with_dirichlet_values(problem, initial_guess), solver_options);
            if (step_result.converged)
                problem.commit_time_step(step_result.state);
            else
                problem.rollback_time_step();
            return step_result;
        } catch (...) {
            problem.rollback_time_step();
            throw;
        }
    };
    const auto run_predicted_step = [&](double target_time, const std::vector<double>& predicted,
                                        const std::vector<double>& committed, bool predictor_used) {
        SolveResult step_result;
        try {
            step_result = run_step(target_time, predicted);
        } catch (const std::domain_error& error) {
            if (!predictor_used) throw;
            step_result.converged = false;
            step_result.failure_category = SolveFailureCategory::physical_domain;
            step_result.failure_message = error.what();
        } catch (const std::overflow_error& error) {
            if (!predictor_used) throw;
            step_result.converged = false;
            step_result.failure_category = SolveFailureCategory::physical_domain;
            step_result.failure_message = error.what();
        }
        if (step_result.converged || !predictor_used) return step_result;
        SolveResult fallback = run_step(target_time, committed);
        merge_attempt(step_result, fallback);
        return step_result;
    };
    while (!reaches_end(problem.committed_time(), options.end_time)) {
        const double controller_time_step = std::min(next_time_step, options.end_time - problem.committed_time());
        double time_step = controller_time_step;
        bool event_truncated = false;
        for (const double event : events) {
            if (reaches_end(problem.committed_time(), event)) continue;
            if (event > options.end_time && !times_equal(event, options.end_time)) break;
            const double event_step = event - problem.committed_time();
            if (event_step < time_step || times_equal(event, problem.committed_time() + time_step)) {
                event_truncated = event_step < time_step && !times_equal(event, problem.committed_time() + time_step);
                time_step = event_step;
            }
            break;
        }
        std::size_t cutbacks = 0;
        for (;;) {
            const double end_time = problem.committed_time() + time_step;
            SolveResult attempt;
            double time_error_estimate = 0.0;
            TransientTimeErrorEstimate time_error_components;
            TransientConservationSummary first_half_conservation;
            int controller_nonlinear_iterations = 0;
            const bool error_control = options.time_error_relative_tolerance > 0.0;
            ProblemStateSnapshot base_state;
            bool base_state_available = false;
            const double base_time = problem.committed_time();
            const std::vector<double> base_solution = problem.committed_solution();
            const bool predictor_used = options.use_linear_time_predictor && base_time > 0.0;
            const std::vector<double> predicted_solution = linear_transient_predictor(
                base_solution, predictor_used ? &predictor_reference_state : nullptr, base_time, 0.0, end_time);
            if (error_control) {
                base_state = problem.capture_state();
                base_state_available = true;
            }
            try {
                if (!error_control) {
                    attempt = run_predicted_step(end_time, predicted_solution, base_solution, predictor_used);
                    controller_nonlinear_iterations = attempt.nonlinear_iterations;
                } else {
                    const SolveResult full_step =
                        run_predicted_step(end_time, predicted_solution, base_solution, predictor_used);
                    ProblemStateSnapshot full_step_state;
                    if (full_step.converged) full_step_state = problem.capture_state();
                    attempt = full_step;
                    controller_nonlinear_iterations = full_step.nonlinear_iterations;
                    if (full_step.converged) {
                        problem.restore_state(base_state);
                        const double half_time = base_time + 0.5 * time_step;
                        const std::vector<double> first_half_prediction = linear_transient_predictor(base_solution,
                            predictor_used ? &predictor_reference_state : nullptr, base_time, 0.0, half_time);
                        const SolveResult first_half =
                            run_predicted_step(half_time, first_half_prediction, base_solution, predictor_used);
                        if (first_half.converged) first_half_conservation = problem.last_conservation_summary();
                        controller_nonlinear_iterations =
                            std::max(controller_nonlinear_iterations, first_half.nonlinear_iterations);
                        merge_attempt(attempt, first_half);
                        if (!first_half.converged) {
                            problem.restore_state(base_state);
                            base_state_available = false;
                        } else {
                            const std::vector<double> first_half_solution = problem.committed_solution();
                            const std::vector<double> second_half_prediction = linear_transient_predictor(
                                first_half_solution, &predictor_reference_state, half_time, 0.0, end_time);
                            const SolveResult second_half =
                                run_predicted_step(end_time, second_half_prediction, first_half_solution, true);
                            controller_nonlinear_iterations =
                                std::max(controller_nonlinear_iterations, second_half.nonlinear_iterations);
                            merge_attempt(attempt, second_half);
                            if (!second_half.converged) {
                                problem.restore_state(base_state);
                                base_state_available = false;
                            } else {
                                time_error_components =
                                    problem.step_doubling_error(full_step_state, problem.capture_state(), options);
                                time_error_estimate = time_error_components.maximum;
                                if (!(time_error_estimate <= 1.0)) {
                                    problem.restore_state(base_state);
                                    base_state_available = false;
                                    attempt.converged = false;
                                    attempt.failure_category = SolveFailureCategory::time_discretization;
                                    attempt.failure_message = "Backward-Euler step-doubling error exceeded one";
                                    ++result.time_error_rejections;
                                } else {
                                    problem.combine_last_half_step_conservation(first_half_conservation);
                                }
                            }
                        }
                    }
                }
            } catch (const std::domain_error& error) {
                if (base_state_available) problem.restore_state(base_state);
                attempt.converged = false;
                attempt.failure_category = SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (const std::overflow_error& error) {
                if (base_state_available) problem.restore_state(base_state);
                attempt.converged = false;
                attempt.failure_category = SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (...) {
                if (base_state_available) problem.restore_state(base_state);
                throw;
            }
            accumulate_timing(result.aggregate_timing, attempt.timing);
            result.total_nonlinear_iterations += attempt.nonlinear_iterations;
            result.total_linear_iterations += attempt.linear_iterations;
            result.last_attempt = std::move(attempt);
            if (result.last_attempt.converged) {
                next_time_step = accepted_next_time_step(options, time_step, controller_time_step, event_truncated,
                    cutbacks, controller_nonlinear_iterations);
                if (error_control) {
                    const double error_limited_step = time_step * step_factor(options, time_error_estimate);
                    next_time_step = std::clamp(std::min(next_time_step, error_limited_step), options.minimum_time_step,
                        options.maximum_time_step);
                }
                result.accepted_steps.push_back(
                    {problem.committed_time(), time_step, next_time_step, problem.committed_load_factor(), cutbacks,
                        result.last_attempt.nonlinear_iterations, result.last_attempt.linear_iterations,
                        time_error_estimate, time_error_components, problem.last_conservation_summary()});
                if (observer != nullptr) observer->accepted_step(problem, result.accepted_steps.back());
                break;
            }
            result.rejected_steps.push_back(
                {problem.committed_time() + time_step, time_step, cutbacks, result.last_attempt.nonlinear_iterations,
                    result.last_attempt.linear_iterations, result.last_attempt.convergence_reason,
                    result.last_attempt.residual_norm, result.last_attempt.failure_category,
                    result.last_attempt.failure_message, time_error_estimate, time_error_components});
            if (cutbacks >= options.maximum_cutbacks_per_step) {
                result.termination_reason = TransientTerminationReason::maximum_cutbacks;
                break;
            }
            double reduction_factor = options.cutback_factor;
            if (result.last_attempt.failure_category == SolveFailureCategory::time_discretization)
                reduction_factor =
                    std::min(0.9, std::max(options.cutback_factor, step_factor(options, time_error_estimate)));
            const double reduced = time_step * reduction_factor;
            const double minimum_tolerance =
                16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, options.minimum_time_step);
            if (time_step <= options.minimum_time_step + minimum_tolerance) {
                result.termination_reason = TransientTerminationReason::minimum_time_step;
                break;
            }
            time_step = std::max(reduced, options.minimum_time_step);
            ++cutbacks;
            ++result.total_cutbacks;
        }
        if (!result.last_attempt.converged) break;
    }
    result.completed = reaches_end(problem.committed_time(), options.end_time);
    if (result.completed) result.termination_reason = TransientTerminationReason::completed;
    result.committed_state = problem.committed_solution();
    result.committed_time = problem.committed_time();
    result.next_time_step = next_time_step;
    result.total_seconds = seconds_since(start);
    return result;
}

const char* transient_termination_reason_name(TransientTerminationReason reason) noexcept {
    switch (reason) {
    case TransientTerminationReason::not_started: return "not_started";
    case TransientTerminationReason::completed: return "completed";
    case TransientTerminationReason::maximum_cutbacks: return "maximum_cutbacks";
    case TransientTerminationReason::minimum_time_step: return "minimum_time_step";
    }
    return "unknown";
}
} // namespace fuelsim
