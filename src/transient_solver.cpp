#include "fuelsim/problem_solver.hpp"

#include "solver_workflow.hpp"
#include "transient_time_control.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

using solver_workflow::SteadyClock;
using solver_workflow::accumulate_timing;
using solver_workflow::elapsed_seconds;
using solver_workflow::initial_guess_with_dirichlet_values;
using solver_workflow::merge_attempt;
using solver_workflow::solve_contact_equilibrium;

class TimeStepTransaction final {
  public:
    TimeStepTransaction(TransientProblem& problem,
                        const TransientStepInput& input)
        : _problem(problem), _committed(false) {
        _problem.begin_time_step(input);
    }

    ~TimeStepTransaction() {
        if (!_committed)
            _problem.rollback_time_step();
    }

    TimeStepTransaction(const TimeStepTransaction&) = delete;
    TimeStepTransaction& operator=(const TimeStepTransaction&) = delete;

    void commit(const std::vector<double>& solution) {
        _problem.commit_time_step(solution);
        _committed = true;
    }

  private:
    TransientProblem& _problem;
    bool _committed;
};

bool reaches_end(double time, double end_time) {
    const double scale = std::max({1.0, std::abs(time), std::abs(end_time)});
    return end_time - time <=
           16.0 * std::numeric_limits<double>::epsilon() * scale;
}

void validate_time_options(const TransientProblem& problem,
                           const TransientTimeOptions& options) {
    if (problem.time_step_active())
        throw std::logic_error(
            "solve_transient cannot start with an active time step");
    const double time_scale = std::max(
        {1.0, std::abs(problem.committed_time()), std::abs(options.end_time)});
    const double time_tolerance =
        16.0 * std::numeric_limits<double>::epsilon() * time_scale;
    if (!std::isfinite(options.end_time) ||
        options.end_time < problem.committed_time() - time_tolerance)
        throw std::invalid_argument(
            "solve_transient end time must not precede committed time");
    if (!std::isfinite(options.initial_time_step) ||
        !std::isfinite(options.minimum_time_step) ||
        !std::isfinite(options.maximum_time_step) ||
        !(options.minimum_time_step > 0.0) ||
        options.initial_time_step < options.minimum_time_step ||
        options.maximum_time_step < options.initial_time_step)
        throw std::invalid_argument("solve_transient requires 0 < minimum <= "
                                    "initial <= maximum time step");
    if (!std::isfinite(options.growth_factor) || options.growth_factor < 1.0)
        throw std::invalid_argument(
            "solve_transient growth factor must be at least one");
    if (!std::isfinite(options.cutback_factor) ||
        !(options.cutback_factor > 0.0 && options.cutback_factor < 1.0))
        throw std::invalid_argument(
            "solve_transient cutback factor must lie between zero and one");
    if (!std::isfinite(options.load_ramp_time) || options.load_ramp_time < 0.0)
        throw std::invalid_argument(
            "solve_transient ramp time must be finite and nonnegative");
    if (options.target_nonlinear_iterations == 0 &&
        options.iteration_window != 0)
        throw std::invalid_argument(
            "solve_transient iteration window requires a target");
    if (options.target_nonlinear_iterations > 0 &&
        options.iteration_window >= options.target_nonlinear_iterations)
        throw std::invalid_argument(
            "solve_transient iteration window must be smaller than target");
    if (!std::isfinite(options.time_error_relative_tolerance) ||
        options.time_error_relative_tolerance < 0.0 ||
        !std::isfinite(options.temperature_time_absolute_tolerance) ||
        !(options.temperature_time_absolute_tolerance > 0.0) ||
        !std::isfinite(options.displacement_time_absolute_tolerance) ||
        !(options.displacement_time_absolute_tolerance > 0.0) ||
        !std::isfinite(
            options.strain_history_time_absolute_tolerance) ||
        !(options.strain_history_time_absolute_tolerance > 0.0) ||
        !std::isfinite(
            options.stress_history_time_absolute_tolerance) ||
        !(options.stress_history_time_absolute_tolerance > 0.0) ||
        !std::isfinite(options.time_error_safety_factor) ||
        !(options.time_error_safety_factor > 0.0 &&
          options.time_error_safety_factor < 1.0))
        throw std::invalid_argument(
            "solve_transient time-error tolerances must be finite and "
            "nonnegative/positive, and safety factor must lie in (0, 1)");
}

double load_factor_at_time(const TransientTimeOptions& options, double time) {
    if (options.load_ramp_time == 0.0)
        return 1.0;
    return std::min(time / options.load_ramp_time, 1.0);
}

double accepted_next_time_step(const TransientTimeOptions& options,
                               double actual_time_step,
                               double controller_time_step,
                               bool event_truncated, std::size_t cutbacks,
                               int nonlinear_iterations) {
    const double base = event_truncated && cutbacks == 0 ? controller_time_step
                                                         : actual_time_step;
    if (options.target_nonlinear_iterations == 0) {
        if (cutbacks > 0)
            return std::clamp(base, options.minimum_time_step,
                              options.maximum_time_step);
        return std::min(options.maximum_time_step,
                        base * options.growth_factor);
    }
    const std::size_t iterations =
        nonlinear_iterations < 0
            ? 0
            : static_cast<std::size_t>(nonlinear_iterations);
    const std::size_t lower =
        options.target_nonlinear_iterations - options.iteration_window;
    const std::size_t upper =
        options.target_nonlinear_iterations >
                std::numeric_limits<std::size_t>::max() -
                    options.iteration_window
            ? std::numeric_limits<std::size_t>::max()
            : options.target_nonlinear_iterations + options.iteration_window;
    if (iterations < lower && cutbacks == 0)
        return std::min(options.maximum_time_step,
                        base * options.growth_factor);
    if (iterations > upper)
        return std::max(options.minimum_time_step,
                        base * options.cutback_factor);
    return std::clamp(base, options.minimum_time_step,
                      options.maximum_time_step);
}

} // namespace

TransientResult solve_transient(TransientProblem& problem,
                                const TransientTimeOptions& options,
                                const SolverOptions& solver_options,
                                TransientStepObserver* observer) {
    validate_time_options(problem, options);
    const SteadyClock::time_point start = SteadyClock::now();
    TransientResult result;
    PetscSolver solver;
    const std::vector<double> events = problem.time_events();
    double next_time_step = options.initial_time_step;
    while (!reaches_end(problem.committed_time(), options.end_time)) {
        const double controller_time_step = std::min(
            next_time_step, options.end_time - problem.committed_time());
        double time_step = controller_time_step;
        bool event_truncated = false;
        for (const double event : events) {
            if (reaches_end(problem.committed_time(), event))
                continue;
            if (event >= options.end_time ||
                reaches_end(event, options.end_time))
                break;
            const double event_step = event - problem.committed_time();
            if (event_step < time_step &&
                !reaches_end(event, problem.committed_time() + time_step)) {
                time_step = event_step;
                event_truncated = true;
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
            const bool error_control =
                options.time_error_relative_tolerance > 0.0;
            TransientCommittedState base_state;
            if (error_control)
                base_state = problem.committed_state();
            try {
                if (!error_control) {
                    TimeStepTransaction transaction(
                        problem,
                        {end_time, load_factor_at_time(options, end_time)});
                    attempt = solve_contact_equilibrium(
                        solver,
                        problem,
                        initial_guess_with_dirichlet_values(
                            problem, problem.committed_solution()),
                        solver_options);
                    controller_nonlinear_iterations =
                        attempt.nonlinear_iterations;
                    if (attempt.converged)
                        transaction.commit(attempt.state);
                } else {
                    SolveResult full_step;
                    TransientCommittedState full_step_state;
                    {
                        TimeStepTransaction transaction(
                            problem,
                            {end_time,
                             load_factor_at_time(options, end_time)});
                        full_step = solve_contact_equilibrium(
                            solver,
                            problem,
                            initial_guess_with_dirichlet_values(
                                problem, problem.committed_solution()),
                            solver_options);
                        if (full_step.converged) {
                            transaction.commit(full_step.state);
                            full_step_state = problem.committed_state();
                        }
                    }
                    attempt = full_step;
                    controller_nonlinear_iterations =
                        full_step.nonlinear_iterations;
                    if (full_step.converged) {
                        problem.restore_committed_state(base_state);
                        const double half_time =
                            base_state.time + 0.5 * time_step;
                        SolveResult first_half;
                        {
                            TimeStepTransaction transaction(
                                problem,
                                {half_time,
                                 load_factor_at_time(options, half_time)});
                            first_half = solve_contact_equilibrium(
                                solver,
                                problem,
                                initial_guess_with_dirichlet_values(
                                    problem, problem.committed_solution()),
                                solver_options);
                            if (first_half.converged)
                                transaction.commit(first_half.state);
                        }
                        if (first_half.converged)
                            first_half_conservation =
                                problem.last_conservation_summary();
                        controller_nonlinear_iterations = std::max(
                            controller_nonlinear_iterations,
                            first_half.nonlinear_iterations);
                        merge_attempt(attempt, first_half);
                        if (!first_half.converged) {
                            problem.restore_committed_state(
                                std::move(base_state));
                        } else {
                            SolveResult second_half;
                            {
                                TimeStepTransaction transaction(
                                    problem,
                                    {end_time,
                                     load_factor_at_time(options, end_time)});
                                second_half = solve_contact_equilibrium(
                                    solver,
                                    problem,
                                    initial_guess_with_dirichlet_values(
                                        problem,
                                        problem.committed_solution()),
                                    solver_options);
                                if (second_half.converged)
                                    transaction.commit(second_half.state);
                            }
                            controller_nonlinear_iterations = std::max(
                                controller_nonlinear_iterations,
                                second_half.nonlinear_iterations);
                            merge_attempt(attempt, second_half);
                            if (!second_half.converged) {
                                problem.restore_committed_state(
                                    std::move(base_state));
                            } else {
                                time_error_components =
                                    time_control::step_doubling_error(
                                    full_step_state,
                                    problem.committed_state(), options);
                                time_error_estimate =
                                    time_error_components.maximum;
                                if (!(time_error_estimate <= 1.0)) {
                                    problem.restore_committed_state(
                                        std::move(base_state));
                                    attempt.converged = false;
                                    attempt.failure_category =
                                        SolveFailureCategory::
                                            time_discretization;
                                    attempt.failure_message =
                                        "Backward-Euler step-doubling error "
                                        "exceeded one";
                                    ++result.time_error_rejections;
                                } else {
                                    TransientCommittedState accepted_state =
                                        problem.committed_state();
                                    accepted_state.conservation =
                                        time_control::combine_half_step_conservation(
                                            first_half_conservation,
                                            accepted_state.conservation);
                                    problem.restore_committed_state(
                                        std::move(accepted_state));
                                }
                            }
                        }
                    }
                }
            } catch (const std::domain_error& error) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                attempt.converged = false;
                attempt.failure_category =
                    SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (const std::overflow_error& error) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                attempt.converged = false;
                attempt.failure_category =
                    SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (...) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                throw;
            }
            accumulate_timing(result.aggregate_timing, attempt.timing);
            result.total_nonlinear_iterations += attempt.nonlinear_iterations;
            result.total_linear_iterations += attempt.linear_iterations;
            result.last_attempt = std::move(attempt);
            if (result.last_attempt.converged) {
                std::vector<RegionInelasticSummary> histories;
                histories.reserve(problem.region_count());
                for (std::size_t region = 0; region < problem.region_count();
                     ++region)
                    histories.push_back(
                        problem.summarize_region_history(region));
                next_time_step = accepted_next_time_step(
                    options, time_step, controller_time_step, event_truncated,
                    cutbacks, controller_nonlinear_iterations);
                if (error_control) {
                    const double error_limited_step =
                        time_step *
                        time_control::step_factor(options, time_error_estimate);
                    next_time_step =
                        std::clamp(std::min(next_time_step, error_limited_step),
                                   options.minimum_time_step,
                                   options.maximum_time_step);
                }
                result.accepted_steps.push_back(
                    {problem.committed_time(), time_step, next_time_step,
                     problem.committed_load_factor(), cutbacks,
                     result.last_attempt.nonlinear_iterations,
                     result.last_attempt.linear_iterations,
                     std::move(histories), time_error_estimate,
                     time_error_components,
                     problem.last_conservation_summary()});
                if (observer != nullptr)
                    observer->accepted_step(problem,
                                            result.accepted_steps.back());
                break;
            }
            result.rejected_steps.push_back(
                {problem.committed_time() + time_step, time_step, cutbacks,
                 result.last_attempt.nonlinear_iterations,
                 result.last_attempt.linear_iterations,
                 result.last_attempt.convergence_reason,
                 result.last_attempt.residual_norm,
                 result.last_attempt.failure_category,
                 result.last_attempt.failure_message, time_error_estimate,
                 time_error_components});
            if (cutbacks >= options.maximum_cutbacks_per_step) {
                result.termination_reason =
                    TransientTerminationReason::maximum_cutbacks;
                break;
            }
            double reduction_factor = options.cutback_factor;
            if (result.last_attempt.failure_category ==
                SolveFailureCategory::time_discretization)
                reduction_factor = std::min(
                    0.9,
                    std::max(options.cutback_factor,
                             time_control::step_factor(options,
                                                       time_error_estimate)));
            const double reduced = time_step * reduction_factor;
            const double minimum_tolerance =
                16.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, options.minimum_time_step);
            if (time_step <= options.minimum_time_step + minimum_tolerance) {
                result.termination_reason =
                    TransientTerminationReason::minimum_time_step;
                break;
            }
            time_step = std::max(reduced, options.minimum_time_step);
            ++cutbacks;
            ++result.total_cutbacks;
        }
        if (!result.last_attempt.converged)
            break;
    }
    result.completed = reaches_end(problem.committed_time(), options.end_time);
    if (result.completed)
        result.termination_reason = TransientTerminationReason::completed;
    result.committed_state = problem.committed_solution();
    result.committed_time = problem.committed_time();
    result.next_time_step = next_time_step;
    result.total_seconds = elapsed_seconds(start);
    return result;
}

const char*
transient_termination_reason_name(TransientTerminationReason reason) noexcept {
    switch (reason) {
    case TransientTerminationReason::not_started:
        return "not_started";
    case TransientTerminationReason::completed:
        return "completed";
    case TransientTerminationReason::maximum_cutbacks:
        return "maximum_cutbacks";
    case TransientTerminationReason::minimum_time_step:
        return "minimum_time_step";
    }
    return "unknown";
}

} // namespace fuelsim
