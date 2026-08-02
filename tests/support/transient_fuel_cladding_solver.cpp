#include "support/transient_fuel_cladding_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

using SteadyClock = std::chrono::steady_clock;

class TimeStepTransaction final {
  public:
    TimeStepTransaction(TransientFuelCladdingProblem& problem,
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

    void commit(const std::vector<double>& converged_solution) {
        _problem.commit_time_step(converged_solution);
        _committed = true;
    }

  private:
    TransientFuelCladdingProblem& _problem;
    bool _committed;
};

double seconds_since(const SteadyClock::time_point& start) {
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

void validate_time_options(const TransientFuelCladdingProblem& problem,
                           const TransientTimeOptions& options) {
    if (problem.time_step_active())
        throw std::logic_error("TransientFuelCladdingTimeStepper cannot start "
                               "with an active problem time step");
    if (!std::isfinite(options.end_time) ||
        !(options.end_time > problem.committed_time()))
        throw std::invalid_argument("TransientFuelCladdingTimeStepper end_time "
                                    "must be finite and greater than the "
                                    "committed time");
    if (!std::isfinite(options.initial_time_step) ||
        !std::isfinite(options.minimum_time_step) ||
        !std::isfinite(options.maximum_time_step) ||
        !(options.minimum_time_step > 0.0) ||
        !(options.initial_time_step >= options.minimum_time_step) ||
        !(options.maximum_time_step >= options.initial_time_step))
        throw std::invalid_argument("TransientFuelCladdingTimeStepper requires "
                                    "0 < minimum_time_step <= "
                                    "initial_time_step <= maximum_time_step");
    if (!std::isfinite(options.growth_factor) ||
        !(options.growth_factor >= 1.0))
        throw std::invalid_argument(
            "TransientFuelCladdingTimeStepper growth_factor must be finite and "
            "at least one");
    if (!std::isfinite(options.cutback_factor) ||
        !(options.cutback_factor > 0.0 && options.cutback_factor < 1.0))
        throw std::invalid_argument(
            "TransientFuelCladdingTimeStepper cutback_factor must lie strictly "
            "between zero and "
            "one");
    if (!std::isfinite(options.heat_source_ramp_time) ||
        !(options.heat_source_ramp_time >= 0.0))
        throw std::invalid_argument("TransientFuelCladdingTimeStepper "
                                    "heat_source_ramp_time must be finite and "
                                    "nonnegative");
}

double heat_source_at_time(const TransientFuelCladdingProblem& problem,
                           const TransientTimeOptions& options, double time) {
    const double target = problem.parameters().steady.volumetric_heat_source;
    if (options.heat_source_ramp_time == 0.0)
        return target;
    return target * std::min(time / options.heat_source_ramp_time, 1.0);
}

bool reaches_end(double time, double end_time) {
    const double scale = std::max({1.0, std::abs(time), std::abs(end_time)});
    return end_time - time <=
           16.0 * std::numeric_limits<double>::epsilon() * scale;
}

} // namespace

TransientFuelCladdingResult TransientFuelCladdingTimeStepper::solve(
    TransientFuelCladdingProblem& problem,
    const TransientTimeOptions& time_options,
    const SolverOptions& solver_options) const {
    validate_time_options(problem, time_options);
    const SteadyClock::time_point total_start = SteadyClock::now();

    TransientFuelCladdingResult result;
    PetscSolver nonlinear_solver;
    double next_time_step = time_options.initial_time_step;

    while (!reaches_end(problem.committed_time(), time_options.end_time)) {
        double time_step = std::min(
            next_time_step, time_options.end_time - problem.committed_time());
        std::size_t cutbacks = 0;

        for (;;) {
            const double end_time = problem.committed_time() + time_step;
            const double heat_source =
                heat_source_at_time(problem, time_options, end_time);

            SolveResult attempt;
            {
                TimeStepTransaction transaction(problem,
                                                {end_time, heat_source});
                attempt = nonlinear_solver.solve(
                    problem, problem.committed_solution(), solver_options);
                accumulate_timing(result.aggregate_timing, attempt.timing);
                result.total_nonlinear_iterations +=
                    attempt.nonlinear_iterations;
                if (attempt.converged)
                    transaction.commit(attempt.state);
            }

            result.last_attempt = std::move(attempt);
            if (result.last_attempt.converged) {
                result.accepted_steps.push_back(
                    {problem.committed_time(), time_step,
                     problem.committed_heat_source(), cutbacks,
                     result.last_attempt.nonlinear_iterations,
                     problem.summarize_fuel_history(),
                     problem.summarize_cladding_history()});
                next_time_step =
                    std::min(time_options.maximum_time_step,
                             time_step * time_options.growth_factor);
                break;
            }

            if (cutbacks >= time_options.maximum_cutbacks_per_step)
                break;
            const double reduced_time_step =
                time_step * time_options.cutback_factor;
            if (reduced_time_step < time_options.minimum_time_step)
                break;
            time_step = reduced_time_step;
            ++cutbacks;
            ++result.total_cutbacks;
        }

        if (!result.last_attempt.converged)
            break;
    }

    result.completed =
        reaches_end(problem.committed_time(), time_options.end_time);
    result.committed_state = problem.committed_solution();
    result.committed_time = problem.committed_time();
    result.total_seconds = seconds_since(total_start);
    return result;
}

} // namespace fuelsim
