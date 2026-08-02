#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

using SteadyClock = std::chrono::steady_clock;

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
}

double load_factor_at_time(const TransientTimeOptions& options, double time) {
    if (options.load_ramp_time == 0.0)
        return 1.0;
    return std::min(time / options.load_ramp_time, 1.0);
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

} // namespace

SteadyResult solve_steady(SteadyProblem& problem, std::size_t load_steps,
                          const SolverOptions& options) {
    if (load_steps == 0)
        throw std::invalid_argument("solve_steady load_steps must be positive");
    const SteadyClock::time_point start = SteadyClock::now();
    SteadyResult result;
    PetscSequentialSolver solver;
    std::vector<double> state = problem.initial_state();
    for (std::size_t step = 1; step <= load_steps; ++step) {
        problem.set_load_factor(static_cast<double>(step) /
                                static_cast<double>(load_steps));
        SolveResult attempt = solver.solve(
            problem, initial_guess_with_dirichlet_values(problem, state),
            options);
        accumulate_timing(result.aggregate_timing, attempt.timing);
        result.total_nonlinear_iterations += attempt.nonlinear_iterations;
        if (!attempt.converged) {
            result.solve = std::move(attempt);
            result.completed_steps = step - 1;
            result.total_seconds = seconds_since(start);
            return result;
        }
        result.completed_steps = step;
        if (step == load_steps)
            result.solve = std::move(attempt);
        else
            state = std::move(attempt.state);
    }
    result.completed = true;
    result.total_seconds = seconds_since(start);
    return result;
}

TransientResult solve_transient(TransientProblem& problem,
                                const TransientTimeOptions& options,
                                const SolverOptions& solver_options,
                                TransientStepObserver* observer) {
    validate_time_options(problem, options);
    const SteadyClock::time_point start = SteadyClock::now();
    TransientResult result;
    PetscSequentialSolver solver;
    double next_time_step = options.initial_time_step;
    while (!reaches_end(problem.committed_time(), options.end_time)) {
        double time_step = std::min(
            next_time_step, options.end_time - problem.committed_time());
        std::size_t cutbacks = 0;
        for (;;) {
            const double end_time = problem.committed_time() + time_step;
            SolveResult attempt;
            {
                TimeStepTransaction transaction(
                    problem,
                    {end_time, load_factor_at_time(options, end_time)});
                attempt =
                    solver.solve(problem,
                                 initial_guess_with_dirichlet_values(
                                     problem, problem.committed_solution()),
                                 solver_options);
                accumulate_timing(result.aggregate_timing, attempt.timing);
                result.total_nonlinear_iterations +=
                    attempt.nonlinear_iterations;
                if (attempt.converged)
                    transaction.commit(attempt.state);
            }
            result.last_attempt = std::move(attempt);
            if (result.last_attempt.converged) {
                std::vector<RegionInelasticSummary> histories;
                histories.reserve(problem.region_count());
                for (std::size_t region = 0; region < problem.region_count();
                     ++region)
                    histories.push_back(
                        problem.summarize_region_history(region));
                result.accepted_steps.push_back(
                    {problem.committed_time(), time_step,
                     problem.committed_load_factor(), cutbacks,
                     result.last_attempt.nonlinear_iterations,
                     std::move(histories)});
                if (observer != nullptr)
                    observer->accepted_step(problem,
                                            result.accepted_steps.back());
                next_time_step = std::min(options.maximum_time_step,
                                          time_step * options.growth_factor);
                break;
            }
            if (cutbacks >= options.maximum_cutbacks_per_step)
                break;
            const double reduced = time_step * options.cutback_factor;
            if (reduced < options.minimum_time_step)
                break;
            time_step = reduced;
            ++cutbacks;
            ++result.total_cutbacks;
        }
        if (!result.last_attempt.converged)
            break;
    }
    result.completed = reaches_end(problem.committed_time(), options.end_time);
    result.committed_state = problem.committed_solution();
    result.committed_time = problem.committed_time();
    result.total_seconds = seconds_since(start);
    return result;
}

} // namespace fuelsim
