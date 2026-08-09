#include "fuelsim/problem_solver.hpp"

#include "solver_workflow.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {

using solver_workflow::SteadyClock;
using solver_workflow::accumulate_timing;
using solver_workflow::elapsed_seconds;
using solver_workflow::initial_guess_with_dirichlet_values;
using solver_workflow::solve_contact_equilibrium;

SteadyResult solve_steady(SteadyProblem& problem,
                          const SteadyLoadOptions& load_options,
                          const SolverOptions& options) {
    if (load_options.load_steps == 0)
        throw std::invalid_argument("solve_steady load_steps must be positive");
    if (!std::isfinite(load_options.cutback_factor) ||
        !(load_options.cutback_factor > 0.0 &&
          load_options.cutback_factor < 1.0))
        throw std::invalid_argument(
            "solve_steady cutback factor must lie between zero and one");
    if (!std::isfinite(load_options.minimum_load_increment) ||
        !(load_options.minimum_load_increment > 0.0) ||
        load_options.minimum_load_increment > 1.0)
        throw std::invalid_argument(
            "solve_steady minimum load increment must lie in (0, 1]");
    const SteadyClock::time_point start = SteadyClock::now();
    SteadyResult result;
    PetscSolver solver;
    std::vector<double> state = problem.initial_state();
    double accepted_load_factor = 0.0;
    for (std::size_t step = 1; step <= load_options.load_steps; ++step) {
        const double target_load_factor =
            static_cast<double>(step) /
            static_cast<double>(load_options.load_steps);
        while (accepted_load_factor < target_load_factor) {
            std::size_t cutbacks = 0;
            double attempted_load_factor = target_load_factor;
            double load_increment =
                attempted_load_factor - accepted_load_factor;
            SolveResult attempt;
            for (;;) {
                const std::vector<std::vector<ContactPointHistory>>
                    committed_contact_histories =
                        problem.committed_contact_histories();
                attempt = SolveResult{};
                try {
                    problem.set_load_factor(attempted_load_factor);
                    attempt = solve_contact_equilibrium(
                        solver,
                        problem,
                        initial_guess_with_dirichlet_values(problem, state),
                        options);
                    if (attempt.converged)
                        problem.commit_contact_state(attempt.state);
                } catch (const std::domain_error& error) {
                    attempt.converged = false;
                    attempt.failure_category =
                        SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                } catch (const std::overflow_error& error) {
                    attempt.converged = false;
                    attempt.failure_category =
                        SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                }
                if (!attempt.converged) {
                    problem.restore_contact_state(
                        state, committed_contact_histories);
                    problem.set_load_factor(accepted_load_factor);
                }
                accumulate_timing(result.aggregate_timing, attempt.timing);
                result.total_nonlinear_iterations +=
                    attempt.nonlinear_iterations;
                result.total_linear_iterations += attempt.linear_iterations;
                if (attempt.converged)
                    break;

                result.rejected_steps.push_back(
                    {attempted_load_factor, load_increment, cutbacks,
                     attempt.failure_category, attempt.failure_message});
                result.solve = attempt;
                if (cutbacks >=
                    load_options.maximum_cutbacks_per_step) {
                    result.total_seconds = elapsed_seconds(start);
                    return result;
                }
                const double tolerance =
                    16.0 * std::numeric_limits<double>::epsilon();
                if (load_increment <=
                    load_options.minimum_load_increment + tolerance) {
                    result.total_seconds = elapsed_seconds(start);
                    return result;
                }
                load_increment = std::max(
                    load_increment * load_options.cutback_factor,
                    load_options.minimum_load_increment);
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
    result.total_seconds = elapsed_seconds(start);
    return result;
}


} // namespace fuelsim
