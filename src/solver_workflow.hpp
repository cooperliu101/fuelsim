#ifndef FUELSIM_SOLVER_WORKFLOW_HPP
#define FUELSIM_SOLVER_WORKFLOW_HPP

#include "fuelsim/problem_solver.hpp"

#include <chrono>
#include <vector>

namespace fuelsim::solver_workflow {

using SteadyClock = std::chrono::steady_clock;

double elapsed_seconds(const SteadyClock::time_point& start);
void accumulate_timing(SolveTiming& total, const SolveTiming& step);
void merge_attempt(SolveResult& aggregate, const SolveResult& addition);

std::vector<double>
initial_guess_with_dirichlet_values(const NonlinearProblem& problem,
                                    const std::vector<double>& state);

SolveResult solve_contact_equilibrium(
    PetscSolver& solver, SteadyProblem& problem,
    const std::vector<double>& initial_guess, const SolverOptions& options);

SolveResult solve_contact_equilibrium(
    PetscSolver& solver, TransientProblem& problem,
    const std::vector<double>& initial_guess, const SolverOptions& options);

} // namespace fuelsim::solver_workflow

#endif
