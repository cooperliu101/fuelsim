#ifndef FUELSIM_PROBLEM_SOLVER_HPP
#define FUELSIM_PROBLEM_SOLVER_HPP

#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

#include <cstddef>
#include <vector>

namespace fuelsim {

struct SteadyResult final {
    SolveResult solve;
    std::size_t completed_steps = 0;
    bool completed = false;
    int total_nonlinear_iterations = 0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
};

SteadyResult solve_steady(SteadyProblem& problem, std::size_t load_steps,
                          const SolverOptions& options = SolverOptions{});

struct TransientTimeOptions final {
    double end_time;
    double initial_time_step;
    double minimum_time_step;
    double maximum_time_step;
    double growth_factor;
    double cutback_factor;
    std::size_t maximum_cutbacks_per_step;
    double load_ramp_time;
};

struct TransientAcceptedStep final {
    double time;
    double time_step;
    double load_factor;
    std::size_t cutbacks;
    int nonlinear_iterations;
    std::vector<RegionInelasticSummary> region_histories;
};

struct TransientResult final {
    SolveResult last_attempt;
    std::vector<double> committed_state;
    std::vector<TransientAcceptedStep> accepted_steps;
    bool completed = false;
    std::size_t total_cutbacks = 0;
    int total_nonlinear_iterations = 0;
    double committed_time = 0.0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
};

class TransientStepObserver {
  public:
    virtual ~TransientStepObserver() = default;
    virtual void accepted_step(const TransientProblem& problem,
                               const TransientAcceptedStep& step) = 0;
};

TransientResult
solve_transient(TransientProblem& problem,
                const TransientTimeOptions& time_options,
                const SolverOptions& solver_options = SolverOptions{},
                TransientStepObserver* observer = nullptr);

} // namespace fuelsim

#endif
