#ifndef FUELSIM_TRANSIENT_FUEL_CLADDING_SOLVER_HPP
#define FUELSIM_TRANSIENT_FUEL_CLADDING_SOLVER_HPP

#include <cstddef>
#include <vector>

#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/transient_fuel_cladding_problem.hpp"

namespace fuelsim {

struct TransientTimeOptions final {
    double end_time;
    double initial_time_step;
    double minimum_time_step;
    double maximum_time_step;
    double growth_factor;
    double cutback_factor;
    std::size_t maximum_cutbacks_per_step;
    double heat_source_ramp_time;
};

struct TransientAcceptedStep final {
    double time;
    double time_step;
    double volumetric_heat_source;
    std::size_t cutbacks;
    int nonlinear_iterations;
    RegionInelasticSummary fuel_history;
    RegionInelasticSummary cladding_history;
};

struct TransientFuelCladdingResult final {
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

class TransientFuelCladdingTimeStepper final {
  public:
    TransientFuelCladdingResult
    solve(TransientFuelCladdingProblem& problem,
          const TransientTimeOptions& time_options,
          const SolverOptions& solver_options = SolverOptions{}) const;
};

} // namespace fuelsim

#endif
