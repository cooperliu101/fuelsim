#ifndef FUELSIM_STEADY_FUEL_CLADDING_SOLVER_HPP
#define FUELSIM_STEADY_FUEL_CLADDING_SOLVER_HPP

#include <cstddef>

#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/steady_fuel_cladding_problem.hpp"

namespace fuelsim {

struct SteadyFuelCladdingLoadResult final {
    SolveResult solve;
    std::size_t completed_steps = 0;
    bool completed = false;
    int total_nonlinear_iterations = 0;
    double problem_setup_seconds = 0.0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
};

class SteadyFuelCladdingLoadStepper final {
  public:
    SteadyFuelCladdingLoadResult
    solve(const SteadyFuelCladdingParameters& target_parameters,
          std::size_t load_steps,
          const SolverOptions& options = SolverOptions{}) const;
    SteadyFuelCladdingLoadResult
    solve(const SteadyFuelCladdingParameters& target_parameters,
          StructuredRzMesh fuel_mesh, StructuredRzMesh cladding_mesh,
          std::size_t load_steps,
          const SolverOptions& options = SolverOptions{}) const;
};

} // namespace fuelsim

#endif
