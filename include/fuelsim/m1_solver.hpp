#ifndef FUELSIM_M1_SOLVER_HPP
#define FUELSIM_M1_SOLVER_HPP

#include <cstddef>

#include "fuelsim/m1_problem.hpp"
#include "fuelsim/petsc_solver.hpp"

namespace fuelsim {

struct M1LoadStepResult final {
    SolveResult solve;
    std::size_t completed_steps = 0;
    bool completed = false;
    int total_nonlinear_iterations = 0;
    double problem_setup_seconds = 0.0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
};

class M1LoadStepper final {
  public:
    M1LoadStepResult
    solve(const M1Parameters& target_parameters, std::size_t load_steps,
          const SolverOptions& options = SolverOptions{}) const;
};

} // namespace fuelsim

#endif
