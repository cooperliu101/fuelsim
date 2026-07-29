#ifndef FUELSIM_M1_SOLVER_HPP
#define FUELSIM_M1_SOLVER_HPP

#include <cstddef>

#include "fuelsim/m1_problem.hpp"
#include "fuelsim/petsc_solver.hpp"

namespace fuelsim {

struct M1LoadStepResult final {
    SolveResult solve;
    std::size_t completed_steps;
    bool completed;
};

class M1LoadStepper final {
  public:
    M1LoadStepResult
    solve(const M1Parameters& target_parameters, std::size_t load_steps,
          const SolverOptions& options = SolverOptions{}) const;
};

} // namespace fuelsim

#endif
