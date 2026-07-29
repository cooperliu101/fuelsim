#ifndef FUELSIM_PETSC_SOLVER_HPP
#define FUELSIM_PETSC_SOLVER_HPP

#include <string>
#include <vector>

#include "fuelsim/nonlinear_problem.hpp"

namespace fuelsim {

class PetscSession final {
  public:
    PetscSession(int& argc, char**& argv, const char* help);
    ~PetscSession();

    PetscSession(const PetscSession&) = delete;
    PetscSession& operator=(const PetscSession&) = delete;

  private:
    bool _owns_initialization;
};

struct SolverOptions final {
    double absolute_tolerance = 1.0e-8;
    double relative_tolerance = 1.0e-10;
    double step_tolerance = 1.0e-12;
    int maximum_iterations = 40;
};

struct SolveResult final {
    std::vector<double> state;
    int nonlinear_iterations;
    double residual_norm;
    int convergence_reason;
    bool converged;
};

class PetscSequentialSolver final {
  public:
    SolveResult solve(const NonlinearProblem& problem,
                      const std::vector<double>& initial_state,
                      const SolverOptions& options = SolverOptions{}) const;
};

std::string petsc_convergence_reason_name(int reason);

} // namespace fuelsim

#endif
