#ifndef FUELSIM_PETSC_SOLVER_HPP
#define FUELSIM_PETSC_SOLVER_HPP

#include <cstddef>
#include <memory>
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

struct SolveTiming final {
    double setup_seconds = 0.0;
    double nonlinear_solve_seconds = 0.0;
    double residual_callback_seconds = 0.0;
    double jacobian_callback_seconds = 0.0;
    double total_seconds = 0.0;
    std::size_t residual_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t workspace_setups = 0;
    std::size_t solve_calls = 0;
};

struct SolveResult final {
    std::vector<double> state;
    int nonlinear_iterations = 0;
    double residual_norm = 0.0;
    int convergence_reason = 0;
    bool converged = false;
    SolveTiming timing;
};

class PetscSequentialSolver final {
  public:
    PetscSequentialSolver();
    ~PetscSequentialSolver();

    PetscSequentialSolver(const PetscSequentialSolver&) = delete;
    PetscSequentialSolver& operator=(const PetscSequentialSolver&) = delete;

    SolveResult solve(const NonlinearProblem& problem,
                      const std::vector<double>& initial_state,
                      const SolverOptions& options = SolverOptions{});

  private:
    class Implementation;
    std::unique_ptr<Implementation> _implementation;
};

std::string petsc_convergence_reason_name(int reason);

} // namespace fuelsim

#endif
