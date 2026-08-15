#pragma once
#include "fuelsim/nonlinear_problem.hpp"
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace fuelsim {
class PetscSession final {
  public:
    PetscSession(int& argc, char**& argv, const char* help);
    ~PetscSession();
    PetscSession(const PetscSession&) = delete;
    PetscSession& operator=(const PetscSession&) = delete;
    int rank() const noexcept;
    int size() const noexcept;
    void collective_root_action(const std::function<void()>& action) const;

  private:
    bool _owns_initialization;
    int _rank, _size;
};

struct SolverOptions final {
    enum class LineSearch {
        backtracking,
        basic,
    };
    enum class LinearSolver {
        automatic,
        direct,
        gmres,
    };
    enum class Preconditioner {
        automatic,
        lu,
        block_jacobi,
        field_split,
        hypre,
    };
    double absolute_tolerance = 1.0e-8, relative_tolerance = 1.0e-10, step_tolerance = 1.0e-12;
    int maximum_iterations = 40;
    LineSearch line_search = LineSearch::basic;
    LinearSolver linear_solver = LinearSolver::automatic;
    Preconditioner preconditioner = Preconditioner::automatic;
    double linear_relative_tolerance = 1.0e-8;
    int maximum_linear_iterations = 500;
    bool backtracking_fallback = true, field_residual_scaling = false;
    double residual_reduction_tolerance = 1.0e-6, temperature_residual_absolute_tolerance = 1.0e-8,
           mechanical_residual_absolute_tolerance = 1.0e-4, temperature_residual_scale = 0.0,
           mechanical_residual_scale = 0.0;
};

struct SolveTiming final {
    double setup_seconds = 0.0, nonlinear_solve_seconds = 0.0, residual_callback_seconds = 0.0,
           jacobian_callback_seconds = 0.0, total_seconds = 0.0;
    std::size_t residual_evaluations = 0, jacobian_evaluations = 0, workspace_setups = 0, solve_calls = 0;
};
enum class SolveFailureCategory {
    none,
    nonlinear_divergence,
    physical_domain,
    residual_verification,
    time_discretization,
    contact_constraint,
};

struct SolveResult final {
    std::vector<double> state;
    int nonlinear_iterations = 0, linear_iterations = 0;
    double residual_norm = 0.0;
    int convergence_reason = 0;
    bool converged = false;
    SolveTiming timing;
    int mpi_rank = 0, mpi_size = 1;
    std::size_t local_contribution_begin = 0, local_contribution_end = 0, global_state_dofs = 0,
                maximum_shadow_state_dofs = 0, total_shadow_state_dofs = 0, total_remote_shadow_state_dofs = 0;
    SolveFailureCategory failure_category = SolveFailureCategory::none;
    std::string failure_message;
    std::size_t nonlinear_attempts = 1, augmented_lagrangian_iterations = 0;
    double maximum_contact_penetration = 0.0;
    bool used_backtracking_fallback = false;
    SolveFailureCategory basic_failure_category = SolveFailureCategory::none;
    std::string basic_failure_message;
    std::vector<std::string> field_names;
    std::vector<double> initial_field_residual_norms, field_residual_reference_norms;
    std::vector<double> final_field_residual_norms, final_scaled_field_residual_norms;
    std::vector<double> field_residual_scalings;
};

class PetscSolver final {
  public:
    PetscSolver();
    ~PetscSolver();
    PetscSolver(const PetscSolver&) = delete;
    PetscSolver& operator=(const PetscSolver&) = delete;
    SolveResult solve(const NonlinearProblem& problem, const std::vector<double>& initial_state,
        const SolverOptions& options = SolverOptions{});

  private:
    SolveResult solve_once(
        const NonlinearProblem& problem, const std::vector<double>& initial_state, const SolverOptions& options);
    class Implementation;
    std::unique_ptr<Implementation> _impl;
};

std::string petsc_convergence_reason_name(int reason);
const char* solve_failure_category_name(SolveFailureCategory category) noexcept;
} // namespace fuelsim
