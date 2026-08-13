#ifndef FUELSIM_PETSC_SOLVER_HPP
#define FUELSIM_PETSC_SOLVER_HPP
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
    int _rank;
    int _size;
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
    double absolute_tolerance = 1.0e-8;
    double relative_tolerance = 1.0e-10;
    double step_tolerance = 1.0e-12;
    int maximum_iterations = 40;
    LineSearch line_search = LineSearch::basic;
    LinearSolver linear_solver = LinearSolver::automatic;
    Preconditioner preconditioner = Preconditioner::automatic;
    double linear_relative_tolerance = 1.0e-8;
    int maximum_linear_iterations = 500;
    bool backtracking_fallback = true;
    bool field_residual_scaling = false;
    double residual_reduction_tolerance = 1.0e-6;
    double temperature_residual_absolute_tolerance = 1.0e-8;
    double mechanical_residual_absolute_tolerance = 1.0e-4;
    double temperature_residual_scale = 0.0;
    double mechanical_residual_scale = 0.0;
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
    int nonlinear_iterations = 0;
    int linear_iterations = 0;
    double residual_norm = 0.0;
    int convergence_reason = 0;
    bool converged = false;
    SolveTiming timing;
    int mpi_rank = 0;
    int mpi_size = 1;
    std::size_t local_contribution_begin = 0;
    std::size_t local_contribution_end = 0;
    std::size_t global_state_dofs = 0;
    std::size_t maximum_shadow_state_dofs = 0;
    std::size_t total_shadow_state_dofs = 0;
    std::size_t total_remote_shadow_state_dofs = 0;
    SolveFailureCategory failure_category = SolveFailureCategory::none;
    std::string failure_message;
    std::size_t nonlinear_attempts = 1;
    std::size_t augmented_lagrangian_iterations = 0;
    double maximum_contact_penetration = 0.0;
    bool used_backtracking_fallback = false;
    SolveFailureCategory basic_failure_category = SolveFailureCategory::none;
    std::string basic_failure_message;
    std::vector<std::string> field_names;
    std::vector<double> initial_field_residual_norms;
    std::vector<double> field_residual_reference_norms;
    std::vector<double> final_field_residual_norms;
    std::vector<double> final_scaled_field_residual_norms;
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
    std::unique_ptr<Implementation> _implementation;
};
std::string petsc_convergence_reason_name(int reason);
const char* solve_failure_category_name(SolveFailureCategory category) noexcept;
} // namespace fuelsim
#endif
