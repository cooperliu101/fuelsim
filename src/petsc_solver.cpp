#include "fuelsim/petsc_solver.hpp"

#include <petsc.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuelsim {
namespace {

using SteadyClock = std::chrono::steady_clock;

double seconds_since(const SteadyClock::time_point& start) {
    return std::chrono::duration<double>(SteadyClock::now() - start).count();
}

void check_petsc(PetscErrorCode code, const char* operation) {
    if (code == PETSC_SUCCESS)
        return;

    const char* text = nullptr;
    const PetscErrorCode message_code = PetscErrorMessage(code, &text, nullptr);
    if (message_code != PETSC_SUCCESS)
        text = nullptr;
    std::string message = operation;
    message += " failed";
    if (text != nullptr) {
        message += ": ";
        message += text;
    }
    throw std::runtime_error(message);
}

PetscInt checked_petsc_int(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
        throw std::length_error("fuelsim DOF index exceeds PetscInt range");
    return static_cast<PetscInt>(value);
}

struct SolverContext final {
    const NonlinearProblem* problem = nullptr;
    bool pattern_locked = false;
    std::vector<double> state_values;
    std::vector<double> residual_values;
    std::vector<PetscInt> constrained_dofs;
    SolveTiming timing;
};

struct PetscObjects final {
    SNES snes = nullptr;
    Vec state = nullptr;
    Vec residual = nullptr;
    Mat jacobian = nullptr;

    ~PetscObjects() {
        if (snes != nullptr) {
            const PetscErrorCode code = SNESDestroy(&snes);
            (void)code;
        }
        if (state != nullptr) {
            const PetscErrorCode code = VecDestroy(&state);
            (void)code;
        }
        if (residual != nullptr) {
            const PetscErrorCode code = VecDestroy(&residual);
            (void)code;
        }
        if (jacobian != nullptr) {
            const PetscErrorCode code = MatDestroy(&jacobian);
            (void)code;
        }
    }
};

PetscErrorCode form_function(SNES, Vec state, Vec residual, void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        const PetscInt count = checked_petsc_int(problem.dof_count());

        const PetscScalar* petsc_state = nullptr;
        PetscCall(VecGetArrayRead(state, &petsc_state));
        context.state_values.resize(problem.dof_count());
        for (PetscInt index = 0; index < count; ++index) {
            context.state_values[static_cast<std::size_t>(index)] =
                PetscRealPart(petsc_state[index]);
        }
        PetscCall(VecRestoreArrayRead(state, &petsc_state));

        problem.assemble_residual(context.state_values,
                                  context.residual_values);
        for (const DirichletCondition& condition :
             problem.dirichlet_conditions()) {
            context.residual_values[condition.dof] =
                context.state_values[condition.dof] - condition.value;
        }

        PetscScalar* petsc_residual = nullptr;
        PetscCall(VecGetArray(residual, &petsc_residual));
        for (PetscInt index = 0; index < count; ++index) {
            petsc_residual[index] =
                context.residual_values[static_cast<std::size_t>(index)];
        }
        PetscCall(VecRestoreArray(residual, &petsc_residual));
        context.timing.residual_callback_seconds += seconds_since(start);
        ++context.timing.residual_evaluations;
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,
                "fuelsim residual assembly failed: %s", error.what());
    } catch (...) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,
                "fuelsim residual assembly failed with unknown error");
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_jacobian(SNES, Vec state, Mat jacobian, Mat preconditioner,
                             void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        if (jacobian != preconditioner)
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
                    "fuelsim requires one matrix for J and P");

        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        const PetscInt count = checked_petsc_int(problem.dof_count());

        const PetscScalar* petsc_state = nullptr;
        PetscCall(VecGetArrayRead(state, &petsc_state));
        context.state_values.resize(problem.dof_count());
        for (PetscInt index = 0; index < count; ++index) {
            context.state_values[static_cast<std::size_t>(index)] =
                PetscRealPart(petsc_state[index]);
        }
        PetscCall(VecRestoreArrayRead(state, &petsc_state));

        PetscCall(MatZeroEntries(jacobian));
        for (std::size_t contribution = 0;
             contribution < problem.contribution_count(); ++contribution) {
            const LocalValues local_state =
                problem.contribution_state(contribution, context.state_values);
            const LocalSystem local_system =
                problem.linearize_contribution(contribution, local_state);
            const LocalDofs size_dofs = problem.contribution_dofs(contribution);
            std::array<PetscInt, local_dof_count> dofs{};
            for (std::size_t local = 0; local < dofs.size(); ++local)
                dofs[local] = checked_petsc_int(size_dofs[local]);

            PetscCall(MatSetValues(
                jacobian, static_cast<PetscInt>(dofs.size()), dofs.data(),
                static_cast<PetscInt>(dofs.size()), dofs.data(),
                local_system.jacobian.data(), ADD_VALUES));
        }

        PetscCall(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));

        if (!context.pattern_locked) {
            PetscCall(
                MatSetOption(jacobian, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
            PetscCall(
                MatSetOption(jacobian, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE));
            context.pattern_locked = true;
        }

        PetscCall(MatZeroRows(
            jacobian, static_cast<PetscInt>(context.constrained_dofs.size()),
            context.constrained_dofs.data(), 1.0, nullptr, nullptr));
        context.timing.jacobian_callback_seconds += seconds_since(start);
        ++context.timing.jacobian_evaluations;
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,
                "fuelsim Jacobian assembly failed: %s", error.what());
    } catch (...) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB,
                "fuelsim Jacobian assembly failed with unknown error");
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

class PetscSequentialSolver::Implementation final {
  public:
    bool prepare(const NonlinearProblem& problem) {
        const PetscInt requested_count = checked_petsc_int(problem.dof_count());
        if (_problem == &problem && _count == requested_count)
            return false;

        _objects = std::make_unique<PetscObjects>();
        _problem = &problem;
        _count = requested_count;
        _context = SolverContext{};
        _context.problem = &problem;
        _context.state_values.resize(problem.dof_count());
        _context.residual_values.resize(problem.dof_count());
        _context.constrained_dofs.reserve(
            problem.dirichlet_conditions().size());
        for (const DirichletCondition& condition :
             problem.dirichlet_conditions()) {
            _context.constrained_dofs.push_back(
                checked_petsc_int(condition.dof));
        }

        check_petsc(VecCreateSeq(PETSC_COMM_SELF, _count, &_objects->state),
                    "VecCreateSeq state");
        check_petsc(VecDuplicate(_objects->state, &_objects->residual),
                    "VecDuplicate residual");
        check_petsc(MatCreateSeqAIJ(PETSC_COMM_SELF, _count, _count, 40,
                                    nullptr, &_objects->jacobian),
                    "MatCreateSeqAIJ");
        check_petsc(MatSetOption(_objects->jacobian,
                                 MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
                    "MatSetOption MAT_NEW_NONZERO_ALLOCATION_ERR");
        check_petsc(SNESCreate(PETSC_COMM_SELF, &_objects->snes), "SNESCreate");
        check_petsc(SNESSetFunction(_objects->snes, _objects->residual,
                                    form_function, &_context),
                    "SNESSetFunction");
        check_petsc(SNESSetJacobian(_objects->snes, _objects->jacobian,
                                    _objects->jacobian, form_jacobian,
                                    &_context),
                    "SNESSetJacobian");
        check_petsc(SNESSetType(_objects->snes, SNESNEWTONLS), "SNESSetType");

        SNESLineSearch line_search = nullptr;
        check_petsc(SNESGetLineSearch(_objects->snes, &line_search),
                    "SNESGetLineSearch");
        check_petsc(SNESLineSearchSetType(line_search, SNESLINESEARCHBASIC),
                    "SNESLineSearchSetType");

        KSP ksp = nullptr;
        PC preconditioner = nullptr;
        check_petsc(SNESGetKSP(_objects->snes, &ksp), "SNESGetKSP");
        check_petsc(KSPSetType(ksp, KSPPREONLY), "KSPSetType");
        check_petsc(KSPGetPC(ksp, &preconditioner), "KSPGetPC");
        check_petsc(PCSetType(preconditioner, PCLU), "PCSetType");
        check_petsc(PCFactorSetReuseOrdering(preconditioner, PETSC_TRUE),
                    "PCFactorSetReuseOrdering");
        check_petsc(PCFactorSetReuseFill(preconditioner, PETSC_TRUE),
                    "PCFactorSetReuseFill");
        return true;
    }

    PetscObjects& objects() {
        return *_objects;
    }

    SolverContext& context() noexcept {
        return _context;
    }

    PetscInt count() const noexcept {
        return _count;
    }

  private:
    const NonlinearProblem* _problem = nullptr;
    PetscInt _count = 0;
    SolverContext _context;
    std::unique_ptr<PetscObjects> _objects;
};

PetscSession::PetscSession(int& argc, char**& argv, const char* help)
    : _owns_initialization(false) {
    PetscBool initialized = PETSC_FALSE;
    check_petsc(PetscInitialized(&initialized), "PetscInitialized");
    if (initialized == PETSC_FALSE) {
        check_petsc(PetscInitialize(&argc, &argv, nullptr, help),
                    "PetscInitialize");
        _owns_initialization = true;
    }
}

PetscSession::~PetscSession() {
    if (!_owns_initialization)
        return;

    PetscBool finalized = PETSC_FALSE;
    if (PetscFinalized(&finalized) == PETSC_SUCCESS &&
        finalized == PETSC_FALSE) {
        const PetscErrorCode code = PetscFinalize();
        (void)code;
    }
}

PetscSequentialSolver::PetscSequentialSolver()
    : _implementation(std::make_unique<Implementation>()) {}

PetscSequentialSolver::~PetscSequentialSolver() = default;

SolveResult
PetscSequentialSolver::solve(const NonlinearProblem& problem,
                             const std::vector<double>& initial_state,
                             const SolverOptions& options) {
    if (initial_state.size() != problem.dof_count())
        throw std::invalid_argument(
            "PetscSequentialSolver initial state size mismatch");
    if (!(options.absolute_tolerance > 0.0) ||
        !(options.relative_tolerance > 0.0) ||
        !(options.step_tolerance > 0.0) || options.maximum_iterations <= 0)
        throw std::invalid_argument(
            "PetscSequentialSolver tolerances and iterations must be positive");

    PetscMPIInt world_size = 0;
    const int mpi_code = MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
    if (mpi_code != MPI_SUCCESS)
        throw std::runtime_error("MPI_Comm_size failed");
    if (world_size != 1)
        throw std::invalid_argument(
            "PetscSequentialSolver supports exactly one MPI rank");

    const SteadyClock::time_point total_start = SteadyClock::now();
    const SteadyClock::time_point setup_start = SteadyClock::now();
    const bool workspace_created = _implementation->prepare(problem);
    PetscObjects& objects = _implementation->objects();
    SolverContext& context = _implementation->context();
    context.timing = SolveTiming{};
    context.timing.workspace_setups = workspace_created ? 1U : 0U;
    context.timing.solve_calls = 1;
    const PetscInt count = _implementation->count();

    PetscScalar* state_array = nullptr;
    check_petsc(VecGetArray(objects.state, &state_array), "VecGetArray state");
    for (PetscInt index = 0; index < count; ++index) {
        state_array[index] = initial_state[static_cast<std::size_t>(index)];
    }
    check_petsc(VecRestoreArray(objects.state, &state_array),
                "VecRestoreArray state");

    check_petsc(SNESSetTolerances(objects.snes, options.absolute_tolerance,
                                  options.relative_tolerance,
                                  options.step_tolerance,
                                  options.maximum_iterations, PETSC_DEFAULT),
                "SNESSetTolerances");
    check_petsc(SNESSetFromOptions(objects.snes), "SNESSetFromOptions");
    context.timing.setup_seconds = seconds_since(setup_start);

    const SteadyClock::time_point solve_start = SteadyClock::now();
    check_petsc(SNESSolve(objects.snes, nullptr, objects.state), "SNESSolve");
    context.timing.nonlinear_solve_seconds = seconds_since(solve_start);

    SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
    PetscInt iterations = 0;
    PetscReal residual_norm = 0.0;
    check_petsc(SNESGetConvergedReason(objects.snes, &reason),
                "SNESGetConvergedReason");
    check_petsc(SNESGetIterationNumber(objects.snes, &iterations),
                "SNESGetIterationNumber");
    check_petsc(SNESGetFunctionNorm(objects.snes, &residual_norm),
                "SNESGetFunctionNorm");

    const PetscScalar* solved_array = nullptr;
    check_petsc(VecGetArrayRead(objects.state, &solved_array),
                "VecGetArrayRead solution");
    std::vector<double> solution(problem.dof_count(), 0.0);
    for (PetscInt index = 0; index < count; ++index) {
        solution[static_cast<std::size_t>(index)] =
            PetscRealPart(solved_array[index]);
    }
    check_petsc(VecRestoreArrayRead(objects.state, &solved_array),
                "VecRestoreArrayRead solution");

    context.timing.total_seconds = seconds_since(total_start);
    SolveResult result;
    result.state = std::move(solution);
    result.nonlinear_iterations = static_cast<int>(iterations);
    result.residual_norm = static_cast<double>(residual_norm);
    result.convergence_reason = static_cast<int>(reason);
    result.converged = reason > 0;
    result.timing = context.timing;
    return result;
}

std::string petsc_convergence_reason_name(int reason) {
    const char* name = SNESConvergedReasons[reason];
    if (name == nullptr)
        return "UNKNOWN";
    return name;
}

} // namespace fuelsim
