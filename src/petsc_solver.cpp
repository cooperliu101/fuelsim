#include "fuelsim/petsc_solver.hpp"

#include <petsc.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuelsim {
namespace {

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
    const NonlinearProblem* problem;
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
    try {
        const SolverContext& context =
            *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        const PetscInt count = checked_petsc_int(problem.dof_count());

        const PetscScalar* petsc_state = nullptr;
        PetscCall(VecGetArrayRead(state, &petsc_state));
        std::vector<double> state_values(problem.dof_count(), 0.0);
        for (PetscInt index = 0; index < count; ++index) {
            state_values[static_cast<std::size_t>(index)] =
                PetscRealPart(petsc_state[index]);
        }
        PetscCall(VecRestoreArrayRead(state, &petsc_state));

        std::vector<double> residual_values;
        problem.assemble_residual(state_values, residual_values);
        for (const DirichletCondition& condition :
             problem.dirichlet_conditions()) {
            residual_values[condition.dof] =
                state_values[condition.dof] - condition.value;
        }

        PetscScalar* petsc_residual = nullptr;
        PetscCall(VecGetArray(residual, &petsc_residual));
        for (PetscInt index = 0; index < count; ++index) {
            petsc_residual[index] =
                residual_values[static_cast<std::size_t>(index)];
        }
        PetscCall(VecRestoreArray(residual, &petsc_residual));
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
    try {
        if (jacobian != preconditioner)
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
                    "fuelsim requires one matrix for J and P");

        const SolverContext& context =
            *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        const PetscInt count = checked_petsc_int(problem.dof_count());

        const PetscScalar* petsc_state = nullptr;
        PetscCall(VecGetArrayRead(state, &petsc_state));
        std::vector<double> state_values(problem.dof_count(), 0.0);
        for (PetscInt index = 0; index < count; ++index) {
            state_values[static_cast<std::size_t>(index)] =
                PetscRealPart(petsc_state[index]);
        }
        PetscCall(VecRestoreArrayRead(state, &petsc_state));

        PetscCall(MatZeroEntries(jacobian));
        for (std::size_t contribution = 0;
             contribution < problem.contribution_count(); ++contribution) {
            const LocalValues local_state =
                problem.contribution_state(contribution, state_values);
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

        std::vector<PetscInt> constrained;
        constrained.reserve(problem.dirichlet_conditions().size());
        for (const DirichletCondition& condition :
             problem.dirichlet_conditions()) {
            constrained.push_back(checked_petsc_int(condition.dof));
        }
        PetscCall(MatZeroRows(jacobian,
                              static_cast<PetscInt>(constrained.size()),
                              constrained.data(), 1.0, nullptr, nullptr));
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

PetscSession::PetscSession(int& argc, char**& argv, const char* help)
    : owns_initialization_(false) {
    PetscBool initialized = PETSC_FALSE;
    check_petsc(PetscInitialized(&initialized), "PetscInitialized");
    if (initialized == PETSC_FALSE) {
        check_petsc(PetscInitialize(&argc, &argv, nullptr, help),
                    "PetscInitialize");
        owns_initialization_ = true;
    }
}

PetscSession::~PetscSession() {
    if (!owns_initialization_)
        return;

    PetscBool finalized = PETSC_FALSE;
    if (PetscFinalized(&finalized) == PETSC_SUCCESS &&
        finalized == PETSC_FALSE) {
        const PetscErrorCode code = PetscFinalize();
        (void)code;
    }
}

SolveResult
PetscSequentialSolver::solve(const NonlinearProblem& problem,
                             const std::vector<double>& initial_state,
                             const SolverOptions& options) const {
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

    const PetscInt count = checked_petsc_int(problem.dof_count());
    PetscObjects objects;
    check_petsc(VecCreateSeq(PETSC_COMM_SELF, count, &objects.state),
                "VecCreateSeq state");
    check_petsc(VecDuplicate(objects.state, &objects.residual),
                "VecDuplicate residual");
    check_petsc(MatCreateSeqAIJ(PETSC_COMM_SELF, count, count, 40, nullptr,
                                &objects.jacobian),
                "MatCreateSeqAIJ");
    check_petsc(MatSetOption(objects.jacobian, MAT_NEW_NONZERO_ALLOCATION_ERR,
                             PETSC_FALSE),
                "MatSetOption");
    check_petsc(SNESCreate(PETSC_COMM_SELF, &objects.snes), "SNESCreate");

    PetscScalar* state_array = nullptr;
    check_petsc(VecGetArray(objects.state, &state_array), "VecGetArray state");
    for (PetscInt index = 0; index < count; ++index) {
        state_array[index] = initial_state[static_cast<std::size_t>(index)];
    }
    check_petsc(VecRestoreArray(objects.state, &state_array),
                "VecRestoreArray state");

    SolverContext context{&problem};
    check_petsc(SNESSetFunction(objects.snes, objects.residual, form_function,
                                &context),
                "SNESSetFunction");
    check_petsc(SNESSetJacobian(objects.snes, objects.jacobian,
                                objects.jacobian, form_jacobian, &context),
                "SNESSetJacobian");
    check_petsc(SNESSetType(objects.snes, SNESNEWTONLS), "SNESSetType");
    SNESLineSearch line_search = nullptr;
    check_petsc(SNESGetLineSearch(objects.snes, &line_search),
                "SNESGetLineSearch");
    check_petsc(SNESLineSearchSetType(line_search, SNESLINESEARCHBASIC),
                "SNESLineSearchSetType");
    check_petsc(SNESSetTolerances(objects.snes, options.absolute_tolerance,
                                  options.relative_tolerance,
                                  options.step_tolerance,
                                  options.maximum_iterations, PETSC_DEFAULT),
                "SNESSetTolerances");

    KSP ksp = nullptr;
    PC preconditioner = nullptr;
    check_petsc(SNESGetKSP(objects.snes, &ksp), "SNESGetKSP");
    check_petsc(KSPSetType(ksp, KSPPREONLY), "KSPSetType");
    check_petsc(KSPGetPC(ksp, &preconditioner), "KSPGetPC");
    check_petsc(PCSetType(preconditioner, PCLU), "PCSetType");
    check_petsc(SNESSetFromOptions(objects.snes), "SNESSetFromOptions");

    check_petsc(SNESSolve(objects.snes, nullptr, objects.state), "SNESSolve");

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

    return {
        std::move(solution),
        static_cast<int>(iterations),
        static_cast<double>(residual_norm),
        static_cast<int>(reason),
        reason > 0,
    };
}

std::string petsc_convergence_reason_name(int reason) {
    const char* name = SNESConvergedReasons[reason];
    if (name == nullptr)
        return "UNKNOWN";
    return name;
}

} // namespace fuelsim
