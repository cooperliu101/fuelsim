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

PetscErrorCode collective_timing(const SolveTiming& local,
                                 SolveTiming& result) {
    PetscFunctionBeginUser;
    std::array<double, 5> local_seconds = {
        local.setup_seconds,
        local.nonlinear_solve_seconds,
        local.residual_callback_seconds,
        local.jacobian_callback_seconds,
        local.total_seconds,
    };
    std::array<double, 5> maximum_seconds{};
    PetscCallMPI(MPIU_Allreduce(
        local_seconds.data(), maximum_seconds.data(),
        static_cast<MPIU_Count>(local_seconds.size()), MPI_DOUBLE, MPI_MAX,
        PETSC_COMM_WORLD));
    std::array<PetscInt64, 4> local_counts = {
        static_cast<PetscInt64>(local.residual_evaluations),
        static_cast<PetscInt64>(local.jacobian_evaluations),
        static_cast<PetscInt64>(local.workspace_setups),
        static_cast<PetscInt64>(local.solve_calls),
    };
    std::array<PetscInt64, 4> maximum_counts{};
    PetscCallMPI(MPIU_Allreduce(
        local_counts.data(), maximum_counts.data(),
        static_cast<MPIU_Count>(local_counts.size()), MPIU_INT64, MPI_MAX,
        PETSC_COMM_WORLD));
    result.setup_seconds = maximum_seconds[0];
    result.nonlinear_solve_seconds = maximum_seconds[1];
    result.residual_callback_seconds = maximum_seconds[2];
    result.jacobian_callback_seconds = maximum_seconds[3];
    result.total_seconds = maximum_seconds[4];
    result.residual_evaluations =
        static_cast<std::size_t>(maximum_counts[0]);
    result.jacobian_evaluations =
        static_cast<std::size_t>(maximum_counts[1]);
    result.workspace_setups = static_cast<std::size_t>(maximum_counts[2]);
    result.solve_calls = static_cast<std::size_t>(maximum_counts[3]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscInt checked_petsc_int(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
        throw std::length_error("fuelsim DOF index exceeds PetscInt range");
    return static_cast<PetscInt>(value);
}

struct SolverContext final {
    const NonlinearProblem* problem = nullptr;
    bool pattern_locked = false;
    PetscMPIInt rank = 0;
    PetscMPIInt size = 1;
    std::size_t contribution_begin = 0;
    std::size_t contribution_end = 0;
    VecScatter state_scatter = nullptr;
    Vec gathered_state = nullptr;
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
    Vec gathered_state = nullptr;
    VecScatter state_scatter = nullptr;

    ~PetscObjects() {
        if (state_scatter != nullptr) {
            const PetscErrorCode code = VecScatterDestroy(&state_scatter);
            (void)code;
        }
        if (gathered_state != nullptr) {
            const PetscErrorCode code = VecDestroy(&gathered_state);
            (void)code;
        }
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

void configure_linear_solver(PetscObjects& objects,
                             const SolverOptions& options,
                             PetscMPIInt world_size) {
    SolverOptions::LinearSolver linear = options.linear_solver;
    if (linear == SolverOptions::LinearSolver::automatic) {
        if (options.preconditioner ==
                SolverOptions::Preconditioner::block_jacobi ||
            options.preconditioner ==
                SolverOptions::Preconditioner::field_split ||
            options.preconditioner == SolverOptions::Preconditioner::hypre)
            linear = SolverOptions::LinearSolver::gmres;
        else
            linear = SolverOptions::LinearSolver::direct;
    }
    SolverOptions::Preconditioner preconditioner_type =
        options.preconditioner;
    if (preconditioner_type == SolverOptions::Preconditioner::automatic)
        preconditioner_type =
            linear == SolverOptions::LinearSolver::direct
                ? SolverOptions::Preconditioner::lu
                : (world_size == 1
                       ? SolverOptions::Preconditioner::block_jacobi
                       : SolverOptions::Preconditioner::field_split);
    if (linear == SolverOptions::LinearSolver::direct &&
        preconditioner_type != SolverOptions::Preconditioner::lu)
        throw std::invalid_argument(
            "direct linear solver requires the LU preconditioner");

    KSP ksp = nullptr;
    PC preconditioner = nullptr;
    check_petsc(SNESGetKSP(objects.snes, &ksp), "SNESGetKSP");
    check_petsc(KSPGetPC(ksp, &preconditioner), "KSPGetPC");
    if (linear == SolverOptions::LinearSolver::direct)
        check_petsc(KSPSetType(ksp, KSPPREONLY), "KSPSetType PREONLY");
    else
        check_petsc(KSPSetType(ksp, KSPGMRES), "KSPSetType GMRES");

    switch (preconditioner_type) {
    case SolverOptions::Preconditioner::automatic:
        throw std::logic_error("automatic preconditioner was not resolved");
    case SolverOptions::Preconditioner::lu:
        check_petsc(PCSetType(preconditioner, PCLU), "PCSetType LU");
        if (world_size > 1)
            check_petsc(PCFactorSetMatSolverType(preconditioner,
                                                MATSOLVERMUMPS),
                        "PCFactorSetMatSolverType MUMPS");
        check_petsc(PCFactorSetReuseOrdering(preconditioner, PETSC_TRUE),
                    "PCFactorSetReuseOrdering");
        check_petsc(PCFactorSetReuseFill(preconditioner, PETSC_TRUE),
                    "PCFactorSetReuseFill");
        break;
    case SolverOptions::Preconditioner::block_jacobi:
        check_petsc(PCSetType(preconditioner, PCBJACOBI),
                    "PCSetType BJACOBI");
        break;
    case SolverOptions::Preconditioner::field_split: {
        PetscInt global_count = 0;
        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        check_petsc(VecGetSize(objects.state, &global_count),
                    "VecGetSize field split");
        if (global_count % 3 != 0)
            throw std::invalid_argument(
                "field_split requires the [T, ur, uz] three-field layout");
        check_petsc(VecGetOwnershipRange(objects.state, &ownership_begin,
                                         &ownership_end),
                    "VecGetOwnershipRange field split");
        const PetscInt node_count = global_count / 3;
        const PetscInt temperature_begin =
            std::max<PetscInt>(ownership_begin, 0);
        const PetscInt temperature_end =
            std::min<PetscInt>(ownership_end, node_count);
        const PetscInt mechanics_begin =
            std::max<PetscInt>(ownership_begin, node_count);
        const PetscInt mechanics_end =
            std::min<PetscInt>(ownership_end, global_count);
        IS temperature = nullptr;
        IS mechanics = nullptr;
        const PetscInt temperature_count =
            std::max<PetscInt>(0, temperature_end - temperature_begin);
        const PetscInt mechanics_count =
            std::max<PetscInt>(0, mechanics_end - mechanics_begin);
        check_petsc(ISCreateStride(PETSC_COMM_WORLD,
                                   temperature_count,
                                   temperature_begin, 1, &temperature),
                    "ISCreateStride temperature");
        check_petsc(ISCreateStride(PETSC_COMM_WORLD,
                                   mechanics_count,
                                   mechanics_begin, 1, &mechanics),
                    "ISCreateStride mechanics");
        check_petsc(PCSetType(preconditioner, PCFIELDSPLIT),
                    "PCSetType FIELDSPLIT");
        check_petsc(PCFieldSplitSetIS(preconditioner, "temperature",
                                     temperature),
                    "PCFieldSplitSetIS temperature");
        check_petsc(PCFieldSplitSetIS(preconditioner, "mechanics", mechanics),
                    "PCFieldSplitSetIS mechanics");
        check_petsc(PCFieldSplitSetType(preconditioner,
                                       PC_COMPOSITE_MULTIPLICATIVE),
                    "PCFieldSplitSetType multiplicative");
        check_petsc(ISDestroy(&temperature), "ISDestroy temperature");
        check_petsc(ISDestroy(&mechanics), "ISDestroy mechanics");
        break;
    }
    case SolverOptions::Preconditioner::hypre:
        check_petsc(PCSetType(preconditioner, PCHYPRE), "PCSetType HYPRE");
        break;
    }
    check_petsc(KSPSetTolerances(ksp, options.linear_relative_tolerance,
                                 PETSC_DEFAULT, PETSC_DEFAULT,
                                 options.maximum_linear_iterations),
                "KSPSetTolerances");
}

PetscErrorCode gather_state(Vec state, SolverContext& context) {
    PetscFunctionBeginUser;
    PetscCall(VecScatterBegin(context.state_scatter, state,
                              context.gathered_state, INSERT_VALUES,
                              SCATTER_FORWARD));
    PetscCall(VecScatterEnd(context.state_scatter, state,
                            context.gathered_state, INSERT_VALUES,
                            SCATTER_FORWARD));
    const PetscScalar* values = nullptr;
    PetscCall(VecGetArrayRead(context.gathered_state, &values));
    for (std::size_t index = 0; index < context.state_values.size(); ++index)
        context.state_values[index] =
            PetscRealPart(values[checked_petsc_int(index)]);
    PetscCall(VecRestoreArrayRead(context.gathered_state, &values));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_function(SNES, Vec state, Vec residual, void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        PetscCall(gather_state(state, context));
        PetscCall(VecSet(residual, 0.0));
        for (std::size_t contribution = context.contribution_begin;
             contribution < context.contribution_end; ++contribution) {
            const LocalDofs size_dofs =
                problem.contribution_dofs(contribution);
            const LocalValues local_state =
                problem.contribution_state(contribution,
                                           context.state_values);
            const LocalResidual local_residual =
                problem.contribution_residual(contribution, local_state);
            std::array<PetscInt, local_dof_count> dofs{};
            for (std::size_t local = 0; local < dofs.size(); ++local)
                dofs[local] = checked_petsc_int(size_dofs[local]);
            PetscCall(VecSetValues(
                residual, static_cast<PetscInt>(dofs.size()), dofs.data(),
                local_residual.data(), ADD_VALUES));
        }
        if (context.rank == 0) {
            problem.assemble_state_independent_residual(
                context.residual_values);
            for (std::size_t index = 0;
                 index < context.residual_values.size(); ++index) {
                if (context.residual_values[index] != 0.0)
                    PetscCall(VecSetValue(
                        residual, checked_petsc_int(index),
                        context.residual_values[index], ADD_VALUES));
            }
        }
        PetscCall(VecAssemblyBegin(residual));
        PetscCall(VecAssemblyEnd(residual));

        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        PetscCall(VecGetOwnershipRange(residual, &ownership_begin,
                                       &ownership_end));
        PetscScalar* local_residual = nullptr;
        PetscCall(VecGetArray(residual, &local_residual));
        for (const DirichletCondition& condition :
             problem.dirichlet_conditions()) {
            const PetscInt dof = checked_petsc_int(condition.dof);
            if (dof >= ownership_begin && dof < ownership_end)
                local_residual[dof - ownership_begin] =
                    context.state_values[condition.dof] - condition.value;
        }
        PetscCall(VecRestoreArray(residual, &local_residual));
        context.timing.residual_callback_seconds += seconds_since(start);
        ++context.timing.residual_evaluations;
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
                "fuelsim residual assembly failed: %s", error.what());
    } catch (...) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
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
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
                    "fuelsim requires one matrix for J and P");

        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        PetscCall(gather_state(state, context));

        PetscCall(MatZeroEntries(jacobian));
        for (std::size_t contribution = context.contribution_begin;
             contribution < context.contribution_end; ++contribution) {
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
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
                "fuelsim Jacobian assembly failed: %s", error.what());
    } catch (...) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
                "fuelsim Jacobian assembly failed with unknown error");
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

class PetscSolver::Implementation final {
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
        _context.rank = PetscGlobalRank;
        _context.size = PetscGlobalSize;
        const std::size_t rank = static_cast<std::size_t>(_context.rank);
        const std::size_t size = static_cast<std::size_t>(_context.size);
        _context.contribution_begin = problem.contribution_count() * rank / size;
        _context.contribution_end =
            problem.contribution_count() * (rank + 1U) / size;
        _context.state_values.resize(problem.dof_count());
        _context.residual_values.resize(problem.dof_count());

        check_petsc(
            VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, _count,
                         &_objects->state),
            "VecCreateMPI state");
        check_petsc(VecDuplicate(_objects->state, &_objects->residual),
                    "VecDuplicate residual");
        check_petsc(VecScatterCreateToAll(
                        _objects->state, &_objects->state_scatter,
                        &_objects->gathered_state),
                    "VecScatterCreateToAll state");
        _context.state_scatter = _objects->state_scatter;
        _context.gathered_state = _objects->gathered_state;

        PetscInt local_count = 0;
        check_petsc(VecGetLocalSize(_objects->state, &local_count),
                    "VecGetLocalSize state");
        const PetscInt diagonal_nonzeros = std::min<PetscInt>(60, local_count);
        const PetscInt off_diagonal_nonzeros =
            std::min<PetscInt>(60, _count - local_count);
        check_petsc(MatCreateAIJ(
                        PETSC_COMM_WORLD, local_count, local_count, _count,
                        _count, diagonal_nonzeros, nullptr,
                        off_diagonal_nonzeros, nullptr,
                        &_objects->jacobian),
                    "MatCreateAIJ");
        check_petsc(MatSetOption(_objects->jacobian,
                                 MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
                    "MatSetOption MAT_NEW_NONZERO_ALLOCATION_ERR");
        check_petsc(SNESCreate(PETSC_COMM_WORLD, &_objects->snes),
                    "SNESCreate");
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

        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        check_petsc(VecGetOwnershipRange(_objects->state, &ownership_begin,
                                         &ownership_end),
                    "VecGetOwnershipRange state");
        for (const DirichletCondition& condition :
             problem.dirichlet_conditions()) {
            const PetscInt dof = checked_petsc_int(condition.dof);
            if (dof >= ownership_begin && dof < ownership_end)
                _context.constrained_dofs.push_back(dof);
        }
        return true;
    }

    PetscObjects& objects() {
        return *_objects;
    }

    SolverContext& context() noexcept {
        return _context;
    }

  private:
    const NonlinearProblem* _problem = nullptr;
    PetscInt _count = 0;
    SolverContext _context;
    std::unique_ptr<PetscObjects> _objects;
};

PetscSession::PetscSession(int& argc, char**& argv, const char* help)
    : _owns_initialization(false), _rank(0), _size(1) {
    PetscBool initialized = PETSC_FALSE;
    check_petsc(PetscInitialized(&initialized), "PetscInitialized");
    if (initialized == PETSC_FALSE) {
        check_petsc(PetscInitialize(&argc, &argv, nullptr, help),
                    "PetscInitialize");
        _owns_initialization = true;
    }
    _rank = static_cast<int>(PetscGlobalRank);
    _size = static_cast<int>(PetscGlobalSize);
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

int PetscSession::rank() const noexcept {
    return _rank;
}

int PetscSession::size() const noexcept {
    return _size;
}

PetscSolver::PetscSolver()
    : _implementation(std::make_unique<Implementation>()) {}

PetscSolver::~PetscSolver() = default;

SolveResult
PetscSolver::solve(const NonlinearProblem& problem,
                   const std::vector<double>& initial_state,
                   const SolverOptions& options) {
    if (initial_state.size() != problem.dof_count())
        throw std::invalid_argument(
            "PetscSolver initial state size mismatch");
    if (!(options.absolute_tolerance > 0.0) ||
        !(options.relative_tolerance > 0.0) ||
        !(options.step_tolerance > 0.0) || options.maximum_iterations <= 0 ||
        !(options.linear_relative_tolerance > 0.0) ||
        options.maximum_linear_iterations <= 0)
        throw std::invalid_argument(
            "PetscSolver tolerances and iterations must be positive");

    const SteadyClock::time_point total_start = SteadyClock::now();
    const SteadyClock::time_point setup_start = SteadyClock::now();
    const bool workspace_created = _implementation->prepare(problem);
    PetscObjects& objects = _implementation->objects();
    SolverContext& context = _implementation->context();
    context.timing = SolveTiming{};
    context.timing.workspace_setups = workspace_created ? 1U : 0U;
    context.timing.solve_calls = 1;
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    check_petsc(VecGetOwnershipRange(objects.state, &ownership_begin,
                                     &ownership_end),
                "VecGetOwnershipRange state");
    PetscScalar* state_array = nullptr;
    check_petsc(VecGetArray(objects.state, &state_array), "VecGetArray state");
    for (PetscInt index = ownership_begin; index < ownership_end; ++index)
        state_array[index - ownership_begin] =
            initial_state[static_cast<std::size_t>(index)];
    check_petsc(VecRestoreArray(objects.state, &state_array),
                "VecRestoreArray state");

    check_petsc(SNESSetTolerances(objects.snes, options.absolute_tolerance,
                                  options.relative_tolerance,
                                  options.step_tolerance,
                                  options.maximum_iterations, PETSC_DEFAULT),
                "SNESSetTolerances");
    configure_linear_solver(objects, options, PetscGlobalSize);
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

    check_petsc(gather_state(objects.state, context), "gather_state solution");
    std::vector<double> solution(problem.dof_count(), 0.0);
    solution = context.state_values;

    context.timing.total_seconds = seconds_since(total_start);
    SolveResult result;
    result.state = std::move(solution);
    result.nonlinear_iterations = static_cast<int>(iterations);
    result.residual_norm = static_cast<double>(residual_norm);
    result.convergence_reason = static_cast<int>(reason);
    result.converged = reason > 0;
    check_petsc(collective_timing(context.timing, result.timing),
                "collective timing reduction");
    result.mpi_rank = static_cast<int>(context.rank);
    result.mpi_size = static_cast<int>(context.size);
    result.local_contribution_begin = context.contribution_begin;
    result.local_contribution_end = context.contribution_end;
    return result;
}

std::string petsc_convergence_reason_name(int reason) {
    const char* name = SNESConvergedReasons[reason];
    if (name == nullptr)
        return "UNKNOWN";
    return name;
}

} // namespace fuelsim
