#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem_solver.hpp"

#include <petsc.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {

using SteadyClock = std::chrono::steady_clock;

double seconds_since(const SteadyClock::time_point& start) {
    return std::chrono::duration<double>(SteadyClock::now() - start).count();
}

void accumulate_timing(SolveTiming& total, const SolveTiming& step) {
    total.setup_seconds += step.setup_seconds;
    total.nonlinear_solve_seconds += step.nonlinear_solve_seconds;
    total.residual_callback_seconds += step.residual_callback_seconds;
    total.jacobian_callback_seconds += step.jacobian_callback_seconds;
    total.total_seconds += step.total_seconds;
    total.residual_evaluations += step.residual_evaluations;
    total.jacobian_evaluations += step.jacobian_evaluations;
    total.workspace_setups += step.workspace_setups;
    total.solve_calls += step.solve_calls;
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

void check_mpi(PetscMPIInt code, const char* operation) {
    if (code == MPI_SUCCESS)
        return;
    std::string message = operation;
    message += " failed with MPI error code ";
    message += std::to_string(code);
    throw std::runtime_error(message);
}

PetscErrorCode collective_timing(const SolveTiming& local, SolveTiming& result) {
    PetscFunctionBeginUser;
    std::array<double, 5> local_seconds = {
        local.setup_seconds,
        local.nonlinear_solve_seconds,
        local.residual_callback_seconds,
        local.jacobian_callback_seconds,
        local.total_seconds,
    };
    std::array<double, 5> maximum_seconds{};
    PetscCallMPI(MPIU_Allreduce(local_seconds.data(), maximum_seconds.data(),
                                static_cast<MPIU_Count>(local_seconds.size()), MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD));
    std::array<PetscInt64, 4> local_counts = {
        static_cast<PetscInt64>(local.residual_evaluations),
        static_cast<PetscInt64>(local.jacobian_evaluations),
        static_cast<PetscInt64>(local.workspace_setups),
        static_cast<PetscInt64>(local.solve_calls),
    };
    std::array<PetscInt64, 4> maximum_counts{};
    PetscCallMPI(MPIU_Allreduce(local_counts.data(), maximum_counts.data(),
                                static_cast<MPIU_Count>(local_counts.size()), MPIU_INT64, MPI_MAX, PETSC_COMM_WORLD));
    result.setup_seconds = maximum_seconds[0];
    result.nonlinear_solve_seconds = maximum_seconds[1];
    result.residual_callback_seconds = maximum_seconds[2];
    result.jacobian_callback_seconds = maximum_seconds[3];
    result.total_seconds = maximum_seconds[4];
    result.residual_evaluations = static_cast<std::size_t>(maximum_counts[0]);
    result.jacobian_evaluations = static_cast<std::size_t>(maximum_counts[1]);
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
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    VecScatter state_scatter = nullptr;
    Vec gathered_state = nullptr;
    std::vector<std::uint32_t> shadow_dofs;
    std::vector<double> state_values;
    std::vector<PetscInt> constrained_dofs;
    std::vector<bool> constrained;
    bool field_residual_scaling = true;
    double residual_scaling_floor = 1.0e-8;
    bool first_residual = true;
    bool thermal_scaling_initialized = false;
    bool mechanics_scaling_initialized = false;
    std::array<double, 3> initial_field_residual_norms{};
    std::array<double, 3> field_residual_reference_norms{};
    std::array<double, 3> latest_unscaled_field_residual_norms{};
    std::array<double, 3> latest_field_residual_norms{};
    std::array<double, 3> field_residual_scalings{{1.0, 1.0, 1.0}};
    double initial_residual_norm = std::numeric_limits<double>::quiet_NaN();
    bool saw_domain_error = false;
    bool last_function_domain_error = false;
    std::string last_domain_error;
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

void configure_linear_solver(PetscObjects& objects, const SolverOptions& options, PetscMPIInt world_size) {
    SolverOptions::LinearSolver linear = options.linear_solver;
    if (linear == SolverOptions::LinearSolver::automatic) {
        if (options.preconditioner == SolverOptions::Preconditioner::block_jacobi ||
            options.preconditioner == SolverOptions::Preconditioner::field_split ||
            options.preconditioner == SolverOptions::Preconditioner::hypre)
            linear = SolverOptions::LinearSolver::gmres;
        else
            linear = SolverOptions::LinearSolver::direct;
    }
    SolverOptions::Preconditioner preconditioner_type = options.preconditioner;
    if (preconditioner_type == SolverOptions::Preconditioner::automatic)
        preconditioner_type = linear == SolverOptions::LinearSolver::direct
                                  ? SolverOptions::Preconditioner::lu
                                  : (world_size == 1 ? SolverOptions::Preconditioner::block_jacobi
                                                     : SolverOptions::Preconditioner::field_split);
    if (linear == SolverOptions::LinearSolver::direct && preconditioner_type != SolverOptions::Preconditioner::lu)
        throw std::invalid_argument("direct linear solver requires the LU preconditioner");

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
            check_petsc(PCFactorSetMatSolverType(preconditioner, MATSOLVERMUMPS), "PCFactorSetMatSolverType MUMPS");
        check_petsc(PCFactorSetReuseOrdering(preconditioner, PETSC_TRUE), "PCFactorSetReuseOrdering");
        check_petsc(PCFactorSetReuseFill(preconditioner, PETSC_TRUE), "PCFactorSetReuseFill");
        break;
    case SolverOptions::Preconditioner::block_jacobi:
        check_petsc(PCSetType(preconditioner, PCBJACOBI), "PCSetType BJACOBI");
        break;
    case SolverOptions::Preconditioner::field_split: {
        PetscBool already_configured = PETSC_FALSE;
        check_petsc(
            PetscObjectTypeCompare(reinterpret_cast<PetscObject>(preconditioner), PCFIELDSPLIT, &already_configured),
            "PetscObjectTypeCompare field split");
        PetscInt global_count = 0;
        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        check_petsc(VecGetSize(objects.state, &global_count), "VecGetSize field split");
        if (global_count % 3 != 0)
            throw std::invalid_argument("field_split requires the [T, ur, uz] three-field layout");
        check_petsc(VecGetOwnershipRange(objects.state, &ownership_begin, &ownership_end),
                    "VecGetOwnershipRange field split");
        const PetscInt node_count = global_count / 3;
        const PetscInt temperature_begin = std::max<PetscInt>(ownership_begin, 0);
        const PetscInt temperature_end = std::min<PetscInt>(ownership_end, node_count);
        const PetscInt mechanics_begin = std::max<PetscInt>(ownership_begin, node_count);
        const PetscInt mechanics_end = std::min<PetscInt>(ownership_end, global_count);
        const PetscInt temperature_count = std::max<PetscInt>(0, temperature_end - temperature_begin);
        const PetscInt mechanics_count = std::max<PetscInt>(0, mechanics_end - mechanics_begin);
        check_petsc(PCSetType(preconditioner, PCFIELDSPLIT), "PCSetType FIELDSPLIT");
        if (already_configured == PETSC_FALSE) {
            IS temperature = nullptr;
            IS mechanics = nullptr;
            try {
                check_petsc(ISCreateStride(PETSC_COMM_WORLD, temperature_count, temperature_begin, 1, &temperature),
                            "ISCreateStride temperature");
                check_petsc(ISCreateStride(PETSC_COMM_WORLD, mechanics_count, mechanics_begin, 1, &mechanics),
                            "ISCreateStride mechanics");
                check_petsc(PCFieldSplitSetIS(preconditioner, "temperature", temperature),
                            "PCFieldSplitSetIS temperature");
                check_petsc(PCFieldSplitSetIS(preconditioner, "mechanics", mechanics), "PCFieldSplitSetIS mechanics");
                check_petsc(ISDestroy(&temperature), "ISDestroy temperature");
                check_petsc(ISDestroy(&mechanics), "ISDestroy mechanics");
            } catch (...) {
                if (temperature != nullptr)
                    (void)ISDestroy(&temperature);
                if (mechanics != nullptr)
                    (void)ISDestroy(&mechanics);
                throw;
            }
        }
        check_petsc(PCFieldSplitSetType(preconditioner, PC_COMPOSITE_MULTIPLICATIVE),
                    "PCFieldSplitSetType multiplicative");
        break;
    }
    case SolverOptions::Preconditioner::hypre:
        check_petsc(PCSetType(preconditioner, PCHYPRE), "PCSetType HYPRE");
        break;
    }
    check_petsc(KSPSetTolerances(ksp, options.linear_relative_tolerance, PETSC_DEFAULT, PETSC_DEFAULT,
                                 options.maximum_linear_iterations),
                "KSPSetTolerances");
}

PetscErrorCode gather_state(Vec state, SolverContext& context) {
    PetscFunctionBeginUser;
    PetscCall(VecScatterBegin(context.state_scatter, state, context.gathered_state, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(context.state_scatter, state, context.gathered_state, INSERT_VALUES, SCATTER_FORWARD));
    const PetscScalar* values = nullptr;
    PetscCall(VecGetArrayRead(context.gathered_state, &values));
    for (std::size_t index = 0; index < context.state_values.size(); ++index)
        context.state_values[index] = PetscRealPart(values[checked_petsc_int(index)]);
    PetscCall(VecRestoreArrayRead(context.gathered_state, &values));
    PetscFunctionReturn(PETSC_SUCCESS);
}

std::vector<double> gather_complete_state(Vec state, std::size_t global_size) {
    VecScatter scatter = nullptr;
    Vec gathered = nullptr;
    check_petsc(VecScatterCreateToAll(state, &scatter, &gathered), "VecScatterCreateToAll final state");
    try {
        check_petsc(VecScatterBegin(scatter, state, gathered, INSERT_VALUES, SCATTER_FORWARD),
                    "VecScatterBegin final state");
        check_petsc(VecScatterEnd(scatter, state, gathered, INSERT_VALUES, SCATTER_FORWARD),
                    "VecScatterEnd final state");
        std::vector<double> result(global_size, 0.0);
        const PetscScalar* values = nullptr;
        check_petsc(VecGetArrayRead(gathered, &values), "VecGetArrayRead final state");
        for (std::size_t index = 0; index < global_size; ++index)
            result[index] = PetscRealPart(values[checked_petsc_int(index)]);
        check_petsc(VecRestoreArrayRead(gathered, &values), "VecRestoreArrayRead final state");
        check_petsc(VecScatterDestroy(&scatter), "VecScatterDestroy final state");
        check_petsc(VecDestroy(&gathered), "VecDestroy final state");
        return result;
    } catch (...) {
        if (scatter != nullptr)
            (void)VecScatterDestroy(&scatter);
        if (gathered != nullptr)
            (void)VecDestroy(&gathered);
        throw;
    }
}

PetscErrorCode synchronize_domain_error(bool local_error, bool& global_error) {
    PetscFunctionBeginUser;
    const PetscMPIInt local = local_error ? 1 : 0;
    PetscMPIInt global = 0;
    PetscCallMPI(MPIU_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD));
    global_error = global != 0;
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode field_norms(Vec vector, std::array<double, 3>& norms) {
    PetscFunctionBeginUser;
    PetscInt global_count = 0;
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    PetscCall(VecGetSize(vector, &global_count));
    if (global_count % 3 != 0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "fuelsim field residual scaling requires [T, ur, uz]");
    PetscCall(VecGetOwnershipRange(vector, &ownership_begin, &ownership_end));
    const PetscScalar* values = nullptr;
    PetscCall(VecGetArrayRead(vector, &values));
    std::array<double, 3> local_squared{};
    const PetscInt node_count = global_count / 3;
    for (PetscInt global = ownership_begin; global < ownership_end; ++global) {
        const std::size_t field = static_cast<std::size_t>(global / node_count);
        const double value = PetscRealPart(values[global - ownership_begin]);
        local_squared[field] += value * value;
    }
    PetscCall(VecRestoreArrayRead(vector, &values));
    std::array<double, 3> global_squared{};
    if (PetscGlobalSize == 1) {
        global_squared = local_squared;
    } else {
        PetscCallMPI(
            MPIU_Allreduce(local_squared.data(), global_squared.data(), 3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD));
    }
    for (std::size_t field = 0; field < norms.size(); ++field)
        norms[field] = std::sqrt(global_squared[field]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode scale_residual(Vec residual, SolverContext& context) {
    PetscFunctionBeginUser;
    PetscCall(field_norms(residual, context.latest_unscaled_field_residual_norms));
    if (context.first_residual) {
        context.initial_field_residual_norms = context.latest_unscaled_field_residual_norms;
        context.first_residual = false;
    }
    if (!context.field_residual_scaling) {
        context.latest_field_residual_norms = context.latest_unscaled_field_residual_norms;
        for (std::size_t field = 0; field < context.field_residual_reference_norms.size(); ++field)
            context.field_residual_reference_norms[field] =
                std::max(context.field_residual_reference_norms[field], context.latest_field_residual_norms[field]);
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    const double thermal_norm = context.latest_unscaled_field_residual_norms[0];
    if (!context.thermal_scaling_initialized && thermal_norm > context.residual_scaling_floor) {
        context.field_residual_scalings[0] = 1.0 / thermal_norm;
        context.initial_field_residual_norms[0] = thermal_norm;
        context.thermal_scaling_initialized = true;
    }
    const double mechanics_norm =
        std::hypot(context.latest_unscaled_field_residual_norms[1], context.latest_unscaled_field_residual_norms[2]);
    if (!context.mechanics_scaling_initialized && mechanics_norm > context.residual_scaling_floor) {
        const double mechanics_scaling = 1.0 / mechanics_norm;
        context.field_residual_scalings[1] = mechanics_scaling;
        context.field_residual_scalings[2] = mechanics_scaling;
        context.initial_field_residual_norms[1] = context.latest_unscaled_field_residual_norms[1];
        context.initial_field_residual_norms[2] = context.latest_unscaled_field_residual_norms[2];
        context.mechanics_scaling_initialized = true;
    }

    PetscInt global_count = 0;
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    PetscCall(VecGetSize(residual, &global_count));
    PetscCall(VecGetOwnershipRange(residual, &ownership_begin, &ownership_end));
    const PetscInt node_count = global_count / 3;
    PetscScalar* values = nullptr;
    PetscCall(VecGetArray(residual, &values));
    for (PetscInt global = ownership_begin; global < ownership_end; ++global) {
        const std::size_t index = static_cast<std::size_t>(global);
        if (!context.constrained[index]) {
            const std::size_t field = static_cast<std::size_t>(global / node_count);
            values[global - ownership_begin] *= context.field_residual_scalings[field];
        }
    }
    PetscCall(VecRestoreArray(residual, &values));
    PetscCall(field_norms(residual, context.latest_field_residual_norms));
    for (std::size_t field = 0; field < context.field_residual_reference_norms.size(); ++field)
        context.field_residual_reference_norms[field] =
            std::max(context.field_residual_reference_norms[field], context.latest_field_residual_norms[field]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_function(SNES snes, Vec state, Vec residual, void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        context.last_function_domain_error = false;
        const NonlinearProblem& problem = *context.problem;
        PetscCall(gather_state(state, context));
        const GlobalStateView state_view(problem.dof_count(), context.shadow_dofs, context.state_values);
        PetscCall(VecSet(residual, 0.0));
        bool local_domain_error = false;
        try {
            problem.validate_local_state(context.contribution_begin, context.contribution_end, state_view);
            for (std::size_t contribution = context.contribution_begin; contribution < context.contribution_end;
                 ++contribution) {
                const LocalDofs size_dofs = problem.contribution_dofs(contribution);
                const LocalValues local_state = problem.contribution_state(contribution, state_view);
                const LocalResidual local_residual = problem.contribution_residual(contribution, local_state);
                std::array<PetscInt, local_dof_count> dofs{};
                for (std::size_t local = 0; local < dofs.size(); ++local)
                    dofs[local] = checked_petsc_int(size_dofs[local]);
                PetscCall(VecSetValues(residual, static_cast<PetscInt>(dofs.size()), dofs.data(), local_residual.data(),
                                       ADD_VALUES));
            }
        } catch (const std::domain_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        } catch (const std::overflow_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        }
        bool global_domain_error = false;
        PetscCall(synchronize_domain_error(local_domain_error, global_domain_error));
        if (global_domain_error && context.last_domain_error.empty())
            context.last_domain_error = "a residual evaluation violated its physical domain on "
                                        "another MPI rank";
        PetscCall(VecAssemblyBegin(residual));
        PetscCall(VecAssemblyEnd(residual));
        if (global_domain_error) {
            context.saw_domain_error = true;
            context.last_function_domain_error = true;
            PetscCall(VecSet(residual, 0.0));
            PetscCall(SNESSetFunctionDomainError(snes));
            context.timing.residual_callback_seconds += seconds_since(start);
            ++context.timing.residual_evaluations;
            PetscFunctionReturn(PETSC_SUCCESS);
        }

        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        PetscCall(VecGetOwnershipRange(residual, &ownership_begin, &ownership_end));
        PetscScalar* local_residual = nullptr;
        PetscCall(VecGetArray(residual, &local_residual));
        for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
            const PetscInt dof = checked_petsc_int(condition.dof);
            if (dof >= ownership_begin && dof < ownership_end)
                local_residual[dof - ownership_begin] = state_view.value(condition.dof) - condition.value;
        }
        PetscCall(VecRestoreArray(residual, &local_residual));
        PetscCall(scale_residual(residual, context));
        if (!std::isfinite(context.initial_residual_norm)) {
            PetscReal norm = 0.0;
            PetscCall(VecNorm(residual, NORM_2, &norm));
            context.initial_residual_norm = static_cast<double>(norm);
        }
        context.timing.residual_callback_seconds += seconds_since(start);
        ++context.timing.residual_evaluations;
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "fuelsim residual assembly failed: %s", error.what());
    } catch (...) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "fuelsim residual assembly failed with unknown error");
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_jacobian(SNES snes, Vec state, Mat jacobian, Mat preconditioner, void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        if (jacobian != preconditioner)
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "fuelsim requires one matrix for J and P");

        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        const NonlinearProblem& problem = *context.problem;
        PetscCall(gather_state(state, context));
        const GlobalStateView state_view(problem.dof_count(), context.shadow_dofs, context.state_values);

        PetscCall(MatZeroEntries(jacobian));
        bool local_domain_error = false;
        try {
            problem.validate_local_state(context.contribution_begin, context.contribution_end, state_view);
            for (std::size_t contribution = context.contribution_begin; contribution < context.contribution_end;
                 ++contribution) {
                const LocalValues local_state = problem.contribution_state(contribution, state_view);
                const LocalSystem local_system = problem.linearize_contribution(contribution, local_state);
                const LocalDofs size_dofs = problem.contribution_dofs(contribution);
                std::array<PetscInt, local_dof_count> dofs{};
                for (std::size_t local = 0; local < dofs.size(); ++local)
                    dofs[local] = checked_petsc_int(size_dofs[local]);

                std::array<double, local_dof_count * local_dof_count> scaled_jacobian = local_system.jacobian;
                if (context.field_residual_scaling) {
                    const std::size_t node_count = problem.dof_count() / 3;
                    for (std::size_t row = 0; row < dofs.size(); ++row) {
                        const std::size_t global_row = size_dofs[row];
                        if (context.constrained[global_row])
                            continue;
                        const std::size_t field = global_row / node_count;
                        for (std::size_t column = 0; column < dofs.size(); ++column)
                            scaled_jacobian[row * dofs.size() + column] *= context.field_residual_scalings[field];
                    }
                }
                PetscCall(MatSetValues(jacobian, static_cast<PetscInt>(dofs.size()), dofs.data(),
                                       static_cast<PetscInt>(dofs.size()), dofs.data(), scaled_jacobian.data(),
                                       ADD_VALUES));
            }
        } catch (const std::domain_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        } catch (const std::overflow_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        }
        bool global_domain_error = false;
        PetscCall(synchronize_domain_error(local_domain_error, global_domain_error));
        if (global_domain_error && context.last_domain_error.empty())
            context.last_domain_error = "a Jacobian evaluation violated its physical domain on "
                                        "another MPI rank";
        PetscCall(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));
        if (global_domain_error) {
            context.saw_domain_error = true;
            PetscCall(MatZeroEntries(jacobian));
            PetscCall(SNESSetJacobianDomainError(snes));
            context.timing.jacobian_callback_seconds += seconds_since(start);
            ++context.timing.jacobian_evaluations;
            PetscFunctionReturn(PETSC_SUCCESS);
        }

        // The first Jacobian assembly inserts zero-valued blocks for every
        // potential contact candidate. Keep that complete pattern locked so
        // dynamic search cannot introduce a new nonzero location later.
        if (!context.pattern_locked) {
            PetscCall(MatSetOption(jacobian, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
            PetscCall(MatSetOption(jacobian, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));
            context.pattern_locked = true;
        }

        PetscCall(MatZeroRows(jacobian, static_cast<PetscInt>(context.constrained_dofs.size()),
                              context.constrained_dofs.data(), 1.0, nullptr, nullptr));
        context.timing.jacobian_callback_seconds += seconds_since(start);
        ++context.timing.jacobian_evaluations;
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "fuelsim Jacobian assembly failed: %s", error.what());
    } catch (...) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "fuelsim Jacobian assembly failed with unknown error");
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
        _context.contribution_end = problem.contribution_count() * (rank + 1U) / size;
        _context.constrained.assign(problem.dof_count(), false);

        check_petsc(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, _count, &_objects->state), "VecCreateMPI state");
        check_petsc(VecDuplicate(_objects->state, &_objects->residual), "VecDuplicate residual");

        PetscInt local_count = 0;
        check_petsc(VecGetLocalSize(_objects->state, &local_count), "VecGetLocalSize state");
        const PetscInt diagonal_nonzeros = std::min<PetscInt>(60, local_count);
        const PetscInt off_diagonal_nonzeros = std::min<PetscInt>(60, _count - local_count);
        check_petsc(MatCreateAIJ(PETSC_COMM_WORLD, local_count, local_count, _count, _count, diagonal_nonzeros, nullptr,
                                 off_diagonal_nonzeros, nullptr, &_objects->jacobian),
                    "MatCreateAIJ");
        check_petsc(MatSetOption(_objects->jacobian, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
                    "MatSetOption MAT_NEW_NONZERO_ALLOCATION_ERR");
        check_petsc(SNESCreate(PETSC_COMM_WORLD, &_objects->snes), "SNESCreate");
        check_petsc(SNESSetFunction(_objects->snes, _objects->residual, form_function, &_context), "SNESSetFunction");
        check_petsc(SNESSetJacobian(_objects->snes, _objects->jacobian, _objects->jacobian, form_jacobian, &_context),
                    "SNESSetJacobian");
        check_petsc(SNESSetType(_objects->snes, SNESNEWTONLS), "SNESSetType");

        SNESLineSearch line_search = nullptr;
        check_petsc(SNESGetLineSearch(_objects->snes, &line_search), "SNESGetLineSearch");
        check_petsc(SNESLineSearchSetType(line_search, SNESLINESEARCHBASIC), "SNESLineSearchSetType");

        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        check_petsc(VecGetOwnershipRange(_objects->state, &ownership_begin, &ownership_end),
                    "VecGetOwnershipRange state");
        _context.ownership_begin = ownership_begin;
        _context.ownership_end = ownership_end;
        for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
            const PetscInt dof = checked_petsc_int(condition.dof);
            if (dof >= ownership_begin && dof < ownership_end)
                _context.constrained_dofs.push_back(dof);
            _context.constrained[condition.dof] = true;
        }
        const std::vector<std::size_t> required_dofs =
            problem.required_state_dofs(_context.contribution_begin, _context.contribution_end);
        _context.shadow_dofs.reserve(required_dofs.size() + _context.constrained_dofs.size());
        for (const std::size_t dof : required_dofs) {
            if (dof > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
                throw std::length_error("fuelsim shadow DOF index exceeds uint32_t range");
            _context.shadow_dofs.push_back(static_cast<std::uint32_t>(dof));
        }
        for (const PetscInt dof : _context.constrained_dofs)
            _context.shadow_dofs.push_back(static_cast<std::uint32_t>(dof));
        std::sort(_context.shadow_dofs.begin(), _context.shadow_dofs.end());
        _context.shadow_dofs.erase(std::unique(_context.shadow_dofs.begin(), _context.shadow_dofs.end()),
                                   _context.shadow_dofs.end());
        _context.state_values.resize(_context.shadow_dofs.size());
        std::vector<PetscInt> shadow_indices;
        shadow_indices.reserve(_context.shadow_dofs.size());
        for (const std::uint32_t dof : _context.shadow_dofs)
            shadow_indices.push_back(checked_petsc_int(static_cast<std::size_t>(dof)));
        const PetscInt shadow_count = checked_petsc_int(_context.shadow_dofs.size());
        IS source_indices = nullptr;
        IS destination_indices = nullptr;
        try {
            check_petsc(ISCreateGeneral(PETSC_COMM_SELF, shadow_count, shadow_indices.data(), PETSC_COPY_VALUES,
                                        &source_indices),
                        "ISCreateGeneral shadow state");
            check_petsc(ISCreateStride(PETSC_COMM_SELF, shadow_count, 0, 1, &destination_indices),
                        "ISCreateStride shadow state");
            check_petsc(VecCreateSeq(PETSC_COMM_SELF, shadow_count, &_objects->gathered_state),
                        "VecCreateSeq shadow state");
            check_petsc(VecScatterCreate(_objects->state, source_indices, _objects->gathered_state, destination_indices,
                                         &_objects->state_scatter),
                        "VecScatterCreate shadow state");
            check_petsc(ISDestroy(&source_indices), "ISDestroy shadow source");
            check_petsc(ISDestroy(&destination_indices), "ISDestroy shadow destination");
        } catch (...) {
            if (source_indices != nullptr)
                (void)ISDestroy(&source_indices);
            if (destination_indices != nullptr)
                (void)ISDestroy(&destination_indices);
            throw;
        }
        _context.state_scatter = _objects->state_scatter;
        _context.gathered_state = _objects->gathered_state;
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
        check_petsc(PetscInitialize(&argc, &argv, nullptr, help), "PetscInitialize");
        _owns_initialization = true;
    }
    _rank = static_cast<int>(PetscGlobalRank);
    _size = static_cast<int>(PetscGlobalSize);
}

PetscSession::~PetscSession() {
    if (!_owns_initialization)
        return;

    PetscBool finalized = PETSC_FALSE;
    if (PetscFinalized(&finalized) == PETSC_SUCCESS && finalized == PETSC_FALSE) {
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

void PetscSession::collective_root_action(const std::function<void()>& action) const {
    if (!action)
        throw std::invalid_argument("PetscSession collective root action must not be empty");
    bool failed = false;
    std::string message;
    if (_rank == 0) {
        try {
            action();
        } catch (const std::exception& error) {
            failed = true;
            message = error.what();
        } catch (...) {
            failed = true;
            message = "unknown root-rank I/O error";
        }
    }
    Vec status = nullptr;
    check_petsc(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 1, &status), "VecCreateMPI(root I/O status)");
    if (_rank == 0 && failed)
        check_petsc(VecSetValue(status, 0, 1.0, INSERT_VALUES), "VecSetValue(root I/O status)");
    check_petsc(VecAssemblyBegin(status), "VecAssemblyBegin(root I/O status)");
    check_petsc(VecAssemblyEnd(status), "VecAssemblyEnd(root I/O status)");
    PetscReal failure_norm = 0.0;
    check_petsc(VecNorm(status, NORM_1, &failure_norm), "VecNorm(root I/O status)");
    check_petsc(VecDestroy(&status), "VecDestroy(root I/O status)");
    if (failure_norm == 0.0)
        return;
    throw std::runtime_error(_rank == 0 ? "collective root-rank I/O failed: " + message
                                        : "collective root-rank I/O failed; see rank 0 for details");
}

PetscSolver::PetscSolver() : _implementation(std::make_unique<Implementation>()) {}

PetscSolver::~PetscSolver() = default;

SolveResult PetscSolver::solve(const NonlinearProblem& problem, const std::vector<double>& initial_state,
                               const SolverOptions& options) {
    SolveResult result = solve_once(problem, initial_state, options);
    if (result.converged || !options.backtracking_fallback || options.line_search != SolverOptions::LineSearch::basic)
        return result;

    SolverOptions fallback_options = options;
    fallback_options.line_search = SolverOptions::LineSearch::backtracking;
    fallback_options.backtracking_fallback = false;
    SolveResult fallback = solve_once(problem, initial_state, fallback_options);
    fallback.nonlinear_iterations += result.nonlinear_iterations;
    fallback.linear_iterations += result.linear_iterations;
    accumulate_timing(fallback.timing, result.timing);
    fallback.nonlinear_attempts = result.nonlinear_attempts + 1;
    fallback.used_backtracking_fallback = true;
    fallback.basic_failure_category = result.failure_category;
    fallback.basic_failure_message = result.failure_message;
    return fallback;
}

SolveResult PetscSolver::solve_once(const NonlinearProblem& problem, const std::vector<double>& initial_state,
                                    const SolverOptions& options) {
    if (initial_state.size() != problem.dof_count())
        throw std::invalid_argument("PetscSolver initial state size mismatch");
    const bool fixed_temperature_scale = options.temperature_residual_scale > 0.0;
    const bool fixed_mechanical_scale = options.mechanical_residual_scale > 0.0;
    const bool residual_scaling = options.field_residual_scaling || fixed_temperature_scale;

    const SteadyClock::time_point total_start = SteadyClock::now();
    const SteadyClock::time_point setup_start = SteadyClock::now();
    const bool workspace_created = _implementation->prepare(problem);
    PetscObjects& objects = _implementation->objects();
    SolverContext& context = _implementation->context();
    context.timing = SolveTiming{};
    context.timing.workspace_setups = workspace_created ? 1U : 0U;
    context.timing.solve_calls = 1;
    context.initial_residual_norm = std::numeric_limits<double>::quiet_NaN();
    context.field_residual_scaling = residual_scaling;
    context.residual_scaling_floor = options.absolute_tolerance;
    context.first_residual = true;
    context.thermal_scaling_initialized = fixed_temperature_scale;
    context.mechanics_scaling_initialized = fixed_mechanical_scale;
    context.initial_field_residual_norms = {};
    context.field_residual_reference_norms = {};
    context.latest_unscaled_field_residual_norms = {};
    context.latest_field_residual_norms = {};
    context.field_residual_scalings = {1.0, 1.0, 1.0};
    if (fixed_temperature_scale) {
        context.field_residual_scalings[0] = 1.0 / options.temperature_residual_scale;
        context.field_residual_scalings[1] = 1.0 / options.mechanical_residual_scale;
        context.field_residual_scalings[2] = 1.0 / options.mechanical_residual_scale;
    }
    context.saw_domain_error = false;
    context.last_function_domain_error = false;
    context.last_domain_error.clear();
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    check_petsc(VecGetOwnershipRange(objects.state, &ownership_begin, &ownership_end), "VecGetOwnershipRange state");
    PetscScalar* state_array = nullptr;
    check_petsc(VecGetArray(objects.state, &state_array), "VecGetArray state");
    for (PetscInt index = ownership_begin; index < ownership_end; ++index)
        state_array[index - ownership_begin] = initial_state[static_cast<std::size_t>(index)];
    check_petsc(VecRestoreArray(objects.state, &state_array), "VecRestoreArray state");

    check_petsc(SNESSetTolerances(objects.snes, options.absolute_tolerance, options.relative_tolerance,
                                  residual_scaling ? 0.0 : options.step_tolerance, options.maximum_iterations,
                                  PETSC_DEFAULT),
                "SNESSetTolerances");
    SNESLineSearch line_search = nullptr;
    check_petsc(SNESGetLineSearch(objects.snes, &line_search), "SNESGetLineSearch");
    check_petsc(SNESLineSearchSetType(line_search, options.line_search == SolverOptions::LineSearch::backtracking
                                                       ? SNESLINESEARCHBT
                                                       : SNESLINESEARCHBASIC),
                "SNESLineSearchSetType");
    configure_linear_solver(objects, options, PetscGlobalSize);
    check_petsc(SNESSetFromOptions(objects.snes), "SNESSetFromOptions");
    context.timing.setup_seconds = seconds_since(setup_start);

    const SteadyClock::time_point solve_start = SteadyClock::now();
    check_petsc(SNESSolve(objects.snes, nullptr, objects.state), "SNESSolve");
    context.timing.nonlinear_solve_seconds = seconds_since(solve_start);

    SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
    PetscInt iterations = 0;
    PetscInt linear_iterations = 0;
    PetscReal residual_norm = 0.0;
    check_petsc(SNESGetConvergedReason(objects.snes, &reason), "SNESGetConvergedReason");
    check_petsc(SNESGetIterationNumber(objects.snes, &iterations), "SNESGetIterationNumber");
    check_petsc(SNESGetLinearSolveIterations(objects.snes, &linear_iterations), "SNESGetLinearSolveIterations");
    check_petsc(SNESGetFunctionNorm(objects.snes, &residual_norm), "SNESGetFunctionNorm");
    const bool requires_explicit_residual_audit = reason == SNES_CONVERGED_SNORM_RELATIVE || reason < 0;
    if (requires_explicit_residual_audit) {
        context.last_function_domain_error = false;
        check_petsc(SNESComputeFunction(objects.snes, objects.state, objects.residual),
                    "SNESComputeFunction final residual");
        check_petsc(VecNorm(objects.residual, NORM_2, &residual_norm), "VecNorm final residual");
    }
    const bool final_domain_error = context.last_function_domain_error;

    std::vector<double> solution = gather_complete_state(objects.state, problem.dof_count());

    context.timing.total_seconds = seconds_since(total_start);
    SolveResult result;
    result.state = std::move(solution);
    result.nonlinear_iterations = static_cast<int>(iterations);
    result.linear_iterations = static_cast<int>(linear_iterations);
    result.residual_norm = static_cast<double>(residual_norm);
    result.convergence_reason = static_cast<int>(reason);
    result.initial_field_residual_norms = context.initial_field_residual_norms;
    result.field_residual_reference_norms = context.field_residual_reference_norms;
    result.final_field_residual_norms = context.latest_unscaled_field_residual_norms;
    result.final_scaled_field_residual_norms = context.latest_field_residual_norms;
    result.field_residual_scalings = context.field_residual_scalings;
    const double configured_residual_threshold =
        std::max(options.absolute_tolerance, options.relative_tolerance * context.initial_residual_norm);
    const double fallback_reduction = options.residual_reduction_tolerance;
    const double numerical_residual_floor = 10.0 * std::sqrt(std::numeric_limits<double>::epsilon());
    const double independently_verified_threshold =
        std::max(numerical_residual_floor, fallback_reduction * context.initial_residual_norm);
    const double physical_absolute_threshold =
        std::hypot(options.temperature_residual_absolute_tolerance * result.field_residual_scalings[0],
                   std::hypot(options.mechanical_residual_absolute_tolerance * result.field_residual_scalings[1],
                              options.mechanical_residual_absolute_tolerance * result.field_residual_scalings[2]));
    const double residual_threshold =
        std::max({configured_residual_threshold, independently_verified_threshold, physical_absolute_threshold});
    const double residual_slack = residual_threshold * (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
    bool fields_verified = true;
    std::array<double, 3> field_thresholds{};
    for (std::size_t field = 0; field < result.final_scaled_field_residual_norms.size(); ++field) {
        const double initial_scaled =
            result.initial_field_residual_norms[field] * result.field_residual_scalings[field];
        const double reference = std::max(initial_scaled, result.field_residual_reference_norms[field]);
        const double configured_field_threshold = std::max(
            field == 0 ? options.temperature_residual_absolute_tolerance * result.field_residual_scalings[field]
                       : options.mechanical_residual_absolute_tolerance * result.field_residual_scalings[field],
            options.relative_tolerance * reference);
        const double independent_field_threshold = std::max(numerical_residual_floor, fallback_reduction * reference);
        const double field_threshold = std::max(configured_field_threshold, independent_field_threshold) *
                                       (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
        field_thresholds[field] = field_threshold;
        fields_verified = std::isfinite(result.final_scaled_field_residual_norms[field]) &&
                          result.final_scaled_field_residual_norms[field] <= field_threshold && fields_verified;
    }
    const bool residual_verified = std::isfinite(result.residual_norm) &&
                                   std::isfinite(context.initial_residual_norm) &&
                                   result.residual_norm <= residual_slack && fields_verified;
    const bool recoverable_stopping_reason =
        reason == SNES_DIVERGED_LINE_SEARCH || reason == SNES_DIVERGED_MAX_IT || reason == SNES_DIVERGED_LOCAL_MIN;
    result.converged = residual_verified && !final_domain_error && (reason > 0 || recoverable_stopping_reason);
    if (!result.converged && reason < 0) {
        result.failure_category = context.saw_domain_error ? SolveFailureCategory::physical_domain
                                                           : SolveFailureCategory::nonlinear_divergence;
        if (!context.last_domain_error.empty()) {
            result.failure_message = context.last_domain_error;
        } else if (context.saw_domain_error) {
            result.failure_message = "a residual or Jacobian evaluation violated its physical domain";
        } else {
            result.failure_message = petsc_convergence_reason_name(result.convergence_reason);
        }
    } else if (!result.converged && !residual_verified) {
        result.failure_category = SolveFailureCategory::residual_verification;
        std::ostringstream message;
        message << "PETSc reported convergence without satisfying residual "
                   "tolerances: global="
                << result.residual_norm << '/' << residual_slack;
        for (std::size_t field = 0; field < field_thresholds.size(); ++field)
            message << ", field" << field << '=' << result.final_scaled_field_residual_norms[field] << '/'
                    << field_thresholds[field] << " (initial=" << result.initial_field_residual_norms[field]
                    << ", reference=" << result.field_residual_reference_norms[field] << ')';
        result.failure_message = message.str();
    }
    check_petsc(collective_timing(context.timing, result.timing), "collective timing reduction");
    result.mpi_rank = static_cast<int>(context.rank);
    result.mpi_size = static_cast<int>(context.size);
    result.local_contribution_begin = context.contribution_begin;
    result.local_contribution_end = context.contribution_end;
    result.global_state_dofs = problem.dof_count();
    const PetscInt64 local_shadow = static_cast<PetscInt64>(context.shadow_dofs.size());
    PetscInt64 maximum_shadow = 0;
    PetscInt64 total_shadow = 0;
    PetscInt64 local_remote_shadow = 0;
    for (const std::uint32_t dof : context.shadow_dofs) {
        const PetscInt petsc_dof = checked_petsc_int(static_cast<std::size_t>(dof));
        if (petsc_dof < context.ownership_begin || petsc_dof >= context.ownership_end)
            ++local_remote_shadow;
    }
    PetscInt64 total_remote_shadow = 0;
    check_mpi(MPIU_Allreduce(&local_shadow, &maximum_shadow, 1, MPIU_INT64, MPI_MAX, PETSC_COMM_WORLD),
              "MPIU_Allreduce maximum shadow state");
    check_mpi(MPIU_Allreduce(&local_shadow, &total_shadow, 1, MPIU_INT64, MPI_SUM, PETSC_COMM_WORLD),
              "MPIU_Allreduce total shadow state");
    check_mpi(MPIU_Allreduce(&local_remote_shadow, &total_remote_shadow, 1, MPIU_INT64, MPI_SUM, PETSC_COMM_WORLD),
              "MPIU_Allreduce remote shadow state");
    result.maximum_shadow_state_dofs = static_cast<std::size_t>(maximum_shadow);
    result.total_shadow_state_dofs = static_cast<std::size_t>(total_shadow);
    result.total_remote_shadow_state_dofs = static_cast<std::size_t>(total_remote_shadow);
    return result;
}

std::string petsc_convergence_reason_name(int reason) {
    const char* name = SNESConvergedReasons[reason];
    if (name == nullptr)
        return "UNKNOWN";
    return name;
}

const char* solve_failure_category_name(SolveFailureCategory category) noexcept {
    switch (category) {
    case SolveFailureCategory::none:
        return "none";
    case SolveFailureCategory::nonlinear_divergence:
        return "nonlinear_divergence";
    case SolveFailureCategory::physical_domain:
        return "physical_domain";
    case SolveFailureCategory::residual_verification:
        return "residual_verification";
    case SolveFailureCategory::time_discretization:
        return "time_discretization";
    case SolveFailureCategory::contact_constraint:
        return "contact_constraint";
    }
    return "unknown";
}

// Shared steady and transient solve helpers.
namespace solver_workflow {

namespace {

void merge_attempt(SolveResult& aggregate, const SolveResult& addition) {
    const int nonlinear_iterations = aggregate.nonlinear_iterations + addition.nonlinear_iterations;
    const int linear_iterations = aggregate.linear_iterations + addition.linear_iterations;
    const std::size_t nonlinear_attempts = aggregate.nonlinear_attempts + addition.nonlinear_attempts;
    const std::size_t augmented_iterations =
        aggregate.augmented_lagrangian_iterations + addition.augmented_lagrangian_iterations;
    const double maximum_penetration =
        std::max(aggregate.maximum_contact_penetration, addition.maximum_contact_penetration);
    SolveTiming timing = aggregate.timing;
    accumulate_timing(timing, addition.timing);
    const bool used_backtracking = aggregate.used_backtracking_fallback || addition.used_backtracking_fallback;
    const SolveFailureCategory basic_failure = aggregate.basic_failure_category != SolveFailureCategory::none
                                                   ? aggregate.basic_failure_category
                                                   : addition.basic_failure_category;
    const std::string basic_message =
        !aggregate.basic_failure_message.empty() ? aggregate.basic_failure_message : addition.basic_failure_message;
    aggregate = addition;
    aggregate.nonlinear_iterations = nonlinear_iterations;
    aggregate.linear_iterations = linear_iterations;
    aggregate.nonlinear_attempts = nonlinear_attempts;
    aggregate.timing = timing;
    aggregate.used_backtracking_fallback = used_backtracking;
    aggregate.basic_failure_category = basic_failure;
    aggregate.basic_failure_message = basic_message;
    aggregate.augmented_lagrangian_iterations = augmented_iterations;
    aggregate.maximum_contact_penetration = maximum_penetration;
}

void mark_augmented_failure(SolveResult& result, const AugmentedContactUpdate& status, std::size_t completed_updates) {
    result.converged = false;
    result.failure_category = SolveFailureCategory::contact_constraint;
    result.failure_message =
        "Augmented contact did not reach penetration tolerance after " + std::to_string(completed_updates) +
        " multiplier updates; maximum constraint violation=" + std::to_string(status.maximum_constraint_violation) +
        ", maximum penetration=" + std::to_string(status.maximum_penetration) +
        ", tolerance=" + std::to_string(status.penetration_tolerance);
}

} // namespace

std::vector<double> initial_guess_with_dirichlet_values(const NonlinearProblem& problem,
                                                        const std::vector<double>& state) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument("Dirichlet initial-guess state size does not match problem");
    std::vector<double> result = state;
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = condition.value;
    return result;
}

SolveResult solve_contact_equilibrium(PetscSolver& solver, SteadyProblem& problem,
                                      const std::vector<double>& initial_guess, const SolverOptions& options) {
    SolveResult result = solver.solve(problem, initial_guess, options);
    if (!problem.uses_augmented_contact())
        return result;
    std::size_t updates = 0;
    while (result.converged) {
        const AugmentedContactUpdate status = problem.update_augmented_contact_multipliers(result.state, updates);
        result.maximum_contact_penetration = status.maximum_penetration;
        result.augmented_lagrangian_iterations = updates;
        if (status.converged)
            return result;
        if (!status.update_allowed) {
            mark_augmented_failure(result, status, updates);
            return result;
        }
        ++updates;
        SolveResult next = solver.solve(problem, initial_guess_with_dirichlet_values(problem, result.state), options);
        merge_attempt(result, next);
    }
    result.augmented_lagrangian_iterations = updates;
    return result;
}

SolveResult solve_contact_equilibrium(PetscSolver& solver, TransientProblem& problem,
                                      const std::vector<double>& initial_guess, const SolverOptions& options) {
    SolveResult result = solver.solve(problem, initial_guess, options);
    if (!problem.uses_augmented_contact())
        return result;
    std::size_t updates = 0;
    while (result.converged) {
        const AugmentedContactUpdate status = problem.update_augmented_contact_multipliers(result.state, updates);
        result.maximum_contact_penetration = status.maximum_penetration;
        result.augmented_lagrangian_iterations = updates;
        if (status.converged)
            return result;
        if (!status.update_allowed) {
            mark_augmented_failure(result, status, updates);
            return result;
        }
        ++updates;
        SolveResult next = solver.solve(problem, initial_guess_with_dirichlet_values(problem, result.state), options);
        merge_attempt(result, next);
    }
    result.augmented_lagrangian_iterations = updates;
    return result;
}

} // namespace solver_workflow

// Steady load stepping.

using solver_workflow::initial_guess_with_dirichlet_values;
using solver_workflow::solve_contact_equilibrium;

SteadyResult solve_steady(SteadyProblem& problem, const SteadyLoadOptions& load_options, const SolverOptions& options) {
    const SteadyClock::time_point start = SteadyClock::now();
    SteadyResult result;
    PetscSolver solver;
    std::vector<double> state = problem.initial_state();
    double accepted_load_factor = 0.0;
    for (std::size_t step = 1; step <= load_options.load_steps; ++step) {
        const double target_load_factor = static_cast<double>(step) / static_cast<double>(load_options.load_steps);
        while (accepted_load_factor < target_load_factor) {
            std::size_t cutbacks = 0;
            double attempted_load_factor = target_load_factor;
            double load_increment = attempted_load_factor - accepted_load_factor;
            SolveResult attempt;
            for (;;) {
                const std::vector<std::vector<ContactPointHistory>> committed_contact_histories =
                    problem.committed_contact_histories();
                attempt = SolveResult{};
                try {
                    problem.set_load_factor(attempted_load_factor);
                    attempt = solve_contact_equilibrium(solver, problem,
                                                        initial_guess_with_dirichlet_values(problem, state), options);
                    if (attempt.converged)
                        problem.commit_contact_state(attempt.state);
                } catch (const std::domain_error& error) {
                    attempt.converged = false;
                    attempt.failure_category = SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                } catch (const std::overflow_error& error) {
                    attempt.converged = false;
                    attempt.failure_category = SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                }
                if (!attempt.converged) {
                    problem.restore_contact_state(state, committed_contact_histories);
                    problem.set_load_factor(accepted_load_factor);
                }
                accumulate_timing(result.aggregate_timing, attempt.timing);
                result.total_nonlinear_iterations += attempt.nonlinear_iterations;
                result.total_linear_iterations += attempt.linear_iterations;
                if (attempt.converged)
                    break;

                result.rejected_steps.push_back({attempted_load_factor, load_increment, cutbacks,
                                                 attempt.failure_category, attempt.failure_message});
                result.solve = attempt;
                if (cutbacks >= load_options.maximum_cutbacks_per_step) {
                    result.total_seconds = seconds_since(start);
                    return result;
                }
                const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
                if (load_increment <= load_options.minimum_load_increment + tolerance) {
                    result.total_seconds = seconds_since(start);
                    return result;
                }
                load_increment =
                    std::max(load_increment * load_options.cutback_factor, load_options.minimum_load_increment);
                attempted_load_factor = accepted_load_factor + load_increment;
                ++cutbacks;
                ++result.total_cutbacks;
            }
            accepted_load_factor = attempted_load_factor;
            state = attempt.state;
            result.solve = std::move(attempt);
        }
        result.completed_steps = step;
    }
    result.completed = true;
    result.total_seconds = seconds_since(start);
    return result;
}

// Backward-Euler step-doubling error control.
namespace time_control {
namespace {

struct ErrorAccumulator final {
    double difference_squared = 0.0;
    double solution_squared = 0.0;
    std::size_t count = 0;
};

void accumulate_error(ErrorAccumulator& accumulator, double full_step, double two_half_steps) {
    const double difference = two_half_steps - full_step;
    accumulator.difference_squared += difference * difference;
    accumulator.solution_squared += two_half_steps * two_half_steps;
    ++accumulator.count;
}

double normalized_error(const ErrorAccumulator& accumulator, double absolute_tolerance, double relative_tolerance) {
    if (accumulator.count == 0)
        return 0.0;
    const double denominator = absolute_tolerance * std::sqrt(static_cast<double>(accumulator.count)) +
                               relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}

} // namespace

TransientConservationSummary combine_half_step_conservation(const TransientConservationSummary& first,
                                                            const TransientConservationSummary& second) {
    TransientConservationSummary result;
    const auto average = [](double left, double right) { return 0.5 * (left + right); };
    result.generated_heat_rate = average(first.generated_heat_rate, second.generated_heat_rate);
    result.stored_heat_rate = average(first.stored_heat_rate, second.stored_heat_rate);
    result.convection_heat_rate = average(first.convection_heat_rate, second.convection_heat_rate);
    result.interface_heat_imbalance = average(first.interface_heat_imbalance, second.interface_heat_imbalance);
    result.dirichlet_heat_input_rate = average(first.dirichlet_heat_input_rate, second.dirichlet_heat_input_rate);
    result.global_thermal_balance = result.stored_heat_rate + result.convection_heat_rate +
                                    result.interface_heat_imbalance - result.generated_heat_rate -
                                    result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.generated_heat_rate) + std::abs(result.stored_heat_rate) +
                                 std::abs(result.convection_heat_rate) + std::abs(result.interface_heat_imbalance) +
                                 std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0 ? std::abs(result.global_thermal_balance) / thermal_scale : 0.0;
    result.unconstrained_thermal_residual_l2 =
        std::max(first.unconstrained_thermal_residual_l2, second.unconstrained_thermal_residual_l2);

    result.internal_mechanical_work_increment =
        first.internal_mechanical_work_increment + second.internal_mechanical_work_increment;
    result.pressure_traction_work_increment =
        first.pressure_traction_work_increment + second.pressure_traction_work_increment;
    result.dirichlet_reaction_work_increment =
        first.dirichlet_reaction_work_increment + second.dirichlet_reaction_work_increment;
    result.contact_work_increment = first.contact_work_increment + second.contact_work_increment;
    result.mechanical_work_balance = result.internal_mechanical_work_increment + result.contact_work_increment -
                                     result.pressure_traction_work_increment - result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) + std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment) + std::abs(result.contact_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0 ? std::abs(result.mechanical_work_balance) / mechanical_scale : 0.0;
    result.unconstrained_mechanical_residual_l2 =
        std::max(first.unconstrained_mechanical_residual_l2, second.unconstrained_mechanical_residual_l2);
    result.elastic_energy_change = first.elastic_energy_change + second.elastic_energy_change;
    result.plastic_dissipation_increment = first.plastic_dissipation_increment + second.plastic_dissipation_increment;
    result.creep_dissipation_increment = first.creep_dissipation_increment + second.creep_dissipation_increment;
    return result;
}

TransientTimeErrorEstimate step_doubling_error(const TransientCommittedState& full_step,
                                               const TransientCommittedState& two_half_steps,
                                               const TransientTimeOptions& options) {
    if (full_step.solution.size() != two_half_steps.solution.size() || full_step.solution.size() % 3 != 0)
        throw std::logic_error("step-doubling states must share the [T, ur, uz] layout");
    if (full_step.material_histories.size() != two_half_steps.material_histories.size() ||
        full_step.material_stresses.size() != two_half_steps.material_stresses.size())
        throw std::logic_error("step-doubling material-state region layouts differ");

    const std::size_t node_count = full_step.solution.size() / 3;
    std::array<ErrorAccumulator, 3> nodal{};
    for (std::size_t field = 0; field < 3; ++field) {
        const std::size_t begin = field * node_count;
        const std::size_t end = begin + node_count;
        for (std::size_t dof = begin; dof < end; ++dof)
            accumulate_error(nodal[field], full_step.solution[dof], two_half_steps.solution[dof]);
    }

    ErrorAccumulator elastic;
    ErrorAccumulator plastic;
    ErrorAccumulator creep;
    ErrorAccumulator equivalent_plastic;
    ErrorAccumulator equivalent_creep;
    ErrorAccumulator stress;
    ErrorAccumulator contact_friction;
    ErrorAccumulator contact_normal_multiplier;
    bool contact_state_mismatch = false;
    for (std::size_t region = 0; region < full_step.material_histories.size(); ++region) {
        const auto& full_history = full_step.material_histories[region];
        const auto& half_history = two_half_steps.material_histories[region];
        const auto& full_stress = full_step.material_stresses[region];
        const auto& half_stress = two_half_steps.material_stresses[region];
        if (full_history.size() != half_history.size() || full_stress.size() != half_stress.size() ||
            full_history.size() != full_stress.size())
            throw std::logic_error("step-doubling material-state element layouts differ");
        for (std::size_t element = 0; element < full_history.size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const MaterialPointState& full_point = full_history[element][q];
                const MaterialPointState& half_point = half_history[element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_error(elastic, full_point.elastic_strain[component],
                                     half_point.elastic_strain[component]);
                    accumulate_error(plastic, full_point.plastic_strain[component],
                                     half_point.plastic_strain[component]);
                    accumulate_error(creep, full_point.creep_strain[component], half_point.creep_strain[component]);
                }
                accumulate_error(equivalent_plastic, full_point.equivalent_plastic_strain,
                                 half_point.equivalent_plastic_strain);
                accumulate_error(equivalent_creep, full_point.equivalent_creep_strain,
                                 half_point.equivalent_creep_strain);

                const AxisymmetricStressValues& full_value = full_stress[element][q];
                const AxisymmetricStressValues& half_value = half_stress[element][q];
                accumulate_error(stress, full_value.rr, half_value.rr);
                accumulate_error(stress, full_value.zz, half_value.zz);
                accumulate_error(stress, full_value.hoop, half_value.hoop);
                accumulate_error(stress, full_value.rz, half_value.rz);
            }
        }
    }
    if (full_step.contact_histories.size() != two_half_steps.contact_histories.size())
        throw std::logic_error("step-doubling contact-history layouts differ");
    for (std::size_t contact = 0; contact < full_step.contact_histories.size(); ++contact) {
        if (full_step.contact_histories[contact].size() != two_half_steps.contact_histories[contact].size())
            throw std::logic_error("step-doubling contact-node history layouts differ");
        for (std::size_t node = 0; node < full_step.contact_histories[contact].size(); ++node) {
            const ContactPointHistory& full = full_step.contact_histories[contact][node];
            const ContactPointHistory& half = two_half_steps.contact_histories[contact][node];
            accumulate_error(contact_friction, full.elastic_tangential_slip, half.elastic_tangential_slip);
            accumulate_error(contact_normal_multiplier, full.normal_multiplier, half.normal_multiplier);
            contact_state_mismatch = contact_state_mismatch || full.sliding != half.sliding;
        }
    }

    TransientTimeErrorEstimate result;
    result.temperature =
        normalized_error(nodal[0], options.temperature_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.radial_displacement =
        normalized_error(nodal[1], options.displacement_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.axial_displacement =
        normalized_error(nodal[2], options.displacement_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.elastic_strain = normalized_error(elastic, options.strain_history_time_absolute_tolerance,
                                             options.time_error_relative_tolerance);
    result.plastic_strain = normalized_error(plastic, options.strain_history_time_absolute_tolerance,
                                             options.time_error_relative_tolerance);
    result.creep_strain =
        normalized_error(creep, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.equivalent_plastic_strain = normalized_error(
        equivalent_plastic, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.equivalent_creep_strain = normalized_error(equivalent_creep, options.strain_history_time_absolute_tolerance,
                                                      options.time_error_relative_tolerance);
    result.stress =
        normalized_error(stress, options.stress_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.contact_friction = contact_state_mismatch
                                  ? std::numeric_limits<double>::infinity()
                                  : normalized_error(contact_friction, options.displacement_time_absolute_tolerance,
                                                     options.time_error_relative_tolerance);
    result.contact_normal_multiplier =
        normalized_error(contact_normal_multiplier, options.stress_history_time_absolute_tolerance,
                         options.time_error_relative_tolerance);
    result.maximum = std::max({result.temperature, result.radial_displacement, result.axial_displacement,
                               result.elastic_strain, result.plastic_strain, result.creep_strain,
                               result.equivalent_plastic_strain, result.equivalent_creep_strain, result.stress,
                               result.contact_friction, result.contact_normal_multiplier});
    return result;
}

double step_factor(const TransientTimeOptions& options, double error) {
    if (!(error > 0.0))
        return options.growth_factor;
    return std::clamp(options.time_error_safety_factor / std::sqrt(error), 0.1, options.growth_factor);
}

} // namespace time_control

// Transient time stepping and state transactions.
namespace {

using solver_workflow::merge_attempt;

class TimeStepTransaction final {
  public:
    TimeStepTransaction(TransientProblem& problem, const TransientStepInput& input)
        : _problem(problem), _committed(false) {
        _problem.begin_time_step(input);
    }

    ~TimeStepTransaction() {
        if (!_committed)
            _problem.rollback_time_step();
    }

    TimeStepTransaction(const TimeStepTransaction&) = delete;
    TimeStepTransaction& operator=(const TimeStepTransaction&) = delete;

    void commit(const std::vector<double>& solution) {
        _problem.commit_time_step(solution);
        _committed = true;
    }

  private:
    TransientProblem& _problem;
    bool _committed;
};

bool reaches_end(double time, double end_time) {
    const double scale = std::max({1.0, std::abs(time), std::abs(end_time)});
    return end_time - time <= 16.0 * std::numeric_limits<double>::epsilon() * scale;
}

void validate_time_options(const TransientProblem& problem, const TransientTimeOptions& options) {
    if (problem.time_step_active())
        throw std::logic_error("solve_transient cannot start with an active time step");
    const double time_scale = std::max({1.0, std::abs(problem.committed_time()), std::abs(options.end_time)});
    const double time_tolerance = 16.0 * std::numeric_limits<double>::epsilon() * time_scale;
    if (!std::isfinite(options.end_time) || options.end_time < problem.committed_time() - time_tolerance)
        throw std::invalid_argument("solve_transient end time must not precede committed time");
}

double load_factor_at_time(const TransientTimeOptions& options, double time) {
    if (options.load_ramp_time == 0.0)
        return 1.0;
    return std::min(time / options.load_ramp_time, 1.0);
}

double accepted_next_time_step(const TransientTimeOptions& options, double actual_time_step,
                               double controller_time_step, bool event_truncated, std::size_t cutbacks,
                               int nonlinear_iterations) {
    const double base = event_truncated && cutbacks == 0 ? controller_time_step : actual_time_step;
    if (options.target_nonlinear_iterations == 0) {
        if (cutbacks > 0)
            return std::clamp(base, options.minimum_time_step, options.maximum_time_step);
        return std::min(options.maximum_time_step, base * options.growth_factor);
    }
    const std::size_t iterations = nonlinear_iterations < 0 ? 0 : static_cast<std::size_t>(nonlinear_iterations);
    const std::size_t lower = options.target_nonlinear_iterations - options.iteration_window;
    const std::size_t upper =
        options.target_nonlinear_iterations > std::numeric_limits<std::size_t>::max() - options.iteration_window
            ? std::numeric_limits<std::size_t>::max()
            : options.target_nonlinear_iterations + options.iteration_window;
    if (iterations < lower && cutbacks == 0)
        return std::min(options.maximum_time_step, base * options.growth_factor);
    if (iterations > upper)
        return std::max(options.minimum_time_step, base * options.cutback_factor);
    return std::clamp(base, options.minimum_time_step, options.maximum_time_step);
}

} // namespace

TransientResult solve_transient(TransientProblem& problem, const TransientTimeOptions& options,
                                const SolverOptions& solver_options, TransientStepObserver* observer) {
    validate_time_options(problem, options);
    const SteadyClock::time_point start = SteadyClock::now();
    TransientResult result;
    PetscSolver solver;
    const std::vector<double> events = problem.time_events();
    double next_time_step = options.initial_time_step;
    const auto run_step = [&](double target_time) {
        TimeStepTransaction transaction(problem, {target_time, load_factor_at_time(options, target_time)});
        SolveResult step_result = solve_contact_equilibrium(
            solver, problem, initial_guess_with_dirichlet_values(problem, problem.committed_solution()),
            solver_options);
        if (step_result.converged)
            transaction.commit(step_result.state);
        return step_result;
    };
    while (!reaches_end(problem.committed_time(), options.end_time)) {
        const double controller_time_step = std::min(next_time_step, options.end_time - problem.committed_time());
        double time_step = controller_time_step;
        bool event_truncated = false;
        for (const double event : events) {
            if (reaches_end(problem.committed_time(), event))
                continue;
            if (event >= options.end_time || reaches_end(event, options.end_time))
                break;
            const double event_step = event - problem.committed_time();
            if (event_step < time_step && !reaches_end(event, problem.committed_time() + time_step)) {
                time_step = event_step;
                event_truncated = true;
            }
            break;
        }
        std::size_t cutbacks = 0;
        for (;;) {
            const double end_time = problem.committed_time() + time_step;
            SolveResult attempt;
            double time_error_estimate = 0.0;
            TransientTimeErrorEstimate time_error_components;
            TransientConservationSummary first_half_conservation;
            int controller_nonlinear_iterations = 0;
            const bool error_control = options.time_error_relative_tolerance > 0.0;
            TransientCommittedState base_state;
            if (error_control)
                base_state = problem.committed_state();
            try {
                if (!error_control) {
                    attempt = run_step(end_time);
                    controller_nonlinear_iterations = attempt.nonlinear_iterations;
                } else {
                    const SolveResult full_step = run_step(end_time);
                    TransientCommittedState full_step_state;
                    if (full_step.converged)
                        full_step_state = problem.committed_state();
                    attempt = full_step;
                    controller_nonlinear_iterations = full_step.nonlinear_iterations;
                    if (full_step.converged) {
                        problem.restore_committed_state(base_state);
                        const double half_time = base_state.time + 0.5 * time_step;
                        const SolveResult first_half = run_step(half_time);
                        if (first_half.converged)
                            first_half_conservation = problem.last_conservation_summary();
                        controller_nonlinear_iterations =
                            std::max(controller_nonlinear_iterations, first_half.nonlinear_iterations);
                        merge_attempt(attempt, first_half);
                        if (!first_half.converged) {
                            problem.restore_committed_state(std::move(base_state));
                        } else {
                            const SolveResult second_half = run_step(end_time);
                            controller_nonlinear_iterations =
                                std::max(controller_nonlinear_iterations, second_half.nonlinear_iterations);
                            merge_attempt(attempt, second_half);
                            if (!second_half.converged) {
                                problem.restore_committed_state(std::move(base_state));
                            } else {
                                time_error_components = time_control::step_doubling_error(
                                    full_step_state, problem.committed_state(), options);
                                time_error_estimate = time_error_components.maximum;
                                if (!(time_error_estimate <= 1.0)) {
                                    problem.restore_committed_state(std::move(base_state));
                                    attempt.converged = false;
                                    attempt.failure_category = SolveFailureCategory::time_discretization;
                                    attempt.failure_message = "Backward-Euler step-doubling error "
                                                              "exceeded one";
                                    ++result.time_error_rejections;
                                } else {
                                    TransientCommittedState accepted_state = problem.committed_state();
                                    accepted_state.conservation = time_control::combine_half_step_conservation(
                                        first_half_conservation, accepted_state.conservation);
                                    problem.restore_committed_state(std::move(accepted_state));
                                }
                            }
                        }
                    }
                }
            } catch (const std::domain_error& error) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                attempt.converged = false;
                attempt.failure_category = SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (const std::overflow_error& error) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                attempt.converged = false;
                attempt.failure_category = SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (...) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                throw;
            }
            accumulate_timing(result.aggregate_timing, attempt.timing);
            result.total_nonlinear_iterations += attempt.nonlinear_iterations;
            result.total_linear_iterations += attempt.linear_iterations;
            result.last_attempt = std::move(attempt);
            if (result.last_attempt.converged) {
                next_time_step = accepted_next_time_step(options, time_step, controller_time_step, event_truncated,
                                                         cutbacks, controller_nonlinear_iterations);
                if (error_control) {
                    const double error_limited_step =
                        time_step * time_control::step_factor(options, time_error_estimate);
                    next_time_step = std::clamp(std::min(next_time_step, error_limited_step), options.minimum_time_step,
                                                options.maximum_time_step);
                }
                result.accepted_steps.push_back(
                    {problem.committed_time(), time_step, next_time_step, problem.committed_load_factor(), cutbacks,
                     result.last_attempt.nonlinear_iterations, result.last_attempt.linear_iterations,
                     time_error_estimate, time_error_components, problem.last_conservation_summary()});
                if (observer != nullptr)
                    observer->accepted_step(problem, result.accepted_steps.back());
                break;
            }
            result.rejected_steps.push_back(
                {problem.committed_time() + time_step, time_step, cutbacks, result.last_attempt.nonlinear_iterations,
                 result.last_attempt.linear_iterations, result.last_attempt.convergence_reason,
                 result.last_attempt.residual_norm, result.last_attempt.failure_category,
                 result.last_attempt.failure_message, time_error_estimate, time_error_components});
            if (cutbacks >= options.maximum_cutbacks_per_step) {
                result.termination_reason = TransientTerminationReason::maximum_cutbacks;
                break;
            }
            double reduction_factor = options.cutback_factor;
            if (result.last_attempt.failure_category == SolveFailureCategory::time_discretization)
                reduction_factor = std::min(
                    0.9, std::max(options.cutback_factor, time_control::step_factor(options, time_error_estimate)));
            const double reduced = time_step * reduction_factor;
            const double minimum_tolerance =
                16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, options.minimum_time_step);
            if (time_step <= options.minimum_time_step + minimum_tolerance) {
                result.termination_reason = TransientTerminationReason::minimum_time_step;
                break;
            }
            time_step = std::max(reduced, options.minimum_time_step);
            ++cutbacks;
            ++result.total_cutbacks;
        }
        if (!result.last_attempt.converged)
            break;
    }
    result.completed = reaches_end(problem.committed_time(), options.end_time);
    if (result.completed)
        result.termination_reason = TransientTerminationReason::completed;
    result.committed_state = problem.committed_solution();
    result.committed_time = problem.committed_time();
    result.next_time_step = next_time_step;
    result.total_seconds = seconds_since(start);
    return result;
}

const char* transient_termination_reason_name(TransientTerminationReason reason) noexcept {
    switch (reason) {
    case TransientTerminationReason::not_started:
        return "not_started";
    case TransientTerminationReason::completed:
        return "completed";
    case TransientTerminationReason::maximum_cutbacks:
        return "maximum_cutbacks";
    case TransientTerminationReason::minimum_time_step:
        return "minimum_time_step";
    }
    return "unknown";
}

} // namespace fuelsim
