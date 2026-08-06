#include "fuelsim/petsc_solver.hpp"

#include <petsc.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <sstream>
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
    double initial_residual_norm =
        std::numeric_limits<double>::quiet_NaN();
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

PetscErrorCode synchronize_domain_error(bool local_error,
                                        bool& global_error) {
    PetscFunctionBeginUser;
    const PetscMPIInt local = local_error ? 1 : 0;
    PetscMPIInt global = 0;
    PetscCallMPI(MPIU_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX,
                                PETSC_COMM_WORLD));
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
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                "fuelsim field residual scaling requires [T, ur, uz]");
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
        PetscCallMPI(MPIU_Allreduce(local_squared.data(),
                                   global_squared.data(), 3, MPI_DOUBLE,
                                   MPI_SUM, PETSC_COMM_WORLD));
    }
    for (std::size_t field = 0; field < norms.size(); ++field)
        norms[field] = std::sqrt(global_squared[field]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode scale_residual(Vec residual, SolverContext& context) {
    PetscFunctionBeginUser;
    PetscCall(
        field_norms(residual, context.latest_unscaled_field_residual_norms));
    if (context.first_residual) {
        context.initial_field_residual_norms =
            context.latest_unscaled_field_residual_norms;
        context.first_residual = false;
    }
    if (!context.field_residual_scaling) {
        context.latest_field_residual_norms =
            context.latest_unscaled_field_residual_norms;
        for (std::size_t field = 0;
             field < context.field_residual_reference_norms.size(); ++field)
            context.field_residual_reference_norms[field] = std::max(
                context.field_residual_reference_norms[field],
                context.latest_field_residual_norms[field]);
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    const double thermal_norm =
        context.latest_unscaled_field_residual_norms[0];
    if (!context.thermal_scaling_initialized &&
        thermal_norm > context.residual_scaling_floor) {
        context.field_residual_scalings[0] = 1.0 / thermal_norm;
        context.initial_field_residual_norms[0] = thermal_norm;
        context.thermal_scaling_initialized = true;
    }
    const double mechanics_norm = std::hypot(
        context.latest_unscaled_field_residual_norms[1],
        context.latest_unscaled_field_residual_norms[2]);
    if (!context.mechanics_scaling_initialized &&
        mechanics_norm > context.residual_scaling_floor) {
        const double mechanics_scaling = 1.0 / mechanics_norm;
        context.field_residual_scalings[1] = mechanics_scaling;
        context.field_residual_scalings[2] = mechanics_scaling;
        context.initial_field_residual_norms[1] =
            context.latest_unscaled_field_residual_norms[1];
        context.initial_field_residual_norms[2] =
            context.latest_unscaled_field_residual_norms[2];
        context.mechanics_scaling_initialized = true;
    }

    PetscInt global_count = 0;
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    PetscCall(VecGetSize(residual, &global_count));
    PetscCall(VecGetOwnershipRange(residual, &ownership_begin,
                                   &ownership_end));
    const PetscInt node_count = global_count / 3;
    PetscScalar* values = nullptr;
    PetscCall(VecGetArray(residual, &values));
    for (PetscInt global = ownership_begin; global < ownership_end; ++global) {
        const std::size_t index = static_cast<std::size_t>(global);
        if (!context.constrained[index]) {
            const std::size_t field =
                static_cast<std::size_t>(global / node_count);
            values[global - ownership_begin] *=
                context.field_residual_scalings[field];
        }
    }
    PetscCall(VecRestoreArray(residual, &values));
    PetscCall(field_norms(residual, context.latest_field_residual_norms));
    for (std::size_t field = 0;
         field < context.field_residual_reference_norms.size(); ++field)
        context.field_residual_reference_norms[field] = std::max(
            context.field_residual_reference_norms[field],
            context.latest_field_residual_norms[field]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_function(SNES snes, Vec state, Vec residual,
                             void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        context.last_function_domain_error = false;
        const NonlinearProblem& problem = *context.problem;
        PetscCall(gather_state(state, context));
        PetscCall(VecSet(residual, 0.0));
        bool local_domain_error = false;
        try {
            problem.validate_state(context.state_values);
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
        } catch (const std::domain_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        } catch (const std::overflow_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        }
        bool global_domain_error = false;
        PetscCall(synchronize_domain_error(local_domain_error,
                                           global_domain_error));
        if (global_domain_error && context.last_domain_error.empty())
            context.last_domain_error =
                "a residual evaluation violated its physical domain on "
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
        PetscCall(scale_residual(residual, context));
        if (!std::isfinite(context.initial_residual_norm)) {
            PetscReal norm = 0.0;
            PetscCall(VecNorm(residual, NORM_2, &norm));
            context.initial_residual_norm = static_cast<double>(norm);
        }
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

PetscErrorCode form_jacobian(SNES snes, Vec state, Mat jacobian,
                             Mat preconditioner,
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
        bool local_domain_error = false;
        try {
            problem.validate_state(context.state_values);
            for (std::size_t contribution = context.contribution_begin;
                 contribution < context.contribution_end; ++contribution) {
                const LocalValues local_state = problem.contribution_state(
                    contribution, context.state_values);
                const LocalSystem local_system =
                    problem.linearize_contribution(contribution, local_state);
                const LocalDofs size_dofs =
                    problem.contribution_dofs(contribution);
                std::array<PetscInt, local_dof_count> dofs{};
                for (std::size_t local = 0; local < dofs.size(); ++local)
                    dofs[local] = checked_petsc_int(size_dofs[local]);

                std::array<double, local_dof_count * local_dof_count>
                    scaled_jacobian = local_system.jacobian;
                if (context.field_residual_scaling) {
                    const std::size_t node_count = problem.dof_count() / 3;
                    for (std::size_t row = 0; row < dofs.size(); ++row) {
                        const std::size_t global_row = size_dofs[row];
                        if (context.constrained[global_row])
                            continue;
                        const std::size_t field = global_row / node_count;
                        for (std::size_t column = 0; column < dofs.size();
                             ++column)
                            scaled_jacobian[row * dofs.size() + column] *=
                                context.field_residual_scalings[field];
                    }
                }
                PetscCall(MatSetValues(
                    jacobian, static_cast<PetscInt>(dofs.size()), dofs.data(),
                    static_cast<PetscInt>(dofs.size()), dofs.data(),
                    scaled_jacobian.data(), ADD_VALUES));
            }
        } catch (const std::domain_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        } catch (const std::overflow_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        }
        bool global_domain_error = false;
        PetscCall(synchronize_domain_error(local_domain_error,
                                           global_domain_error));
        if (global_domain_error && context.last_domain_error.empty())
            context.last_domain_error =
                "a Jacobian evaluation violated its physical domain on "
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

        if (!context.pattern_locked) {
            PetscCall(
                MatSetOption(jacobian, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
            PetscCall(MatSetOption(jacobian, MAT_NEW_NONZERO_LOCATION_ERR,
                                   PETSC_TRUE));
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
        _context.constrained.assign(problem.dof_count(), false);

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
            _context.constrained[condition.dof] = true;
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

void PetscSession::collective_root_action(
    const std::function<void()>& action) const {
    if (!action)
        throw std::invalid_argument(
            "PetscSession collective root action must not be empty");
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
    check_petsc(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 1, &status),
                "VecCreateMPI(root I/O status)");
    if (_rank == 0 && failed)
        check_petsc(VecSetValue(status, 0, 1.0, INSERT_VALUES),
                    "VecSetValue(root I/O status)");
    check_petsc(VecAssemblyBegin(status),
                "VecAssemblyBegin(root I/O status)");
    check_petsc(VecAssemblyEnd(status), "VecAssemblyEnd(root I/O status)");
    PetscReal failure_norm = 0.0;
    check_petsc(VecNorm(status, NORM_1, &failure_norm),
                "VecNorm(root I/O status)");
    check_petsc(VecDestroy(&status), "VecDestroy(root I/O status)");
    if (failure_norm == 0.0)
        return;
    throw std::runtime_error(
        _rank == 0
            ? "collective root-rank I/O failed: " + message
            : "collective root-rank I/O failed; see rank 0 for details");
}

PetscSolver::PetscSolver()
    : _implementation(std::make_unique<Implementation>()) {}

PetscSolver::~PetscSolver() = default;

SolveResult
PetscSolver::solve(const NonlinearProblem& problem,
                   const std::vector<double>& initial_state,
                   const SolverOptions& options) {
    SolveResult result = solve_once(problem, initial_state, options);
    if (result.converged || !options.backtracking_fallback ||
        options.line_search != SolverOptions::LineSearch::basic)
        return result;

    SolverOptions fallback_options = options;
    fallback_options.line_search = SolverOptions::LineSearch::backtracking;
    fallback_options.backtracking_fallback = false;
    SolveResult fallback = solve_once(problem, initial_state, fallback_options);
    fallback.nonlinear_iterations += result.nonlinear_iterations;
    fallback.linear_iterations += result.linear_iterations;
    fallback.timing.setup_seconds += result.timing.setup_seconds;
    fallback.timing.nonlinear_solve_seconds +=
        result.timing.nonlinear_solve_seconds;
    fallback.timing.residual_callback_seconds +=
        result.timing.residual_callback_seconds;
    fallback.timing.jacobian_callback_seconds +=
        result.timing.jacobian_callback_seconds;
    fallback.timing.total_seconds += result.timing.total_seconds;
    fallback.timing.residual_evaluations += result.timing.residual_evaluations;
    fallback.timing.jacobian_evaluations += result.timing.jacobian_evaluations;
    fallback.timing.workspace_setups += result.timing.workspace_setups;
    fallback.timing.solve_calls += result.timing.solve_calls;
    fallback.nonlinear_attempts = result.nonlinear_attempts + 1;
    fallback.used_backtracking_fallback = true;
    fallback.basic_failure_category = result.failure_category;
    fallback.basic_failure_message = result.failure_message;
    return fallback;
}

SolveResult
PetscSolver::solve_once(const NonlinearProblem& problem,
                        const std::vector<double>& initial_state,
                        const SolverOptions& options) {
    if (initial_state.size() != problem.dof_count())
        throw std::invalid_argument(
            "PetscSolver initial state size mismatch");
    if (!(options.absolute_tolerance > 0.0) ||
        !(options.relative_tolerance > 0.0) ||
        !(options.step_tolerance > 0.0) || options.maximum_iterations <= 0 ||
        !(options.linear_relative_tolerance > 0.0) ||
        !(options.residual_reduction_tolerance > 0.0) ||
        !(options.temperature_residual_absolute_tolerance > 0.0) ||
        !(options.mechanical_residual_absolute_tolerance > 0.0) ||
        !std::isfinite(options.temperature_residual_scale) ||
        !std::isfinite(options.mechanical_residual_scale) ||
        options.temperature_residual_scale < 0.0 ||
        options.mechanical_residual_scale < 0.0 ||
        options.maximum_linear_iterations <= 0)
        throw std::invalid_argument(
            "PetscSolver tolerances and iterations must be positive");
    const bool fixed_temperature_scale =
        options.temperature_residual_scale > 0.0;
    const bool fixed_mechanical_scale =
        options.mechanical_residual_scale > 0.0;
    if (fixed_temperature_scale != fixed_mechanical_scale)
        throw std::invalid_argument(
            "PetscSolver fixed temperature and mechanical residual scales "
            "must both be zero or positive");
    if (fixed_temperature_scale && options.field_residual_scaling)
        throw std::invalid_argument(
            "PetscSolver fixed residual scales cannot be combined with "
            "automatic field residual scaling");
    const bool residual_scaling =
        options.field_residual_scaling || fixed_temperature_scale;

    const SteadyClock::time_point total_start = SteadyClock::now();
    const SteadyClock::time_point setup_start = SteadyClock::now();
    const bool workspace_created = _implementation->prepare(problem);
    PetscObjects& objects = _implementation->objects();
    SolverContext& context = _implementation->context();
    context.timing = SolveTiming{};
    context.timing.workspace_setups = workspace_created ? 1U : 0U;
    context.timing.solve_calls = 1;
    context.initial_residual_norm =
        std::numeric_limits<double>::quiet_NaN();
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
        context.field_residual_scalings[0] =
            1.0 / options.temperature_residual_scale;
        context.field_residual_scalings[1] =
            1.0 / options.mechanical_residual_scale;
        context.field_residual_scalings[2] =
            1.0 / options.mechanical_residual_scale;
    }
    context.saw_domain_error = false;
    context.last_function_domain_error = false;
    context.last_domain_error.clear();
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
                                  residual_scaling ? 0.0
                                                   : options.step_tolerance,
                                  options.maximum_iterations, PETSC_DEFAULT),
                "SNESSetTolerances");
    SNESLineSearch line_search = nullptr;
    check_petsc(SNESGetLineSearch(objects.snes, &line_search),
                "SNESGetLineSearch");
    check_petsc(
        SNESLineSearchSetType(
            line_search,
            options.line_search == SolverOptions::LineSearch::backtracking
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
    check_petsc(SNESGetConvergedReason(objects.snes, &reason),
                "SNESGetConvergedReason");
    check_petsc(SNESGetIterationNumber(objects.snes, &iterations),
                "SNESGetIterationNumber");
    check_petsc(SNESGetLinearSolveIterations(objects.snes,
                                             &linear_iterations),
                "SNESGetLinearSolveIterations");
    check_petsc(SNESGetFunctionNorm(objects.snes, &residual_norm),
                "SNESGetFunctionNorm");
    const bool requires_explicit_residual_audit =
        reason == SNES_CONVERGED_SNORM_RELATIVE || reason < 0;
    if (requires_explicit_residual_audit) {
        context.last_function_domain_error = false;
        check_petsc(SNESComputeFunction(objects.snes, objects.state,
                                        objects.residual),
                    "SNESComputeFunction final residual");
        check_petsc(VecNorm(objects.residual, NORM_2, &residual_norm),
                    "VecNorm final residual");
    }
    const bool final_domain_error = context.last_function_domain_error;

    check_petsc(gather_state(objects.state, context), "gather_state solution");
    std::vector<double> solution(problem.dof_count(), 0.0);
    solution = context.state_values;

    context.timing.total_seconds = seconds_since(total_start);
    SolveResult result;
    result.state = std::move(solution);
    result.nonlinear_iterations = static_cast<int>(iterations);
    result.linear_iterations = static_cast<int>(linear_iterations);
    result.residual_norm = static_cast<double>(residual_norm);
    result.convergence_reason = static_cast<int>(reason);
    result.initial_field_residual_norms =
        context.initial_field_residual_norms;
    result.field_residual_reference_norms =
        context.field_residual_reference_norms;
    result.final_field_residual_norms =
        context.latest_unscaled_field_residual_norms;
    result.final_scaled_field_residual_norms =
        context.latest_field_residual_norms;
    result.field_residual_scalings = context.field_residual_scalings;
    const double configured_residual_threshold = std::max(
        options.absolute_tolerance,
        options.relative_tolerance * context.initial_residual_norm);
    const double fallback_reduction =
        options.residual_reduction_tolerance;
    const double numerical_residual_floor =
        10.0 * std::sqrt(std::numeric_limits<double>::epsilon());
    const double independently_verified_threshold =
        std::max(numerical_residual_floor,
                 fallback_reduction * context.initial_residual_norm);
    const double physical_absolute_threshold = std::hypot(
        options.temperature_residual_absolute_tolerance *
            result.field_residual_scalings[0],
        std::hypot(options.mechanical_residual_absolute_tolerance *
                       result.field_residual_scalings[1],
                   options.mechanical_residual_absolute_tolerance *
                       result.field_residual_scalings[2]));
    const double residual_threshold = std::max(
        {configured_residual_threshold, independently_verified_threshold,
         physical_absolute_threshold});
    const double residual_slack =
        residual_threshold *
        (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
    bool fields_verified = true;
    std::array<double, 3> field_thresholds{};
    for (std::size_t field = 0;
         field < result.final_scaled_field_residual_norms.size(); ++field) {
        const double initial_scaled =
            result.initial_field_residual_norms[field] *
            result.field_residual_scalings[field];
        const double reference = std::max(
            initial_scaled, result.field_residual_reference_norms[field]);
        const double configured_field_threshold =
            std::max(field == 0
                         ? options.temperature_residual_absolute_tolerance *
                               result.field_residual_scalings[field]
                         : options.mechanical_residual_absolute_tolerance *
                               result.field_residual_scalings[field],
                     options.relative_tolerance * reference);
        const double independent_field_threshold =
            std::max(numerical_residual_floor,
                     fallback_reduction * reference);
        const double field_threshold =
            std::max(configured_field_threshold,
                     independent_field_threshold) *
            (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
        field_thresholds[field] = field_threshold;
        fields_verified =
            std::isfinite(result.final_scaled_field_residual_norms[field]) &&
            result.final_scaled_field_residual_norms[field] <=
                field_threshold &&
            fields_verified;
    }
    const bool residual_verified =
        std::isfinite(result.residual_norm) &&
        std::isfinite(context.initial_residual_norm) &&
        result.residual_norm <= residual_slack && fields_verified;
    const bool recoverable_stopping_reason =
        reason == SNES_DIVERGED_LINE_SEARCH ||
        reason == SNES_DIVERGED_MAX_IT ||
        reason == SNES_DIVERGED_LOCAL_MIN;
    result.converged = residual_verified && !final_domain_error &&
                       (reason > 0 || recoverable_stopping_reason);
    if (!result.converged && reason < 0) {
        result.failure_category = context.saw_domain_error
                                      ? SolveFailureCategory::physical_domain
                                      : SolveFailureCategory::nonlinear_divergence;
        if (!context.last_domain_error.empty()) {
            result.failure_message = context.last_domain_error;
        } else if (context.saw_domain_error) {
            result.failure_message =
                "a residual or Jacobian evaluation violated its physical domain";
        } else {
            result.failure_message =
                petsc_convergence_reason_name(result.convergence_reason);
        }
    } else if (!result.converged && !residual_verified) {
        result.failure_category = SolveFailureCategory::residual_verification;
        std::ostringstream message;
        message << "PETSc reported convergence without satisfying residual "
                   "tolerances: global="
                << result.residual_norm << '/' << residual_slack;
        for (std::size_t field = 0; field < field_thresholds.size(); ++field)
            message << ", field" << field << '='
                    << result.final_scaled_field_residual_norms[field] << '/'
                    << field_thresholds[field] << " (initial="
                    << result.initial_field_residual_norms[field]
                    << ", reference="
                    << result.field_residual_reference_norms[field] << ')';
        result.failure_message = message.str();
    }
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
    }
    return "unknown";
}

} // namespace fuelsim
