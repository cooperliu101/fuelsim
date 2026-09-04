#include "fuelsim/core/spatial_definition.hpp"
#include "fuelsim/solver/petsc_solver.hpp"
#include "solver_detail.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <petsc.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {
using SteadyClock = solver_detail::Clock;
using solver_detail::accumulate_timing;
using solver_detail::seconds_since;

void check_petsc(PetscErrorCode code, const char* operation) {
    if (code == PETSC_SUCCESS) return;
    const char* text = nullptr;
    const PetscErrorCode message_code = PetscErrorMessage(code, &text, nullptr);
    if (message_code != PETSC_SUCCESS) text = nullptr;
    std::string message = operation;
    message += " failed";
    if (text != nullptr) {
        message += ": ";
        message += text;
    }
    throw std::runtime_error(message);
}

void check_mpi(PetscMPIInt code, const char* operation) {
    if (code == MPI_SUCCESS) return;
    std::string message = operation;
    message += " failed with MPI error code ";
    message += std::to_string(code);
    throw std::runtime_error(message);
}

struct MemorySnapshot final {
    PetscInt64 resident_bytes = 0;
    PetscInt64 maximum_resident_bytes = 0;
};

PetscInt64 checked_memory_bytes(PetscLogDouble value) {
    if (!std::isfinite(value) || value <= 0.0) return 0;
    const PetscLogDouble maximum = static_cast<PetscLogDouble>(std::numeric_limits<PetscInt64>::max());
    if (value >= maximum) return std::numeric_limits<PetscInt64>::max();
    return static_cast<PetscInt64>(value + 0.5);
}

MemorySnapshot memory_snapshot() {
    PetscLogDouble resident = 0.0, maximum = 0.0;
    check_petsc(PetscMemoryGetCurrentUsage(&resident), "PetscMemoryGetCurrentUsage");
    check_petsc(PetscMemoryGetMaximumUsage(&maximum), "PetscMemoryGetMaximumUsage");
    return {checked_memory_bytes(resident), checked_memory_bytes(maximum)};
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
    const std::array<double, 2> local_assembly = {
        local.maximum_residual_assembly_seconds, local.maximum_jacobian_assembly_seconds};
    std::array<double, 2> minimum_assembly{}, maximum_assembly{};
    PetscCallMPI(MPIU_Allreduce(local_assembly.data(), minimum_assembly.data(),
        static_cast<MPIU_Count>(local_assembly.size()), MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPIU_Allreduce(local_assembly.data(), maximum_assembly.data(),
        static_cast<MPIU_Count>(local_assembly.size()), MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD));
    std::array<PetscInt64, 4> local_counts = {
        static_cast<PetscInt64>(local.residual_evaluations),
        static_cast<PetscInt64>(local.jacobian_evaluations),
        static_cast<PetscInt64>(local.workspace_setups),
        static_cast<PetscInt64>(local.solve_calls),
    };
    std::array<PetscInt64, 4> maximum_counts{};
    PetscCallMPI(MPIU_Allreduce(local_counts.data(), maximum_counts.data(),
        static_cast<MPIU_Count>(local_counts.size()), MPIU_INT64, MPI_MAX, PETSC_COMM_WORLD));
    const std::array<PetscInt64, 4> local_resident = {static_cast<PetscInt64>(local.initial_resident_bytes),
        static_cast<PetscInt64>(local.setup_resident_bytes), static_cast<PetscInt64>(local.solve_resident_bytes),
        static_cast<PetscInt64>(local.final_resident_bytes)};
    std::array<PetscInt64, 4> maximum_resident{};
    PetscCallMPI(MPIU_Allreduce(local_resident.data(), maximum_resident.data(),
        static_cast<MPIU_Count>(local_resident.size()), MPIU_INT64, MPI_MAX, PETSC_COMM_WORLD));
    const std::array<PetscInt64, 1> local_peak = {static_cast<PetscInt64>(local.maximum_peak_resident_bytes)};
    std::array<PetscInt64, 1> minimum_peak{}, maximum_peak{}, total_peak{};
    PetscCallMPI(MPIU_Allreduce(local_peak.data(), minimum_peak.data(), 1, MPIU_INT64, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPIU_Allreduce(local_peak.data(), maximum_peak.data(), 1, MPIU_INT64, MPI_MAX, PETSC_COMM_WORLD));
    PetscCallMPI(MPIU_Allreduce(local_peak.data(), total_peak.data(), 1, MPIU_INT64, MPI_SUM, PETSC_COMM_WORLD));
    result.setup_seconds = maximum_seconds[0];
    result.nonlinear_solve_seconds = maximum_seconds[1];
    result.residual_callback_seconds = maximum_seconds[2];
    result.jacobian_callback_seconds = maximum_seconds[3];
    result.minimum_residual_assembly_seconds = minimum_assembly[0];
    result.maximum_residual_assembly_seconds = maximum_assembly[0];
    result.minimum_jacobian_assembly_seconds = minimum_assembly[1];
    result.maximum_jacobian_assembly_seconds = maximum_assembly[1];
    result.local_residual_assembly_seconds = local_assembly[0];
    result.local_jacobian_assembly_seconds = local_assembly[1];
    result.total_seconds = maximum_seconds[4];
    result.initial_resident_bytes = static_cast<std::size_t>(maximum_resident[0]);
    result.setup_resident_bytes = static_cast<std::size_t>(maximum_resident[1]);
    result.solve_resident_bytes = static_cast<std::size_t>(maximum_resident[2]);
    result.final_resident_bytes = static_cast<std::size_t>(maximum_resident[3]);
    result.minimum_peak_resident_bytes = static_cast<std::size_t>(minimum_peak[0]);
    result.maximum_peak_resident_bytes = static_cast<std::size_t>(maximum_peak[0]);
    result.total_peak_resident_bytes = static_cast<std::size_t>(total_peak[0]);
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

struct MatrixInsertionWorkspace final {
    std::vector<PetscInt> rows, columns;
    std::vector<unsigned char> grouped_rows;
    std::vector<double> values;
};

PetscErrorCode insert_pattern_blocks(Mat matrix, const std::vector<PetscInt>& dofs,
    const std::vector<unsigned char>& pattern, const std::vector<double>& dense_values, InsertMode mode,
    PetscInt ownership_begin, PetscInt ownership_end, MatrixInsertionWorkspace& workspace) {
    PetscFunctionBeginUser;
    const std::size_t count = dofs.size();
    PetscCheck(pattern.size() == count * count, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
        "Local Jacobian pattern size does not match its DOFs");
    PetscCheck(dense_values.size() == count * count, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
        "Local Jacobian value count does not match its DOFs");
    workspace.grouped_rows.assign(count, 0U);
    for (std::size_t row = 0; row < count; ++row) {
        if (workspace.grouped_rows[row] != 0U) continue;
        if (dofs[row] < ownership_begin || dofs[row] >= ownership_end) {
            workspace.grouped_rows[row] = 1U;
            continue;
        }
        workspace.columns.clear();
        for (std::size_t column = 0; column < count; ++column)
            if (pattern[row * count + column] != 0U) workspace.columns.push_back(dofs[column]);
        workspace.grouped_rows[row] = 1U;
        if (workspace.columns.empty()) continue;
        workspace.rows.assign(1, dofs[row]);
        for (std::size_t candidate = row + 1U; candidate < count; ++candidate) {
            if (workspace.grouped_rows[candidate] != 0U || dofs[candidate] < ownership_begin ||
                dofs[candidate] >= ownership_end)
                continue;
            const auto first = pattern.begin() + static_cast<std::ptrdiff_t>(row * count);
            const auto candidate_first = pattern.begin() + static_cast<std::ptrdiff_t>(candidate * count);
            if (std::equal(first, first + static_cast<std::ptrdiff_t>(count), candidate_first)) {
                workspace.rows.push_back(dofs[candidate]);
                workspace.grouped_rows[candidate] = 1U;
            }
        }
        workspace.values.clear();
        workspace.values.reserve(workspace.rows.size() * workspace.columns.size());
        for (std::size_t local_row = 0; local_row < count; ++local_row) {
            if (workspace.grouped_rows[local_row] == 0U || dofs[local_row] < ownership_begin ||
                dofs[local_row] >= ownership_end)
                continue;
            if (!std::equal(pattern.begin() + static_cast<std::ptrdiff_t>(row * count),
                    pattern.begin() + static_cast<std::ptrdiff_t>((row + 1U) * count),
                    pattern.begin() + static_cast<std::ptrdiff_t>(local_row * count)))
                continue;
            for (std::size_t column = 0; column < count; ++column)
                if (pattern[row * count + column] != 0U)
                    workspace.values.push_back(dense_values[local_row * count + column]);
        }
        PetscCall(MatSetValues(matrix, static_cast<PetscInt>(workspace.rows.size()), workspace.rows.data(),
            static_cast<PetscInt>(workspace.columns.size()), workspace.columns.data(), workspace.values.data(), mode));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

struct SolverContext final {
    const NonlinearProblem* problem = nullptr;
    bool pattern_locked = false;
    PetscMPIInt rank = 0;
    PetscMPIInt size = 1;
    std::size_t contribution_begin = 0, contribution_end = 0;
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
    bool field_residual_convergence = false;
    double relative_tolerance = 1.0e-10, residual_reduction_tolerance = 1.0e-6,
           temperature_residual_absolute_tolerance = 1.0e-8, mechanical_residual_absolute_tolerance = 1.0e-4;
    bool first_residual = true, thermal_scaling_initialized = false, mechanics_scaling_initialized = false;
    std::vector<double> initial_field_residual_norms, field_residual_reference_norms;
    std::vector<double> latest_unscaled_field_residual_norms, latest_field_residual_norms;
    std::vector<double> field_residual_scalings, local_field_squared_norms;
    std::vector<double> global_field_squared_norms;
    std::vector<std::size_t> dof_fields;
    std::vector<PetscInt> problem_to_petsc;
    std::vector<std::size_t> petsc_to_problem;
    std::vector<PetscInt> petsc_contribution_dofs;
    std::vector<unsigned char> contribution_jacobian_pattern;
    ContributionWorkspace contribution_workspace;

    struct FixedContributionMetadata final {
        std::vector<std::size_t> dofs;
        std::vector<PetscInt> petsc_dofs;
        std::vector<unsigned char> jacobian_pattern;
    };

    std::vector<FixedContributionMetadata> fixed_contributions;
    MatrixInsertionWorkspace matrix_insertion;
    double initial_residual_norm = std::numeric_limits<double>::quiet_NaN();
    bool saw_domain_error = false, last_function_domain_error = false;
    std::string last_domain_error;
    SolveTiming timing;
};

PetscErrorCode assemble_contributions(SolverContext& context, Vec residual, Mat jacobian) {
    PetscFunctionBeginUser;
    const bool linearize = jacobian != nullptr;
    context.problem->validate_local_state(context.contribution_begin, context.contribution_end, context.state_values);
    for (std::size_t entry = context.contribution_begin; entry < context.contribution_end; ++entry) {
        const SolverContext::FixedContributionMetadata* fixed =
            context.fixed_contributions.empty() ? nullptr
                                                : &context.fixed_contributions.at(entry - context.contribution_begin);
        if (fixed == nullptr) context.problem->contribution_dofs(entry, context.contribution_workspace.dofs);
        const std::vector<std::size_t>& dofs = fixed == nullptr ? context.contribution_workspace.dofs : fixed->dofs;
        const std::size_t local_count = dofs.size();
        context.contribution_workspace.resize(local_count, linearize);
        for (std::size_t local = 0; local < local_count; ++local)
            context.contribution_workspace.state[local] = context.state_values[dofs[local]];
        if (fixed == nullptr) {
            context.petsc_contribution_dofs.resize(local_count);
            for (std::size_t local = 0; local < local_count; ++local)
                context.petsc_contribution_dofs[local] = context.problem_to_petsc[dofs[local]];
        }
        std::vector<double>* local_jacobian = linearize ? &context.contribution_workspace.jacobian : nullptr;
        context.problem->compute_contribution(
            entry, context.contribution_workspace.state, context.contribution_workspace.residual, local_jacobian);
        if (context.contribution_workspace.residual.size() != local_count ||
            (linearize && context.contribution_workspace.jacobian.size() != local_count * local_count))
            throw std::logic_error("NonlinearProblem contribution output has the wrong size");
        const PetscInt petsc_local_count = checked_petsc_int(local_count);
        const std::vector<PetscInt>& petsc_dof_vector =
            fixed == nullptr ? context.petsc_contribution_dofs : fixed->petsc_dofs;
        const PetscInt* petsc_dofs = petsc_dof_vector.data();
        if (!linearize) {
            PetscCall(VecSetValues(
                residual, petsc_local_count, petsc_dofs, context.contribution_workspace.residual.data(), ADD_VALUES));
            continue;
        }
        std::vector<double>& scaled_jacobian = context.contribution_workspace.jacobian;
        if (context.field_residual_scaling) {
            for (std::size_t row = 0; row < local_count; ++row) {
                const std::size_t global_row = dofs[row];
                if (context.constrained[global_row]) continue;
                const std::size_t field = context.dof_fields[global_row];
                for (std::size_t column = 0; column < local_count; ++column)
                    scaled_jacobian[row * local_count + column] *= context.field_residual_scalings[field];
            }
        }
        if (fixed == nullptr)
            context.problem->contribution_jacobian_pattern(entry, context.contribution_jacobian_pattern);
        const std::vector<unsigned char>& jacobian_pattern =
            fixed == nullptr ? context.contribution_jacobian_pattern : fixed->jacobian_pattern;
        PetscCall(insert_pattern_blocks(jacobian, petsc_dof_vector, jacobian_pattern, scaled_jacobian, ADD_VALUES, 0,
            checked_petsc_int(context.problem->dof_count()), context.matrix_insertion));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

struct PetscObjects final {
    SNES snes = nullptr;
    Vec state = nullptr;
    Vec residual = nullptr;
    Mat jacobian = nullptr;
    Vec gathered_state = nullptr;
    VecScatter state_scatter = nullptr;

    ~PetscObjects() {
        (void)VecScatterDestroy(&state_scatter);
        (void)VecDestroy(&gathered_state);
        (void)SNESDestroy(&snes);
        (void)VecDestroy(&state);
        (void)VecDestroy(&residual);
        (void)MatDestroy(&jacobian);
    }
};

void configure_linear_solver(PetscObjects& objects, const NonlinearProblem& problem, const SolverOptions& options,
    PetscMPIInt world_size, const SolverContext& context) {
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
    case SolverOptions::Preconditioner::automatic: throw std::logic_error("automatic preconditioner was not resolved");
    case SolverOptions::Preconditioner::lu: {
        check_petsc(PCSetType(preconditioner, PCLU), "PCSetType LU");
        std::array<char, 64> requested_solver{};
        PetscBool solver_set = PETSC_FALSE;
        check_petsc(PetscOptionsGetString(nullptr, nullptr, "-pc_factor_mat_solver_type", requested_solver.data(),
                        requested_solver.size(), &solver_set),
            "PetscOptionsGetString factor matrix solver");
        const bool default_mumps =
            world_size > 1 || options.direct_factorization == SolverOptions::DirectFactorization::mumps;
        if (solver_set == PETSC_FALSE && default_mumps)
            check_petsc(PCFactorSetMatSolverType(preconditioner, MATSOLVERMUMPS), "PCFactorSetMatSolverType MUMPS");
        const bool uses_mumps = solver_set == PETSC_TRUE
                                    ? std::string(requested_solver.data()) == std::string(MATSOLVERMUMPS)
                                    : default_mumps;
        if (uses_mumps) {
            // Preserve SCOTCH as the robust default. A case may select PORD after a controlled same-workload audit,
            // while an explicit PETSc command-line setting always wins.
            PetscBool ordering_set = PETSC_FALSE;
            check_petsc(PetscOptionsHasName(nullptr, nullptr, "-mat_mumps_icntl_7", &ordering_set),
                "PetscOptionsHasName MUMPS ordering");
            if (ordering_set == PETSC_FALSE) {
                const char* ordering = options.mumps_ordering == SolverOptions::MumpsOrdering::pord ? "4" : "3";
                check_petsc(PetscOptionsSetValue(nullptr, "-mat_mumps_icntl_7", ordering),
                    "PetscOptionsSetValue MUMPS ordering");
            }
            if (world_size >= 4) {
                PetscBool memory_relaxation_set = PETSC_FALSE;
                check_petsc(PetscOptionsHasName(nullptr, nullptr, "-mat_mumps_icntl_14", &memory_relaxation_set),
                    "PetscOptionsHasName MUMPS memory relaxation");
                if (memory_relaxation_set == PETSC_FALSE)
                    check_petsc(PetscOptionsSetValue(nullptr, "-mat_mumps_icntl_14", "100"),
                        "PetscOptionsSetValue MUMPS memory relaxation");
            }
        }
        check_petsc(PCFactorSetReuseOrdering(preconditioner, PETSC_TRUE), "PCFactorSetReuseOrdering");
        check_petsc(PCFactorSetReuseFill(preconditioner, PETSC_TRUE), "PCFactorSetReuseFill");
        break;
    }
    case SolverOptions::Preconditioner::block_jacobi:
        check_petsc(PCSetType(preconditioner, PCBJACOBI), "PCSetType BJACOBI");
        break;
    case SolverOptions::Preconditioner::field_split: {
        PetscBool already_configured = PETSC_FALSE;
        check_petsc(
            PetscObjectTypeCompare(reinterpret_cast<PetscObject>(preconditioner), PCFIELDSPLIT, &already_configured),
            "PetscObjectTypeCompare field split");
        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        check_petsc(
            VecGetOwnershipRange(objects.state, &ownership_begin, &ownership_end), "VecGetOwnershipRange field split");
        std::vector<PetscInt> thermal_indices, mechanical_indices;
        bool has_thermal_field = false, has_mechanical_field = false;
        for (const FieldDescriptor& field : problem.field_layout()) {
            has_thermal_field = has_thermal_field || field.category == FieldCategory::thermal;
            has_mechanical_field = has_mechanical_field || field.category == FieldCategory::mechanical;
        }
        for (PetscInt dof = ownership_begin; dof < ownership_end; ++dof) {
            const std::size_t problem_dof = context.petsc_to_problem[static_cast<std::size_t>(dof)];
            std::vector<PetscInt>& indices =
                problem.field_layout()[context.dof_fields[problem_dof]].category == FieldCategory::thermal
                    ? thermal_indices
                    : mechanical_indices;
            indices.push_back(dof);
        }
        if (!has_thermal_field || !has_mechanical_field)
            throw std::invalid_argument("field_split requires thermal and mechanical field metadata");
        check_petsc(PCSetType(preconditioner, PCFIELDSPLIT), "PCSetType FIELDSPLIT");
        if (already_configured == PETSC_FALSE) {
            IS temperature = nullptr;
            IS mechanics = nullptr;
            try {
                check_petsc(ISCreateGeneral(PETSC_COMM_WORLD, checked_petsc_int(thermal_indices.size()),
                                thermal_indices.data(), PETSC_COPY_VALUES, &temperature),
                    "ISCreateGeneral temperature");
                check_petsc(ISCreateGeneral(PETSC_COMM_WORLD, checked_petsc_int(mechanical_indices.size()),
                                mechanical_indices.data(), PETSC_COPY_VALUES, &mechanics),
                    "ISCreateGeneral mechanics");
                check_petsc(
                    PCFieldSplitSetIS(preconditioner, "temperature", temperature), "PCFieldSplitSetIS temperature");
                check_petsc(PCFieldSplitSetIS(preconditioner, "mechanics", mechanics), "PCFieldSplitSetIS mechanics");
                check_petsc(ISDestroy(&temperature), "ISDestroy temperature");
                check_petsc(ISDestroy(&mechanics), "ISDestroy mechanics");
            } catch (...) {
                if (temperature != nullptr) (void)ISDestroy(&temperature);
                if (mechanics != nullptr) (void)ISDestroy(&mechanics);
                throw;
            }
        }
        check_petsc(
            PCFieldSplitSetType(preconditioner, PC_COMPOSITE_MULTIPLICATIVE), "PCFieldSplitSetType multiplicative");
        if (world_size >= 4) {
            PetscBool mechanics_fill_set = PETSC_FALSE;
            check_petsc(PetscOptionsHasName(
                            nullptr, nullptr, "-fieldsplit_mechanics_sub_pc_factor_levels", &mechanics_fill_set),
                "PetscOptionsHasName mechanics field-split ILU fill level");
            if (mechanics_fill_set == PETSC_FALSE)
                check_petsc(PetscOptionsSetValue(nullptr, "-fieldsplit_mechanics_sub_pc_factor_levels", "1"),
                    "PetscOptionsSetValue mechanics field-split ILU fill level");
        }
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
    for (std::size_t index = 0; index < context.shadow_dofs.size(); ++index)
        context.state_values[context.shadow_dofs[index]] = PetscRealPart(values[checked_petsc_int(index)]);
    PetscCall(VecRestoreArrayRead(context.gathered_state, &values));
    PetscFunctionReturn(PETSC_SUCCESS);
}

std::vector<double> gather_complete_state(
    Vec state, std::size_t global_size, const std::vector<std::size_t>& petsc_to_problem) {
    VecScatter scatter = nullptr;
    Vec gathered = nullptr;
    check_petsc(VecScatterCreateToAll(state, &scatter, &gathered), "VecScatterCreateToAll final state");
    try {
        check_petsc(
            VecScatterBegin(scatter, state, gathered, INSERT_VALUES, SCATTER_FORWARD), "VecScatterBegin final state");
        check_petsc(
            VecScatterEnd(scatter, state, gathered, INSERT_VALUES, SCATTER_FORWARD), "VecScatterEnd final state");
        std::vector<double> result(global_size, 0.0);
        const PetscScalar* values = nullptr;
        check_petsc(VecGetArrayRead(gathered, &values), "VecGetArrayRead final state");
        for (std::size_t index = 0; index < global_size; ++index)
            result[petsc_to_problem[index]] = PetscRealPart(values[checked_petsc_int(index)]);
        check_petsc(VecRestoreArrayRead(gathered, &values), "VecRestoreArrayRead final state");
        check_petsc(VecScatterDestroy(&scatter), "VecScatterDestroy final state");
        check_petsc(VecDestroy(&gathered), "VecDestroy final state");
        return result;
    } catch (...) {
        if (scatter != nullptr) (void)VecScatterDestroy(&scatter);
        if (gathered != nullptr) (void)VecDestroy(&gathered);
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

PetscErrorCode field_norms(Vec vector, SolverContext& context, std::vector<double>& norms) {
    PetscFunctionBeginUser;
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    PetscCall(VecGetOwnershipRange(vector, &ownership_begin, &ownership_end));
    const PetscScalar* values = nullptr;
    PetscCall(VecGetArrayRead(vector, &values));
    context.local_field_squared_norms.assign(context.problem->field_layout().size(), 0.0);
    for (PetscInt global = ownership_begin; global < ownership_end; ++global) {
        const std::size_t problem_dof = context.petsc_to_problem[static_cast<std::size_t>(global)];
        const std::size_t field = context.dof_fields[problem_dof];
        const double value = PetscRealPart(values[global - ownership_begin]);
        context.local_field_squared_norms[field] += value * value;
    }
    PetscCall(VecRestoreArrayRead(vector, &values));
    context.global_field_squared_norms.assign(context.problem->field_layout().size(), 0.0);
    if (PetscGlobalSize == 1) {
        context.global_field_squared_norms = context.local_field_squared_norms;
    } else {
        PetscCallMPI(MPIU_Allreduce(context.local_field_squared_norms.data(), context.global_field_squared_norms.data(),
            static_cast<MPIU_Count>(context.local_field_squared_norms.size()), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD));
    }
    norms.resize(context.problem->field_layout().size());
    for (std::size_t field = 0; field < norms.size(); ++field)
        norms[field] = std::sqrt(context.global_field_squared_norms[field]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode scale_residual(Vec residual, SolverContext& context) {
    PetscFunctionBeginUser;
    PetscCall(field_norms(residual, context, context.latest_unscaled_field_residual_norms));
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
    double thermal_norm = 0.0, mechanics_norm = 0.0;
    for (std::size_t field = 0; field < context.problem->field_layout().size(); ++field)
        if (context.problem->field_layout()[field].category == FieldCategory::thermal)
            thermal_norm = std::hypot(thermal_norm, context.latest_unscaled_field_residual_norms[field]);
        else
            mechanics_norm = std::hypot(mechanics_norm, context.latest_unscaled_field_residual_norms[field]);
    if (!context.thermal_scaling_initialized && thermal_norm > context.residual_scaling_floor) {
        for (std::size_t field = 0; field < context.problem->field_layout().size(); ++field) {
            if (context.problem->field_layout()[field].category == FieldCategory::thermal) {
                context.field_residual_scalings[field] = 1.0 / thermal_norm;
                context.initial_field_residual_norms[field] = context.latest_unscaled_field_residual_norms[field];
            }
        }
        context.thermal_scaling_initialized = true;
    }
    if (!context.mechanics_scaling_initialized && mechanics_norm > context.residual_scaling_floor) {
        const double mechanics_scaling = 1.0 / mechanics_norm;
        for (std::size_t field = 0; field < context.problem->field_layout().size(); ++field) {
            if (context.problem->field_layout()[field].category == FieldCategory::mechanical) {
                context.field_residual_scalings[field] = mechanics_scaling;
                context.initial_field_residual_norms[field] = context.latest_unscaled_field_residual_norms[field];
            }
        }
        context.mechanics_scaling_initialized = true;
    }
    PetscInt ownership_begin = 0;
    PetscInt ownership_end = 0;
    PetscCall(VecGetOwnershipRange(residual, &ownership_begin, &ownership_end));
    PetscScalar* values = nullptr;
    PetscCall(VecGetArray(residual, &values));
    for (PetscInt global = ownership_begin; global < ownership_end; ++global) {
        const std::size_t problem_dof = context.petsc_to_problem[static_cast<std::size_t>(global)];
        if (!context.constrained[problem_dof]) {
            const std::size_t field = context.dof_fields[problem_dof];
            values[global - ownership_begin] *= context.field_residual_scalings[field];
        }
    }
    PetscCall(VecRestoreArray(residual, &values));
    PetscCall(field_norms(residual, context, context.latest_field_residual_norms));
    for (std::size_t field = 0; field < context.field_residual_reference_norms.size(); ++field)
        context.field_residual_reference_norms[field] =
            std::max(context.field_residual_reference_norms[field], context.latest_field_residual_norms[field]);
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode assemble_callback(SNES snes, Vec state, Vec residual, Mat jacobian, void* raw_context) {
    PetscFunctionBeginUser;
    const SteadyClock::time_point start = SteadyClock::now();
    try {
        SolverContext& context = *static_cast<SolverContext*>(raw_context);
        const bool linearize = jacobian != nullptr;
        double& callback_seconds =
            linearize ? context.timing.jacobian_callback_seconds : context.timing.residual_callback_seconds;
        std::size_t& callback_evaluations =
            linearize ? context.timing.jacobian_evaluations : context.timing.residual_evaluations;
        if (!linearize) context.last_function_domain_error = false;
        const NonlinearProblem& problem = *context.problem;
        PetscCall(gather_state(state, context));
        if (linearize)
            PetscCall(MatZeroEntries(jacobian));
        else
            PetscCall(VecSet(residual, 0.0));
        bool local_domain_error = false;
        const SteadyClock::time_point assembly_start = SteadyClock::now();
        try {
            PetscCall(assemble_contributions(context, residual, jacobian));
        } catch (const std::domain_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        } catch (const std::overflow_error& error) {
            local_domain_error = true;
            context.last_domain_error = error.what();
        }
        const double local_assembly_seconds = seconds_since(assembly_start);
        if (linearize)
            context.timing.maximum_jacobian_assembly_seconds += local_assembly_seconds;
        else
            context.timing.maximum_residual_assembly_seconds += local_assembly_seconds;
        bool global_domain_error = false;
        PetscCall(synchronize_domain_error(local_domain_error, global_domain_error));
        if (global_domain_error && context.last_domain_error.empty())
            context.last_domain_error = linearize
                                            ? "a Jacobian evaluation violated its physical domain on another MPI rank"
                                            : "a residual evaluation violated its physical domain on another MPI rank";
        if (linearize) {
            PetscCall(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));
        } else {
            PetscCall(VecAssemblyBegin(residual));
            PetscCall(VecAssemblyEnd(residual));
        }
        if (global_domain_error) {
            context.saw_domain_error = true;
            if (linearize) {
                PetscCall(MatZeroEntries(jacobian));
                PetscCall(SNESSetJacobianDomainError(snes));
            } else {
                context.last_function_domain_error = true;
                PetscCall(VecSet(residual, 0.0));
                PetscCall(SNESSetFunctionDomainError(snes));
            }
            callback_seconds += seconds_since(start);
            ++callback_evaluations;
            PetscFunctionReturn(PETSC_SUCCESS);
        }
        if (linearize) {
            if (!context.pattern_locked) {
                PetscCall(MatSetOption(jacobian, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
                if (!context.problem->jacobian_sparsity_is_state_dependent())
                    PetscCall(MatSetOption(jacobian, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));
                context.pattern_locked = true;
            }
            PetscCall(MatZeroRows(jacobian, static_cast<PetscInt>(context.constrained_dofs.size()),
                context.constrained_dofs.data(), 1.0, nullptr, nullptr));
        } else {
            PetscInt ownership_begin = 0;
            PetscInt ownership_end = 0;
            PetscCall(VecGetOwnershipRange(residual, &ownership_begin, &ownership_end));
            PetscScalar* local_residual = nullptr;
            PetscCall(VecGetArray(residual, &local_residual));
            for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
                const PetscInt dof = context.problem_to_petsc[condition.dof];
                if (dof >= ownership_begin && dof < ownership_end)
                    local_residual[dof - ownership_begin] = context.state_values[condition.dof] - condition.value;
            }
            PetscCall(VecRestoreArray(residual, &local_residual));
            PetscCall(scale_residual(residual, context));
            if (!std::isfinite(context.initial_residual_norm)) {
                PetscReal norm = 0.0;
                PetscCall(VecNorm(residual, NORM_2, &norm));
                context.initial_residual_norm = static_cast<double>(norm);
            }
        }
        callback_seconds += seconds_since(start);
        ++callback_evaluations;
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "fuelsim nonlinear assembly failed: %s", error.what());
    } catch (...) { SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "fuelsim nonlinear assembly failed with unknown error"); }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_function(SNES snes, Vec state, Vec residual, void* raw_context) {
    return assemble_callback(snes, state, residual, nullptr, raw_context);
}

PetscErrorCode form_jacobian(SNES snes, Vec state, Mat jacobian, Mat preconditioner, void* raw_context) {
    PetscFunctionBeginUser;
    PetscCheck(jacobian == preconditioner, PETSC_COMM_WORLD, PETSC_ERR_SUP, "fuelsim requires one matrix for J and P");
    PetscCall(assemble_callback(snes, state, nullptr, jacobian, raw_context));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode field_residual_convergence_test(SNES snes, PetscInt iteration, PetscReal state_norm, PetscReal step_norm,
    PetscReal residual_norm, SNESConvergedReason* reason, void* raw_context) {
    PetscFunctionBeginUser;
    PetscCall(SNESConvergedDefault(snes, iteration, state_norm, step_norm, residual_norm, reason, nullptr));
    if (*reason != SNES_CONVERGED_ITERATING) PetscFunctionReturn(PETSC_SUCCESS);
    const SolverContext& context = *static_cast<const SolverContext*>(raw_context);
    if (!context.field_residual_convergence || context.last_function_domain_error || iteration == 0)
        PetscFunctionReturn(PETSC_SUCCESS);
    constexpr double numerical_residual_floor = 10.0 * std::sqrt(std::numeric_limits<double>::epsilon());
    bool converged = true;
    for (std::size_t field = 0; field < context.latest_field_residual_norms.size(); ++field) {
        const double initial_scaled =
            context.initial_field_residual_norms[field] * context.field_residual_scalings[field];
        const double reference = std::max(initial_scaled, context.field_residual_reference_norms[field]);
        const double physical_tolerance = context.problem->field_layout()[field].category == FieldCategory::thermal
                                              ? context.temperature_residual_absolute_tolerance
                                              : context.mechanical_residual_absolute_tolerance;
        const double configured_threshold = std::max(
            physical_tolerance * context.field_residual_scalings[field], context.relative_tolerance * reference);
        const double independent_threshold =
            std::max(numerical_residual_floor, context.residual_reduction_tolerance * reference);
        const double threshold = std::max(configured_threshold, independent_threshold) *
                                 (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
        converged = std::isfinite(context.latest_field_residual_norms[field]) &&
                    context.latest_field_residual_norms[field] <= threshold && converged;
    }
    if (converged) *reason = SNES_CONVERGED_FNORM_ABS;
    PetscFunctionReturn(PETSC_SUCCESS);
}
} // namespace

class PetscSolver::Implementation final {
  public:
    bool prepare(const NonlinearProblem& problem) {
        const PetscInt requested_count = checked_petsc_int(problem.dof_count());
        const std::shared_ptr<const void> requested_identity = problem.discretization_identity();
        if (_problem_identity == requested_identity) {
            _context.problem = &problem;
            return false;
        }
        problem.validate_discretization();
        _problem_identity.reset();
        _objects = std::make_unique<PetscObjects>();
        _count = requested_count;
        _context = SolverContext{};
        _context.problem = &problem;
        _context.rank = PetscGlobalRank;
        _context.size = PetscGlobalSize;
        const std::size_t rank = static_cast<std::size_t>(_context.rank),
                          size = static_cast<std::size_t>(_context.size);
        const auto contribution_partition = problem.contribution_partition(rank, size);
        _context.contribution_begin = contribution_partition.first;
        _context.contribution_end = contribution_partition.second;
        _context.constrained.assign(problem.dof_count(), false);
        const std::size_t field_count = problem.field_layout().size();
        _context.initial_field_residual_norms.resize(field_count);
        _context.field_residual_reference_norms.resize(field_count);
        _context.latest_unscaled_field_residual_norms.resize(field_count);
        _context.latest_field_residual_norms.resize(field_count);
        _context.field_residual_scalings.resize(field_count, 1.0);
        _context.local_field_squared_norms.resize(field_count);
        _context.global_field_squared_norms.resize(field_count);
        _context.dof_fields.resize(problem.dof_count());
        for (std::size_t field = 0; field < field_count; ++field) {
            const FieldDescriptor& descriptor = problem.field_layout()[field];
            std::fill(_context.dof_fields.begin() + static_cast<std::ptrdiff_t>(descriptor.begin),
                _context.dof_fields.begin() + static_cast<std::ptrdiff_t>(descriptor.end), field);
        }
        _context.problem_to_petsc.resize(problem.dof_count());
        _context.petsc_to_problem.resize(problem.dof_count());
        bool interleave_fields = field_count > 1U;
        const std::size_t field_size = field_count == 0U ? 0U : problem.field_layout().front().end;
        for (std::size_t field = 0; field < field_count; ++field) {
            const FieldDescriptor& descriptor = problem.field_layout()[field];
            interleave_fields = interleave_fields && descriptor.begin == field * field_size &&
                                descriptor.end == (field + 1U) * field_size;
        }
        for (std::size_t problem_dof = 0; problem_dof < problem.dof_count(); ++problem_dof) {
            const std::size_t petsc_dof =
                interleave_fields ? (problem_dof % field_size) * field_count + problem_dof / field_size : problem_dof;
            _context.problem_to_petsc[problem_dof] = checked_petsc_int(petsc_dof);
            _context.petsc_to_problem[petsc_dof] = problem_dof;
        }
        std::size_t maximum_local_dofs = 0;
        std::vector<std::size_t> contribution_dofs;
        const bool fixed_contributions = problem.contribution_metadata_is_fixed();
        if (fixed_contributions)
            _context.fixed_contributions.resize(_context.contribution_end - _context.contribution_begin);
        for (std::size_t entry = _context.contribution_begin; entry < _context.contribution_end; ++entry) {
            problem.contribution_dofs(entry, contribution_dofs);
            maximum_local_dofs = std::max(maximum_local_dofs, contribution_dofs.size());
            for (const std::size_t dof : contribution_dofs)
                if (dof >= problem.dof_count())
                    throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
            if (fixed_contributions) {
                SolverContext::FixedContributionMetadata& fixed =
                    _context.fixed_contributions.at(entry - _context.contribution_begin);
                fixed.dofs = contribution_dofs;
                fixed.petsc_dofs.resize(contribution_dofs.size());
                for (std::size_t local = 0; local < contribution_dofs.size(); ++local)
                    fixed.petsc_dofs[local] = _context.problem_to_petsc[contribution_dofs[local]];
                problem.contribution_jacobian_pattern(entry, fixed.jacobian_pattern);
            }
        }
        _context.contribution_workspace.reserve(maximum_local_dofs);
        _context.petsc_contribution_dofs.reserve(maximum_local_dofs);
        _context.contribution_jacobian_pattern.reserve(maximum_local_dofs * maximum_local_dofs);
        _context.matrix_insertion.rows.reserve(maximum_local_dofs);
        _context.matrix_insertion.columns.reserve(maximum_local_dofs);
        _context.matrix_insertion.grouped_rows.reserve(maximum_local_dofs);
        _context.matrix_insertion.values.reserve(maximum_local_dofs * maximum_local_dofs);
        const PetscInt requested_local_count =
            interleave_fields
                ? checked_petsc_int(field_count * (field_size * (rank + 1U) / size - field_size * rank / size))
                : PETSC_DECIDE;
        check_petsc(
            VecCreateMPI(PETSC_COMM_WORLD, requested_local_count, _count, &_objects->state), "VecCreateMPI state");
        check_petsc(VecDuplicate(_objects->state, &_objects->residual), "VecDuplicate residual");
        PetscInt local_count = 0;
        check_petsc(VecGetLocalSize(_objects->state, &local_count), "VecGetLocalSize state");
        PetscInt ownership_begin = 0;
        PetscInt ownership_end = 0;
        check_petsc(
            VecGetOwnershipRange(_objects->state, &ownership_begin, &ownership_end), "VecGetOwnershipRange state");
        const std::size_t sparsity_count = problem.sparsity_contribution_count();
        std::vector<PetscInt> sparsity_dofs;
        std::vector<double> sparsity_zeros;
        std::vector<unsigned char> sparsity_pattern;
        MatrixInsertionWorkspace sparsity_insertion;
        Mat preallocator = nullptr;
        try {
            check_petsc(MatCreate(PETSC_COMM_WORLD, &preallocator), "MatCreate sparsity preallocator");
            check_petsc(MatSetSizes(preallocator, local_count, local_count, _count, _count),
                "MatSetSizes sparsity preallocator");
            check_petsc(MatSetType(preallocator, MATPREALLOCATOR), "MatSetType sparsity preallocator");
            check_petsc(MatSetUp(preallocator), "MatSetUp sparsity preallocator");
            check_petsc(MatSetOption(preallocator, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE),
                "MatSetOption sparsity preallocator MAT_NO_OFF_PROC_ENTRIES");
            // Every rank sees the same graph, but inserts only its owned rows. This keeps
            // preallocation communication-free without changing contribution ownership.
            for (std::size_t entry = 0; entry < sparsity_count; ++entry) {
                problem.sparsity_contribution_dofs(entry, contribution_dofs);
                sparsity_dofs.resize(contribution_dofs.size());
                for (std::size_t local = 0; local < contribution_dofs.size(); ++local)
                    sparsity_dofs[local] = _context.problem_to_petsc[contribution_dofs[local]];
                sparsity_zeros.assign(contribution_dofs.size() * contribution_dofs.size(), 0.0);
                problem.sparsity_contribution_jacobian_pattern(entry, sparsity_pattern);
                check_petsc(insert_pattern_blocks(preallocator, sparsity_dofs, sparsity_pattern, sparsity_zeros,
                                INSERT_VALUES, ownership_begin, ownership_end, sparsity_insertion),
                    "insert_pattern_blocks sparsity preallocator");
            }
            check_petsc(MatAssemblyBegin(preallocator, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin sparsity preallocator");
            check_petsc(MatAssemblyEnd(preallocator, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd sparsity preallocator");
            check_petsc(MatCreate(PETSC_COMM_WORLD, &_objects->jacobian), "MatCreate Jacobian");
            check_petsc(
                MatSetSizes(_objects->jacobian, local_count, local_count, _count, _count), "MatSetSizes Jacobian");
            check_petsc(MatSetType(_objects->jacobian, MATAIJ), "MatSetType Jacobian");
            check_petsc(
                MatPreallocatorPreallocate(preallocator,
                    problem.jacobian_sparsity_is_state_dependent() ? PETSC_FALSE : PETSC_TRUE, _objects->jacobian),
                "MatPreallocatorPreallocate Jacobian");
            check_petsc(MatDestroy(&preallocator), "MatDestroy sparsity preallocator");
        } catch (...) {
            if (preallocator != nullptr) (void)MatDestroy(&preallocator);
            throw;
        }
        check_petsc(MatSetOption(_objects->jacobian, MAT_NEW_NONZERO_ALLOCATION_ERR,
                        problem.jacobian_sparsity_is_state_dependent() ? PETSC_FALSE : PETSC_TRUE),
            "MatSetOption MAT_NEW_NONZERO_ALLOCATION_ERR");
        check_petsc(SNESCreate(PETSC_COMM_WORLD, &_objects->snes), "SNESCreate");
        check_petsc(SNESSetFunction(_objects->snes, _objects->residual, form_function, &_context), "SNESSetFunction");
        check_petsc(SNESSetJacobian(_objects->snes, _objects->jacobian, _objects->jacobian, form_jacobian, &_context),
            "SNESSetJacobian");
        check_petsc(SNESSetType(_objects->snes, SNESNEWTONLS), "SNESSetType");
        SNESLineSearch line_search = nullptr;
        check_petsc(SNESGetLineSearch(_objects->snes, &line_search), "SNESGetLineSearch");
        check_petsc(SNESLineSearchSetType(line_search, SNESLINESEARCHBASIC), "SNESLineSearchSetType");
        _context.ownership_begin = ownership_begin;
        _context.ownership_end = ownership_end;
        for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
            if (condition.dof >= problem.dof_count())
                throw std::out_of_range("Dirichlet condition DOF is out of range");
            const PetscInt dof = _context.problem_to_petsc[condition.dof];
            if (dof >= ownership_begin && dof < ownership_end) _context.constrained_dofs.push_back(dof);
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
        for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
            const PetscInt dof = _context.problem_to_petsc[condition.dof];
            if (dof >= ownership_begin && dof < ownership_end)
                _context.shadow_dofs.push_back(static_cast<std::uint32_t>(condition.dof));
        }
        std::sort(_context.shadow_dofs.begin(), _context.shadow_dofs.end());
        _context.shadow_dofs.erase(
            std::unique(_context.shadow_dofs.begin(), _context.shadow_dofs.end()), _context.shadow_dofs.end());
        _context.state_values.resize(problem.dof_count(), std::numeric_limits<double>::quiet_NaN());
        std::vector<PetscInt> shadow_indices;
        shadow_indices.reserve(_context.shadow_dofs.size());
        for (const std::uint32_t dof : _context.shadow_dofs)
            shadow_indices.push_back(_context.problem_to_petsc[static_cast<std::size_t>(dof)]);
        const PetscInt shadow_count = checked_petsc_int(_context.shadow_dofs.size());
        IS source_indices = nullptr;
        IS destination_indices = nullptr;
        try {
            check_petsc(ISCreateGeneral(
                            PETSC_COMM_SELF, shadow_count, shadow_indices.data(), PETSC_COPY_VALUES, &source_indices),
                "ISCreateGeneral shadow state");
            check_petsc(ISCreateStride(PETSC_COMM_SELF, shadow_count, 0, 1, &destination_indices),
                "ISCreateStride shadow state");
            check_petsc(
                VecCreateSeq(PETSC_COMM_SELF, shadow_count, &_objects->gathered_state), "VecCreateSeq shadow state");
            check_petsc(VecScatterCreate(_objects->state, source_indices, _objects->gathered_state, destination_indices,
                            &_objects->state_scatter),
                "VecScatterCreate shadow state");
            check_petsc(ISDestroy(&source_indices), "ISDestroy shadow source");
            check_petsc(ISDestroy(&destination_indices), "ISDestroy shadow destination");
        } catch (...) {
            if (source_indices != nullptr) (void)ISDestroy(&source_indices);
            if (destination_indices != nullptr) (void)ISDestroy(&destination_indices);
            throw;
        }
        _context.state_scatter = _objects->state_scatter;
        _context.gathered_state = _objects->gathered_state;
        _problem_identity = requested_identity;
        return true;
    }

    PetscObjects& objects() { return *_objects; }

    SolverContext& context() noexcept { return _context; }

  private:
    std::shared_ptr<const void> _problem_identity;
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
    check_petsc(PetscMemorySetGetMaximumUsage(), "PetscMemorySetGetMaximumUsage");
    _rank = static_cast<int>(PetscGlobalRank);
    _size = static_cast<int>(PetscGlobalSize);
}

PetscSession::~PetscSession() {
    if (!_owns_initialization) return;
    PetscBool finalized = PETSC_FALSE;
    if (PetscFinalized(&finalized) == PETSC_SUCCESS && finalized == PETSC_FALSE) {
        const PetscErrorCode code = PetscFinalize();
        (void)code;
    }
}

int PetscSession::rank() const noexcept { return _rank; }

int PetscSession::size() const noexcept { return _size; }

void PetscSession::collective_root_action(const std::function<void()>& action) const {
    if (!action) throw std::invalid_argument("PetscSession collective root action must not be empty");
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
    const PetscInt local_failure[2] = {failed ? 1 : 0, failed ? 1 : 0};
    PetscInt collective_failure[2]{};
    check_petsc(PetscGlobalMinMaxInt(PETSC_COMM_WORLD, local_failure, collective_failure),
        "PetscGlobalMinMaxInt root I/O status");
    if (collective_failure[1] == 0) return;
    throw std::runtime_error(_rank == 0 ? "collective root-rank I/O failed: " + message
                                        : "collective root-rank I/O failed; see rank 0 for details");
}

PetscSolver::PetscSolver() : _impl(std::make_unique<Implementation>()) {}

PetscSolver::~PetscSolver() = default;

SolveResult PetscSolver::solve(
    const NonlinearProblem& problem, const std::vector<double>& initial_state, const SolverOptions& options) {
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

SolveResult PetscSolver::solve_once(
    const NonlinearProblem& problem, const std::vector<double>& initial_state, const SolverOptions& options) {
    if (initial_state.size() != problem.dof_count())
        throw std::invalid_argument("PetscSolver initial state size mismatch");
    const bool fixed_temperature_scale = options.temperature_residual_scale > 0.0,
               fixed_mechanical_scale = options.mechanical_residual_scale > 0.0;
    if (fixed_temperature_scale != fixed_mechanical_scale)
        throw std::invalid_argument("fixed residual scaling requires paired temperature and mechanical scales");
    if (options.field_residual_scaling && fixed_temperature_scale)
        throw std::invalid_argument("automatic and fixed residual scaling are mutually exclusive");
    const bool residual_scaling = options.field_residual_scaling || fixed_temperature_scale;
    const SteadyClock::time_point total_start = SteadyClock::now();
    const SteadyClock::time_point setup_start = SteadyClock::now();
    const MemorySnapshot initial_memory = memory_snapshot();
    const bool workspace_created = _impl->prepare(problem);
    PetscObjects& objects = _impl->objects();
    SolverContext& context = _impl->context();
    context.timing = SolveTiming{};
    context.timing.workspace_setups = workspace_created ? 1U : 0U;
    context.timing.solve_calls = 1;
    context.timing.initial_resident_bytes = static_cast<std::size_t>(initial_memory.resident_bytes);
    context.timing.maximum_peak_resident_bytes = static_cast<std::size_t>(initial_memory.maximum_resident_bytes);
    context.initial_residual_norm = std::numeric_limits<double>::quiet_NaN();
    context.field_residual_scaling = residual_scaling;
    context.residual_scaling_floor = options.absolute_tolerance;
    context.field_residual_convergence = options.field_residual_convergence;
    context.relative_tolerance = options.relative_tolerance;
    context.residual_reduction_tolerance = options.residual_reduction_tolerance;
    context.temperature_residual_absolute_tolerance = options.temperature_residual_absolute_tolerance;
    context.mechanical_residual_absolute_tolerance = options.mechanical_residual_absolute_tolerance;
    context.first_residual = true;
    context.thermal_scaling_initialized = fixed_temperature_scale;
    context.mechanics_scaling_initialized = fixed_mechanical_scale;
    const std::size_t field_count = problem.field_layout().size();
    context.initial_field_residual_norms.assign(field_count, 0.0);
    context.field_residual_reference_norms.assign(field_count, 0.0);
    context.latest_unscaled_field_residual_norms.assign(field_count, 0.0);
    context.latest_field_residual_norms.assign(field_count, 0.0);
    context.field_residual_scalings.assign(field_count, 1.0);
    if (fixed_temperature_scale) {
        for (std::size_t field = 0; field < field_count; ++field) {
            context.field_residual_scalings[field] = problem.field_layout()[field].category == FieldCategory::thermal
                                                         ? 1.0 / options.temperature_residual_scale
                                                         : 1.0 / options.mechanical_residual_scale;
        }
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
        state_array[index - ownership_begin] = initial_state[context.petsc_to_problem[static_cast<std::size_t>(index)]];
    check_petsc(VecRestoreArray(objects.state, &state_array), "VecRestoreArray state");
    check_petsc(SNESSetTolerances(objects.snes, options.absolute_tolerance, options.relative_tolerance,
                    residual_scaling ? 0.0 : options.step_tolerance, options.maximum_iterations, PETSC_DEFAULT),
        "SNESSetTolerances");
    check_petsc(SNESSetLagJacobian(objects.snes, options.jacobian_lag), "SNESSetLagJacobian");
    SNESLineSearch line_search = nullptr;
    check_petsc(SNESGetLineSearch(objects.snes, &line_search), "SNESGetLineSearch");
    check_petsc(
        SNESLineSearchSetType(line_search,
            options.line_search == SolverOptions::LineSearch::backtracking ? SNESLINESEARCHBT : SNESLINESEARCHBASIC),
        "SNESLineSearchSetType");
    configure_linear_solver(objects, problem, options, PetscGlobalSize, context);
    check_petsc(SNESSetFromOptions(objects.snes), "SNESSetFromOptions");
    check_petsc(SNESSetConvergenceTest(objects.snes,
                    options.field_residual_convergence ? field_residual_convergence_test : SNESConvergedDefault,
                    options.field_residual_convergence ? static_cast<void*>(&context) : nullptr, nullptr),
        "SNESSetConvergenceTest");
    const MemorySnapshot setup_memory = memory_snapshot();
    context.timing.setup_resident_bytes = static_cast<std::size_t>(setup_memory.resident_bytes);
    context.timing.maximum_peak_resident_bytes = std::max(
        context.timing.maximum_peak_resident_bytes, static_cast<std::size_t>(setup_memory.maximum_resident_bytes));
    context.timing.setup_seconds = seconds_since(setup_start);
    const SteadyClock::time_point solve_start = SteadyClock::now();
    check_petsc(SNESSolve(objects.snes, nullptr, objects.state), "SNESSolve");
    const MemorySnapshot solve_memory = memory_snapshot();
    context.timing.solve_resident_bytes = static_cast<std::size_t>(solve_memory.resident_bytes);
    context.timing.maximum_peak_resident_bytes = std::max(
        context.timing.maximum_peak_resident_bytes, static_cast<std::size_t>(solve_memory.maximum_resident_bytes));
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
        check_petsc(
            SNESComputeFunction(objects.snes, objects.state, objects.residual), "SNESComputeFunction final residual");
        check_petsc(VecNorm(objects.residual, NORM_2, &residual_norm), "VecNorm final residual");
    }
    const bool final_domain_error = context.last_function_domain_error;
    std::vector<double> solution = gather_complete_state(objects.state, problem.dof_count(), context.petsc_to_problem);
    const MemorySnapshot final_memory = memory_snapshot();
    context.timing.final_resident_bytes = static_cast<std::size_t>(final_memory.resident_bytes);
    context.timing.maximum_peak_resident_bytes = std::max(
        context.timing.maximum_peak_resident_bytes, static_cast<std::size_t>(final_memory.maximum_resident_bytes));
    SolveResult result;
    result.state = std::move(solution);
    result.nonlinear_iterations = static_cast<int>(iterations);
    result.linear_iterations = static_cast<int>(linear_iterations);
    result.residual_norm = static_cast<double>(residual_norm);
    result.convergence_reason = static_cast<int>(reason);
    result.field_names.reserve(problem.field_layout().size());
    for (const FieldDescriptor& field : problem.field_layout()) result.field_names.push_back(field.name);
    result.initial_field_residual_norms = context.initial_field_residual_norms;
    result.field_residual_reference_norms = context.field_residual_reference_norms;
    result.final_field_residual_norms = context.latest_unscaled_field_residual_norms;
    result.final_scaled_field_residual_norms = context.latest_field_residual_norms;
    result.field_residual_scalings = context.field_residual_scalings;
    context.timing.total_seconds = seconds_since(total_start);
    const double configured_residual_threshold =
        std::max(options.absolute_tolerance, options.relative_tolerance * context.initial_residual_norm);
    const double fallback_reduction = options.residual_reduction_tolerance,
                 numerical_residual_floor = 10.0 * std::sqrt(std::numeric_limits<double>::epsilon());
    const double independently_verified_threshold =
        std::max(numerical_residual_floor, fallback_reduction * context.initial_residual_norm);
    double physical_absolute_threshold = 0.0;
    for (std::size_t field = 0; field < problem.field_layout().size(); ++field) {
        const double tolerance = problem.field_layout()[field].category == FieldCategory::thermal
                                     ? options.temperature_residual_absolute_tolerance
                                     : options.mechanical_residual_absolute_tolerance;
        physical_absolute_threshold =
            std::hypot(physical_absolute_threshold, tolerance * result.field_residual_scalings[field]);
    }
    const double residual_threshold =
        std::max({configured_residual_threshold, independently_verified_threshold, physical_absolute_threshold});
    const double residual_slack = residual_threshold * (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
    bool fields_verified = true;
    std::vector<double> field_thresholds(result.final_scaled_field_residual_norms.size());
    for (std::size_t field = 0; field < result.final_scaled_field_residual_norms.size(); ++field) {
        const double initial_scaled =
            result.initial_field_residual_norms[field] * result.field_residual_scalings[field];
        const double reference = std::max(initial_scaled, result.field_residual_reference_norms[field]);
        const double configured_field_threshold =
            std::max(problem.field_layout()[field].category == FieldCategory::thermal
                         ? options.temperature_residual_absolute_tolerance * result.field_residual_scalings[field]
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
        if (!context.last_domain_error.empty())
            result.failure_message = context.last_domain_error;
        else if (context.saw_domain_error)
            result.failure_message = "a residual or Jacobian evaluation violated its physical domain";
        else
            result.failure_message = petsc_convergence_reason_name(result.convergence_reason);
    } else if (!result.converged && !residual_verified) {
        result.failure_category = SolveFailureCategory::residual_verification;
        std::ostringstream message;
        message << "PETSc reported convergence without satisfying residual tolerances: global=" << result.residual_norm
                << '/' << residual_slack;
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
        const PetscInt petsc_dof = context.problem_to_petsc[static_cast<std::size_t>(dof)];
        if (petsc_dof < context.ownership_begin || petsc_dof >= context.ownership_end) ++local_remote_shadow;
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
    if (name == nullptr) return "UNKNOWN";
    return name;
}

const char* solve_failure_category_name(SolveFailureCategory category) noexcept {
    switch (category) {
    case SolveFailureCategory::none: return "none";
    case SolveFailureCategory::nonlinear_divergence: return "nonlinear_divergence";
    case SolveFailureCategory::physical_domain: return "physical_domain";
    case SolveFailureCategory::residual_verification: return "residual_verification";
    case SolveFailureCategory::time_discretization: return "time_discretization";
    case SolveFailureCategory::contact_constraint: return "contact_constraint";
    }
    return "unknown";
}

} // namespace fuelsim
