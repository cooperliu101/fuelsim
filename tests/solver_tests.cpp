#include "core/cax4_evaluation.hpp"
#include "core/nonlinear_problem.hpp"
#include "core/steady_problem.hpp"
#include "solver/petsc_solver.hpp"
#include "solver/solve_workflows.hpp"
#include "support/jacobian_check.hpp"
#include "support/mesh_fixture.hpp"
#include "support/rz_problem_access.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double metric_relative_error(double actual, double expected) {
    return std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-30);
}

fuelsim::ThermoelasticProperties constant_material(double conductivity, double thermal_expansion) {
    return fuelsim::test::thermoelastic(0.0, conductivity, 75.0e9, 0.3, thermal_expansion, 600.0);
}

fuelsim::BoundaryConditionDefinition
dirichlet(const std::string& name, const std::string& boundary, fuelsim::Field field, double value) {
    fuelsim::BoundaryConditionDefinition condition{};
    condition.name = name;
    condition.type = fuelsim::BoundaryConditionType::dirichlet;
    condition.boundary = boundary;
    condition.field = field;
    condition.value = value;
    return condition;
}

fuelsim::BoundaryConditionDefinition pressure(const std::string& name, const std::string& boundary, double value) {
    fuelsim::BoundaryConditionDefinition condition{};
    condition.name = name;
    condition.type = fuelsim::BoundaryConditionType::pressure;
    condition.boundary = boundary;
    condition.field = fuelsim::Field::radial_displacement;
    condition.value = value;
    return condition;
}

fuelsim::SpatialDefinition single_region_definition(const fuelsim::ThermoelasticProperties& material,
    double heat_source,
    double initial_temperature,
    double outer_temperature,
    double inner_radius,
    double inner_pressure,
    double outer_pressure) {
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", material, heat_source, initial_temperature});
    if (inner_radius == 0.0)
        definition.boundary_conditions.push_back(
            dirichlet("inner_radial", "solid_inner", fuelsim::Field::radial_displacement, 0.0));
    definition.boundary_conditions.push_back(
        dirichlet("bottom_axial", "solid_bottom", fuelsim::Field::axial_displacement, 0.0));
    definition.boundary_conditions.push_back(
        dirichlet("outer_temperature", "solid_outer", fuelsim::Field::temperature, outer_temperature));
    if (inner_pressure > 0.0)
        definition.boundary_conditions.push_back(pressure("inner_pressure", "solid_inner", inner_pressure));
    if (outer_pressure > 0.0)
        definition.boundary_conditions.push_back(pressure("outer_pressure", "solid_outer", outer_pressure));
    return definition;
}

class TwelveDofProblem : public fuelsim::NonlinearProblem {
  public:
    std::size_t dof_count() const noexcept override { return fuelsim::cax4_local_dof_count; }

    std::size_t contribution_count() const noexcept override { return 1; }

    const std::vector<fuelsim::FieldDescriptor>& field_layout() const noexcept override { return _fields; }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override {
        validate_contribution(index);
        dofs.resize(fuelsim::cax4_local_dof_count);
        for (std::size_t dof = 0; dof < dofs.size(); ++dof)
            dofs[dof] = dof;
    }

    const std::vector<fuelsim::DirichletCondition>& dirichlet_conditions() const noexcept override {
        return _conditions;
    }

  protected:
    static void validate_contribution(std::size_t index) {
        if (index != 0)
            throw std::out_of_range("TwelveDofProblem contribution index");
    }

  private:
    std::vector<fuelsim::FieldDescriptor> _fields = {
        {"temperature", 0, 4, fuelsim::FieldCategory::thermal},
        {"radial", 4, 8, fuelsim::FieldCategory::mechanical},
        {"axial", 8, 12, fuelsim::FieldCategory::mechanical},
    };
    std::vector<fuelsim::DirichletCondition> _conditions;
};

class ChangingCouplingProblem final : public TwelveDofProblem {
  public:
    void set_coupling(double value) { _coupling = value; }

    bool jacobian_sparsity_is_state_dependent() const noexcept override { return true; }

    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const override {
        validate_contribution(index);
        residual.resize(state.size());
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            residual[dof] = 1.0e-20 * (state[dof] - 1.0);
        // The uniformly scaled system is well-conditioned. Zero filtering
        // must preserve its tiny coefficients, including a remote coupling.
        residual[0] += 0.5 * _coupling * (state.back() * state.back() - 1.0);
        if (!jacobian)
            return;
        jacobian->assign(state.size() * state.size(), 0.0);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            (*jacobian)[dof * state.size() + dof] = 1.0e-20;
        (*jacobian)[state.size() - 1] = _coupling * state.back();
    }

  private:
    double _coupling = 0.0;
};

bool test_changing_direct_coupling() {
    ChangingCouplingProblem problem;
    fuelsim::PetscSolver solver;
    fuelsim::SolverOptions options;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.direct_factorization = fuelsim::SolverOptions::DirectFactorization::mumps;
    options.backtracking_fallback = false;
    options.absolute_tolerance = 1.0e-30;
    bool passed = true;
    for (const double coupling : {0.0, 1.0e-20, 0.0, -2.0e-20}) {
        problem.set_coupling(coupling);
        const auto result = solver.solve(problem, std::vector<double>(problem.dof_count(), 0.0), options);
        const int expected_iterations = coupling == 0.0 ? 1 : 2;
        if (!result.converged || result.nonlinear_iterations != expected_iterations)
            std::cerr << "changing_coupling=" << coupling << " iterations=" << result.nonlinear_iterations
                      << " reason=" << result.convergence_reason << " state0=" << result.state[0]
                      << " state_last=" << result.state.back() << '\n';
        passed = check(result.converged && result.nonlinear_iterations == expected_iterations,
                     "direct factorization updates tiny couplings within Newton iteration and across solves")
                 && passed;
        for (std::size_t dof = 0; dof < result.state.size(); ++dof) {
            const double expected = 1.0;
            passed = check(std::abs(result.state[dof] / expected - 1.0) < 1.0e-12,
                         "direct factorization preserves the exact coupled solution")
                     && passed;
        }
    }
    return passed;
}

class LogDomainProblem final : public TwelveDofProblem {
  public:
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const override {
        validate_contribution(index);
        residual.resize(state.size());
        for (std::size_t dof = 0; dof < state.size(); ++dof) {
            if (!(state[dof] > 0.0))
                throw std::domain_error("log-domain Newton iterate must remain positive");
            residual[dof] = std::log(state[dof]) + 10.0;
        }
        if (jacobian == nullptr)
            return;
        jacobian->assign(state.size() * state.size(), 0.0);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            (*jacobian)[dof * state.size() + dof] = 1.0 / state[dof];
    }
};

class StagnatingProblem final : public TwelveDofProblem {
  public:
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const override {
        validate_contribution(index);
        residual.assign(state.size(), 1.0);
        if (jacobian == nullptr)
            return;
        jacobian->assign(state.size() * state.size(), 0.0);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            (*jacobian)[dof * state.size() + dof] = 1.0e20;
    }
};

class FieldStagnatingProblem final : public TwelveDofProblem {
  public:
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const override {
        validate_contribution(index);
        residual.assign(state.size(), 0.0);
        for (std::size_t dof = 0; dof < 4; ++dof)
            residual[dof] = 1.1;
        if (jacobian == nullptr)
            return;
        jacobian->assign(state.size() * state.size(), 0.0);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            (*jacobian)[dof * state.size() + dof] = 1.0e20;
    }
};

class QuadraticProblem final : public TwelveDofProblem {
  public:
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const override {
        validate_contribution(index);
        residual.resize(state.size());
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            residual[dof] = state[dof] * state[dof] - 2.0;
        if (jacobian == nullptr)
            return;
        jacobian->assign(state.size() * state.size(), 0.0);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            (*jacobian)[dof * state.size() + dof] = 2.0 * state[dof];
    }
};

class RuntimeLayoutProblem final : public fuelsim::NonlinearProblem {
  public:
    std::size_t dof_count() const noexcept override { return 32; }

    std::size_t contribution_count() const noexcept override { return 2; }

    const std::vector<fuelsim::FieldDescriptor>& field_layout() const noexcept override { return _fields; }

    std::size_t cax4_local_dof_count(std::size_t index) const {
        if (index == 0)
            return 32;
        if (index == 1)
            return _narrow_dofs.size();
        throw std::out_of_range("RuntimeLayoutProblem contribution index");
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override {
        ++_dof_mapping_calls;
        if (index == 0) {
            dofs.resize(dof_count());
            for (std::size_t dof = 0; dof < dofs.size(); ++dof)
                dofs[dof] = dof;
            return;
        }
        if (index == 1) {
            dofs = _narrow_dofs;
            return;
        }
        throw std::out_of_range("RuntimeLayoutProblem contribution index");
    }

    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const override {
        const std::size_t local_count = cax4_local_dof_count(index);
        if (state.size() != local_count)
            throw std::invalid_argument("RuntimeLayoutProblem contribution state size");
        if (jacobian == nullptr)
            ++_residual_calls.at(index);
        else
            ++_system_calls.at(index);
        compute_residual_values(index, state, residual);
        if (jacobian == nullptr)
            return;
        jacobian->resize(local_count * local_count);
        for (std::size_t row = 0; row < local_count; ++row)
            for (std::size_t column = 0; column < local_count; ++column)
                (*jacobian)[row * local_count + column] = contribution_coefficient(index, row, column);
    }

    const std::vector<fuelsim::DirichletCondition>& dirichlet_conditions() const noexcept override {
        return _conditions;
    }

    void reset_callback_counts() const noexcept {
        _residual_calls = {};
        _system_calls = {};
    }

    std::size_t residual_call_count(std::size_t contribution) const { return _residual_calls.at(contribution); }

    std::size_t system_call_count(std::size_t contribution) const { return _system_calls.at(contribution); }

    std::size_t dof_mapping_call_count() const noexcept { return _dof_mapping_calls; }

    double target_value(std::size_t dof) const {
        if (dof >= dof_count())
            throw std::out_of_range("RuntimeLayoutProblem target DOF");
        return 0.5 + 0.025 * static_cast<double>(dof);
    }

    double contribution_coefficient(std::size_t index, std::size_t row, std::size_t column) const {
        const std::size_t local_count = cax4_local_dof_count(index);
        if (row >= local_count || column >= local_count)
            throw std::out_of_range("RuntimeLayoutProblem coefficient index");
        if (index == 0) {
            if (row == column)
                return 4.0 + 0.01 * static_cast<double>(row + 1);
            if (column == (row + 1) % local_count)
                return 0.2;
            if (column == (row + 11) % local_count)
                return -0.04;
            return 0.0;
        }
        if (row == column)
            return 0.7 + 0.02 * static_cast<double>(row);
        return 0.005 * static_cast<double>((row + 1) * (column + 1));
    }

  private:
    void
    compute_residual_values(std::size_t index, const std::vector<double>& state, std::vector<double>& residual) const {
        const std::size_t local_count = cax4_local_dof_count(index);
        if (state.size() != local_count)
            throw std::invalid_argument("RuntimeLayoutProblem contribution state size");
        residual.assign(local_count, 0.0);
        for (std::size_t row = 0; row < local_count; ++row)
            for (std::size_t column = 0; column < local_count; ++column)
                residual[row] += contribution_coefficient(index, row, column)
                                 * (state[column] - target_value(index == 0 ? column : _narrow_dofs[column]));
    }

    std::vector<fuelsim::FieldDescriptor> _fields = {
        {"displacement_x", 0, 8, fuelsim::FieldCategory::mechanical},
        {"temperature", 8, 16, fuelsim::FieldCategory::thermal},
        {"displacement_y", 16, 24, fuelsim::FieldCategory::mechanical},
        {"displacement_z", 24, 32, fuelsim::FieldCategory::mechanical},
    };
    std::vector<std::size_t> _narrow_dofs = {0, 9, 14, 18, 23, 27, 31};
    std::vector<fuelsim::DirichletCondition> _conditions;
    mutable std::array<std::size_t, 2> _residual_calls{};
    mutable std::array<std::size_t, 2> _system_calls{};
    mutable std::size_t _dof_mapping_calls = 0;
};

bool test_runtime_contribution_layout() {
    RuntimeLayoutProblem problem;
    const std::vector<fuelsim::FieldDescriptor>& fields = problem.field_layout();
    bool passed = check(fields.size() == 4 && fields[0].name == "displacement_x"
                            && fields[0].category == fuelsim::FieldCategory::mechanical
                            && fields[1].name == "temperature" && fields[1].begin == 8 && fields[1].end == 16
                            && fields[1].category == fuelsim::FieldCategory::thermal
                            && fields[2].category == fuelsim::FieldCategory::mechanical
                            && fields[3].category == fuelsim::FieldCategory::mechanical,
        "runtime field metadata supports a non-leading thermal field and three mechanical fields");
    passed = check(problem.cax4_local_dof_count(0) == 32 && problem.cax4_local_dof_count(1) == 7,
                 "runtime contributions report a wide contribution followed by a narrow contribution")
             && passed;
    std::vector<double> state(problem.dof_count());
    std::vector<double> direction(problem.dof_count());
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        const int centered_index = static_cast<int>(dof % 5) - 2;
        state[dof] = problem.target_value(dof) + 0.03 * static_cast<double>(centered_index);
        direction[dof] = 0.2 + 0.015 * static_cast<double>(dof % 9);
    }
    fuelsim::ContributionWorkspace workspace;
    problem.evaluate_contribution(0, state, workspace, true);
    passed =
        check(workspace.dofs.size() == 32 && workspace.residual.size() == 32 && workspace.jacobian.size() == 32 * 32
                  && workspace.jacobian[3 * 32 + 4] == 0.2 && workspace.jacobian[5 * 32 + 16] == -0.04,
            "wide contribution exposes a 32 by 32 row-major Jacobian")
        && passed;
    problem.evaluate_contribution(1, state, workspace, true);
    passed = check(workspace.dofs.size() == 7 && workspace.residual.size() == 7 && workspace.jacobian.size() == 7 * 7
                       && std::abs(workspace.jacobian[2 * 7 + 5] - 0.09) < 1.0e-15,
                 "reused contribution workspace shrinks to a 7 by 7 row-major Jacobian")
             && passed;
    std::vector<double> expected_residual(problem.dof_count(), 0.0);
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        for (std::size_t row = 0; row < dofs.size(); ++row) {
            double value = 0.0;
            for (std::size_t column = 0; column < dofs.size(); ++column)
                value += problem.contribution_coefficient(contribution, row, column)
                         * (state[dofs[column]] - problem.target_value(dofs[column]));
            expected_residual[dofs[row]] += value;
        }
    }
    std::vector<double> assembled_residual;
    assembled_residual = fuelsim::test::assembled_residual(problem, state);
    double maximum_assembly_difference = 0.0;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof)
        maximum_assembly_difference =
            std::max(maximum_assembly_difference, std::abs(assembled_residual[dof] - expected_residual[dof]));
    passed = check(maximum_assembly_difference < 1.0e-14,
                 "runtime residual assembly handles consecutive 32 and 7 DOF contributions")
             && passed;
    const fuelsim::test::DirectionalJacobianCheck directional =
        fuelsim::test::check_directional_jacobian(problem, state, direction, 1.0e-6);
    passed = check(directional.difference.l2.size() == fields.size()
                       && directional.difference.maximum_absolute.size() == fields.size(),
                 "directional Jacobian diagnostics return one result per runtime field")
             && passed;
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const double reference = directional.finite_difference_directional_derivative.l2[field];
        passed = check(reference > 0.0 && directional.difference.l2[field] < 1.0e-8 * (1.0 + reference)
                           && directional.difference.maximum_absolute[field] < 1.0e-8 * (1.0 + reference),
                     "runtime row-major Jacobian matches the centered directional difference for " + fields[field].name)
                 && passed;
    }
    std::vector<double> initial(problem.dof_count(), 0.0);
    fuelsim::SolverOptions options;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::gmres;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::field_split;
    options.backtracking_fallback = false;
    options.linear_relative_tolerance = 1.0e-12;
    options.maximum_linear_iterations = 200;
    options.temperature_residual_scale = 200.0;
    options.mechanical_residual_scale = 20.0;
    std::vector<double> expected_initial_residual;
    expected_initial_residual = fuelsim::test::assembled_residual(problem, initial);
    std::vector<double> expected_initial_field_norms(fields.size(), 0.0);
    for (std::size_t field = 0; field < fields.size(); ++field)
        for (std::size_t dof = fields[field].begin; dof < fields[field].end; ++dof)
            expected_initial_field_norms[field] =
                std::hypot(expected_initial_field_norms[field], expected_initial_residual[dof]);
    fuelsim::PetscSolver solver;
    problem.reset_callback_counts();
    const fuelsim::SolveResult result = solver.solve(problem, initial, options);
    passed = check(result.converged && result.nonlinear_iterations == 1,
                 "runtime layout linear system converges in one exact Newton update")
             && passed;
    const bool field_result_sizes = result.field_names.size() == fields.size()
                                    && result.initial_field_residual_norms.size() == fields.size()
                                    && result.field_residual_reference_norms.size() == fields.size()
                                    && result.final_field_residual_norms.size() == fields.size()
                                    && result.final_scaled_field_residual_norms.size() == fields.size()
                                    && result.field_residual_scalings.size() == fields.size();
    passed = check(field_result_sizes, "PETSc solve reports all four runtime fields") && passed;
    if (field_result_sizes) {
        passed = check(result.field_names[0] == "displacement_x" && result.field_names[1] == "temperature"
                           && result.field_names[2] == "displacement_y" && result.field_names[3] == "displacement_z",
                     "PETSc solve preserves runtime field names and order")
                 && passed;
        passed = check(std::abs(result.field_residual_scalings[0] - 0.05) < 1.0e-15
                           && std::abs(result.field_residual_scalings[1] - 0.005) < 1.0e-15
                           && std::abs(result.field_residual_scalings[2] - 0.05) < 1.0e-15
                           && std::abs(result.field_residual_scalings[3] - 0.05) < 1.0e-15,
                     "fixed scaling follows thermal and mechanical categories instead of field index")
                 && passed;
        for (std::size_t field = 0; field < fields.size(); ++field) {
            passed = check(std::abs(result.initial_field_residual_norms[field] - expected_initial_field_norms[field])
                               < 1.0e-13 * (1.0 + expected_initial_field_norms[field]),
                         "PETSc initial residual includes every runtime contribution for " + fields[field].name)
                     && passed;
        }
    }
    if (result.state.size() == problem.dof_count()) {
        double maximum_solution_error = 0.0;
        for (std::size_t dof = 0; dof < problem.dof_count(); ++dof)
            maximum_solution_error =
                std::max(maximum_solution_error, std::abs(result.state[dof] - problem.target_value(dof)));
        passed =
            check(maximum_solution_error < 1.0e-10, "runtime layout solve reaches the known linear root") && passed;
    } else {
        passed = check(false, "runtime layout solve returns all 32 global DOFs") && passed;
    }
    const auto expected_partition = problem.contribution_partition(static_cast<std::size_t>(result.mpi_rank),
        static_cast<std::size_t>(result.mpi_size));
    const std::size_t expected_begin = expected_partition.first;
    const std::size_t expected_end = expected_partition.second;
    passed = check(result.local_contribution_begin == expected_begin && result.local_contribution_end == expected_end,
                 "runtime layout solve preserves the exact per-rank contribution interval")
             && passed;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        const bool owned = contribution >= expected_begin && contribution < expected_end;
        passed = check((problem.residual_call_count(contribution) > 0) == owned
                           && (problem.system_call_count(contribution) > 0) == owned,
                     "PETSc evaluates each runtime contribution only on its assigned rank")
                 && passed;
    }
    const std::size_t mapping_calls_after_setup = problem.dof_mapping_call_count();
    problem.reset_callback_counts();
    const fuelsim::SolveResult reused = solver.solve(problem, initial, options);
    passed = check(reused.converged && reused.timing.workspace_setups == 0
                       && problem.dof_mapping_call_count() > mapping_calls_after_setup,
                 "reused PETSc workspace refreshes dynamic contribution mappings without rebuilding")
             && passed;
    if (result.mpi_size == 2 && result.local_contribution_begin < result.local_contribution_end) {
        const std::size_t expected_local_width = result.mpi_rank == 0 ? 32 : 7;
        passed = check(problem.cax4_local_dof_count(result.local_contribution_begin) == expected_local_width,
                     "two-rank solve assigns the wide and narrow contributions to different ranks")
                 && passed;
    }
    fuelsim::SolverOptions automatic_options = options;
    automatic_options.field_residual_scaling = true;
    automatic_options.temperature_residual_scale = 0.0;
    automatic_options.mechanical_residual_scale = 0.0;
    fuelsim::PetscSolver automatic_solver;
    problem.reset_callback_counts();
    const fuelsim::SolveResult automatic = automatic_solver.solve(problem, initial, automatic_options);
    passed = check(automatic.converged && automatic.initial_field_residual_norms.size() == fields.size()
                       && automatic.field_residual_scalings.size() == fields.size(),
                 "automatic scaling solves the four-field runtime problem")
             && passed;
    if (automatic.initial_field_residual_norms.size() == fields.size()
        && automatic.field_residual_scalings.size() == fields.size()) {
        const double thermal_scale = 1.0 / expected_initial_field_norms[1];
        double mechanical_norm = 0.0;
        for (const std::size_t field : {0U, 2U, 3U})
            mechanical_norm = std::hypot(mechanical_norm, expected_initial_field_norms[field]);
        const double mechanical_scale = 1.0 / mechanical_norm;
        bool initial_norms_match = true;
        for (std::size_t field = 0; field < fields.size(); ++field) {
            initial_norms_match =
                initial_norms_match
                && std::abs(automatic.initial_field_residual_norms[field] - expected_initial_field_norms[field])
                       < 1.0e-13 * (1.0 + expected_initial_field_norms[field]);
        }
        passed =
            check(initial_norms_match
                      && std::abs(automatic.field_residual_scalings[1] - thermal_scale) < 1.0e-14 * thermal_scale
                      && std::abs(automatic.field_residual_scalings[0] - mechanical_scale) < 1.0e-14 * mechanical_scale
                      && std::abs(automatic.field_residual_scalings[2] - mechanical_scale) < 1.0e-14 * mechanical_scale
                      && std::abs(automatic.field_residual_scalings[3] - mechanical_scale) < 1.0e-14 * mechanical_scale,
                "automatic scaling uses independently assembled thermal and three-field mechanical norms")
            && passed;
    }
    RuntimeLayoutProblem changed_problem;
    passed = check(problem.discretization_identity() != changed_problem.discretization_identity(),
                 "distinct nonlinear problems have distinct PETSc workspace identities")
             && passed;
    fuelsim::SolverOptions guard_options = options;
    fuelsim::PetscSolver guard_solver;
    const fuelsim::SolveResult guard_baseline = guard_solver.solve(changed_problem, initial, guard_options);
    RuntimeLayoutProblem replacement_problem;
    const fuelsim::SolveResult replacement = guard_solver.solve(replacement_problem, initial, guard_options);
    passed = check(guard_baseline.converged && replacement.converged && replacement.timing.workspace_setups == 1,
                 "PETSc workspace rebuilds for a distinct discretization identity")
             && passed;
    return passed;
}

bool test_global_newton_safeguards() {
    LogDomainProblem domain_problem;
    fuelsim::PetscSolver failing_domain_solver;
    const std::vector<double> initial(fuelsim::cax4_local_dof_count, 1.0);
    fuelsim::SolverOptions failing_domain_options;
    failing_domain_options.backtracking_fallback = false;
    const fuelsim::SolveResult domain_failure =
        failing_domain_solver.solve(domain_problem, initial, failing_domain_options);
    bool passed = check(!domain_failure.converged
                            && domain_failure.failure_category == fuelsim::SolveFailureCategory::physical_domain
                            && !domain_failure.failure_message.empty(),
        "domain failure category and message are available on every rank");
    fuelsim::PetscSolver domain_solver;
    fuelsim::SolverOptions domain_options;
    const fuelsim::SolveResult domain = domain_solver.solve(domain_problem, initial, domain_options);
    if (!domain.converged)
        std::cerr << "log-domain failure: category=" << fuelsim::solve_failure_category_name(domain.failure_category)
                  << " reason=" << domain.convergence_reason << " residual=" << domain.residual_norm
                  << " message=" << domain.failure_message << '\n';
    passed = check(domain.converged && domain.used_backtracking_fallback && domain.nonlinear_attempts == 2
                       && domain.linear_iterations > 0
                       && domain.initial_failure_category == fuelsim::SolveFailureCategory::physical_domain,
                 "BASIC failure automatically retries with backtracking "
                 "from the original state and reports KSP work")
             && passed;
    const double target = std::exp(-10.0);
    fuelsim::PetscSolver critical_point_solver;
    auto critical_point_options = domain_options;
    critical_point_options.line_search = fuelsim::SolverOptions::LineSearch::critical_point;
    const auto critical_point = critical_point_solver.solve(domain_problem, initial, critical_point_options);
    passed = check(critical_point.converged && critical_point.used_backtracking_fallback
                       && critical_point.nonlinear_attempts == 2
                       && critical_point.initial_failure_category == fuelsim::SolveFailureCategory::physical_domain,
                 "Critical-point domain failure retries from the original state using backtracking")
             && passed;
    for (const double value : critical_point.state)
        passed = check(std::abs(value - target) < 1.0e-10 * target,
                     "Critical-point fallback reaches the positive logarithmic root")
                 && passed;
    for (double value : domain.state)
        passed =
            check(std::abs(value - target) < 1.0e-10 * target, "backtracking reaches the positive logarithmic root")
            && passed;
    StagnatingProblem stagnating_problem;
    fuelsim::PetscSolver stagnating_solver;
    fuelsim::SolverOptions options;
    options.line_search = fuelsim::SolverOptions::LineSearch::basic;
    options.backtracking_fallback = false;
    options.field_residual_scaling = false;
    options.step_tolerance = 1.0e-8;
    const fuelsim::SolveResult stagnating = stagnating_solver.solve(stagnating_problem, initial, options);
    passed = check(stagnating.convergence_reason > 0 && !stagnating.converged
                       && stagnating.failure_category == fuelsim::SolveFailureCategory::residual_verification,
                 "positive step-stagnation reason fails residual review")
             && passed;
    FieldStagnatingProblem field_problem;
    fuelsim::PetscSolver field_solver;
    fuelsim::SolverOptions field_options = options;
    field_options.temperature_residual_absolute_tolerance = 2.0;
    field_options.mechanical_residual_absolute_tolerance = 2.0;
    const fuelsim::SolveResult field_failure = field_solver.solve(field_problem, initial, field_options);
    passed = check(!field_failure.converged && field_failure.residual_norm < std::sqrt(12.0)
                       && field_failure.final_scaled_field_residual_norms[0]
                              > field_options.temperature_residual_absolute_tolerance
                       && field_failure.failure_category == fuelsim::SolveFailureCategory::residual_verification
                       && field_failure.failure_message.find("field0=") != std::string::npos,
                 "field residual audit rejects a state that passes the "
                 "combined physical absolute scale")
             && passed;
    QuadraticProblem quadratic_problem;
    fuelsim::PetscSolver quadratic_solver;
    fuelsim::SolverOptions quadratic_options;
    quadratic_options.maximum_iterations = 1;
    quadratic_options.backtracking_fallback = false;
    quadratic_options.field_residual_scaling = false;
    quadratic_options.absolute_tolerance = 1.0e-14;
    quadratic_options.relative_tolerance = 1.0e-14;
    quadratic_options.residual_reduction_tolerance = 0.3;
    const fuelsim::SolveResult rescued = quadratic_solver.solve(quadratic_problem, initial, quadratic_options);
    passed = check(rescued.convergence_reason < 0 && rescued.converged
                       && rescued.failure_category == fuelsim::SolveFailureCategory::none
                       && rescued.failure_message.empty() && rescued.residual_norm < 0.3 * std::sqrt(12.0),
                 "audited residual reduction rescues a PETSc maximum-"
                 "iteration reason at an acceptable state")
             && passed;
    fuelsim::PetscSolver field_convergence_solver;
    fuelsim::SolverOptions field_convergence_options;
    field_convergence_options.field_residual_scaling = false;
    field_convergence_options.field_residual_convergence = true;
    field_convergence_options.absolute_tolerance = 1.0e-14;
    field_convergence_options.relative_tolerance = 1.0e-14;
    field_convergence_options.residual_reduction_tolerance = 1.0e-14;
    field_convergence_options.temperature_residual_absolute_tolerance = 0.6;
    field_convergence_options.mechanical_residual_absolute_tolerance = 0.6;
    const fuelsim::SolveResult field_converged =
        field_convergence_solver.solve(quadratic_problem, initial, field_convergence_options);
    passed = check(field_converged.converged && field_converged.nonlinear_iterations == 1
                       && field_converged.residual_norm > field_convergence_options.absolute_tolerance,
                 "field residual convergence stops when every physical field meets its configured tolerance")
             && passed;
    return passed;
}

bool test_linear_load_predictor() {
    const auto mesh = fuelsim::test::make_disconnected_annular_mesh({{1, "solid", 0.0, 0.004, 0.01, 3, 2}});
    const auto definition =
        single_region_definition(constant_material(4.0, 1.0e-5), 2.0e8, 600.0, 600.0, 0.0, 0.0, 0.0);
    fuelsim::SteadyProblem baseline_problem(definition, mesh), predicted_problem(definition, mesh);
    fuelsim::SteadyLoadOptions loading;
    loading.load_steps = 5;
    fuelsim::SolverOptions options;
    options.line_search = fuelsim::SolverOptions::LineSearch::basic;
    const auto baseline = fuelsim::solve_steady(baseline_problem, loading, options);
    loading.use_linear_load_predictor = true;
    const auto predicted = fuelsim::solve_steady(predicted_problem, loading, options);
    bool fields_equal = true;
    for (const auto& field : baseline_problem.field_layout()) {
        double difference = 0.0, scale = 0.0;
        for (std::size_t i = field.begin; i < field.end; ++i) {
            difference = std::hypot(difference, predicted.solve.state[i] - baseline.solve.state[i]);
            scale = std::hypot(scale, baseline.solve.state[i]);
        }
        const double absolute = field.category == fuelsim::FieldCategory::thermal ? 1.0e-9 : 1.0e-12;
        fields_equal = fields_equal && difference < absolute + 1.0e-10 * scale;
    }
    return check(baseline.completed && predicted.completed && predicted.load_predictor_attempts == 3
                     && predicted.load_predictor_fallbacks == 0 && predicted.total_cutbacks == 0
                     && predicted.total_nonlinear_iterations < baseline.total_nonlinear_iterations && fields_equal,
        "linear load prediction uses two accepted equilibria and reduces iterations without changing the solution");
}

bool test_thermal_cylinder() {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double conductivity = 4.0;
    constexpr double heat_source = 2.0e8;
    constexpr double outer_temperature = 600.0;
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::test::make_disconnected_annular_mesh({{1, "solid", 0.0, radius, length, 32, 2}});
    fuelsim::SteadyProblem problem(single_region_definition(constant_material(conductivity, 0.0),
                                       heat_source,
                                       outer_temperature,
                                       outer_temperature,
                                       0.0,
                                       0.0,
                                       0.0),
        mesh);
    fuelsim::PetscSolver solver;
    fuelsim::SolverOptions scaled_options;
    scaled_options.field_residual_scaling = true;
    const fuelsim::SolveResult result = solver.solve(problem, problem.initial_state(), scaled_options);
    bool passed = check(result.converged, "thermal cylinder SNES converged");
    const auto expected_partition = problem.contribution_partition(static_cast<std::size_t>(result.mpi_rank),
        static_cast<std::size_t>(result.mpi_size));
    const std::size_t expected_begin = expected_partition.first;
    const std::size_t expected_end = expected_partition.second;
    passed = check(result.local_contribution_begin == expected_begin && result.local_contribution_end == expected_end,
                 "PETSc rank owns its exact nonoverlapping contribution range")
             && passed;
    if (result.mpi_size > 1)
        passed = check(result.local_contribution_end - result.local_contribution_begin < problem.contribution_count(),
                     "MPI rank does not repeat the complete model assembly")
                 && passed;
    double maximum_scaled_error = 0.0;
    const double center_rise = heat_source * radius * radius / (4.0 * conductivity);
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node) {
        const double r = fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes()[node].r;
        const double expected = outer_temperature + heat_source * (radius * radius - r * r) / (4.0 * conductivity);
        const double actual =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, node)];
        maximum_scaled_error = std::max(maximum_scaled_error, std::abs(actual - expected) / center_rise);
    }
    passed = check(maximum_scaled_error < 1.0e-3,
                 "thermal cylinder temperature error is below 0.1%; actual=" + std::to_string(maximum_scaled_error))
             && passed;
    std::cout << "thermal_cylinder_maximum_scaled_error=" << maximum_scaled_error << '\n';
    return passed;
}

double thermal_cylinder_error(std::size_t radial_elements) {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double conductivity = 4.0;
    constexpr double heat_source = 2.0e8;
    constexpr double outer_temperature = 600.0;
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::test::make_disconnected_annular_mesh({{1, "solid", 0.0, radius, length, radial_elements, 2}});
    fuelsim::SteadyProblem problem(single_region_definition(constant_material(conductivity, 0.0),
                                       heat_source,
                                       outer_temperature,
                                       outer_temperature,
                                       0.0,
                                       0.0,
                                       0.0),
        mesh);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result = solver.solve(problem, problem.initial_state());
    if (!result.converged)
        throw std::runtime_error("thermal mesh-convergence solve did not converge");
    const double center_rise = heat_source * radius * radius / (4.0 * conductivity);
    double maximum_error = 0.0;
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node) {
        const double r = fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes()[node].r;
        const double expected = outer_temperature + heat_source * (radius * radius - r * r) / (4.0 * conductivity);
        const double actual =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, node)];
        maximum_error = std::max(maximum_error, std::abs(actual - expected) / center_rise);
    }
    return maximum_error;
}

bool test_thermal_mesh_convergence() {
    const std::array<double, 3> errors = {thermal_cylinder_error(8),
        thermal_cylinder_error(16),
        thermal_cylinder_error(32)};
    const double first_ratio = errors[0] / errors[1];
    const double second_ratio = errors[1] / errors[2];
    const double first_order = std::log2(first_ratio);
    const double second_order = std::log2(second_ratio);
    std::cout << "thermal_mesh_convergence_errors=" << errors[0] << ',' << errors[1] << ',' << errors[2] << '\n';
    std::cout << "thermal_mesh_convergence_ratios=" << first_ratio << ',' << second_ratio << '\n';
    std::cout << "thermal_mesh_convergence_orders=" << first_order << ',' << second_order << '\n';
    return check(first_order > 1.7 && second_order > 1.7 && second_order > first_order,
        "successive radial mesh refinement approaches second-order "
        "thermal convergence");
}

bool test_free_thermal_expansion() {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double alpha = 1.0e-5;
    constexpr double temperature = 700.0;
    constexpr double temperature_change = 100.0;
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::test::make_disconnected_annular_mesh({{1, "solid", 0.0, radius, length, 8, 4}});
    const fuelsim::SteadyProblem problem(
        single_region_definition(constant_material(4.0, alpha), 0.0, temperature, temperature, 0.0, 0.0, 0.0),
        mesh);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result = solver.solve(problem, problem.initial_state());
    bool passed = check(result.converged, "free expansion SNES converged");
    double maximum_temperature_error = 0.0;
    double maximum_displacement_error = 0.0;
    const double displacement_scale = alpha * temperature_change * length;
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node) {
        const fuelsim::RzPoint& point = fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes()[node];
        const double actual_temperature =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, node)];
        const double actual_radial =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, node)];
        const double actual_axial =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, node)];
        maximum_temperature_error = std::max(maximum_temperature_error, std::abs(actual_temperature - temperature));
        maximum_displacement_error =
            std::max(maximum_displacement_error, std::abs(actual_radial - alpha * temperature_change * point.r));
        maximum_displacement_error =
            std::max(maximum_displacement_error, std::abs(actual_axial - alpha * temperature_change * point.z));
    }
    double maximum_stress = 0.0;
    for (std::size_t element = 0; element < fuelsim::rz::ProblemAccess::region_element_count(problem, 0); ++element) {
        const fuelsim::Cax4LocalValues state =
            fuelsim::rz::ProblemAccess::contribution_state(problem, element, result.state);
        const auto stresses =
            fuelsim::compute_cax4_thermoelastic_stress(fuelsim::rz::ProblemAccess::region_kernel_data(problem, 0),
                fuelsim::rz::ProblemAccess::region_element_geometry(problem, 0, element),
                state);
        for (const fuelsim::AxisymmetricStressValues& stress : stresses) {
            maximum_stress = std::max(maximum_stress, std::abs(stress.rr));
            maximum_stress = std::max(maximum_stress, std::abs(stress.zz));
            maximum_stress = std::max(maximum_stress, std::abs(stress.hoop));
            maximum_stress = std::max(maximum_stress, std::abs(stress.rz));
        }
    }
    passed = check(maximum_temperature_error < 1.0e-9, "free expansion temperature is uniform") && passed;
    passed = check(maximum_displacement_error / displacement_scale < 1.0e-8,
                 "free expansion displacement matches analytic field")
             && passed;
    passed = check(maximum_stress < 100.0, "free expansion stress is below 100 Pa") && passed;
    std::cout << "free_expansion_maximum_temperature_error=" << maximum_temperature_error << '\n';
    std::cout << "free_expansion_maximum_displacement_relative_error="
              << maximum_displacement_error / displacement_scale << '\n';
    std::cout << "free_expansion_maximum_stress=" << maximum_stress << '\n';
    return passed;
}

bool test_lame_open_ended_cylinder() {
    constexpr double inner_radius = 0.004;
    constexpr double outer_radius = 0.005;
    constexpr double length = 0.01;
    constexpr double pressure = 1.0e6;
    constexpr double young_modulus = 75.0e9;
    constexpr double poisson_ratio = 0.3;
    const fuelsim::ThermoelasticProperties material =
        fuelsim::test::thermoelastic(0.0, 4.0, young_modulus, poisson_ratio, 0.0, 600.0);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::test::make_disconnected_annular_mesh({{1, "solid", inner_radius, outer_radius, length, 48, 2}});
    const fuelsim::SteadyProblem problem(
        single_region_definition(material, 0.0, 600.0, 600.0, inner_radius, pressure, 0.0),
        mesh);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result = solver.solve(problem, problem.initial_state());
    bool passed = check(result.converged, "Lame cylinder SNES converged");
    const double A =
        pressure * inner_radius * inner_radius / (outer_radius * outer_radius - inner_radius * inner_radius);
    std::cout << "lame_initial_field_residuals=" << result.initial_field_residual_norms[0] << ','
              << result.initial_field_residual_norms[1] << ',' << result.initial_field_residual_norms[2] << '\n';
    std::cout << "lame_final_scaled_field_residuals=" << result.final_scaled_field_residual_norms[0] << ','
              << result.final_scaled_field_residual_norms[1] << ',' << result.final_scaled_field_residual_norms[2]
              << '\n';
    const double B = pressure * inner_radius * inner_radius * outer_radius * outer_radius
                     / (outer_radius * outer_radius - inner_radius * inner_radius);
    const double axial_strain = -2.0 * poisson_ratio * A / young_modulus;
    double maximum_radial_relative_error = 0.0;
    double maximum_axial_relative_error = 0.0;
    const double radial_scale =
        ((1.0 - poisson_ratio) * A * inner_radius + (1.0 + poisson_ratio) * B / inner_radius) / young_modulus;
    const double axial_scale = std::abs(axial_strain * length);
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node) {
        const fuelsim::RzPoint& point = fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes()[node];
        const double expected_radial =
            ((1.0 - poisson_ratio) * A * point.r + (1.0 + poisson_ratio) * B / point.r) / young_modulus;
        const double expected_axial = axial_strain * point.z;
        const double actual_radial =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, node)];
        const double actual_axial =
            result.state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, node)];
        maximum_radial_relative_error =
            std::max(maximum_radial_relative_error, std::abs(actual_radial - expected_radial) / radial_scale);
        if (axial_scale > 0.0) {
            maximum_axial_relative_error =
                std::max(maximum_axial_relative_error, std::abs(actual_axial - expected_axial) / axial_scale);
        }
    }
    passed = check(maximum_radial_relative_error < 3.0e-3, "Lame radial displacement error is below 0.3%") && passed;
    passed = check(maximum_axial_relative_error < 3.0e-3, "Lame axial displacement error is below 0.3%") && passed;
    std::cout << "lame_maximum_radial_relative_error=" << maximum_radial_relative_error << '\n';
    std::cout << "lame_maximum_axial_relative_error=" << maximum_axial_relative_error << '\n';
    return passed;
}

bool test_m1_open_gap_analytic_thermal() {
    constexpr double fuel_radius = 0.004;
    constexpr double cladding_inner_radius = 0.0041;
    constexpr double cladding_outer_radius = 0.0046;
    constexpr double length = 0.010;
    constexpr double fuel_conductivity = 4.0;
    constexpr double cladding_conductivity = 16.0;
    constexpr double gap_conductivity = 0.4;
    constexpr double heat_source = 1.0e8;
    constexpr double outer_temperature = 600.0;
    constexpr std::size_t fuel_radial_elements = 32;
    constexpr std::size_t cladding_radial_elements = 8;
    constexpr std::size_t axial_elements = 2;
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::test::make_disconnected_annular_mesh({{1,
                                                                                                   "fuel",
                                                                                                   0.0,
                                                                                                   fuel_radius,
                                                                                                   length,
                                                                                                   fuel_radial_elements,
                                                                                                   axial_elements},
        {2, "clad", cladding_inner_radius, cladding_outer_radius, length, cladding_radial_elements, axial_elements}});
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back(
        {"fuel", "fuel", constant_material(fuel_conductivity, 0.0), heat_source, outer_temperature});
    definition.regions.push_back(
        {"clad", "clad", constant_material(cladding_conductivity, 0.0), 0.0, outer_temperature});
    definition.contacts.push_back(
        {"fuel_clad", "clad_inner", "fuel_outer", true, true, gap_conductivity, 1.0e-6, 1.0e14});
    definition.boundary_conditions.push_back(
        dirichlet("fuel_axis", "fuel_inner", fuelsim::Field::radial_displacement, 0.0));
    definition.boundary_conditions.push_back(
        dirichlet("fuel_bottom", "fuel_bottom", fuelsim::Field::axial_displacement, 0.0));
    definition.boundary_conditions.push_back(
        dirichlet("clad_bottom", "clad_bottom", fuelsim::Field::axial_displacement, 0.0));
    definition.boundary_conditions.push_back(
        dirichlet("clad_temperature", "clad_outer", fuelsim::Field::temperature, outer_temperature));
    const fuelsim::SteadyProblem problem(std::move(definition), mesh);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result = solver.solve(problem, problem.initial_state());
    const double cladding_rise = heat_source * fuel_radius * fuel_radius / (2.0 * cladding_conductivity)
                                 * std::log(cladding_outer_radius / cladding_inner_radius);
    const double gap_rise =
        heat_source * fuel_radius * (cladding_inner_radius - fuel_radius) / (2.0 * gap_conductivity);
    const double fuel_rise = heat_source * fuel_radius * fuel_radius / (4.0 * fuel_conductivity);
    const double expected_cladding_inner = outer_temperature + cladding_rise;
    const double expected_fuel_surface = expected_cladding_inner + gap_rise;
    const double expected_center = expected_fuel_surface + fuel_rise;
    const std::size_t axial_mid = axial_elements / 2;
    const std::size_t fuel_center_local = fuelsim::test::annular_node_id(fuel_radial_elements, 0, axial_mid);
    const std::size_t fuel_surface_local =
        fuelsim::test::annular_node_id(fuel_radial_elements, fuel_radial_elements, axial_mid);
    const std::size_t cladding_inner_local = fuelsim::test::annular_node_id(cladding_radial_elements, 0, axial_mid);
    const auto& dofs = fuelsim::rz::ProblemAccess::dof_map(problem);
    const double actual_center = result.state[dofs.dof(fuelsim::Field::temperature,
        fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) + fuel_center_local)];
    const double actual_fuel_surface = result.state[dofs.dof(fuelsim::Field::temperature,
        fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) + fuel_surface_local)];
    const double actual_cladding_inner = result.state[dofs.dof(fuelsim::Field::temperature,
        fuelsim::rz::ProblemAccess::region_node_offset(problem, 1) + cladding_inner_local)];
    const double temperature_scale = expected_center - outer_temperature;
    const double center_error = std::abs(actual_center - expected_center) / temperature_scale;
    const double fuel_surface_error = std::abs(actual_fuel_surface - expected_fuel_surface) / temperature_scale;
    const double cladding_inner_error = std::abs(actual_cladding_inner - expected_cladding_inner) / temperature_scale;
    const fuelsim::InterfaceSummary interface =
        fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, result.state);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double expected_heat_rate = heat_source * pi * fuel_radius * fuel_radius * length;
    const double heat_balance_error = metric_relative_error(interface.total_heat_rate, expected_heat_rate);
    bool passed = check(result.converged, "M1 open-gap thermal SNES converged");
    passed = check(center_error < 1.0e-3, "M1 analytic center temperature error is below 0.1%") && passed;
    passed = check(fuel_surface_error < 1.0e-3,
                 "M1 analytic fuel-surface temperature error is below "
                 "0.1%")
             && passed;
    passed = check(cladding_inner_error < 1.0e-3,
                 "M1 analytic cladding-inner temperature error is below "
                 "0.1%")
             && passed;
    passed = check(interface.minimum_gap > 0.0 && interface.maximum_contact_pressure == 0.0,
                 "M1 analytic thermal case remains out of contact")
             && passed;
    passed = check(heat_balance_error < 1.0e-9, "M1 interface heat rate balances generated power") && passed;
    std::cout << "m1_analytic_center_temperature_scaled_error=" << center_error << '\n';
    std::cout << "m1_analytic_fuel_surface_scaled_error=" << fuel_surface_error << '\n';
    std::cout << "m1_analytic_cladding_inner_scaled_error=" << cladding_inner_error << '\n';
    std::cout << "m1_analytic_heat_balance_relative_error=" << heat_balance_error << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    try {
        std::string selected_case = "all";
        if (argc >= 3 && std::string(argv[1]) == "--case") {
            selected_case = argv[2];
            for (int index = 3; index < argc; ++index)
                argv[index - 2] = argv[index];
            argc -= 2;
            argv[argc] = nullptr;
        }
        if (selected_case != "all" && selected_case != "runtime-layout"
            && selected_case != "thermal-mesh-convergence") {
            std::cerr << "Usage: fuelsim_solver_tests "
                         "[--case runtime-layout|thermal-mesh-convergence] [PETSc options]\n";
            return 2;
        }
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M0 and M1 numerical acceptance tests\n");
        bool passed = true;
        passed = test_runtime_contribution_layout() && passed;
        passed = test_changing_direct_coupling() && passed;
        if (selected_case == "runtime-layout") {
            if (!passed)
                return 1;
            std::cout << "[PASS] fuelsim runtime-layout solver tests\n";
            return 0;
        }
        if (selected_case == "thermal-mesh-convergence") {
            passed = test_thermal_mesh_convergence() && passed;
            if (!passed)
                return 1;
            std::cout << "[PASS] fuelsim thermal mesh-convergence study\n";
            return 0;
        }
        passed = test_global_newton_safeguards() && passed;
        passed = test_linear_load_predictor() && passed;
        passed = test_thermal_cylinder() && passed;
        passed = test_free_thermal_expansion() && passed;
        passed = test_lame_open_ended_cylinder() && passed;
        passed = test_m1_open_gap_analytic_thermal() && passed;
        if (!passed)
            return 1;
        std::cout << "[PASS] fuelsim M0 and M1 solver acceptance tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] solver tests raised: " << error.what() << '\n';
        return 1;
    }
}
