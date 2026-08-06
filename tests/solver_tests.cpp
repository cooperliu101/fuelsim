#include "fuelsim/petsc_solver.hpp"
#include "support/steady_fuel_cladding_problem.hpp"
#include "support/steady_single_region_problem.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
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

bool test_shadow_state_view() {
    const std::vector<std::uint32_t> dofs = {1U, 4U, 7U};
    const std::vector<double> values = {2.0, 5.0, 8.0};
    const fuelsim::GlobalStateView state(9, dofs, values);
    bool missing_rejected = false;
    try {
        (void)state.value(3);
    } catch (const std::out_of_range&) {
        missing_rejected = true;
    }
    return check(state.global_size() == 9 && state.local_size() == 3 &&
                     state.contains(1) && !state.contains(3) &&
                     state.value(4) == 5.0 && missing_rejected,
                 "shadow state exposes only declared global DOFs");
}

fuelsim::ThermoelasticProperties constant_material(double conductivity,
                                                   double thermal_expansion) {
    return {
        0.0, conductivity, 75.0e9, 0.3, thermal_expansion, 600.0,
    };
}

class LogDomainProblem final : public fuelsim::NonlinearProblem {
  public:
    std::size_t dof_count() const noexcept override {
        return fuelsim::local_dof_count;
    }
    std::size_t contribution_count() const noexcept override {
        return 1;
    }
    fuelsim::LocalDofs contribution_dofs(std::size_t index) const override {
        if (index != 0)
            throw std::out_of_range("LogDomainProblem contribution index");
        fuelsim::LocalDofs dofs{};
        for (std::size_t dof = 0; dof < dofs.size(); ++dof)
            dofs[dof] = dof;
        return dofs;
    }
    fuelsim::LocalResidual
    contribution_residual(std::size_t,
                          const fuelsim::LocalValues& state) const override {
        fuelsim::LocalResidual residual{};
        for (std::size_t dof = 0; dof < state.size(); ++dof) {
            if (!(state[dof] > 0.0))
                throw std::domain_error(
                    "log-domain Newton iterate must remain positive");
            residual[dof] = std::log(state[dof]) + 10.0;
        }
        return residual;
    }
    fuelsim::LocalSystem
    linearize_contribution(std::size_t index,
                           const fuelsim::LocalValues& state) const override {
        fuelsim::LocalSystem system{};
        system.residual = contribution_residual(index, state);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            system.jacobian[dof * state.size() + dof] = 1.0 / state[dof];
        return system;
    }
    const std::vector<fuelsim::DirichletCondition>&
    dirichlet_conditions() const noexcept override {
        return _conditions;
    }

  protected:
    void add_state_independent_residual(
        std::vector<double>&) const override {}

  private:
    std::vector<fuelsim::DirichletCondition> _conditions;
};

class StagnatingProblem final : public fuelsim::NonlinearProblem {
  public:
    std::size_t dof_count() const noexcept override {
        return fuelsim::local_dof_count;
    }
    std::size_t contribution_count() const noexcept override {
        return 1;
    }
    fuelsim::LocalDofs contribution_dofs(std::size_t) const override {
        fuelsim::LocalDofs dofs{};
        for (std::size_t dof = 0; dof < dofs.size(); ++dof)
            dofs[dof] = dof;
        return dofs;
    }
    fuelsim::LocalResidual
    contribution_residual(std::size_t,
                          const fuelsim::LocalValues&) const override {
        fuelsim::LocalResidual residual{};
        residual.fill(1.0);
        return residual;
    }
    fuelsim::LocalSystem
    linearize_contribution(std::size_t index,
                           const fuelsim::LocalValues& state) const override {
        fuelsim::LocalSystem system{};
        system.residual = contribution_residual(index, state);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            system.jacobian[dof * state.size() + dof] = 1.0e20;
        return system;
    }
    const std::vector<fuelsim::DirichletCondition>&
    dirichlet_conditions() const noexcept override {
        return _conditions;
    }

  protected:
    void add_state_independent_residual(
        std::vector<double>&) const override {}

  private:
    std::vector<fuelsim::DirichletCondition> _conditions;
};

class FieldStagnatingProblem final : public fuelsim::NonlinearProblem {
  public:
    std::size_t dof_count() const noexcept override {
        return fuelsim::local_dof_count;
    }
    std::size_t contribution_count() const noexcept override {
        return 1;
    }
    fuelsim::LocalDofs contribution_dofs(std::size_t) const override {
        fuelsim::LocalDofs dofs{};
        for (std::size_t dof = 0; dof < dofs.size(); ++dof)
            dofs[dof] = dof;
        return dofs;
    }
    fuelsim::LocalResidual
    contribution_residual(std::size_t,
                          const fuelsim::LocalValues&) const override {
        fuelsim::LocalResidual residual{};
        for (std::size_t dof = 0; dof < 4; ++dof)
            residual[dof] = 1.1;
        return residual;
    }
    fuelsim::LocalSystem
    linearize_contribution(std::size_t index,
                           const fuelsim::LocalValues& state) const override {
        fuelsim::LocalSystem system{};
        system.residual = contribution_residual(index, state);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            system.jacobian[dof * state.size() + dof] = 1.0e20;
        return system;
    }
    const std::vector<fuelsim::DirichletCondition>&
    dirichlet_conditions() const noexcept override {
        return _conditions;
    }

  protected:
    void add_state_independent_residual(
        std::vector<double>&) const override {}

  private:
    std::vector<fuelsim::DirichletCondition> _conditions;
};

class QuadraticProblem final : public fuelsim::NonlinearProblem {
  public:
    std::size_t dof_count() const noexcept override {
        return fuelsim::local_dof_count;
    }
    std::size_t contribution_count() const noexcept override {
        return 1;
    }
    fuelsim::LocalDofs contribution_dofs(std::size_t) const override {
        fuelsim::LocalDofs dofs{};
        for (std::size_t dof = 0; dof < dofs.size(); ++dof)
            dofs[dof] = dof;
        return dofs;
    }
    fuelsim::LocalResidual
    contribution_residual(std::size_t,
                          const fuelsim::LocalValues& state) const override {
        fuelsim::LocalResidual residual{};
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            residual[dof] = state[dof] * state[dof] - 2.0;
        return residual;
    }
    fuelsim::LocalSystem
    linearize_contribution(std::size_t index,
                           const fuelsim::LocalValues& state) const override {
        fuelsim::LocalSystem system{};
        system.residual = contribution_residual(index, state);
        for (std::size_t dof = 0; dof < state.size(); ++dof)
            system.jacobian[dof * state.size() + dof] = 2.0 * state[dof];
        return system;
    }
    const std::vector<fuelsim::DirichletCondition>&
    dirichlet_conditions() const noexcept override {
        return _conditions;
    }

  protected:
    void add_state_independent_residual(
        std::vector<double>&) const override {}

  private:
    std::vector<fuelsim::DirichletCondition> _conditions;
};

bool test_global_newton_safeguards() {
    LogDomainProblem domain_problem;
    fuelsim::PetscSolver failing_domain_solver;
    const std::vector<double> initial(fuelsim::local_dof_count, 1.0);
    fuelsim::SolverOptions failing_domain_options;
    failing_domain_options.backtracking_fallback = false;
    const fuelsim::SolveResult domain_failure = failing_domain_solver.solve(
        domain_problem, initial, failing_domain_options);
    bool passed = check(
        !domain_failure.converged &&
            domain_failure.failure_category ==
                fuelsim::SolveFailureCategory::physical_domain &&
            !domain_failure.failure_message.empty(),
        "domain failure category and message are available on every rank");

    fuelsim::PetscSolver domain_solver;
    fuelsim::SolverOptions domain_options;
    const fuelsim::SolveResult domain = domain_solver.solve(
        domain_problem, initial, domain_options);
    if (!domain.converged)
        std::cerr << "log-domain failure: category="
                  << fuelsim::solve_failure_category_name(
                         domain.failure_category)
                  << " reason=" << domain.convergence_reason
                  << " residual=" << domain.residual_norm
                  << " message=" << domain.failure_message << '\n';
    passed = check(domain.converged &&
                       domain.used_backtracking_fallback &&
                       domain.nonlinear_attempts == 2 &&
                       domain.linear_iterations > 0 &&
                       domain.basic_failure_category ==
                           fuelsim::SolveFailureCategory::physical_domain,
                   "BASIC failure automatically retries with backtracking "
                   "from the original state and reports KSP work") &&
             passed;
    const double target = std::exp(-10.0);
    for (double value : domain.state)
        passed = check(std::abs(value - target) < 1.0e-10 * target,
                       "backtracking reaches the positive logarithmic root") &&
                 passed;

    StagnatingProblem stagnating_problem;
    fuelsim::PetscSolver stagnating_solver;
    fuelsim::SolverOptions options;
    options.line_search = fuelsim::SolverOptions::LineSearch::basic;
    options.backtracking_fallback = false;
    options.field_residual_scaling = false;
    options.step_tolerance = 1.0e-8;
    const fuelsim::SolveResult stagnating =
        stagnating_solver.solve(stagnating_problem, initial, options);
    passed = check(stagnating.convergence_reason > 0 &&
                       !stagnating.converged &&
                       stagnating.failure_category ==
                           fuelsim::SolveFailureCategory::residual_verification,
                   "positive step-stagnation reason fails residual review") &&
             passed;

    FieldStagnatingProblem field_problem;
    fuelsim::PetscSolver field_solver;
    fuelsim::SolverOptions field_options = options;
    field_options.temperature_residual_absolute_tolerance = 2.0;
    field_options.mechanical_residual_absolute_tolerance = 2.0;
    const fuelsim::SolveResult field_failure =
        field_solver.solve(field_problem, initial, field_options);
    passed = check(
                 !field_failure.converged &&
                     field_failure.residual_norm < std::sqrt(12.0) &&
                     field_failure.final_scaled_field_residual_norms[0] >
                         field_options.temperature_residual_absolute_tolerance &&
                     field_failure.failure_category ==
                         fuelsim::SolveFailureCategory::residual_verification &&
                     field_failure.failure_message.find("field0=") !=
                         std::string::npos,
                 "field residual audit rejects a state that passes the "
                 "combined physical absolute scale") &&
             passed;

    QuadraticProblem quadratic_problem;
    fuelsim::PetscSolver quadratic_solver;
    fuelsim::SolverOptions quadratic_options;
    quadratic_options.maximum_iterations = 1;
    quadratic_options.backtracking_fallback = false;
    quadratic_options.field_residual_scaling = false;
    quadratic_options.absolute_tolerance = 1.0e-14;
    quadratic_options.relative_tolerance = 1.0e-14;
    quadratic_options.residual_reduction_tolerance = 0.3;
    const fuelsim::SolveResult rescued = quadratic_solver.solve(
        quadratic_problem, initial, quadratic_options);
    passed = check(rescued.convergence_reason < 0 && rescued.converged &&
                       rescued.failure_category ==
                           fuelsim::SolveFailureCategory::none &&
                       rescued.failure_message.empty() &&
                       rescued.residual_norm <
                           0.3 * std::sqrt(12.0),
                   "audited residual reduction rescues a PETSc maximum-"
                   "iteration reason at an acceptable state") &&
             passed;
    return passed;
}

bool test_thermal_cylinder() {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double conductivity = 4.0;
    constexpr double heat_source = 2.0e8;
    constexpr double outer_temperature = 600.0;

    const fuelsim::SteadySingleRegionParameters parameters = {
        0.0,
        radius,
        length,
        32,
        2,
        constant_material(conductivity, 0.0),
        heat_source,
        outer_temperature,
        outer_temperature,
        0.0,
        0.0,
    };

    fuelsim::SteadySingleRegionProblem problem(parameters);
    fuelsim::PetscSolver solver;
    fuelsim::SolverOptions scaled_options;
    scaled_options.field_residual_scaling = true;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state(), scaled_options);

    bool passed = check(result.converged, "thermal cylinder SNES converged");
    const std::size_t expected_begin =
        problem.contribution_count() *
        static_cast<std::size_t>(result.mpi_rank) /
        static_cast<std::size_t>(result.mpi_size);
    const std::size_t expected_end =
        problem.contribution_count() *
        static_cast<std::size_t>(result.mpi_rank + 1) /
        static_cast<std::size_t>(result.mpi_size);
    passed =
        check(result.local_contribution_begin == expected_begin &&
                  result.local_contribution_end == expected_end,
              "PETSc rank owns its exact nonoverlapping contribution range") &&
        passed;
    if (result.mpi_size > 1)
        passed =
            check(result.local_contribution_end -
                          result.local_contribution_begin <
                      problem.contribution_count(),
                  "MPI rank does not repeat the complete model assembly") &&
            passed;
    double maximum_scaled_error = 0.0;
    const double center_rise =
        heat_source * radius * radius / (4.0 * conductivity);
    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const double r = problem.mesh().nodes()[node].r;
        const double expected =
            outer_temperature +
            heat_source * (radius * radius - r * r) / (4.0 * conductivity);
        const double actual = result.state[problem.dof_map().temperature(node)];
        maximum_scaled_error = std::max(
            maximum_scaled_error, std::abs(actual - expected) / center_rise);
    }
    passed = check(maximum_scaled_error < 1.0e-3,
                   "thermal cylinder temperature error is below 0.1%; actual=" +
                       std::to_string(maximum_scaled_error)) &&
             passed;
    std::cout << "thermal_cylinder_maximum_scaled_error="
              << maximum_scaled_error << '\n';
    return passed;
}

double thermal_cylinder_error(std::size_t radial_elements) {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double conductivity = 4.0;
    constexpr double heat_source = 2.0e8;
    constexpr double outer_temperature = 600.0;
    const fuelsim::SteadySingleRegionParameters parameters = {
        0.0,
        radius,
        length,
        radial_elements,
        2,
        constant_material(conductivity, 0.0),
        heat_source,
        outer_temperature,
        outer_temperature,
        0.0,
        0.0,
    };
    fuelsim::SteadySingleRegionProblem problem(parameters);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());
    if (!result.converged)
        throw std::runtime_error(
            "thermal mesh-convergence solve did not converge");
    const double center_rise =
        heat_source * radius * radius / (4.0 * conductivity);
    double maximum_error = 0.0;
    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const double r = problem.mesh().nodes()[node].r;
        const double expected =
            outer_temperature +
            heat_source * (radius * radius - r * r) /
                (4.0 * conductivity);
        const double actual =
            result.state[problem.dof_map().temperature(node)];
        maximum_error = std::max(maximum_error,
                                 std::abs(actual - expected) / center_rise);
    }
    return maximum_error;
}

bool test_thermal_mesh_convergence() {
    const std::array<double, 3> errors = {
        thermal_cylinder_error(8), thermal_cylinder_error(16),
        thermal_cylinder_error(32)};
    const double first_ratio = errors[0] / errors[1];
    const double second_ratio = errors[1] / errors[2];
    const double first_order = std::log2(first_ratio);
    const double second_order = std::log2(second_ratio);
    std::cout << "thermal_mesh_convergence_errors=" << errors[0] << ','
              << errors[1] << ',' << errors[2] << '\n';
    std::cout << "thermal_mesh_convergence_ratios=" << first_ratio << ','
              << second_ratio << '\n';
    std::cout << "thermal_mesh_convergence_orders=" << first_order << ','
              << second_order << '\n';
    return check(first_order > 1.7 && second_order > 1.7 &&
                     second_order > first_order,
                 "successive radial mesh refinement approaches second-order "
                 "thermal convergence");
}

bool test_free_thermal_expansion() {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double alpha = 1.0e-5;
    constexpr double temperature = 700.0;
    constexpr double temperature_change = 100.0;

    const fuelsim::SteadySingleRegionParameters parameters = {
        0.0, radius,      length, 8,   4,   constant_material(4.0, alpha),
        0.0, temperature, 600.0,  0.0, 0.0,
    };

    const fuelsim::SteadySingleRegionProblem problem(parameters);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    bool passed = check(result.converged, "free expansion SNES converged");
    double maximum_temperature_error = 0.0;
    double maximum_displacement_error = 0.0;
    const double displacement_scale = alpha * temperature_change * length;
    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const fuelsim::RzPoint& point = problem.mesh().nodes()[node];
        const double actual_temperature =
            result.state[problem.dof_map().temperature(node)];
        const double actual_radial =
            result.state[problem.dof_map().radial_displacement(node)];
        const double actual_axial =
            result.state[problem.dof_map().axial_displacement(node)];
        maximum_temperature_error =
            std::max(maximum_temperature_error,
                     std::abs(actual_temperature - temperature));
        maximum_displacement_error = std::max(
            maximum_displacement_error,
            std::abs(actual_radial - alpha * temperature_change * point.r));
        maximum_displacement_error = std::max(
            maximum_displacement_error,
            std::abs(actual_axial - alpha * temperature_change * point.z));
    }

    double maximum_stress = 0.0;
    for (std::size_t element = 0; element < problem.element_count();
         ++element) {
        const fuelsim::LocalValues state =
            problem.element_state(element, result.state);
        const auto stresses = problem.kernel().stress_values(
            problem.element_geometry(element), state);
        for (const fuelsim::AxisymmetricStressValues& stress : stresses) {
            maximum_stress = std::max(maximum_stress, std::abs(stress.rr));
            maximum_stress = std::max(maximum_stress, std::abs(stress.zz));
            maximum_stress = std::max(maximum_stress, std::abs(stress.hoop));
            maximum_stress = std::max(maximum_stress, std::abs(stress.rz));
        }
    }

    passed = check(maximum_temperature_error < 1.0e-9,
                   "free expansion temperature is uniform") &&
             passed;
    passed = check(maximum_displacement_error / displacement_scale < 1.0e-8,
                   "free expansion displacement matches analytic field") &&
             passed;
    passed = check(maximum_stress < 100.0,
                   "free expansion stress is below 100 Pa") &&
             passed;
    std::cout << "free_expansion_maximum_temperature_error="
              << maximum_temperature_error << '\n';
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

    fuelsim::ThermoelasticProperties material = constant_material(4.0, 0.0);
    material.young_modulus = young_modulus;
    material.poisson_ratio = poisson_ratio;

    const fuelsim::SteadySingleRegionParameters parameters = {
        inner_radius, outer_radius, length, 48,       2,   material,
        0.0,          600.0,        600.0,  pressure, 0.0,
    };

    const fuelsim::SteadySingleRegionProblem problem(parameters);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    bool passed = check(result.converged, "Lame cylinder SNES converged");
    const double A =
        pressure * inner_radius * inner_radius /
        (outer_radius * outer_radius - inner_radius * inner_radius);
    std::cout << "lame_initial_field_residuals="
              << result.initial_field_residual_norms[0] << ','
              << result.initial_field_residual_norms[1] << ','
              << result.initial_field_residual_norms[2] << '\n';
    std::cout << "lame_final_scaled_field_residuals="
              << result.final_scaled_field_residual_norms[0] << ','
              << result.final_scaled_field_residual_norms[1] << ','
              << result.final_scaled_field_residual_norms[2] << '\n';
    const double B =
        pressure * inner_radius * inner_radius * outer_radius * outer_radius /
        (outer_radius * outer_radius - inner_radius * inner_radius);
    const double axial_strain = -2.0 * poisson_ratio * A / young_modulus;

    double maximum_radial_relative_error = 0.0;
    double maximum_axial_relative_error = 0.0;
    const double radial_scale = ((1.0 - poisson_ratio) * A * inner_radius +
                                 (1.0 + poisson_ratio) * B / inner_radius) /
                                young_modulus;
    const double axial_scale = std::abs(axial_strain * length);

    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const fuelsim::RzPoint& point = problem.mesh().nodes()[node];
        const double expected_radial = ((1.0 - poisson_ratio) * A * point.r +
                                        (1.0 + poisson_ratio) * B / point.r) /
                                       young_modulus;
        const double expected_axial = axial_strain * point.z;
        const double actual_radial =
            result.state[problem.dof_map().radial_displacement(node)];
        const double actual_axial =
            result.state[problem.dof_map().axial_displacement(node)];
        maximum_radial_relative_error =
            std::max(maximum_radial_relative_error,
                     std::abs(actual_radial - expected_radial) / radial_scale);
        if (axial_scale > 0.0) {
            maximum_axial_relative_error =
                std::max(maximum_axial_relative_error,
                         std::abs(actual_axial - expected_axial) / axial_scale);
        }
    }

    passed = check(maximum_radial_relative_error < 3.0e-3,
                   "Lame radial displacement error is below 0.3%") &&
             passed;
    passed = check(maximum_axial_relative_error < 3.0e-3,
                   "Lame axial displacement error is below 0.3%") &&
             passed;
    std::cout << "lame_maximum_radial_relative_error="
              << maximum_radial_relative_error << '\n';
    std::cout << "lame_maximum_axial_relative_error="
              << maximum_axial_relative_error << '\n';
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

    const fuelsim::SteadyFuelCladdingParameters parameters = {
        fuel_radius,
        cladding_inner_radius,
        cladding_outer_radius,
        length,
        length,
        32,
        8,
        2,
        constant_material(fuel_conductivity, 0.0),
        constant_material(cladding_conductivity, 0.0),
        heat_source,
        outer_temperature,
        outer_temperature,
        gap_conductivity,
        1.0e-6,
        1.0e14,
    };
    const fuelsim::SteadyFuelCladdingProblem problem(parameters);
    fuelsim::PetscSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    const double cladding_rise =
        heat_source * fuel_radius * fuel_radius /
        (2.0 * cladding_conductivity) *
        std::log(cladding_outer_radius / cladding_inner_radius);
    const double gap_rise = heat_source * fuel_radius *
                            (cladding_inner_radius - fuel_radius) /
                            (2.0 * gap_conductivity);
    const double fuel_rise =
        heat_source * fuel_radius * fuel_radius / (4.0 * fuel_conductivity);
    const double expected_cladding_inner = outer_temperature + cladding_rise;
    const double expected_fuel_surface = expected_cladding_inner + gap_rise;
    const double expected_center = expected_fuel_surface + fuel_rise;

    const std::size_t axial_mid = problem.fuel_mesh().axial_elements() / 2;
    const std::size_t fuel_center_local =
        problem.fuel_mesh().node_id(0, axial_mid);
    const std::size_t fuel_surface_local = problem.fuel_mesh().node_id(
        problem.fuel_mesh().radial_elements(), axial_mid);
    const std::size_t cladding_inner_local =
        problem.cladding_mesh().node_id(0, axial_mid);
    const fuelsim::DofMap& dofs = problem.dof_map();

    const double actual_center = result.state[dofs.temperature(
        problem.fuel_global_node(fuel_center_local))];
    const double actual_fuel_surface = result.state[dofs.temperature(
        problem.fuel_global_node(fuel_surface_local))];
    const double actual_cladding_inner = result.state[dofs.temperature(
        problem.cladding_global_node(cladding_inner_local))];
    const double temperature_scale = expected_center - outer_temperature;
    const double center_error =
        std::abs(actual_center - expected_center) / temperature_scale;
    const double fuel_surface_error =
        std::abs(actual_fuel_surface - expected_fuel_surface) /
        temperature_scale;
    const double cladding_inner_error =
        std::abs(actual_cladding_inner - expected_cladding_inner) /
        temperature_scale;

    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.state);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double expected_heat_rate =
        heat_source * pi * fuel_radius * fuel_radius * length;
    const double heat_balance_error =
        metric_relative_error(interface.total_heat_rate, expected_heat_rate);

    bool passed = check(result.converged, "M1 open-gap thermal SNES converged");
    passed = check(center_error < 1.0e-3,
                   "M1 analytic center temperature error is below 0.1%") &&
             passed;
    passed = check(fuel_surface_error < 1.0e-3,
                   "M1 analytic fuel-surface temperature error is below "
                   "0.1%") &&
             passed;
    passed = check(cladding_inner_error < 1.0e-3,
                   "M1 analytic cladding-inner temperature error is below "
                   "0.1%") &&
             passed;
    passed = check(interface.minimum_gap > 0.0 &&
                       interface.maximum_contact_pressure == 0.0,
                   "M1 analytic thermal case remains out of contact") &&
             passed;
    passed = check(heat_balance_error < 1.0e-9,
                   "M1 interface heat rate balances generated power") &&
             passed;

    std::cout << "m1_analytic_center_temperature_scaled_error=" << center_error
              << '\n';
    std::cout << "m1_analytic_fuel_surface_scaled_error=" << fuel_surface_error
              << '\n';
    std::cout << "m1_analytic_cladding_inner_scaled_error="
              << cladding_inner_error << '\n';
    std::cout << "m1_analytic_heat_balance_relative_error="
              << heat_balance_error << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M0 and M1 numerical acceptance tests\n");

        bool passed = true;
        passed = test_shadow_state_view() && passed;
        passed = test_global_newton_safeguards() && passed;
        passed = test_thermal_cylinder() && passed;
        passed = test_thermal_mesh_convergence() && passed;
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
