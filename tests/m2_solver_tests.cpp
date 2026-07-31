#include "fuelsim/m2_problem.hpp"
#include "fuelsim/m2_solver.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/quad4_rz_transient.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
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

fuelsim::ThermoelasticProperties constant_material(double conductivity,
                                                   double young_modulus,
                                                   double poisson_ratio,
                                                   double thermal_expansion) {
    return {
        0.0,           conductivity,      young_modulus,
        poisson_ratio, thermal_expansion, 600.0,
    };
}

fuelsim::TransientInelasticProperties elastic_transient_properties() {
    return {
        10000.0,         300.0,      fuelsim::InelasticBehavior::elastic,
        {0.0, 1.0, 1.0}, {1.0, 0.0},
    };
}

class SingleElementTransientProblem final : public fuelsim::NonlinearProblem {
  public:
    SingleElementTransientProblem()
        : _geometry(fuelsim::make_quad4_rz_geometry({{
              {0.0, 0.0},
              {0.004, 0.0},
              {0.004, 0.010},
              {0.0, 0.010},
          }})),
          _kernel(fuelsim::IsotropicInelasticMaterial(
                      constant_material(2.0, 2.0e11, 0.3, 0.0),
                      elastic_transient_properties()),
                  3.0e6),
          _committed_temperature{600.0, 600.0, 600.0, 600.0},
          _committed_material{}, _dofs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
          _dirichlet_conditions{
              {4, 0.0}, {5, 0.0}, {6, 0.0},  {7, 0.0},
              {8, 0.0}, {9, 0.0}, {10, 0.0}, {11, 0.0},
          },
          _time_step(1.0) {}

    std::vector<double> initial_state() const {
        std::vector<double> state(fuelsim::local_dof_count, 0.0);
        for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node)
            state[node] = _committed_temperature[node];
        return state;
    }

    void commit(const std::vector<double>& converged_state) {
        if (converged_state.size() != dof_count())
            throw std::invalid_argument(
                "SingleElementTransientProblem committed state size "
                "mismatch");

        fuelsim::LocalValues local_state{};
        std::copy(converged_state.begin(), converged_state.end(),
                  local_state.begin());
        _committed_material = _kernel.trial_state_values(
            _geometry, local_state, _committed_material, _time_step);
        for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node)
            _committed_temperature[node] = converged_state[node];
    }

    const fuelsim::Quad4TemperatureHistory&
    committed_temperature() const noexcept {
        return _committed_temperature;
    }

    std::size_t dof_count() const noexcept override {
        return fuelsim::local_dof_count;
    }

    std::size_t contribution_count() const noexcept override {
        return 1;
    }

    fuelsim::LocalDofs
    contribution_dofs(std::size_t contribution_index) const override {
        if (contribution_index != 0)
            throw std::out_of_range(
                "SingleElementTransientProblem contribution is out of "
                "range");
        return _dofs;
    }

    fuelsim::LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const fuelsim::LocalValues& state) const override {
        if (contribution_index != 0)
            throw std::out_of_range(
                "SingleElementTransientProblem contribution is out of "
                "range");
        return _kernel.residual(_geometry, state, _committed_temperature,
                                _committed_material, _time_step);
    }

    fuelsim::LocalSystem
    linearize_contribution(std::size_t contribution_index,
                           const fuelsim::LocalValues& state) const override {
        if (contribution_index != 0)
            throw std::out_of_range(
                "SingleElementTransientProblem contribution is out of "
                "range");
        return _kernel.linearize(_geometry, state, _committed_temperature,
                                 _committed_material, _time_step);
    }

    const std::vector<fuelsim::DirichletCondition>&
    dirichlet_conditions() const noexcept override {
        return _dirichlet_conditions;
    }

  protected:
    void add_state_independent_residual(
        std::vector<double>& residual) const override {
        (void)residual;
    }

  private:
    fuelsim::Quad4RzGeometry _geometry;
    fuelsim::Quad4RzTransientKernel _kernel;
    fuelsim::Quad4TemperatureHistory _committed_temperature;
    fuelsim::Quad4MaterialHistory _committed_material;
    fuelsim::LocalDofs _dofs;
    std::vector<fuelsim::DirichletCondition> _dirichlet_conditions;
    double _time_step;
};

bool test_single_element_backward_euler_heat_source() {
    constexpr std::size_t step_count = 10;
    constexpr double expected_temperature = 610.0;
    // verification/moose/m21_transient_heat_rz_out.csv at t = 10 s.
    constexpr double moose_snapshot_temperature = 610.0;

    SingleElementTransientProblem problem;
    fuelsim::PetscSequentialSolver solver;
    std::vector<double> state = problem.initial_state();
    std::size_t workspace_setups = 0;
    std::size_t solve_calls = 0;
    bool passed = true;

    for (std::size_t step = 0; step < step_count; ++step) {
        const fuelsim::SolveResult result = solver.solve(problem, state);
        passed = check(result.converged, "single-element backward-Euler step " +
                                             std::to_string(step + 1) +
                                             " converged") &&
                 passed;
        workspace_setups += result.timing.workspace_setups;
        solve_calls += result.timing.solve_calls;
        if (!result.converged)
            break;

        problem.commit(result.state);
        state = result.state;
    }

    double maximum_temperature_error = 0.0;
    double average_temperature = 0.0;
    for (double temperature : problem.committed_temperature()) {
        maximum_temperature_error =
            std::max(maximum_temperature_error,
                     std::abs(temperature - expected_temperature));
        average_temperature += temperature;
    }
    average_temperature /= static_cast<double>(fuelsim::quad4_node_count);
    const double moose_snapshot_relative_error =
        std::abs(average_temperature - moose_snapshot_temperature) /
        moose_snapshot_temperature;

    passed = check(maximum_temperature_error < 1.0e-9,
                   "single-element temperature reaches 610 K") &&
             passed;
    passed = check(moose_snapshot_relative_error < 1.0e-3,
                   "single-element temperature differs from the MOOSE "
                   "snapshot by less than 0.1%") &&
             passed;
    passed = check(workspace_setups == 1,
                   "ten backward-Euler steps create one PETSc workspace") &&
             passed;
    passed = check(solve_calls == step_count,
                   "one PETSc solve is issued per backward-Euler step") &&
             passed;

    std::cout << "m2_single_element_temperature=" << average_temperature
              << '\n';
    std::cout << "m2_single_element_moose_relative_error="
              << moose_snapshot_relative_error << '\n';
    std::cout << "m2_single_element_workspace_setups=" << workspace_setups
              << '\n';
    std::cout << "m2_single_element_solve_calls=" << solve_calls << '\n';
    return passed;
}

fuelsim::M2Parameters zero_source_m2_parameters() {
    const fuelsim::ThermoelasticProperties fuel =
        constant_material(4.0, 2.0e11, 0.3, 1.0e-5);
    const fuelsim::ThermoelasticProperties cladding =
        constant_material(16.0, 7.5e10, 0.3, 5.0e-6);
    const fuelsim::TransientInelasticProperties transient =
        elastic_transient_properties();

    return {
        {
            0.0040,
            0.0041,
            0.0046,
            0.010,
            0.01002,
            1,
            1,
            1,
            fuel,
            cladding,
            0.0,
            600.0,
            600.0,
            0.4,
            1.0e-6,
            1.0e14,
        },
        transient,
        transient,
    };
}

bool material_history_is_zero(const fuelsim::Quad4MaterialHistory& history) {
    for (const fuelsim::MaterialPointState& point : history) {
        if (point.equivalent_plastic_strain != 0.0 ||
            point.equivalent_creep_strain != 0.0)
            return false;
        for (std::size_t component = 0; component < 4; ++component) {
            if (point.plastic_strain[component] != 0.0 ||
                point.creep_strain[component] != 0.0)
                return false;
        }
    }
    return true;
}

bool m2_histories_are_zero(const fuelsim::M2Problem& problem) {
    for (std::size_t element = 0;
         element < problem.base_problem().fuel_element_count(); ++element) {
        if (!material_history_is_zero(problem.fuel_material_history(element)))
            return false;
    }
    for (std::size_t element = 0;
         element < problem.base_problem().cladding_element_count(); ++element) {
        if (!material_history_is_zero(
                problem.cladding_material_history(element)))
            return false;
    }
    return true;
}

bool test_m2_zero_source_history_and_interface() {
    fuelsim::M2Problem problem(zero_source_m2_parameters());
    const std::vector<double> initial_state = problem.committed_solution();
    const fuelsim::M2TimeOptions time_options = {
        2.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0, 0.0,
    };

    const fuelsim::M2TimeStepper time_stepper;
    const fuelsim::M2TransientResult result =
        time_stepper.solve(problem, time_options);

    bool passed = check(result.completed && result.last_attempt.converged,
                        "M2 zero-source transient completes two steps");
    passed = check(result.accepted_steps.size() == 2,
                   "M2 zero-source transient commits two accepted steps") &&
             passed;
    passed =
        check(result.committed_time == 2.0 && problem.committed_time() == 2.0,
              "M2 zero-source transient commits time 2 s") &&
        passed;
    passed = check(!problem.time_step_active(),
                   "M2 has no active transaction after completion") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "M2 time stepper creates one PETSc workspace") &&
             passed;
    passed = check(result.aggregate_timing.solve_calls == 2,
                   "M2 time stepper reuses the workspace for two solves") &&
             passed;

    if (result.accepted_steps.size() == 2) {
        passed = check(result.accepted_steps[0].time == 1.0 &&
                           result.accepted_steps[1].time == 2.0 &&
                           result.accepted_steps[0].time_step == 1.0 &&
                           result.accepted_steps[1].time_step == 1.0,
                       "M2 accepted-step times and increments are committed") &&
                 passed;
        passed =
            check(result.accepted_steps[0].volumetric_heat_source == 0.0 &&
                      result.accepted_steps[1].volumetric_heat_source == 0.0,
                  "M2 accepted steps retain the zero heat source") &&
            passed;
    }

    double maximum_state_drift = 0.0;
    for (std::size_t dof = 0; dof < initial_state.size(); ++dof) {
        maximum_state_drift =
            std::max(maximum_state_drift, std::abs(result.committed_state[dof] -
                                                   initial_state[dof]));
    }
    passed = check(maximum_state_drift < 1.0e-12,
                   "M2 zero-source committed solution does not drift") &&
             passed;
    passed = check(m2_histories_are_zero(problem),
                   "M2 zero-source material history does not drift") &&
             passed;

    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.committed_state);
    const std::size_t expected_projected_nodes =
        problem.parameters().base.axial_elements + 1;
    passed =
        check(interface.projected_contact_nodes == expected_projected_nodes,
              "M2 preserves every fuel-surface contact projection") &&
        passed;
    passed = check(std::isfinite(interface.minimum_gap) &&
                       interface.minimum_gap > 0.0 &&
                       std::isfinite(interface.minimum_contact_gap) &&
                       interface.minimum_contact_gap > 0.0 &&
                       interface.active_contact_nodes == 0,
                   "M2 zero-source interface remains a finite open gap") &&
             passed;
    passed = check(std::abs(interface.total_heat_rate) < 1.0e-12,
                   "M2 zero-source interface has zero heat transfer") &&
             passed;

    std::cout << "m2_zero_source_committed_time=" << result.committed_time
              << '\n';
    std::cout << "m2_zero_source_accepted_steps="
              << result.accepted_steps.size() << '\n';
    std::cout << "m2_zero_source_workspace_setups="
              << result.aggregate_timing.workspace_setups << '\n';
    std::cout << "m2_zero_source_maximum_state_drift=" << maximum_state_drift
              << '\n';
    std::cout << "m2_zero_source_projected_contact_nodes="
              << interface.projected_contact_nodes << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M2 transient solver acceptance tests\n");

        bool passed = true;
        passed = test_single_element_backward_euler_heat_source() && passed;
        passed = test_m2_zero_source_history_and_interface() && passed;
        if (!passed)
            return 1;

        std::cout << "[PASS] fuelsim M2 transient solver acceptance tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M2 solver tests raised: " << error.what() << '\n';
        return 1;
    }
}
