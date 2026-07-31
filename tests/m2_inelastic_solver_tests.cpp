#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/quad4_rz_transient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double radius = 1.0e-3;
constexpr double height = 1.0e-3;
constexpr double temperature = 600.0;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double moose_relative_tolerance = 1.0e-3;

enum class LoadingMode {
    prescribed_top_displacement,
    axial_traction,
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double relative_error(double actual, double expected) {
    if (!std::isfinite(actual) || !std::isfinite(expected) || expected == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - expected) / std::abs(expected);
}

fuelsim::ThermoelasticProperties thermoelastic_properties() {
    return {
        0.0, 1.0, 2.0e11, 0.3, 0.0, temperature,
    };
}

fuelsim::TransientInelasticProperties j2_properties() {
    return {
        1.0,
        1.0,
        fuelsim::InelasticBehavior::j2_plasticity,
        {0.0, 1.0, 1.0},
        {2.0e8, 2.0e9},
    };
}

fuelsim::TransientInelasticProperties norton_properties() {
    return {
        1.0,
        1.0,
        fuelsim::InelasticBehavior::norton_creep,
        {1.0e-30, 1.0, 3.0},
        {1.0, 0.0},
    };
}

bool histories_equal(const fuelsim::Quad4MaterialHistory& lhs,
                     const fuelsim::Quad4MaterialHistory& rhs) {
    for (std::size_t point = 0; point < lhs.size(); ++point) {
        if (lhs[point].equivalent_plastic_strain !=
                rhs[point].equivalent_plastic_strain ||
            lhs[point].equivalent_creep_strain !=
                rhs[point].equivalent_creep_strain)
            return false;
        for (std::size_t component = 0; component < 4; ++component) {
            if (lhs[point].plastic_strain[component] !=
                    rhs[point].plastic_strain[component] ||
                lhs[point].creep_strain[component] !=
                    rhs[point].creep_strain[component])
                return false;
        }
    }
    return true;
}

class SingleElementInelasticProblem final : public fuelsim::NonlinearProblem {
  public:
    SingleElementInelasticProblem(
        fuelsim::TransientInelasticProperties properties,
        LoadingMode loading_mode, double time_step)
        : _geometry(fuelsim::make_quad4_rz_geometry({{
              {0.0, 0.0},
              {radius, 0.0},
              {radius, height},
              {0.0, height},
          }})),
          _kernel(fuelsim::IsotropicInelasticMaterial(
                      thermoelastic_properties(), properties),
                  0.0),
          _committed_temperature{
              temperature,
              temperature,
              temperature,
              temperature,
          },
          _committed_material{}, _last_stress{},
          _dofs{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
          _dirichlet_conditions{
              {4, 0.0},
              {7, 0.0},
              {8, 0.0},
              {9, 0.0},
          },
          _loading_mode(loading_mode), _time_step(time_step),
          _traction(loading_mode == LoadingMode::axial_traction ? 1.0e8 : 0.0),
          _committed_steps(0) {
        if (!std::isfinite(_time_step) || !(_time_step > 0.0))
            throw std::invalid_argument(
                "SingleElementInelasticProblem time step must be finite and "
                "positive");
        if (_loading_mode == LoadingMode::prescribed_top_displacement) {
            _dirichlet_conditions.push_back({10, 0.0});
            _dirichlet_conditions.push_back({11, 0.0});
        }
    }

    std::vector<double> initial_state() const {
        std::vector<double> state(fuelsim::local_dof_count, 0.0);
        for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node)
            state[node] = _committed_temperature[node];
        return state;
    }

    void set_step_end_time(double end_time) {
        if (!std::isfinite(end_time) || !(end_time > 0.0))
            throw std::invalid_argument(
                "SingleElementInelasticProblem end time must be finite and "
                "positive");
        if (_loading_mode != LoadingMode::prescribed_top_displacement)
            return;

        const double top_displacement = height * 0.002 * end_time;
        _dirichlet_conditions[4].value = top_displacement;
        _dirichlet_conditions[5].value = top_displacement;
    }

    void apply_dirichlet_values(std::vector<double>& state) const {
        if (state.size() != dof_count())
            throw std::invalid_argument(
                "SingleElementInelasticProblem state size mismatch");
        for (const fuelsim::DirichletCondition& condition :
             _dirichlet_conditions)
            state[condition.dof] = condition.value;
    }

    void commit(const std::vector<double>& converged_state) {
        if (converged_state.size() != dof_count())
            throw std::invalid_argument(
                "SingleElementInelasticProblem committed state size "
                "mismatch");

        fuelsim::LocalValues local_state{};
        std::copy(converged_state.begin(), converged_state.end(),
                  local_state.begin());
        const std::array<fuelsim::AxisymmetricStressValues, 4> staged_stress =
            _kernel.stress_values(_geometry, local_state, _committed_material,
                                  _time_step);
        const fuelsim::Quad4MaterialHistory staged_material =
            _kernel.trial_state_values(_geometry, local_state,
                                       _committed_material, _time_step);

        for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node) {
            const double value = converged_state[node];
            if (!std::isfinite(value) || !(value > 0.0))
                throw std::domain_error(
                    "SingleElementInelasticProblem committed temperature "
                    "must be finite and positive");
            _committed_temperature[node] = value;
        }
        _last_stress = staged_stress;
        _committed_material = staged_material;
        ++_committed_steps;
    }

    const fuelsim::Quad4MaterialHistory& committed_material() const noexcept {
        return _committed_material;
    }

    const fuelsim::Quad4TemperatureHistory&
    committed_temperature() const noexcept {
        return _committed_temperature;
    }

    const std::array<fuelsim::AxisymmetricStressValues, 4>&
    last_stress() const noexcept {
        return _last_stress;
    }

    std::size_t committed_steps() const noexcept {
        return _committed_steps;
    }

    std::size_t dof_count() const noexcept override {
        return fuelsim::local_dof_count;
    }

    std::size_t contribution_count() const noexcept override {
        return 1;
    }

    fuelsim::LocalDofs
    contribution_dofs(std::size_t contribution_index) const override {
        require_single_contribution(contribution_index);
        return _dofs;
    }

    fuelsim::LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const fuelsim::LocalValues& state) const override {
        require_single_contribution(contribution_index);
        return _kernel.residual(_geometry, state, _committed_temperature,
                                _committed_material, _time_step);
    }

    fuelsim::LocalSystem
    linearize_contribution(std::size_t contribution_index,
                           const fuelsim::LocalValues& state) const override {
        require_single_contribution(contribution_index);
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
        if (_loading_mode != LoadingMode::axial_traction)
            return;
        if (residual.size() != dof_count())
            throw std::logic_error(
                "SingleElementInelasticProblem residual size mismatch");

        const double top_area = pi * radius * radius;
        residual[10] -= _traction * 2.0 * top_area / 3.0;
        residual[11] -= _traction * top_area / 3.0;
    }

  private:
    static void require_single_contribution(std::size_t index) {
        if (index != 0)
            throw std::out_of_range(
                "SingleElementInelasticProblem contribution is out of "
                "range");
    }

    fuelsim::Quad4RzGeometry _geometry;
    fuelsim::Quad4RzTransientKernel _kernel;
    fuelsim::Quad4TemperatureHistory _committed_temperature;
    fuelsim::Quad4MaterialHistory _committed_material;
    std::array<fuelsim::AxisymmetricStressValues, 4> _last_stress;
    fuelsim::LocalDofs _dofs;
    std::vector<fuelsim::DirichletCondition> _dirichlet_conditions;
    LoadingMode _loading_mode;
    double _time_step;
    double _traction;
    std::size_t _committed_steps;
};

double average_axial_stress(
    const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses) {
    double sum = 0.0;
    for (const fuelsim::AxisymmetricStressValues& stress : stresses)
        sum += stress.zz;
    return sum / static_cast<double>(stresses.size());
}

double average_equivalent_plastic_strain(
    const fuelsim::Quad4MaterialHistory& history) {
    double sum = 0.0;
    for (const fuelsim::MaterialPointState& point : history)
        sum += point.equivalent_plastic_strain;
    return sum / static_cast<double>(history.size());
}

double
average_equivalent_creep_strain(const fuelsim::Quad4MaterialHistory& history) {
    double sum = 0.0;
    for (const fuelsim::MaterialPointState& point : history)
        sum += point.equivalent_creep_strain;
    return sum / static_cast<double>(history.size());
}

double maximum_inelastic_trace(const fuelsim::Quad4MaterialHistory& history,
                               LoadingMode loading_mode) {
    double maximum = 0.0;
    for (const fuelsim::MaterialPointState& point : history) {
        const std::array<double, 4>& strain =
            loading_mode == LoadingMode::prescribed_top_displacement
                ? point.plastic_strain
                : point.creep_strain;
        maximum =
            std::max(maximum, std::abs(strain[0] + strain[1] + strain[2]));
    }
    return maximum;
}

bool temperatures_are_committed(
    const fuelsim::Quad4TemperatureHistory& history) {
    for (double value : history) {
        if (std::abs(value - temperature) > 1.0e-12)
            return false;
    }
    return true;
}

struct LoadPathResult final {
    std::vector<double> state;
    std::size_t workspace_setups;
    std::size_t solve_calls;
    bool converged;
    bool callbacks_preserved_history;
};

LoadPathResult solve_load_path(SingleElementInelasticProblem& problem,
                               std::size_t step_count, double time_step) {
    fuelsim::PetscSequentialSolver solver;
    std::vector<double> state = problem.initial_state();
    std::size_t workspace_setups = 0;
    std::size_t solve_calls = 0;
    bool converged = true;
    bool callbacks_preserved_history = true;
    const fuelsim::SolverOptions options = {
        1.0e-10,
        1.0e-10,
        1.0e-12,
        40,
    };

    for (std::size_t step = 0; step < step_count; ++step) {
        const double end_time = static_cast<double>(step + 1) * time_step;
        problem.set_step_end_time(end_time);
        problem.apply_dirichlet_values(state);
        const fuelsim::Quad4MaterialHistory history_before =
            problem.committed_material();
        const fuelsim::SolveResult result =
            solver.solve(problem, state, options);
        workspace_setups += result.timing.workspace_setups;
        solve_calls += result.timing.solve_calls;
        callbacks_preserved_history =
            histories_equal(history_before, problem.committed_material()) &&
            callbacks_preserved_history;
        if (!result.converged) {
            std::cerr << "[FAIL] load step " << step + 1 << " PETSc reason="
                      << fuelsim::petsc_convergence_reason_name(
                             result.convergence_reason)
                      << ", iterations=" << result.nonlinear_iterations
                      << ", residual_norm=" << result.residual_norm << '\n';
            converged = false;
            state = result.state;
            break;
        }

        problem.commit(result.state);
        state = result.state;
    }

    return {
        state,
        workspace_setups,
        solve_calls,
        converged,
        callbacks_preserved_history,
    };
}

bool test_j2_moose_comparison() {
    constexpr std::size_t step_count = 10;
    constexpr double time_step = 0.1;
    constexpr double moose_axial_stress = 201980198.0198;
    constexpr double moose_equivalent_plastic = 0.000990099009901;
    constexpr double moose_top_displacement = 2.0e-6;

    SingleElementInelasticProblem problem(
        j2_properties(), LoadingMode::prescribed_top_displacement, time_step);
    const LoadPathResult result =
        solve_load_path(problem, step_count, time_step);

    const double axial_stress = average_axial_stress(problem.last_stress());
    const double equivalent_plastic =
        average_equivalent_plastic_strain(problem.committed_material());
    const double top_displacement = 0.5 * (result.state[10] + result.state[11]);
    const double top_displacement_spread =
        std::abs(result.state[10] - result.state[11]);
    const double stress_error =
        relative_error(axial_stress, moose_axial_stress);
    const double plastic_error =
        relative_error(equivalent_plastic, moose_equivalent_plastic);
    const double displacement_error =
        relative_error(top_displacement, moose_top_displacement);
    const double maximum_trace = maximum_inelastic_trace(
        problem.committed_material(), LoadingMode::prescribed_top_displacement);

    bool passed =
        check(result.converged, "J2 ten-step PETSc load path converged");
    passed = check(result.workspace_setups == 1,
                   "J2 load path creates one reusable PETSc workspace") &&
             passed;
    passed = check(result.solve_calls == step_count,
                   "J2 load path issues ten PETSc solve calls") &&
             passed;
    passed =
        check(result.callbacks_preserved_history,
              "J2 residual and Jacobian callbacks do not commit history") &&
        passed;
    passed = check(problem.committed_steps() == step_count,
                   "J2 history is committed after every accepted step") &&
             passed;
    passed = check(temperatures_are_committed(problem.committed_temperature()),
                   "J2 committed temperature remains 600 K") &&
             passed;
    passed = check(stress_error < moose_relative_tolerance,
                   "J2 final axial stress matches MOOSE within 0.1%") &&
             passed;
    passed = check(plastic_error < moose_relative_tolerance,
                   "J2 final equivalent plastic strain matches MOOSE within "
                   "0.1%") &&
             passed;
    passed = check(displacement_error < moose_relative_tolerance,
                   "J2 final top displacement matches MOOSE within 0.1%") &&
             passed;
    passed = check(top_displacement_spread < 1.0e-12,
                   "J2 top-side axial displacement is uniform") &&
             passed;
    passed = check(maximum_trace < 1.0e-12,
                   "J2 committed plastic strain is trace-free") &&
             passed;

    std::cout << "m22_j2_axial_stress=" << axial_stress << '\n';
    std::cout << "m22_j2_equivalent_plastic=" << equivalent_plastic << '\n';
    std::cout << "m22_j2_top_displacement=" << top_displacement << '\n';
    std::cout << "m22_j2_stress_relative_error=" << stress_error << '\n';
    std::cout << "m22_j2_plastic_relative_error=" << plastic_error << '\n';
    std::cout << "m22_j2_displacement_relative_error=" << displacement_error
              << '\n';
    std::cout << "m22_j2_top_displacement_spread=" << top_displacement_spread
              << '\n';
    std::cout << "m22_j2_maximum_plastic_trace=" << maximum_trace << '\n';
    std::cout << "m22_j2_workspace_setups=" << result.workspace_setups << '\n';
    return passed;
}

bool test_norton_moose_comparison() {
    constexpr std::size_t step_count = 10;
    constexpr double time_step = 10.0;
    constexpr double moose_axial_stress = 99998007.620195;
    constexpr double moose_equivalent_creep = 9.9991036503199e-5;
    constexpr double moose_top_displacement = 5.9997808603447e-7;

    SingleElementInelasticProblem problem(
        norton_properties(), LoadingMode::axial_traction, time_step);
    const LoadPathResult result =
        solve_load_path(problem, step_count, time_step);

    const double axial_stress = average_axial_stress(problem.last_stress());
    const double equivalent_creep =
        average_equivalent_creep_strain(problem.committed_material());
    const double top_displacement = 0.5 * (result.state[10] + result.state[11]);
    const double top_displacement_spread =
        std::abs(result.state[10] - result.state[11]);
    const double stress_error =
        relative_error(axial_stress, moose_axial_stress);
    const double creep_error =
        relative_error(equivalent_creep, moose_equivalent_creep);
    const double displacement_error =
        relative_error(top_displacement, moose_top_displacement);
    const double maximum_trace = maximum_inelastic_trace(
        problem.committed_material(), LoadingMode::axial_traction);

    bool passed =
        check(result.converged, "Norton ten-step PETSc load path converged");
    passed = check(result.workspace_setups == 1,
                   "Norton load path creates one reusable PETSc workspace") &&
             passed;
    passed = check(result.solve_calls == step_count,
                   "Norton load path issues ten PETSc solve calls") &&
             passed;
    passed = check(result.callbacks_preserved_history,
                   "Norton residual and Jacobian callbacks do not commit "
                   "history") &&
             passed;
    passed = check(problem.committed_steps() == step_count,
                   "Norton history is committed after every accepted step") &&
             passed;
    passed = check(temperatures_are_committed(problem.committed_temperature()),
                   "Norton committed temperature remains 600 K") &&
             passed;
    passed = check(stress_error < moose_relative_tolerance,
                   "Norton final axial stress matches MOOSE within 0.1%") &&
             passed;
    passed = check(creep_error < moose_relative_tolerance,
                   "Norton final equivalent creep strain matches MOOSE within "
                   "0.1%") &&
             passed;
    passed = check(displacement_error < moose_relative_tolerance,
                   "Norton final top displacement matches MOOSE within 0.1%") &&
             passed;
    passed = check(top_displacement_spread < 1.0e-12,
                   "Norton top-side axial displacement is uniform") &&
             passed;
    passed = check(maximum_trace < 1.0e-12,
                   "Norton committed creep strain is trace-free") &&
             passed;

    std::cout << "m22_norton_axial_stress=" << axial_stress << '\n';
    std::cout << "m22_norton_equivalent_creep=" << equivalent_creep << '\n';
    std::cout << "m22_norton_top_displacement=" << top_displacement << '\n';
    std::cout << "m22_norton_stress_relative_error=" << stress_error << '\n';
    std::cout << "m22_norton_creep_relative_error=" << creep_error << '\n';
    std::cout << "m22_norton_displacement_relative_error=" << displacement_error
              << '\n';
    std::cout << "m22_norton_top_displacement_spread="
              << top_displacement_spread << '\n';
    std::cout << "m22_norton_maximum_creep_trace=" << maximum_trace << '\n';
    std::cout << "m22_norton_workspace_setups=" << result.workspace_setups
              << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M2.2 inelastic single-element solver tests\n");

        bool passed = true;
        passed = test_j2_moose_comparison() && passed;
        passed = test_norton_moose_comparison() && passed;
        if (!passed)
            return 1;

        std::cout << "[PASS] fuelsim M2.2 inelastic solver acceptance tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M2.2 inelastic solver tests raised: "
                  << error.what() << '\n';
        return 1;
    }
}
