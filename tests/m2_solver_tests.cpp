#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/quad4_rz_transient.hpp"
#include "support/transient_fuel_cladding_problem.hpp"
#include "support/transient_fuel_cladding_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
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

class ImportedTransientHeatProblem final : public fuelsim::NonlinearProblem {
  public:
    explicit ImportedTransientHeatProblem(fuelsim::StructuredRzMesh mesh)
        : _mesh(std::move(mesh)), _dof_map(_mesh.nodes().size()),
          _kernel(fuelsim::IsotropicInelasticMaterial(
                      constant_material(2.0, 2.0e11, 0.3, 0.0),
                      elastic_transient_properties()),
                  3.0e6),
          _committed_temperature(_mesh.elements().size()),
          _committed_material(_mesh.elements().size()),
          _committed_solution(_dof_map.dof_count(), 0.0), _time_step(1.0) {
        _geometries.reserve(_mesh.elements().size());
        for (std::size_t element = 0; element < _mesh.elements().size();
             ++element) {
            fuelsim::Quad4Coordinates coordinates{};
            for (std::size_t node = 0; node < coordinates.size(); ++node)
                coordinates[node] =
                    _mesh.nodes().at(_mesh.elements()[element].nodes[node]);
            _geometries.push_back(fuelsim::make_quad4_rz_geometry(coordinates));
            _committed_temperature[element].fill(600.0);
        }
        for (std::size_t node = 0; node < _mesh.nodes().size(); ++node) {
            _committed_solution[_dof_map.temperature(node)] = 600.0;
            _dirichlet_conditions.push_back(
                {_dof_map.radial_displacement(node), 0.0});
            _dirichlet_conditions.push_back(
                {_dof_map.axial_displacement(node), 0.0});
        }
    }

    std::vector<double> initial_state() const {
        return _committed_solution;
    }

    void commit(const std::vector<double>& converged_state) {
        if (converged_state.size() != dof_count())
            throw std::invalid_argument(
                "ImportedTransientHeatProblem committed state size "
                "mismatch");

        std::vector<fuelsim::Quad4MaterialHistory> staged_material(
            _committed_material.size());
        std::vector<fuelsim::Quad4TemperatureHistory> staged_temperature(
            _committed_temperature.size());
        for (std::size_t element = 0; element < contribution_count();
             ++element) {
            const fuelsim::LocalValues local_state =
                contribution_state(element, converged_state);
            staged_material[element] = _kernel.trial_state_values(
                _geometries[element], local_state, _committed_material[element],
                _time_step);
            for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node)
                staged_temperature[element][node] = local_state[node];
        }
        _committed_material.swap(staged_material);
        _committed_temperature.swap(staged_temperature);
        _committed_solution = converged_state;
    }

    double committed_temperature(std::size_t node) const {
        return _committed_solution.at(_dof_map.temperature(node));
    }

    std::size_t node_count() const noexcept {
        return _mesh.nodes().size();
    }

    std::size_t dof_count() const noexcept override {
        return _dof_map.dof_count();
    }

    std::size_t contribution_count() const noexcept override {
        return _mesh.elements().size();
    }

    fuelsim::LocalDofs
    contribution_dofs(std::size_t contribution_index) const override {
        return _dof_map.element_dofs(_mesh.elements().at(contribution_index));
    }

    fuelsim::LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const fuelsim::LocalValues& state) const override {
        return _kernel.residual(_geometries.at(contribution_index), state,
                                _committed_temperature.at(contribution_index),
                                _committed_material.at(contribution_index),
                                _time_step);
    }

    fuelsim::LocalSystem
    linearize_contribution(std::size_t contribution_index,
                           const fuelsim::LocalValues& state) const override {
        return _kernel.linearize(_geometries.at(contribution_index), state,
                                 _committed_temperature.at(contribution_index),
                                 _committed_material.at(contribution_index),
                                 _time_step);
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
    fuelsim::StructuredRzMesh _mesh;
    fuelsim::DofMap _dof_map;
    std::vector<fuelsim::Quad4RzGeometry> _geometries;
    fuelsim::Quad4RzTransientKernel _kernel;
    std::vector<fuelsim::Quad4TemperatureHistory> _committed_temperature;
    std::vector<fuelsim::Quad4MaterialHistory> _committed_material;
    std::vector<double> _committed_solution;
    std::vector<fuelsim::DirichletCondition> _dirichlet_conditions;
    double _time_step;
};

bool test_moose_mesh_backward_euler_heat_source(const std::string& mesh_path) {
    constexpr std::size_t step_count = 10;
    constexpr double expected_temperature = 610.0;
    // verification/moose/m21_transient_heat_rz_out.csv at t = 10 s.
    constexpr double moose_snapshot_temperature = 610.0;

    const fuelsim::UnstructuredQuad4Mesh imported =
        fuelsim::ExodusMeshIo::read_quad4(mesh_path);
    fuelsim::StructuredRzMesh mesh =
        fuelsim::StructuredRzMesh::from_unstructured_block(
            imported, 0, {"left", "right", "bottom", "top"});
    ImportedTransientHeatProblem problem(std::move(mesh));
    fuelsim::PetscSequentialSolver solver;
    std::vector<double> state = problem.initial_state();
    std::size_t workspace_setups = 0;
    std::size_t solve_calls = 0;
    bool passed =
        check(problem.node_count() == 15 && problem.contribution_count() == 8,
              "M2.1 uses all 8 elements from the MOOSE Exodus mesh");

    for (std::size_t step = 0; step < step_count; ++step) {
        const fuelsim::SolveResult result = solver.solve(problem, state);
        passed = check(result.converged,
                       "M2.1 MOOSE-mesh backward-Euler step " +
                           std::to_string(step + 1) + " converged") &&
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
    for (std::size_t node = 0; node < problem.node_count(); ++node) {
        const double temperature = problem.committed_temperature(node);
        maximum_temperature_error =
            std::max(maximum_temperature_error,
                     std::abs(temperature - expected_temperature));
        average_temperature += temperature;
    }
    average_temperature /= static_cast<double>(problem.node_count());
    const double moose_snapshot_relative_error =
        std::abs(average_temperature - moose_snapshot_temperature) /
        moose_snapshot_temperature;

    passed = check(maximum_temperature_error < 1.0e-9,
                   "M2.1 MOOSE-mesh temperature reaches 610 K") &&
             passed;
    passed = check(moose_snapshot_relative_error < 1.0e-3,
                   "M2.1 temperature differs from the MOOSE "
                   "snapshot by less than 0.1%") &&
             passed;
    passed = check(workspace_setups == 1,
                   "ten backward-Euler steps create one PETSc workspace") &&
             passed;
    passed = check(solve_calls == step_count,
                   "one PETSc solve is issued per backward-Euler step") &&
             passed;

    std::cout << "m21_moose_mesh_average_temperature=" << average_temperature
              << '\n';
    std::cout << "m21_moose_mesh_relative_error="
              << moose_snapshot_relative_error << '\n';
    std::cout << "m21_moose_mesh_workspace_setups=" << workspace_setups << '\n';
    std::cout << "m21_moose_mesh_solve_calls=" << solve_calls << '\n';
    return passed;
}

fuelsim::TransientFuelCladdingParameters zero_source_m2_parameters() {
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

bool m2_histories_are_zero(
    const fuelsim::TransientFuelCladdingProblem& problem) {
    for (std::size_t element = 0;
         element < problem.steady_problem().fuel_element_count(); ++element) {
        if (!material_history_is_zero(problem.fuel_material_history(element)))
            return false;
    }
    for (std::size_t element = 0;
         element < problem.steady_problem().cladding_element_count();
         ++element) {
        if (!material_history_is_zero(
                problem.cladding_material_history(element)))
            return false;
    }
    return true;
}

bool test_m2_zero_source_history_and_interface() {
    fuelsim::TransientFuelCladdingProblem problem(zero_source_m2_parameters());
    const std::vector<double> initial_state = problem.committed_solution();
    const fuelsim::TransientTimeOptions time_options = {
        2.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0, 0.0,
    };

    const fuelsim::TransientFuelCladdingTimeStepper time_stepper;
    const fuelsim::TransientFuelCladdingResult result =
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
        problem.parameters().steady.axial_elements + 1;
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
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m2_solver_tests <m21_mesh.e>\n";
        return 2;
    }

    try {
        const std::string m21_mesh_path = argv[1];
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M2 transient solver acceptance tests\n");

        bool passed = true;
        passed =
            test_moose_mesh_backward_euler_heat_source(m21_mesh_path) && passed;
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
