#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/m2_problem.hpp"
#include "fuelsim/quad4_rz_transient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double scaled_error(double actual, double expected) {
    return std::abs(actual - expected) /
           (1.0 + std::max(std::abs(actual), std::abs(expected)));
}

fuelsim::ThermoelasticProperties
simple_thermoelastic(double young_modulus = 200.0) {
    return {
        0.0, 1.0, young_modulus, 0.25, 0.0, 600.0,
    };
}

fuelsim::TransientInelasticProperties
elastic_properties(double density = 10.0, double specific_heat = 20.0) {
    return {
        density,         specific_heat, fuelsim::InelasticBehavior::elastic,
        {0.0, 1.0, 1.0}, {1.0, 0.0},
    };
}

fuelsim::TransientInelasticProperties
creep_properties(double coefficient, double reference_stress, double exponent) {
    return {
        10.0,
        20.0,
        fuelsim::InelasticBehavior::norton_creep,
        {coefficient, reference_stress, exponent},
        {1.0, 0.0},
    };
}

fuelsim::TransientInelasticProperties
plastic_properties(double yield_stress, double hardening_modulus) {
    return {
        10.0,
        20.0,
        fuelsim::InelasticBehavior::j2_plasticity,
        {0.0, 1.0, 1.0},
        {yield_stress, hardening_modulus},
    };
}

double equivalent_stress(const fuelsim::AxisymmetricStress& stress) {
    const double mean =
        (stress.rr.value() + stress.zz.value() + stress.hoop.value()) / 3.0;
    const double rr = stress.rr.value() - mean;
    const double zz = stress.zz.value() - mean;
    const double hoop = stress.hoop.value() - mean;
    const double rz = stress.rz.value();
    return std::sqrt(1.5 * (rr * rr + zz * zz + hoop * hoop + 2.0 * rz * rz));
}

double inelastic_trace(const std::array<double, 4>& strain) {
    return strain[0] + strain[1] + strain[2];
}

bool same_state(const fuelsim::MaterialPointState& lhs,
                const fuelsim::MaterialPointState& rhs) {
    return lhs.plastic_strain == rhs.plastic_strain &&
           lhs.creep_strain == rhs.creep_strain &&
           lhs.equivalent_plastic_strain == rhs.equivalent_plastic_strain &&
           lhs.equivalent_creep_strain == rhs.equivalent_creep_strain;
}

bool test_j2_plasticity_material_point() {
    const fuelsim::IsotropicInelasticMaterial material(
        simple_thermoelastic(), plastic_properties(20.0, 40.0));
    const fuelsim::MaterialPointState committed{};
    const fuelsim::InelasticStressResponse response =
        material.response(0.2, -0.1, -0.1, 0.0, 600.0, 1.0, committed);
    const fuelsim::MaterialPointState state =
        fuelsim::IsotropicInelasticMaterial::state_values(response.trial_state);

    constexpr double expected_increment = 0.1;
    constexpr double expected_equivalent_stress = 24.0;
    bool passed =
        check(scaled_error(equivalent_stress(response.stress),
                           expected_equivalent_stress) < 1.0e-13,
              "J2 linear-hardening return reaches the analytic yield stress");
    passed = check(std::abs(state.equivalent_plastic_strain -
                            expected_increment) < 1.0e-14,
                   "J2 equivalent plastic increment matches the closed form") &&
             passed;
    passed = check(std::abs(state.plastic_strain[0] - 0.1) < 1.0e-14 &&
                       std::abs(state.plastic_strain[1] + 0.05) < 1.0e-14 &&
                       std::abs(state.plastic_strain[2] + 0.05) < 1.0e-14 &&
                       state.plastic_strain[3] == 0.0,
                   "J2 plastic flow follows the axisymmetric deviatoric "
                   "direction") &&
             passed;
    passed = check(std::abs(inelastic_trace(state.plastic_strain)) < 1.0e-14,
                   "J2 plastic strain is trace free") &&
             passed;
    passed = check(same_state(committed, fuelsim::MaterialPointState{}),
                   "J2 trial evaluation does not mutate committed state") &&
             passed;

    const fuelsim::InelasticStressResponse unloading = material.response(
        state.plastic_strain[0], state.plastic_strain[1],
        state.plastic_strain[2], state.plastic_strain[3], 600.0, 1.0, state);
    const fuelsim::MaterialPointState unloaded_state =
        fuelsim::IsotropicInelasticMaterial::state_values(
            unloading.trial_state);
    passed = check(same_state(state, unloaded_state),
                   "J2 unloading does not add plastic strain") &&
             passed;
    passed = check(equivalent_stress(unloading.stress) < 1.0e-13,
                   "J2 unloading to the plastic strain gives zero stress") &&
             passed;

    const adlite::Scalar active_rr = adlite::Scalar::independent(0.2, 0, 1);
    const fuelsim::InelasticStressResponse active =
        material.response(active_rr, -0.1, -0.1, 0.0, 600.0, 1.0, committed);
    constexpr double perturbation = 1.0e-6;
    const double plus = material
                            .response(0.2 + perturbation, -0.1, -0.1, 0.0,
                                      600.0, 1.0, committed)
                            .stress.rr.value();
    const double minus = material
                             .response(0.2 - perturbation, -0.1, -0.1, 0.0,
                                       600.0, 1.0, committed)
                             .stress.rr.value();
    const double finite_difference = (plus - minus) / (2.0 * perturbation);
    const double derivative_error =
        scaled_error(active.stress.rr.derivative(0), finite_difference);
    passed = check(derivative_error < 1.0e-9,
                   "active J2 AD tangent matches centered finite difference") &&
             passed;

    std::cout << "m22_j2_equivalent_stress="
              << equivalent_stress(response.stress) << '\n';
    std::cout << "m22_j2_equivalent_plastic_strain="
              << state.equivalent_plastic_strain << '\n';
    std::cout << "m22_j2_material_ad_scaled_error=" << derivative_error << '\n';
    return passed;
}

bool test_norton_creep_material_point() {
    constexpr double coefficient = 0.5;
    constexpr double reference_stress = 1.0;
    constexpr double time_step = 0.01;
    constexpr double shear_modulus = 80.0;
    constexpr double trial_stress = 48.0;
    const fuelsim::IsotropicInelasticMaterial material(
        simple_thermoelastic(),
        creep_properties(coefficient, reference_stress, 1.0));
    const fuelsim::MaterialPointState committed{};
    const fuelsim::InelasticStressResponse response =
        material.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState state =
        fuelsim::IsotropicInelasticMaterial::state_values(response.trial_state);

    const double expected_stress =
        trial_stress / (1.0 + 3.0 * shear_modulus * time_step * coefficient /
                                  reference_stress);
    const double expected_increment =
        (trial_stress - expected_stress) / (3.0 * shear_modulus);
    bool passed = check(scaled_error(equivalent_stress(response.stress),
                                     expected_stress) < 1.0e-13,
                        "linear Norton backward-Euler stress matches the "
                        "analytic root");
    passed = check(std::abs(state.equivalent_creep_strain -
                            expected_increment) < 1.0e-14,
                   "linear Norton equivalent creep increment matches the "
                   "analytic root") &&
             passed;
    passed = check(std::abs(inelastic_trace(state.creep_strain)) < 1.0e-14,
                   "Norton creep strain is trace free") &&
             passed;
    passed = check(same_state(committed, fuelsim::MaterialPointState{}),
                   "Norton trial evaluation does not mutate committed state") &&
             passed;

    const adlite::Scalar active_rr = adlite::Scalar::independent(0.2, 0, 1);
    const fuelsim::InelasticStressResponse active = material.response(
        active_rr, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    constexpr double perturbation = 1.0e-6;
    const double plus = material
                            .response(0.2 + perturbation, -0.1, -0.1, 0.0,
                                      600.0, time_step, committed)
                            .stress.rr.value();
    const double minus = material
                             .response(0.2 - perturbation, -0.1, -0.1, 0.0,
                                       600.0, time_step, committed)
                             .stress.rr.value();
    const double finite_difference = (plus - minus) / (2.0 * perturbation);
    const double derivative_error =
        scaled_error(active.stress.rr.derivative(0), finite_difference);
    passed =
        check(derivative_error < 1.0e-9,
              "active Norton AD tangent matches centered finite difference") &&
        passed;

    const adlite::Scalar zero_active_rr =
        adlite::Scalar::independent(0.0, 0, 1);
    const fuelsim::InelasticStressResponse zero_active = material.response(
        zero_active_rr, 0.0, 0.0, 0.0, 600.0, time_step, committed);
    constexpr double zero_perturbation = 1.0e-8;
    const double zero_plus = material
                                 .response(zero_perturbation, 0.0, 0.0, 0.0,
                                           600.0, time_step, committed)
                                 .stress.rr.value();
    const double zero_minus = material
                                  .response(-zero_perturbation, 0.0, 0.0, 0.0,
                                            600.0, time_step, committed)
                                  .stress.rr.value();
    const double zero_finite_difference =
        (zero_plus - zero_minus) / (2.0 * zero_perturbation);
    const double zero_derivative_error = scaled_error(
        zero_active.stress.rr.derivative(0), zero_finite_difference);
    passed = check(zero_derivative_error < 1.0e-9,
                   "linear Norton zero-stress AD tangent matches centered "
                   "finite difference") &&
             passed;

    const fuelsim::IsotropicInelasticMaterial cubic_material(
        simple_thermoelastic(), creep_properties(1.0e-4, 1.0, 3.0));
    const fuelsim::InelasticStressResponse cubic =
        cubic_material.response(0.2, -0.1, -0.1, 0.0, 600.0, 0.01, committed);
    const double cubic_stress = equivalent_stress(cubic.stress);
    const double cubic_residual =
        cubic_stress +
        3.0 * shear_modulus * 0.01 * 1.0e-4 * std::pow(cubic_stress, 3.0) -
        trial_stress;
    passed = check(std::abs(cubic_residual) < 1.0e-11 * (1.0 + trial_stress),
                   "nonlinear Norton safeguarded Newton satisfies its local "
                   "equation") &&
             passed;

    const fuelsim::ThermoelasticProperties small_scale_elastic = {
        0.0, 1.0, 1.0, 0.0, 0.0, 600.0,
    };
    const fuelsim::IsotropicInelasticMaterial extreme_scale_material(
        small_scale_elastic,
        creep_properties(1.0, std::numeric_limits<double>::max(), 2.0));
    const fuelsim::InelasticStressResponse extreme_scale =
        extreme_scale_material.response(1.0e-16, -0.5e-16, -0.5e-16, 0.0, 600.0,
                                        1.0, committed);
    const double extreme_scale_stress = equivalent_stress(extreme_scale.stress);
    passed =
        check(extreme_scale_stress > 0.0 &&
                  std::abs(extreme_scale_stress - 1.5e-16) / 1.5e-16 < 1.0e-12,
              "Norton extreme stress normalization does not erase a "
              "nonzero deviatoric stress") &&
        passed;

    constexpr double large_coefficient = 5.0e307;
    const fuelsim::IsotropicInelasticMaterial large_beta_material(
        small_scale_elastic, creep_properties(large_coefficient, 1.0, 2.0));
    const double one_direction = 1.0;
    const double negative_half_direction = -0.5;
    const adlite::Scalar large_beta_rr =
        adlite::Scalar::seeded(2.0 / 3.0, &one_direction, 1);
    const adlite::Scalar large_beta_zz =
        adlite::Scalar::seeded(-1.0 / 3.0, &negative_half_direction, 1);
    const adlite::Scalar large_beta_hoop =
        adlite::Scalar::seeded(-1.0 / 3.0, &negative_half_direction, 1);
    const fuelsim::InelasticStressResponse large_beta =
        large_beta_material.response(large_beta_rr, large_beta_zz,
                                     large_beta_hoop, 0.0, 600.0, 1.0,
                                     committed);
    const double large_beta_stress = equivalent_stress(large_beta.stress);
    const double beta = 1.5 * large_coefficient;
    const double expected_large_beta_derivative =
        1.0 / (1.0 + 2.0 * (beta * large_beta_stress));
    const double actual_large_beta_derivative =
        large_beta.stress.rr.derivative(0);
    passed = check(actual_large_beta_derivative > 0.0 &&
                       std::isfinite(actual_large_beta_derivative) &&
                       std::abs(actual_large_beta_derivative -
                                expected_large_beta_derivative) /
                               expected_large_beta_derivative <
                           1.0e-12,
                   "Norton large dimensionless coefficient keeps a finite "
                   "consistent tangent") &&
             passed;

    std::cout << "m22_norton_equivalent_stress="
              << equivalent_stress(response.stress) << '\n';
    std::cout << "m22_norton_equivalent_creep_strain="
              << state.equivalent_creep_strain << '\n';
    std::cout << "m22_norton_material_ad_scaled_error=" << derivative_error
              << '\n';
    std::cout << "m22_norton_zero_stress_ad_scaled_error="
              << zero_derivative_error << '\n';
    std::cout << "m22_norton_cubic_local_residual=" << cubic_residual << '\n';
    std::cout << "m22_norton_extreme_scale_equivalent_stress="
              << extreme_scale_stress << '\n';
    std::cout << "m22_norton_large_beta_tangent="
              << actual_large_beta_derivative << '\n';
    return passed;
}

fuelsim::Quad4RzGeometry test_geometry() {
    const fuelsim::Quad4Coordinates coordinates = {{
        {1.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.0},
        {1.0, 1.0},
    }};
    return fuelsim::make_quad4_rz_geometry(coordinates);
}

bool test_transient_element() {
    const fuelsim::Quad4RzGeometry geometry = test_geometry();
    const fuelsim::Quad4RzTransientKernel kernel(
        fuelsim::IsotropicInelasticMaterial(simple_thermoelastic(),
                                            elastic_properties()),
        100.0);
    const fuelsim::Quad4TemperatureHistory old_temperature = {
        600.0,
        600.0,
        600.0,
        600.0,
    };
    const fuelsim::Quad4MaterialHistory history{};

    fuelsim::LocalValues uniform_state = {
        605.0, 605.0, 605.0, 605.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    const fuelsim::LocalResidual balanced = kernel.residual(
        geometry, uniform_state, old_temperature, history, 10.0);
    bool passed = true;
    for (std::size_t row = 0; row < fuelsim::quad4_node_count; ++row)
        passed =
            check(std::abs(balanced[row]) < 1.0e-10,
                  "uniform transient heat source balances heat capacity") &&
            passed;

    const fuelsim::LocalValues state = {
        605.0, 607.0, 604.0, 603.0, -0.05, -0.10,
        -0.10, -0.05, 0.0,   0.0,   0.10,  0.10,
    };
    const fuelsim::LocalValues direction = {
        0.7, -0.4, 0.3, -0.6, 0.2, -0.4, 0.5, -0.1, -0.3, 0.6, -0.2, 0.4,
    };
    const fuelsim::LocalSystem system =
        kernel.linearize(geometry, state, old_temperature, history, 2.0);
    constexpr double perturbation = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += perturbation * direction[dof];
        minus[dof] -= perturbation * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        kernel.residual(geometry, plus, old_temperature, history, 2.0);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus, old_temperature, history, 2.0);

    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column) {
            ad_direction +=
                system.jacobian[row * fuelsim::quad4_local_dof_count + column] *
                direction[column];
        }
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
        maximum_jacobian_error =
            std::max(maximum_jacobian_error,
                     scaled_error(ad_direction, finite_difference));
    }
    passed = check(maximum_jacobian_error < 1.0e-7,
                   "transient Quad4 AD Jacobian matches centered finite "
                   "difference") &&
             passed;

    const fuelsim::LocalSystem slow =
        kernel.linearize(geometry, state, old_temperature, history, 5.0);
    double maximum_capacity_error = 0.0;
    constexpr double heat_capacity = 200.0;
    for (std::size_t row = 0; row < fuelsim::quad4_node_count; ++row) {
        for (std::size_t column = 0; column < fuelsim::quad4_node_count;
             ++column) {
            double mass = 0.0;
            for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
                mass += point.weighted_measure * point.shape[row] *
                        point.shape[column];
            }
            const double expected =
                heat_capacity * (1.0 / 2.0 - 1.0 / 5.0) * mass;
            const double actual =
                system.jacobian[row * fuelsim::quad4_local_dof_count + column] -
                slow.jacobian[row * fuelsim::quad4_local_dof_count + column];
            maximum_capacity_error = std::max(maximum_capacity_error,
                                              scaled_error(actual, expected));
        }
    }
    passed = check(maximum_capacity_error < 1.0e-13,
                   "transient Jacobian contains the consistent capacity "
                   "matrix divided by dt") &&
             passed;

    const fuelsim::Quad4MaterialHistory trial =
        kernel.trial_state_values(geometry, state, history, 2.0);
    for (std::size_t q = 0; q < trial.size(); ++q)
        passed = check(same_state(trial[q], history[q]),
                       "elastic transient element leaves material history "
                       "unchanged") &&
                 passed;

    std::cout << "m21_element_jacobian_maximum_scaled_error="
              << maximum_jacobian_error << '\n';
    std::cout << "m21_capacity_matrix_maximum_scaled_error="
              << maximum_capacity_error << '\n';
    return passed;
}

fuelsim::M2Parameters transaction_parameters() {
    const fuelsim::M1Parameters base = {
        1.0,
        1.1,
        1.2,
        1.0,
        1.02,
        1,
        1,
        1,
        simple_thermoelastic(),
        simple_thermoelastic(),
        100.0,
        600.0,
        600.0,
        1.0,
        0.01,
        1.0e4,
    };
    return {
        base,
        creep_properties(1.0e-4, 1.0, 1.0),
        plastic_properties(1.0, 10.0),
    };
}

bool test_problem_history_transaction() {
    fuelsim::M2Problem problem(transaction_parameters());
    const std::vector<double> initial_solution = problem.committed_solution();
    const fuelsim::MaterialPointState initial_history =
        problem.fuel_material_history(0)[0];

    problem.begin_time_step({1.0, 50.0});
    std::vector<double> trial_solution = initial_solution;
    const fuelsim::M1Problem& base = problem.base_problem();
    const fuelsim::DofMap& dofs = base.dof_map();
    for (std::size_t local_node = 0;
         local_node < base.fuel_mesh().nodes().size(); ++local_node) {
        const fuelsim::RzPoint& point = base.fuel_mesh().nodes()[local_node];
        const std::size_t global_node = base.fuel_global_node(local_node);
        trial_solution[dofs.radial_displacement(global_node)] = -0.05 * point.r;
        trial_solution[dofs.axial_displacement(global_node)] = 0.10 * point.z;
    }
    const fuelsim::LocalValues local_trial =
        problem.contribution_state(0, trial_solution);
    (void)problem.contribution_residual(0, local_trial);
    (void)problem.linearize_contribution(0, local_trial);

    bool passed =
        check(same_state(problem.fuel_material_history(0)[0], initial_history),
              "Newton residual and Jacobian callbacks do not mutate "
              "committed history");
    problem.commit_time_step(trial_solution);
    const double committed_creep =
        problem.summarize_fuel_history().maximum_equivalent_creep_strain;
    passed = check(problem.committed_time() == 1.0 &&
                       problem.committed_heat_source() == 50.0 &&
                       !problem.time_step_active(),
                   "accepted M2 time step commits time, load, and phase") &&
             passed;
    passed = check(committed_creep > 0.0,
                   "accepted M2 time step commits nonzero creep history") &&
             passed;
    passed = check(problem.committed_solution() == trial_solution,
                   "accepted M2 time step commits the converged nodal state") &&
             passed;

    problem.begin_time_step({2.0, 75.0});
    problem.rollback_time_step();
    passed = check(problem.committed_time() == 1.0 &&
                       problem.committed_heat_source() == 50.0 &&
                       problem.fuel_kernel().volumetric_heat_source() == 50.0,
                   "rollback restores the committed time and heat load") &&
             passed;
    passed =
        check(
            problem.summarize_fuel_history().maximum_equivalent_creep_strain ==
                    committed_creep &&
                problem.committed_solution() == trial_solution,
            "rollback preserves committed nodal and material history") &&
        passed;

    problem.begin_time_step({2.0, 75.0});
    std::vector<double> invalid_solution = trial_solution;
    invalid_solution[base.dof_map().temperature(base.fuel_global_node(0))] =
        -100.0;
    bool invalid_commit_threw = false;
    try {
        problem.commit_time_step(invalid_solution);
    } catch (const std::domain_error&) {
        invalid_commit_threw = true;
    }
    passed = check(invalid_commit_threw && problem.time_step_active(),
                   "M2 rejects a nonpositive nodal temperature before "
                   "changing committed state") &&
             passed;
    problem.rollback_time_step();
    passed =
        check(problem.committed_time() == 1.0 &&
                  problem.committed_solution() == trial_solution &&
                  problem.summarize_fuel_history()
                          .maximum_equivalent_creep_strain == committed_creep,
              "failed commit retains the previous complete transaction") &&
        passed;

    bool inactive_threw = false;
    try {
        (void)problem.contribution_residual(0, local_trial);
    } catch (const std::logic_error&) {
        inactive_threw = true;
    }
    passed =
        check(inactive_threw,
              "M2 residual evaluation without an active step is rejected") &&
        passed;

    std::cout << "m21_committed_time=" << problem.committed_time() << '\n';
    std::cout << "m22_committed_maximum_creep_strain=" << committed_creep
              << '\n';
    return passed;
}

} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    try {
        bool passed = true;
        passed = test_j2_plasticity_material_point() && passed;
        passed = test_norton_creep_material_point() && passed;
        passed = test_transient_element() && passed;
        passed = test_problem_history_transaction() && passed;
        if (!passed)
            return 1;
        std::cout << "[PASS] fuelsim M2 transient and inelastic core tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] unexpected exception: " << error.what() << '\n';
        return 1;
    }
}
