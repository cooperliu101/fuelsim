#include "support/c3d8_common_tests.hpp"

namespace {
using namespace fuelsim::test::c3d8;

bool test_reduced_integration_inelastic_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
    fuelsim::Hex8LocalValues committed_state{}, state{};
    for (std::size_t node = 0; node < 8; ++node) {
        committed_state[node] = 300.0;
        state[node] = 302.0 + 0.35 * static_cast<double>(node);
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[8 + node] = 0.20 * point.x + 0.02 * point.y + 0.01 * point.z;
        state[16 + node] = 0.02 * point.x - 0.04 * point.y - 0.015 * point.z;
        state[24 + node] = 0.01 * point.x - 0.015 * point.y - 0.03 * point.z;
    }
    const fuelsim::CartesianMaterialHistory committed_material(1);
    bool passed = true;
    for (const std::array<bool, 2> branch : {std::array<bool, 2>{false, true}, {true, false}, {true, true}}) {
        const fuelsim::CartesianTestData data{
            fuelsim::IsotropicThermoelasticMaterial(inelastic_properties(branch[0], branch[1])),
            3.0,
            1.0,
            fuelsim::StrainFormulation::small,
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        fuelsim::Hex8LocalJacobian jacobian{};
        (void)
            fuelsim::compute_c3d8_transient(data, geometry, state, committed_state, committed_material, 1.0, &jacobian);

        std::array<double, 32> direction{};
        for (std::size_t dof = 0; dof < direction.size(); ++dof)
            direction[dof] = std::sin(0.31 * static_cast<double>(dof + 1));
        std::array<double, 3> block_errors{};
        for (std::size_t block = 0; block < block_errors.size(); ++block) {
            const bool thermal_columns = block != 2;
            const bool thermal_rows = block == 0;
            const double step = thermal_columns ? 1.0e-4 : 1.0e-7;
            fuelsim::Hex8LocalValues plus = state, minus = state;
            for (std::size_t column = 0; column < direction.size(); ++column) {
                if ((column < 8) != thermal_columns)
                    continue;
                plus[column] += step * direction[column];
                minus[column] -= step * direction[column];
            }
            const fuelsim::Hex8LocalResidual plus_residual =
                fuelsim::compute_c3d8_transient(data, geometry, plus, committed_state, committed_material, 1.0);
            const fuelsim::Hex8LocalResidual minus_residual =
                fuelsim::compute_c3d8_transient(data, geometry, minus, committed_state, committed_material, 1.0);
            double difference_squared = 0.0, reference_squared = 0.0;
            for (std::size_t row = 0; row < 32; ++row) {
                if ((row < 8) != thermal_rows)
                    continue;
                double analytic = 0.0;
                for (std::size_t column = 0; column < 32; ++column) {
                    if ((column < 8) != thermal_columns)
                        continue;
                    analytic += jacobian[row * 32 + column] * direction[column];
                }
                const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
                const double difference = analytic - numerical;
                difference_squared += difference * difference;
                reference_squared += numerical * numerical;
            }
            block_errors[block] = std::sqrt(difference_squared / reference_squared);
        }
        double thermal_displacement_maximum = 0.0;
        for (std::size_t row = 0; row < 8; ++row)
            for (std::size_t column = 8; column < 32; ++column)
                thermal_displacement_maximum =
                    std::max(thermal_displacement_maximum, std::abs(jacobian[row * 32 + column]));
        const fuelsim::CartesianMaterialHistory update =
            fuelsim::compute_c3d8_transient_update(data, geometry, state, committed_state, committed_material, 1.0);
        const fuelsim::CartesianMaterialPointState& point = update.front();
        const double plastic_trace = point.plastic_strain[0] + point.plastic_strain[1] + point.plastic_strain[2];
        const double creep_trace = point.creep_strain[0] + point.creep_strain[1] + point.creep_strain[2];
        std::cout << "hex8_c3d8rt_inelastic_ktt_relative_error=" << block_errors[0] << '\n'
                  << "hex8_c3d8rt_inelastic_kut_relative_error=" << block_errors[1] << '\n'
                  << "hex8_c3d8rt_inelastic_kuu_relative_error=" << block_errors[2] << '\n';
        passed =
            check(*std::max_element(block_errors.begin(), block_errors.end()) < 2.0e-6,
                "C3D8RT plastic, creep, and coupled Jacobian blocks match centered directional differences")
            && check(thermal_displacement_maximum == 0.0,
                "small-strain C3D8RT thermal residual has no displacement coupling")
            && check(update.size() == 1
                         && (branch[1] ? point.equivalent_plastic_strain > 0.0 : point.equivalent_plastic_strain == 0.0)
                         && (branch[0] ? point.equivalent_creep_strain > 0.0 : point.equivalent_creep_strain == 0.0)
                         && std::abs(plastic_trace) < 2.0e-15 && std::abs(creep_trace) < 2.0e-15,
                "C3D8RT commits exactly one traceless plastic-creep material point with the requested branches")
            && check(committed_material.front().equivalent_plastic_strain == 0.0
                         && committed_material.front().equivalent_creep_strain == 0.0,
                "C3D8RT residual, Jacobian, and trial update do not mutate committed history")
            && passed;
    }
    const fuelsim::CartesianTestData finite_data{fuelsim::IsotropicThermoelasticMaterial(properties()),
        3.0,
        1.0,
        fuelsim::StrainFormulation::finite,
        fuelsim::Hex8ElementFormulation::c3d8rt,
        300.0};
    fuelsim::Hex8LocalJacobian finite_jacobian{};
    const fuelsim::Hex8LocalResidual finite_residual = fuelsim::compute_c3d8_transient(finite_data,
        geometry,
        state,
        committed_state,
        committed_material,
        1.0,
        &finite_jacobian);
    const fuelsim::Hex8LocalResidual finite_residual_without_jacobian =
        fuelsim::compute_c3d8_transient(finite_data, geometry, state, committed_state, committed_material, 1.0);
    const fuelsim::CartesianMaterialHistory finite_update =
        fuelsim::compute_c3d8_transient_update(finite_data, geometry, state, committed_state, committed_material, 1.0);
    double thermal_displacement_maximum = 0.0;
    for (std::size_t row = 0; row < 8; ++row)
        for (std::size_t column = 8; column < 32; ++column)
            thermal_displacement_maximum =
                std::max(thermal_displacement_maximum, std::abs(finite_jacobian[row * 32 + column]));
    passed = check(std::all_of(finite_residual.begin(),
                       finite_residual.end(),
                       [](double value) { return std::isfinite(value); })
                       && thermal_displacement_maximum > 0.0 && finite_update.size() == 1,
                 "finite-strain C3D8RT evaluates current-geometry thermal-mechanical coupling and one material point")
             && check(finite_residual == finite_residual_without_jacobian,
                 "finite-strain C3D8RT Jacobian and residual-only paths produce identical residual values")
             && passed;
    const auto rejects_domain = [&](const fuelsim::Hex8LocalValues& trial, const fuelsim::Hex8LocalValues& committed) {
        try {
            (void)fuelsim::compute_c3d8_transient(finite_data, geometry, trial, committed, committed_material, 1.0);
            return false;
        } catch (const std::domain_error&) {
            return true;
        }
    };
    fuelsim::Hex8LocalValues invalid_current = state;
    for (std::size_t node = 0; node < 8; ++node)
        invalid_current[8 + node] = -2.0 * coordinates[node].x;
    fuelsim::Hex8LocalValues invalid_midpoint = committed_state;
    for (std::size_t node = 0; node < 8; ++node) {
        invalid_midpoint[8 + node] = -2.0 * coordinates[node].x;
        invalid_midpoint[16 + node] = -2.0 * coordinates[node].y;
    }
    fuelsim::Hex8LocalValues invalid_committed = committed_state, expanded_current = committed_state;
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t component = 0; component < 3; ++component) {
            const double coordinate = component == 0   ? coordinates[node].x
                                      : component == 1 ? coordinates[node].y
                                                       : coordinates[node].z;
            invalid_committed[8 * (component + 1) + node] = -2.0 * coordinate;
            expanded_current[8 * (component + 1) + node] = 2.0 * coordinate;
        }
    passed = check(rejects_domain(invalid_current, committed_state) && rejects_domain(invalid_midpoint, committed_state)
                       && rejects_domain(expanded_current, invalid_committed),
                 "finite-strain C3D8RT rejects nonpositive current, midpoint, and committed configurations")
             && passed;
    return passed;
}

bool test_reduced_integration_thermoelastic_capacity_gate() {
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
    fuelsim::Hex8LocalValues committed_state{}, state{};
    for (std::size_t node = 0; node < 8; ++node) {
        committed_state[node] = 300.0;
        state[node] = 310.0;
    }
    std::array<fuelsim::Hex8LocalResidual, 2> steady{}, transient{};
    std::array<fuelsim::Hex8LocalJacobian, 2> jacobian{};
    const std::array<fuelsim::StrainFormulation, 2> formulations = {fuelsim::StrainFormulation::small,
        fuelsim::StrainFormulation::finite};
    bool passed = true;
    for (std::size_t formulation = 0; formulation < formulations.size(); ++formulation) {
        const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(capacity_properties()),
            0.0,
            1.0,
            formulations[formulation],
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        steady[formulation] = fuelsim::compute_c3d8_thermoelastic(data, geometry, state);
        const fuelsim::elements::C3d8Input input{data.material,
            geometry,
            state,
            committed_state,
            nullptr,
            2.0,
            data.time,
            0.0,
            formulations[formulation],
            false,
            300.0};
        const auto without_capacity = fuelsim::elements::evaluate_c3d8rt(input, {true, true, false, false});
        passed = check(without_capacity.residual == steady[formulation],
                     "C3D8RT model request can disable capacity while retaining the supplied committed configuration")
                 && passed;
        passed = check(without_capacity.history.empty(), "Unrequested C3D8RT history stays empty") && passed;

        transient[formulation] =
            fuelsim::compute_c3d8_thermoelastic(data, geometry, state, &committed_state, 2.0, &jacobian[formulation]);
        for (std::size_t node = 0; node < 8; ++node)
            passed = check(transient[formulation][node] - steady[formulation][node] > 0.0,
                         "C3D8RT thermoelastic local interface adds heat capacity when a committed state is supplied")
                     && passed;
    }
    for (std::size_t node = 0; node < 8; ++node) {
        passed = check(near(transient[0][node] - steady[0][node], transient[1][node] - steady[1][node], 1.0e-13),
                     "small- and finite-strain C3D8RT use the same committed-state heat-capacity gate")
                 && passed;
        for (std::size_t column = 0; column < 8; ++column)
            passed = check(near(jacobian[0][node * 32 + column], jacobian[1][node * 32 + column], 1.0e-13),
                         "undeformed small- and finite-strain C3D8RT thermal Jacobians include the same capacity term")
                     && passed;
    }
    return passed;
}

bool test_reduced_integration_hourglass_energy() {
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = 400.0;
        const double mode = (node == 0 || node == 2 || node == 5 || node == 7) ? 1.0 : -1.0;
        state[8 + node] = 0.03 * mode + 0.01 * geometry.reduced_point.gradient[node][0];
        state[16 + node] = -0.02 * mode;
        state[24 + node] = 0.015 * mode;
    }
    const fuelsim::ThermoelasticProperties fixed = fuelsim::test::thermoelastic(0.0, 1.0, 200.0, 0.25, 0.0, 300.0, 0.0);
    const fuelsim::ThermoelasticProperties varying_initial =
        fuelsim::test::thermoelastic(0.0, 1.0, 100.0, 0.25, 0.0, 300.0, 1.0);
    bool passed = true;
    for (const fuelsim::StrainFormulation formulation :
        {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        const fuelsim::CartesianTestData fixed_data{fuelsim::IsotropicThermoelasticMaterial(fixed),
            0.0,
            0.0,
            formulation,
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        const fuelsim::CartesianTestData varying_data{fuelsim::IsotropicThermoelasticMaterial(varying_initial),
            0.0,
            0.0,
            formulation,
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        const fuelsim::Hex8LocalResidual fixed_residual =
            fuelsim::compute_c3d8_thermoelastic(fixed_data, geometry, state);
        const fuelsim::Hex8LocalResidual varying_residual =
            fuelsim::compute_c3d8_thermoelastic(varying_data, geometry, state);
        double maximum_error = 0.0, maximum_scale = 0.0;
        constexpr double step = 1.0e-7;
        for (std::size_t column = 8; column < state.size(); ++column) {
            fuelsim::Hex8LocalValues plus = state, minus = state;
            plus[column] += step;
            minus[column] -= step;
            const double plus_energy = fuelsim::elements::c3d8rt_hourglass_energy({fixed_data.material,
                                           geometry,
                                           plus,
                                           plus,
                                           nullptr,
                                           0.0,
                                           fixed_data.time,
                                           fixed_data.volumetric_heat_source,
                                           fixed_data.strain_formulation,
                                           false,
                                           fixed_data.initial_temperature})
                                       - fuelsim::elements::c3d8rt_hourglass_energy({varying_data.material,
                                           geometry,
                                           plus,
                                           plus,
                                           nullptr,
                                           0.0,
                                           varying_data.time,
                                           varying_data.volumetric_heat_source,
                                           varying_data.strain_formulation,
                                           false,
                                           varying_data.initial_temperature});
            const double minus_energy = fuelsim::elements::c3d8rt_hourglass_energy({fixed_data.material,
                                            geometry,
                                            minus,
                                            minus,
                                            nullptr,
                                            0.0,
                                            fixed_data.time,
                                            fixed_data.volumetric_heat_source,
                                            fixed_data.strain_formulation,
                                            false,
                                            fixed_data.initial_temperature})
                                        - fuelsim::elements::c3d8rt_hourglass_energy({varying_data.material,
                                            geometry,
                                            minus,
                                            minus,
                                            nullptr,
                                            0.0,
                                            varying_data.time,
                                            varying_data.volumetric_heat_source,
                                            varying_data.strain_formulation,
                                            false,
                                            varying_data.initial_temperature});
            const double numerical = (plus_energy - minus_energy) / (2.0 * step);
            const double analytic = fixed_residual[column] - varying_residual[column];
            maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
            maximum_scale = std::max({maximum_scale, std::abs(analytic), std::abs(numerical)});
        }
        const double relative_error = maximum_error / maximum_scale;
        std::cout << "hex8_c3d8rt_hourglass_energy_gradient_relative_error=" << relative_error << '\n';
        passed = check(relative_error < 2.0e-8,
                     "C3D8RT mechanical hourglass energy gradient equals the residual hourglass force")
                 && passed;
    }
    return passed;
}

int run_c3d8rt_tests() {
    bool passed = test_history_geometry(true);
    passed = test_reduced_integration_inelastic_jacobian() && passed;
    passed = test_reduced_integration_thermoelastic_capacity_gate() && passed;
    passed = test_reduced_integration_hourglass_energy() && passed;
    return passed ? 0 : 1;
}
} // namespace

int main() {
    return run_c3d8rt_tests();
}
