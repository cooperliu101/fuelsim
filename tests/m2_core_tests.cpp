#include "fuelsim/core/rz_quad4.hpp"
#include "fuelsim/core/transient_problem.hpp"
#include "material.hpp"
#include "support/material_factory.hpp"
#include "support/mesh_fixture.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace fuelsim {}

#include <utility>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double scaled_error(double actual, double expected) {
    return std::abs(actual - expected) / (1.0 + std::max(std::abs(actual), std::abs(expected)));
}

fuelsim::ThermoelasticProperties
simple_thermoelastic(double young_modulus = 200.0, double density = 10.0, double specific_heat = 20.0) {
    return fuelsim::test::thermoelastic(0.0,
        1.0,
        young_modulus,
        0.25,
        0.0,
        600.0,
        0.0,
        0.0,
        0.0,
        density,
        specific_heat);
}

fuelsim::ThermoelasticProperties elastic_properties(double density = 10.0, double specific_heat = 20.0) {
    return simple_thermoelastic(200.0, density, specific_heat);
}

fuelsim::ThermoelasticProperties creep_properties(double coefficient, double reference_stress, double exponent) {
    return fuelsim::test::with_norton(simple_thermoelastic(), coefficient, reference_stress, exponent);
}

fuelsim::ThermoelasticProperties plastic_properties(double yield_stress, double hardening_modulus) {
    return fuelsim::test::with_plasticity(simple_thermoelastic(), yield_stress, hardening_modulus);
}

fuelsim::ThermoelasticProperties coupled_properties(fuelsim::ThermoelasticProperties material,
    double coefficient,
    double reference_stress,
    double exponent,
    double yield_stress,
    double hardening_modulus) {
    return fuelsim::test::with_plasticity(
        fuelsim::test::with_norton(std::move(material), coefficient, reference_stress, exponent),
        yield_stress,
        hardening_modulus);
}

fuelsim::ThermoelasticProperties coupled_properties(double coefficient,
    double reference_stress,
    double exponent,
    double yield_stress,
    double hardening_modulus) {
    return coupled_properties(simple_thermoelastic(),
        coefficient,
        reference_stress,
        exponent,
        yield_stress,
        hardening_modulus);
}

double equivalent_stress(const fuelsim::AxisymmetricStress& stress) {
    const double mean = (stress.rr.value() + stress.zz.value() + stress.hoop.value()) / 3.0;
    const double rr = stress.rr.value() - mean;
    const double zz = stress.zz.value() - mean;
    const double hoop = stress.hoop.value() - mean;
    const double rz = stress.rz.value();
    return std::sqrt(1.5 * (rr * rr + zz * zz + hoop * hoop + 2.0 * rz * rz));
}

double inelastic_trace(const std::array<double, 4>& strain) {
    return strain[0] + strain[1] + strain[2];
}

bool same_state(const fuelsim::MaterialPointState& lhs, const fuelsim::MaterialPointState& rhs) {
    return lhs.elastic_strain == rhs.elastic_strain && lhs.plastic_strain == rhs.plastic_strain
           && lhs.creep_strain == rhs.creep_strain && lhs.equivalent_plastic_strain == rhs.equivalent_plastic_strain
           && lhs.equivalent_creep_strain == rhs.equivalent_creep_strain;
}

bool same_inelastic_state(const fuelsim::MaterialPointState& lhs, const fuelsim::MaterialPointState& rhs) {
    return lhs.plastic_strain == rhs.plastic_strain && lhs.creep_strain == rhs.creep_strain
           && lhs.equivalent_plastic_strain == rhs.equivalent_plastic_strain
           && lhs.equivalent_creep_strain == rhs.equivalent_creep_strain;
}

bool test_builtin_material_parameter_order() {
    const fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>();
    functions->name = "ordered_builtin_test";
    functions->thermal = registry.bind_thermal("inverse_temperature_thermophysical",
        {{"specific_heat", 20.0},
            {"conductivity_constant", 2.0},
            {"density", 10.0},
            {"conductivity_inverse_temperature", 300.0}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"poisson_ratio_temperature_coefficient", 1.0e-3},
            {"reference_temperature", 100.0},
            {"young_modulus", 1000.0},
            {"young_modulus_temperature_coefficient", -1.0},
            {"poisson_ratio", 0.2}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "linear_temperature_isotropic_thermal_expansion",
        {{"reference_temperature", 100.0},
            {"thermal_expansion_temperature_coefficient", 1.0e-3},
            {"thermal_expansion", 1.0e-2}}));
    const fuelsim::IsotropicThermoelasticMaterial material({std::move(functions), 1000.0});
    const adlite::Scalar temperature = adlite::Scalar::independent(110.0, 0, 1);
    const adlite::Scalar conductivity = material.conductivity(temperature);
    const adlite::Scalar heat_capacity = material.heat_capacity(temperature);
    const fuelsim::ActiveThermoelasticProperties elasticity = material.active_properties(temperature);
    const fuelsim::AxisymmetricStrain eigenstrain = material.eigenstrain_rz(temperature);
    constexpr double young_modulus = 990.0, poisson_ratio = 0.21;
    const double expected_shear = young_modulus / (2.0 * (1.0 + poisson_ratio));
    const double expected_lame = young_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    constexpr double perturbation = 1.0e-5;
    const fuelsim::ActiveThermoelasticProperties plus = material.active_properties(110.0 + perturbation);
    const fuelsim::ActiveThermoelasticProperties minus = material.active_properties(110.0 - perturbation);
    bool passed = check(scaled_error(conductivity.value(), 300.0 / 110.0 + 2.0) < 1.0e-14
                            && scaled_error(conductivity.derivative(0), -300.0 / (110.0 * 110.0)) < 1.0e-14
                            && heat_capacity.value() == 200.0 && !heat_capacity.is_active(),
        "built-in thermal functions use registry schema order after named binding");
    passed = check(scaled_error(elasticity.shear_modulus.value(), expected_shear) < 1.0e-14
                       && scaled_error(elasticity.lame_lambda.value(), expected_lame) < 1.0e-14
                       && scaled_error(elasticity.shear_modulus.derivative(0),
                              (plus.shear_modulus.value() - minus.shear_modulus.value()) / (2.0 * perturbation))
                              < 1.0e-8
                       && scaled_error(elasticity.lame_lambda.derivative(0),
                              (plus.lame_lambda.value() - minus.lame_lambda.value()) / (2.0 * perturbation))
                              < 1.0e-8,
                 "built-in elasticity uses registry schema order after named binding")
             && passed;
    passed = check(scaled_error(eigenstrain.rr.value(), 0.2) < 1.0e-14
                       && scaled_error(eigenstrain.rr.derivative(0), 0.03) < 1.0e-14
                       && eigenstrain.rr.value() == eigenstrain.zz.value()
                       && eigenstrain.rr.value() == eigenstrain.hoop.value() && eigenstrain.rz.value() == 0.0,
                 "built-in eigenstrain uses registry schema order after named binding")
             && passed;
    auto linear_functions = std::make_shared<fuelsim::MaterialFunctionSet>(material.functions());
    linear_functions->name = "linear_temperature_thermophysical_test";
    linear_functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
        {{"specific_heat_temperature_coefficient", 0.3},
            {"density", 10.0},
            {"reference_temperature", 100.0},
            {"conductivity", 4.0},
            {"specific_heat", 20.0},
            {"density_temperature_coefficient", 0.2},
            {"conductivity_temperature_coefficient", 0.1}});
    const fuelsim::IsotropicThermoelasticMaterial linear_material({std::move(linear_functions), 1000.0});
    const adlite::Scalar linear_conductivity = linear_material.conductivity(temperature);
    const adlite::Scalar linear_capacity = linear_material.heat_capacity(temperature);
    passed =
        check(linear_conductivity.value() == 5.0 && linear_conductivity.derivative(0) == 0.1
                  && linear_capacity.value() == 276.0 && scaled_error(linear_capacity.derivative(0), 8.2) < 1.0e-14,
            "linear-temperature conductivity, density, and specific heat preserve exact temperature derivatives")
        && passed;
    return passed;
}

fuelsim::ThermalPropertyEvaluator custom_thermal_properties(const fuelsim::MaterialParameters& named) {
    const double conductivity_offset = named.value("conductivity_offset");
    const double conductivity_slope = named.value("conductivity_slope");
    const double density = named.value("density");
    const double specific_heat_offset = named.value("specific_heat_offset");
    const double specific_heat_slope = named.value("specific_heat_slope");
    return [=](const fuelsim::ThermoelasticFunctionInput& input, fuelsim::ThermalPropertyOutput& output) {
        output.conductivity = conductivity_offset + conductivity_slope * input.temperature;
        output.density = density;
        output.specific_heat = specific_heat_offset + specific_heat_slope * input.temperature;
    };
}

fuelsim::ElasticPropertyEvaluator custom_elastic_properties(const fuelsim::MaterialParameters& named) {
    const double young_modulus = named.value("young_modulus");
    const double young_modulus_temperature_coefficient = named.value("young_modulus_temperature_coefficient");
    const double reference_temperature = named.value("reference_temperature");
    const double poisson_ratio = named.value("poisson_ratio");
    return [=](const fuelsim::ThermoelasticFunctionInput& input, fuelsim::ElasticPropertyOutput& output) {
        output.young_modulus =
            young_modulus + young_modulus_temperature_coefficient * (input.temperature - reference_temperature);
        output.poisson_ratio = poisson_ratio;
    };
}

fuelsim::EigenstrainEvaluator custom_eigenstrain(const fuelsim::MaterialParameters& named) {
    const double coefficient = named.value("coefficient");
    const double reference_temperature = named.value("reference_temperature");
    return [=](const fuelsim::ThermoelasticFunctionInput& input, fuelsim::SymmetricTensor3& output) {
        const adlite::Scalar value = coefficient * (input.temperature - reference_temperature);
        output = {value, value, value, 0.0, 0.0, 0.0};
    };
}

fuelsim::CreepRateEvaluator custom_creep_rate(const fuelsim::MaterialParameters& named) {
    const double coefficient = named.value("coefficient");
    const double reference_stress = named.value("reference_stress");
    return [=](const fuelsim::CreepRateInput& input) {
        return coefficient * input.equivalent_stress / reference_stress;
    };
}

fuelsim::PlasticFlowStressEvaluator custom_flow_stress(const fuelsim::MaterialParameters& named) {
    const double yield_stress = named.value("yield_stress");
    const double hardening_modulus = named.value("hardening_modulus");
    return [=](const fuelsim::PlasticFlowStressInput& input) {
        return yield_stress + hardening_modulus * input.equivalent_plastic_strain;
    };
}

bool test_registered_material_functions() {
    fuelsim::MaterialFunctionRegistry registry;
    registry.add_thermal("custom_thermal",
        {{"conductivity_offset", "W/(m*K)"},
            {"conductivity_slope", "W/(m*K^2)"},
            {"density", "kg/m^3"},
            {"specific_heat_offset", "J/(kg*K)"},
            {"specific_heat_slope", "J/(kg*K^2)"}},
        &custom_thermal_properties);
    registry.add_elasticity("custom_elastic",
        {{"young_modulus", "Pa"},
            {"poisson_ratio", "1"},
            {"reference_temperature", "K"},
            {"young_modulus_temperature_coefficient", "Pa/K"}},
        &custom_elastic_properties);
    registry.add_eigenstrain("custom_eigenstrain",
        {{"coefficient", "1/K"}, {"reference_temperature", "K"}},
        &custom_eigenstrain);
    registry.add_creep("custom_creep", {{"coefficient", "1/s"}, {"reference_stress", "Pa"}}, &custom_creep_rate);
    registry.add_plasticity("custom_plasticity",
        {{"yield_stress", "Pa"}, {"hardening_modulus", "Pa"}},
        &custom_flow_stress);
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>();
    functions->name = "registered_test";
    functions->thermal = registry.bind_thermal("custom_thermal",
        {{"specific_heat_slope", 0.1},
            {"density", 10.0},
            {"conductivity_slope", 0.01},
            {"specific_heat_offset", 20.0},
            {"conductivity_offset", 1.0}});
    functions->elasticity = registry.bind_elasticity("custom_elastic",
        {{"poisson_ratio", 0.25},
            {"young_modulus_temperature_coefficient", -0.1},
            {"reference_temperature", 600.0},
            {"young_modulus", 1000.0}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "custom_eigenstrain",
        {{"reference_temperature", 600.0}, {"coefficient", 1.0e-5}}));
    functions->creep = registry.bind_creep("custom_creep", {{"reference_stress", 10.0}, {"coefficient", 1.0e-3}});
    functions->plasticity =
        registry.bind_plasticity("custom_plasticity", {{"hardening_modulus", 50.0}, {"yield_stress", 10.0}});
    const fuelsim::IsotropicThermoelasticMaterial material({functions, 1000.0});
    const fuelsim::MaterialPointState committed{};
    const adlite::Scalar active_strain = adlite::Scalar::independent(0.05, 0, 1);
    const fuelsim::InelasticStressResponse active =
        material.response(active_strain, -0.01, -0.01, 0.005, 650.0, 2.0, committed, {4.0, 0.004, 0.0, 0.01});
    constexpr double perturbation = 1.0e-7;
    const double plus =
        material.response(0.05 + perturbation, -0.01, -0.01, 0.005, 650.0, 2.0, committed, {4.0, 0.004, 0.0, 0.01})
            .stress.rr.value();
    const double minus =
        material.response(0.05 - perturbation, -0.01, -0.01, 0.005, 650.0, 2.0, committed, {4.0, 0.004, 0.0, 0.01})
            .stress.rr.value();
    const double centered = (plus - minus) / (2.0 * perturbation);
    const adlite::Scalar active_temperature = adlite::Scalar::independent(650.0, 0, 1);
    const adlite::Scalar conductivity = material.conductivity(active_temperature, {4.0, 0.004, 0.0, 0.01});
    const adlite::Scalar heat_capacity = material.heat_capacity(active_temperature, {4.0, 0.004, 0.0, 0.01});
    bool passed = check(functions->thermal.parameters.value("density") == 10.0
                            && functions->thermal.parameters.value("conductivity_offset") == 1.0,
        "registered material parameters are bound by name and stored in schema order");
    passed = check(scaled_error(active.stress.rr.derivative(0), centered) < 1.0e-7,
                 "registered creep-plastic material AD tangent matches a centered difference")
             && passed;
    passed =
        check(std::abs(conductivity.value() - 7.5) < 1.0e-14 && std::abs(conductivity.derivative(0) - 0.01) < 1.0e-14
                  && std::abs(heat_capacity.value() - 850.0) < 1.0e-12
                  && std::abs(heat_capacity.derivative(0) - 1.0) < 1.0e-14,
            "registered thermal properties preserve temperature derivatives")
        && passed;
    passed = check(same_state(committed, fuelsim::MaterialPointState{}),
                 "registered material evaluation does not mutate committed history")
             && passed;
    fuelsim::MaterialFunctionSet changed_functions = *functions;
    changed_functions.creep =
        registry.bind_creep("custom_creep", {{"reference_stress", 10.0}, {"coefficient", 1.1e-3}});
    passed = check(changed_functions.signature() != functions->signature(),
                 "registered material signature includes named parameter values")
             && passed;
    bool missing_parameter_threw = false;
    try {
        (void)registry.bind_creep("custom_creep", {{"coefficient", 1.0}});
    } catch (const std::invalid_argument& error) {
        missing_parameter_threw = std::string(error.what()).find("reference_stress") != std::string::npos;
    }
    passed = check(missing_parameter_threw, "registered material functions reject a missing named parameter") && passed;
    bool empty_unit_threw = false;
    try {
        registry.add_creep("invalid_unit", {{"coefficient", ""}}, &custom_creep_rate);
    } catch (const std::invalid_argument& error) {
        empty_unit_threw = std::string(error.what()).find("empty SI unit") != std::string::npos;
    }
    passed = check(empty_unit_threw, "registered material parameters require an SI unit description") && passed;
    return passed;
}

bool test_objective_incremental_history_rotation() {
    const fuelsim::IsotropicThermoelasticMaterial material(elastic_properties());
    fuelsim::MaterialPointState committed;
    committed.elastic_strain = {0.020, -0.012, -0.008, 0.006};
    committed.plastic_strain = {0.030, -0.018, -0.012, 0.004};
    committed.creep_strain = {-0.016, 0.010, 0.006, -0.003};
    committed.equivalent_plastic_strain = 0.041;
    committed.equivalent_creep_strain = 0.027;
    constexpr double angle = 0.37;
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const fuelsim::AxisymmetricRotation rotation = {cosine, -sine, sine, cosine, 1.0};
    const std::array<double, 4> total = {
        committed.elastic_strain[0] + committed.plastic_strain[0] + committed.creep_strain[0],
        committed.elastic_strain[1] + committed.plastic_strain[1] + committed.creep_strain[1],
        committed.elastic_strain[2] + committed.plastic_strain[2] + committed.creep_strain[2],
        committed.elastic_strain[3] + committed.plastic_strain[3] + committed.creep_strain[3],
    };
    const fuelsim::InelasticStressResponse unrotated =
        material.response(total[0], total[1], total[2], total[3], 600.0, 1.0, committed);
    const fuelsim::InelasticStressResponse rotated =
        material.incremental_response(0.0, 0.0, 0.0, 0.0, rotation, 600.0, 600.0, 1.0, committed);
    const fuelsim::AxisymmetricStress expected_stress = fuelsim::rotate_axisymmetric_tensor(unrotated.stress, rotation);
    const auto rotate_values = [&](const std::array<double, 4>& values) {
        const fuelsim::AxisymmetricStress tensor = {values[0], values[1], values[2], values[3]};
        const fuelsim::AxisymmetricStress value = fuelsim::rotate_axisymmetric_tensor(tensor, rotation);
        return std::array<double, 4>{value.rr.value(), value.zz.value(), value.hoop.value(), value.rz.value()};
    };
    const fuelsim::MaterialPointState state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(rotated.trial_state);
    const std::array<double, 4> expected_elastic = rotate_values(committed.elastic_strain);
    const std::array<double, 4> expected_plastic = rotate_values(committed.plastic_strain);
    const std::array<double, 4> expected_creep = rotate_values(committed.creep_strain);
    double maximum_error = 0.0;
    for (std::size_t component = 0; component < 4; ++component) {
        maximum_error = std::max({maximum_error,
            std::abs(state.elastic_strain[component] - expected_elastic[component]),
            std::abs(state.plastic_strain[component] - expected_plastic[component]),
            std::abs(state.creep_strain[component] - expected_creep[component])});
    }
    maximum_error = std::max({maximum_error,
        scaled_error(rotated.stress.rr.value(), expected_stress.rr.value()),
        scaled_error(rotated.stress.zz.value(), expected_stress.zz.value()),
        scaled_error(rotated.stress.hoop.value(), expected_stress.hoop.value()),
        scaled_error(rotated.stress.rz.value(), expected_stress.rz.value())});
    bool passed = check(maximum_error < 1.0e-14,
        "incremental finite strain objectively rotates stress and all "
        "tensor histories");
    passed = check(std::abs(inelastic_trace(state.plastic_strain)) < 1.0e-14
                       && std::abs(inelastic_trace(state.creep_strain)) < 1.0e-14
                       && state.equivalent_plastic_strain == committed.equivalent_plastic_strain
                       && state.equivalent_creep_strain == committed.equivalent_creep_strain,
                 "objective rotation preserves trace-free histories and "
                 "equivalent scalars")
             && passed;
    std::cout << "m41_objective_history_rotation_maximum_error=" << maximum_error << '\n';
    return passed;
}
enum class TestInelasticBehavior { elastic, creep, plastic, coupled };

double temperature_tangent_error(TestInelasticBehavior behavior) {
    fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 1.0, 200.0, 0.25, 1.0e-5, 600.0, -0.08, 1.0e-5, 2.0e-8, 10.0, 20.0);
    if (behavior == TestInelasticBehavior::creep || behavior == TestInelasticBehavior::coupled)
        properties = fuelsim::test::with_norton(std::move(properties), 0.5, 1.0, 1.0, 600.0, 1.0e-4, 2.0e-3, 1.0e-4);
    if (behavior == TestInelasticBehavior::plastic || behavior == TestInelasticBehavior::coupled)
        properties = fuelsim::test::with_plasticity(std::move(properties), 20.0, 40.0, 600.0, -1.0e-2, -2.0e-2);
    const fuelsim::IsotropicThermoelasticMaterial material(std::move(properties));
    const fuelsim::MaterialPointState committed{};
    constexpr double temperature = 610.0;
    constexpr double time_step = 0.01;
    const adlite::Scalar active_temperature = adlite::Scalar::independent(temperature, 0, 1);
    const fuelsim::InelasticStressResponse active =
        material.response(0.2, -0.1, -0.1, 0.02, active_temperature, time_step, committed);
    constexpr double perturbation = 1.0e-4;
    const double plus =
        material.response(0.2, -0.1, -0.1, 0.02, temperature + perturbation, time_step, committed).stress.rr.value();
    const double minus =
        material.response(0.2, -0.1, -0.1, 0.02, temperature - perturbation, time_step, committed).stress.rr.value();
    return scaled_error(active.stress.rr.derivative(0), (plus - minus) / (2.0 * perturbation));
}

bool test_temperature_active_inelastic_properties() {
    const double elastic_error = temperature_tangent_error(TestInelasticBehavior::elastic);
    const double plastic_error = temperature_tangent_error(TestInelasticBehavior::plastic);
    const double creep_error = temperature_tangent_error(TestInelasticBehavior::creep);
    const double coupled_error = temperature_tangent_error(TestInelasticBehavior::coupled);
    const double maximum_error = std::max({elastic_error, plastic_error, creep_error, coupled_error});
    std::cout << "m40_active_elastic_temperature_tangent_error=" << elastic_error << '\n';
    std::cout << "m40_active_plastic_temperature_tangent_error=" << plastic_error << '\n';
    std::cout << "m40_active_creep_temperature_tangent_error=" << creep_error << '\n';
    std::cout << "m40_active_coupled_temperature_tangent_error=" << coupled_error << '\n';
    std::cout << "m40_active_material_temperature_tangent_error=" << maximum_error << '\n';
    return check(maximum_error < 1.0e-7,
        "elastic, plastic, creep, and coupled temperature-active "
        "tangents match centered differences");
}

bool test_j2_plasticity_material_point() {
    const fuelsim::IsotropicThermoelasticMaterial material(plastic_properties(20.0, 40.0));
    const fuelsim::MaterialPointState committed{};
    const fuelsim::InelasticStressResponse response = material.response(0.2, -0.1, -0.1, 0.0, 600.0, 1.0, committed);
    const fuelsim::MaterialPointState state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(response.trial_state);
    constexpr double expected_increment = 0.1;
    constexpr double expected_equivalent_stress = 24.0;
    bool passed = check(scaled_error(equivalent_stress(response.stress), expected_equivalent_stress) < 1.0e-13,
        "J2 linear-hardening return reaches the analytic yield stress");
    passed = check(std::abs(state.equivalent_plastic_strain - expected_increment) < 1.0e-14,
                 "J2 equivalent plastic increment matches the closed form")
             && passed;
    passed =
        check(std::abs(state.plastic_strain[0] - 0.1) < 1.0e-14 && std::abs(state.plastic_strain[1] + 0.05) < 1.0e-14
                  && std::abs(state.plastic_strain[2] + 0.05) < 1.0e-14 && state.plastic_strain[3] == 0.0,
            "J2 plastic flow follows the axisymmetric deviatoric "
            "direction")
        && passed;
    passed =
        check(std::abs(inelastic_trace(state.plastic_strain)) < 1.0e-14, "J2 plastic strain is trace free") && passed;
    passed = check(same_state(committed, fuelsim::MaterialPointState{}),
                 "J2 trial evaluation does not mutate committed state")
             && passed;
    const fuelsim::InelasticStressResponse unloading = material.response(state.plastic_strain[0],
        state.plastic_strain[1],
        state.plastic_strain[2],
        state.plastic_strain[3],
        600.0,
        1.0,
        state);
    const fuelsim::MaterialPointState unloaded_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(unloading.trial_state);
    passed = check(same_inelastic_state(state, unloaded_state), "J2 unloading does not add plastic strain") && passed;
    passed =
        check(equivalent_stress(unloading.stress) < 1.0e-13, "J2 unloading to the plastic strain gives zero stress")
        && passed;
    const adlite::Scalar active_rr = adlite::Scalar::independent(0.2, 0, 1);
    const fuelsim::InelasticStressResponse active =
        material.response(active_rr, -0.1, -0.1, 0.0, 600.0, 1.0, committed);
    constexpr double perturbation = 1.0e-6;
    const double plus = material.response(0.2 + perturbation, -0.1, -0.1, 0.0, 600.0, 1.0, committed).stress.rr.value();
    const double minus =
        material.response(0.2 - perturbation, -0.1, -0.1, 0.0, 600.0, 1.0, committed).stress.rr.value();
    const double finite_difference = (plus - minus) / (2.0 * perturbation);
    const double derivative_error = scaled_error(active.stress.rr.derivative(0), finite_difference);
    passed = check(derivative_error < 1.0e-9, "active J2 AD tangent matches centered finite difference") && passed;
    std::cout << "m22_j2_equivalent_stress=" << equivalent_stress(response.stress) << '\n';
    std::cout << "m22_j2_equivalent_plastic_strain=" << state.equivalent_plastic_strain << '\n';
    std::cout << "m22_j2_material_ad_scaled_error=" << derivative_error << '\n';
    return passed;
}

bool test_norton_creep_material_point() {
    constexpr double coefficient = 0.5;
    constexpr double reference_stress = 1.0;
    constexpr double time_step = 0.01;
    constexpr double shear_modulus = 80.0;
    constexpr double trial_stress = 48.0;
    const fuelsim::IsotropicThermoelasticMaterial material(creep_properties(coefficient, reference_stress, 1.0));
    const fuelsim::MaterialPointState committed{};
    const fuelsim::InelasticStressResponse response =
        material.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(response.trial_state);
    const double expected_stress =
        trial_stress / (1.0 + 3.0 * shear_modulus * time_step * coefficient / reference_stress);
    const double expected_increment = (trial_stress - expected_stress) / (3.0 * shear_modulus);
    bool passed = check(scaled_error(equivalent_stress(response.stress), expected_stress) < 1.0e-13,
        "linear Norton backward-Euler stress matches the "
        "analytic root");
    passed = check(std::abs(state.equivalent_creep_strain - expected_increment) < 1.0e-14,
                 "linear Norton equivalent creep increment matches the "
                 "analytic root")
             && passed;
    passed =
        check(std::abs(inelastic_trace(state.creep_strain)) < 1.0e-14, "Norton creep strain is trace free") && passed;
    passed = check(same_state(committed, fuelsim::MaterialPointState{}),
                 "Norton trial evaluation does not mutate committed state")
             && passed;
    const adlite::Scalar active_rr = adlite::Scalar::independent(0.2, 0, 1);
    const fuelsim::InelasticStressResponse active =
        material.response(active_rr, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState active_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(active.trial_state);
    passed = check(active.stress.rr.value() == response.stress.rr.value() && same_state(active_state, state),
                 "active Norton evaluation preserves the passive primal response exactly")
             && passed;
    constexpr double perturbation = 1.0e-6;
    const double plus =
        material.response(0.2 + perturbation, -0.1, -0.1, 0.0, 600.0, time_step, committed).stress.rr.value();
    const double minus =
        material.response(0.2 - perturbation, -0.1, -0.1, 0.0, 600.0, time_step, committed).stress.rr.value();
    const double finite_difference = (plus - minus) / (2.0 * perturbation);
    const double derivative_error = scaled_error(active.stress.rr.derivative(0), finite_difference);
    passed = check(derivative_error < 1.0e-9, "active Norton AD tangent matches centered finite difference") && passed;
    const adlite::Scalar zero_active_rr = adlite::Scalar::independent(0.0, 0, 1);
    const fuelsim::InelasticStressResponse zero_active =
        material.response(zero_active_rr, 0.0, 0.0, 0.0, 600.0, time_step, committed);
    constexpr double zero_perturbation = 1.0e-8;
    const double zero_plus =
        material.response(zero_perturbation, 0.0, 0.0, 0.0, 600.0, time_step, committed).stress.rr.value();
    const double zero_minus =
        material.response(-zero_perturbation, 0.0, 0.0, 0.0, 600.0, time_step, committed).stress.rr.value();
    const double zero_finite_difference = (zero_plus - zero_minus) / (2.0 * zero_perturbation);
    const double zero_derivative_error = scaled_error(zero_active.stress.rr.derivative(0), zero_finite_difference);
    passed = check(zero_derivative_error < 1.0e-9,
                 "linear Norton zero-stress AD tangent matches centered "
                 "finite difference")
             && passed;
    const fuelsim::IsotropicThermoelasticMaterial cubic_material(creep_properties(1.0e-4, 1.0, 3.0));
    const fuelsim::InelasticStressResponse cubic =
        cubic_material.response(0.2, -0.1, -0.1, 0.0, 600.0, 0.01, committed);
    const double cubic_stress = equivalent_stress(cubic.stress);
    const double cubic_residual =
        cubic_stress + 3.0 * shear_modulus * 0.01 * 1.0e-4 * std::pow(cubic_stress, 3.0) - trial_stress;
    passed = check(std::abs(cubic_residual) < 1.0e-11 * (1.0 + trial_stress),
                 "nonlinear Norton safeguarded Newton satisfies its local "
                 "equation")
             && passed;
    const fuelsim::ThermoelasticProperties small_scale_elastic =
        fuelsim::test::thermoelastic(0.0, 1.0, 1.0, 0.0, 0.0, 600.0);
    const fuelsim::IsotropicThermoelasticMaterial extreme_scale_material(
        fuelsim::test::with_norton(small_scale_elastic, 1.0, std::numeric_limits<double>::max(), 2.0));
    const fuelsim::InelasticStressResponse extreme_scale =
        extreme_scale_material.response(1.0e-16, -0.5e-16, -0.5e-16, 0.0, 600.0, 1.0, committed);
    const double extreme_scale_stress = equivalent_stress(extreme_scale.stress);
    passed = check(extreme_scale_stress > 0.0 && std::abs(extreme_scale_stress - 1.5e-16) / 1.5e-16 < 1.0e-12,
                 "Norton extreme stress normalization does not erase a "
                 "nonzero deviatoric stress")
             && passed;
    constexpr double large_coefficient = 5.0e307;
    const fuelsim::IsotropicThermoelasticMaterial large_beta_material(
        fuelsim::test::with_norton(small_scale_elastic, large_coefficient, 1.0, 2.0));
    const double one_direction = 1.0;
    const double negative_half_direction = -0.5;
    const adlite::Scalar large_beta_rr = adlite::Scalar::seeded(2.0 / 3.0, &one_direction, 1);
    const adlite::Scalar large_beta_zz = adlite::Scalar::seeded(-1.0 / 3.0, &negative_half_direction, 1);
    const adlite::Scalar large_beta_hoop = adlite::Scalar::seeded(-1.0 / 3.0, &negative_half_direction, 1);
    const fuelsim::InelasticStressResponse large_beta =
        large_beta_material.response(large_beta_rr, large_beta_zz, large_beta_hoop, 0.0, 600.0, 1.0, committed);
    const double large_beta_stress = equivalent_stress(large_beta.stress);
    const double beta = 1.5 * large_coefficient;
    const double expected_large_beta_derivative = 1.0 / (1.0 + 2.0 * (beta * large_beta_stress));
    const double actual_large_beta_derivative = large_beta.stress.rr.derivative(0);
    passed = check(actual_large_beta_derivative > 0.0 && std::isfinite(actual_large_beta_derivative)
                       && std::abs(actual_large_beta_derivative - expected_large_beta_derivative)
                                  / expected_large_beta_derivative
                              < 1.0e-12,
                 "Norton large dimensionless coefficient keeps a finite "
                 "consistent tangent")
             && passed;
    const fuelsim::IsotropicThermoelasticMaterial subnormal_material(
        fuelsim::test::with_norton(small_scale_elastic, 1.0e-300, 1.0e-300, 2.0));
    const fuelsim::InelasticStressResponse subnormal = subnormal_material.response((2.0 / 3.0) * 1.0e-300,
        (-1.0 / 3.0) * 1.0e-300,
        (-1.0 / 3.0) * 1.0e-300,
        0.0,
        600.0,
        1.0,
        committed);
    const double subnormal_stress = 1.5 * subnormal.stress.rr.value();
    const double expected_subnormal_ratio = (-1.0 + std::sqrt(7.0)) / 3.0;
    passed = check(subnormal_stress > 0.0
                       && std::abs(subnormal_stress / (expected_subnormal_ratio * 1.0e-300) - 1.0) < 1.0e-12,
                 "Norton AD hypot norm preserves a subnormal-scale trial "
                 "stress")
             && passed;
    constexpr double huge_trial_stress = 1.0e300;
    constexpr double huge_reference_stress = 1.0e100;
    constexpr double huge_creep_coefficient = 1.0e250;
    const fuelsim::IsotropicThermoelasticMaterial logarithmic_root_material(
        fuelsim::test::with_norton(small_scale_elastic, huge_creep_coefficient, huge_reference_stress, 2.0));
    const fuelsim::InelasticStressResponse logarithmic_root =
        logarithmic_root_material.response((2.0 / 3.0) * huge_trial_stress,
            (-1.0 / 3.0) * huge_trial_stress,
            (-1.0 / 3.0) * huge_trial_stress,
            0.0,
            600.0,
            1.0,
            committed);
    const double logarithmic_root_stress = 1.5 * logarithmic_root.stress.rr.value();
    const double expected_logarithmic_root =
        std::sqrt(huge_trial_stress / (1.5 * huge_creep_coefficient / (huge_reference_stress * huge_reference_stress)));
    passed = check(std::isfinite(logarithmic_root_stress) && logarithmic_root_stress > 0.0
                       && std::abs(logarithmic_root_stress / expected_logarithmic_root - 1.0) < 1.0e-12,
                 "Norton logarithmic root remains finite when its "
                 "dimensionless coefficient exceeds double range")
             && passed;
    const fuelsim::ThermoelasticProperties overflow_three_g_elastic =
        fuelsim::test::thermoelastic(0.0, 1.0, 1.6e308, 0.0, 0.0, 600.0);
    constexpr double overflow_three_g_young_modulus = 1.6e308;
    const double overflow_three_g_coefficient = (2.0 / 3.0) / overflow_three_g_young_modulus;
    const fuelsim::IsotropicThermoelasticMaterial overflow_three_g_material(
        fuelsim::test::with_norton(overflow_three_g_elastic, overflow_three_g_coefficient, 1.0, 1.0));
    constexpr double overflow_three_g_trial_stress = 48.0;
    const double overflow_three_g_strain_scale = overflow_three_g_trial_stress / overflow_three_g_young_modulus;
    const fuelsim::InelasticStressResponse overflow_three_g =
        overflow_three_g_material.response((2.0 / 3.0) * overflow_three_g_strain_scale,
            (-1.0 / 3.0) * overflow_three_g_strain_scale,
            (-1.0 / 3.0) * overflow_three_g_strain_scale,
            0.0,
            600.0,
            1.0,
            committed);
    const fuelsim::MaterialPointState overflow_three_g_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(overflow_three_g.trial_state);
    const double overflow_three_g_stress = equivalent_stress(overflow_three_g.stress);
    const double expected_overflow_three_g_creep = overflow_three_g_coefficient * overflow_three_g_stress;
    passed =
        check(std::abs(overflow_three_g_stress - 24.0) < 1.0e-12 && expected_overflow_three_g_creep > 0.0
                  && std::abs(overflow_three_g_state.equivalent_creep_strain / expected_overflow_three_g_creep - 1.0)
                         < 1.0e-12,
            "Norton stress and creep history remain consistent when 3G "
            "exceeds double range")
        && passed;
    std::cout << "m22_norton_equivalent_stress=" << equivalent_stress(response.stress) << '\n';
    std::cout << "m22_norton_equivalent_creep_strain=" << state.equivalent_creep_strain << '\n';
    std::cout << "m22_norton_material_ad_scaled_error=" << derivative_error << '\n';
    std::cout << "m22_norton_zero_stress_ad_scaled_error=" << zero_derivative_error << '\n';
    std::cout << "m22_norton_cubic_local_residual=" << cubic_residual << '\n';
    std::cout << "m22_norton_extreme_scale_equivalent_stress=" << extreme_scale_stress << '\n';
    std::cout << "m22_norton_large_beta_tangent=" << actual_large_beta_derivative << '\n';
    std::cout << "m22_norton_subnormal_equivalent_stress=" << subnormal_stress << '\n';
    std::cout << "m22_norton_logarithmic_root_equivalent_stress=" << logarithmic_root_stress << '\n';
    return passed;
}

bool test_coupled_plastic_creep_material_point() {
    constexpr double coefficient = 0.02;
    constexpr double reference_stress = 10.0;
    constexpr double exponent = 2.0;
    constexpr double yield_stress = 20.0;
    constexpr double hardening = 40.0;
    constexpr double time_step = 0.1;
    constexpr double shear_modulus = 80.0;
    constexpr double trial_stress = 48.0;
    const fuelsim::IsotropicThermoelasticMaterial material(
        coupled_properties(coefficient, reference_stress, exponent, yield_stress, hardening));
    const fuelsim::MaterialPointState committed{};
    const fuelsim::InelasticStressResponse response =
        material.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(response.trial_state);
    const double stress = equivalent_stress(response.stress);
    const double expected_creep = time_step * coefficient * std::pow(stress / reference_stress, exponent);
    const double stress_balance =
        stress + 3.0 * shear_modulus * (state.equivalent_plastic_strain + state.equivalent_creep_strain) - trial_stress;
    const double yield_residual = stress - (yield_stress + hardening * state.equivalent_plastic_strain);
    bool passed = check(state.equivalent_plastic_strain > 0.0 && state.equivalent_creep_strain > 0.0,
        "coupled update activates plastic and creep strain at one "
        "material point");
    passed = check(std::abs(state.equivalent_creep_strain - expected_creep) < 1.0e-13,
                 "coupled creep increment satisfies the backward-Euler Norton "
                 "law")
             && passed;
    passed = check(std::abs(stress_balance) < 1.0e-12, "coupled update satisfies deviatoric stress balance") && passed;
    passed = check(std::abs(yield_residual) < 1.0e-12, "coupled update ends on the hardened yield surface") && passed;
    passed = check(std::abs(inelastic_trace(state.plastic_strain)) < 1.0e-14
                       && std::abs(inelastic_trace(state.creep_strain)) < 1.0e-14,
                 "coupled plastic and creep strain increments are trace free")
             && passed;
    passed = check(same_state(committed, fuelsim::MaterialPointState{}),
                 "coupled trial evaluation does not mutate committed "
                 "state")
             && passed;
    const adlite::Scalar active_rr = adlite::Scalar::independent(0.2, 0, 1);
    const fuelsim::InelasticStressResponse active =
        material.response(active_rr, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState active_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(active.trial_state);
    passed = check(active.stress.rr.value() == response.stress.rr.value() && same_state(active_state, state),
                 "active coupled evaluation preserves the passive primal response exactly")
             && passed;
    constexpr double perturbation = 1.0e-6;
    const double plus =
        material.response(0.2 + perturbation, -0.1, -0.1, 0.0, 600.0, time_step, committed).stress.rr.value();
    const double minus =
        material.response(0.2 - perturbation, -0.1, -0.1, 0.0, 600.0, time_step, committed).stress.rr.value();
    const double finite_difference = (plus - minus) / (2.0 * perturbation);
    const double derivative_error = scaled_error(active.stress.rr.derivative(0), finite_difference);
    passed = check(derivative_error < 1.0e-9,
                 "coupled material AD tangent matches centered finite "
                 "difference")
             && passed;
    const fuelsim::IsotropicThermoelasticMaterial creep_limit(
        coupled_properties(coefficient, reference_stress, exponent, 100.0, hardening));
    const fuelsim::MaterialPointState creep_limit_state = fuelsim::IsotropicThermoelasticMaterial::state_values(
        creep_limit.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed).trial_state);
    passed =
        check(creep_limit_state.equivalent_plastic_strain == 0.0 && creep_limit_state.equivalent_creep_strain > 0.0,
            "strong creep below yield leaves the coupled plastic branch "
            "inactive")
        && passed;
    const fuelsim::IsotropicThermoelasticMaterial plastic_limit(
        coupled_properties(0.0, reference_stress, exponent, yield_stress, hardening));
    const fuelsim::MaterialPointState plastic_limit_state = fuelsim::IsotropicThermoelasticMaterial::state_values(
        plastic_limit.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed).trial_state);
    const double expected_plastic = (trial_stress - yield_stress) / (3.0 * shear_modulus + hardening);
    passed = check(std::abs(plastic_limit_state.equivalent_plastic_strain - expected_plastic) < 1.0e-14
                       && plastic_limit_state.equivalent_creep_strain == 0.0,
                 "zero Norton coefficient reduces the coupled update to J2 "
                 "plasticity")
             && passed;
    const fuelsim::IsotropicThermoelasticMaterial perfect_plastic(
        coupled_properties(coefficient, reference_stress, exponent, yield_stress, 0.0));
    const fuelsim::InelasticStressResponse perfect_response =
        perfect_plastic.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState perfect_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(perfect_response.trial_state);
    passed = check(std::abs(equivalent_stress(perfect_response.stress) - yield_stress) < 1.0e-13
                       && perfect_state.equivalent_plastic_strain > 0.0 && perfect_state.equivalent_creep_strain > 0.0,
                 "coupled perfect plasticity keeps stress on the fixed yield "
                 "surface")
             && passed;
    const fuelsim::IsotropicThermoelasticMaterial tiny_hardening(coupled_properties(coefficient,
        reference_stress,
        exponent,
        yield_stress,
        std::numeric_limits<double>::denorm_min()));
    const fuelsim::InelasticStressResponse tiny_hardening_response =
        tiny_hardening.response(0.2, -0.1, -0.1, 0.0, 600.0, time_step, committed);
    const fuelsim::MaterialPointState tiny_hardening_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(tiny_hardening_response.trial_state);
    passed =
        check(std::abs(equivalent_stress(tiny_hardening_response.stress) - equivalent_stress(perfect_response.stress))
                      < 1.0e-13
                  && std::abs(tiny_hardening_state.equivalent_plastic_strain - perfect_state.equivalent_plastic_strain)
                         < 1.0e-13
                  && std::abs(tiny_hardening_state.equivalent_creep_strain - perfect_state.equivalent_creep_strain)
                         < 1.0e-13,
            "vanishing hardening approaches the coupled perfect-plastic "
            "solution without cancellation")
        && passed;
    const fuelsim::ThermoelasticProperties large_modulus_elastic =
        fuelsim::test::thermoelastic(0.0, 1.0, 6.0e307, 0.0, 0.0, 600.0);
    const fuelsim::IsotropicThermoelasticMaterial large_modulus_material(
        coupled_properties(large_modulus_elastic, 0.0, reference_stress, exponent, yield_stress, 1.0e308));
    constexpr double large_modulus_trial_stress = 48.0;
    constexpr double large_young_modulus = 6.0e307;
    const double strain_scale = large_modulus_trial_stress / large_young_modulus;
    const fuelsim::InelasticStressResponse large_modulus_response =
        large_modulus_material.response((2.0 / 3.0) * strain_scale,
            (-1.0 / 3.0) * strain_scale,
            (-1.0 / 3.0) * strain_scale,
            0.0,
            600.0,
            time_step,
            committed);
    const fuelsim::MaterialPointState large_modulus_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(large_modulus_response.trial_state);
    const double expected_large_modulus_stress =
        yield_stress + (10.0 / 19.0) * (large_modulus_trial_stress - yield_stress);
    const double expected_large_modulus_plastic = ((large_modulus_trial_stress - yield_stress) / 1.9) * 1.0e-308;
    passed = check(std::abs(equivalent_stress(large_modulus_response.stress) - expected_large_modulus_stress) < 1.0e-12
                       && std::abs(large_modulus_state.equivalent_plastic_strain / expected_large_modulus_plastic - 1.0)
                              < 1.0e-12
                       && large_modulus_state.equivalent_creep_strain == 0.0,
                 "scaled coupled weights preserve a finite solution when "
                 "3G plus hardening exceeds double range")
             && passed;
    const fuelsim::ThermoelasticProperties overflow_three_g_elastic =
        fuelsim::test::thermoelastic(0.0, 1.0, 1.6e308, 0.0, 0.0, 600.0);
    const fuelsim::IsotropicThermoelasticMaterial overflow_three_g_material(
        coupled_properties(overflow_three_g_elastic, 0.0, reference_stress, exponent, yield_stress, 1.0e308));
    constexpr double overflow_three_g_trial_stress = 48.0;
    constexpr double overflow_three_g_young_modulus = 1.6e308;
    const double overflow_three_g_strain_scale = overflow_three_g_trial_stress / overflow_three_g_young_modulus;
    const fuelsim::InelasticStressResponse overflow_three_g_response =
        overflow_three_g_material.response((2.0 / 3.0) * overflow_three_g_strain_scale,
            (-1.0 / 3.0) * overflow_three_g_strain_scale,
            (-1.0 / 3.0) * overflow_three_g_strain_scale,
            0.0,
            600.0,
            time_step,
            committed);
    const fuelsim::MaterialPointState overflow_three_g_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(overflow_three_g_response.trial_state);
    const double expected_overflow_three_g_stress =
        yield_stress + (1.0 / 3.4) * (overflow_three_g_trial_stress - yield_stress);
    const double expected_overflow_three_g_plastic = ((overflow_three_g_trial_stress - yield_stress) / 3.4) * 1.0e-308;
    passed =
        check(
            std::isfinite(equivalent_stress(overflow_three_g_response.stress))
                && std::isfinite(overflow_three_g_state.equivalent_plastic_strain)
                && std::abs(equivalent_stress(overflow_three_g_response.stress) - expected_overflow_three_g_stress)
                       < 1.0e-12
                && std::abs(overflow_three_g_state.equivalent_plastic_strain / expected_overflow_three_g_plastic - 1.0)
                       < 1.0e-12
                && overflow_three_g_state.equivalent_creep_strain == 0.0,
            "logarithmic coupled coefficient preserves a finite solution "
            "when 3G exceeds double range")
        && passed;
    const fuelsim::ThermoelasticProperties small_modulus_elastic =
        fuelsim::test::thermoelastic(0.0, 1.0, 6.0e-309, 0.0, 0.0, 600.0);
    const fuelsim::IsotropicThermoelasticMaterial small_modulus_material(coupled_properties(small_modulus_elastic,
        0.0,
        reference_stress,
        exponent,
        std::numeric_limits<double>::denorm_min(),
        0.0));
    constexpr double small_modulus_trial_stress = 1.0e-308;
    constexpr double small_young_modulus = 6.0e-309;
    const double small_modulus_strain_scale = small_modulus_trial_stress / small_young_modulus;
    const fuelsim::InelasticStressResponse small_modulus_response =
        small_modulus_material.response((2.0 / 3.0) * small_modulus_strain_scale,
            (-1.0 / 3.0) * small_modulus_strain_scale,
            (-1.0 / 3.0) * small_modulus_strain_scale,
            0.0,
            600.0,
            time_step,
            committed);
    const fuelsim::MaterialPointState small_modulus_state =
        fuelsim::IsotropicThermoelasticMaterial::state_values(small_modulus_response.trial_state);
    const double expected_small_modulus_plastic =
        (small_modulus_trial_stress - std::numeric_limits<double>::denorm_min()) / (1.5 * small_young_modulus);
    passed = check(std::isfinite(small_modulus_state.equivalent_plastic_strain)
                       && std::abs(small_modulus_state.equivalent_plastic_strain / expected_small_modulus_plastic - 1.0)
                              < 1.0e-12
                       && small_modulus_state.equivalent_creep_strain == 0.0,
                 "scaled inverse shear modulus remains finite when 1/G exceeds "
                 "double range")
             && passed;
    std::cout << "m22_coupled_equivalent_stress=" << stress << '\n';
    std::cout << "m22_coupled_equivalent_plastic_strain=" << state.equivalent_plastic_strain << '\n';
    std::cout << "m22_coupled_equivalent_creep_strain=" << state.equivalent_creep_strain << '\n';
    std::cout << "m22_coupled_stress_balance=" << stress_balance << '\n';
    std::cout << "m22_coupled_yield_residual=" << yield_residual << '\n';
    std::cout << "m22_coupled_material_ad_scaled_error=" << derivative_error << '\n';
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
    const fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(elastic_properties()),
        100.0,
        0.0,
        fuelsim::StrainFormulation::small};
    const fuelsim::LocalValues old_temperature = {
        600.0,
        600.0,
        600.0,
        600.0,
    };
    const fuelsim::Quad4MaterialHistory history{};
    fuelsim::LocalValues uniform_state = {
        605.0,
        605.0,
        605.0,
        605.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    const fuelsim::LocalResidual balanced =
        fuelsim::compute_quad4_rz_transient(data, geometry, uniform_state, old_temperature, history, 10.0);
    bool passed = true;
    for (std::size_t row = 0; row < fuelsim::quad4_node_count; ++row)
        passed =
            check(std::abs(balanced[row]) < 1.0e-10, "uniform transient heat source balances heat capacity") && passed;
    const fuelsim::LocalResidual without_capacity = fuelsim::compute_quad4_rz_transient(data,
        geometry,
        uniform_state,
        old_temperature,
        history,
        10.0,
        nullptr,
        false);
    double maximum_without_capacity = 0.0;
    for (std::size_t row = 0; row < fuelsim::quad4_node_count; ++row)
        maximum_without_capacity = std::max(maximum_without_capacity, std::abs(without_capacity[row]));
    passed = check(maximum_without_capacity > 1.0e-6,
                 "transient thermal time term can be disabled independently of heat conduction")
             && passed;
    const fuelsim::LocalValues state = {
        605.0,
        607.0,
        604.0,
        603.0,
        -0.05,
        -0.10,
        -0.10,
        -0.05,
        0.0,
        0.0,
        0.10,
        0.10,
    };
    const fuelsim::LocalValues direction = {
        0.7,
        -0.4,
        0.3,
        -0.6,
        0.2,
        -0.4,
        0.5,
        -0.1,
        -0.3,
        0.6,
        -0.2,
        0.4,
    };
    const fuelsim::rz::LocalLinearization system =
        fuelsim::rz::linearize_quad4_rz_transient(data, geometry, state, old_temperature, history, 2.0);
    constexpr double perturbation = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += perturbation * direction[dof];
        minus[dof] -= perturbation * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        fuelsim::compute_quad4_rz_transient(data, geometry, plus, old_temperature, history, 2.0);
    const fuelsim::LocalResidual minus_residual =
        fuelsim::compute_quad4_rz_transient(data, geometry, minus, old_temperature, history, 2.0);
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::quad4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
        maximum_jacobian_error = std::max(maximum_jacobian_error, scaled_error(ad_direction, finite_difference));
    }
    passed = check(maximum_jacobian_error < 1.0e-7,
                 "transient Quad4 AD Jacobian matches centered finite "
                 "difference")
             && passed;
    const fuelsim::rz::LocalLinearization slow =
        fuelsim::rz::linearize_quad4_rz_transient(data, geometry, state, old_temperature, history, 5.0);
    double maximum_capacity_error = 0.0;
    constexpr double heat_capacity = 200.0;
    for (std::size_t row = 0; row < fuelsim::quad4_node_count; ++row) {
        for (std::size_t column = 0; column < fuelsim::quad4_node_count; ++column) {
            double mass = 0.0;
            for (const fuelsim::RzQuadraturePoint& point : geometry.points)
                mass += point.weighted_measure * point.shape[row] * point.shape[column];
            const double expected = heat_capacity * (1.0 / 2.0 - 1.0 / 5.0) * mass;
            const double actual = system.jacobian[row * fuelsim::quad4_local_dof_count + column]
                                  - slow.jacobian[row * fuelsim::quad4_local_dof_count + column];
            maximum_capacity_error = std::max(maximum_capacity_error, scaled_error(actual, expected));
        }
    }
    passed = check(maximum_capacity_error < 1.0e-13,
                 "transient Jacobian contains the consistent capacity "
                 "matrix divided by dt")
             && passed;
    const fuelsim::Quad4MaterialHistory trial =
        fuelsim::compute_quad4_rz_transient_update(data, geometry, state, old_temperature, history, 2.0);
    for (std::size_t q = 0; q < trial.size(); ++q)
        passed = check(same_inelastic_state(trial[q], history[q]),
                     "elastic transient element leaves inelastic history "
                     "unchanged")
                 && passed;
    std::cout << "m21_disabled_capacity_residual_maximum=" << maximum_without_capacity << '\n';
    std::cout << "m21_element_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
    std::cout << "m21_capacity_matrix_maximum_scaled_error=" << maximum_capacity_error << '\n';
    return passed;
}

bool test_coupled_transient_element_jacobian(fuelsim::StrainFormulation strain_formulation,
    const std::string& formulation_name) {
    const fuelsim::Quad4RzGeometry geometry = test_geometry();
    const fuelsim::ThermoelasticProperties properties = coupled_properties(0.02, 10.0, 2.0, 20.0, 40.0);
    const fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(properties), 0.0, 1.25, strain_formulation};
    fuelsim::LocalValues committed_state = {
        600.0,
        600.0,
        600.0,
        600.0,
    };
    fuelsim::Quad4MaterialHistory history{};
    fuelsim::LocalValues state = {
        600.0,
        600.0,
        600.0,
        600.0,
        -0.1,
        -0.2,
        -0.2,
        -0.1,
        0.0,
        0.0,
        0.2,
        0.2,
    };
    if (strain_formulation == fuelsim::StrainFormulation::finite) {
        const std::array<fuelsim::RzPoint, 4> coordinates = {{
            {1.0, 0.0},
            {2.0, 0.0},
            {2.0, 1.0},
            {1.0, 1.0},
        }};
        for (std::size_t node = 0; node < coordinates.size(); ++node) {
            const double r = coordinates[node].r;
            const double z = coordinates[node].z;
            committed_state[4 + node] = 0.02 * r + 0.03 * z;
            committed_state[8 + node] = -0.01 * r + 0.01 * z;
            state[4 + node] = 0.05 * r + 0.12 * z;
            state[8 + node] = -0.06 * r + 0.03 * z;
        }
        for (fuelsim::MaterialPointState& point : history) {
            point.elastic_strain = {0.020, -0.012, -0.008, 0.006};
            point.plastic_strain = {0.012, -0.007, -0.005, 0.004};
            point.creep_strain = {-0.008, 0.005, 0.003, -0.002};
            point.equivalent_plastic_strain = 0.02;
            point.equivalent_creep_strain = 0.01;
        }
    }
    const fuelsim::LocalValues direction = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.2,
        -0.4,
        0.5,
        -0.1,
        -0.3,
        0.6,
        -0.2,
        0.4,
    };
    constexpr double time_step = 0.1;
    const fuelsim::rz::LocalLinearization system =
        fuelsim::rz::linearize_quad4_rz_transient(data, geometry, state, committed_state, history, time_step);
    constexpr double perturbation = 1.0e-6;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += perturbation * direction[dof];
        minus[dof] -= perturbation * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        fuelsim::compute_quad4_rz_transient(data, geometry, plus, committed_state, history, time_step);
    const fuelsim::LocalResidual minus_residual =
        fuelsim::compute_quad4_rz_transient(data, geometry, minus, committed_state, history, time_step);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::quad4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
        maximum_error = std::max(maximum_error, scaled_error(ad_direction, finite_difference));
    }
    const fuelsim::Quad4MaterialHistory trial =
        fuelsim::compute_quad4_rz_transient_update(data, geometry, state, committed_state, history, time_step);
    bool both_histories_active = true;
    for (const fuelsim::MaterialPointState& point : trial) {
        both_histories_active =
            point.equivalent_plastic_strain > 0.0 && point.equivalent_creep_strain > 0.0 && both_histories_active;
    }
    bool passed = check(maximum_error < 1.0e-7,
        formulation_name
            + " coupled transient Quad4 AD Jacobian matches centered "
              "finite difference");
    passed = check(both_histories_active,
                 formulation_name
                     + " coupled transient Quad4 activates both histories at "
                       "every quadrature point")
             && passed;
    std::cout << formulation_name << "_coupled_element_jacobian_maximum_scaled_error=" << maximum_error << '\n';
    return passed;
}

bool test_coupled_transient_element_jacobians() {
    bool passed = test_coupled_transient_element_jacobian(fuelsim::StrainFormulation::small, "m22_small_strain");
    passed = test_coupled_transient_element_jacobian(fuelsim::StrainFormulation::finite, "m41_finite_strain") && passed;
    return passed;
}

fuelsim::SpatialDefinition transaction_definition() {
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", coupled_properties(1.0e-4, 1.0, 1.0, 1.0, 10.0), 100.0, 600.0});
    definition.boundary_conditions.push_back({"axis_radial",
        fuelsim::BoundaryConditionType::dirichlet,
        "solid_inner",
        fuelsim::Field::radial_displacement,
        0.0});
    definition.boundary_conditions.push_back({"bottom_axial",
        fuelsim::BoundaryConditionType::dirichlet,
        "solid_bottom",
        fuelsim::Field::axial_displacement,
        0.0});
    definition.boundary_conditions.push_back({"outer_temperature",
        fuelsim::BoundaryConditionType::dirichlet,
        "solid_outer",
        fuelsim::Field::temperature,
        600.0});
    return definition;
}

bool test_problem_history_transaction() {
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::test::make_disconnected_annular_mesh({{1, "solid", 0.0, 1.0, 1.0, 1, 1}});
    fuelsim::TransientProblem problem(transaction_definition(), mesh);
    const std::vector<double> initial_solution = problem.committed_solution();
    const fuelsim::MaterialPointState initial_history = fuelsim::rz::ProblemAccess::material_history(problem, 0, 0)[0];
    problem.begin_time_step({1.0, 0.5});
    std::vector<double> trial_solution = initial_solution;
    const auto& dofs = fuelsim::rz::ProblemAccess::dof_map(problem);
    for (std::size_t local_node = 0; local_node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size();
        ++local_node) {
        const fuelsim::RzPoint& point = fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes()[local_node];
        const std::size_t global_node = fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) + local_node;
        trial_solution[dofs.dof(fuelsim::Field::radial_displacement, global_node)] = -0.05 * point.r;
        trial_solution[dofs.dof(fuelsim::Field::axial_displacement, global_node)] = 0.10 * point.z;
    }
    const fuelsim::LocalValues local_trial = fuelsim::rz::ProblemAccess::contribution_state(problem, 0, trial_solution);
    (void)fuelsim::rz::ProblemAccess::contribution_residual(problem, 0, local_trial);
    (void)fuelsim::rz::ProblemAccess::linearize_contribution(problem, 0, local_trial);
    bool passed = check(same_state(fuelsim::rz::ProblemAccess::material_history(problem, 0, 0)[0], initial_history),
        "Newton residual and Jacobian callbacks do not mutate "
        "committed history");
    problem.commit_time_step(trial_solution);
    const fuelsim::RegionStateSummary committed_history =
        fuelsim::rz::ProblemAccess::summarize_region_history(problem, 0);
    const double committed_plastic = committed_history.maximum_equivalent_plastic_strain;
    const double committed_creep = committed_history.maximum_equivalent_creep_strain;
    passed =
        check(problem.committed_time() == 1.0 && problem.committed_load_factor() == 0.5 && !problem.time_step_active(),
            "accepted M2 time step commits time, load, and phase")
        && passed;
    passed = check(committed_plastic > 0.0 && committed_creep > 0.0,
                 "accepted M2 time step commits coupled plastic and creep "
                 "history")
             && passed;
    passed =
        check(problem.committed_solution() == trial_solution, "accepted M2 time step commits the converged nodal state")
        && passed;
    problem.begin_time_step({2.0, 0.75});
    problem.rollback_time_step();
    passed = check(problem.committed_time() == 1.0 && problem.committed_load_factor() == 0.5
                       && fuelsim::rz::ProblemAccess::region_kernel_data(problem, 0).volumetric_heat_source == 50.0,
                 "rollback restores the committed time and heat load")
             && passed;
    passed =
        check(fuelsim::rz::ProblemAccess::summarize_region_history(problem, 0).maximum_equivalent_creep_strain
                      == committed_creep
                  && fuelsim::rz::ProblemAccess::summarize_region_history(problem, 0).maximum_equivalent_plastic_strain
                         == committed_plastic
                  && problem.committed_solution() == trial_solution,
            "rollback preserves committed nodal and material history")
        && passed;
    problem.begin_time_step({2.0, 0.75});
    std::vector<double> invalid_solution = trial_solution;
    invalid_solution[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature,
        fuelsim::rz::ProblemAccess::region_node_offset(problem, 0))] = -100.0;
    bool invalid_commit_threw = false;
    try {
        problem.commit_time_step(invalid_solution);
    } catch (const std::domain_error&) {
        invalid_commit_threw = true;
    }
    passed = check(invalid_commit_threw && problem.time_step_active(),
                 "M2 rejects a nonpositive nodal temperature before "
                 "changing committed state")
             && passed;
    problem.rollback_time_step();
    passed =
        check(problem.committed_time() == 1.0 && problem.committed_solution() == trial_solution
                  && fuelsim::rz::ProblemAccess::summarize_region_history(problem, 0).maximum_equivalent_creep_strain
                         == committed_creep
                  && fuelsim::rz::ProblemAccess::summarize_region_history(problem, 0).maximum_equivalent_plastic_strain
                         == committed_plastic,
            "failed commit retains the previous complete transaction")
        && passed;
    bool inactive_threw = false;
    try {
        (void)fuelsim::rz::ProblemAccess::contribution_residual(problem, 0, local_trial);
    } catch (const std::logic_error&) {
        inactive_threw = true;
    }
    passed = check(inactive_threw, "M2 residual evaluation without an active step is rejected") && passed;
    std::cout << "m21_committed_time=" << problem.committed_time() << '\n';
    std::cout << "m22_committed_maximum_plastic_strain=" << committed_plastic << '\n';
    std::cout << "m22_committed_maximum_creep_strain=" << committed_creep << '\n';
    return passed;
}
} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    try {
        bool passed = true;
        passed = test_builtin_material_parameter_order() && passed;
        passed = test_registered_material_functions() && passed;
        passed = test_objective_incremental_history_rotation() && passed;
        passed = test_j2_plasticity_material_point() && passed;
        passed = test_norton_creep_material_point() && passed;
        passed = test_coupled_plastic_creep_material_point() && passed;
        passed = test_temperature_active_inelastic_properties() && passed;
        passed = test_transient_element() && passed;
        passed = test_coupled_transient_element_jacobians() && passed;
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
