#pragma once
#include "fuelsim/core/material.hpp"
#include <memory>
#include <utility>

namespace fuelsim::test {
inline ThermoelasticProperties thermoelastic(double conductivity_inverse_temperature, double conductivity_offset,
    double young_modulus, double poisson_ratio, double thermal_expansion, double reference_temperature,
    double young_modulus_temperature_coefficient = 0.0, double poisson_ratio_temperature_coefficient = 0.0,
    double thermal_expansion_temperature_coefficient = 0.0, double density = 1.0, double specific_heat = 1.0) {
    MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "test_material";
    functions->thermal = registry.bind_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", conductivity_inverse_temperature},
            {"conductivity_constant", conductivity_offset}, {"density", density}, {"specific_heat", specific_heat}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", young_modulus}, {"poisson_ratio", poisson_ratio},
            {"reference_temperature", reference_temperature},
            {"young_modulus_temperature_coefficient", young_modulus_temperature_coefficient},
            {"poisson_ratio_temperature_coefficient", poisson_ratio_temperature_coefficient}});
    functions->eigenstrains.push_back(
        registry.bind_eigenstrain("thermal", "linear_temperature_isotropic_thermal_expansion",
            {{"thermal_expansion", thermal_expansion}, {"reference_temperature", reference_temperature},
                {"thermal_expansion_temperature_coefficient", thermal_expansion_temperature_coefficient}}));
    return {std::move(functions), young_modulus};
}

inline ThermoelasticProperties with_norton(ThermoelasticProperties material, double coefficient,
    double reference_stress, double stress_exponent, double reference_temperature = 600.0,
    double coefficient_temperature_coefficient = 0.0, double reference_stress_temperature_coefficient = 0.0,
    double stress_exponent_temperature_coefficient = 0.0) {
    MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>(*material.functions);
    functions->creep = registry.bind_creep("linear_temperature_norton",
        {{"coefficient", coefficient}, {"reference_stress", reference_stress}, {"stress_exponent", stress_exponent},
            {"reference_temperature", reference_temperature},
            {"coefficient_temperature_coefficient", coefficient_temperature_coefficient},
            {"reference_stress_temperature_coefficient", reference_stress_temperature_coefficient},
            {"stress_exponent_temperature_coefficient", stress_exponent_temperature_coefficient}});
    material.functions = std::move(functions);
    return material;
}

inline ThermoelasticProperties with_plasticity(ThermoelasticProperties material, double yield_stress,
    double hardening_modulus, double reference_temperature = 600.0, double yield_stress_temperature_coefficient = 0.0,
    double hardening_temperature_coefficient = 0.0) {
    MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>(*material.functions);
    functions->plasticity = registry.bind_plasticity("linear_temperature_isotropic_hardening",
        {{"yield_stress", yield_stress}, {"hardening_modulus", hardening_modulus},
            {"reference_temperature", reference_temperature},
            {"yield_stress_temperature_coefficient", yield_stress_temperature_coefficient},
            {"hardening_temperature_coefficient", hardening_temperature_coefficient}});
    material.functions = std::move(functions);
    return material;
}
} // namespace fuelsim::test
