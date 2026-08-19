#include "fuelsim/core/material_functions.hpp"
#include "detail/fnv_hash.hpp"
#include "fuelsim/core/material.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
MaterialParameters::MaterialParameters(std::vector<MaterialParameterValue> values) : _values(std::move(values)) {
    for (std::size_t index = 0; index < _values.size(); ++index) {
        if (_values[index].name.empty()) throw std::invalid_argument("Material parameter name must not be empty");
        if (!std::isfinite(_values[index].value))
            throw std::invalid_argument("Material parameter '" + _values[index].name + "' must be finite");
        for (std::size_t previous = 0; previous < index; ++previous)
            if (_values[previous].name == _values[index].name)
                throw std::invalid_argument("Duplicate material parameter '" + _values[index].name + "'");
    }
}

double MaterialParameters::value(const std::string& name) const {
    const auto found = std::find_if(_values.begin(), _values.end(),
        [&name](const MaterialParameterValue& parameter) { return parameter.name == name; });
    if (found == _values.end()) throw std::out_of_range("Unknown material parameter '" + name + "'");
    return found->value;
}

namespace {
void material_hash_string(std::uint64_t& hash, const std::string& value) noexcept {
    detail::fnv1a_bytes(hash, value.data(), value.size());
    const unsigned char terminator = 0;
    detail::fnv1a_bytes(hash, &terminator, 1);
}

void material_hash_u32(std::uint64_t& hash, std::uint32_t value) noexcept {
    detail::fnv1a_bytes(hash, &value, sizeof(value));
}

void material_hash_double(std::uint64_t& hash, double value) noexcept {
    detail::fnv1a_bytes(hash, &value, sizeof(value));
}

void hash_instance(std::uint64_t& hash, const std::string& name, std::uint32_t version,
    const MaterialParameters& parameters) noexcept {
    material_hash_string(hash, name);
    material_hash_u32(hash, version);
    for (const MaterialParameterValue& parameter : parameters.values()) {
        material_hash_string(hash, parameter.name);
        material_hash_double(hash, parameter.value);
    }
}

void validate_registration(
    const std::string& name, const std::vector<MaterialParameterDefinition>& parameters, bool has_function) {
    if (name.empty()) throw std::invalid_argument("Material function name must not be empty");
    if (!has_function) throw std::invalid_argument("Material function '" + name + "' must not be null");
    for (std::size_t index = 0; index < parameters.size(); ++index) {
        if (parameters[index].name.empty())
            throw std::invalid_argument("Material function '" + name + "' has an empty parameter name");
        if (parameters[index].unit.empty())
            throw std::invalid_argument(
                "Material function '" + name + "' parameter '" + parameters[index].name + "' has an empty SI unit");
        for (std::size_t previous = 0; previous < index; ++previous)
            if (parameters[previous].name == parameters[index].name)
                throw std::invalid_argument(
                    "Material function '" + name + "' repeats parameter '" + parameters[index].name + "'");
    }
}

void validate_bound_values(const std::string& function_name,
    const std::vector<MaterialParameterDefinition>& definitions, const std::vector<MaterialParameterValue>& values) {
    for (const MaterialParameterValue& value : values) {
        const bool known = std::any_of(definitions.begin(), definitions.end(),
            [&](const MaterialParameterDefinition& definition) { return definition.name == value.name; });
        if (!known)
            throw std::invalid_argument(
                "Material function '" + function_name + "' has unknown key '" + value.name + "'");
    }
    for (const MaterialParameterDefinition& definition : definitions) {
        const auto found = std::find_if(values.begin(), values.end(),
            [&](const MaterialParameterValue& value) { return value.name == definition.name; });
        if (found == values.end())
            throw std::invalid_argument(
                "Material function '" + function_name + "' is missing required key '" + definition.name + "'");
    }
    if (values.size() != definitions.size())
        throw std::invalid_argument("Material function '" + function_name + "' parameter count does not match");
}

std::vector<MaterialParameterValue> ordered_values(const std::string& function_name,
    const std::vector<MaterialParameterDefinition>& definitions, const std::vector<MaterialParameterValue>& values) {
    validate_bound_values(function_name, definitions, values);
    std::vector<MaterialParameterValue> result;
    result.reserve(definitions.size());
    for (const MaterialParameterDefinition& definition : definitions) {
        const auto found = std::find_if(values.begin(), values.end(),
            [&](const MaterialParameterValue& value) { return value.name == definition.name; });
        result.push_back(*found);
    }
    return result;
}

// Binding preserves each strict named registration's declaration order for these built-in evaluations.
void constant_thermophysical(const ThermoelasticFunctionInput& input, ThermalPropertyOutput& output) {
    output.conductivity = input.parameters->values()[0].value;
    output.density = input.parameters->values()[1].value;
    output.specific_heat = input.parameters->values()[2].value;
}

void inverse_temperature_thermophysical(const ThermoelasticFunctionInput& input, ThermalPropertyOutput& output) {
    output.conductivity = input.parameters->values()[0].value / input.temperature + input.parameters->values()[1].value;
    output.density = input.parameters->values()[2].value;
    output.specific_heat = input.parameters->values()[3].value;
}

void constant_isotropic_elasticity(const ThermoelasticFunctionInput& input, ElasticPropertyOutput& output) {
    output.young_modulus = input.parameters->values()[0].value;
    output.poisson_ratio = input.parameters->values()[1].value;
}

void linear_temperature_isotropic_elasticity(const ThermoelasticFunctionInput& input, ElasticPropertyOutput& output) {
    const adlite::Scalar temperature_change = input.temperature - input.parameters->values()[2].value;
    output.young_modulus =
        input.parameters->values()[0].value + input.parameters->values()[3].value * temperature_change;
    output.poisson_ratio =
        input.parameters->values()[1].value + input.parameters->values()[4].value * temperature_change;
}

void isotropic_thermal_expansion(const ThermoelasticFunctionInput& input, SymmetricTensor3& output) {
    const adlite::Scalar value =
        input.parameters->values()[0].value * (input.temperature - input.parameters->values()[1].value);
    output = {value, value, value, 0.0, 0.0, 0.0};
}

void linear_temperature_isotropic_thermal_expansion(const ThermoelasticFunctionInput& input, SymmetricTensor3& output) {
    const adlite::Scalar temperature_change = input.temperature - input.parameters->values()[1].value;
    const adlite::Scalar coefficient =
        input.parameters->values()[0].value + input.parameters->values()[2].value * temperature_change;
    const adlite::Scalar value = coefficient * temperature_change;
    output = {value, value, value, 0.0, 0.0, 0.0};
}

adlite::Scalar norton_creep_rate(const CreepRateInput& input) {
    if (!(input.equivalent_stress.value() > 0.0)) return 0.0;
    return input.parameters->value("coefficient") *
           adlite::pow(input.equivalent_stress / input.parameters->value("reference_stress"),
               input.parameters->value("stress_exponent"));
}

adlite::Scalar linear_temperature_norton_creep_rate(const CreepRateInput& input) {
    const adlite::Scalar temperature_change = input.temperature - input.parameters->value("reference_temperature");
    const adlite::Scalar coefficient =
        input.parameters->value("coefficient") +
        input.parameters->value("coefficient_temperature_coefficient") * temperature_change;
    const adlite::Scalar reference_stress =
        input.parameters->value("reference_stress") +
        input.parameters->value("reference_stress_temperature_coefficient") * temperature_change;
    const adlite::Scalar exponent =
        input.parameters->value("stress_exponent") +
        input.parameters->value("stress_exponent_temperature_coefficient") * temperature_change;
    if (!(input.equivalent_stress.value() > 0.0)) return 0.0;
    return coefficient * adlite::pow(input.equivalent_stress / reference_stress, exponent);
}

adlite::Scalar linear_isotropic_flow_stress(const PlasticFlowStressInput& input) {
    return input.parameters->value("yield_stress") +
           input.parameters->value("hardening_modulus") * input.equivalent_plastic_strain;
}

adlite::Scalar linear_temperature_isotropic_flow_stress(const PlasticFlowStressInput& input) {
    const adlite::Scalar temperature_change = input.temperature - input.parameters->value("reference_temperature");
    const adlite::Scalar yield_stress =
        input.parameters->value("yield_stress") +
        input.parameters->value("yield_stress_temperature_coefficient") * temperature_change;
    const adlite::Scalar hardening = input.parameters->value("hardening_modulus") +
                                     input.parameters->value("hardening_temperature_coefficient") * temperature_change;
    return yield_stress + hardening * input.equivalent_plastic_strain;
}
} // namespace

std::uint64_t MaterialFunctionSet::signature() const noexcept {
    std::uint64_t hash = detail::fnv1a_offset;
    material_hash_string(hash, name);
    hash_instance(hash, thermal.name, thermal.version, thermal.parameters);
    hash_instance(hash, elasticity.name, elasticity.version, elasticity.parameters);
    for (const EigenstrainFunctionInstance& eigenstrain : eigenstrains) {
        material_hash_string(hash, eigenstrain.instance_name);
        hash_instance(hash, eigenstrain.name, eigenstrain.version, eigenstrain.parameters);
    }
    hash_instance(hash, creep.name, creep.version, creep.parameters);
    hash_instance(hash, plasticity.name, plasticity.version, plasticity.parameters);
    return hash;
}

void MaterialFunctionRegistry::add_registration(std::string name, std::vector<MaterialParameterDefinition> parameters,
    std::uint32_t version, Category category, Function function, bool has_function, const char* category_name) {
    validate_registration(name, parameters, has_function);
    if (version == 0)
        throw std::invalid_argument(std::string(category_name) + " material function version must be positive");
    if (std::any_of(_registrations.begin(), _registrations.end(),
            [&name, category](const Registration& entry) { return entry.category == category && entry.name == name; }))
        throw std::invalid_argument("Duplicate " + std::string(category_name) + " material function '" + name + "'");
    _registrations.push_back({std::move(name), version, std::move(parameters), category, function});
}

const MaterialFunctionRegistry::Registration& MaterialFunctionRegistry::find_registration(
    Category category, const std::string& name, const char* category_name) const {
    const auto found = std::find_if(_registrations.begin(), _registrations.end(),
        [&name, category](const Registration& entry) { return entry.category == category && entry.name == name; });
    if (found == _registrations.end())
        throw std::invalid_argument("Unknown " + std::string(category_name) + " material function '" + name + "'");
    return *found;
}

void MaterialFunctionRegistry::add_thermal(std::string name, std::vector<MaterialParameterDefinition> parameters,
    ThermalPropertyFunction function, std::uint32_t version) {
    add_registration(
        std::move(name), std::move(parameters), version, Category::thermal, function, function != nullptr, "thermal");
}

void MaterialFunctionRegistry::add_elasticity(std::string name, std::vector<MaterialParameterDefinition> parameters,
    ElasticPropertyFunction function, std::uint32_t version) {
    add_registration(std::move(name), std::move(parameters), version, Category::elasticity, function,
        function != nullptr, "elasticity");
}

void MaterialFunctionRegistry::add_eigenstrain(std::string name, std::vector<MaterialParameterDefinition> parameters,
    EigenstrainFunction function, std::uint32_t version) {
    add_registration(std::move(name), std::move(parameters), version, Category::eigenstrain, function,
        function != nullptr, "eigenstrain");
}

void MaterialFunctionRegistry::add_creep(std::string name, std::vector<MaterialParameterDefinition> parameters,
    CreepRateFunction function, std::uint32_t version) {
    add_registration(
        std::move(name), std::move(parameters), version, Category::creep, function, function != nullptr, "creep");
}

void MaterialFunctionRegistry::add_plasticity(std::string name, std::vector<MaterialParameterDefinition> parameters,
    PlasticFlowStressFunction function, std::uint32_t version) {
    add_registration(std::move(name), std::move(parameters), version, Category::plasticity, function,
        function != nullptr, "plasticity");
}

ThermalFunctionInstance MaterialFunctionRegistry::bind_thermal(
    const std::string& name, std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::thermal, name, "thermal");
    return {name, entry.version, MaterialParameters(ordered_values(name, entry.parameters, values)),
        std::get<ThermalPropertyFunction>(entry.function)};
}

ElasticFunctionInstance MaterialFunctionRegistry::bind_elasticity(
    const std::string& name, std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::elasticity, name, "elasticity");
    return {name, entry.version, MaterialParameters(ordered_values(name, entry.parameters, values)),
        std::get<ElasticPropertyFunction>(entry.function)};
}

EigenstrainFunctionInstance MaterialFunctionRegistry::bind_eigenstrain(
    const std::string& instance_name, const std::string& name, std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::eigenstrain, name, "eigenstrain");
    return {instance_name, name, entry.version, MaterialParameters(ordered_values(name, entry.parameters, values)),
        std::get<EigenstrainFunction>(entry.function)};
}

CreepFunctionInstance MaterialFunctionRegistry::bind_creep(
    const std::string& name, std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::creep, name, "creep");
    return {name, entry.version, MaterialParameters(ordered_values(name, entry.parameters, values)),
        std::get<CreepRateFunction>(entry.function)};
}

PlasticFunctionInstance MaterialFunctionRegistry::bind_plasticity(
    const std::string& name, std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::plasticity, name, "plasticity");
    return {name, entry.version, MaterialParameters(ordered_values(name, entry.parameters, values)),
        std::get<PlasticFlowStressFunction>(entry.function)};
}

MaterialFunctionRegistry make_builtin_material_function_registry() {
    MaterialFunctionRegistry registry;
    registry.add_thermal("constant_thermophysical",
        {{"conductivity", "W/(m*K)"}, {"density", "kg/m^3"}, {"specific_heat", "J/(kg*K)"}}, &constant_thermophysical);
    registry.add_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", "W/m"}, {"conductivity_constant", "W/(m*K)"}, {"density", "kg/m^3"},
            {"specific_heat", "J/(kg*K)"}},
        &inverse_temperature_thermophysical);
    registry.add_elasticity(
        "constant_isotropic", {{"young_modulus", "Pa"}, {"poisson_ratio", "1"}}, &constant_isotropic_elasticity);
    registry.add_elasticity("linear_temperature_isotropic",
        {{"young_modulus", "Pa"}, {"poisson_ratio", "1"}, {"reference_temperature", "K"},
            {"young_modulus_temperature_coefficient", "Pa/K"}, {"poisson_ratio_temperature_coefficient", "1/K"}},
        &linear_temperature_isotropic_elasticity);
    registry.add_eigenstrain("isotropic_thermal_expansion",
        {{"thermal_expansion", "1/K"}, {"reference_temperature", "K"}}, &isotropic_thermal_expansion);
    registry.add_eigenstrain("linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", "1/K"}, {"reference_temperature", "K"},
            {"thermal_expansion_temperature_coefficient", "1/K^2"}},
        &linear_temperature_isotropic_thermal_expansion);
    registry.add_creep(
        "norton", {{"coefficient", "1/s"}, {"reference_stress", "Pa"}, {"stress_exponent", "1"}}, &norton_creep_rate);
    registry.add_creep("linear_temperature_norton",
        {{"coefficient", "1/s"}, {"reference_stress", "Pa"}, {"stress_exponent", "1"}, {"reference_temperature", "K"},
            {"coefficient_temperature_coefficient", "1/(s*K)"}, {"reference_stress_temperature_coefficient", "Pa/K"},
            {"stress_exponent_temperature_coefficient", "1/K"}},
        &linear_temperature_norton_creep_rate);
    registry.add_plasticity("linear_isotropic_hardening", {{"yield_stress", "Pa"}, {"hardening_modulus", "Pa"}},
        &linear_isotropic_flow_stress);
    registry.add_plasticity("linear_temperature_isotropic_hardening",
        {{"yield_stress", "Pa"}, {"hardening_modulus", "Pa"}, {"reference_temperature", "K"},
            {"yield_stress_temperature_coefficient", "Pa/K"}, {"hardening_temperature_coefficient", "Pa/K"}},
        &linear_temperature_isotropic_flow_stress);
    return registry;
}

} // namespace fuelsim
