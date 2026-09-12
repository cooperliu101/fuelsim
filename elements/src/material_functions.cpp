#include "material_functions.hpp"
#include "fnv_hash.hpp"
#include "material.hpp"
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
        if (_values[index].name.empty())
            throw std::invalid_argument("Material parameter name must not be empty");
        if (!std::isfinite(_values[index].value))
            throw std::invalid_argument("Material parameter '" + _values[index].name + "' must be finite");
        for (std::size_t previous = 0; previous < index; ++previous)
            if (_values[previous].name == _values[index].name)
                throw std::invalid_argument("Duplicate material parameter '" + _values[index].name + "'");
    }
}

double MaterialParameters::value(const std::string& name) const {
    const auto found = std::find_if(_values.begin(), _values.end(), [&name](const MaterialParameterValue& parameter) {
        return parameter.name == name;
    });
    if (found == _values.end())
        throw std::out_of_range("Unknown material parameter '" + name + "'");
    return found->value;
}

namespace {
void material_hash_string(std::uint64_t& hash, const std::string& value) noexcept {
    hashing::fnv1a_bytes(hash, value.data(), value.size());
    const unsigned char terminator = 0;
    hashing::fnv1a_bytes(hash, &terminator, 1);
}

void material_hash_u32(std::uint64_t& hash, std::uint32_t value) noexcept {
    hashing::fnv1a_bytes(hash, &value, sizeof(value));
}

void material_hash_double(std::uint64_t& hash, double value) noexcept {
    hashing::fnv1a_bytes(hash, &value, sizeof(value));
}

void hash_instance(std::uint64_t& hash,
    const std::string& name,
    std::uint32_t version,
    const MaterialParameters& parameters) noexcept {
    material_hash_string(hash, name);
    material_hash_u32(hash, version);
    for (const MaterialParameterValue& parameter : parameters.values()) {
        material_hash_string(hash, parameter.name);
        material_hash_double(hash, parameter.value);
    }
}

void validate_registration(const std::string& name,
    const std::vector<MaterialParameterDefinition>& parameters,
    bool has_function) {
    if (name.empty())
        throw std::invalid_argument("Material function name must not be empty");
    if (!has_function)
        throw std::invalid_argument("Material function '" + name + "' must not be null");
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

std::vector<MaterialParameterValue> ordered_values(const std::string& function_name,
    const std::vector<MaterialParameterDefinition>& definitions,
    const std::vector<MaterialParameterValue>& values) {
    for (const MaterialParameterValue& value : values) {
        const bool known = std::any_of(definitions.begin(),
            definitions.end(),
            [&](const MaterialParameterDefinition& definition) { return definition.name == value.name; });
        if (!known)
            throw std::invalid_argument(
                "Material function '" + function_name + "' has unknown key '" + value.name + "'");
    }
    std::vector<MaterialParameterValue> result;
    result.reserve(definitions.size());
    for (const MaterialParameterDefinition& definition : definitions) {
        const auto found = std::find_if(values.begin(), values.end(), [&](const MaterialParameterValue& value) {
            return value.name == definition.name;
        });
        if (found == values.end())
            throw std::invalid_argument(
                "Material function '" + function_name + "' is missing required key '" + definition.name + "'");
        result.push_back(*found);
    }
    if (values.size() != definitions.size())
        throw std::invalid_argument("Material function '" + function_name + "' parameter count does not match");
    return result;
}

ThermalPropertyEvaluator constant_thermophysical(const MaterialParameters& named) {
    const double conductivity = named.value("conductivity");
    const double density = named.value("density");
    const double specific_heat = named.value("specific_heat");
    return [conductivity, density, specific_heat](const ThermoelasticFunctionInput&, ThermalPropertyOutput& output) {
        output.conductivity = conductivity;
        output.density = density;
        output.specific_heat = specific_heat;
    };
}

ThermalPropertyEvaluator inverse_temperature_thermophysical(const MaterialParameters& named) {
    const double conductivity_inverse_temperature = named.value("conductivity_inverse_temperature");
    const double conductivity_constant = named.value("conductivity_constant");
    const double density = named.value("density");
    const double specific_heat = named.value("specific_heat");
    return [conductivity_inverse_temperature, conductivity_constant, density, specific_heat](
               const ThermoelasticFunctionInput& input,
               ThermalPropertyOutput& output) {
        output.conductivity = conductivity_inverse_temperature / input.temperature + conductivity_constant;
        output.density = density;
        output.specific_heat = specific_heat;
    };
}

ThermalPropertyEvaluator linear_temperature_thermophysical(const MaterialParameters& named) {
    const double conductivity = named.value("conductivity");
    const double density = named.value("density");
    const double specific_heat = named.value("specific_heat");
    const double reference_temperature = named.value("reference_temperature");
    const double conductivity_temperature_coefficient = named.value("conductivity_temperature_coefficient");
    const double density_temperature_coefficient = named.value("density_temperature_coefficient");
    const double specific_heat_temperature_coefficient = named.value("specific_heat_temperature_coefficient");
    return [conductivity,
               density,
               specific_heat,
               reference_temperature,
               conductivity_temperature_coefficient,
               density_temperature_coefficient,
               specific_heat_temperature_coefficient](const ThermoelasticFunctionInput& input,
               ThermalPropertyOutput& output) {
        const adlite::Scalar temperature_change = input.temperature - reference_temperature;
        output.conductivity = conductivity + conductivity_temperature_coefficient * temperature_change;
        output.density = density + density_temperature_coefficient * temperature_change;
        output.specific_heat = specific_heat + specific_heat_temperature_coefficient * temperature_change;
    };
}

ElasticPropertyEvaluator constant_isotropic_elasticity(const MaterialParameters& named) {
    const double young_modulus = named.value("young_modulus");
    const double poisson_ratio = named.value("poisson_ratio");
    return [young_modulus, poisson_ratio](const ThermoelasticFunctionInput&, ElasticPropertyOutput& output) {
        output.young_modulus = young_modulus;
        output.poisson_ratio = poisson_ratio;
    };
}

ElasticPropertyEvaluator linear_temperature_isotropic_elasticity(const MaterialParameters& named) {
    const double young_modulus = named.value("young_modulus");
    const double poisson_ratio = named.value("poisson_ratio");
    const double reference_temperature = named.value("reference_temperature");
    const double young_modulus_temperature_coefficient = named.value("young_modulus_temperature_coefficient");
    const double poisson_ratio_temperature_coefficient = named.value("poisson_ratio_temperature_coefficient");
    return [young_modulus,
               poisson_ratio,
               reference_temperature,
               young_modulus_temperature_coefficient,
               poisson_ratio_temperature_coefficient](const ThermoelasticFunctionInput& input,
               ElasticPropertyOutput& output) {
        const adlite::Scalar temperature_change = input.temperature - reference_temperature;
        output.young_modulus = young_modulus + young_modulus_temperature_coefficient * temperature_change;
        output.poisson_ratio = poisson_ratio + poisson_ratio_temperature_coefficient * temperature_change;
    };
}

EigenstrainEvaluator isotropic_thermal_expansion(const MaterialParameters& named) {
    const double thermal_expansion = named.value("thermal_expansion");
    const double reference_temperature = named.value("reference_temperature");
    return
        [thermal_expansion, reference_temperature](const ThermoelasticFunctionInput& input, SymmetricTensor3& output) {
            const adlite::Scalar value = thermal_expansion * (input.temperature - reference_temperature);
            output = {value, value, value, 0.0, 0.0, 0.0};
        };
}

EigenstrainEvaluator linear_temperature_isotropic_thermal_expansion(const MaterialParameters& named) {
    const double thermal_expansion = named.value("thermal_expansion");
    const double reference_temperature = named.value("reference_temperature");
    const double thermal_expansion_temperature_coefficient = named.value("thermal_expansion_temperature_coefficient");
    return [thermal_expansion, reference_temperature, thermal_expansion_temperature_coefficient](
               const ThermoelasticFunctionInput& input,
               SymmetricTensor3& output) {
        const adlite::Scalar temperature_change = input.temperature - reference_temperature;
        const adlite::Scalar coefficient =
            thermal_expansion + thermal_expansion_temperature_coefficient * temperature_change;
        const adlite::Scalar value = coefficient * temperature_change;
        output = {value, value, value, 0.0, 0.0, 0.0};
    };
}

CreepRateEvaluator norton_creep_rate(const MaterialParameters& named) {
    const double coefficient = named.value("coefficient");
    const double reference_stress = named.value("reference_stress");
    const double stress_exponent = named.value("stress_exponent");
    return [coefficient, reference_stress, stress_exponent](const CreepRateInput& input) {
        if (!(input.equivalent_stress.value() > 0.0))
            return adlite::Scalar(0.0);
        return coefficient * adlite::pow(input.equivalent_stress / reference_stress, stress_exponent);
    };
}

CreepRateEvaluator linear_temperature_norton_creep_rate(const MaterialParameters& named) {
    const double coefficient = named.value("coefficient");
    const double reference_stress = named.value("reference_stress");
    const double stress_exponent = named.value("stress_exponent");
    const double reference_temperature = named.value("reference_temperature");
    const double coefficient_temperature_coefficient = named.value("coefficient_temperature_coefficient");
    const double reference_stress_temperature_coefficient = named.value("reference_stress_temperature_coefficient");
    const double stress_exponent_temperature_coefficient = named.value("stress_exponent_temperature_coefficient");
    return [coefficient,
               reference_stress,
               stress_exponent,
               reference_temperature,
               coefficient_temperature_coefficient,
               reference_stress_temperature_coefficient,
               stress_exponent_temperature_coefficient](const CreepRateInput& input) {
        const adlite::Scalar temperature_change = input.temperature - reference_temperature;
        const adlite::Scalar active_coefficient =
            coefficient + coefficient_temperature_coefficient * temperature_change;
        const adlite::Scalar active_reference_stress =
            reference_stress + reference_stress_temperature_coefficient * temperature_change;
        const adlite::Scalar active_stress_exponent =
            stress_exponent + stress_exponent_temperature_coefficient * temperature_change;
        if (!(input.equivalent_stress.value() > 0.0))
            return adlite::Scalar(0.0);
        return active_coefficient
               * adlite::pow(input.equivalent_stress / active_reference_stress, active_stress_exponent);
    };
}

PlasticFlowStressEvaluator linear_isotropic_flow_stress(const MaterialParameters& named) {
    const double yield_stress = named.value("yield_stress");
    const double hardening_modulus = named.value("hardening_modulus");
    return [yield_stress, hardening_modulus](const PlasticFlowStressInput& input) {
        return yield_stress + hardening_modulus * input.equivalent_plastic_strain;
    };
}

PlasticFlowStressEvaluator linear_temperature_isotropic_flow_stress(const MaterialParameters& named) {
    const double yield_stress = named.value("yield_stress");
    const double hardening_modulus = named.value("hardening_modulus");
    const double reference_temperature = named.value("reference_temperature");
    const double yield_stress_temperature_coefficient = named.value("yield_stress_temperature_coefficient");
    const double hardening_temperature_coefficient = named.value("hardening_temperature_coefficient");
    return [yield_stress,
               hardening_modulus,
               reference_temperature,
               yield_stress_temperature_coefficient,
               hardening_temperature_coefficient](const PlasticFlowStressInput& input) {
        const adlite::Scalar temperature_change = input.temperature - reference_temperature;
        const adlite::Scalar active_yield_stress =
            yield_stress + yield_stress_temperature_coefficient * temperature_change;
        const adlite::Scalar active_hardening =
            hardening_modulus + hardening_temperature_coefficient * temperature_change;
        return active_yield_stress + active_hardening * input.equivalent_plastic_strain;
    };
}
} // namespace

std::uint64_t MaterialFunctionSet::signature() const noexcept {
    std::uint64_t hash = hashing::fnv1a_offset;
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

void MaterialFunctionRegistry::add_registration(std::string name,
    std::vector<MaterialParameterDefinition> parameters,
    std::uint32_t version,
    Category category,
    Function function,
    bool has_function,
    const char* category_name) {
    validate_registration(name, parameters, has_function);
    if (version == 0)
        throw std::invalid_argument(std::string(category_name) + " material function version must be positive");
    if (std::any_of(_registrations.begin(), _registrations.end(), [&name, category](const Registration& entry) {
            return entry.category == category && entry.name == name;
        }))
        throw std::invalid_argument("Duplicate " + std::string(category_name) + " material function '" + name + "'");
    _registrations.push_back({std::move(name), version, std::move(parameters), category, function});
}

const MaterialFunctionRegistry::Registration& MaterialFunctionRegistry::find_registration(Category category,
    const std::string& name,
    const char* category_name) const {
    const auto found = std::find_if(_registrations.begin(),
        _registrations.end(),
        [&name, category](const Registration& entry) { return entry.category == category && entry.name == name; });
    if (found == _registrations.end())
        throw std::invalid_argument("Unknown " + std::string(category_name) + " material function '" + name + "'");
    return *found;
}

void MaterialFunctionRegistry::add_thermal(std::string name,
    std::vector<MaterialParameterDefinition> parameters,
    ThermalPropertyBinder function,
    std::uint32_t version) {
    add_registration(std::move(name),
        std::move(parameters),
        version,
        Category::thermal,
        function,
        function != nullptr,
        "thermal");
}

void MaterialFunctionRegistry::add_elasticity(std::string name,
    std::vector<MaterialParameterDefinition> parameters,
    ElasticPropertyBinder function,
    std::uint32_t version) {
    add_registration(std::move(name),
        std::move(parameters),
        version,
        Category::elasticity,
        function,
        function != nullptr,
        "elasticity");
}

void MaterialFunctionRegistry::add_eigenstrain(std::string name,
    std::vector<MaterialParameterDefinition> parameters,
    EigenstrainBinder function,
    std::uint32_t version) {
    add_registration(std::move(name),
        std::move(parameters),
        version,
        Category::eigenstrain,
        function,
        function != nullptr,
        "eigenstrain");
}

void MaterialFunctionRegistry::add_creep(std::string name,
    std::vector<MaterialParameterDefinition> parameters,
    CreepRateBinder function,
    std::uint32_t version) {
    add_registration(std::move(name),
        std::move(parameters),
        version,
        Category::creep,
        function,
        function != nullptr,
        "creep");
}

void MaterialFunctionRegistry::add_plasticity(std::string name,
    std::vector<MaterialParameterDefinition> parameters,
    PlasticFlowStressBinder function,
    std::uint32_t version) {
    add_registration(std::move(name),
        std::move(parameters),
        version,
        Category::plasticity,
        function,
        function != nullptr,
        "plasticity");
}

ThermalFunctionInstance MaterialFunctionRegistry::bind_thermal(const std::string& name,
    std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::thermal, name, "thermal");
    MaterialParameters bound(ordered_values(name, entry.parameters, values));
    return {name, entry.version, bound, std::get<ThermalPropertyBinder>(entry.function)(bound)};
}

ElasticFunctionInstance MaterialFunctionRegistry::bind_elasticity(const std::string& name,
    std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::elasticity, name, "elasticity");
    MaterialParameters bound(ordered_values(name, entry.parameters, values));
    return {name, entry.version, bound, std::get<ElasticPropertyBinder>(entry.function)(bound)};
}

EigenstrainFunctionInstance MaterialFunctionRegistry::bind_eigenstrain(const std::string& instance_name,
    const std::string& name,
    std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::eigenstrain, name, "eigenstrain");
    MaterialParameters bound(ordered_values(name, entry.parameters, values));
    return {instance_name, name, entry.version, bound, std::get<EigenstrainBinder>(entry.function)(bound)};
}

CreepFunctionInstance MaterialFunctionRegistry::bind_creep(const std::string& name,
    std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::creep, name, "creep");
    MaterialParameters bound(ordered_values(name, entry.parameters, values));
    CreepBuiltinParameters builtin;
    if (name == "norton" || name == "linear_temperature_norton") {
        builtin.kind = name == "norton" ? CreepBuiltinParameters::Kind::norton
                                        : CreepBuiltinParameters::Kind::linear_temperature_norton;
        builtin.coefficient = bound.value("coefficient");
        builtin.reference_stress = bound.value("reference_stress");
        builtin.stress_exponent = bound.value("stress_exponent");
        if (name == "linear_temperature_norton") {
            builtin.reference_temperature = bound.value("reference_temperature");
            builtin.coefficient_temperature_coefficient = bound.value("coefficient_temperature_coefficient");
            builtin.reference_stress_temperature_coefficient = bound.value("reference_stress_temperature_coefficient");
            builtin.stress_exponent_temperature_coefficient = bound.value("stress_exponent_temperature_coefficient");
        }
    }
    return {name, entry.version, bound, std::get<CreepRateBinder>(entry.function)(bound), builtin};
}

PlasticFunctionInstance MaterialFunctionRegistry::bind_plasticity(const std::string& name,
    std::vector<MaterialParameterValue> values) const {
    const Registration& entry = find_registration(Category::plasticity, name, "plasticity");
    MaterialParameters bound(ordered_values(name, entry.parameters, values));
    PlasticBuiltinParameters builtin;
    if (name == "linear_isotropic_hardening" || name == "linear_temperature_isotropic_hardening") {
        builtin.kind = name == "linear_isotropic_hardening"
                           ? PlasticBuiltinParameters::Kind::linear_isotropic_hardening
                           : PlasticBuiltinParameters::Kind::linear_temperature_isotropic_hardening;
        builtin.yield_stress = bound.value("yield_stress");
        builtin.hardening_modulus = bound.value("hardening_modulus");
        if (name == "linear_temperature_isotropic_hardening") {
            builtin.reference_temperature = bound.value("reference_temperature");
            builtin.yield_stress_temperature_coefficient = bound.value("yield_stress_temperature_coefficient");
            builtin.hardening_temperature_coefficient = bound.value("hardening_temperature_coefficient");
        }
    }
    return {name, entry.version, bound, std::get<PlasticFlowStressBinder>(entry.function)(bound), builtin};
}

MaterialFunctionRegistry make_builtin_material_function_registry() {
    MaterialFunctionRegistry registry;
    registry.add_thermal("constant_thermophysical",
        {{"conductivity", "W/(m*K)"}, {"density", "kg/m^3"}, {"specific_heat", "J/(kg*K)"}},
        &constant_thermophysical);
    registry.add_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", "W/m"},
            {"conductivity_constant", "W/(m*K)"},
            {"density", "kg/m^3"},
            {"specific_heat", "J/(kg*K)"}},
        &inverse_temperature_thermophysical);
    registry.add_thermal("linear_temperature_thermophysical",
        {{"conductivity", "W/(m*K)"},
            {"density", "kg/m^3"},
            {"specific_heat", "J/(kg*K)"},
            {"reference_temperature", "K"},
            {"conductivity_temperature_coefficient", "W/(m*K^2)"},
            {"density_temperature_coefficient", "kg/(m^3*K)"},
            {"specific_heat_temperature_coefficient", "J/(kg*K^2)"}},
        &linear_temperature_thermophysical);
    registry.add_elasticity("constant_isotropic",
        {{"young_modulus", "Pa"}, {"poisson_ratio", "1"}},
        &constant_isotropic_elasticity);
    registry.add_elasticity("linear_temperature_isotropic",
        {{"young_modulus", "Pa"},
            {"poisson_ratio", "1"},
            {"reference_temperature", "K"},
            {"young_modulus_temperature_coefficient", "Pa/K"},
            {"poisson_ratio_temperature_coefficient", "1/K"}},
        &linear_temperature_isotropic_elasticity);
    registry.add_eigenstrain("isotropic_thermal_expansion",
        {{"thermal_expansion", "1/K"}, {"reference_temperature", "K"}},
        &isotropic_thermal_expansion);
    registry.add_eigenstrain("linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", "1/K"},
            {"reference_temperature", "K"},
            {"thermal_expansion_temperature_coefficient", "1/K^2"}},
        &linear_temperature_isotropic_thermal_expansion);
    registry.add_creep("norton",
        {{"coefficient", "1/s"}, {"reference_stress", "Pa"}, {"stress_exponent", "1"}},
        &norton_creep_rate);
    registry.add_creep("linear_temperature_norton",
        {{"coefficient", "1/s"},
            {"reference_stress", "Pa"},
            {"stress_exponent", "1"},
            {"reference_temperature", "K"},
            {"coefficient_temperature_coefficient", "1/(s*K)"},
            {"reference_stress_temperature_coefficient", "Pa/K"},
            {"stress_exponent_temperature_coefficient", "1/K"}},
        &linear_temperature_norton_creep_rate);
    registry.add_plasticity("linear_isotropic_hardening",
        {{"yield_stress", "Pa"}, {"hardening_modulus", "Pa"}},
        &linear_isotropic_flow_stress);
    registry.add_plasticity("linear_temperature_isotropic_hardening",
        {{"yield_stress", "Pa"},
            {"hardening_modulus", "Pa"},
            {"reference_temperature", "K"},
            {"yield_stress_temperature_coefficient", "Pa/K"},
            {"hardening_temperature_coefficient", "Pa/K"}},
        &linear_temperature_isotropic_flow_stress);
    return registry;
}

} // namespace fuelsim
