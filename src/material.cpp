#include "fuelsim/material.hpp"
#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/material_functions.hpp"

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
        for (std::size_t previous = 0; previous < index; ++previous) {
            if (_values[previous].name == _values[index].name)
                throw std::invalid_argument("Duplicate material parameter '" + _values[index].name + "'");
        }
    }
}

std::size_t MaterialParameters::size() const noexcept {
    return _values.size();
}

double MaterialParameters::value(const std::string& name) const {
    const auto found = std::find_if(_values.begin(), _values.end(), [&name](const MaterialParameterValue& parameter) {
        return parameter.name == name;
    });
    if (found == _values.end())
        throw std::out_of_range("Unknown material parameter '" + name + "'");
    return found->value;
}

double MaterialParameters::value(std::size_t index) const {
    return _values.at(index).value;
}

const std::vector<MaterialParameterValue>& MaterialParameters::values() const noexcept {
    return _values;
}

namespace {

constexpr std::uint64_t material_signature_offset = 14695981039346656037ULL;
constexpr std::uint64_t material_signature_prime = 1099511628211ULL;

void material_hash_bytes(std::uint64_t& hash, const void* data, std::size_t size) noexcept {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= bytes[index];
        hash *= material_signature_prime;
    }
}

void material_hash_string(std::uint64_t& hash, const std::string& value) noexcept {
    material_hash_bytes(hash, value.data(), value.size());
    const unsigned char terminator = 0;
    material_hash_bytes(hash, &terminator, 1);
}

void material_hash_u32(std::uint64_t& hash, std::uint32_t value) noexcept {
    material_hash_bytes(hash, &value, sizeof(value));
}

void material_hash_double(std::uint64_t& hash, double value) noexcept {
    material_hash_bytes(hash, &value, sizeof(value));
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

void validate_registration(const std::string& name, const std::vector<MaterialParameterDefinition>& parameters,
                           bool has_function) {
    if (name.empty())
        throw std::invalid_argument("Material function name must not be empty");
    if (!has_function)
        throw std::invalid_argument("Material function '" + name + "' must not be null");
    for (std::size_t index = 0; index < parameters.size(); ++index) {
        if (parameters[index].name.empty())
            throw std::invalid_argument("Material function '" + name + "' has an empty parameter name");
        if (parameters[index].unit.empty())
            throw std::invalid_argument("Material function '" + name + "' parameter '" + parameters[index].name +
                                        "' has an empty SI unit");
        for (std::size_t previous = 0; previous < index; ++previous) {
            if (parameters[previous].name == parameters[index].name)
                throw std::invalid_argument("Material function '" + name + "' repeats parameter '" +
                                            parameters[index].name + "'");
        }
    }
}

void validate_bound_values(const std::string& function_name,
                           const std::vector<MaterialParameterDefinition>& definitions,
                           const std::vector<MaterialParameterValue>& values) {
    for (std::size_t index = 0; index < values.size(); ++index) {
        const MaterialParameterValue& value = values[index];
        const bool known =
            std::any_of(definitions.begin(), definitions.end(),
                        [&](const MaterialParameterDefinition& definition) { return definition.name == value.name; });
        if (!known)
            throw std::invalid_argument("Material function '" + function_name + "' does not accept parameter '" +
                                        value.name + "'");
        for (std::size_t previous = 0; previous < index; ++previous) {
            if (values[previous].name == value.name)
                throw std::invalid_argument("Material function '" + function_name + "' repeats parameter '" +
                                            value.name + "'");
        }
    }
    for (const MaterialParameterDefinition& definition : definitions) {
        const auto found = std::find_if(values.begin(), values.end(), [&](const MaterialParameterValue& value) {
            return value.name == definition.name;
        });
        if (found == values.end())
            throw std::invalid_argument("Material function '" + function_name + "' is missing parameter '" +
                                        definition.name + "'");
    }
    if (values.size() != definitions.size())
        throw std::invalid_argument("Material function '" + function_name + "' parameter count does not match");
}

std::vector<MaterialParameterValue> ordered_values(const std::string& function_name,
                                                   const std::vector<MaterialParameterDefinition>& definitions,
                                                   const std::vector<MaterialParameterValue>& values) {
    validate_bound_values(function_name, definitions, values);
    std::vector<MaterialParameterValue> result;
    result.reserve(definitions.size());
    for (const MaterialParameterDefinition& definition : definitions) {
        const auto found = std::find_if(values.begin(), values.end(), [&](const MaterialParameterValue& value) {
            return value.name == definition.name;
        });
        result.push_back(*found);
    }
    return result;
}

void constant_thermophysical(const ThermalPropertyInput& input, ThermalPropertyOutput& output) {
    output.conductivity = input.parameters->value("conductivity");
    output.density = input.parameters->value("density");
    output.specific_heat = input.parameters->value("specific_heat");
}

void inverse_temperature_thermophysical(const ThermalPropertyInput& input, ThermalPropertyOutput& output) {
    output.conductivity = input.parameters->value("conductivity_inverse_temperature") / input.temperature +
                          input.parameters->value("conductivity_constant");
    output.density = input.parameters->value("density");
    output.specific_heat = input.parameters->value("specific_heat");
}

void constant_isotropic_elasticity(const ElasticPropertyInput& input, ElasticPropertyOutput& output) {
    output.young_modulus = input.parameters->value("young_modulus");
    output.poisson_ratio = input.parameters->value("poisson_ratio");
}

void linear_temperature_isotropic_elasticity(const ElasticPropertyInput& input, ElasticPropertyOutput& output) {
    const adlite::Scalar temperature_change = input.temperature - input.parameters->value("reference_temperature");
    output.young_modulus = input.parameters->value("young_modulus") +
                           input.parameters->value("young_modulus_temperature_coefficient") * temperature_change;
    output.poisson_ratio = input.parameters->value("poisson_ratio") +
                           input.parameters->value("poisson_ratio_temperature_coefficient") * temperature_change;
}

void isotropic_thermal_expansion(const EigenstrainInput& input, AxisymmetricStrain& output) {
    const adlite::Scalar value = input.parameters->value("thermal_expansion") *
                                 (input.temperature - input.parameters->value("reference_temperature"));
    output = {value, value, value, 0.0};
}

void linear_temperature_isotropic_thermal_expansion(const EigenstrainInput& input, AxisymmetricStrain& output) {
    const adlite::Scalar temperature_change = input.temperature - input.parameters->value("reference_temperature");
    const adlite::Scalar coefficient =
        input.parameters->value("thermal_expansion") +
        input.parameters->value("thermal_expansion_temperature_coefficient") * temperature_change;
    const adlite::Scalar value = coefficient * temperature_change;
    output = {value, value, value, 0.0};
}

adlite::Scalar norton_creep_rate(const CreepRateInput& input) {
    if (!(input.equivalent_stress.value() > 0.0))
        return 0.0;
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
    if (!(input.equivalent_stress.value() > 0.0))
        return 0.0;
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

bool MaterialFunctionSet::has_creep() const noexcept {
    return creep.function != nullptr;
}

bool MaterialFunctionSet::has_plasticity() const noexcept {
    return plasticity.function != nullptr;
}

std::uint64_t MaterialFunctionSet::signature() const noexcept {
    std::uint64_t hash = material_signature_offset;
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

void MaterialFunctionRegistry::add_thermal(std::string name, std::vector<MaterialParameterDefinition> parameters,
                                           ThermalPropertyFunction function, std::uint32_t version) {
    validate_registration(name, parameters, function != nullptr);
    if (version == 0)
        throw std::invalid_argument("Thermal material function version must be positive");
    if (std::any_of(_thermal.begin(), _thermal.end(),
                    [&](const ThermalRegistration& entry) { return entry.name == name; }))
        throw std::invalid_argument("Duplicate thermal material function '" + name + "'");
    _thermal.push_back({std::move(name), version, std::move(parameters), function});
}

void MaterialFunctionRegistry::add_elasticity(std::string name, std::vector<MaterialParameterDefinition> parameters,
                                              ElasticPropertyFunction function, std::uint32_t version) {
    validate_registration(name, parameters, function != nullptr);
    if (version == 0)
        throw std::invalid_argument("Elasticity material function version must be positive");
    if (std::any_of(_elasticity.begin(), _elasticity.end(),
                    [&](const ElasticRegistration& entry) { return entry.name == name; }))
        throw std::invalid_argument("Duplicate elasticity material function '" + name + "'");
    _elasticity.push_back({std::move(name), version, std::move(parameters), function});
}

void MaterialFunctionRegistry::add_eigenstrain(std::string name, std::vector<MaterialParameterDefinition> parameters,
                                               EigenstrainFunction function, std::uint32_t version) {
    validate_registration(name, parameters, function != nullptr);
    if (version == 0)
        throw std::invalid_argument("Eigenstrain material function version must be positive");
    if (std::any_of(_eigenstrain.begin(), _eigenstrain.end(),
                    [&](const EigenstrainRegistration& entry) { return entry.name == name; }))
        throw std::invalid_argument("Duplicate eigenstrain material function '" + name + "'");
    _eigenstrain.push_back({std::move(name), version, std::move(parameters), function});
}

void MaterialFunctionRegistry::add_creep(std::string name, std::vector<MaterialParameterDefinition> parameters,
                                         CreepRateFunction function, std::uint32_t version) {
    validate_registration(name, parameters, function != nullptr);
    if (version == 0)
        throw std::invalid_argument("Creep material function version must be positive");
    if (std::any_of(_creep.begin(), _creep.end(), [&](const CreepRegistration& entry) { return entry.name == name; }))
        throw std::invalid_argument("Duplicate creep material function '" + name + "'");
    _creep.push_back({std::move(name), version, std::move(parameters), function});
}

void MaterialFunctionRegistry::add_plasticity(std::string name, std::vector<MaterialParameterDefinition> parameters,
                                              PlasticFlowStressFunction function, std::uint32_t version) {
    validate_registration(name, parameters, function != nullptr);
    if (version == 0)
        throw std::invalid_argument("Plasticity material function version must be positive");
    if (std::any_of(_plasticity.begin(), _plasticity.end(),
                    [&](const PlasticRegistration& entry) { return entry.name == name; }))
        throw std::invalid_argument("Duplicate plasticity material function '" + name + "'");
    _plasticity.push_back({std::move(name), version, std::move(parameters), function});
}

const std::vector<MaterialParameterDefinition>&
MaterialFunctionRegistry::thermal_parameters(const std::string& name) const {
    const auto found = std::find_if(_thermal.begin(), _thermal.end(),
                                    [&](const ThermalRegistration& entry) { return entry.name == name; });
    if (found == _thermal.end())
        throw std::invalid_argument("Unknown thermal material function '" + name + "'");
    return found->parameters;
}

const std::vector<MaterialParameterDefinition>&
MaterialFunctionRegistry::elasticity_parameters(const std::string& name) const {
    const auto found = std::find_if(_elasticity.begin(), _elasticity.end(),
                                    [&](const ElasticRegistration& entry) { return entry.name == name; });
    if (found == _elasticity.end())
        throw std::invalid_argument("Unknown elasticity material function '" + name + "'");
    return found->parameters;
}

const std::vector<MaterialParameterDefinition>&
MaterialFunctionRegistry::eigenstrain_parameters(const std::string& name) const {
    const auto found = std::find_if(_eigenstrain.begin(), _eigenstrain.end(),
                                    [&](const EigenstrainRegistration& entry) { return entry.name == name; });
    if (found == _eigenstrain.end())
        throw std::invalid_argument("Unknown eigenstrain material function '" + name + "'");
    return found->parameters;
}

const std::vector<MaterialParameterDefinition>&
MaterialFunctionRegistry::creep_parameters(const std::string& name) const {
    const auto found =
        std::find_if(_creep.begin(), _creep.end(), [&](const CreepRegistration& entry) { return entry.name == name; });
    if (found == _creep.end())
        throw std::invalid_argument("Unknown creep material function '" + name + "'");
    return found->parameters;
}

const std::vector<MaterialParameterDefinition>&
MaterialFunctionRegistry::plasticity_parameters(const std::string& name) const {
    const auto found = std::find_if(_plasticity.begin(), _plasticity.end(),
                                    [&](const PlasticRegistration& entry) { return entry.name == name; });
    if (found == _plasticity.end())
        throw std::invalid_argument("Unknown plasticity material function '" + name + "'");
    return found->parameters;
}

ThermalFunctionInstance MaterialFunctionRegistry::bind_thermal(const std::string& name,
                                                               std::vector<MaterialParameterValue> values) const {
    const auto found = std::find_if(_thermal.begin(), _thermal.end(),
                                    [&](const ThermalRegistration& entry) { return entry.name == name; });
    if (found == _thermal.end())
        throw std::invalid_argument("Unknown thermal material function '" + name + "'");
    return {name, found->version, MaterialParameters(ordered_values(name, found->parameters, values)), found->function};
}

ElasticFunctionInstance MaterialFunctionRegistry::bind_elasticity(const std::string& name,
                                                                  std::vector<MaterialParameterValue> values) const {
    const auto found = std::find_if(_elasticity.begin(), _elasticity.end(),
                                    [&](const ElasticRegistration& entry) { return entry.name == name; });
    if (found == _elasticity.end())
        throw std::invalid_argument("Unknown elasticity material function '" + name + "'");
    return {name, found->version, MaterialParameters(ordered_values(name, found->parameters, values)), found->function};
}

EigenstrainFunctionInstance
MaterialFunctionRegistry::bind_eigenstrain(const std::string& instance_name, const std::string& name,
                                           std::vector<MaterialParameterValue> values) const {
    const auto found = std::find_if(_eigenstrain.begin(), _eigenstrain.end(),
                                    [&](const EigenstrainRegistration& entry) { return entry.name == name; });
    if (found == _eigenstrain.end())
        throw std::invalid_argument("Unknown eigenstrain material function '" + name + "'");
    return {instance_name, name, found->version, MaterialParameters(ordered_values(name, found->parameters, values)),
            found->function};
}

CreepFunctionInstance MaterialFunctionRegistry::bind_creep(const std::string& name,
                                                           std::vector<MaterialParameterValue> values) const {
    const auto found =
        std::find_if(_creep.begin(), _creep.end(), [&](const CreepRegistration& entry) { return entry.name == name; });
    if (found == _creep.end())
        throw std::invalid_argument("Unknown creep material function '" + name + "'");
    return {name, found->version, MaterialParameters(ordered_values(name, found->parameters, values)), found->function};
}

PlasticFunctionInstance MaterialFunctionRegistry::bind_plasticity(const std::string& name,
                                                                  std::vector<MaterialParameterValue> values) const {
    const auto found = std::find_if(_plasticity.begin(), _plasticity.end(),
                                    [&](const PlasticRegistration& entry) { return entry.name == name; });
    if (found == _plasticity.end())
        throw std::invalid_argument("Unknown plasticity material function '" + name + "'");
    return {name, found->version, MaterialParameters(ordered_values(name, found->parameters, values)), found->function};
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
    registry.add_elasticity("constant_isotropic", {{"young_modulus", "Pa"}, {"poisson_ratio", "1"}},
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
    registry.add_creep("norton", {{"coefficient", "1/s"}, {"reference_stress", "Pa"}, {"stress_exponent", "1"}},
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
    registry.add_plasticity("linear_isotropic_hardening", {{"yield_stress", "Pa"}, {"hardening_modulus", "Pa"}},
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

AxisymmetricStress rotate_axisymmetric_tensor(const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation) {
    return {
        rotation.rr * rotation.rr * tensor.rr + rotation.rz * rotation.rz * tensor.zz +
            2.0 * rotation.rr * rotation.rz * tensor.rz,
        rotation.zr * rotation.zr * tensor.rr + rotation.zz * rotation.zz * tensor.zz +
            2.0 * rotation.zr * rotation.zz * tensor.rz,
        rotation.hoop * rotation.hoop * tensor.hoop,
        rotation.rr * rotation.zr * tensor.rr + rotation.rz * rotation.zz * tensor.zz +
            (rotation.rr * rotation.zz + rotation.rz * rotation.zr) * tensor.rz,
    };
}

IsotropicThermoelasticMaterial::IsotropicThermoelasticMaterial(ThermoelasticProperties properties)
    : _properties(properties), _lame_lambda(0.0), _shear_modulus(0.0), _temperature_dependent(false) {
    if (_properties.functions) {
        if (_properties.functions->thermal.function == nullptr)
            throw std::invalid_argument("Material function set requires a thermal function");
        if (_properties.functions->elasticity.function == nullptr)
            throw std::invalid_argument("Material function set requires an elasticity function");
        return;
    }
    _lame_lambda = _properties.young_modulus * _properties.poisson_ratio /
                   ((1.0 + _properties.poisson_ratio) * (1.0 - 2.0 * _properties.poisson_ratio));
    _shear_modulus = _properties.young_modulus / (2.0 * (1.0 + _properties.poisson_ratio));
    _temperature_dependent = _properties.young_modulus_temperature_coefficient != 0.0 ||
                             _properties.poisson_ratio_temperature_coefficient != 0.0 ||
                             _properties.thermal_expansion_temperature_coefficient != 0.0;
}

const ThermoelasticProperties& IsotropicThermoelasticMaterial::properties() const noexcept {
    return _properties;
}

adlite::Scalar IsotropicThermoelasticMaterial::conductivity(const adlite::Scalar& temperature, double time,
                                                            double radius, double axial_coordinate) const {
    if (!std::isfinite(temperature.value()) || !(temperature.value() > 0.0))
        throw std::domain_error("Thermoelastic material temperature must be finite and "
                                "positive");
    if (_properties.functions) {
        ThermalPropertyOutput output{};
        const ThermalFunctionInstance& instance = _properties.functions->thermal;
        instance.function({temperature, time, radius, axial_coordinate, &instance.parameters}, output);
        if (!std::isfinite(output.conductivity.value()) || !(output.conductivity.value() > 0.0))
            throw std::domain_error("Thermal material function conductivity must be finite and positive");
        return output.conductivity;
    }
    return _properties.conductivity_inverse_temperature / temperature + _properties.conductivity_offset;
}

adlite::Scalar IsotropicThermoelasticMaterial::heat_capacity(const adlite::Scalar& temperature, double time,
                                                             double radius, double axial_coordinate) const {
    if (_properties.functions) {
        ThermalPropertyOutput output{};
        const ThermalFunctionInstance& instance = _properties.functions->thermal;
        instance.function({temperature, time, radius, axial_coordinate, &instance.parameters}, output);
        if (!std::isfinite(output.density.value()) || !(output.density.value() > 0.0) ||
            !std::isfinite(output.specific_heat.value()) || !(output.specific_heat.value() > 0.0))
            throw std::domain_error("Thermal material function density and specific_heat must be finite and positive");
        return output.density * output.specific_heat;
    }
    throw std::logic_error("Thermoelastic material heat capacity is only available for a material function set");
}

ActiveThermoelasticProperties IsotropicThermoelasticMaterial::active_properties(const adlite::Scalar& temperature,
                                                                                double time, double radius,
                                                                                double axial_coordinate) const {
    if (!std::isfinite(temperature.value()))
        throw std::domain_error("Thermoelastic material temperature must be finite");
    if (_properties.functions) {
        ElasticPropertyOutput output{};
        const ElasticFunctionInstance& instance = _properties.functions->elasticity;
        instance.function({temperature, time, radius, axial_coordinate, &instance.parameters}, output);
        if (!std::isfinite(output.young_modulus.value()) || !(output.young_modulus.value() > 0.0))
            throw std::domain_error("Elasticity material function young_modulus must be finite and positive");
        if (!std::isfinite(output.poisson_ratio.value()) ||
            !(output.poisson_ratio.value() > -1.0 && output.poisson_ratio.value() < 0.5))
            throw std::domain_error("Elasticity material function poisson_ratio must lie between -1 and 0.5");
        const adlite::Scalar lame_lambda = output.young_modulus * output.poisson_ratio /
                                           ((1.0 + output.poisson_ratio) * (1.0 - 2.0 * output.poisson_ratio));
        const adlite::Scalar shear_modulus = output.young_modulus / (2.0 * (1.0 + output.poisson_ratio));
        return {output.young_modulus, output.poisson_ratio, 0.0, lame_lambda, shear_modulus};
    }
    if (!_temperature_dependent)
        return {_properties.young_modulus, _properties.poisson_ratio, _properties.thermal_expansion, _lame_lambda,
                _shear_modulus};
    const adlite::Scalar temperature_change = temperature - _properties.reference_temperature;
    const adlite::Scalar young_modulus =
        _properties.young_modulus + _properties.young_modulus_temperature_coefficient * temperature_change;
    const adlite::Scalar poisson_ratio =
        _properties.poisson_ratio + _properties.poisson_ratio_temperature_coefficient * temperature_change;
    const adlite::Scalar thermal_expansion =
        _properties.thermal_expansion + _properties.thermal_expansion_temperature_coefficient * temperature_change;
    if (!std::isfinite(young_modulus.value()) || !(young_modulus.value() > 0.0))
        throw std::domain_error("Thermoelastic material active young_modulus must be positive");
    if (!std::isfinite(poisson_ratio.value()) || !(poisson_ratio.value() > -1.0 && poisson_ratio.value() < 0.5))
        throw std::domain_error("Thermoelastic material active poisson_ratio must lie between "
                                "-1 and 0.5");
    if (!std::isfinite(thermal_expansion.value()))
        throw std::domain_error("Thermoelastic material active thermal_expansion must be "
                                "finite");
    const adlite::Scalar lame_lambda =
        young_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const adlite::Scalar shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio));
    return {young_modulus, poisson_ratio, thermal_expansion, lame_lambda, shear_modulus};
}

AxisymmetricStrain IsotropicThermoelasticMaterial::eigenstrain(const adlite::Scalar& temperature, double time,
                                                               double radius, double axial_coordinate) const {
    if (!_properties.functions) {
        const ActiveThermoelasticProperties active = active_properties(temperature);
        const adlite::Scalar thermal_strain =
            active.thermal_expansion * (temperature - _properties.reference_temperature);
        return {thermal_strain, thermal_strain, thermal_strain, 0.0};
    }
    AxisymmetricStrain result{};
    for (const EigenstrainFunctionInstance& instance : _properties.functions->eigenstrains) {
        AxisymmetricStrain value{};
        instance.function({temperature, time, radius, axial_coordinate, &instance.parameters}, value);
        if (!std::isfinite(value.rr.value()) || !std::isfinite(value.zz.value()) ||
            !std::isfinite(value.hoop.value()) || !std::isfinite(value.rz.value()))
            throw std::domain_error("Eigenstrain material function output must be finite");
        result.rr += value.rr;
        result.zz += value.zz;
        result.hoop += value.hoop;
        result.rz += value.rz;
    }
    return result;
}

namespace {

AxisymmetricStress hooke_stress(const adlite::Scalar& lame_lambda, const adlite::Scalar& shear_modulus,
                                const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
                                const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz) {
    const adlite::Scalar trace = strain_rr + strain_zz + strain_hoop;
    return {
        lame_lambda * trace + 2.0 * shear_modulus * strain_rr,
        lame_lambda * trace + 2.0 * shear_modulus * strain_zz,
        lame_lambda * trace + 2.0 * shear_modulus * strain_hoop,
        2.0 * shear_modulus * strain_rz,
    };
}

} // namespace

AxisymmetricStress IsotropicThermoelasticMaterial::stress(const adlite::Scalar& strain_rr,
                                                          const adlite::Scalar& strain_zz,
                                                          const adlite::Scalar& strain_hoop,
                                                          const adlite::Scalar& strain_rz,
                                                          const adlite::Scalar& temperature, double time, double radius,
                                                          double axial_coordinate) const {
    ActiveThermoelasticProperties active{_properties.young_modulus, _properties.poisson_ratio,
                                         _properties.thermal_expansion, _lame_lambda, _shear_modulus};
    if (_temperature_dependent || _properties.functions)
        active = active_properties(temperature, time, radius, axial_coordinate);
    const AxisymmetricStrain imposed = eigenstrain(temperature, time, radius, axial_coordinate);

    return hooke_stress(active.lame_lambda, active.shear_modulus, strain_rr - imposed.rr, strain_zz - imposed.zz,
                        strain_hoop - imposed.hoop, strain_rz - imposed.rz);
}

// Inelastic constitutive updates.
namespace {

constexpr std::size_t component_count = 4;
constexpr int maximum_creep_iterations = 100;
const double log_denorm_min = std::log(std::numeric_limits<double>::denorm_min());

struct NortonRoot final {
    double equivalent_stress;
    double trial_stress_derivative;
    double relaxed_fraction;
};

struct CoupledUpdate final {
    double equivalent_stress;
    double stress_derivative;
    double plastic_increment;
    double plastic_increment_derivative;
    double creep_increment;
    double creep_increment_derivative;
};

void validate_material_point_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < component_count; ++component) {
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component]) ||
            !std::isfinite(state.creep_strain[component]))
            throw std::domain_error("Material point strain histories must be finite");
    }

    if (!std::isfinite(state.equivalent_plastic_strain) || !(state.equivalent_plastic_strain >= 0.0))
        throw std::domain_error("Material point equivalent plastic strain must be finite and "
                                "nonnegative");
    if (!std::isfinite(state.equivalent_creep_strain) || !(state.equivalent_creep_strain >= 0.0))
        throw std::domain_error("Material point equivalent creep strain must be finite and "
                                "nonnegative");
}

MaterialPointTrialState passive_trial_state(const MaterialPointState& committed) {
    MaterialPointTrialState trial;
    std::copy(committed.elastic_strain.begin(), committed.elastic_strain.end(), trial.elastic_strain.begin());
    std::copy(committed.plastic_strain.begin(), committed.plastic_strain.end(), trial.plastic_strain.begin());
    std::copy(committed.creep_strain.begin(), committed.creep_strain.end(), trial.creep_strain.begin());
    trial.equivalent_plastic_strain = committed.equivalent_plastic_strain;
    trial.equivalent_creep_strain = committed.equivalent_creep_strain;
    return trial;
}

AxisymmetricStress tensor(const std::array<adlite::Scalar, component_count>& components) {
    return {components[0], components[1], components[2], components[3]};
}

std::array<adlite::Scalar, component_count> components(const AxisymmetricStress& tensor) {
    return {tensor.rr, tensor.zz, tensor.hoop, tensor.rz};
}

double log_add_exp(double first, double second) {
    const double maximum = std::max(first, second);
    if (std::isinf(maximum))
        return maximum;
    return maximum + std::log1p(std::exp(std::min(first, second) - maximum));
}

double second_exponential_weight(double first, double second) {
    if (second >= first)
        return 1.0 / (1.0 + std::exp(first - second));
    const double ratio = std::exp(second - first);
    return ratio / (1.0 + ratio);
}

NortonRoot solve_power_law_equivalent_stress(double driving_stress, double log_stress_coefficient, double time_step,
                                             const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0 ||
        log_stress_coefficient == -std::numeric_limits<double>::infinity())
        return {driving_stress, 1.0, 0.0};
    if (!(driving_stress > 0.0)) {
        if (creep.stress_exponent > 1.0)
            return {driving_stress, 1.0, 0.0};

        const double log_linear_coefficient = log_stress_coefficient + std::log(time_step) +
                                              std::log(creep.coefficient) - std::log(creep.reference_stress);
        if (std::isnan(log_linear_coefficient))
            throw std::overflow_error("Norton creep linear coefficient is not a number");
        if (log_linear_coefficient == -std::numeric_limits<double>::infinity())
            return {driving_stress, 1.0, 0.0};
        if (log_linear_coefficient == std::numeric_limits<double>::infinity())
            return {driving_stress, 0.0, 0.0};
        const double log_denominator = log_add_exp(0.0, log_linear_coefficient);
        return {
            driving_stress,
            std::exp(-log_denominator),
            0.0,
        };
    }

    const double log_driving_stress = std::log(driving_stress);
    const double log_coefficient = log_stress_coefficient + std::log(time_step) + std::log(creep.coefficient) -
                                   log_driving_stress +
                                   creep.stress_exponent * (log_driving_stress - std::log(creep.reference_stress));
    if (std::isnan(log_coefficient))
        throw std::overflow_error("Norton creep dimensionless local coefficient is not a number");
    if (log_coefficient == -std::numeric_limits<double>::infinity())
        return {driving_stress, 1.0, 0.0};
    if (log_coefficient == std::numeric_limits<double>::infinity())
        throw std::overflow_error("Norton creep logarithmic local coefficient overflowed");

    double log_stress_ratio = 0.0;
    if (creep.stress_exponent == 1.0) {
        log_stress_ratio = log_add_exp(0.0, log_coefficient);
    } else {
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
        const double initial_residual = log_add_exp(0.0, log_coefficient);
        if (initial_residual == 0.0)
            return {driving_stress, 1.0, 0.0};

        constexpr double log_two = 0.693147180559945309417232121458176568;
        double lower = 0.0;
        double upper = std::max(log_two, (log_coefficient + log_two) / creep.stress_exponent);
        double current = log_coefficient > 0.0 ? log_coefficient / creep.stress_exponent : std::exp(log_coefficient);
        if (!std::isfinite(upper) || !(upper > lower))
            throw std::overflow_error("Norton creep logarithmic root bracket is invalid");
        if (!std::isfinite(current) || !(current > lower) || !(current < upper))
            current = 0.5 * (lower + upper);

        bool converged = false;
        for (int iteration = 0; iteration < maximum_creep_iterations; ++iteration) {
            const double first = -current;
            const double second = log_coefficient - creep.stress_exponent * current;
            const double residual = log_add_exp(first, second);
            if (!std::isfinite(residual))
                throw std::overflow_error("Norton creep logarithmic residual is not finite");
            if (std::fabs(residual) <= tolerance) {
                converged = true;
                break;
            }

            if (residual > 0.0)
                lower = current;
            else
                upper = current;
            if (upper - lower <= tolerance * std::max(1.0, std::fabs(current))) {
                current = 0.5 * (lower + upper);
                converged = true;
                break;
            }

            const double weight = second_exponential_weight(first, second);
            const double derivative = -(1.0 + (creep.stress_exponent - 1.0) * weight);
            double candidate = std::isfinite(derivative) ? current - residual / derivative : 0.5 * (lower + upper);
            if (!std::isfinite(candidate) || !(candidate > lower) || !(candidate < upper) || candidate == current)
                candidate = 0.5 * (lower + upper);
            current = candidate;
        }
        if (!converged)
            throw std::domain_error("Norton creep logarithmic local Newton solve did not "
                                    "converge");
        log_stress_ratio = current;
    }

    const double stress_ratio = std::exp(-log_stress_ratio);
    const double relaxed_fraction = -std::expm1(-log_stress_ratio);
    const double log_equivalent_stress = log_driving_stress - log_stress_ratio;
    const double equivalent_stress = log_equivalent_stress <= log_denorm_min ? 0.0 : std::exp(log_equivalent_stress);
    const double derivative_denominator = stress_ratio + creep.stress_exponent * relaxed_fraction;
    const double driving_stress_derivative =
        derivative_denominator == 0.0 ? 0.0 : stress_ratio / derivative_denominator;
    if (!std::isfinite(equivalent_stress) || !std::isfinite(driving_stress_derivative) ||
        !(driving_stress_derivative >= 0.0) || !std::isfinite(relaxed_fraction) || !(relaxed_fraction >= 0.0) ||
        !(relaxed_fraction <= 1.0))
        throw std::overflow_error("Norton creep logarithmic root or derivative is not finite");

    return {
        equivalent_stress,
        driving_stress_derivative,
        relaxed_fraction,
    };
}

NortonRoot solve_norton_equivalent_stress(double trial_stress, double shear_modulus, double time_step,
                                          const NortonCreepProperties& creep) {
    return solve_power_law_equivalent_stress(trial_stress, std::log(3.0) + std::log(shear_modulus), time_step, creep);
}

double evaluate_creep_increment(double equivalent_stress, double time_step, const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0)
        return 0.0;
    if (!(equivalent_stress > 0.0))
        return 0.0;

    const double log_increment =
        std::log(time_step) + std::log(creep.coefficient) +
        creep.stress_exponent * (std::log(equivalent_stress) - std::log(creep.reference_stress));
    if (std::isnan(log_increment))
        throw std::overflow_error("Norton creep increment logarithm is not a number");

    const double maximum_log = std::log(std::numeric_limits<double>::max());
    if (log_increment <= log_denorm_min)
        return 0.0;
    if (!std::isfinite(log_increment) || log_increment >= maximum_log)
        throw std::overflow_error("Norton creep increment overflowed");

    const double increment = std::exp(log_increment);
    if (!std::isfinite(increment) || !(increment >= 0.0))
        throw std::overflow_error("Norton creep increment is not finite");
    return increment;
}

CoupledUpdate solve_coupled_update(double trial_stress, double shear_modulus, double time_step,
                                   const NortonCreepProperties& creep, const J2PlasticityProperties& plasticity,
                                   double committed_equivalent_plastic_strain) {
    const double current_yield_stress =
        plasticity.yield_stress + plasticity.isotropic_hardening_modulus * committed_equivalent_plastic_strain;
    if (!std::isfinite(current_yield_stress))
        throw std::overflow_error("Coupled plastic-creep current yield stress is not finite");

    const NortonRoot creep_only = solve_norton_equivalent_stress(trial_stress, shear_modulus, time_step, creep);
    const double inverse_three_shear_modulus = (1.0 / 3.0) / shear_modulus;
    if (!(creep_only.equivalent_stress > current_yield_stress)) {
        const double creep_increment = trial_stress * creep_only.relaxed_fraction * inverse_three_shear_modulus;
        const double creep_increment_derivative =
            (1.0 - creep_only.trial_stress_derivative) * inverse_three_shear_modulus;
        return {
            creep_only.equivalent_stress, creep_only.trial_stress_derivative, 0.0, 0.0, creep_increment,
            creep_increment_derivative,
        };
    }

    const double hardening = plasticity.isotropic_hardening_modulus;
    if (hardening == 0.0) {
        const double creep_increment = evaluate_creep_increment(current_yield_stress, time_step, creep);
        const double plastic_increment =
            (trial_stress - current_yield_stress) * inverse_three_shear_modulus - creep_increment;
        if (!std::isfinite(plastic_increment) || !(plastic_increment > 0.0))
            throw std::domain_error("Coupled perfect-plastic update produced a nonpositive "
                                    "plastic increment");
        return {
            current_yield_stress, 0.0, plastic_increment, inverse_three_shear_modulus, creep_increment, 0.0,
        };
    }

    double hardening_weight = 0.0;
    double inverse_denominator = 0.0;
    double log_root_stress_coefficient = 0.0;
    if (hardening <= shear_modulus) {
        const double ratio = hardening / shear_modulus;
        const double scaled_denominator = 3.0 + ratio;
        hardening_weight = ratio / scaled_denominator;
        inverse_denominator = (1.0 / scaled_denominator) / shear_modulus;
        log_root_stress_coefficient = std::log(hardening) + std::log(3.0 / scaled_denominator);
    } else {
        const double ratio = shear_modulus / hardening;
        const double scaled_denominator = 1.0 + 3.0 * ratio;
        hardening_weight = 1.0 / scaled_denominator;
        inverse_denominator = (1.0 / scaled_denominator) / hardening;
        log_root_stress_coefficient = std::log(shear_modulus) + std::log(3.0 / scaled_denominator);
    }
    const double root_driving_stress = current_yield_stress + hardening_weight * (trial_stress - current_yield_stress);
    if (!std::isfinite(inverse_denominator) || !(inverse_denominator > 0.0) || !std::isfinite(root_driving_stress) ||
        !(root_driving_stress > 0.0) || !std::isfinite(log_root_stress_coefficient))
        throw std::overflow_error("Coupled plastic-creep transformed root is invalid");

    const NortonRoot root =
        solve_power_law_equivalent_stress(root_driving_stress, log_root_stress_coefficient, time_step, creep);
    const double creep_increment = evaluate_creep_increment(root.equivalent_stress, time_step, creep);
    const double plastic_increment =
        (trial_stress - current_yield_stress - 3.0 * (shear_modulus * creep_increment)) * inverse_denominator;
    const double plastic_increment_derivative = root.trial_stress_derivative * inverse_denominator;
    const double stress_derivative = hardening_weight * root.trial_stress_derivative;
    const double creep_increment_derivative = (1.0 - root.trial_stress_derivative) * inverse_three_shear_modulus;
    if (!std::isfinite(root.equivalent_stress) || !std::isfinite(stress_derivative) ||
        !std::isfinite(plastic_increment) || !(plastic_increment > 0.0) ||
        !std::isfinite(plastic_increment_derivative) || !std::isfinite(creep_increment) ||
        !std::isfinite(creep_increment_derivative))
        throw std::overflow_error("Coupled plastic-creep update or derivative is not finite");

    return {
        root.equivalent_stress,       stress_derivative, plastic_increment,
        plastic_increment_derivative, creep_increment,   creep_increment_derivative,
    };
}

double zero_stress_scale_shear_derivative(double scale, double time_step, const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0 || creep.stress_exponent > 1.0)
        return 0.0;
    return -3.0 * time_step * creep.coefficient / creep.reference_stress * scale * scale;
}

struct NortonDirectPartials final {
    double creep_stress_derivative;
    double direct_coefficient;
    double direct_reference;
    double direct_exponent;
};

NortonDirectPartials norton_direct_partials(double equivalent_stress, double creep_increment,
                                            const NortonCreepProperties& creep) {
    return {
        equivalent_stress > 0.0 ? creep.stress_exponent * creep_increment / equivalent_stress : 0.0,
        creep.coefficient > 0.0 ? creep_increment / creep.coefficient : 0.0,
        -creep.stress_exponent * creep_increment / creep.reference_stress,
        equivalent_stress > 0.0 ? creep_increment * std::log(equivalent_stress / creep.reference_stress) : 0.0,
    };
}

void fill_creep_partials(double* creep, std::size_t count, double creep_stress_derivative, const double* stress,
                         const double* direct) {
    for (std::size_t parameter = 0; parameter < count; ++parameter)
        creep[parameter] = creep_stress_derivative * stress[parameter] + direct[parameter];
}

struct NortonPartials final {
    std::array<double, 5> stress{};
    std::array<double, 5> creep{};
};

NortonPartials norton_partials(double equivalent_trial_stress, double shear_modulus, double time_step,
                               const NortonCreepProperties& creep, const NortonRoot& root, double creep_increment) {
    NortonPartials result;
    const double stress = root.equivalent_stress;
    const double stress_derivative = root.trial_stress_derivative;
    const NortonDirectPartials direct = norton_direct_partials(stress, creep_increment, creep);
    result.stress = {stress_derivative, -3.0 * creep_increment * stress_derivative,
                     -3.0 * shear_modulus * direct.direct_coefficient * stress_derivative,
                     -3.0 * shear_modulus * direct.direct_reference * stress_derivative,
                     -3.0 * shear_modulus * direct.direct_exponent * stress_derivative};
    const std::array<double, 5> creep_direct = {0.0, 0.0, direct.direct_coefficient, direct.direct_reference,
                                                direct.direct_exponent};
    fill_creep_partials(result.creep.data(), result.creep.size(), direct.creep_stress_derivative, result.stress.data(),
                        creep_direct.data());
    if (time_step == 0.0 || equivalent_trial_stress == 0.0)
        result.creep.fill(0.0);
    return result;
}

adlite::Scalar compose_with_derivatives(double value, const adlite::Scalar* inputs, std::size_t input_count,
                                        const double* derivatives) {
    return adlite::compose(value, inputs, derivatives, input_count);
}

struct CoupledPartials final {
    std::array<double, 7> stress{};
    std::array<double, 7> plastic{};
    std::array<double, 7> creep{};
};

CoupledPartials coupled_partials(double equivalent_trial_stress, double shear_modulus,
                                 const NortonCreepProperties& creep, const J2PlasticityProperties& plasticity,
                                 double committed_equivalent_plastic_strain, const CoupledUpdate& update) {
    CoupledPartials result;
    const double stress = update.equivalent_stress;
    const double plastic_increment = update.plastic_increment;
    const double creep_increment = update.creep_increment;
    const NortonDirectPartials direct = norton_direct_partials(stress, creep_increment, creep);
    const std::array<double, 7> direct_creep = {
        0.0, 0.0, direct.direct_coefficient, direct.direct_reference, direct.direct_exponent, 0.0, 0.0};
    const double hardening = plasticity.isotropic_hardening_modulus;
    if (!(plastic_increment > 0.0)) {
        const NortonRoot root{stress, update.stress_derivative, 0.0};
        const NortonPartials norton =
            norton_partials(equivalent_trial_stress, shear_modulus, 1.0, creep, root, creep_increment);
        for (std::size_t parameter = 0; parameter < 5; ++parameter) {
            result.stress[parameter] = norton.stress[parameter];
            result.creep[parameter] = norton.creep[parameter];
        }
        return result;
    }
    if (hardening == 0.0) {
        result.stress[5] = 1.0;
        fill_creep_partials(result.creep.data(), result.creep.size(), direct.creep_stress_derivative,
                            result.stress.data(), direct_creep.data());
        const double inverse_three_shear = 1.0 / (3.0 * shear_modulus);
        result.plastic[0] = inverse_three_shear;
        result.plastic[1] =
            -(equivalent_trial_stress - plasticity.yield_stress) / (3.0 * shear_modulus * shear_modulus) -
            result.creep[1];
        for (std::size_t parameter = 2; parameter < result.plastic.size(); ++parameter)
            result.plastic[parameter] = -result.creep[parameter];
        result.plastic[5] -= inverse_three_shear;
        return result;
    }

    const double stress_derivative = update.stress_derivative;
    result.stress[0] = stress_derivative;
    result.stress[1] = -3.0 * (plastic_increment + creep_increment) * stress_derivative;
    result.stress[2] = -3.0 * shear_modulus * direct.direct_coefficient * stress_derivative;
    result.stress[3] = -3.0 * shear_modulus * direct.direct_reference * stress_derivative;
    result.stress[4] = -3.0 * shear_modulus * direct.direct_exponent * stress_derivative;
    result.stress[5] = 3.0 * shear_modulus / hardening * stress_derivative;
    result.stress[6] =
        3.0 * shear_modulus * (stress - plasticity.yield_stress) / (hardening * hardening) * stress_derivative;
    fill_creep_partials(result.creep.data(), result.creep.size(), direct.creep_stress_derivative, result.stress.data(),
                        direct_creep.data());
    for (std::size_t parameter = 0; parameter < result.plastic.size(); ++parameter) {
        const double yield_derivative = parameter == 5 ? 1.0 : 0.0;
        const double hardening_derivative = parameter == 6 ? 1.0 : 0.0;
        result.plastic[parameter] = (result.stress[parameter] - yield_derivative -
                                     (committed_equivalent_plastic_strain + plastic_increment) * hardening_derivative) /
                                    hardening;
    }
    return result;
}

AxisymmetricStress returned_stress(const adlite::Scalar& mean_stress,
                                   const std::array<adlite::Scalar, component_count>& deviatoric_trial,
                                   const adlite::Scalar& scale) {
    return {
        mean_stress + scale * deviatoric_trial[0],
        mean_stress + scale * deviatoric_trial[1],
        mean_stress + scale * deviatoric_trial[2],
        scale * deviatoric_trial[3],
    };
}

void update_flow_strain(std::array<adlite::Scalar, component_count>& trial,
                        const std::array<double, component_count>& committed,
                        const std::array<adlite::Scalar, component_count>& deviatoric_trial,
                        const adlite::Scalar& equivalent_trial_stress, const adlite::Scalar& increment) {
    for (std::size_t component = 0; component < component_count; ++component) {
        const adlite::Scalar flow_direction = 1.5 * deviatoric_trial[component] / equivalent_trial_stress;
        trial[component] = committed[component] + increment * flow_direction;
    }
}

void update_relaxed_creep_strain(MaterialPointTrialState& trial, const MaterialPointState& committed,
                                 const std::array<adlite::Scalar, component_count>& deviatoric_trial,
                                 const adlite::Scalar& stress_scale, const adlite::Scalar& shear_modulus) {
    for (std::size_t component = 0; component < component_count; ++component)
        trial.creep_strain[component] = committed.creep_strain[component] +
                                        (1.0 - stress_scale) * deviatoric_trial[component] / (2.0 * shear_modulus);
}

adlite::Scalar function_creep_rate(const CreepFunctionInstance& function, const adlite::Scalar& equivalent_stress,
                                   const adlite::Scalar& temperature, const adlite::Scalar& equivalent_creep_strain,
                                   double time, double radius, double axial_coordinate) {
    const adlite::Scalar rate = function.function({equivalent_stress, temperature, equivalent_creep_strain, time,
                                                   radius, axial_coordinate, &function.parameters});
    if (std::isnan(rate.value()) || rate.value() < 0.0)
        throw std::domain_error("Creep material function rate must be nonnegative and must not be NaN");
    return rate;
}

adlite::Scalar function_flow_stress(const PlasticFunctionInstance& function,
                                    const adlite::Scalar& equivalent_plastic_strain, const adlite::Scalar& temperature,
                                    double time, double radius, double axial_coordinate) {
    const adlite::Scalar stress = function.function(
        {equivalent_plastic_strain, temperature, time, radius, axial_coordinate, &function.parameters});
    if (!std::isfinite(stress.value()) || !(stress.value() > 0.0))
        throw std::domain_error("Plasticity material function flow stress must be finite and positive");
    return stress;
}

struct FunctionCreepRoot final {
    adlite::Scalar equivalent_stress;
    adlite::Scalar creep_increment;
    double zero_stress_scale;
};

double passive_creep_rate_and_derivative(const CreepFunctionInstance& function, double equivalent_stress,
                                         double temperature, double equivalent_creep_strain, double time, double radius,
                                         double axial_coordinate, double& derivative) {
    const adlite::Scalar active_stress = adlite::Scalar::independent(equivalent_stress, 0, 1);
    const adlite::Scalar rate = function_creep_rate(function, active_stress, temperature, equivalent_creep_strain, time,
                                                    radius, axial_coordinate);
    derivative = rate.is_active() ? rate.derivative(0) : 0.0;
    return rate.value();
}

FunctionCreepRoot solve_function_creep(const CreepFunctionInstance& function,
                                       const adlite::Scalar& equivalent_trial_stress,
                                       const adlite::Scalar& shear_modulus, const adlite::Scalar& temperature,
                                       double equivalent_creep_strain, double time, double time_step, double radius,
                                       double axial_coordinate) {
    double zero_derivative = 0.0;
    const double zero_rate = passive_creep_rate_and_derivative(
        function, 0.0, temperature.value(), equivalent_creep_strain, time, radius, axial_coordinate, zero_derivative);
    if (zero_rate != 0.0)
        throw std::domain_error("Creep material function rate must be zero at zero equivalent stress");
    if (!std::isfinite(zero_derivative) || zero_derivative < 0.0)
        throw std::domain_error("Creep material function must be nondecreasing at zero equivalent stress");
    const double zero_scale = 1.0 / (1.0 + 3.0 * shear_modulus.value() * time_step * zero_derivative);
    if (!(equivalent_trial_stress.value() > 0.0))
        return {0.0, 0.0, zero_scale};

    const double trial = equivalent_trial_stress.value();
    const double shear = shear_modulus.value();
    double lower = 0.0;
    double upper = trial;
    double current = trial;
    double derivative = 1.0;
    bool converged = false;
    constexpr int maximum_function_iterations = 256;
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * trial;
    for (int iteration = 0; iteration < maximum_function_iterations; ++iteration) {
        double rate_derivative = 0.0;
        const double rate =
            passive_creep_rate_and_derivative(function, current, temperature.value(), equivalent_creep_strain, time,
                                              radius, axial_coordinate, rate_derivative);
        if (std::isfinite(rate_derivative) && rate_derivative < 0.0)
            throw std::domain_error("Creep material function must be nondecreasing in equivalent stress");
        const double residual = current + 3.0 * shear * time_step * rate - trial;
        derivative = 1.0 + 3.0 * shear * time_step * rate_derivative;
        if (std::isfinite(residual) && std::isfinite(derivative) && derivative > 0.0 &&
            std::fabs(residual) <= tolerance) {
            converged = true;
            break;
        }
        if (!std::isfinite(residual) || residual > 0.0)
            upper = current;
        else
            lower = current;
        if (upper - lower <= tolerance) {
            current = 0.5 * (lower + upper);
            double final_rate_derivative = 0.0;
            (void)passive_creep_rate_and_derivative(function, current, temperature.value(), equivalent_creep_strain,
                                                    time, radius, axial_coordinate, final_rate_derivative);
            derivative = 1.0 + 3.0 * shear * time_step * final_rate_derivative;
            converged = std::isfinite(derivative) && derivative > 0.0;
            break;
        }
        double candidate = 0.5 * (lower + upper);
        if (std::isfinite(residual) && std::isfinite(derivative) && derivative > 0.0) {
            const double newton = current - residual / derivative;
            if (std::isfinite(newton) && newton > lower && newton < upper)
                candidate = newton;
        }
        current = candidate;
    }
    if (!converged)
        throw std::domain_error("Registered creep material function local solve did not converge");

    const adlite::Scalar passive_stress(current);
    const adlite::Scalar active_rate = function_creep_rate(function, passive_stress, temperature,
                                                           equivalent_creep_strain, time, radius, axial_coordinate);
    const adlite::Scalar residual =
        passive_stress + 3.0 * shear_modulus * time_step * active_rate - equivalent_trial_stress;
    const adlite::Scalar returned_stress = passive_stress - residual / derivative;
    const adlite::Scalar increment =
        time_step * function_creep_rate(function, returned_stress, temperature, equivalent_creep_strain, time, radius,
                                        axial_coordinate);
    return {returned_stress, increment, zero_scale};
}

struct FunctionCoupledUpdate final {
    adlite::Scalar equivalent_stress;
    adlite::Scalar plastic_increment;
    adlite::Scalar creep_increment;
};

FunctionCoupledUpdate solve_function_coupled(const MaterialFunctionSet& functions,
                                             const adlite::Scalar& equivalent_trial_stress,
                                             const adlite::Scalar& shear_modulus, const adlite::Scalar& temperature,
                                             double equivalent_plastic_strain, double equivalent_creep_strain,
                                             double time, double time_step, double radius, double axial_coordinate) {
    const adlite::Scalar flow_at_committed = function_flow_stress(functions.plasticity, equivalent_plastic_strain,
                                                                  temperature, time, radius, axial_coordinate);
    FunctionCreepRoot creep_only{equivalent_trial_stress, 0.0, 1.0};
    if (functions.has_creep())
        creep_only = solve_function_creep(functions.creep, equivalent_trial_stress, shear_modulus, temperature,
                                          equivalent_creep_strain, time, time_step, radius, axial_coordinate);
    if (!(creep_only.equivalent_stress.value() > flow_at_committed.value()))
        return {creep_only.equivalent_stress, 0.0, creep_only.creep_increment};

    const double trial = equivalent_trial_stress.value();
    const double shear = shear_modulus.value();
    double lower = 0.0;
    double upper = trial / (3.0 * shear);
    double current = 0.5 * upper;
    double derivative = 0.0;
    bool converged = false;
    constexpr int maximum_function_iterations = 256;
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * trial;
    for (int iteration = 0; iteration < maximum_function_iterations; ++iteration) {
        const adlite::Scalar active_increment = adlite::Scalar::independent(current, 0, 1);
        const adlite::Scalar flow =
            function_flow_stress(functions.plasticity, equivalent_plastic_strain + active_increment,
                                 temperature.value(), time, radius, axial_coordinate);
        if (flow.is_active() && flow.derivative(0) < 0.0)
            throw std::domain_error(
                "Plasticity material function flow stress must be nondecreasing in equivalent plastic strain");
        adlite::Scalar creep_increment(0.0);
        if (functions.has_creep())
            creep_increment = time_step * function_creep_rate(functions.creep, flow, temperature.value(),
                                                              equivalent_creep_strain, time, radius, axial_coordinate);
        const adlite::Scalar residual = flow + 3.0 * shear * (active_increment + creep_increment) - trial;
        derivative = residual.derivative(0);
        if (std::isfinite(residual.value()) && std::isfinite(derivative) && derivative > 0.0 &&
            std::fabs(residual.value()) <= tolerance) {
            converged = true;
            break;
        }
        if (!std::isfinite(residual.value()) || residual.value() > 0.0)
            upper = current;
        else
            lower = current;
        if (upper - lower <= tolerance / (3.0 * shear)) {
            current = 0.5 * (lower + upper);
            const adlite::Scalar final_increment = adlite::Scalar::independent(current, 0, 1);
            const adlite::Scalar final_flow =
                function_flow_stress(functions.plasticity, equivalent_plastic_strain + final_increment,
                                     temperature.value(), time, radius, axial_coordinate);
            adlite::Scalar final_creep(0.0);
            if (functions.has_creep())
                final_creep = time_step * function_creep_rate(functions.creep, final_flow, temperature.value(),
                                                              equivalent_creep_strain, time, radius, axial_coordinate);
            derivative = (final_flow + 3.0 * shear * (final_increment + final_creep) - trial).derivative(0);
            converged = std::isfinite(derivative) && derivative > 0.0;
            break;
        }
        double candidate = 0.5 * (lower + upper);
        if (std::isfinite(residual.value()) && std::isfinite(derivative) && derivative > 0.0) {
            const double newton = current - residual.value() / derivative;
            if (std::isfinite(newton) && newton > lower && newton < upper)
                candidate = newton;
        }
        current = candidate;
    }
    if (!converged)
        throw std::domain_error("Registered plasticity material function local solve did not converge");

    const adlite::Scalar passive_increment(current);
    const adlite::Scalar active_flow =
        function_flow_stress(functions.plasticity, equivalent_plastic_strain + passive_increment, temperature, time,
                             radius, axial_coordinate);
    adlite::Scalar active_creep(0.0);
    if (functions.has_creep())
        active_creep = time_step * function_creep_rate(functions.creep, active_flow, temperature,
                                                       equivalent_creep_strain, time, radius, axial_coordinate);
    const adlite::Scalar residual =
        active_flow + 3.0 * shear_modulus * (passive_increment + active_creep) - equivalent_trial_stress;
    const adlite::Scalar plastic_increment = passive_increment - residual / derivative;
    const adlite::Scalar returned_stress =
        function_flow_stress(functions.plasticity, equivalent_plastic_strain + plastic_increment, temperature, time,
                             radius, axial_coordinate);
    adlite::Scalar creep_increment(0.0);
    if (functions.has_creep())
        creep_increment = time_step * function_creep_rate(functions.creep, returned_stress, temperature,
                                                          equivalent_creep_strain, time, radius, axial_coordinate);
    return {returned_stress, plastic_increment, creep_increment};
}

} // namespace

IsotropicInelasticMaterial::IsotropicInelasticMaterial(ThermoelasticProperties thermoelastic_properties,
                                                       TransientInelasticProperties properties)
    : _thermoelastic_material(thermoelastic_properties), _properties(properties) {}

const TransientInelasticProperties& IsotropicInelasticMaterial::properties() const noexcept {
    return _properties;
}

adlite::Scalar IsotropicInelasticMaterial::conductivity(const adlite::Scalar& temperature, double time, double radius,
                                                        double axial_coordinate) const {
    return _thermoelastic_material.conductivity(temperature, time, radius, axial_coordinate);
}

adlite::Scalar IsotropicInelasticMaterial::heat_capacity(const adlite::Scalar& temperature, double time, double radius,
                                                         double axial_coordinate) const {
    if (_thermoelastic_material.properties().functions)
        return _thermoelastic_material.heat_capacity(temperature, time, radius, axial_coordinate);
    return _properties.density * _properties.specific_heat;
}

ActiveNortonCreepProperties
IsotropicInelasticMaterial::active_creep_properties(const adlite::Scalar& temperature) const {
    const double reference_temperature = _thermoelastic_material.properties().reference_temperature;
    const adlite::Scalar temperature_change = temperature - reference_temperature;
    ActiveNortonCreepProperties active{
        _properties.creep.coefficient + _properties.creep.coefficient_temperature_coefficient * temperature_change,
        _properties.creep.reference_stress +
            _properties.creep.reference_stress_temperature_coefficient * temperature_change,
        _properties.creep.stress_exponent +
            _properties.creep.stress_exponent_temperature_coefficient * temperature_change};
    if (!std::isfinite(active.coefficient.value()) || !(active.coefficient.value() >= 0.0) ||
        !std::isfinite(active.reference_stress.value()) || !(active.reference_stress.value() > 0.0) ||
        !std::isfinite(active.stress_exponent.value()) || !(active.stress_exponent.value() >= 1.0))
        throw std::domain_error("Active Norton creep properties violate their physical domain");
    return active;
}

ActiveJ2PlasticityProperties
IsotropicInelasticMaterial::active_plasticity_properties(const adlite::Scalar& temperature) const {
    const double reference_temperature = _thermoelastic_material.properties().reference_temperature;
    const adlite::Scalar temperature_change = temperature - reference_temperature;
    ActiveJ2PlasticityProperties active{
        _properties.plasticity.yield_stress +
            _properties.plasticity.yield_stress_temperature_coefficient * temperature_change,
        _properties.plasticity.isotropic_hardening_modulus +
            _properties.plasticity.hardening_temperature_coefficient * temperature_change};
    if (!std::isfinite(active.yield_stress.value()) || !(active.yield_stress.value() > 0.0) ||
        !std::isfinite(active.isotropic_hardening_modulus.value()) ||
        !(active.isotropic_hardening_modulus.value() >= 0.0))
        throw std::domain_error("Active J2 plasticity properties violate their physical domain");
    return active;
}

InelasticStressResponse IsotropicInelasticMaterial::response(
    const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz, const adlite::Scalar& strain_hoop,
    const adlite::Scalar& strain_rz, const adlite::Scalar& temperature, double time_step,
    const MaterialPointState& committed, double time, double radius, double axial_coordinate) const {
    if (!std::isfinite(strain_rr.value()) || !std::isfinite(strain_zz.value()) || !std::isfinite(strain_hoop.value()) ||
        !std::isfinite(strain_rz.value()) || !std::isfinite(temperature.value()))
        throw std::domain_error("Inelastic material strain and temperature inputs must be "
                                "finite");
    if (!std::isfinite(time_step) || !(time_step >= 0.0))
        throw std::domain_error("Inelastic material time_step must be finite and nonnegative");
    validate_material_point_state(committed);

    MaterialPointTrialState trial_state = passive_trial_state(committed);
    const ActiveThermoelasticProperties active =
        _thermoelastic_material.active_properties(temperature, time, radius, axial_coordinate);
    const AxisymmetricStrain imposed = _thermoelastic_material.eigenstrain(temperature, time, radius, axial_coordinate);
    const std::array<adlite::Scalar, component_count> elastic_strain = {
        strain_rr - imposed.rr - committed.plastic_strain[0] - committed.creep_strain[0],
        strain_zz - imposed.zz - committed.plastic_strain[1] - committed.creep_strain[1],
        strain_hoop - imposed.hoop - committed.plastic_strain[2] - committed.creep_strain[2],
        strain_rz - imposed.rz - committed.plastic_strain[3] - committed.creep_strain[3],
    };
    const AxisymmetricStress stress_trial = hooke_stress(active.lame_lambda, active.shear_modulus, elastic_strain[0],
                                                         elastic_strain[1], elastic_strain[2], elastic_strain[3]);
    const adlite::Scalar mean_stress = (stress_trial.rr + stress_trial.zz + stress_trial.hoop) / 3.0;
    const std::array<adlite::Scalar, component_count> deviatoric_trial = {
        stress_trial.rr - mean_stress,
        stress_trial.zz - mean_stress,
        stress_trial.hoop - mean_stress,
        stress_trial.rz,
    };
    const bool zero_deviatoric_stress = deviatoric_trial[0].value() == 0.0 && deviatoric_trial[1].value() == 0.0 &&
                                        deviatoric_trial[2].value() == 0.0 && deviatoric_trial[3].value() == 0.0;
    constexpr double square_root_two = 1.414213562373095048801688724209698079;
    constexpr double square_root_three_halves = 1.224744871391589049098642037352945695;
    const adlite::Scalar equivalent_trial_stress =
        zero_deviatoric_stress
            ? adlite::Scalar(0.0)
            : square_root_three_halves *
                  adlite::hypot(adlite::hypot(deviatoric_trial[0], deviatoric_trial[1]),
                                adlite::hypot(deviatoric_trial[2], square_root_two * deviatoric_trial[3]));
    if (!std::isfinite(equivalent_trial_stress.value()))
        throw std::overflow_error("Inelastic material trial equivalent stress is not finite");

    const auto finalize = [&](const AxisymmetricStress& stress) {
        for (std::size_t component = 0; component < component_count; ++component) {
            trial_state.elastic_strain[component] =
                elastic_strain[component] -
                (trial_state.plastic_strain[component] - committed.plastic_strain[component]) -
                (trial_state.creep_strain[component] - committed.creep_strain[component]);
        }
        return InelasticStressResponse{stress, trial_state};
    };

    const auto zero_stress_creep_response = [&](double stress_scale_derivative,
                                                const NortonCreepProperties& creep_values) {
        const adlite::Scalar stress_scale =
            adlite::compose(stress_scale_derivative, active.shear_modulus,
                            zero_stress_scale_shear_derivative(stress_scale_derivative, time_step, creep_values));
        const AxisymmetricStress stress = returned_stress(mean_stress, deviatoric_trial, stress_scale);
        update_relaxed_creep_strain(trial_state, committed, deviatoric_trial, stress_scale, active.shear_modulus);
        return finalize(stress);
    };

    const std::shared_ptr<const MaterialFunctionSet>& functions = _thermoelastic_material.properties().functions;
    const bool registered_creep_uses_builtin_integrator = !functions || !functions->has_creep() ||
                                                          functions->creep.name == "norton" ||
                                                          functions->creep.name == "linear_temperature_norton";
    const bool registered_plasticity_uses_builtin_integrator =
        !functions || !functions->has_plasticity() || functions->plasticity.name == "linear_isotropic_hardening" ||
        functions->plasticity.name == "linear_temperature_isotropic_hardening";
    bool registered_builtin_reference_temperatures_match = true;
    if (functions && functions->creep.name == "linear_temperature_norton" &&
        functions->plasticity.name == "linear_temperature_isotropic_hardening") {
        registered_builtin_reference_temperatures_match =
            functions->creep.parameters.value("reference_temperature") ==
            functions->plasticity.parameters.value("reference_temperature");
    }
    if (functions && !(registered_creep_uses_builtin_integrator && registered_plasticity_uses_builtin_integrator &&
                       registered_builtin_reference_temperatures_match)) {
        if (!functions->has_creep() && !functions->has_plasticity())
            return finalize(stress_trial);
        if (equivalent_trial_stress.value() == 0.0) {
            if (!functions->has_creep())
                return finalize(stress_trial);
            const FunctionCreepRoot update =
                solve_function_creep(functions->creep, equivalent_trial_stress, active.shear_modulus, temperature,
                                     committed.equivalent_creep_strain, time, time_step, radius, axial_coordinate);
            const adlite::Scalar stress_scale(update.zero_stress_scale);
            const AxisymmetricStress stress = returned_stress(mean_stress, deviatoric_trial, stress_scale);
            update_relaxed_creep_strain(trial_state, committed, deviatoric_trial, stress_scale, active.shear_modulus);
            return finalize(stress);
        }

        FunctionCoupledUpdate update{};
        if (functions->has_plasticity()) {
            update = solve_function_coupled(*functions, equivalent_trial_stress, active.shear_modulus, temperature,
                                            committed.equivalent_plastic_strain, committed.equivalent_creep_strain,
                                            time, time_step, radius, axial_coordinate);
        } else {
            const FunctionCreepRoot creep =
                solve_function_creep(functions->creep, equivalent_trial_stress, active.shear_modulus, temperature,
                                     committed.equivalent_creep_strain, time, time_step, radius, axial_coordinate);
            update = {creep.equivalent_stress, 0.0, creep.creep_increment};
        }
        const adlite::Scalar stress_scale = update.equivalent_stress / equivalent_trial_stress;
        const AxisymmetricStress stress = returned_stress(mean_stress, deviatoric_trial, stress_scale);
        update_flow_strain(trial_state.plastic_strain, committed.plastic_strain, deviatoric_trial,
                           equivalent_trial_stress, update.plastic_increment);
        update_flow_strain(trial_state.creep_strain, committed.creep_strain, deviatoric_trial, equivalent_trial_stress,
                           update.creep_increment);
        trial_state.equivalent_plastic_strain = committed.equivalent_plastic_strain + update.plastic_increment;
        trial_state.equivalent_creep_strain = committed.equivalent_creep_strain + update.creep_increment;
        return finalize(stress);
    }

    if (_properties.behavior == InelasticBehavior::elastic)
        return finalize(stress_trial);

    if (_properties.behavior == InelasticBehavior::j2_plasticity) {
        const ActiveJ2PlasticityProperties plasticity = active_plasticity_properties(temperature);
        if (equivalent_trial_stress.value() == 0.0)
            return finalize(stress_trial);

        const adlite::Scalar hardening = plasticity.isotropic_hardening_modulus;
        const adlite::Scalar current_yield_stress =
            plasticity.yield_stress + hardening * committed.equivalent_plastic_strain;
        if (!std::isfinite(current_yield_stress.value()))
            throw std::overflow_error("J2 plasticity current yield stress is not finite");

        const double yield_function = equivalent_trial_stress.value() - current_yield_stress.value();
        if (!(yield_function > 0.0))
            return finalize(stress_trial);

        const adlite::Scalar denominator = 3.0 * active.shear_modulus + hardening;
        const adlite::Scalar plastic_increment = (equivalent_trial_stress - current_yield_stress) / denominator;
        const adlite::Scalar returned_equivalent_stress =
            equivalent_trial_stress - 3.0 * active.shear_modulus * plastic_increment;
        const adlite::Scalar stress_scale = returned_equivalent_stress / equivalent_trial_stress;
        const AxisymmetricStress stress = returned_stress(mean_stress, deviatoric_trial, stress_scale);

        update_flow_strain(trial_state.plastic_strain, committed.plastic_strain, deviatoric_trial,
                           equivalent_trial_stress, plastic_increment);
        trial_state.equivalent_plastic_strain = committed.equivalent_plastic_strain + plastic_increment;
        return finalize(stress);
    }

    if (_properties.behavior == InelasticBehavior::norton_creep_j2_plasticity) {
        const ActiveNortonCreepProperties creep = active_creep_properties(temperature);
        const ActiveJ2PlasticityProperties plasticity = active_plasticity_properties(temperature);
        const NortonCreepProperties creep_values = {creep.coefficient.value(), creep.reference_stress.value(),
                                                    creep.stress_exponent.value()};
        const J2PlasticityProperties plasticity_values = {plasticity.yield_stress.value(),
                                                          plasticity.isotropic_hardening_modulus.value()};
        const CoupledUpdate update =
            solve_coupled_update(equivalent_trial_stress.value(), active.shear_modulus.value(), time_step, creep_values,
                                 plasticity_values, committed.equivalent_plastic_strain);
        if (equivalent_trial_stress.value() == 0.0)
            return zero_stress_creep_response(update.stress_derivative, creep_values);

        const CoupledPartials partials =
            coupled_partials(equivalent_trial_stress.value(), active.shear_modulus.value(), creep_values,
                             plasticity_values, committed.equivalent_plastic_strain, update);
        const std::array<adlite::Scalar, 7> coupled_inputs = {equivalent_trial_stress,
                                                              active.shear_modulus,
                                                              creep.coefficient,
                                                              creep.reference_stress,
                                                              creep.stress_exponent,
                                                              plasticity.yield_stress,
                                                              plasticity.isotropic_hardening_modulus};
        const adlite::Scalar returned_equivalent_stress = compose_with_derivatives(
            update.equivalent_stress, coupled_inputs.data(), coupled_inputs.size(), partials.stress.data());
        const adlite::Scalar plastic_increment = compose_with_derivatives(
            update.plastic_increment, coupled_inputs.data(), coupled_inputs.size(), partials.plastic.data());
        const adlite::Scalar creep_increment = compose_with_derivatives(update.creep_increment, coupled_inputs.data(),
                                                                        coupled_inputs.size(), partials.creep.data());
        const adlite::Scalar stress_scale = returned_equivalent_stress / equivalent_trial_stress;
        const AxisymmetricStress stress = returned_stress(mean_stress, deviatoric_trial, stress_scale);

        update_flow_strain(trial_state.plastic_strain, committed.plastic_strain, deviatoric_trial,
                           equivalent_trial_stress, plastic_increment);
        update_flow_strain(trial_state.creep_strain, committed.creep_strain, deviatoric_trial, equivalent_trial_stress,
                           creep_increment);
        trial_state.equivalent_plastic_strain = committed.equivalent_plastic_strain + plastic_increment;
        trial_state.equivalent_creep_strain = committed.equivalent_creep_strain + creep_increment;
        return finalize(stress);
    }

    const ActiveNortonCreepProperties creep = active_creep_properties(temperature);
    const NortonCreepProperties creep_values = {creep.coefficient.value(), creep.reference_stress.value(),
                                                creep.stress_exponent.value()};
    if (time_step == 0.0 || creep.coefficient.value() == 0.0)
        return finalize(stress_trial);

    const NortonRoot root = solve_norton_equivalent_stress(equivalent_trial_stress.value(),
                                                           active.shear_modulus.value(), time_step, creep_values);
    if (equivalent_trial_stress.value() == 0.0)
        return zero_stress_creep_response(root.trial_stress_derivative, creep_values);

    const double evaluated_creep = evaluate_creep_increment(root.equivalent_stress, time_step, creep_values);
    const NortonPartials partials = norton_partials(equivalent_trial_stress.value(), active.shear_modulus.value(),
                                                    time_step, creep_values, root, evaluated_creep);
    const std::array<adlite::Scalar, 5> norton_inputs = {equivalent_trial_stress, active.shear_modulus,
                                                         creep.coefficient, creep.reference_stress,
                                                         creep.stress_exponent};
    const adlite::Scalar returned_equivalent_stress = compose_with_derivatives(
        root.equivalent_stress, norton_inputs.data(), norton_inputs.size(), partials.stress.data());
    const adlite::Scalar creep_increment =
        compose_with_derivatives(evaluated_creep, norton_inputs.data(), norton_inputs.size(), partials.creep.data());
    const adlite::Scalar stress_scale = returned_equivalent_stress / equivalent_trial_stress;
    const AxisymmetricStress stress = returned_stress(mean_stress, deviatoric_trial, stress_scale);

    update_flow_strain(trial_state.creep_strain, committed.creep_strain, deviatoric_trial, equivalent_trial_stress,
                       creep_increment);
    trial_state.equivalent_creep_strain = committed.equivalent_creep_strain + creep_increment;
    return finalize(stress);
}

InelasticStressResponse IsotropicInelasticMaterial::incremental_response(
    const adlite::Scalar& strain_increment_rr, const adlite::Scalar& strain_increment_zz,
    const adlite::Scalar& strain_increment_hoop, const adlite::Scalar& strain_increment_rz,
    const AxisymmetricRotation& rotation, const adlite::Scalar& temperature, double committed_temperature,
    double time_step, const MaterialPointState& committed, double time, double radius, double axial_coordinate) const {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental material committed temperature must be finite and "
                                "positive");

    const adlite::Scalar old_temperature(committed_temperature);
    const AxisymmetricStrain old_imposed =
        _thermoelastic_material.eigenstrain(old_temperature, time - time_step, radius, axial_coordinate);
    const std::array<adlite::Scalar, component_count> synthetic_total = {
        committed.elastic_strain[0] + strain_increment_rr + old_imposed.rr + committed.plastic_strain[0] +
            committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment_zz + old_imposed.zz + committed.plastic_strain[1] +
            committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment_hoop + old_imposed.hoop + committed.plastic_strain[2] +
            committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment_rz + old_imposed.rz + committed.plastic_strain[3] +
            committed.creep_strain[3],
    };

    InelasticStressResponse result =
        response(synthetic_total[0], synthetic_total[1], synthetic_total[2], synthetic_total[3], temperature, time_step,
                 committed, time, radius, axial_coordinate);
    result.stress = rotate_axisymmetric_tensor(result.stress, rotation);
    result.trial_state.elastic_strain =
        components(rotate_axisymmetric_tensor(tensor(result.trial_state.elastic_strain), rotation));
    result.trial_state.plastic_strain =
        components(rotate_axisymmetric_tensor(tensor(result.trial_state.plastic_strain), rotation));
    result.trial_state.creep_strain =
        components(rotate_axisymmetric_tensor(tensor(result.trial_state.creep_strain), rotation));
    return result;
}

MaterialPointState IsotropicInelasticMaterial::state_values(const MaterialPointTrialState& trial_state) {
    MaterialPointState state;
    for (std::size_t component = 0; component < component_count; ++component) {
        state.elastic_strain[component] = trial_state.elastic_strain[component].value();
        state.plastic_strain[component] = trial_state.plastic_strain[component].value();
        state.creep_strain[component] = trial_state.creep_strain[component].value();
    }
    state.equivalent_plastic_strain = trial_state.equivalent_plastic_strain.value();
    state.equivalent_creep_strain = trial_state.equivalent_creep_strain.value();
    validate_material_point_state(state);
    return state;
}

} // namespace fuelsim
