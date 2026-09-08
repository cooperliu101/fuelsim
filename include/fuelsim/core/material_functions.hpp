#pragma once
#include <adlite/adlite.hpp>
#include <cstdint>
#include <functional>
#include <string>
#include <variant>
#include <vector>

namespace fuelsim {
struct MaterialParameterDefinition final {
    std::string name, unit;
};

struct MaterialParameterValue final {
    std::string name;
    double value;
};

struct MaterialFunctionContext final {
    double time = 0.0, x = 0.0, y = 0.0, z = 0.0;
};

class MaterialParameters final {
  public:
    MaterialParameters() = default;
    explicit MaterialParameters(std::vector<MaterialParameterValue> values);
    double value(const std::string& name) const;

    const std::vector<MaterialParameterValue>& values() const noexcept { return _values; }

  private:
    std::vector<MaterialParameterValue> _values;
};

struct ThermoelasticFunctionInput final {
    adlite::Scalar temperature;
    MaterialFunctionContext context;
};

struct ThermalPropertyOutput final {
    adlite::Scalar conductivity, density, specific_heat;
};

struct ElasticPropertyOutput final {
    adlite::Scalar young_modulus, poisson_ratio;
};

struct AxisymmetricStrain final {
    adlite::Scalar rr, zz, hoop, rz;
};

struct SymmetricTensor3 final {
    adlite::Scalar xx, yy, zz, xy, yz, xz;
};

struct CreepRateInput final {
    adlite::Scalar equivalent_stress, temperature, equivalent_creep_strain;
    MaterialFunctionContext context;
};

struct PlasticFlowStressInput final {
    adlite::Scalar equivalent_plastic_strain, temperature;
    MaterialFunctionContext context;
};

using ThermalPropertyEvaluator = std::function<void(const ThermoelasticFunctionInput&, ThermalPropertyOutput&)>;
using ThermalPropertyBinder = ThermalPropertyEvaluator (*)(const MaterialParameters&);
using ElasticPropertyEvaluator = std::function<void(const ThermoelasticFunctionInput&, ElasticPropertyOutput&)>;
using ElasticPropertyBinder = ElasticPropertyEvaluator (*)(const MaterialParameters&);
using EigenstrainEvaluator = std::function<void(const ThermoelasticFunctionInput&, SymmetricTensor3&)>;
using EigenstrainBinder = EigenstrainEvaluator (*)(const MaterialParameters&);
using CreepRateEvaluator = std::function<adlite::Scalar(const CreepRateInput&)>;
using CreepRateBinder = CreepRateEvaluator (*)(const MaterialParameters&);
using PlasticFlowStressEvaluator = std::function<adlite::Scalar(const PlasticFlowStressInput&)>;
using PlasticFlowStressBinder = PlasticFlowStressEvaluator (*)(const MaterialParameters&);

struct CreepBuiltinParameters final {
    enum class Kind { custom, norton, linear_temperature_norton };
    Kind kind = Kind::custom;
    double coefficient = 0.0, reference_stress = 0.0, stress_exponent = 0.0;
    double reference_temperature = 0.0;
    double coefficient_temperature_coefficient = 0.0;
    double reference_stress_temperature_coefficient = 0.0;
    double stress_exponent_temperature_coefficient = 0.0;
};

struct PlasticBuiltinParameters final {
    enum class Kind { custom, linear_isotropic_hardening, linear_temperature_isotropic_hardening };
    Kind kind = Kind::custom;
    double yield_stress = 0.0, hardening_modulus = 0.0, reference_temperature = 0.0;
    double yield_stress_temperature_coefficient = 0.0;
    double hardening_temperature_coefficient = 0.0;
};

struct ThermalFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    ThermalPropertyEvaluator function;
};

struct ElasticFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    ElasticPropertyEvaluator function;
};

struct EigenstrainFunctionInstance final {
    std::string instance_name, name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    EigenstrainEvaluator function;
};

struct CreepFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    CreepRateEvaluator function;
    CreepBuiltinParameters builtin;
};

struct PlasticFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    PlasticFlowStressEvaluator function;
    PlasticBuiltinParameters builtin;
};

struct MaterialFunctionSet final {
    std::string name;
    ThermalFunctionInstance thermal;
    ElasticFunctionInstance elasticity;
    std::vector<EigenstrainFunctionInstance> eigenstrains;
    CreepFunctionInstance creep;
    PlasticFunctionInstance plasticity;

    bool has_creep() const noexcept { return creep.function != nullptr; }

    bool has_plasticity() const noexcept { return plasticity.function != nullptr; }

    std::uint64_t signature() const noexcept;
};

class MaterialFunctionRegistry final {
  public:
    void add_thermal(std::string name,
        std::vector<MaterialParameterDefinition> parameters,
        ThermalPropertyBinder function,
        std::uint32_t version = 1);
    void add_elasticity(std::string name,
        std::vector<MaterialParameterDefinition> parameters,
        ElasticPropertyBinder function,
        std::uint32_t version = 1);
    void add_eigenstrain(std::string name,
        std::vector<MaterialParameterDefinition> parameters,
        EigenstrainBinder function,
        std::uint32_t version = 1);
    void add_creep(std::string name,
        std::vector<MaterialParameterDefinition> parameters,
        CreepRateBinder function,
        std::uint32_t version = 1);
    void add_plasticity(std::string name,
        std::vector<MaterialParameterDefinition> parameters,
        PlasticFlowStressBinder function,
        std::uint32_t version = 1);
    ThermalFunctionInstance bind_thermal(const std::string& name, std::vector<MaterialParameterValue> values) const;
    ElasticFunctionInstance bind_elasticity(const std::string& name, std::vector<MaterialParameterValue> values) const;
    EigenstrainFunctionInstance bind_eigenstrain(const std::string& instance_name,
        const std::string& name,
        std::vector<MaterialParameterValue> values) const;
    CreepFunctionInstance bind_creep(const std::string& name, std::vector<MaterialParameterValue> values) const;
    PlasticFunctionInstance bind_plasticity(const std::string& name, std::vector<MaterialParameterValue> values) const;

  private:
    enum class Category { thermal, elasticity, eigenstrain, creep, plasticity };
    using Function = std::variant<ThermalPropertyBinder,
        ElasticPropertyBinder,
        EigenstrainBinder,
        CreepRateBinder,
        PlasticFlowStressBinder>;

    struct Registration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        Category category;
        Function function;
    };

    void add_registration(std::string name,
        std::vector<MaterialParameterDefinition> parameters,
        std::uint32_t version,
        Category category,
        Function function,
        bool has_function,
        const char* category_name);
    const Registration& find_registration(Category category, const std::string& name, const char* category_name) const;
    std::vector<Registration> _registrations;
};

MaterialFunctionRegistry make_builtin_material_function_registry();
} // namespace fuelsim
