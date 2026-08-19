#pragma once
#include <adlite/adlite.hpp>
#include <cstdint>
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
    const MaterialParameters* parameters;
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
    const MaterialParameters* parameters;
};

struct PlasticFlowStressInput final {
    adlite::Scalar equivalent_plastic_strain, temperature;
    MaterialFunctionContext context;
    const MaterialParameters* parameters;
};

using ThermalPropertyFunction = void (*)(const ThermoelasticFunctionInput&, ThermalPropertyOutput&);
using ElasticPropertyFunction = void (*)(const ThermoelasticFunctionInput&, ElasticPropertyOutput&);
using EigenstrainFunction = void (*)(const ThermoelasticFunctionInput&, SymmetricTensor3&);
using CreepRateFunction = adlite::Scalar (*)(const CreepRateInput&);
using PlasticFlowStressFunction = adlite::Scalar (*)(const PlasticFlowStressInput&);

struct ThermalFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    ThermalPropertyFunction function = nullptr;
};

struct ElasticFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    ElasticPropertyFunction function = nullptr;
};

struct EigenstrainFunctionInstance final {
    std::string instance_name, name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    EigenstrainFunction function = nullptr;
};

struct CreepFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    CreepRateFunction function = nullptr;
};

struct PlasticFunctionInstance final {
    std::string name;
    std::uint32_t version = 0;
    MaterialParameters parameters;
    PlasticFlowStressFunction function = nullptr;
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
    void add_thermal(std::string name, std::vector<MaterialParameterDefinition> parameters,
        ThermalPropertyFunction function, std::uint32_t version = 1);
    void add_elasticity(std::string name, std::vector<MaterialParameterDefinition> parameters,
        ElasticPropertyFunction function, std::uint32_t version = 1);
    void add_eigenstrain(std::string name, std::vector<MaterialParameterDefinition> parameters,
        EigenstrainFunction function, std::uint32_t version = 1);
    void add_creep(std::string name, std::vector<MaterialParameterDefinition> parameters, CreepRateFunction function,
        std::uint32_t version = 1);
    void add_plasticity(std::string name, std::vector<MaterialParameterDefinition> parameters,
        PlasticFlowStressFunction function, std::uint32_t version = 1);
    ThermalFunctionInstance bind_thermal(const std::string& name, std::vector<MaterialParameterValue> values) const;
    ElasticFunctionInstance bind_elasticity(const std::string& name, std::vector<MaterialParameterValue> values) const;
    EigenstrainFunctionInstance bind_eigenstrain(
        const std::string& instance_name, const std::string& name, std::vector<MaterialParameterValue> values) const;
    CreepFunctionInstance bind_creep(const std::string& name, std::vector<MaterialParameterValue> values) const;
    PlasticFunctionInstance bind_plasticity(const std::string& name, std::vector<MaterialParameterValue> values) const;

  private:
    enum class Category { thermal, elasticity, eigenstrain, creep, plasticity };
    using Function = std::variant<ThermalPropertyFunction, ElasticPropertyFunction, EigenstrainFunction,
        CreepRateFunction, PlasticFlowStressFunction>;

    struct Registration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        Category category;
        Function function;
    };

    void add_registration(std::string name, std::vector<MaterialParameterDefinition> parameters, std::uint32_t version,
        Category category, Function function, bool has_function, const char* category_name);
    const Registration& find_registration(Category category, const std::string& name, const char* category_name) const;
    std::vector<Registration> _registrations;
};

MaterialFunctionRegistry make_builtin_material_function_registry();
} // namespace fuelsim
