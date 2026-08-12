#ifndef FUELSIM_MATERIAL_FUNCTIONS_HPP
#define FUELSIM_MATERIAL_FUNCTIONS_HPP

#include <adlite/adlite.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim {

struct MaterialParameterDefinition final {
    std::string name;
    std::string unit;
};

struct MaterialParameterValue final {
    std::string name;
    double value;
};

class MaterialParameters final {
  public:
    MaterialParameters() = default;
    explicit MaterialParameters(std::vector<MaterialParameterValue> values);

    std::size_t size() const noexcept;
    double value(const std::string& name) const;
    double value(std::size_t index) const;
    const std::vector<MaterialParameterValue>& values() const noexcept;

  private:
    std::vector<MaterialParameterValue> _values;
};

struct ThermalPropertyInput final {
    adlite::Scalar temperature;
    double time;
    double x;
    double y;
    double z;
    const MaterialParameters* parameters;
};

struct ThermalPropertyOutput final {
    adlite::Scalar conductivity;
    adlite::Scalar density;
    adlite::Scalar specific_heat;
};

struct ElasticPropertyInput final {
    adlite::Scalar temperature;
    double time;
    double x;
    double y;
    double z;
    const MaterialParameters* parameters;
};

struct ElasticPropertyOutput final {
    adlite::Scalar young_modulus;
    adlite::Scalar poisson_ratio;
};

struct EigenstrainInput final {
    adlite::Scalar temperature;
    double time;
    double x;
    double y;
    double z;
    const MaterialParameters* parameters;
};

struct AxisymmetricStrain final {
    adlite::Scalar rr;
    adlite::Scalar zz;
    adlite::Scalar hoop;
    adlite::Scalar rz;
};

struct SymmetricTensor3 final {
    adlite::Scalar xx;
    adlite::Scalar yy;
    adlite::Scalar zz;
    adlite::Scalar xy;
    adlite::Scalar yz;
    adlite::Scalar xz;
};

struct CreepRateInput final {
    adlite::Scalar equivalent_stress;
    adlite::Scalar temperature;
    adlite::Scalar equivalent_creep_strain;
    double time;
    double x;
    double y;
    double z;
    const MaterialParameters* parameters;
};

struct PlasticFlowStressInput final {
    adlite::Scalar equivalent_plastic_strain;
    adlite::Scalar temperature;
    double time;
    double x;
    double y;
    double z;
    const MaterialParameters* parameters;
};

using ThermalPropertyFunction = void (*)(const ThermalPropertyInput&, ThermalPropertyOutput&);
using ElasticPropertyFunction = void (*)(const ElasticPropertyInput&, ElasticPropertyOutput&);
using EigenstrainFunction = void (*)(const EigenstrainInput&, SymmetricTensor3&);
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
    std::string instance_name;
    std::string name;
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

    bool has_creep() const noexcept;
    bool has_plasticity() const noexcept;
    std::uint64_t signature() const noexcept;
};

class MaterialFunctionRegistry final {
  public:
    void add_thermal(std::string name, std::vector<MaterialParameterDefinition> parameters, ThermalPropertyFunction function, std::uint32_t version = 1);
    void add_elasticity(std::string name, std::vector<MaterialParameterDefinition> parameters, ElasticPropertyFunction function, std::uint32_t version = 1);
    void add_eigenstrain(std::string name, std::vector<MaterialParameterDefinition> parameters, EigenstrainFunction function, std::uint32_t version = 1);
    void add_creep(std::string name, std::vector<MaterialParameterDefinition> parameters, CreepRateFunction function, std::uint32_t version = 1);
    void add_plasticity(std::string name, std::vector<MaterialParameterDefinition> parameters, PlasticFlowStressFunction function, std::uint32_t version = 1);

    const std::vector<MaterialParameterDefinition>& thermal_parameters(const std::string& name) const;
    const std::vector<MaterialParameterDefinition>& elasticity_parameters(const std::string& name) const;
    const std::vector<MaterialParameterDefinition>& eigenstrain_parameters(const std::string& name) const;
    const std::vector<MaterialParameterDefinition>& creep_parameters(const std::string& name) const;
    const std::vector<MaterialParameterDefinition>& plasticity_parameters(const std::string& name) const;

    ThermalFunctionInstance bind_thermal(const std::string& name, std::vector<MaterialParameterValue> values) const;
    ElasticFunctionInstance bind_elasticity(const std::string& name, std::vector<MaterialParameterValue> values) const;
    EigenstrainFunctionInstance bind_eigenstrain(const std::string& instance_name, const std::string& name, std::vector<MaterialParameterValue> values) const;
    CreepFunctionInstance bind_creep(const std::string& name, std::vector<MaterialParameterValue> values) const;
    PlasticFunctionInstance bind_plasticity(const std::string& name, std::vector<MaterialParameterValue> values) const;

  private:
    struct ThermalRegistration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        ThermalPropertyFunction function;
    };
    struct ElasticRegistration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        ElasticPropertyFunction function;
    };
    struct EigenstrainRegistration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        EigenstrainFunction function;
    };
    struct CreepRegistration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        CreepRateFunction function;
    };
    struct PlasticRegistration final {
        std::string name;
        std::uint32_t version;
        std::vector<MaterialParameterDefinition> parameters;
        PlasticFlowStressFunction function;
    };

    std::vector<ThermalRegistration> _thermal;
    std::vector<ElasticRegistration> _elasticity;
    std::vector<EigenstrainRegistration> _eigenstrain;
    std::vector<CreepRegistration> _creep;
    std::vector<PlasticRegistration> _plasticity;
};

MaterialFunctionRegistry make_builtin_material_function_registry();

} // namespace fuelsim

#endif
