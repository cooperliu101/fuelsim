#pragma once
#include "material_functions.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <memory>

namespace fuelsim {
struct ThermoelasticProperties final {
    std::shared_ptr<const MaterialFunctionSet> functions;
    double reference_young_modulus;
};

struct ActiveThermoelasticProperties final {
    adlite::Scalar lame_lambda, shear_modulus;
};

struct AxisymmetricStress final {
    adlite::Scalar rr, zz, hoop, rz;
};

struct AxisymmetricStressValues final {
    double rr, zz, hoop, rz;
};

struct SymmetricTensor3Values final {
    double xx, yy, zz, xy, yz, xz;
};

struct CartesianMaterialPointState final {
    std::array<double, 6> elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0;
    SymmetricTensor3Values stress{};
};

struct CartesianInelasticStressResponse final {
    SymmetricTensor3 stress;
    CartesianMaterialPointState trial_state;
};

struct CartesianRotation final {
    adlite::Scalar xx{1.0}, xy{0.0}, xz{0.0};
    adlite::Scalar yx{0.0}, yy{1.0}, yz{0.0};
    adlite::Scalar zx{0.0}, zy{0.0}, zz{1.0};
};

struct AxisymmetricRotation final {
    adlite::Scalar rr{1.0}, rz{0.0}, zr{0.0}, zz{1.0}, hoop{1.0};
};

struct MaterialPointState final {
    std::array<double, 4> elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0;
    AxisymmetricStressValues stress{};
};

struct MaterialPointTrialState final {
    std::array<adlite::Scalar, 4> elastic_strain{}, plastic_strain{}, creep_strain{};
    adlite::Scalar equivalent_plastic_strain{0.0}, equivalent_creep_strain{0.0};
};

struct InelasticStressResponse final {
    AxisymmetricStress stress;
    MaterialPointTrialState trial_state;
};

AxisymmetricStress rotate_axisymmetric_tensor(const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation);
SymmetricTensor3 rotate_cartesian_tensor(const SymmetricTensor3& tensor, const CartesianRotation& rotation);

class IsotropicThermoelasticMaterial final {
  public:
    explicit IsotropicThermoelasticMaterial(ThermoelasticProperties properties);

    const MaterialFunctionSet& functions() const noexcept { return *_properties.functions; }

    adlite::Scalar conductivity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    adlite::Scalar heat_capacity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    ActiveThermoelasticProperties active_properties(const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    AxisymmetricStrain eigenstrain_rz(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    SymmetricTensor3 eigenstrain(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    AxisymmetricStress stress(const adlite::Scalar& strain_rr,
        const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop,
        const adlite::Scalar& strain_rz,
        const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    SymmetricTensor3 stress(const SymmetricTensor3& strain,
        const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    SymmetricTensor3Values
    stress_values(const SymmetricTensor3Values& strain, double temperature, MaterialFunctionContext context = {}) const;
    InelasticStressResponse response(const adlite::Scalar& strain_rr,
        const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop,
        const adlite::Scalar& strain_rz,
        const adlite::Scalar& temperature,
        double time_step,
        const MaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianInelasticStressResponse response(const SymmetricTensor3& strain,
        const adlite::Scalar& temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianInelasticStressResponse incremental_response(const SymmetricTensor3& strain_increment,
        const CartesianRotation& rotation,
        const adlite::Scalar& temperature,
        double committed_temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianMaterialPointState response_values(const SymmetricTensor3Values& strain,
        double temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianMaterialPointState incremental_response_values(const SymmetricTensor3Values& strain_increment,
        const CartesianRotation& rotation,
        double temperature,
        double committed_temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    InelasticStressResponse incremental_response(const adlite::Scalar& strain_increment_rr,
        const adlite::Scalar& strain_increment_zz,
        const adlite::Scalar& strain_increment_hoop,
        const adlite::Scalar& strain_increment_rz,
        const AxisymmetricRotation& rotation,
        const adlite::Scalar& temperature,
        double committed_temperature,
        double time_step,
        const MaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    static MaterialPointState state_values(const MaterialPointTrialState& trial_state);

  private:
    ThermoelasticProperties _properties;
};
} // namespace fuelsim

namespace fuelsim {
struct AxisymmetricStressTangent final {
    AxisymmetricStressValues stress;
    std::array<std::array<double, 4>, 4> tangent{};
    std::array<double, 4> thermal{};
};

AxisymmetricStressTangent evaluate_axisymmetric_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 4>& fed_strain,
    double temperature,
    double time_step,
    const MaterialPointState* committed_material,
    MaterialFunctionContext context);
} // namespace fuelsim
