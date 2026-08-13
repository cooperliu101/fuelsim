#pragma once
#include "fuelsim/material_functions.hpp"
#include <adlite/adlite.hpp>
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
struct AxisymmetricRotation final {
    adlite::Scalar rr{1.0}, rz{0.0}, zr{0.0}, zz{1.0}, hoop{1.0};
};
AxisymmetricStress rotate_axisymmetric_tensor(const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation);
class IsotropicThermoelasticMaterial final {
  public:
    explicit IsotropicThermoelasticMaterial(ThermoelasticProperties properties);
    const MaterialFunctionSet& functions() const noexcept { return *_properties.functions; }
    adlite::Scalar conductivity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    adlite::Scalar heat_capacity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    ActiveThermoelasticProperties active_properties(
        const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    AxisymmetricStrain eigenstrain_rz(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    SymmetricTensor3 eigenstrain(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    AxisymmetricStress stress(const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz, const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    SymmetricTensor3 stress(
        const SymmetricTensor3& strain, const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;

  private:
    ThermoelasticProperties _properties;
};
} // namespace fuelsim
