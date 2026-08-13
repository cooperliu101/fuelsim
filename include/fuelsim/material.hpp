#pragma once
#include "fuelsim/material_functions.hpp"
#include <adlite/adlite.hpp>
#include <memory>
namespace fuelsim {
struct ThermoelasticProperties final {
    double conductivity_inverse_temperature, conductivity_offset, young_modulus, poisson_ratio, thermal_expansion,
        reference_temperature;
    double young_modulus_temperature_coefficient = 0.0, poisson_ratio_temperature_coefficient = 0.0,
           thermal_expansion_temperature_coefficient = 0.0;
    std::shared_ptr<const MaterialFunctionSet> functions{};
};
struct ActiveThermoelasticProperties final {
    adlite::Scalar young_modulus, poisson_ratio, thermal_expansion, lame_lambda, shear_modulus;
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
    const ThermoelasticProperties& properties() const noexcept { return _properties; }
    adlite::Scalar conductivity(
        const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    adlite::Scalar heat_capacity(
        const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    ActiveThermoelasticProperties active_properties(
        const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    AxisymmetricStrain eigenstrain(
        const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    AxisymmetricStress stress(const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz, const adlite::Scalar& temperature,
        double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    adlite::Scalar conductivity_cartesian(
        const adlite::Scalar& temperature, double time, double x, double y, double z) const;
    adlite::Scalar heat_capacity_cartesian(
        const adlite::Scalar& temperature, double time, double x, double y, double z) const;
    ActiveThermoelasticProperties active_properties_cartesian(
        const adlite::Scalar& temperature, double time, double x, double y, double z) const;
    SymmetricTensor3 eigenstrain_cartesian(
        const adlite::Scalar& temperature, double time, double x, double y, double z) const;
    SymmetricTensor3 stress_cartesian(const SymmetricTensor3& strain, const adlite::Scalar& temperature, double time,
        double x, double y, double z) const;

  private:
    ThermoelasticProperties _properties;
    double _lame_lambda, _shear_modulus;
    bool _temperature_dependent;
};
} // namespace fuelsim
