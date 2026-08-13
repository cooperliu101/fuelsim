#ifndef FUELSIM_MATERIAL_HPP
#define FUELSIM_MATERIAL_HPP

#include "fuelsim/material_functions.hpp"

#include <adlite/adlite.hpp>

#include <memory>

namespace fuelsim {

struct ThermoelasticProperties final {
    double conductivity_inverse_temperature;
    double conductivity_offset;
    double young_modulus;
    double poisson_ratio;
    double thermal_expansion;
    double reference_temperature;
    double young_modulus_temperature_coefficient = 0.0;
    double poisson_ratio_temperature_coefficient = 0.0;
    double thermal_expansion_temperature_coefficient = 0.0;
    std::shared_ptr<const MaterialFunctionSet> functions{};
};

struct ActiveThermoelasticProperties final {
    adlite::Scalar young_modulus;
    adlite::Scalar poisson_ratio;
    adlite::Scalar thermal_expansion;
    adlite::Scalar lame_lambda;
    adlite::Scalar shear_modulus;
};

struct AxisymmetricStress final {
    adlite::Scalar rr;
    adlite::Scalar zz;
    adlite::Scalar hoop;
    adlite::Scalar rz;
};

struct AxisymmetricStressValues final {
    double rr;
    double zz;
    double hoop;
    double rz;
};

struct SymmetricTensor3Values final {
    double xx;
    double yy;
    double zz;
    double xy;
    double yz;
    double xz;
};

struct AxisymmetricRotation final {
    adlite::Scalar rr{1.0};
    adlite::Scalar rz{0.0};
    adlite::Scalar zr{0.0};
    adlite::Scalar zz{1.0};
    adlite::Scalar hoop{1.0};
};

AxisymmetricStress rotate_axisymmetric_tensor(const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation);

class IsotropicThermoelasticMaterial final {
  public:
    explicit IsotropicThermoelasticMaterial(ThermoelasticProperties properties);

    const ThermoelasticProperties& properties() const noexcept;

    adlite::Scalar conductivity(const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0,
                                double axial_coordinate = 0.0) const;

    adlite::Scalar heat_capacity(const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0,
                                 double axial_coordinate = 0.0) const;

    ActiveThermoelasticProperties active_properties(const adlite::Scalar& temperature, double time = 0.0,
                                                    double radius = 0.0, double axial_coordinate = 0.0) const;

    AxisymmetricStrain eigenstrain(const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0,
                                   double axial_coordinate = 0.0) const;

    AxisymmetricStress stress(const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
                              const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz,
                              const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0,
                              double axial_coordinate = 0.0) const;

    adlite::Scalar conductivity_cartesian(const adlite::Scalar& temperature, double time, double x, double y,
                                          double z) const;
    adlite::Scalar heat_capacity_cartesian(const adlite::Scalar& temperature, double time, double x, double y,
                                           double z) const;
    ActiveThermoelasticProperties active_properties_cartesian(const adlite::Scalar& temperature, double time, double x,
                                                              double y, double z) const;
    SymmetricTensor3 eigenstrain_cartesian(const adlite::Scalar& temperature, double time, double x, double y,
                                           double z) const;
    SymmetricTensor3 stress_cartesian(const SymmetricTensor3& strain, const adlite::Scalar& temperature, double time,
                                      double x, double y, double z) const;

  private:
    ThermoelasticProperties _properties;
    double _lame_lambda;
    double _shear_modulus;
    bool _temperature_dependent;
};

} // namespace fuelsim

#endif
