#ifndef FUELSIM_MATERIAL_HPP
#define FUELSIM_MATERIAL_HPP

#include <adlite/adlite.hpp>

namespace fuelsim {

struct ThermoelasticProperties final {
    double conductivity_inverse_temperature;
    double conductivity_offset;
    double young_modulus;
    double poisson_ratio;
    double thermal_expansion;
    double reference_temperature;
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

class IsotropicThermoelasticMaterial final {
  public:
    explicit IsotropicThermoelasticMaterial(ThermoelasticProperties properties);

    const ThermoelasticProperties& properties() const noexcept;

    adlite::Scalar conductivity(const adlite::Scalar& temperature) const;

    AxisymmetricStress stress(const adlite::Scalar& strain_rr,
                              const adlite::Scalar& strain_zz,
                              const adlite::Scalar& strain_hoop,
                              const adlite::Scalar& strain_rz,
                              const adlite::Scalar& temperature) const;

  private:
    ThermoelasticProperties properties_;
    double lame_lambda_;
    double shear_modulus_;
};

} // namespace fuelsim

#endif
