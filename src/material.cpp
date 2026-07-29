#include "fuelsim/material.hpp"

#include <cmath>
#include <stdexcept>

namespace fuelsim {

IsotropicThermoelasticMaterial::IsotropicThermoelasticMaterial(
    ThermoelasticProperties properties)
    : properties_(properties), lame_lambda_(0.0), shear_modulus_(0.0) {
    if (!std::isfinite(properties_.conductivity_inverse_temperature) ||
        !(properties_.conductivity_inverse_temperature >= 0.0) ||
        !std::isfinite(properties_.conductivity_offset) ||
        !(properties_.conductivity_offset >= 0.0) ||
        !(properties_.conductivity_inverse_temperature > 0.0 ||
          properties_.conductivity_offset > 0.0))
        throw std::invalid_argument(
            "Thermoelastic material conductivity coefficients must be "
            "finite, nonnegative, and not both zero");
    if (!std::isfinite(properties_.young_modulus) ||
        !(properties_.young_modulus > 0.0))
        throw std::invalid_argument(
            "Thermoelastic material young_modulus must be positive");
    if (!std::isfinite(properties_.poisson_ratio) ||
        !(properties_.poisson_ratio > -1.0 && properties_.poisson_ratio < 0.5))
        throw std::invalid_argument(
            "Thermoelastic material poisson_ratio must lie between -1 and "
            "0.5");
    if (!std::isfinite(properties_.thermal_expansion) ||
        !std::isfinite(properties_.reference_temperature))
        throw std::invalid_argument(
            "Thermoelastic material thermal properties must be finite");

    const double nu = properties_.poisson_ratio;
    const double E = properties_.young_modulus;
    lame_lambda_ = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    shear_modulus_ = E / (2.0 * (1.0 + nu));
}

const ThermoelasticProperties&
IsotropicThermoelasticMaterial::properties() const noexcept {
    return properties_;
}

adlite::Scalar IsotropicThermoelasticMaterial::conductivity(
    const adlite::Scalar& temperature) const {
    if (!std::isfinite(temperature.value()) || !(temperature.value() > 0.0))
        throw std::domain_error(
            "Thermoelastic material temperature must be finite and "
            "positive");
    return properties_.conductivity_inverse_temperature / temperature +
           properties_.conductivity_offset;
}

AxisymmetricStress IsotropicThermoelasticMaterial::stress(
    const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz,
    const adlite::Scalar& temperature) const {
    const adlite::Scalar thermal_strain =
        properties_.thermal_expansion *
        (temperature - properties_.reference_temperature);

    const adlite::Scalar elastic_rr = strain_rr - thermal_strain;
    const adlite::Scalar elastic_zz = strain_zz - thermal_strain;
    const adlite::Scalar elastic_hoop = strain_hoop - thermal_strain;
    const adlite::Scalar trace = elastic_rr + elastic_zz + elastic_hoop;

    return {
        lame_lambda_ * trace + 2.0 * shear_modulus_ * elastic_rr,
        lame_lambda_ * trace + 2.0 * shear_modulus_ * elastic_zz,
        lame_lambda_ * trace + 2.0 * shear_modulus_ * elastic_hoop,
        2.0 * shear_modulus_ * strain_rz,
    };
}

} // namespace fuelsim
