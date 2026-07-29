#include "fuelsim/material.hpp"

#include <cmath>
#include <stdexcept>

namespace fuelsim {

IsotropicThermoelasticMaterial::IsotropicThermoelasticMaterial(
    ThermoelasticProperties properties)
    : _properties(properties), _lame_lambda(0.0), _shear_modulus(0.0) {
    if (!std::isfinite(_properties.conductivity_inverse_temperature) ||
        !(_properties.conductivity_inverse_temperature >= 0.0) ||
        !std::isfinite(_properties.conductivity_offset) ||
        !(_properties.conductivity_offset >= 0.0) ||
        !(_properties.conductivity_inverse_temperature > 0.0 ||
          _properties.conductivity_offset > 0.0))
        throw std::invalid_argument(
            "Thermoelastic material conductivity coefficients must be "
            "finite, nonnegative, and not both zero");
    if (!std::isfinite(_properties.young_modulus) ||
        !(_properties.young_modulus > 0.0))
        throw std::invalid_argument(
            "Thermoelastic material young_modulus must be positive");
    if (!std::isfinite(_properties.poisson_ratio) ||
        !(_properties.poisson_ratio > -1.0 && _properties.poisson_ratio < 0.5))
        throw std::invalid_argument(
            "Thermoelastic material poisson_ratio must lie between -1 and "
            "0.5");
    if (!std::isfinite(_properties.thermal_expansion) ||
        !std::isfinite(_properties.reference_temperature))
        throw std::invalid_argument(
            "Thermoelastic material thermal properties must be finite");

    const double nu = _properties.poisson_ratio;
    const double E = _properties.young_modulus;
    _lame_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    _shear_modulus = E / (2.0 * (1.0 + nu));
}

const ThermoelasticProperties&
IsotropicThermoelasticMaterial::properties() const noexcept {
    return _properties;
}

adlite::Scalar IsotropicThermoelasticMaterial::conductivity(
    const adlite::Scalar& temperature) const {
    if (!std::isfinite(temperature.value()) || !(temperature.value() > 0.0))
        throw std::domain_error(
            "Thermoelastic material temperature must be finite and "
            "positive");
    return _properties.conductivity_inverse_temperature / temperature +
           _properties.conductivity_offset;
}

AxisymmetricStress IsotropicThermoelasticMaterial::stress(
    const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz,
    const adlite::Scalar& temperature) const {
    const adlite::Scalar thermal_strain =
        _properties.thermal_expansion *
        (temperature - _properties.reference_temperature);

    const adlite::Scalar elastic_rr = strain_rr - thermal_strain;
    const adlite::Scalar elastic_zz = strain_zz - thermal_strain;
    const adlite::Scalar elastic_hoop = strain_hoop - thermal_strain;
    const adlite::Scalar trace = elastic_rr + elastic_zz + elastic_hoop;

    return {
        _lame_lambda * trace + 2.0 * _shear_modulus * elastic_rr,
        _lame_lambda * trace + 2.0 * _shear_modulus * elastic_zz,
        _lame_lambda * trace + 2.0 * _shear_modulus * elastic_hoop,
        2.0 * _shear_modulus * strain_rz,
    };
}

} // namespace fuelsim
