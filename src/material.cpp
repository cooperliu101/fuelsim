#include "fuelsim/material.hpp"

#include <cmath>
#include <stdexcept>

namespace fuelsim {

AxisymmetricStress rotate_axisymmetric_tensor(
    const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation) {
    return {
        rotation.rr * rotation.rr * tensor.rr +
            rotation.rz * rotation.rz * tensor.zz +
            2.0 * rotation.rr * rotation.rz * tensor.rz,
        rotation.zr * rotation.zr * tensor.rr +
            rotation.zz * rotation.zz * tensor.zz +
            2.0 * rotation.zr * rotation.zz * tensor.rz,
        rotation.hoop * rotation.hoop * tensor.hoop,
        rotation.rr * rotation.zr * tensor.rr +
            rotation.rz * rotation.zz * tensor.zz +
            (rotation.rr * rotation.zz + rotation.rz * rotation.zr) *
                tensor.rz,
    };
}

IsotropicThermoelasticMaterial::IsotropicThermoelasticMaterial(
    ThermoelasticProperties properties)
    : _properties(properties), _lame_lambda(0.0), _shear_modulus(0.0),
      _temperature_dependent(false) {
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
        !std::isfinite(_properties.young_modulus_temperature_coefficient) ||
        !std::isfinite(_properties.poisson_ratio_temperature_coefficient) ||
        !std::isfinite(
            _properties.thermal_expansion_temperature_coefficient) ||
        !std::isfinite(_properties.reference_temperature))
        throw std::invalid_argument(
            "Thermoelastic material thermal properties must be finite");
    _lame_lambda =
        _properties.young_modulus * _properties.poisson_ratio /
        ((1.0 + _properties.poisson_ratio) *
         (1.0 - 2.0 * _properties.poisson_ratio));
    _shear_modulus =
        _properties.young_modulus /
        (2.0 * (1.0 + _properties.poisson_ratio));
    _temperature_dependent =
        _properties.young_modulus_temperature_coefficient != 0.0 ||
        _properties.poisson_ratio_temperature_coefficient != 0.0 ||
        _properties.thermal_expansion_temperature_coefficient != 0.0;
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

ActiveThermoelasticProperties
IsotropicThermoelasticMaterial::active_properties(
    const adlite::Scalar& temperature) const {
    if (!std::isfinite(temperature.value()))
        throw std::domain_error(
            "Thermoelastic material temperature must be finite");
    if (!_temperature_dependent)
        return {_properties.young_modulus, _properties.poisson_ratio,
                _properties.thermal_expansion, _lame_lambda,
                _shear_modulus};
    const adlite::Scalar temperature_change =
        temperature - _properties.reference_temperature;
    const adlite::Scalar young_modulus =
        _properties.young_modulus +
        _properties.young_modulus_temperature_coefficient *
            temperature_change;
    const adlite::Scalar poisson_ratio =
        _properties.poisson_ratio +
        _properties.poisson_ratio_temperature_coefficient *
            temperature_change;
    const adlite::Scalar thermal_expansion =
        _properties.thermal_expansion +
        _properties.thermal_expansion_temperature_coefficient *
            temperature_change;
    if (!std::isfinite(young_modulus.value()) ||
        !(young_modulus.value() > 0.0))
        throw std::domain_error(
            "Thermoelastic material active young_modulus must be positive");
    if (!std::isfinite(poisson_ratio.value()) ||
        !(poisson_ratio.value() > -1.0 && poisson_ratio.value() < 0.5))
        throw std::domain_error(
            "Thermoelastic material active poisson_ratio must lie between "
            "-1 and 0.5");
    if (!std::isfinite(thermal_expansion.value()))
        throw std::domain_error(
            "Thermoelastic material active thermal_expansion must be "
            "finite");
    const adlite::Scalar lame_lambda =
        young_modulus * poisson_ratio /
        ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const adlite::Scalar shear_modulus =
        young_modulus / (2.0 * (1.0 + poisson_ratio));
    return {young_modulus, poisson_ratio, thermal_expansion, lame_lambda,
            shear_modulus};
}

AxisymmetricStress IsotropicThermoelasticMaterial::stress(
    const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz,
    const adlite::Scalar& temperature) const {
    if (!_temperature_dependent) {
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
    const ActiveThermoelasticProperties active = active_properties(temperature);
    const adlite::Scalar thermal_strain =
        active.thermal_expansion *
        (temperature - _properties.reference_temperature);

    const adlite::Scalar elastic_rr = strain_rr - thermal_strain;
    const adlite::Scalar elastic_zz = strain_zz - thermal_strain;
    const adlite::Scalar elastic_hoop = strain_hoop - thermal_strain;
    const adlite::Scalar trace = elastic_rr + elastic_zz + elastic_hoop;

    return {
        active.lame_lambda * trace + 2.0 * active.shear_modulus * elastic_rr,
        active.lame_lambda * trace + 2.0 * active.shear_modulus * elastic_zz,
        active.lame_lambda * trace +
            2.0 * active.shear_modulus * elastic_hoop,
        2.0 * active.shear_modulus * strain_rz,
    };
}

} // namespace fuelsim
