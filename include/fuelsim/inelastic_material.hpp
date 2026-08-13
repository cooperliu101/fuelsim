#pragma once
#include "fuelsim/material.hpp"
#include <array>
namespace fuelsim {
enum class InelasticBehavior {
    elastic,
    norton_creep,
    j2_plasticity,
    norton_creep_j2_plasticity,
};
struct NortonCreepProperties final {
    double coefficient, reference_stress, stress_exponent;
    double coefficient_temperature_coefficient = 0.0, reference_stress_temperature_coefficient = 0.0,
           stress_exponent_temperature_coefficient = 0.0;
};
struct J2PlasticityProperties final {
    double yield_stress, isotropic_hardening_modulus;
    double yield_stress_temperature_coefficient = 0.0, hardening_temperature_coefficient = 0.0;
};
struct ActiveNortonCreepProperties final {
    adlite::Scalar coefficient, reference_stress, stress_exponent;
};
struct ActiveJ2PlasticityProperties final {
    adlite::Scalar yield_stress, isotropic_hardening_modulus;
};
struct TransientInelasticProperties final {
    double density, specific_heat;
    InelasticBehavior behavior;
    NortonCreepProperties creep;
    J2PlasticityProperties plasticity;
};
struct MaterialPointState final {
    std::array<double, 4> elastic_strain{};
    std::array<double, 4> plastic_strain{};
    std::array<double, 4> creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0;
};
struct MaterialPointTrialState final {
    std::array<adlite::Scalar, 4> elastic_strain{};
    std::array<adlite::Scalar, 4> plastic_strain{};
    std::array<adlite::Scalar, 4> creep_strain{};
    adlite::Scalar equivalent_plastic_strain{0.0}, equivalent_creep_strain{0.0};
};
struct InelasticStressResponse final {
    AxisymmetricStress stress;
    MaterialPointTrialState trial_state;
};
class IsotropicInelasticMaterial final {
  public:
    IsotropicInelasticMaterial(
        ThermoelasticProperties thermoelastic_properties, TransientInelasticProperties properties);
    const TransientInelasticProperties& properties() const noexcept { return _properties; }
    adlite::Scalar conductivity(
        const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    adlite::Scalar heat_capacity(
        const adlite::Scalar& temperature, double time = 0.0, double radius = 0.0, double axial_coordinate = 0.0) const;
    InelasticStressResponse response(const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz, const adlite::Scalar& temperature,
        double time_step, const MaterialPointState& committed, double time = 0.0, double radius = 0.0,
        double axial_coordinate = 0.0) const;
    InelasticStressResponse incremental_response(const adlite::Scalar& strain_increment_rr,
        const adlite::Scalar& strain_increment_zz, const adlite::Scalar& strain_increment_hoop,
        const adlite::Scalar& strain_increment_rz, const AxisymmetricRotation& rotation,
        const adlite::Scalar& temperature, double committed_temperature, double time_step,
        const MaterialPointState& committed, double time = 0.0, double radius = 0.0,
        double axial_coordinate = 0.0) const;
    static MaterialPointState state_values(const MaterialPointTrialState& trial_state);

  private:
    ActiveNortonCreepProperties active_creep_properties(const adlite::Scalar& temperature) const;
    ActiveJ2PlasticityProperties active_plasticity_properties(const adlite::Scalar& temperature) const;
    IsotropicThermoelasticMaterial _thermoelastic_material;
    TransientInelasticProperties _properties;
};
} // namespace fuelsim
