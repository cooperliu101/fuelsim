#pragma once
#include "fuelsim/material.hpp"
#include <array>
namespace fuelsim {
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
    explicit IsotropicInelasticMaterial(ThermoelasticProperties thermoelastic_properties);
    adlite::Scalar conductivity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    adlite::Scalar heat_capacity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    InelasticStressResponse response(const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz, const adlite::Scalar& temperature,
        double time_step, const MaterialPointState& committed, MaterialFunctionContext context = {}) const;
    InelasticStressResponse incremental_response(const adlite::Scalar& strain_increment_rr,
        const adlite::Scalar& strain_increment_zz, const adlite::Scalar& strain_increment_hoop,
        const adlite::Scalar& strain_increment_rz, const AxisymmetricRotation& rotation,
        const adlite::Scalar& temperature, double committed_temperature, double time_step,
        const MaterialPointState& committed, MaterialFunctionContext context = {}) const;
    static MaterialPointState state_values(const MaterialPointTrialState& trial_state);

  private:
    IsotropicThermoelasticMaterial _thermoelastic_material;
};
} // namespace fuelsim
