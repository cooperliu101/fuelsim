#ifndef FUELSIM_QUAD4_RZ_THERMOELASTIC_HPP
#define FUELSIM_QUAD4_RZ_THERMOELASTIC_HPP

#include "fuelsim/material.hpp"
#include "fuelsim/quad4_rz_kinematics.hpp"

#include <array>

namespace fuelsim {

class Quad4RzThermoelasticKernel final {
  public:
    Quad4RzThermoelasticKernel(IsotropicThermoelasticMaterial material, double volumetric_heat_source,
                               StrainFormulation strain_formulation);

    double volumetric_heat_source() const noexcept;
    void set_volumetric_heat_source(double volumetric_heat_source) noexcept;
    void set_time(double time) noexcept;

    LocalResidual residual(const Quad4RzGeometry& geometry, const LocalValues& state) const;

    LocalSystem linearize(const Quad4RzGeometry& geometry, const LocalValues& state) const;

    std::array<AxisymmetricStressValues, 4> stress_values(const Quad4RzGeometry& geometry,
                                                          const LocalValues& state) const;

  private:
    void residual_ad(const Quad4RzGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;

    IsotropicThermoelasticMaterial _material;
    double _volumetric_heat_source;
    double _time;
    StrainFormulation _strain_formulation;
};

} // namespace fuelsim

#endif
