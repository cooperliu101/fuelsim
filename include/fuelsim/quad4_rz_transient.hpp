#ifndef FUELSIM_QUAD4_RZ_TRANSIENT_HPP
#define FUELSIM_QUAD4_RZ_TRANSIENT_HPP

#include <array>

#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/quad4_rz.hpp"

namespace fuelsim {

using Quad4TemperatureHistory = std::array<double, quad4_node_count>;
using Quad4MaterialHistory = std::array<MaterialPointState, 4>;

class Quad4RzTransientKernel final {
  public:
    Quad4RzTransientKernel(IsotropicInelasticMaterial material,
                           double volumetric_heat_source);

    double volumetric_heat_source() const noexcept;
    void set_volumetric_heat_source(double volumetric_heat_source);

    LocalResidual residual(const Quad4RzGeometry& geometry,
                           const LocalValues& current_state,
                           const Quad4TemperatureHistory& committed_temperature,
                           const Quad4MaterialHistory& committed_material,
                           double time_step) const;

    LocalSystem linearize(const Quad4RzGeometry& geometry,
                          const LocalValues& current_state,
                          const Quad4TemperatureHistory& committed_temperature,
                          const Quad4MaterialHistory& committed_material,
                          double time_step) const;

    Quad4MaterialHistory trial_state_values(
        const Quad4RzGeometry& geometry, const LocalValues& converged_state,
        const Quad4MaterialHistory& committed_material, double time_step) const;

    std::array<AxisymmetricStressValues, 4>
    stress_values(const Quad4RzGeometry& geometry, const LocalValues& state,
                  const Quad4MaterialHistory& committed_material,
                  double time_step) const;

  private:
    void residual_ad(const Quad4RzGeometry& geometry,
                     const LocalAdValues& current_state,
                     const Quad4TemperatureHistory& committed_temperature,
                     const Quad4MaterialHistory& committed_material,
                     double time_step, LocalAdValues& residual) const;

    IsotropicInelasticMaterial _material;
    double _volumetric_heat_source;
};

} // namespace fuelsim

#endif
