#include "fuelsim/quad4_rz_thermoelastic.hpp"

#include "quad4_rz_assembly.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

struct ThermoelasticPointResponse final {
    adlite::Scalar temperature;
    adlite::Scalar gradient_temperature_r;
    adlite::Scalar gradient_temperature_z;
    AxisymmetricKinematics kinematics;
    AxisymmetricStress stress;
};

ThermoelasticPointResponse
point_response(const RzQuadraturePoint& point, const LocalAdValues& state,
               const IsotropicThermoelasticMaterial& material,
               StrainFormulation strain_formulation) {
    const adlite::Scalar temperature =
        quad4_rz_detail::interpolate(point.shape, state, 0);
    const AxisymmetricKinematics kinematics =
        evaluate_axisymmetric_kinematics(point, state, strain_formulation);
    AxisymmetricStress stress = material.stress(
        kinematics.strain_rr, kinematics.strain_zz, kinematics.strain_hoop,
        kinematics.strain_rz, temperature);
    if (strain_formulation == StrainFormulation::finite)
        stress = rotate_axisymmetric_tensor(stress, kinematics.rotation);
    return {
        temperature,
        quad4_rz_detail::interpolate(point.gradient_r, state, 0),
        quad4_rz_detail::interpolate(point.gradient_z, state, 0),
        kinematics,
        stress,
    };
}

} // namespace

Quad4RzThermoelasticKernel::Quad4RzThermoelasticKernel(
    IsotropicThermoelasticMaterial material, double volumetric_heat_source,
    StrainFormulation strain_formulation)
    : _material(material), _volumetric_heat_source(0.0),
      _strain_formulation(strain_formulation) {
    set_volumetric_heat_source(volumetric_heat_source);
}

double Quad4RzThermoelasticKernel::volumetric_heat_source() const noexcept {
    return _volumetric_heat_source;
}

void Quad4RzThermoelasticKernel::set_volumetric_heat_source(
    double volumetric_heat_source) {
    if (!std::isfinite(volumetric_heat_source) ||
        !(volumetric_heat_source >= 0.0))
        throw std::invalid_argument(
            "Quad4RzThermoelasticKernel volumetric heat source must be finite "
            "and nonnegative");
    _volumetric_heat_source = volumetric_heat_source;
}

LocalResidual
Quad4RzThermoelasticKernel::residual(const Quad4RzGeometry& geometry,
                                     const LocalValues& state) const {
    const LocalAdValues passive_state = quad4_rz_detail::passive_state(state);

    LocalAdValues passive_residual{};
    residual_ad(geometry, passive_state, passive_residual);

    LocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = passive_residual[row].value();
    return result;
}

LocalSystem
Quad4RzThermoelasticKernel::linearize(const Quad4RzGeometry& geometry,
                                      const LocalValues& state) const {
    LocalAdValues active_state{};
    adlite::seed_identity(state.data(), state.size(), active_state.data());

    LocalAdValues active_residual{};
    residual_ad(geometry, active_state, active_residual);

    LocalSystem result{};
    adlite::extract_jacobian(active_residual.data(), active_residual.size(),
                             active_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

std::array<AxisymmetricStressValues, 4>
Quad4RzThermoelasticKernel::stress_values(const Quad4RzGeometry& geometry,
                                          const LocalValues& state) const {
    const LocalAdValues passive_state = quad4_rz_detail::passive_state(state);

    std::array<AxisymmetricStressValues, 4> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const AxisymmetricStress& stress =
            point_response(geometry.points[q], passive_state, _material,
                           _strain_formulation)
                .stress;
        result[q] = {stress.rr.value(), stress.zz.value(), stress.hoop.value(),
                     stress.rz.value()};
    }
    return result;
}

void Quad4RzThermoelasticKernel::residual_ad(const Quad4RzGeometry& geometry,
                                             const LocalAdValues& state,
                                             LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));

    for (const RzQuadraturePoint& point : geometry.points) {
        const ThermoelasticPointResponse response =
            point_response(point, state, _material, _strain_formulation);
        const adlite::Scalar conductivity =
            _material.conductivity(response.temperature);
        quad4_rz_detail::add_steady_point_residual(
            point, response.gradient_temperature_r,
            response.gradient_temperature_z, response.kinematics, conductivity,
            _volumetric_heat_source, response.stress, residual);
    }
}

} // namespace fuelsim
