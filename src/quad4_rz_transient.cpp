#include "fuelsim/quad4_rz_transient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

struct PointFields final {
    adlite::Scalar temperature;
    adlite::Scalar gradient_temperature_r;
    adlite::Scalar gradient_temperature_z;
    AxisymmetricKinematics kinematics;
};

adlite::Scalar interpolate(const std::array<double, quad4_node_count>& values,
                           const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node)
        result += values[node] * state[offset + node];
    return result;
}

double interpolate_committed_temperature(
    const std::array<double, quad4_node_count>& shape,
    const LocalValues& committed_state) {
    double result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node)
        result += shape[node] * committed_state[node];
    return result;
}

PointFields point_fields(const RzQuadraturePoint& point,
                         const LocalAdValues& state,
                         const LocalValues& committed_state,
                         StrainFormulation strain_formulation) {
    const adlite::Scalar temperature = interpolate(point.shape, state, 0);

    return {
        temperature,
        interpolate(point.gradient_r, state, 0),
        interpolate(point.gradient_z, state, 0),
        evaluate_axisymmetric_incremental_kinematics(
            point, state, committed_state, strain_formulation),
    };
}

void validate_time_step(double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument(
            "Quad4RzTransientKernel time_step must be finite and positive");
}

void validate_committed_state(const LocalValues& committed_state) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        const double temperature = committed_state[node];
        if (!std::isfinite(temperature) || !(temperature > 0.0))
            throw std::invalid_argument(
                "Quad4RzTransientKernel committed temperatures must be finite "
                "and positive");
    }
}

double volumetric_heat_capacity(const IsotropicInelasticMaterial& material) {
    const TransientInelasticProperties& properties = material.properties();
    const double value = properties.density * properties.specific_heat;
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error(
            "Quad4RzTransientKernel volumetric heat capacity must be finite "
            "and positive");
    return value;
}

} // namespace

Quad4RzTransientKernel::Quad4RzTransientKernel(
    IsotropicInelasticMaterial material, double volumetric_heat_source,
    StrainFormulation strain_formulation)
    : _material(material), _volumetric_heat_source(0.0),
      _strain_formulation(strain_formulation) {
    (void)volumetric_heat_capacity(_material);
    set_volumetric_heat_source(volumetric_heat_source);
}

double Quad4RzTransientKernel::volumetric_heat_source() const noexcept {
    return _volumetric_heat_source;
}

void Quad4RzTransientKernel::set_volumetric_heat_source(
    double volumetric_heat_source) {
    if (!std::isfinite(volumetric_heat_source) ||
        !(volumetric_heat_source >= 0.0))
        throw std::invalid_argument(
            "Quad4RzTransientKernel volumetric heat source must be finite and "
            "nonnegative");
    _volumetric_heat_source = volumetric_heat_source;
}

LocalResidual Quad4RzTransientKernel::residual(
    const Quad4RzGeometry& geometry, const LocalValues& current_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step) const {
    validate_time_step(time_step);
    validate_committed_state(committed_state);

    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < current_state.size(); ++dof)
        passive_state[dof] = current_state[dof];

    LocalAdValues passive_residual{};
    residual_ad(geometry, passive_state, committed_state,
                committed_material, time_step, passive_residual);

    LocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = passive_residual[row].value();
    return result;
}

LocalSystem Quad4RzTransientKernel::linearize(
    const Quad4RzGeometry& geometry, const LocalValues& current_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step) const {
    validate_time_step(time_step);
    validate_committed_state(committed_state);

    LocalAdValues active_state{};
    adlite::seed_identity(current_state.data(), current_state.size(),
                          active_state.data());

    LocalAdValues active_residual{};
    residual_ad(geometry, active_state, committed_state,
                committed_material, time_step, active_residual);

    LocalSystem result{};
    adlite::extract_jacobian(active_residual.data(), active_residual.size(),
                             active_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

Quad4MaterialHistory Quad4RzTransientKernel::trial_state_values(
    const Quad4RzGeometry& geometry, const LocalValues& converged_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step) const {
    validate_time_step(time_step);

    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < converged_state.size(); ++dof)
        passive_state[dof] = converged_state[dof];

    Quad4MaterialHistory result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const PointFields fields =
            point_fields(geometry.points[q], passive_state, committed_state,
                         _strain_formulation);
        if (!std::isfinite(fields.temperature.value()) ||
            !(fields.temperature.value() > 0.0))
            throw std::domain_error(
                "Quad4RzTransientKernel trial temperature must be finite and "
                "positive");

        const double old_temperature = interpolate_committed_temperature(
            geometry.points[q].shape, committed_state);
        const auto response =
            _strain_formulation == StrainFormulation::finite
                ? _material.incremental_response(
                      fields.kinematics.strain_rr,
                      fields.kinematics.strain_zz,
                      fields.kinematics.strain_hoop,
                      fields.kinematics.strain_rz,
                      fields.kinematics.rotation, fields.temperature,
                      old_temperature, time_step, committed_material[q])
                : _material.response(
                      fields.kinematics.strain_rr,
                      fields.kinematics.strain_zz,
                      fields.kinematics.strain_hoop,
                      fields.kinematics.strain_rz, fields.temperature,
                      time_step, committed_material[q]);
        result[q] =
            IsotropicInelasticMaterial::state_values(response.trial_state);
    }
    return result;
}

std::array<AxisymmetricStressValues, 4> Quad4RzTransientKernel::stress_values(
    const Quad4RzGeometry& geometry, const LocalValues& state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step) const {
    validate_time_step(time_step);

    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];

    std::array<AxisymmetricStressValues, 4> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const PointFields fields =
            point_fields(geometry.points[q], passive_state, committed_state,
                         _strain_formulation);
        const double old_temperature = interpolate_committed_temperature(
            geometry.points[q].shape, committed_state);
        const InelasticStressResponse response =
            _strain_formulation == StrainFormulation::finite
                ? _material.incremental_response(
                      fields.kinematics.strain_rr,
                      fields.kinematics.strain_zz,
                      fields.kinematics.strain_hoop,
                      fields.kinematics.strain_rz,
                      fields.kinematics.rotation, fields.temperature,
                      old_temperature, time_step, committed_material[q])
                : _material.response(
                      fields.kinematics.strain_rr,
                      fields.kinematics.strain_zz,
                      fields.kinematics.strain_hoop,
                      fields.kinematics.strain_rz, fields.temperature,
                      time_step, committed_material[q]);
        result[q] = {
            response.stress.rr.value(),
            response.stress.zz.value(),
            response.stress.hoop.value(),
            response.stress.rz.value(),
        };
    }
    return result;
}

void Quad4RzTransientKernel::residual_ad(
    const Quad4RzGeometry& geometry, const LocalAdValues& current_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step,
    LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));
    const double heat_capacity = volumetric_heat_capacity(_material);

    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const RzQuadraturePoint& point = geometry.points[q];
        const PointFields fields =
            point_fields(point, current_state, committed_state,
                         _strain_formulation);
        const double old_temperature = interpolate_committed_temperature(
            point.shape, committed_state);
        const adlite::Scalar temperature_rate =
            (fields.temperature - old_temperature) / time_step;
        const adlite::Scalar conductivity =
            _material.conductivity(fields.temperature);
        const auto response =
            _strain_formulation == StrainFormulation::finite
                ? _material.incremental_response(
                      fields.kinematics.strain_rr,
                      fields.kinematics.strain_zz,
                      fields.kinematics.strain_hoop,
                      fields.kinematics.strain_rz,
                      fields.kinematics.rotation, fields.temperature,
                      old_temperature, time_step, committed_material[q])
                : _material.response(
                      fields.kinematics.strain_rr,
                      fields.kinematics.strain_zz,
                      fields.kinematics.strain_hoop,
                      fields.kinematics.strain_rz, fields.temperature,
                      time_step, committed_material[q]);
        const auto& stress = response.stress;

        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            residual[node] +=
                point.weighted_measure *
                (heat_capacity * point.shape[node] * temperature_rate +
                 conductivity *
                     (point.gradient_r[node] * fields.gradient_temperature_r +
                      point.gradient_z[node] * fields.gradient_temperature_z) -
                 _volumetric_heat_source * point.shape[node]);

            residual[4 + node] +=
                fields.kinematics.weighted_measure *
                (stress.rr * fields.kinematics.gradient_r[node] +
                 stress.hoop * point.shape[node] /
                     fields.kinematics.radius +
                 stress.rz * fields.kinematics.gradient_z[node]);

            residual[8 + node] +=
                fields.kinematics.weighted_measure *
                (stress.zz * fields.kinematics.gradient_z[node] +
                 stress.rz * fields.kinematics.gradient_r[node]);
        }
    }
}

} // namespace fuelsim
