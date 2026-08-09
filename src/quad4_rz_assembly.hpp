#ifndef FUELSIM_QUAD4_RZ_ASSEMBLY_HPP
#define FUELSIM_QUAD4_RZ_ASSEMBLY_HPP

#include "fuelsim/quad4_rz_kinematics.hpp"

#include <array>
#include <cstddef>

namespace fuelsim::quad4_rz_detail {

inline adlite::Scalar
interpolate(const std::array<double, quad4_node_count>& coefficients,
            const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node)
        result += coefficients[node] * state[offset + node];
    return result;
}

inline double
interpolate(const std::array<double, quad4_node_count>& coefficients,
            const LocalValues& state, std::size_t offset) {
    double result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node)
        result += coefficients[node] * state[offset + node];
    return result;
}

inline LocalAdValues passive_state(const LocalValues& state) {
    LocalAdValues result{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        result[dof] = state[dof];
    return result;
}

inline void add_mechanical_point_residual(
    const RzQuadraturePoint& point, const AxisymmetricKinematics& kinematics,
    const AxisymmetricStress& stress, LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        residual[4 + node] +=
            kinematics.weighted_measure *
            (stress.rr * kinematics.gradient_r[node] +
             stress.hoop * point.shape[node] / kinematics.radius +
             stress.rz * kinematics.gradient_z[node]);

        residual[8 + node] += kinematics.weighted_measure *
                              (stress.zz * kinematics.gradient_z[node] +
                               stress.rz * kinematics.gradient_r[node]);
    }
}

inline void add_steady_point_residual(
    const RzQuadraturePoint& point,
    const adlite::Scalar& gradient_temperature_r,
    const adlite::Scalar& gradient_temperature_z,
    const AxisymmetricKinematics& kinematics,
    const adlite::Scalar& conductivity, double volumetric_heat_source,
    const AxisymmetricStress& stress, LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        residual[node] +=
            point.weighted_measure *
            (conductivity * (point.gradient_r[node] * gradient_temperature_r +
                             point.gradient_z[node] * gradient_temperature_z) -
             volumetric_heat_source * point.shape[node]);
    }
    add_mechanical_point_residual(point, kinematics, stress, residual);
}

inline void add_transient_point_residual(
    const RzQuadraturePoint& point,
    const adlite::Scalar& gradient_temperature_r,
    const adlite::Scalar& gradient_temperature_z,
    const AxisymmetricKinematics& kinematics, double heat_capacity,
    const adlite::Scalar& temperature_rate, const adlite::Scalar& conductivity,
    double volumetric_heat_source, const AxisymmetricStress& stress,
    LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        residual[node] +=
            point.weighted_measure *
            (heat_capacity * point.shape[node] * temperature_rate +
             conductivity * (point.gradient_r[node] * gradient_temperature_r +
                             point.gradient_z[node] * gradient_temperature_z) -
             volumetric_heat_source * point.shape[node]);
    }
    add_mechanical_point_residual(point, kinematics, stress, residual);
}

} // namespace fuelsim::quad4_rz_detail

#endif
