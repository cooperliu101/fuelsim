#pragma once
#include "c3d20_types.hpp"
#include "matrix3.hpp"

namespace fuelsim::c3d20_detail {
using cartesian_detail::ActiveMatrix3;
using cartesian_detail::Matrix3;

struct Hex20KinematicsValues final {
    SymmetricTensor3Values strain_increment{};
    CartesianRotation rotation{};
    std::array<std::array<double, 3>, 20> current_gradient{};
    double current_weighted_measure = 0.0;
};

struct Hex20Kinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, 20> current_gradient;
    Matrix3 current_inverse_values{}, midpoint_inverse_values{};
    adlite::Scalar current_weighted_measure{0.0};
};

struct Hex20SourceMeasureValues final {
    Matrix3 measure_derivative{};
    double weighted_measure = 0.0;
};

struct Hex20FiniteThermalKinematicsValues final {
    std::array<std::array<double, 3>, 8> midpoint_temperature_gradient{};
    double current_weighted_measure = 0.0;
};

ActiveMatrix3 displacement_gradient(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state);

Hex20KinematicsValues evaluate_kinematics_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation);

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point,
    const ActiveMatrix3& gradient,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation);

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalAdValues& state,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation);

Hex20SourceMeasureValues evaluate_source_measure_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state);

std::array<std::array<double, 3>, 8> temperature_shape_gradients(const Hex20MechanicalQuadraturePoint& point,
    const Matrix3& inverse_map);

Hex20FiniteThermalKinematicsValues evaluate_finite_thermal_kinematics_values(
    const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state);
} // namespace fuelsim::c3d20_detail
