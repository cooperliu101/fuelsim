#pragma once
#include "fuelsim/elements/coordinates.hpp"
#include "fuelsim/elements/kinematics.hpp"
#include "fuelsim/elements/material.hpp"
#include <array>
#include <cstddef>

namespace fuelsim {
constexpr std::size_t local_dof_count = 12;
constexpr std::size_t local_jacobian_size = local_dof_count * local_dof_count;
using LocalValues = std::array<double, local_dof_count>;
using LocalResidual = std::array<double, local_dof_count>;
using LocalJacobian = std::array<double, local_jacobian_size>;
using LocalAdValues = std::array<adlite::Scalar, local_dof_count>;
constexpr std::size_t quad4_node_count = 4;
constexpr std::size_t quad4_local_dof_count = local_dof_count;
using Quad4Coordinates = std::array<RzPoint, quad4_node_count>;

struct RzQuadraturePoint final {
    std::array<double, quad4_node_count> shape, gradient_r, gradient_z;
    double radius, axial_coordinate, weighted_measure;
};

struct Quad4RzGeometry final {
    std::array<RzQuadraturePoint, 4> points;
    Quad4Coordinates coordinates{};
};

Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates);

struct AxisymmetricKinematics final {
    std::array<adlite::Scalar, quad4_node_count> gradient_r, gradient_z;
    adlite::Scalar radius, weighted_measure, strain_rr, strain_zz, strain_hoop, strain_rz;
    AxisymmetricRotation rotation;
    adlite::Scalar midpoint_weighted_measure = 0.0;
};

AxisymmetricKinematics evaluate_axisymmetric_kinematics(const RzQuadraturePoint& point,
    const LocalAdValues& state,
    StrainFormulation strain_formulation);
AxisymmetricKinematics evaluate_axisymmetric_incremental_kinematics(const RzQuadraturePoint& point,
    const LocalAdValues& current_state,
    const LocalValues& committed_state,
    StrainFormulation strain_formulation);

using Quad4MaterialHistory = std::array<MaterialPointState, 4>;
AxisymmetricKinematics evaluate_axisymmetric_kinematics_from_point(const RzQuadraturePoint& point,
    const adlite::Scalar& radial_displacement,
    const adlite::Scalar& displacement_gradient_rr,
    const adlite::Scalar& displacement_gradient_rz,
    const adlite::Scalar& displacement_gradient_zr,
    const adlite::Scalar& displacement_gradient_zz,
    const LocalValues& committed_state,
    StrainFormulation strain_formulation,
    bool cax4t = false);

struct AxisymmetricStressTangent final {
    AxisymmetricStressValues stress;
    std::array<std::array<double, 4>, 4> tangent{};
    std::array<double, 4> thermal{};
};

AxisymmetricStressTangent evaluate_axisymmetric_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 4>& fed_strain,
    double temperature,
    double time_step,
    const MaterialPointState* committed_material,
    MaterialFunctionContext context);
} // namespace fuelsim
