#pragma once
#include "fuelsim/core/kinematics.hpp"
#include "fuelsim/core/material.hpp"
#include "fuelsim/core/mesh.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>

namespace fuelsim {
constexpr std::size_t local_dof_count = 12;
constexpr std::size_t local_jacobian_size = local_dof_count * local_dof_count;
using LocalDofs = std::array<std::size_t, local_dof_count>;
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

struct Quad4RzData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source = 0.0, time = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    RzElementFormulation element_formulation = RzElementFormulation::quad4;
    double initial_temperature = 600.0;
};

LocalResidual compute_quad4_rz_thermoelastic(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian = nullptr);
std::array<AxisymmetricStressValues, 4> compute_quad4_rz_thermoelastic_stress(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state);
using Quad4MaterialHistory = std::array<MaterialPointState, 4>;
LocalResidual compute_quad4_rz_transient(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& current_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material,
    double time_step,
    LocalJacobian* jacobian = nullptr,
    bool include_thermal_time_term = true);
Quad4MaterialHistory compute_quad4_rz_transient_update(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& converged_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material,
    double time_step);

struct Line2RzBoundaryGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};
enum class TractionComponent { radial, axial };
Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes);
enum class Line2RzBoundaryKind { pressure, traction, convection };

struct Line2RzBoundaryData final {
    Line2RzBoundaryKind kind;
    TractionComponent component;
    double load, ambient;
    bool use_displaced_geometry;
};

LocalResidual compute_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian = nullptr);
} // namespace fuelsim
