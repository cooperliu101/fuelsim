#pragma once
#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/mesh.hpp"
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
struct LocalSystem final {
    LocalResidual residual;
    LocalJacobian jacobian;
};
constexpr std::size_t quad4_node_count = 4;
constexpr std::size_t quad4_local_dof_count = local_dof_count;
constexpr std::size_t quad4_jacobian_size = local_jacobian_size;
using Quad4Coordinates = std::array<RzPoint, quad4_node_count>;
struct RzQuadraturePoint final {
    std::array<double, quad4_node_count> shape, gradient_r, gradient_z;
    double radius, axial_coordinate, weighted_measure;
};
struct Quad4RzGeometry final {
    std::array<RzQuadraturePoint, 4> points;
};
Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates);
enum class StrainFormulation { small, finite };
struct AxisymmetricKinematics final {
    std::array<adlite::Scalar, quad4_node_count> gradient_r, gradient_z;
    adlite::Scalar radius, weighted_measure, strain_rr, strain_zz, strain_hoop, strain_rz;
    AxisymmetricRotation rotation;
};
AxisymmetricKinematics evaluate_axisymmetric_kinematics(
    const RzQuadraturePoint& point, const LocalAdValues& state, StrainFormulation strain_formulation);
AxisymmetricKinematics evaluate_axisymmetric_incremental_kinematics(const RzQuadraturePoint& point,
    const LocalAdValues& current_state, const LocalValues& committed_state, StrainFormulation strain_formulation);
struct Quad4RzThermoelasticData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
    StrainFormulation strain_formulation;
};
LocalResidual compute_quad4_rz_thermoelastic_residual(
    const Quad4RzThermoelasticData& data, const Quad4RzGeometry& geometry, const LocalValues& state);
LocalSystem compute_quad4_rz_thermoelastic_system(
    const Quad4RzThermoelasticData& data, const Quad4RzGeometry& geometry, const LocalValues& state);
std::array<AxisymmetricStressValues, 4> compute_quad4_rz_thermoelastic_stress(
    const Quad4RzThermoelasticData& data, const Quad4RzGeometry& geometry, const LocalValues& state);
using Quad4MaterialHistory = std::array<MaterialPointState, 4>;
struct Quad4RzTransientData final {
    IsotropicInelasticMaterial material;
    double volumetric_heat_source = 0.0, time = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
};
double compute_quad4_rz_transient_heat_capacity(
    const Quad4RzTransientData& data, double temperature, double radius, double axial_coordinate);
LocalResidual compute_quad4_rz_transient_residual(const Quad4RzTransientData& data, const Quad4RzGeometry& geometry,
    const LocalValues& current_state, const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step);
LocalSystem compute_quad4_rz_transient_system(const Quad4RzTransientData& data, const Quad4RzGeometry& geometry,
    const LocalValues& current_state, const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step);
Quad4MaterialHistory compute_quad4_rz_transient_trial_state(const Quad4RzTransientData& data,
    const Quad4RzGeometry& geometry, const LocalValues& converged_state, const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step);
std::array<AxisymmetricStressValues, 4> compute_quad4_rz_transient_stress(const Quad4RzTransientData& data,
    const Quad4RzGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step);
struct Line2RzBoundaryGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};
enum class TractionComponent { radial, axial };
Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(
    const std::array<RzPoint, 2>& coordinates, const std::array<std::size_t, 2>& local_nodes);
enum class Line2RzBoundaryKind { pressure, traction, convection };
struct Line2RzBoundaryData final {
    Line2RzBoundaryKind kind;
    TractionComponent component;
    double load, ambient;
    bool use_displaced_geometry;
};
Line2RzBoundaryData make_line2_rz_pressure_data(double pressure, bool use_displaced_geometry);
Line2RzBoundaryData make_line2_rz_traction_data(
    TractionComponent component, double traction, bool use_displaced_geometry);
Line2RzBoundaryData make_line2_rz_convection_data(double heat_transfer_coefficient, double ambient_temperature);
LocalResidual compute_line2_rz_boundary_residual(
    const Line2RzBoundaryData& data, const Line2RzBoundaryGeometry& geometry, const LocalValues& state);
LocalSystem compute_line2_rz_boundary_system(
    const Line2RzBoundaryData& data, const Line2RzBoundaryGeometry& geometry, const LocalValues& state);
} // namespace fuelsim
