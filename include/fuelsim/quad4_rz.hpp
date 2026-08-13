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
class Quad4RzTransientKernel final {
  public:
    Quad4RzTransientKernel(
        IsotropicInelasticMaterial material, double volumetric_heat_source, StrainFormulation strain_formulation);
    double volumetric_heat_source() const noexcept { return _volumetric_heat_source; }
    double heat_capacity(double temperature, double radius, double axial_coordinate) const;
    void set_volumetric_heat_source(double value) noexcept { _volumetric_heat_source = value; }
    void set_time(double value) noexcept { _time = value; }
    LocalResidual residual(const Quad4RzGeometry& geometry, const LocalValues& current_state,
        const LocalValues& committed_state, const Quad4MaterialHistory& committed_material, double time_step) const;
    LocalSystem linearize(const Quad4RzGeometry& geometry, const LocalValues& current_state,
        const LocalValues& committed_state, const Quad4MaterialHistory& committed_material, double time_step) const;
    Quad4MaterialHistory trial_state_values(const Quad4RzGeometry& geometry, const LocalValues& converged_state,
        const LocalValues& committed_state, const Quad4MaterialHistory& committed_material, double time_step) const;
    std::array<AxisymmetricStressValues, 4> stress_values(const Quad4RzGeometry& geometry, const LocalValues& state,
        const LocalValues& committed_state, const Quad4MaterialHistory& committed_material, double time_step) const;

  private:
    void residual_ad(const Quad4RzGeometry& geometry, const LocalAdValues& current_state,
        const LocalValues& committed_state, const Quad4MaterialHistory& committed_material, double time_step,
        LocalAdValues& residual) const;
    IsotropicInelasticMaterial _material;
    double _volumetric_heat_source, _time;
    StrainFormulation _strain_formulation;
};
struct Line2RzBoundaryGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};
enum class TractionComponent { radial, axial };
Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(
    const std::array<RzPoint, 2>& coordinates, const std::array<std::size_t, 2>& local_nodes);
class Line2RzBoundaryKernel final {
  public:
    Line2RzBoundaryKernel(double pressure, bool use_displaced_geometry);
    Line2RzBoundaryKernel(TractionComponent component, double traction, bool use_displaced_geometry);
    Line2RzBoundaryKernel(double heat_transfer_coefficient, double ambient_temperature);
    void set_load(double value) noexcept { _load = value; }
    void set_convection(double coefficient, double ambient) noexcept {
        _load = coefficient;
        _ambient = ambient;
    }
    LocalResidual residual(const Line2RzBoundaryGeometry& geometry, const LocalValues& state) const;
    LocalSystem linearize(const Line2RzBoundaryGeometry& geometry, const LocalValues& state) const;

  private:
    void residual_ad(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    void mechanical_residual(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    void convection_residual(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    enum class Kind { pressure, traction, convection };
    Kind _kind;
    TractionComponent _component;
    double _load, _ambient;
    bool _use_displaced_geometry;
};
} // namespace fuelsim
