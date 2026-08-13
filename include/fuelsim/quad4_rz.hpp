#pragma once
#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/mesh.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <variant>
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
class Quad4RzThermoelasticKernel final {
  public:
    Quad4RzThermoelasticKernel(
        IsotropicThermoelasticMaterial material, double volumetric_heat_source, StrainFormulation strain_formulation);
    double volumetric_heat_source() const noexcept { return _volumetric_heat_source; }
    void set_volumetric_heat_source(double value) noexcept { _volumetric_heat_source = value; }
    void set_time(double value) noexcept { _time = value; }
    LocalResidual residual(const Quad4RzGeometry& geometry, const LocalValues& state) const;
    LocalSystem linearize(const Quad4RzGeometry& geometry, const LocalValues& state) const;
    std::array<AxisymmetricStressValues, 4> stress_values(
        const Quad4RzGeometry& geometry, const LocalValues& state) const;

  private:
    void residual_ad(const Quad4RzGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    IsotropicThermoelasticMaterial _material;
    double _volumetric_heat_source, _time;
    StrainFormulation _strain_formulation;
};
using Quad4MaterialHistory = std::array<MaterialPointState, 4>;
class Quad4RzTransientKernel final {
  public:
    Quad4RzTransientKernel(
        IsotropicInelasticMaterial material, double volumetric_heat_source, StrainFormulation strain_formulation);
    double volumetric_heat_source() const noexcept { return _volumetric_heat_source; }
    const TransientInelasticProperties& properties() const noexcept { return _material.properties(); }
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
struct ConvectionProperties final {
    double heat_transfer_coefficient, ambient_temperature;
};
struct PressureProperties final {
    double pressure;
    bool use_displaced_geometry;
};
enum class TractionComponent { radial, axial };
struct TractionProperties final {
    TractionComponent component;
    double traction;
    bool use_displaced_geometry;
};
Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(
    const std::array<RzPoint, 2>& coordinates, const std::array<std::size_t, 2>& local_nodes);
class Line2RzBoundaryKernel final {
  public:
    explicit Line2RzBoundaryKernel(PressureProperties properties) : _properties(properties) {}
    explicit Line2RzBoundaryKernel(TractionProperties properties) : _properties(properties) {}
    explicit Line2RzBoundaryKernel(ConvectionProperties properties) : _properties(properties) {}
    const PressureProperties& pressure_properties() const { return std::get<PressureProperties>(_properties); }
    const TractionProperties& traction_properties() const { return std::get<TractionProperties>(_properties); }
    void set_properties(PressureProperties value) noexcept { _properties = value; }
    void set_properties(TractionProperties value) noexcept { _properties = value; }
    void set_properties(ConvectionProperties value) noexcept { _properties = value; }
    LocalResidual residual(const Line2RzBoundaryGeometry& geometry, const LocalValues& state) const;
    LocalSystem linearize(const Line2RzBoundaryGeometry& geometry, const LocalValues& state) const;

  private:
    void residual_ad(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    void pressure_residual(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    void traction_residual(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    void convection_residual(
        const Line2RzBoundaryGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
    std::variant<PressureProperties, TractionProperties, ConvectionProperties> _properties;
};
} // namespace fuelsim
