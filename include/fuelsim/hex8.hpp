#pragma once
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
namespace fuelsim {
inline constexpr std::size_t hex8_node_count = 8;
inline constexpr std::size_t hex8_local_dof_count = 32;
inline constexpr std::size_t quad4_face_node_count = 4;
inline constexpr std::size_t quad4_face_local_dof_count = 16;
using Hex8LocalDofs = std::array<std::size_t, hex8_local_dof_count>;
using Hex8LocalValues = std::array<double, hex8_local_dof_count>;
using Hex8LocalResidual = std::array<double, hex8_local_dof_count>;
using Hex8LocalJacobian = std::array<double, hex8_local_dof_count * hex8_local_dof_count>;
using Hex8LocalAdValues = std::array<adlite::Scalar, hex8_local_dof_count>;
struct Hex8LocalSystem final {
    Hex8LocalResidual residual;
    Hex8LocalJacobian jacobian;
};
using Quad4FaceLocalDofs = std::array<std::size_t, quad4_face_local_dof_count>;
using Quad4FaceLocalValues = std::array<double, quad4_face_local_dof_count>;
using Quad4FaceLocalResidual = std::array<double, quad4_face_local_dof_count>;
using Quad4FaceLocalJacobian = std::array<double, quad4_face_local_dof_count * quad4_face_local_dof_count>;
using Quad4FaceLocalAdValues = std::array<adlite::Scalar, quad4_face_local_dof_count>;
struct Quad4FaceLocalSystem final {
    Quad4FaceLocalResidual residual;
    Quad4FaceLocalJacobian jacobian;
};
using Hex8Coordinates = std::array<CartesianPoint3, hex8_node_count>;
using Quad4FaceCoordinates = std::array<CartesianPoint3, quad4_face_node_count>;
struct Hex8QuadraturePoint final {
    std::array<double, hex8_node_count> shape;
    std::array<std::array<double, 3>, hex8_node_count> gradient;
    CartesianPoint3 position;
    double weighted_measure;
};
struct Hex8Geometry final {
    std::array<Hex8QuadraturePoint, 8> points;
};
struct Quad4FaceQuadraturePoint final {
    std::array<double, quad4_face_node_count> shape;
    CartesianPoint3 outward_area_vector;
    double weighted_measure;
};
struct Quad4FaceGeometry final {
    std::array<Quad4FaceQuadraturePoint, 4> points;
};
Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates);
Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates);
class Hex8ThermoelasticKernel final {
  public:
    Hex8ThermoelasticKernel(IsotropicThermoelasticMaterial material, double volumetric_heat_source);
    double heat_capacity(double temperature, double x, double y, double z) const;
    void set_volumetric_heat_source(double value) noexcept { _volumetric_heat_source = value; }
    void set_time(double value) noexcept { _time = value; }
    Hex8LocalResidual residual(const Hex8Geometry& geometry, const Hex8LocalValues& state) const;
    Hex8LocalSystem linearize(const Hex8Geometry& geometry, const Hex8LocalValues& state) const;
    Hex8LocalResidual residual(const Hex8Geometry& geometry, const Hex8LocalValues& current_state,
        const Hex8LocalValues& committed_state, double time_step) const;
    Hex8LocalSystem linearize(const Hex8Geometry& geometry, const Hex8LocalValues& current_state,
        const Hex8LocalValues& committed_state, double time_step) const;
    std::array<SymmetricTensor3Values, 8> stress_values(
        const Hex8Geometry& geometry, const Hex8LocalValues& state) const;

  private:
    void residual_ad(const Hex8Geometry& geometry, const Hex8LocalAdValues& state,
        const Hex8LocalValues* committed_state, double time_step, Hex8LocalAdValues& residual) const;
    IsotropicThermoelasticMaterial _material;
    double _volumetric_heat_source, _time;
};
enum class CartesianTractionComponent { x, y, z };
class Quad4FaceBoundaryKernel final {
  public:
    explicit Quad4FaceBoundaryKernel(double pressure);
    Quad4FaceBoundaryKernel(CartesianTractionComponent component, double traction);
    Quad4FaceBoundaryKernel(double heat_transfer_coefficient, double ambient_temperature);
    void set_load(double value) noexcept { _load = value; }
    void set_convection(double coefficient, double ambient) noexcept {
        _load = coefficient;
        _ambient_temperature = ambient;
    }
    Quad4FaceLocalResidual residual(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;
    Quad4FaceLocalSystem linearize(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;

  private:
    enum class Kind { pressure, traction, convection };
    void residual_ad(
        const Quad4FaceGeometry& geometry, const Quad4FaceLocalAdValues& state, Quad4FaceLocalAdValues& residual) const;
    Kind _kind;
    CartesianTractionComponent _component;
    double _load, _ambient_temperature;
};
} // namespace fuelsim
