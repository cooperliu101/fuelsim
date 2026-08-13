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
struct Hex8ThermoelasticData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
};
double compute_hex8_heat_capacity(const Hex8ThermoelasticData& data, double temperature, double x, double y, double z);
Hex8LocalResidual compute_hex8_residual(
    const Hex8ThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);
Hex8LocalSystem compute_hex8_system(
    const Hex8ThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);
Hex8LocalResidual compute_hex8_transient_residual(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& current_state, const Hex8LocalValues& committed_state, double time_step);
Hex8LocalSystem compute_hex8_transient_system(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& current_state, const Hex8LocalValues& committed_state, double time_step);
std::array<SymmetricTensor3Values, 8> compute_hex8_stress(
    const Hex8ThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);
enum class CartesianTractionComponent { x, y, z };
enum class Quad4FaceBoundaryKind { pressure, traction, convection };
struct Quad4FaceBoundaryData final {
    Quad4FaceBoundaryKind kind;
    CartesianTractionComponent component;
    double load, ambient_temperature;
};
Quad4FaceBoundaryData make_quad4_face_pressure_data(double pressure);
Quad4FaceBoundaryData make_quad4_face_traction_data(CartesianTractionComponent component, double traction);
Quad4FaceBoundaryData make_quad4_face_convection_data(double heat_transfer_coefficient, double ambient_temperature);
Quad4FaceLocalResidual compute_quad4_face_boundary_residual(
    const Quad4FaceBoundaryData& data, const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state);
Quad4FaceLocalSystem compute_quad4_face_boundary_system(
    const Quad4FaceBoundaryData& data, const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state);
} // namespace fuelsim
