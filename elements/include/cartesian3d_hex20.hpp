#pragma once
#include "cartesian3d_hex8.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {
inline constexpr std::size_t hex20_temperature_node_count = 8;
inline constexpr std::size_t hex20_displacement_node_count = 20;
inline constexpr std::size_t hex20_local_dof_count = 68;
inline constexpr std::size_t hex20_thermal_quadrature_point_count = 8;
inline constexpr std::size_t hex20_mechanical_quadrature_point_count = 27;
inline constexpr std::size_t quad8_face_temperature_node_count = 4;
inline constexpr std::size_t quad8_face_displacement_node_count = 8;
inline constexpr std::size_t quad8_face_local_dof_count = 28;
inline constexpr std::size_t quad8_face_thermal_quadrature_point_count = 4;
inline constexpr std::size_t quad8_face_mechanical_quadrature_point_count = 9;
using Hex20LocalDofs = std::array<std::size_t, hex20_local_dof_count>;
using Hex20LocalValues = std::array<double, hex20_local_dof_count>;
using Hex20LocalResidual = std::array<double, hex20_local_dof_count>;
using Hex20LocalJacobian = std::array<double, hex20_local_dof_count * hex20_local_dof_count>;
using Hex20LocalAdValues = std::array<adlite::Scalar, hex20_local_dof_count>;
using Hex20Coordinates = std::array<CartesianPoint3, hex20_displacement_node_count>;
using Quad8FaceCoordinates = std::array<CartesianPoint3, quad8_face_displacement_node_count>;
using Quad8FaceLocalDofs = std::array<std::size_t, quad8_face_local_dof_count>;
using Quad8FaceLocalValues = std::array<double, quad8_face_local_dof_count>;
using Quad8FaceLocalResidual = std::array<double, quad8_face_local_dof_count>;
using Quad8FaceLocalJacobian = std::array<double, quad8_face_local_dof_count * quad8_face_local_dof_count>;
using Quad8FaceLocalAdValues = std::array<adlite::Scalar, quad8_face_local_dof_count>;

struct Hex20ThermalQuadraturePoint final {
    std::array<double, hex20_temperature_node_count> temperature_shape;
    std::array<std::array<double, 3>, hex20_temperature_node_count> temperature_gradient;
    CartesianPoint3 position;
    double weighted_measure;
};

struct Hex20MechanicalQuadraturePoint final {
    std::array<double, hex20_displacement_node_count> displacement_shape;
    std::array<std::array<double, 3>, hex20_displacement_node_count> displacement_gradient;
    std::array<double, hex20_temperature_node_count> temperature_shape;
    std::array<std::array<double, 3>, hex20_temperature_node_count> temperature_gradient;
    std::array<std::array<double, 3>, hex20_temperature_node_count> source_displacement_gradient;
    CartesianPoint3 position;
    double weighted_measure;
    double source_weighted_measure;
};

struct Hex20Geometry final {
    std::array<Hex20ThermalQuadraturePoint, hex20_thermal_quadrature_point_count> thermal_points;
    std::vector<Hex20MechanicalQuadraturePoint> mechanical_points;
};

struct Quad8FaceThermalQuadraturePoint final {
    std::array<double, quad8_face_temperature_node_count> temperature_shape;
    std::array<double, quad8_face_displacement_node_count> displacement_shape;
    std::array<double, quad8_face_displacement_node_count> derivative_xi, derivative_eta;
    double quadrature_weight;
    double weighted_measure;
};

struct Quad8FaceMechanicalQuadraturePoint final {
    std::array<double, quad8_face_displacement_node_count> displacement_shape;
    std::array<double, quad8_face_displacement_node_count> derivative_xi, derivative_eta;
    CartesianPoint3 tangent_xi, tangent_eta;
    double quadrature_weight;
};

struct Quad8FaceGeometry final {
    std::array<Quad8FaceThermalQuadraturePoint, quad8_face_thermal_quadrature_point_count> thermal_points;
    std::array<Quad8FaceMechanicalQuadraturePoint, quad8_face_mechanical_quadrature_point_count> mechanical_points;
};

Hex20Geometry make_hex20_geometry(const Hex20Coordinates& coordinates,
    Hex20ElementFormulation formulation = Hex20ElementFormulation::c3d20t);
Quad8FaceGeometry make_quad8_face_geometry(const Quad8FaceCoordinates& coordinates);
Quad8FaceMechanicalQuadraturePoint make_quad8_face_mechanical_point(const Quad8FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight);
void validate_hex20_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state);
Hex20LocalResidual compute_hex20_thermoelastic(const CartesianThermoelasticData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues* committed_state = nullptr,
    double time_step = 0.0,
    Hex20LocalJacobian* jacobian = nullptr);
Hex20LocalResidual compute_hex20_transient(const CartesianThermoelasticData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex20LocalJacobian* jacobian = nullptr,
    bool include_thermal_time_term = true);
CartesianMaterialHistory compute_hex20_transient_update(const CartesianThermoelasticData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step);
std::vector<SymmetricTensor3Values> compute_hex20_stress(const CartesianThermoelasticData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state);
Quad8FaceLocalResidual compute_quad8_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad8FaceGeometry& geometry,
    const Quad8FaceLocalValues& state,
    Quad8FaceLocalJacobian* jacobian = nullptr);
} // namespace fuelsim
