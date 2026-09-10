#pragma once
#include "boundary_types.hpp"
#include "coordinates.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>

namespace fuelsim {
inline constexpr std::size_t quad8_face_temperature_node_count = 4;

inline constexpr std::size_t quad8_face_displacement_node_count = 8;

inline constexpr std::size_t quad8_face_local_dof_count = 28;

inline constexpr std::size_t quad8_face_thermal_quadrature_point_count = 4;

inline constexpr std::size_t quad8_face_mechanical_quadrature_point_count = 9;

using Quad8FaceCoordinates = std::array<CartesianPoint3, quad8_face_displacement_node_count>;

using Quad8FaceLocalDofs = std::array<std::size_t, quad8_face_local_dof_count>;

using Quad8FaceLocalValues = std::array<double, quad8_face_local_dof_count>;

using Quad8FaceLocalResidual = std::array<double, quad8_face_local_dof_count>;

using Quad8FaceLocalJacobian = std::array<double, quad8_face_local_dof_count * quad8_face_local_dof_count>;

using Quad8FaceLocalAdValues = std::array<adlite::Scalar, quad8_face_local_dof_count>;

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

Quad8FaceGeometry make_quad8_face_geometry(const Quad8FaceCoordinates& coordinates);

Quad8FaceMechanicalQuadraturePoint make_quad8_face_mechanical_point(const Quad8FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight);

Quad8FaceLocalResidual compute_quad8_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad8FaceGeometry& geometry,
    const Quad8FaceLocalValues& state,
    Quad8FaceLocalJacobian* jacobian = nullptr);
} // namespace fuelsim
