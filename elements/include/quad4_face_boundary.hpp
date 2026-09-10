#pragma once
#include "boundary_types.hpp"
#include "coordinates.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>

namespace fuelsim {
inline constexpr std::size_t quad4_face_node_count = 4;

inline constexpr std::size_t quad4_face_local_dof_count = 16;

using Quad4FaceLocalDofs = std::array<std::size_t, quad4_face_local_dof_count>;

using Quad4FaceLocalValues = std::array<double, quad4_face_local_dof_count>;

using Quad4FaceLocalResidual = std::array<double, quad4_face_local_dof_count>;

using Quad4FaceLocalJacobian = std::array<double, quad4_face_local_dof_count * quad4_face_local_dof_count>;

using Quad4FaceLocalAdValues = std::array<adlite::Scalar, quad4_face_local_dof_count>;

using Quad4FaceCoordinates = std::array<CartesianPoint3, quad4_face_node_count>;

struct Quad4FaceQuadraturePoint final {
    std::array<double, quad4_face_node_count> shape;
    std::array<double, quad4_face_node_count> derivative_xi, derivative_eta;
    CartesianPoint3 tangent_xi, tangent_eta;
    std::array<double, quad4_face_node_count> normal_derivative_xi, normal_derivative_eta;
    CartesianPoint3 normal_tangent_xi, normal_tangent_eta;
    double weighted_measure;
};

struct Quad4FaceGeometry final {
    std::array<Quad4FaceQuadraturePoint, 4> points;
    std::array<Quad4FaceQuadraturePoint, 4> thermal_points;
};

Quad4FaceQuadraturePoint make_quad4_face_quadrature_point(const Quad4FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight);

Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates);

Quad4FaceLocalResidual compute_quad4_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad4FaceGeometry& geometry,
    const Quad4FaceLocalValues& state,
    Quad4FaceLocalJacobian* jacobian = nullptr);
} // namespace fuelsim
