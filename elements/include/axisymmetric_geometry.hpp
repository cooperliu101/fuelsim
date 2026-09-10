#pragma once
#include "coordinates.hpp"
#include "element_types.hpp"
#include "material.hpp"
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

using Quad4MaterialHistory = std::array<MaterialPointState, 4>;

} // namespace fuelsim
