#ifndef FUELSIM_QUAD4_RZ_HPP
#define FUELSIM_QUAD4_RZ_HPP
#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"
#include <array>
#include <cstddef>
namespace fuelsim {
constexpr std::size_t quad4_node_count = 4;
constexpr std::size_t quad4_local_dof_count = local_dof_count;
constexpr std::size_t quad4_jacobian_size = local_jacobian_size;
using Quad4Coordinates = std::array<RzPoint, quad4_node_count>;
struct RzQuadraturePoint final {
    std::array<double, quad4_node_count> shape;
    std::array<double, quad4_node_count> gradient_r;
    std::array<double, quad4_node_count> gradient_z;
    double radius;
    double axial_coordinate;
    double weighted_measure;
};
struct Quad4RzGeometry final {
    std::array<RzQuadraturePoint, 4> points;
};
Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates);
} // namespace fuelsim
#endif
