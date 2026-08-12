#ifndef FUELSIM_HEX8_HPP
#define FUELSIM_HEX8_HPP

#include "fuelsim/hex8_local_system.hpp"
#include "fuelsim/mesh.hpp"

#include <array>

namespace fuelsim {

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

} // namespace fuelsim

#endif
