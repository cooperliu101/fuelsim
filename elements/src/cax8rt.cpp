#include "cax8rt.hpp"
#include "cax8t.hpp"

namespace fuelsim::elements {
Cax8Result evaluate_cax8rt(const Cax8Input& input, ElementRequest request) {
    return evaluate_cax8t(input, request, Cax8Quadrature::reduced);
}

Quad8RzGeometry make_cax8rt_geometry(const Quad8RzCoordinates& coordinates) {
    return make_cax8t_geometry(coordinates, Cax8Quadrature::reduced);
}
} // namespace fuelsim::elements
