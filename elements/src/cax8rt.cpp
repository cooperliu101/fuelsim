#include "cax8rt.hpp"
#include "cax_common.hpp"
#include <stdexcept>

namespace fuelsim::elements {
Cax8Result evaluate_cax8rt(const Cax8Input& input, ElementRequest request) {
    if (input.geometry.point_count != 4)
        throw std::invalid_argument("CAX8RT requires its model-specific material quadrature");
    return cax8_detail::evaluate(input, request);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Quad8RzGeometry make_cax8rt_geometry(const Quad8RzCoordinates& coordinates) {
    return cax8_detail::make_quad8_rz_geometry(coordinates, 2);
}
} // namespace fuelsim::elements
