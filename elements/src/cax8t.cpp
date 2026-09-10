#include "cax8t.hpp"
#include "detail/cax8_assembly.hpp"
#include "detail/cax8_geometry.hpp"
#include <stdexcept>

namespace fuelsim::elements {
Cax8Result evaluate_cax8t(const Cax8Input& input, ElementRequest request) {
    if (input.geometry.point_count != 9)
        throw std::invalid_argument("CAX8T requires its model-specific material quadrature");
    return cax8_detail::evaluate(input, request);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates) {
    return cax8_detail::make_quad8_rz_geometry(coordinates, 3);
}
} // namespace fuelsim::elements
