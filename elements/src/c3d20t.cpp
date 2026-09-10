#include "c3d20t.hpp"
#include "c3d_common.hpp"
#include <stdexcept>

namespace fuelsim::elements {
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request) {
    if (input.geometry.mechanical_points.size() != 27)
        throw std::invalid_argument("C3D20T requires its model-specific material quadrature");
    return c3d20_detail::evaluate(input, request);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates) {
    return c3d20_detail::make_hex20_geometry(coordinates, 3);
}
} // namespace fuelsim::elements
