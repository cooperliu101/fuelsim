#include "c3d20rt.hpp"
#include "detail/c3d20_assembly.hpp"
#include "detail/c3d20_geometry.hpp"
#include <stdexcept>

namespace fuelsim::elements {
C3d20Result evaluate_c3d20rt(const C3d20Input& input, ElementRequest request) {
    if (input.geometry.mechanical_points.size() != 8)
        throw std::invalid_argument("C3D20RT requires its model-specific material quadrature");
    return c3d20_detail::evaluate(input, request);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Hex20Geometry make_c3d20rt_geometry(const Hex20Coordinates& coordinates) {
    return c3d20_detail::make_hex20_geometry(coordinates, 2);
}
} // namespace fuelsim::elements
