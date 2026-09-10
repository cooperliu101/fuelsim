#include "c3d20rt.hpp"
#include "c3d_common.hpp"
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

namespace fuelsim::elements {
void validate_c3d20rt_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    c3d20_detail::validate_hex20_deformation(point, state);
}
} // namespace fuelsim::elements
