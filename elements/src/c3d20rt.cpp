#include "c3d20rt.hpp"
#include "c3d20t.hpp"

namespace fuelsim::elements {
C3d20Result evaluate_c3d20rt(const C3d20Input& input, ElementRequest request) {
    return evaluate_c3d20t(input, request, C3d20Quadrature::reduced);
}

Hex20Geometry make_c3d20rt_geometry(const Hex20Coordinates& coordinates) {
    return make_c3d20t_geometry(coordinates, C3d20Quadrature::reduced);
}

void validate_c3d20rt_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    validate_c3d20t_deformation(point, state);
}
} // namespace fuelsim::elements
