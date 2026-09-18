#include "c3d20rt.hpp"
#include "c3d20t.hpp"
#include <stdexcept>

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

namespace fuelsim::elements {
std::vector<double> c3d20rt_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const CartesianMaterialHistory& history,
    double time) {
    return c3d20t_creep_rates(material, geometry, state, history, time, C3d20Quadrature::reduced);
}
} // namespace fuelsim::elements
