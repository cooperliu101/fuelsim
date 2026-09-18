#include "dc3d20.hpp"
#include "thermal_common.hpp"
#include <stdexcept>

namespace fuelsim::elements {
ThermalGeometry make_dc3d20_geometry(const std::array<CartesianPoint3, 20>& coordinates) {
    return thermal_detail::volume_geometry({coordinates.begin(), coordinates.end()}, ThermalElement::dc3d20);
}

ThermalResult evaluate_dc3d20(const ThermalInput& input, bool jacobian) {
    if (input.geometry.element != ThermalElement::dc3d20 || input.geometry.node_count != 20)
        throw std::invalid_argument("DC3D20 temperature layout mismatch");
    return thermal_detail::integrate(input, jacobian);
}
} // namespace fuelsim::elements
