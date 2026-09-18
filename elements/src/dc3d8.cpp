#include "dc3d8.hpp"
#include "thermal_common.hpp"
#include <stdexcept>

namespace fuelsim::elements {
ThermalGeometry make_dc3d8_geometry(const std::array<CartesianPoint3, 8>& coordinates) {
    return thermal_detail::volume_geometry({coordinates.begin(), coordinates.end()}, ThermalElement::dc3d8);
}

ThermalResult evaluate_dc3d8(const ThermalInput& input, bool jacobian) {
    if (input.geometry.element != ThermalElement::dc3d8 || input.geometry.node_count != 8)
        throw std::invalid_argument("DC3D8 temperature layout mismatch");
    return thermal_detail::integrate(input, jacobian);
}
} // namespace fuelsim::elements
