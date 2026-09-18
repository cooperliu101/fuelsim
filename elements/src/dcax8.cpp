#include "dcax8.hpp"
#include "thermal_common.hpp"
#include <stdexcept>

namespace fuelsim::elements {
ThermalGeometry make_dcax8_geometry(const std::array<CartesianPoint3, 8>& coordinates) {
    return thermal_detail::volume_geometry({coordinates.begin(), coordinates.end()}, ThermalElement::dcax8);
}

ThermalResult evaluate_dcax8(const ThermalInput& input, bool jacobian) {
    if (input.geometry.element != ThermalElement::dcax8 || input.geometry.node_count != 8)
        throw std::invalid_argument("DCAX8 temperature layout mismatch");
    return thermal_detail::integrate(input, jacobian);
}
} // namespace fuelsim::elements
