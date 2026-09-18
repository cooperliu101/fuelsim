#include "dcax4.hpp"
#include "thermal_common.hpp"
#include <stdexcept>

namespace fuelsim::elements {
ThermalGeometry make_dcax4_geometry(const std::array<CartesianPoint3, 4>& coordinates) {
    return thermal_detail::volume_geometry({coordinates.begin(), coordinates.end()}, ThermalElement::dcax4);
}

ThermalResult evaluate_dcax4(const ThermalInput& input, bool jacobian) {
    if (input.geometry.element != ThermalElement::dcax4 || input.geometry.node_count != 4)
        throw std::invalid_argument("DCAX4 temperature layout mismatch");
    return thermal_detail::integrate(input, jacobian);
}
} // namespace fuelsim::elements
