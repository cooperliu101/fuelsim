#pragma once
#include "thermal_types.hpp"

namespace fuelsim::elements {
ThermalGeometry make_dcax4_geometry(const std::array<CartesianPoint3, 4>& coordinates);
ThermalResult evaluate_dcax4(const ThermalInput& input, bool jacobian = false);
} // namespace fuelsim::elements
