#pragma once
#include "thermal_types.hpp"

namespace fuelsim::elements {
ThermalGeometry make_dcax8_geometry(const std::array<CartesianPoint3, 8>& coordinates);
ThermalResult evaluate_dcax8(const ThermalInput& input, bool jacobian = false);
} // namespace fuelsim::elements
