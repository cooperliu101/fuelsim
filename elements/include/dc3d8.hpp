#pragma once
#include "thermal_types.hpp"

namespace fuelsim::elements {
ThermalGeometry make_dc3d8_geometry(const std::array<CartesianPoint3, 8>& coordinates);
ThermalResult evaluate_dc3d8(const ThermalInput& input, bool jacobian = false);
} // namespace fuelsim::elements
