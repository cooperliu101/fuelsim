#pragma once
#include "thermal_types.hpp"

namespace fuelsim::elements {
ThermalGeometry make_dc3d20_geometry(const std::array<CartesianPoint3, 20>& coordinates);
ThermalResult evaluate_dc3d20(const ThermalInput& input, bool jacobian = false);
} // namespace fuelsim::elements
