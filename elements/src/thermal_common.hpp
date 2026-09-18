#pragma once
#include "thermal_types.hpp"

namespace fuelsim::elements::thermal_detail {
ThermalResult integrate(const ThermalInput& input, bool jacobian);
ThermalGeometry volume_geometry(const std::vector<CartesianPoint3>& coordinates, ThermalElement element);
} // namespace fuelsim::elements::thermal_detail
