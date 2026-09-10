#pragma once
#include "c3d20_types.hpp"

namespace fuelsim::elements {
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request = {});
Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates);
// Checks the current deformation at one supplied reference integration point.
void validate_c3d20t_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state);
} // namespace fuelsim::elements
