#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::elements {

C3d8Result evaluate_c3d8rt(const C3d8Input& input, ElementRequest request = {});
double c3d8rt_hourglass_energy(const C3d8Input& input);
Hex8Geometry make_c3d8rt_geometry(const Hex8Coordinates& coordinates);
// Checks the current deformation at one supplied reference integration point.
void validate_c3d8rt_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
} // namespace fuelsim::elements
