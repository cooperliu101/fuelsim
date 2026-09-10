#pragma once
#include "cax4_types.hpp"

namespace fuelsim::elements {
Cax4Result evaluate_cax4rt(const Cax4Input& input, ElementRequest request = {});
double cax4rt_hourglass_energy(const Cax4Input& input);
Quad4RzGeometry make_cax4rt_geometry(const Quad4Coordinates& coordinates);
} // namespace fuelsim::elements
