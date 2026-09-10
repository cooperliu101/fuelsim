#pragma once
#include "cax4_types.hpp"
#include "element_types.hpp"

namespace fuelsim::elements {
Cax4Result evaluate_cax4t(const Cax4Input& input, ElementRequest request = {});
Quad4RzGeometry make_cax4t_geometry(const Quad4Coordinates& coordinates);
} // namespace fuelsim::elements
