#pragma once
#include "c3d20_types.hpp"

namespace fuelsim::elements {
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request = {});
}

namespace fuelsim::elements {
Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates);
}
