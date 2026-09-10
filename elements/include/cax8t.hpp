#pragma once
#include "cax8_types.hpp"

namespace fuelsim::elements {
Cax8Result evaluate_cax8t(const Cax8Input& input, ElementRequest request = {});
}

namespace fuelsim::elements {
Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates);
}
