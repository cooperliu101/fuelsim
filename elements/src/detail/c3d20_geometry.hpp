#pragma once
#include "c3d20_types.hpp"

namespace fuelsim::c3d20_detail {
Hex20Geometry make_hex20_geometry(const Hex20Coordinates& coordinates, std::size_t order);
} // namespace fuelsim::c3d20_detail
