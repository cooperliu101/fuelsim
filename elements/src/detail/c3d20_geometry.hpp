#pragma once
#include "c3d20_types.hpp"

namespace fuelsim {
Hex20Geometry make_hex20_geometry(const Hex20Coordinates& coordinates, std::size_t order);
} // namespace fuelsim
