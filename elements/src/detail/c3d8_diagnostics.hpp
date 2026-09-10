#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::c3d8_detail {
elements::C3d8Diagnostics diagnose_hex8(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    StrainFormulation formulation,
    bool reduced);
} // namespace fuelsim::c3d8_detail
