#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::elements {
C3d8Diagnostics diagnose_c3d8t(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);
C3d8Result evaluate_c3d8t(const C3d8Input& input, ElementRequest request = {});
} // namespace fuelsim::elements
