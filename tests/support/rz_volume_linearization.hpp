#pragma once
#include "core/cax4_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "support/test_support.hpp"

namespace fuelsim::rz {
inline LocalLinearization linearize_cax4_thermoelastic(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_cax4_thermoelastic(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_cax4_transient(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed_state,
    const Quad4MaterialHistory& history,
    double time_step) {
    LocalLinearization result{};
    result.residual =
        compute_cax4_transient(data, geometry, state, committed_state, history, time_step, &result.jacobian);
    return result;
}

} // namespace fuelsim::rz
