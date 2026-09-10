#pragma once
#include "contact.hpp"
#include "rz_quad4.hpp"

namespace fuelsim::rz {
struct LocalLinearization final {
    LocalResidual residual;
    LocalJacobian jacobian;
};

inline LocalLinearization
linearize_quad4_rz_thermoelastic(const Quad4RzData& data, const Quad4RzGeometry& geometry, const LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_quad4_rz_thermoelastic(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_quad4_rz_transient(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& history,
    double time_step) {
    LocalLinearization result{};
    result.residual =
        compute_quad4_rz_transient(data, geometry, state, committed_state, history, time_step, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_line2_rz_boundary(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_line2_rz_gap_heat(const GapHeatProperties& properties,
    const Line2RzHeatPointGeometry& geometry,
    const LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_line2_rz_gap_heat(properties, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry,
    const LocalValues& state,
    const LocalValues& committed_state,
    const ContactPointHistory& history) {
    LocalLinearization result{};
    result.residual =
        compute_node_to_line_rz_contact(properties, geometry, state, committed_state, history, &result.jacobian);
    return result;
}

} // namespace fuelsim::rz
