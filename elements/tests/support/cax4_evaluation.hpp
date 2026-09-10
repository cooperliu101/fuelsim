#pragma once
#include "support/element_test_data.hpp"
// Test-only model selection; the standalone library has no model dispatcher.
#include "axisymmetric_types.hpp"
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include <stdexcept>

namespace fuelsim {
inline elements::Cax4Result evaluate_cax4(const AxisymmetricTestData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory* history,
    double time_step,
    bool jacobian,
    bool thermal_time) {
    const elements::Cax4Input input{data.material,
        geometry,
        state,
        committed,
        history,
        time_step,
        data.time,
        data.volumetric_heat_source,
        data.strain_formulation,
        thermal_time,
        data.initial_temperature};
    switch (data.element_formulation) {
    case RzElementFormulation::cax4t:
        return elements::evaluate_cax4t(input, {true, jacobian, true, false});
    case RzElementFormulation::cax4rt:
        return elements::evaluate_cax4rt(input, {true, jacobian, true, false});
    default:
        throw std::invalid_argument("Four-node axisymmetric region requires CAX4T or CAX4RT");
    }
}

inline Cax4LocalResidual compute_cax4_thermoelastic(const AxisymmetricTestData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    Cax4LocalJacobian* jacobian = nullptr) {
    const auto result = evaluate_cax4(data, geometry, state, {}, nullptr, 0, jacobian != nullptr, false);
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline std::array<AxisymmetricStressValues, 4> compute_cax4_thermoelastic_stress(const AxisymmetricTestData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state) {
    const auto result = evaluate_cax4(data, geometry, state, {}, nullptr, 0, false, false);
    std::array<AxisymmetricStressValues, 4> stress{};
    for (std::size_t q = 0; q < stress.size(); ++q)
        stress[q] = result.history[q].stress;
    return stress;
}

inline Cax4LocalResidual compute_cax4_transient(const AxisymmetricTestData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory& history,
    double time_step,
    Cax4LocalJacobian* jacobian = nullptr,
    bool thermal_time = true) {
    const auto result =
        evaluate_cax4(data, geometry, state, committed, &history, time_step, jacobian != nullptr, thermal_time);
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline Quad4MaterialHistory compute_cax4_transient_update(const AxisymmetricTestData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory& history,
    double time_step) {
    return evaluate_cax4(data, geometry, state, committed, &history, time_step, false, false).history;
}
} // namespace fuelsim

#include "support/rz_linearization.hpp"

namespace fuelsim::rz {
inline LocalLinearization linearize_cax4_thermoelastic(const AxisymmetricTestData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_cax4_thermoelastic(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_cax4_transient(const AxisymmetricTestData& data,
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
