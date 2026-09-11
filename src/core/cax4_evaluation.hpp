#pragma once
#include "cax4_types.hpp"
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "core/element_region_data.hpp"
#include "material_types.hpp"
#include <stdexcept>

namespace fuelsim {
inline elements::Cax4Result evaluate_cax4(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory* history,
    double time_step,
    bool thermal_time,
    elements::ElementRequest request) {
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
        return elements::evaluate_cax4t(input, request);
    case RzElementFormulation::cax4rt:
        return elements::evaluate_cax4rt(input, request);
    default:
        throw std::invalid_argument("Four-node axisymmetric region requires CAX4T or CAX4RT");
    }
}

inline Cax4LocalResidual compute_cax4_thermoelastic(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    Cax4LocalJacobian* jacobian = nullptr) {
    const auto result =
        evaluate_cax4(data, geometry, state, {}, nullptr, 0, false, {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline std::array<AxisymmetricStressValues, 4> compute_cax4_thermoelastic_stress(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state) {
    return evaluate_cax4(data, geometry, state, {}, nullptr, 0, false, {false, false, false, true}).stress;
}

inline Cax4LocalResidual compute_cax4_transient(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory& history,
    double time_step,
    Cax4LocalJacobian* jacobian = nullptr,
    bool thermal_time = true) {
    const auto result = evaluate_cax4(data,
        geometry,
        state,
        committed,
        &history,
        time_step,
        thermal_time,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline Quad4MaterialHistory compute_cax4_transient_update(const AxisymmetricRegionData& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory& history,
    double time_step) {
    return evaluate_cax4(data, geometry, state, committed, &history, time_step, false, {false, false, true, false})
        .history;
}
} // namespace fuelsim
