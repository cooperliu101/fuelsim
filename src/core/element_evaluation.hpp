#pragma once
#include "c3d20rt.hpp"
#include "c3d20t.hpp"
#include "c3d8rt.hpp"
#include "c3d8t.hpp"
#include "core/element_region_data.hpp"

namespace fuelsim {
inline elements::C3d8Result evaluate_c3d8(const CartesianRegionData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    const CartesianMaterialHistory* history,
    double time_step,
    bool thermal_time,
    elements::ElementRequest request) {
    const elements::C3d8Input input{data.material,
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
    return data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt ? elements::evaluate_c3d8rt(input, request)
                                                                           : elements::evaluate_c3d8t(input, request);
}

inline Hex8LocalValues compute_c3d8_thermoelastic(const CartesianRegionData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed = nullptr,
    double time_step = 0,
    Hex8LocalJacobian* jacobian = nullptr) {
    const auto result = evaluate_c3d8(data,
        geometry,
        state,
        committed ? *committed : Hex8LocalValues{},
        nullptr,
        time_step,
        time_step > 0,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline Hex8LocalValues compute_c3d8_transient(const CartesianRegionData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step,
    Hex8LocalJacobian* jacobian = nullptr,
    bool thermal_time = true) {
    const auto result = evaluate_c3d8(data,
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

inline elements::C3d8Result compute_c3d8_transient_update(const CartesianRegionData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step) {
    return evaluate_c3d8(data, geometry, state, committed, &history, time_step, false, {false, false, true, false});
}

inline std::array<SymmetricTensor3Values, 8>
compute_c3d8_stress(const CartesianRegionData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    return evaluate_c3d8(data, geometry, state, {}, nullptr, 0, false, {false, false, false, true}).stress;
}

inline elements::C3d20Result evaluate_c3d20(const CartesianRegionData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed,
    const CartesianMaterialHistory* history,
    double time_step,
    bool thermal_time,
    elements::ElementRequest request) {
    const elements::C3d20Input input{data.material,
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
    return geometry.mechanical_points.size() == 8 ? elements::evaluate_c3d20rt(input, request)
                                                  : elements::evaluate_c3d20t(input, request);
}

inline Hex20LocalValues compute_c3d20_thermoelastic(const CartesianRegionData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues* committed = nullptr,
    double time_step = 0,
    Hex20LocalJacobian* jacobian = nullptr) {
    const auto result = evaluate_c3d20(data,
        geometry,
        state,
        committed ? *committed : Hex20LocalValues{},
        nullptr,
        time_step,
        time_step > 0,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline Hex20LocalValues compute_c3d20_transient(const CartesianRegionData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step,
    Hex20LocalJacobian* jacobian = nullptr,
    bool thermal_time = true) {
    const auto result = evaluate_c3d20(data,
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

inline CartesianMaterialHistory compute_c3d20_transient_update(const CartesianRegionData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step) {
    return evaluate_c3d20(data, geometry, state, committed, &history, time_step, false, {false, false, true, false})
        .history;
}

inline std::vector<SymmetricTensor3Values>
compute_c3d20_stress(const CartesianRegionData& data, const Hex20Geometry& geometry, const Hex20LocalValues& state) {
    return evaluate_c3d20(data, geometry, state, {}, nullptr, 0, false, {false, false, false, true}).stress;
}
} // namespace fuelsim

namespace fuelsim {
inline Hex20Geometry make_c3d20_geometry(const Hex20Coordinates& coordinates,
    Hex20ElementFormulation model = Hex20ElementFormulation::c3d20t) {
    return model == Hex20ElementFormulation::c3d20rt ? elements::make_c3d20rt_geometry(coordinates)
                                                     : elements::make_c3d20t_geometry(coordinates);
}
} // namespace fuelsim
