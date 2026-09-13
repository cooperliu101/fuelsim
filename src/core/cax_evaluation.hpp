#pragma once
#include "cax4_types.hpp"
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "cax8rt.hpp"
#include "cax8t.hpp"
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

inline elements::Cax8Result compute_cax8(const AxisymmetricRegionData& data,
    const Quad8RzGeometry& geometry,
    const Quad8RzValues& state,
    const Quad8RzValues& committed,
    const Quad8MaterialHistory* history,
    double time_step,
    bool thermal_time,
    elements::ElementRequest request) {
    const elements::Cax8Input input{data.material,
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
    case RzElementFormulation::cax8t:
        return elements::evaluate_cax8t(input, request);
    case RzElementFormulation::cax8rt:
        return elements::evaluate_cax8rt(input, request);
    default:
        throw std::invalid_argument("Eight-node axisymmetric region requires CAX8T or CAX8RT");
    }
}

inline Quad8RzGeometry make_cax8_geometry(const Quad8RzCoordinates& coordinates, RzElementFormulation model) {
    switch (model) {
    case RzElementFormulation::cax8t:
        return elements::make_cax8t_geometry(coordinates);
    case RzElementFormulation::cax8rt:
        return elements::make_cax8rt_geometry(coordinates);
    default:
        throw std::invalid_argument("Eight-node axisymmetric geometry requires CAX8T or CAX8RT");
    }
}

} // namespace fuelsim
