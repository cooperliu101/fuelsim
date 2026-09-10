#pragma once
#include "cax8_types.hpp"

namespace fuelsim {
elements::Cax8Result compute_quad8_rz(const AxisymmetricElementData& data,
    const Quad8RzGeometry& geometry,
    const Quad8RzValues& state,
    const Quad8RzValues& committed,
    const Quad8MaterialHistory* history,
    double time_step,
    bool jacobian,
    bool thermal_time = true);
}

namespace fuelsim {
Quad8RzGeometry make_quad8_rz_geometry(const Quad8RzCoordinates& coordinates, std::size_t order);
}
