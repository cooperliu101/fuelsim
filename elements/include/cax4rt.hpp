#pragma once
#include "rz_quad4.hpp"

namespace fuelsim::rz {
struct Cax4rtResult final {
    LocalResidual residual{};
    LocalJacobian jacobian{};
    Quad4MaterialHistory history{};
};

Cax4rtResult compute_cax4rt(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state,
    const LocalValues& committed,
    const Quad4MaterialHistory* history,
    double time_step,
    bool jacobian,
    bool thermal_time);
std::array<double, 2> cax4rt_thermal_rates(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state,
    const LocalValues& committed,
    double time_step,
    bool thermal_time);
double cax4rt_hourglass_energy(const Quad4RzData& data, const Quad4RzGeometry& geometry, const LocalValues& state);
} // namespace fuelsim::rz
