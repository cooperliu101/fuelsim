#pragma once
#include "ad_local_system.hpp"
#include "cax4_types.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim::quad4_rz_detail {
inline Cax4LocalAdValues ad_state(const Cax4LocalValues& state, bool active = false) {
    Cax4LocalAdValues result{};
    if (active)
        ad_local_system::make_active(state.data(), state.size(), result.data());
    else
        ad_local_system::make_passive(state.data(), state.size(), result.data());
    return result;
}

inline Cax4LocalResidual
values(const Cax4LocalAdValues& state, const Cax4LocalAdValues& residual, Cax4LocalJacobian* jacobian = nullptr) {
    Cax4LocalResidual result{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), result.data(), jacobian->data());
    return result;
}

inline bool finite_point(const RzPoint& point) {
    return std::isfinite(point.r) && std::isfinite(point.z);
}

inline void validate_line(const std::array<RzPoint, 2>& coordinates, const char* name) {
    for (const RzPoint& point : coordinates) {
        if (!finite_point(point))
            throw std::invalid_argument(std::string(name) + " requires finite coordinates");
        if (!(point.r >= 0.0))
            throw std::invalid_argument(std::string(name) + " requires nonnegative radii");
    }
    const double dr = coordinates[1].r - coordinates[0].r, dz = coordinates[1].z - coordinates[0].z;
    if (!(std::hypot(dr, dz) > 0.0))
        throw std::invalid_argument(std::string(name) + " requires a nonzero line length");
}

} // namespace fuelsim::quad4_rz_detail
