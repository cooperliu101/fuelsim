#pragma once
#include "coordinates.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim::line2_rz_detail {
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

} // namespace fuelsim::line2_rz_detail
