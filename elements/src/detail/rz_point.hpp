#pragma once
#include "cax4_types.hpp"

namespace fuelsim::quad4_rz_detail {
inline adlite::Scalar interpolate(const std::array<double, quad4_node_count>& coefficients,
    const Cax4LocalAdValues& state,
    std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node)
        result += coefficients[node] * state[offset + node];
    return result;
}

inline double interpolate(const std::array<double, quad4_node_count>& coefficients,
    const Cax4LocalValues& state,
    std::size_t offset) {
    double result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node)
        result += coefficients[node] * state[offset + node];
    return result;
}

} // namespace fuelsim::quad4_rz_detail
