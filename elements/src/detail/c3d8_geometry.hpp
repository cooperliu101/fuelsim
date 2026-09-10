#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::c3d8_detail {
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {{{{-1.0, -1.0, -1.0}},
    {{1.0, -1.0, -1.0}},
    {{1.0, 1.0, -1.0}},
    {{-1.0, 1.0, -1.0}},
    {{-1.0, -1.0, 1.0}},
    {{1.0, -1.0, 1.0}},
    {{1.0, 1.0, 1.0}},
    {{-1.0, 1.0, 1.0}}}};
// This permutation is its own inverse: node order to Gauss order and back.
constexpr std::array<std::size_t, hex8_node_count> hex8_node_gauss_permutation = {0, 1, 3, 2, 4, 5, 7, 6};

} // namespace fuelsim::c3d8_detail
