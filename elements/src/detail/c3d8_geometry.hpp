#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::element_detail {
constexpr double gauss = 0.577350269189625764509148780501957456;
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {{{{-1.0, -1.0, -1.0}},
    {{1.0, -1.0, -1.0}},
    {{1.0, 1.0, -1.0}},
    {{-1.0, 1.0, -1.0}},
    {{-1.0, -1.0, 1.0}},
    {{1.0, -1.0, 1.0}},
    {{1.0, 1.0, 1.0}},
    {{-1.0, 1.0, 1.0}}}};
constexpr std::array<std::size_t, hex8_node_count> hex8_node_to_gauss = {0, 1, 3, 2, 4, 5, 7, 6};

} // namespace fuelsim::element_detail
