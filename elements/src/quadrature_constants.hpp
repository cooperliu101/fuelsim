#pragma once
#include <array>

namespace fuelsim::quadrature {
inline constexpr double gauss2 = 0.577350269189625764509148780501957456;
inline constexpr double gauss3 = 0.774596669241483377035853079956479922;
inline constexpr std::array<double, 2> gauss2_points = {-gauss2, gauss2};
inline constexpr std::array<double, 2> gauss2_weights = {1.0, 1.0};
inline constexpr std::array<double, 3> gauss3_points = {-gauss3, 0.0, gauss3};
inline constexpr std::array<double, 3> gauss3_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
} // namespace fuelsim::quadrature
