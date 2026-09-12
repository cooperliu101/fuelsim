#pragma once
#include <adlite/adlite.hpp>
#include <array>

namespace fuelsim::quad8_face_detail {
struct Quad8ShapeValues final {
    std::array<adlite::Scalar, 8> shape{}, derivative_xi{}, derivative_eta{};
    std::array<adlite::Scalar, 8> second_xi{}, second_xi_eta{}, second_eta{};
};

struct DoubleQuad8ShapeValues final {
    std::array<double, 8> shape{}, derivative_xi{}, derivative_eta{};
    std::array<double, 8> second_xi{}, second_xi_eta{}, second_eta{};
};

void quad8_shape(const adlite::Scalar& xi, const adlite::Scalar& eta, Quad8ShapeValues& result);
void double_quad8_shape(double xi, double eta, DoubleQuad8ShapeValues& result);
} // namespace fuelsim::quad8_face_detail
