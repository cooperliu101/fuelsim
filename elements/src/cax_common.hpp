#pragma once
#include "cax4_types.hpp"
#include "cax8_types.hpp"

namespace fuelsim::cax8_detail {
Quad8RzGeometry make_quad8_rz_geometry(const Quad8RzCoordinates& coordinates, std::size_t order);
} // namespace fuelsim::cax8_detail

namespace fuelsim::cax8_detail {
elements::Cax8Result evaluate(const elements::Cax8Input& input, elements::ElementRequest request);
} // namespace fuelsim::cax8_detail
