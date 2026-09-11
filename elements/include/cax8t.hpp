#pragma once
#include "cax8_types.hpp"

namespace fuelsim::elements {
enum class Cax8Quadrature { full, reduced };

Cax8Result evaluate_cax8t(const Cax8Input& input, ElementRequest request = {});
Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates);
// Shared algorithm owned by CAX8T; the reduced model selects its quadrature explicitly.
Cax8Result evaluate_cax8t(const Cax8Input& input, ElementRequest request, Cax8Quadrature quadrature);
Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates, Cax8Quadrature quadrature);
} // namespace fuelsim::elements
