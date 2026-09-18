#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::elements {

C3d8Result evaluate_c3d8t(const C3d8Input& input, ElementRequest request = {});
Hex8Geometry make_c3d8t_geometry(const Hex8Coordinates& coordinates);
std::array<double, 8> c3d8t_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const CartesianMaterialHistory& history,
    double time);
// Checks the current deformation at one supplied reference integration point.
void validate_c3d8t_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
} // namespace fuelsim::elements
