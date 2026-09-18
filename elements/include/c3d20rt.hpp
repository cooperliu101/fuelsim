#pragma once
#include "c3d20_types.hpp"
#include <vector>

namespace fuelsim::elements {
C3d20Result evaluate_c3d20rt(const C3d20Input& input, ElementRequest request = {});
Hex20Geometry make_c3d20rt_geometry(const Hex20Coordinates& coordinates);
// Checks the current deformation at one supplied reference integration point.
void validate_c3d20rt_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state);
std::vector<double> c3d20rt_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const CartesianMaterialHistory& history,
    double time);
} // namespace fuelsim::elements
