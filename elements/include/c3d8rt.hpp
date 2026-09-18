#pragma once
#include "c3d8_types.hpp"
#include <vector>

namespace fuelsim::elements {

C3d8Result evaluate_c3d8rt(const C3d8Input& input, ElementRequest request = {});
double c3d8rt_hourglass_energy(const C3d8Input& input);
Hex8Geometry make_c3d8rt_geometry(const Hex8Coordinates& coordinates);
// Checks the current deformation at one supplied reference integration point.
void validate_c3d8rt_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
std::vector<double> c3d8rt_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const CartesianMaterialHistory& history,
    double time,
    StrainFormulation formulation);
} // namespace fuelsim::elements
