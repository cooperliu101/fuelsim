#pragma once
#include "c3d20_types.hpp"
#include <vector>

namespace fuelsim::elements {
enum class C3d20Quadrature { full, reduced };

C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request = {});
Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates);
// Checks the current deformation at one supplied reference integration point.
void validate_c3d20t_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state);
// Shared algorithm owned by C3D20T; the reduced model selects its quadrature explicitly.
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request, C3d20Quadrature quadrature);
Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates, C3d20Quadrature quadrature);
std::vector<double> c3d20t_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const CartesianMaterialHistory& history,
    double time,
    C3d20Quadrature quadrature = C3d20Quadrature::full);
} // namespace fuelsim::elements
