#pragma once
#include "cax4_types.hpp"
#include <vector>

namespace fuelsim::elements {
Cax4Result evaluate_cax4rt(const Cax4Input& input, ElementRequest request = {});
double cax4rt_hourglass_energy(const Cax4Input& input);
Quad4RzGeometry make_cax4rt_geometry(const Quad4Coordinates& coordinates);
std::vector<double> cax4rt_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Quad4MaterialHistory& history,
    double time,
    StrainFormulation formulation);
} // namespace fuelsim::elements
