#pragma once
#include "cax4_types.hpp"
#include "element_types.hpp"
#include <vector>

namespace fuelsim::elements {
Cax4Result evaluate_cax4t(const Cax4Input& input, ElementRequest request = {});
Quad4RzGeometry make_cax4t_geometry(const Quad4Coordinates& coordinates);
std::vector<double> cax4t_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Quad4MaterialHistory& history,
    double time);
} // namespace fuelsim::elements
