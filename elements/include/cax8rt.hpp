#pragma once
#include "cax8_types.hpp"
#include <vector>

namespace fuelsim::elements {
Cax8Result evaluate_cax8rt(const Cax8Input& input, ElementRequest request = {});
Quad8RzGeometry make_cax8rt_geometry(const Quad8RzCoordinates& coordinates);
std::vector<double> cax8rt_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Quad8RzGeometry& geometry,
    const Quad8RzValues& state,
    const Quad8MaterialHistory& history,
    double time);
} // namespace fuelsim::elements
