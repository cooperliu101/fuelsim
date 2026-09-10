#pragma once
#include "element_types.hpp"
#include "material.hpp"
#include <vector>

namespace fuelsim {
struct CartesianThermoelasticData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
    StrainFormulation strain_formulation = StrainFormulation::small;
    Hex8ElementFormulation hex8_element_formulation = Hex8ElementFormulation::c3d8t;
    double initial_temperature = 0.0;
};

using CartesianMaterialHistory = std::vector<CartesianMaterialPointState>;
} // namespace fuelsim
