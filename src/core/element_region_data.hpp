#pragma once
#include "element_types.hpp"
#include "material.hpp"

namespace fuelsim {
// Host-owned material, loading and model selection for each region.
struct AxisymmetricRegionData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source = 0.0, time = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    RzElementFormulation element_formulation = RzElementFormulation::cax4t;
    double initial_temperature = 600.0;
};

struct CartesianRegionData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
    StrainFormulation strain_formulation = StrainFormulation::small;
    Hex8ElementFormulation hex8_element_formulation = Hex8ElementFormulation::c3d8t;
    double initial_temperature = 0.0;
};
} // namespace fuelsim
