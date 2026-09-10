#pragma once
#include "axisymmetric_geometry.hpp"

namespace fuelsim {
using LocalDofs = std::array<std::size_t, local_dof_count>;

struct AxisymmetricElementData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source = 0.0, time = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    RzElementFormulation element_formulation = RzElementFormulation::cax4t;
    double initial_temperature = 600.0;
};
} // namespace fuelsim
