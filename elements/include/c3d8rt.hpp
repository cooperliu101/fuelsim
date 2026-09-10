#pragma once
#include "c3d8_types.hpp"

namespace fuelsim::elements {
C3d8Result evaluate_c3d8rt(const C3d8Input& input, ElementRequest request = {});
}

namespace fuelsim {
double compute_c3d8_mechanical_hourglass_energy(const CartesianThermoelasticData&,
    const Hex8Geometry&,
    const Hex8LocalValues&);
}
