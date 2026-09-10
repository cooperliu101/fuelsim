#pragma once
#include "axisymmetric_types.hpp"
#include "cax4_types.hpp"
#include "element_types.hpp"

namespace fuelsim::elements {
Cax4Result evaluate_cax4rt(const Cax4Input& input, ElementRequest request = {});
}

namespace fuelsim::elements {
double
cax4rt_hourglass_energy(const AxisymmetricElementData& data, const Quad4RzGeometry& geometry, const LocalValues& state);
}
