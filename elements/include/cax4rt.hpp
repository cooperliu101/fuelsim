#pragma once
#include "cax4_types.hpp"

namespace fuelsim::elements {
Cax4Result evaluate_cax4rt(const Cax4Input& input, ElementRequest request = {});
double cax4rt_hourglass_energy(const Cax4Input& input);
} // namespace fuelsim::elements
