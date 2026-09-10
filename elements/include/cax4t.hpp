#pragma once
#include "cax4_types.hpp"
#include "element_types.hpp"

namespace fuelsim::elements {
Cax4Result evaluate_cax4t(const Cax4Input& input, ElementRequest request = {});
}
