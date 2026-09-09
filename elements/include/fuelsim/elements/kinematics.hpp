#pragma once

namespace fuelsim {
enum class StrainFormulation {
    small,
    finite,
};

enum class Hex20ElementFormulation { c3d20t, c3d20rt };

enum class Hex8ElementFormulation {
    c3d8t,
    c3d8rt,
};
enum class RzElementFormulation { quad4, cax4t, cax4rt, cax8t, cax8rt };
} // namespace fuelsim
