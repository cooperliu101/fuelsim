#pragma once

namespace fuelsim {
enum class StrainFormulation {
    small,
    finite,
};

enum class Hex8ElementFormulation {
    c3d8t,
    c3d8rt,
};
enum class RzElementFormulation { quad4, cax4t };
} // namespace fuelsim
