#pragma once
#include "c3d20t.hpp"
#include "core/cross_section.hpp"
#include <utility>

namespace fuelsim {
// A sparse row of u_s = T a. Indices are global on input to the constructor,
// then local to this transformed element. No penalty or interface spring.
using SolidDisplacementRow = std::vector<std::pair<std::size_t, double>>;

struct ModalSolidElement final {
    std::size_t layer = 0, region = 0;
    Hex20Geometry geometry;
    std::vector<std::size_t> dofs;
    std::array<SolidDisplacementRow, 60> displacement;
};

struct ModalSolidResponse final {
    std::vector<double> residual, jacobian;
    std::vector<SectionStrain> strain, stress;
    std::vector<CartesianPoint3> displacement;
    double energy = 0.0;
};

ModalSolidElement make_modal_solid_element(const Hex20Coordinates& coordinates,
    std::array<SolidDisplacementRow, 60> displacement,
    std::size_t layer,
    std::size_t region);

ModalSolidResponse evaluate_modal_solid(const ModalSolidElement& element,
    const SectionRegion& region,
    const std::vector<double>& state,
    bool jacobian,
    bool output);
} // namespace fuelsim
