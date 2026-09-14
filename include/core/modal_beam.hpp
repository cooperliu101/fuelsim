#pragma once
#include "core/cross_section.hpp"

namespace fuelsim {
struct SectionMode final {
    enum class Kind { extension, bending_x_displacement, bending_y_displacement, torsion, distortion };
    Kind kind;
    // Coefficients of q, q', q''; each field-major [ux(:),uy(:),uz(:)].
    std::array<std::vector<double>, 3> coefficient;
};

struct ReducedSectionBasis final {
    std::vector<SectionMode> modes;
};

struct ModalBeamPoint final {
    double z;
    double weight;
    std::array<std::array<double, 6>, 4> shape;
};

struct ModalBeamElement final {
    // Quintic Hermite interpolation, local order per mode [qL,q'L,q''L,qR,q'R,q''R].
    double lower, upper;
    std::vector<ModalBeamPoint> points;
};

struct ModalBeamResponse final {
    std::vector<double> residual, jacobian;
    double energy = 0.0;
    // Axial-point major, then section-point major. No material-state modal reduction.
    std::vector<SectionStrain> strain, stress;
};

ModalBeamElement make_modal_beam_element(double lower, double upper);
SectionKinematics modal_section_kinematics(const CrossSection& section,
    const SectionMode& mode,
    std::size_t point,
    const std::array<double, 4>& amplitude);
ModalBeamResponse evaluate_modal_beam(const CrossSection& section,
    const ReducedSectionBasis& basis,
    const ModalBeamElement& element,
    const std::vector<double>& state,
    bool include_jacobian = true);
} // namespace fuelsim
