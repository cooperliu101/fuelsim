#pragma once
#include "core/cross_section.hpp"

namespace fuelsim {
// Complete Exodus section node sets bound the interior modal interval. The
// selected outer intervals use native C3D20T displacement unknowns locally.
struct ModalSolidEnds final {
    std::string lower_interface, upper_interface;
};

struct ModalEndRegion final {
    std::size_t mode_count = 0;
    // Complete axial interface node sets in the input Exodus mesh.
    std::string lower_interface, upper_interface;
};

struct SectionMode final {
    enum class Kind {
        extension,
        bending_x_displacement,
        bending_y_displacement,
        torsion,
        distortion,
        poisson_relaxation,
        shear_x,
        shear_y,
        axial_warping,
        shear_free_distortion,
        transverse_corrector,
        width_enrichment
    };
    Kind kind;
    // Coefficients of q, q', q''; each field-major [ux(:),uy(:),uz(:)].
    std::array<std::vector<double>, 3> coefficient;
    // Eigenvalues of physical vector reflection: -1 or +1, 0 if unavailable.
    std::array<int, 2> reflection_parity{};
};

struct ReducedSectionBasis final {
    std::vector<SectionMode> modes;
    // Actual discrete section/material symmetries, not assumed geometry.
    std::array<std::vector<std::size_t>, 2> reflected_nodes;
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
    std::vector<CartesianPoint3> displacement;
};

struct ModalSectionLinearization final {
    std::size_t mode_count = 0;
    // Row-major in [mode, derivative order 0..3]. Only a tangent, never a
    // replacement for axial material-point strain, stress or history evaluation.
    std::vector<double> tangent;
};

// Fixed geometry and basis only. Material response and history are never cached.
struct ModalSectionKinematics final {
    std::size_t mode_count = 0;
    // Section-point major, then mode; coefficients of q through q'''.
    std::vector<std::array<SectionStrain, 4>> strain;
    std::vector<std::array<std::array<double, 3>, 4>> displacement;
};

ModalSectionKinematics sample_modal_section(const CrossSection& section, const ReducedSectionBasis& basis);

std::array<std::array<double, 6>, 4> modal_beam_shape(double lower, double upper, double z);
SectionKinematics modal_section_point_kinematics(const CrossSection& section,
    const SectionMode& mode,
    const SectionPoint& point,
    const std::array<double, 4>& amplitude);
ModalBeamElement make_modal_beam_element(double lower, double upper);
SectionKinematics modal_section_kinematics(const CrossSection& section,
    const SectionMode& mode,
    std::size_t point,
    const std::array<double, 4>& amplitude);
ModalBeamResponse evaluate_modal_beam(const CrossSection& section,
    const ReducedSectionBasis& basis,
    const ModalBeamElement& element,
    const std::vector<double>& state,
    bool include_jacobian = true,
    const ModalSectionKinematics* kinematics = nullptr);
// Fixed, z-independent linear material only. Integrate the actual section
// material tangents, then contract with axial Hermite derivative products.
ModalSectionLinearization linearize_modal_section(const CrossSection& section, const ReducedSectionBasis& basis);
std::vector<double> modal_beam_linear_stiffness(const ModalSectionLinearization& section,
    const ModalBeamElement& element);
} // namespace fuelsim
