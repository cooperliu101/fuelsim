#pragma once
#include "solver/section_warping.hpp"

namespace fuelsim {
struct SectionDistortionMode final {
    std::vector<double> transverse;
    SectionWarping warping;
    double eigenvalue = 0.0;
    double distortion_energy = 0.0;
    double section_norm = 0.0;
    double shear_norm = 0.0;
    double eigen_relative_residual = 0.0;
    double classic_projection = 0.0;
    double orthogonality_error = 0.0;
};

struct SectionDistortionSpectrum final {
    std::vector<SectionDistortionMode> modes;
    std::vector<SectionDistortionMode> shear_free_modes;
    std::size_t shear_kernel_dimension = 0;
};

// Offline projected dense solve, deliberately bounded to 512 transverse DOFs.
// Uses sparse warping solves and a private symmetric section eigensolver, no new dependency.
SectionDistortionSpectrum
build_section_distortion_modes(const CrossSection& section, std::size_t count, std::size_t shear_free_count = 0);
} // namespace fuelsim
