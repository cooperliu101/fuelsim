#pragma once
#include "core/cross_section.hpp"
#include <memory>

namespace fuelsim {
struct SectionWarping final {
    std::vector<double> axial;
    std::vector<double> condensed_force;
    double stiffness = 0.0;
    double relative_residual = 0.0;
};

// Owns a section copy and reuses one sparse factorization for all transverse fields.
class SectionWarpingSolver final {
  public:
    explicit SectionWarpingSolver(CrossSection section);
    ~SectionWarpingSolver();
    SectionWarpingSolver(const SectionWarpingSolver&) = delete;
    SectionWarpingSolver& operator=(const SectionWarpingSolver&) = delete;
    SectionWarping solve(const std::vector<double>& transverse) const;

  private:
    class Impl;
    std::unique_ptr<Impl> _impl;
};

struct SectionTorsionMode final {
    CartesianPoint3 origin;
    std::vector<double> transverse;
    SectionWarping warping;
};

SectionTorsionMode build_section_torsion_mode(const CrossSection& section);
SectionKinematics section_torsion_kinematics(const CrossSection& section,
    const SectionTorsionMode& mode,
    std::size_t point,
    double angle,
    double twist,
    double twist_derivative);
} // namespace fuelsim
