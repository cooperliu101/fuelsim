#pragma once
#include "core/cross_section.hpp"
#include <memory>

namespace fuelsim {
struct SectionTransverseRelaxation final {
    std::vector<double> transverse;
    double relative_residual = 0.0;
};

// Minimize full 3D elastic energy over transverse displacement for a prescribed
// nodal axial strain field. Reuses the section distortion factorization.
class SectionRelaxationSolver final {
  public:
    explicit SectionRelaxationSolver(CrossSection section);
    ~SectionRelaxationSolver();
    SectionRelaxationSolver(const SectionRelaxationSolver&) = delete;
    SectionRelaxationSolver& operator=(const SectionRelaxationSolver&) = delete;
    SectionTransverseRelaxation solve(const std::vector<double>& axial_strain) const;
    // KD*p = integral(Nperp^T*Cs*source), projected off the global force space.
    // This is the transverse axial-variation term of 3D elasticity, not a vibration mode.
    SectionTransverseRelaxation solve_transverse_corrector(const std::vector<double>& source) const;

  private:
    class Impl;
    std::unique_ptr<Impl> _impl;
};

// Requires an initialized PETSc session. Offline preprocessing is local to each caller.
ClassicSectionBasis build_classic_section_basis(const CrossSection& section);
} // namespace fuelsim
