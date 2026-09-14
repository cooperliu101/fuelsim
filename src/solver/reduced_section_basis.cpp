#include "solver/reduced_section_basis.hpp"
#include "solver/section_distortion.hpp"
#include "solver/section_modes.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
ReducedSectionBasis
build_reduced_section_basis(const CrossSection& section, bool torsion, std::size_t enrichment_modes) {
    const auto classic = build_classic_section_basis(section);
    const auto n = section.nodes().size();
    ReducedSectionBasis result;
    for (std::size_t i = 0; i < 3; ++i) {
        SectionMode mode;
        mode.kind = static_cast<SectionMode::Kind>(i);
        for (auto& coefficient : mode.coefficient)
            coefficient.resize(3 * n);
        for (std::size_t node = 0; node < n; ++node) {
            if (i == 0)
                mode.coefficient[0][2 * n + node] = 1.0;
            else {
                mode.coefficient[0][(i - 1) * n + node] = 1.0;
                mode.coefficient[1][2 * n + node] = i == 1 ? -(section.nodes()[node].x - classic.origin.x)
                                                           : -(section.nodes()[node].y - classic.origin.y);
            }
            const auto order = i == 0 ? 1 : 2;
            mode.coefficient[static_cast<std::size_t>(order)][node] = classic.modes[i].correction[node];
            mode.coefficient[static_cast<std::size_t>(order)][n + node] = classic.modes[i].correction[n + node];
        }
        result.modes.push_back(std::move(mode));
    }
    if (torsion) {
        const auto source = build_section_torsion_mode(section);
        SectionMode mode;
        mode.kind = SectionMode::Kind::torsion;
        for (auto& coefficient : mode.coefficient)
            coefficient.resize(3 * n);
        for (std::size_t i = 0; i < 2 * n; ++i)
            mode.coefficient[0][i] = source.transverse[i];
        for (std::size_t i = 0; i < n; ++i)
            mode.coefficient[1][2 * n + i] = source.warping.axial[i];
        result.modes.push_back(std::move(mode));
    }
    if (enrichment_modes == 0)
        return result;
    std::size_t added = 0;
    // Independent shear rotations retain the exact Euler-Bernoulli modes while
    // allowing uz = -x*w' + x*sx (and the y counterpart). Setting sx=sy=0
    // recovers the original bending kinematics; varying them tests physical
    // section moment balance independently of transverse displacement.
    for (std::size_t c = 0; c < 2 && added < enrichment_modes; ++c) {
        SectionMode mode;
        mode.kind = c == 0 ? SectionMode::Kind::shear_x : SectionMode::Kind::shear_y;
        for (auto& coefficient : mode.coefficient)
            coefficient.resize(3 * n);
        for (std::size_t node = 0; node < n; ++node)
            mode.coefficient[0][2 * n + node] =
                c == 0 ? section.nodes()[node].x - classic.origin.x : section.nodes()[node].y - classic.origin.y;
        result.modes.push_back(std::move(mode));
        ++added;
    }
    if (added == enrichment_modes)
        return result;
    // The classical modes prescribe lateral Poisson contraction through axial
    // derivatives. Independent amplitudes of these same transverse fields are
    // essential at a physical clamp: otherwise its nodal constraints force the
    // classical axial strain/curvature to vanish. They are distortion directions,
    // not additional rigid/global modes, and their warping is solved variationally.
    std::vector<std::vector<double>> independent;
    const auto inner = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double value = 0.0;
        for (const auto& point : section.points())
            for (std::size_t c = 0; c < 2; ++c) {
                double left = 0.0, right = 0.0;
                for (std::size_t i = 0; i < 8; ++i) {
                    left += point.shape[i] * a[c * n + point.nodes[i]];
                    right += point.shape[i] * b[c * n + point.nodes[i]];
                }
                value += point.weight * left * right;
            }
        return value;
    };
    SectionWarpingSolver warping(section);
    const auto append = [&](const std::vector<double>& transverse, SectionMode::Kind kind) {
        const double norm = inner(transverse, transverse);
        if (norm == 0.0)
            return false;
        auto remainder = transverse;
        for (int pass = 0; pass < 2; ++pass)
            for (const auto& previous : independent) {
                const double projection = inner(previous, remainder);
                for (std::size_t i = 0; i < 2 * n; ++i)
                    remainder[i] -= projection * previous[i];
            }
        const double independent_norm = inner(remainder, remainder);
        if (independent_norm <= 1.0e-16 * norm)
            return false;
        for (double& value : remainder)
            value /= std::sqrt(independent_norm);
        independent.push_back(std::move(remainder));
        // Preserve the original shape (and its pencil eigenvalue if applicable).
        // Orthogonalization above is solely a rank check; it does not silently
        // rotate a finite eigenmode away from its defining generalized equation.
        std::vector<double> normalized = transverse;
        for (double& value : normalized)
            value /= std::sqrt(norm);
        const auto axial = warping.solve(normalized);
        SectionMode mode;
        mode.kind = kind;
        for (auto& coefficient : mode.coefficient)
            coefficient.resize(3 * n);
        for (std::size_t i = 0; i < 2 * n; ++i)
            mode.coefficient[0][i] = normalized[i];
        for (std::size_t i = 0; i < n; ++i)
            mode.coefficient[1][2 * n + i] = axial.axial[i];
        result.modes.push_back(std::move(mode));
        ++added;
        return true;
    };
    for (const auto& mode : classic.modes)
        if (added < enrichment_modes)
            (void)append(mode.correction, SectionMode::Kind::poisson_relaxation);
    if (added < enrichment_modes) {
        const auto remaining = enrichment_modes - added;
        const auto spectrum = build_section_distortion_modes(section, remaining);
        for (const auto& mode : spectrum.modes)
            if (!append(mode.transverse, SectionMode::Kind::distortion))
                throw std::runtime_error("Selected section enrichment is linearly dependent");
    }
    return result;
}
} // namespace fuelsim
