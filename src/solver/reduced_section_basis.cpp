#include "solver/reduced_section_basis.hpp"
#include "solver/section_distortion.hpp"
#include "solver/section_modes.hpp"

namespace fuelsim {
ReducedSectionBasis
build_reduced_section_basis(const CrossSection& section, bool torsion, std::size_t distortion_modes) {
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
    if (distortion_modes > 0) {
        const auto spectrum = build_section_distortion_modes(section, distortion_modes);
        for (const auto& source : spectrum.modes) {
            SectionMode mode;
            mode.kind = SectionMode::Kind::distortion;
            for (auto& coefficient : mode.coefficient)
                coefficient.resize(3 * n);
            for (std::size_t i = 0; i < 2 * n; ++i)
                mode.coefficient[0][i] = source.transverse[i];
            for (std::size_t i = 0; i < n; ++i)
                mode.coefficient[1][2 * n + i] = source.warping.axial[i];
            result.modes.push_back(std::move(mode));
        }
    }
    return result;
}
} // namespace fuelsim
