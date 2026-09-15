#include "solver/reduced_section_basis.hpp"
#include "solver/section_distortion.hpp"
#include "solver/section_modes.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim {
namespace {
using ReflectionMaps = std::array<std::vector<std::size_t>, 2>;

ReflectionMaps section_reflections(const CrossSection& section) {
    // Recognize actual discrete symmetries, including quadrature and material,
    // rather than assuming a rectangular mesh or symmetric material layering.
    // A missing symmetry simply leaves the original candidate space unchanged.
    ReflectionMaps result;
    const auto n = section.nodes().size();
    const auto center = section.centroid();
    for (std::size_t axis = 0; axis < 2; ++axis) {
        const double origin = axis == 0 ? center.x : center.y;
        double scale = 0.0;
        for (const auto& node : section.nodes())
            scale = std::max({scale, std::abs(node.x), std::abs(node.y)});
        const double tolerance = 64 * std::numeric_limits<double>::epsilon() * scale;
        auto& map = result[axis];
        map.resize(n, n);
        bool valid = true;
        for (std::size_t i = 0; i < n && valid; ++i) {
            auto reflected = section.nodes()[i];
            (axis == 0 ? reflected.x : reflected.y) = 2 * origin - (axis == 0 ? reflected.x : reflected.y);
            for (std::size_t j = 0; j < n; ++j)
                if (std::abs(reflected.x - section.nodes()[j].x) <= tolerance
                    && std::abs(reflected.y - section.nodes()[j].y) <= tolerance) {
                    if (map[i] != n) {
                        valid = false;
                        break;
                    }
                    map[i] = j;
                }
            valid = valid && map[i] != n;
        }
        for (std::size_t i = 0; i < n && valid; ++i)
            valid = map[map[i]] == i;
        for (const auto& point : section.points()) {
            if (!valid)
                break;
            auto position = point.position;
            (axis == 0 ? position.x : position.y) = 2 * origin - (axis == 0 ? position.x : position.y);
            const auto found = std::find_if(section.points().begin(), section.points().end(), [&](const auto& other) {
                return other.region == point.region && std::abs(other.position.x - position.x) <= tolerance
                       && std::abs(other.position.y - position.y) <= tolerance;
            });
            if (found == section.points().end()) {
                valid = false;
                break;
            }
            const double roundoff = 256 * std::numeric_limits<double>::epsilon();
            valid = std::abs(point.weight - found->weight) <= roundoff * point.weight;
            for (std::size_t i = 0; i < 8 && valid; ++i) {
                const auto node = std::find(found->nodes.begin(), found->nodes.end(), map[point.nodes[i]]);
                if (node == found->nodes.end()) {
                    valid = false;
                    break;
                }
                const auto j = static_cast<std::size_t>(node - found->nodes.begin());
                valid = std::abs(point.shape[i] - found->shape[j]) <= roundoff;
                for (std::size_t c = 0; c < 2; ++c) {
                    const double expected = (c == axis ? -1.0 : 1.0) * point.gradient[i][c];
                    const double gradient_scale = std::hypot(point.gradient[i][0], point.gradient[i][1]);
                    valid = valid && std::abs(expected - found->gradient[j][c]) <= roundoff * gradient_scale;
                }
            }
            const auto& region = section.regions()[point.region];
            const auto left =
                region.material.active_properties(region.temperature, material_context(0.0, point.position));
            const auto right =
                region.material.active_properties(region.temperature, material_context(0.0, found->position));
            valid = valid && left.lame_lambda.value() == right.lame_lambda.value()
                    && left.shear_modulus.value() == right.shear_modulus.value();
        }
        if (!valid)
            map.clear();
    }
    return result;
}

struct SymmetryComponent final {
    std::vector<double> field;
    std::array<int, 2> parity{};
};

void project_symmetry(std::vector<double>& field, const ReflectionMaps& maps, const std::array<int, 2>& parity) {
    for (std::size_t axis = 0; axis < 2; ++axis) {
        if (maps[axis].empty())
            continue;
        const auto n = maps[axis].size();
        const auto original = field;
        for (std::size_t c = 0; c < field.size() / n; ++c) {
            // A transverse vector changes its normal component's sign under
            // reflection. A scalar axial displacement does not.
            const double sign = field.size() == 2 * n && c == axis ? -1.0 : 1.0;
            for (std::size_t i = 0; i < n; ++i)
                field[c * n + i] = 0.5 * (original[c * n + i] + parity[axis] * sign * original[c * n + maps[axis][i]]);
        }
    }
}

std::vector<SymmetryComponent> symmetry_components(const std::vector<double>& field, const ReflectionMaps& maps) {
    std::vector<SymmetryComponent> result;
    for (int x : {1, -1})
        for (int y : {1, -1}) {
            if ((maps[0].empty() && x == -1) || (maps[1].empty() && y == -1))
                continue;
            SymmetryComponent part{field, {x, y}};
            project_symmetry(part.field, maps, part.parity);
            result.push_back(std::move(part));
        }
    return result;
}

void retain_mode(ReducedSectionBasis& basis, SectionMode mode) {
    // Auxiliary solves have roundoff in otherwise exact symmetry sectors.
    // Apply the known vector reflection to every kinematic coefficient, so
    // stiffness and physical constraints use the same discrete symmetry.
    for (std::size_t axis = 0; axis < 2; ++axis) {
        const auto& map = basis.reflected_nodes[axis];
        if (map.empty()) {
            mode.reflection_parity[axis] = 0;
            continue;
        }
        const auto n = map.size();
        for (auto& field : mode.coefficient) {
            const auto original = field;
            for (std::size_t c = 0; c < 3; ++c) {
                const double sign = mode.reflection_parity[axis] * (c == axis ? -1.0 : 1.0);
                for (std::size_t i = 0; i < n; ++i)
                    field[c * n + i] = 0.5 * (original[c * n + i] + sign * original[c * n + map[i]]);
            }
        }
    }
    basis.modes.push_back(std::move(mode));
}
} // namespace

ReducedSectionBasis
build_reduced_section_basis(const CrossSection& section, bool torsion, std::size_t enrichment_modes) {
    const auto n = section.nodes().size();
    const std::size_t classic_count = torsion ? 4 : 3;
    if (enrichment_modes > 3 * n - classic_count)
        throw std::invalid_argument("Requested modes exceed the independent three-component section space");
    const auto classic = build_classic_section_basis(section);
    ReducedSectionBasis result;
    result.reflected_nodes = section_reflections(section);
    for (std::size_t i = 0; i < 3; ++i) {
        SectionMode mode;
        mode.kind = static_cast<SectionMode::Kind>(i);
        mode.reflection_parity = {i == 1 ? -1 : 1, i == 2 ? -1 : 1};
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
        retain_mode(result, std::move(mode));
    }
    if (torsion) {
        const auto source = build_section_torsion_mode(section);
        SectionMode mode;
        mode.kind = SectionMode::Kind::torsion;
        mode.reflection_parity = {-1, -1};
        for (auto& coefficient : mode.coefficient)
            coefficient.resize(3 * n);
        for (std::size_t i = 0; i < 2 * n; ++i)
            mode.coefficient[0][i] = source.transverse[i];
        for (std::size_t i = 0; i < n; ++i)
            mode.coefficient[1][2 * n + i] = source.warping.axial[i];
        retain_mode(result, std::move(mode));
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
        mode.reflection_parity = {c == 0 ? -1 : 1, c == 1 ? -1 : 1};
        for (auto& coefficient : mode.coefficient)
            coefficient.resize(3 * n);
        for (std::size_t node = 0; node < n; ++node)
            mode.coefficient[0][2 * n + node] =
                c == 0 ? section.nodes()[node].x - classic.origin.x : section.nodes()[node].y - classic.origin.y;
        retain_mode(result, std::move(mode));
        ++added;
    }
    if (added == enrichment_modes)
        return result;
    const auto& reflections = result.reflected_nodes;
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
    for (std::size_t global = 0; global < 3; ++global) {
        std::vector<double> field(2 * n);
        for (std::size_t i = 0; i < n; ++i) {
            if (global == 0)
                field[i] = 1.0;
            else if (global == 1)
                field[n + i] = 1.0;
            else {
                field[i] = -(section.nodes()[i].y - section.centroid().y);
                field[n + i] = section.nodes()[i].x - section.centroid().x;
            }
        }
        const double norm = std::sqrt(inner(field, field));
        for (double& value : field)
            value /= norm;
        independent.push_back(std::move(field));
    }
    SectionWarpingSolver warping(section);
    const auto append_component =
        [&](const std::vector<double>& transverse, const std::array<int, 2>& parity, SectionMode::Kind kind) {
            const double norm = inner(transverse, transverse);
            if (norm == 0.0)
                return false;
            auto remainder = transverse;
            for (int pass = 0; pass < 2; ++pass) {
                for (const auto& previous : independent) {
                    const double projection = inner(previous, remainder);
                    for (std::size_t i = 0; i < 2 * n; ++i)
                        remainder[i] -= projection * previous[i];
                }
                project_symmetry(remainder, reflections, parity);
            }
            const double independent_norm = inner(remainder, remainder);
            if (independent_norm <= 1.0e-16 * norm)
                return false;
            for (double& value : remainder)
                value /= std::sqrt(independent_norm);
            independent.push_back(remainder);
            // Stable coordinates in the nested enrichment span. The offline
            // eigensolver retains the original eigenvectors and their diagnostics;
            // these production Ritz directions are explicitly orthogonal mixtures,
            // not individual eigenvectors. Warping is linear in the transverse
            // field, so recomputing it preserves the complete coupled span.
            const auto& normalized = remainder;
            const auto axial = warping.solve(normalized);
            SectionMode mode;
            mode.kind = kind;
            mode.reflection_parity = parity;
            for (auto& coefficient : mode.coefficient)
                coefficient.resize(3 * n);
            for (std::size_t i = 0; i < 2 * n; ++i)
                mode.coefficient[0][i] = normalized[i];
            for (std::size_t i = 0; i < n; ++i)
                mode.coefficient[1][2 * n + i] = axial.axial[i];
            retain_mode(result, std::move(mode));
            ++added;
            return true;
        };
    const auto append = [&](const std::vector<double>& transverse, SectionMode::Kind kind) {
        auto parts = symmetry_components(transverse, reflections);
        std::stable_sort(parts.begin(), parts.end(), [&](const auto& a, const auto& b) {
            return inner(a.field, a.field) > inner(b.field, b.field);
        });
        const double initial = inner(transverse, transverse);
        bool accepted = false;
        for (const auto& part : parts) {
            if (added == enrichment_modes)
                break;
            if (inner(part.field, part.field) > 1e-16 * initial)
                accepted = append_component(part.field, part.parity, kind) || accepted;
        }
        return accepted;
    };
    for (const auto& mode : classic.modes)
        if (added < enrichment_modes)
            (void)append(mode.correction, SectionMode::Kind::poisson_relaxation);
    if (added == enrichment_modes)
        return result;
    // Minimizing transverse shear fixes omega for a slowly varying transverse
    // field, but does not make uz=omega*q' a complete axial displacement space.
    // Give the computed warping shapes independent amplitudes as well. Their
    // energy includes G*|grad(omega)|^2*r^2 and Czz*omega^2*r'^2; the full 3D
    // material call retains their coupling to every transverse strain. No
    // prescribed polynomial warping or load-dependent fitted shape is used.
    const auto axial_inner = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double value = 0.0;
        for (const auto& point : section.points()) {
            double left = 0.0, right = 0.0;
            for (std::size_t i = 0; i < 8; ++i) {
                left += point.shape[i] * a[point.nodes[i]];
                right += point.shape[i] * b[point.nodes[i]];
            }
            value += point.weight * left * right;
        }
        return value;
    };
    std::vector<std::vector<double>> axial_basis;
    const auto axial_direction = [&](std::vector<double> field, const std::array<int, 2>& parity) {
        const double initial = axial_inner(field, field);
        for (int pass = 0; pass < 2; ++pass) {
            for (const auto& previous : axial_basis) {
                const double projection = axial_inner(previous, field);
                for (std::size_t i = 0; i < n; ++i)
                    field[i] -= projection * previous[i];
            }
            project_symmetry(field, reflections, parity);
        }
        const double norm = axial_inner(field, field);
        if (norm <= 1.0e-16 * initial || norm == 0.0)
            return false;
        for (double& value : field)
            value /= std::sqrt(norm);
        const auto pivot =
            std::max_element(field.begin(), field.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
        if (*pivot < 0.0)
            for (double& value : field)
                value = -value;
        axial_basis.push_back(std::move(field));
        return true;
    };
    for (std::size_t global = 0; global < 3; ++global) {
        std::vector<double> field(n, 1.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (global == 1)
                field[i] = section.nodes()[i].x - classic.origin.x;
            if (global == 2)
                field[i] = section.nodes()[i].y - classic.origin.y;
        }
        const std::array<int, 2> parity{global == 1 ? -1 : 1, global == 2 ? -1 : 1};
        if (!axial_direction(std::move(field), parity))
            throw std::runtime_error("Classical axial section space is linearly dependent");
    }
    const auto append_axial = [&](const std::vector<double>& field) {
        auto parts = symmetry_components(field, reflections);
        std::stable_sort(parts.begin(), parts.end(), [&](const auto& a, const auto& b) {
            return axial_inner(a.field, a.field) > axial_inner(b.field, b.field);
        });
        const double initial = axial_inner(field, field);
        for (const auto& part : parts) {
            if (added == enrichment_modes)
                break;
            if (axial_inner(part.field, part.field) <= 1e-16 * initial || !axial_direction(part.field, part.parity))
                continue;
            SectionMode mode;
            mode.kind = SectionMode::Kind::axial_warping;
            mode.reflection_parity = part.parity;
            for (auto& coefficient : mode.coefficient)
                coefficient.resize(3 * n);
            for (std::size_t i = 0; i < n; ++i)
                mode.coefficient[0][2 * n + i] = axial_basis.back()[i];
            retain_mode(result, std::move(mode));
            ++added;
        }
    };
    const auto initial_modes = result.modes.size();
    for (std::size_t mode = 3; mode < initial_modes; ++mode) {
        const auto& field = result.modes[mode].coefficient[1];
        append_axial(std::vector<double>(field.begin() + static_cast<std::ptrdiff_t>(2 * n), field.end()));
    }
    // A clamp does not in general preserve the ratio of the x and y Poisson
    // corrections. Split the already solved classical fields into Cartesian
    // components and let them relax independently. The full field is already
    // retained, so only one component needs adding; the other is in its span.
    // The rigid-space projection above is essential on nonsymmetric sections.
    for (const auto& classic_mode : classic.modes) {
        if (added == enrichment_modes)
            break;
        std::vector<double> transverse(2 * n);
        std::copy_n(classic_mode.correction.begin(), n, transverse.begin());
        if (append(transverse, SectionMode::Kind::poisson_relaxation))
            append_axial(warping.solve(transverse).axial);
    }
    if (added < enrichment_modes) {
        const auto remaining = enrichment_modes - added;
        const auto spectrum = build_section_distortion_modes(section, (remaining + 1) / 2, remaining);
        const auto directions = std::max(spectrum.modes.size(), spectrum.shear_free_modes.size());
        // Block equilibrium corrections seeded by the two physical bending
        // axial fields. Advance each family recursively instead of exhausting
        // all soft width fields before reaching higher thickness response.
        std::array<std::vector<double>, 2> axial_seeds{axial_basis[1], axial_basis[2]};
        std::vector<std::vector<double>> transverse_seeds;
        for (std::size_t j = 3; j < independent.size(); ++j)
            transverse_seeds.push_back(independent[j]);
        SectionRelaxationSolver relaxation(section);
        for (std::size_t i = 0; i < directions && added < enrichment_modes; ++i) {
            if (!transverse_seeds.empty()) {
                auto& source = transverse_seeds[i % transverse_seeds.size()];
                const auto correction = relaxation.solve_transverse_corrector(source);
                if (append(correction.transverse, SectionMode::Kind::transverse_corrector)) {
                    source = independent.back();
                    append_axial(warping.solve(source).axial);
                }
            }
            if (added == enrichment_modes)
                break;
            auto& seed = axial_seeds[i % 2];
            const auto correction = warping.solve_axial_corrector(seed);
            const auto previous_axial_count = axial_basis.size();
            append_axial(correction.axial);
            if (axial_basis.size() > previous_axial_count)
                seed = axial_basis.back();
            else {
                seed = correction.axial;
                const double norm = std::sqrt(axial_inner(seed, seed));
                for (double& value : seed)
                    value /= norm;
            }
            if (added == enrichment_modes)
                break;
            // The same retained axial strain drives its transverse Poisson
            // equilibrium: KD*p=-integral(Bperp^T*Cperp,z*psi). Both directions
            // enter the full 3D material response with independent amplitudes.
            const auto transverse = relaxation.solve(seed).transverse;
            if (append(transverse, SectionMode::Kind::poisson_relaxation))
                append_axial(warping.solve(transverse).axial);
            if (added == enrichment_modes)
                break;
            if (i < spectrum.modes.size()) {
                const auto& mode = spectrum.modes[i];
                if (append(mode.transverse, SectionMode::Kind::distortion))
                    append_axial(mode.warping.axial);
            }
            if (i < spectrum.shear_free_modes.size() && added < enrichment_modes) {
                const auto& free = spectrum.shear_free_modes[i];
                if (append(free.transverse, SectionMode::Kind::shear_free_distortion))
                    append_axial(free.warping.axial);
            }
        }
        if (added != enrichment_modes)
            throw std::invalid_argument("Requested enrichment exceeds the independent section spectrum");
    }
    return result;
}

ReducedSectionBasis make_mixed_bending_basis(ReducedSectionBasis basis) {
    for (std::size_t direction = 0; direction < 2; ++direction) {
        const auto bending_kind =
            direction == 0 ? SectionMode::Kind::bending_x_displacement : SectionMode::Kind::bending_y_displacement;
        const auto rotation_kind = direction == 0 ? SectionMode::Kind::shear_x : SectionMode::Kind::shear_y;
        auto bending = std::find_if(basis.modes.begin(), basis.modes.end(), [bending_kind](const SectionMode& mode) {
            return mode.kind == bending_kind;
        });
        auto rotation = std::find_if(basis.modes.begin(), basis.modes.end(), [rotation_kind](const SectionMode& mode) {
            return mode.kind == rotation_kind;
        });
        if (bending == basis.modes.end() || rotation == basis.modes.end())
            throw std::invalid_argument("Mixed bending requires both classical displacement and rotation modes");
        // ux=w, uz=X*r, with gamma_xz=w'+r and epsilon_zz=X*r'.
        // The Euler-Bernoulli subspace is r=-w'. Its lateral Poisson field
        // P*w'' becomes -P*r'. Retain every such field, but avoid a third
        // derivative of transverse displacement in the production unknowns.
        // Finite shear is governed by the full 3D material, without a fitted
        // shear-correction factor. The analogous y direction is independent.
        for (std::size_t i = 0; i < bending->coefficient[2].size(); ++i) {
            rotation->coefficient[1][i] = -bending->coefficient[2][i];
            bending->coefficient[1][i] = 0.0;
            bending->coefficient[2][i] = 0.0;
        }
    }
    return basis;
}

ReducedSectionBasis make_independent_section_basis(const CrossSection& section, ReducedSectionBasis basis) {
    const auto n = section.nodes().size();
    const auto inner = [&](const std::vector<double>& left, const std::vector<double>& right) {
        double result = 0.0;
        for (const auto& point : section.points())
            for (std::size_t c = 0; c < 3; ++c) {
                double a = 0.0, b = 0.0;
                for (std::size_t i = 0; i < 8; ++i) {
                    a += point.shape[i] * left[c * n + point.nodes[i]];
                    b += point.shape[i] * right[c * n + point.nodes[i]];
                }
                result += point.weight * a * b;
            }
        return result;
    };
    std::vector<std::vector<double>> displacement;
    const auto remove = [&](std::vector<double>& field) {
        for (int pass = 0; pass < 2; ++pass)
            for (const auto& previous : displacement) {
                const double projection = inner(previous, field);
                for (std::size_t i = 0; i < 3 * n; ++i)
                    field[i] -= projection * previous[i];
            }
    };
    for (const auto& mode : basis.modes) {
        auto field = mode.coefficient[0];
        const double initial = inner(field, field);
        remove(field);
        const double norm = inner(field, field);
        if (!(norm > 1e-16 * initial))
            throw std::invalid_argument("Independent modal displacement fields are linearly dependent");
        for (double& value : field)
            value /= std::sqrt(norm);
        displacement.push_back(std::move(field));
    }
    // If Phi1_j=sum_i Phi0_i*T_ij, the physical displacement coefficient
    // is a_i=q_i+T_ij*q'_j. Discretize that coefficient independently instead
    // of retaining derivative ties between already independent fields. The
    // latter creates artificial higher-order constraints in a polynomial FE
    // space, which is not closed under differentiation with shared C2 jets.
    // This is a mixed kinematic choice, not an algebraic change of the old FE
    // unknowns. Never drop a derivative field not fully represented by Phi0.
    for (auto& mode : basis.modes)
        for (std::size_t order = 1; order < 3; ++order) {
            auto& field = mode.coefficient[order];
            const double initial = inner(field, field);
            if (initial == 0.0)
                continue;
            auto remainder = field;
            remove(remainder);
            if (inner(remainder, remainder) <= 1e-16 * initial)
                std::fill(field.begin(), field.end(), 0.0);
        }
    return basis;
}
} // namespace fuelsim
