#include "core/modal_beam.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
std::array<std::array<double, 6>, 4> modal_beam_shape(double lower, double upper, double z) {
    const double length = upper - lower;
    if (!std::isfinite(length) || length <= 0.0 || !std::isfinite(z) || z < lower || z > upper)
        throw std::invalid_argument("Invalid modal beam interpolation coordinates");
    const double s = (z - lower) / length;
    std::array<std::array<double, 6>, 4> shape{};
    const std::array<std::array<double, 6>, 6> polynomial{{{1.0, 0.0, 0.0, -10.0, 15.0, -6.0},
        {0.0, length, 0.0, -6.0 * length, 8.0 * length, -3.0 * length},
        {0.0, 0.0, length * length / 2.0, -1.5 * length * length, 1.5 * length * length, -0.5 * length * length},
        {0.0, 0.0, 0.0, 10.0, -15.0, 6.0},
        {0.0, 0.0, 0.0, -4.0 * length, 7.0 * length, -3.0 * length},
        {0.0, 0.0, 0.0, 0.5 * length * length, -length * length, 0.5 * length * length}}};
    for (std::size_t node = 0; node < 6; ++node)
        for (std::size_t order = 0; order < 4; ++order)
            for (std::size_t power = order; power < 6; ++power) {
                double coefficient = polynomial[node][power];
                for (std::size_t d = 0; d < order; ++d)
                    coefficient *= static_cast<double>(power - d) / length;
                shape[order][node] += coefficient * std::pow(s, static_cast<int>(power - order));
            }
    // These are nodal jets by definition. Do not recover their endpoint
    // Kronecker values by cancellation of large inverse-length terms; tiny
    // spurious entries otherwise contaminate physical constraint rank tests
    // on a graded axial mesh. The third derivative remains the polynomial one.
    if (z == lower || z == upper) {
        const std::size_t offset = z == lower ? 0 : 3;
        for (std::size_t order = 0; order < 3; ++order) {
            shape[order].fill(0.0);
            shape[order][offset + order] = 1.0;
        }
    }
    return shape;
}

ModalBeamElement make_modal_beam_element(double lower, double upper) {
    if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower)
        throw std::invalid_argument("Modal beam requires finite increasing axial coordinates");
    ModalBeamElement element{lower, upper, {}};
    const double length = upper - lower;
    // Six-point Gauss integrates the degree-ten energy of a quintic distortion amplitude.
    constexpr std::array<double, 6> locations{-0.93246951420315202781,
        -0.66120938646626451366,
        -0.23861918608319690863,
        0.23861918608319690863,
        0.66120938646626451366,
        0.93246951420315202781};
    constexpr std::array<double, 6> weights{0.17132449237917034504,
        0.36076157304813860757,
        0.46791393457269104739,
        0.46791393457269104739,
        0.36076157304813860757,
        0.17132449237917034504};
    for (std::size_t q = 0; q < 6; ++q) {
        ModalBeamPoint point{};
        const double s = 0.5 * (1.0 + locations[q]);
        point.z = lower + s * length;
        point.weight = 0.5 * length * weights[q];
        point.shape = modal_beam_shape(lower, upper, point.z);
        element.points.push_back(point);
    }
    return element;
}

SectionKinematics modal_section_kinematics(const CrossSection& section,
    const SectionMode& mode,
    std::size_t index,
    const std::array<double, 4>& amplitude) {
    return modal_section_point_kinematics(section, mode, section.points().at(index), amplitude);
}

SectionKinematics modal_section_point_kinematics(const CrossSection& section,
    const SectionMode& mode,
    const SectionPoint& point,
    const std::array<double, 4>& amplitude) {
    const auto n = section.nodes().size();
    SectionKinematics result;
    for (double value : amplitude)
        if (!std::isfinite(value))
            throw std::invalid_argument("Modal amplitude derivatives must be finite");
    for (std::size_t order = 0; order < 3; ++order) {
        if (mode.coefficient[order].size() != 3 * n)
            throw std::invalid_argument("Section mode coefficient dimension mismatch");
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t i = 0; i < 8; ++i) {
                const double value = mode.coefficient[order][c * n + point.nodes[i]];
                if (!std::isfinite(value))
                    throw std::invalid_argument("Section mode coefficients must be finite");
                result.displacement[c] += point.shape[i] * value * amplitude[order];
                result.gradient[3 * c] += point.gradient[i][0] * value * amplitude[order];
                result.gradient[3 * c + 1] += point.gradient[i][1] * value * amplitude[order];
                result.gradient[3 * c + 2] += point.shape[i] * value * amplitude[order + 1];
            }
    }
    return result;
}

ModalSectionLinearization linearize_modal_section(const CrossSection& section, const ReducedSectionBasis& basis) {
    if (basis.modes.empty())
        throw std::invalid_argument("Cannot linearize an empty section basis");
    ModalSectionLinearization result;
    result.mode_count = basis.modes.size();
    const auto size = 4 * result.mode_count;
    result.tangent.resize(size * size);
    constexpr std::array<double, 6> metric{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    for (std::size_t q = 0; q < section.points().size(); ++q) {
        const auto& point = section.points()[q];
        const auto& region = section.regions()[point.region];
        const auto material = evaluate_stress_tangent(region.material,
            {},
            region.temperature,
            0.0,
            nullptr,
            material_context(0.0, point.position));
        std::vector<SectionStrain> b(size), cb(size);
        for (std::size_t i = 0; i < size; ++i) {
            std::array<double, 4> amplitude{};
            amplitude[i % 4] = 1.0;
            const auto g = modal_section_kinematics(section, basis.modes[i / 4], q, amplitude).gradient;
            b[i] = {g[0], g[4], g[8], 0.5 * (g[1] + g[3]), 0.5 * (g[5] + g[7]), 0.5 * (g[2] + g[6])};
            for (std::size_t a = 0; a < 6; ++a)
                for (std::size_t c = 0; c < 6; ++c)
                    cb[i][a] += material.tangent[a][c] * b[i][c];
        }
        for (std::size_t i = 0; i < size; ++i)
            for (std::size_t j = 0; j < size; ++j)
                for (std::size_t c = 0; c < 6; ++c)
                    result.tangent[size * i + j] += point.weight * metric[c] * b[i][c] * cb[j][c];
    }
    return result;
}

std::vector<double> modal_beam_linear_stiffness(const ModalSectionLinearization& section,
    const ModalBeamElement& element) {
    const auto n = section.mode_count, count = 6 * n;
    if (n == 0 || section.tangent.size() != 16 * n * n || element.points.size() != 6)
        throw std::invalid_argument("Invalid section tangent or axial integration rule");
    std::array<double, 16 * 36> products{};
    for (const auto& point : element.points)
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t b = 0; b < 4; ++b)
                for (std::size_t i = 0; i < 6; ++i)
                    for (std::size_t j = 0; j < 6; ++j)
                        products[36 * (4 * a + b) + 6 * i + j] += point.weight * point.shape[a][i] * point.shape[b][j];
    std::vector<double> stiffness(count * count);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t a = 0; a < 4; ++a)
                for (std::size_t b = 0; b < 4; ++b) {
                    const double value = section.tangent[4 * n * (4 * i + a) + 4 * j + b];
                    if (value == 0.0)
                        continue;
                    for (std::size_t left = 0; left < 6; ++left)
                        for (std::size_t right = 0; right < 6; ++right)
                            stiffness[count * (6 * i + left) + 6 * j + right] +=
                                value * products[36 * (4 * a + b) + 6 * left + right];
                }
    return stiffness;
}

ModalSectionKinematics sample_modal_section(const CrossSection& section, const ReducedSectionBasis& basis) {
    ModalSectionKinematics result;
    result.mode_count = basis.modes.size();
    result.strain.resize(section.points().size() * result.mode_count);
    result.displacement.resize(result.strain.size());
    for (std::size_t q = 0; q < section.points().size(); ++q)
        for (std::size_t mode = 0; mode < result.mode_count; ++mode)
            for (std::size_t order = 0; order < 4; ++order) {
                std::array<double, 4> amplitude{};
                amplitude[order] = 1.0;
                const auto k = modal_section_kinematics(section, basis.modes[mode], q, amplitude);
                const auto& g = k.gradient;
                result.displacement[q * result.mode_count + mode][order] = k.displacement;
                result.strain[q * result.mode_count + mode][order] =
                    {g[0], g[4], g[8], 0.5 * (g[1] + g[3]), 0.5 * (g[5] + g[7]), 0.5 * (g[2] + g[6])};
            }
    return result;
}

ModalBeamResponse evaluate_modal_beam(const CrossSection& section,
    const ReducedSectionBasis& basis,
    const ModalBeamElement& element,
    const std::vector<double>& state,
    bool include_jacobian,
    const ModalSectionKinematics* kinematics) {
    const auto count = 6 * basis.modes.size();
    if (basis.modes.empty() || state.size() != count || element.points.size() != 6)
        throw std::invalid_argument("Modal beam state or quadrature dimension mismatch");
    for (double value : state)
        if (!std::isfinite(value))
            throw std::invalid_argument("Modal beam state must be finite");
    ModalSectionKinematics sampled;
    if (!kinematics) {
        sampled = sample_modal_section(section, basis);
        kinematics = &sampled;
    }
    if (kinematics->mode_count != basis.modes.size()
        || kinematics->strain.size() != basis.modes.size() * section.points().size()
        || kinematics->displacement.size() != kinematics->strain.size())
        throw std::invalid_argument("Modal section kinematics do not match the fixed section and basis");
    ModalBeamResponse result;
    result.residual.resize(count);
    if (include_jacobian)
        result.jacobian.resize(count * count);
    result.strain.resize(6 * section.points().size());
    result.stress.resize(6 * section.points().size());
    result.displacement.resize(6 * section.points().size());
    std::vector<std::vector<std::array<double, 4>>> amplitudes(6,
        std::vector<std::array<double, 4>>(basis.modes.size()));
    const double length = element.upper - element.lower;
    for (std::size_t mode = 0; mode < basis.modes.size(); ++mode) {
        const auto first = 6 * mode;
        // Evaluate the Hermite polynomial about its left Taylor jet. Removing
        // the already represented constant/linear/quadratic part before taking
        // derivatives avoids subtracting large rigid displacements at every
        // material point of a slender, finely subdivided beam.
        const double value = state[first], slope = state[first + 1], curvature = state[first + 2];
        const std::array<double, 3> remainder{
            std::fma(-0.5 * length * length, curvature, std::fma(-length, slope, state[first + 3] - value)),
            std::fma(-length, curvature, state[first + 4] - slope),
            state[first + 5] - curvature};
        for (std::size_t z = 0; z < element.points.size(); ++z) {
            const double t = element.points[z].z - element.lower;
            auto& amplitude = amplitudes[z][mode];
            amplitude = {std::fma(0.5 * t * t, curvature, std::fma(t, slope, value)),
                std::fma(t, curvature, slope),
                curvature,
                0.0};
            for (std::size_t order = 0; order < 4; ++order)
                for (std::size_t j = 0; j < 3; ++j)
                    amplitude[order] = std::fma(element.points[z].shape[order][j + 3], remainder[j], amplitude[order]);
        }
    }
    constexpr std::array<double, 6> metric{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    for (std::size_t q = 0; q < section.points().size(); ++q) {
        const auto& point = section.points()[q];
        // Geometry-only strain coefficients of q, q', q'', q'''. Reuse them
        // across the axial quadrature and six Hermite functions. Constitutive
        // evaluation and full material-point output remain inside both loops.
        const auto* section_strain = kinematics->strain.data() + q * basis.modes.size();
        const auto* section_displacement = kinematics->displacement.data() + q * basis.modes.size();
        for (std::size_t z = 0; z < element.points.size(); ++z) {
            const auto& axial = element.points[z];
            std::vector<SectionStrain> b(count), cb(count);
            SectionStrain strain{};
            std::array<double, 3> displacement{};
            for (std::size_t mode = 0; mode < basis.modes.size(); ++mode)
                for (std::size_t order = 0; order < 4; ++order) {
                    for (std::size_t c = 0; c < 6; ++c)
                        strain[c] += section_strain[mode][order][c] * amplitudes[z][mode][order];
                    for (std::size_t c = 0; c < 3; ++c)
                        displacement[c] += section_displacement[mode][order][c] * amplitudes[z][mode][order];
                }
            for (std::size_t i = 0; i < count; ++i) {
                for (std::size_t order = 0; order < 4; ++order)
                    for (std::size_t c = 0; c < 6; ++c)
                        b[i][c] += section_strain[i / 6][order][c] * axial.shape[order][i % 6];
            }
            const auto& region = section.regions()[point.region];
            auto position = point.position;
            position.z = axial.z;
            const auto response = evaluate_stress_tangent(region.material,
                strain,
                region.temperature,
                0.0,
                nullptr,
                material_context(0.0, position));
            const SectionStrain stress{response.stress.xx,
                response.stress.yy,
                response.stress.zz,
                response.stress.xy,
                response.stress.yz,
                response.stress.xz};
            result.strain[z * section.points().size() + q] = strain;
            result.stress[z * section.points().size() + q] = stress;
            result.displacement[z * section.points().size() + q] = {displacement[0], displacement[1], displacement[2]};
            const double weight = point.weight * axial.weight;
            for (std::size_t c = 0; c < 6; ++c)
                result.energy += 0.5 * weight * metric[c] * strain[c] * stress[c];
            for (std::size_t i = 0; i < count; ++i)
                for (std::size_t c = 0; c < 6; ++c) {
                    result.residual[i] += weight * metric[c] * b[i][c] * stress[c];
                    if (include_jacobian)
                        for (std::size_t d = 0; d < 6; ++d)
                            cb[i][c] += response.tangent[c][d] * b[i][d];
                }
            if (include_jacobian)
                for (std::size_t i = 0; i < count; ++i)
                    for (std::size_t j = 0; j < count; ++j)
                        for (std::size_t c = 0; c < 6; ++c)
                            result.jacobian[count * i + j] += weight * metric[c] * b[i][c] * cb[j][c];
        }
    }
    return result;
}
} // namespace fuelsim
