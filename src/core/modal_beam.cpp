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

ModalBeamResponse evaluate_modal_beam(const CrossSection& section,
    const ReducedSectionBasis& basis,
    const ModalBeamElement& element,
    const std::vector<double>& state,
    bool include_jacobian) {
    const auto count = 6 * basis.modes.size();
    if (basis.modes.empty() || state.size() != count || element.points.size() != 6)
        throw std::invalid_argument("Modal beam state or quadrature dimension mismatch");
    for (double value : state)
        if (!std::isfinite(value))
            throw std::invalid_argument("Modal beam state must be finite");
    ModalBeamResponse result;
    result.residual.resize(count);
    if (include_jacobian)
        result.jacobian.resize(count * count);
    constexpr std::array<double, 6> metric{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    for (const auto& axial : element.points)
        for (std::size_t q = 0; q < section.points().size(); ++q) {
            const auto& point = section.points()[q];
            std::vector<SectionStrain> b(count), cb(count);
            SectionStrain strain{};
            for (std::size_t i = 0; i < count; ++i) {
                std::array<double, 4> amplitude{};
                for (std::size_t order = 0; order < 4; ++order)
                    amplitude[order] = axial.shape[order][i % 6];
                const auto k = modal_section_kinematics(section, basis.modes[i / 6], q, amplitude);
                const auto& g = k.gradient;
                b[i] = {g[0], g[4], g[8], 0.5 * (g[1] + g[3]), 0.5 * (g[5] + g[7]), 0.5 * (g[2] + g[6])};
                for (std::size_t c = 0; c < 6; ++c)
                    strain[c] += b[i][c] * state[i];
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
            result.strain.push_back(strain);
            result.stress.push_back(stress);
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
    return result;
}
} // namespace fuelsim
