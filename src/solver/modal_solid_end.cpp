#include "modal_solid_end.hpp"
#include <algorithm>
#include <stdexcept>

namespace fuelsim {
ModalSolidElement make_modal_solid_element(const Hex20Coordinates& coordinates,
    std::array<SolidDisplacementRow, 60> displacement,
    std::size_t layer,
    std::size_t region) {
    ModalSolidElement result;
    result.layer = layer;
    result.region = region;
    result.geometry = elements::make_c3d20t_geometry(coordinates);
    for (const auto& row : displacement)
        for (const auto& entry : row)
            result.dofs.push_back(entry.first);
    std::sort(result.dofs.begin(), result.dofs.end());
    result.dofs.erase(std::unique(result.dofs.begin(), result.dofs.end()), result.dofs.end());
    for (auto& row : displacement) {
        if (row.empty())
            throw std::invalid_argument("Solid end displacement has no retained or local representation");
        for (auto& entry : row)
            entry.first = static_cast<std::size_t>(
                std::lower_bound(result.dofs.begin(), result.dofs.end(), entry.first) - result.dofs.begin());
    }
    result.displacement = std::move(displacement);
    return result;
}

ModalSolidResponse evaluate_modal_solid(const ModalSolidElement& element,
    const SectionRegion& region,
    const std::vector<double>& state,
    bool jacobian,
    bool output) {
    const auto n = element.dofs.size();
    if (state.size() != n)
        throw std::invalid_argument("Solid end element state size mismatch");
    Hex20LocalValues nodal{}, initial{};
    std::fill_n(nodal.begin(), 8, region.temperature);
    std::fill_n(initial.begin(), 8, region.temperature);
    for (std::size_t i = 0; i < 60; ++i)
        for (const auto& entry : element.displacement[i])
            nodal[8 + i] += entry.second * state[entry.first];
    elements::C3d20Input input{region.material, element.geometry, nodal, initial};
    input.initial_temperature = region.temperature;
    const auto native = elements::evaluate_c3d20t(input, {true, jacobian, false, output});
    ModalSolidResponse result;
    result.residual.resize(n);
    if (jacobian)
        result.jacobian.resize(n * n);
    // delta u_s = T delta a gives r_a = T^T r_s and K_a = T^T K_s T.
    // The native material-point residual remains the source of virtual work;
    // evaluating K_a*a would lose accuracy for slender, point-supported plates.
    for (std::size_t i = 0; i < 60; ++i)
        for (const auto& a : element.displacement[i]) {
            result.residual[a.first] += a.second * native.residual[8 + i];
            if (jacobian)
                for (std::size_t j = 0; j < 60; ++j)
                    for (const auto& b : element.displacement[j])
                        result.jacobian[n * a.first + b.first] +=
                            a.second * native.jacobian[68 * (8 + i) + 8 + j] * b.second;
        }
    if (output) {
        for (const auto& stress : native.stress)
            result.stress.push_back({stress.xx, stress.yy, stress.zz, stress.xy, stress.yz, stress.xz});
        for (std::size_t q = 0; q < element.geometry.mechanical_points.size(); ++q) {
            const auto& point = element.geometry.mechanical_points[q];
            std::array<double, 9> gradient{};
            std::array<double, 3> u{};
            for (std::size_t c = 0; c < 3; ++c)
                for (std::size_t i = 0; i < 20; ++i) {
                    u[c] += point.displacement_shape[i] * nodal[8 + 20 * c + i];
                    for (std::size_t d = 0; d < 3; ++d)
                        gradient[3 * c + d] += point.displacement_gradient[i][d] * nodal[8 + 20 * c + i];
                }
            const SectionStrain strain{gradient[0],
                gradient[4],
                gradient[8],
                0.5 * (gradient[1] + gradient[3]),
                0.5 * (gradient[5] + gradient[7]),
                0.5 * (gradient[2] + gradient[6])};
            result.strain.push_back(strain);
            result.displacement.push_back({u[0], u[1], u[2]});
            for (std::size_t c = 0; c < 6; ++c)
                result.energy += 0.5 * point.weighted_measure * (c < 3 ? 1.0 : 2.0) * strain[c] * result.stress[q][c];
        }
    }
    return result;
}
} // namespace fuelsim
