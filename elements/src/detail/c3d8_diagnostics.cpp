#include "c3d8_diagnostics.hpp"
#include "c3d8_geometry.hpp"
#include "c3d8_kinematics.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::c3d8_detail {
elements::C3d8Diagnostics diagnose_hex8(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    StrainFormulation formulation,
    bool reduced) {
    elements::C3d8Diagnostics result;
    result.material_point_count = reduced ? 1 : 8;
    Hex8LocalAdValues active{}, old{};
    for (std::size_t i = 0; i < state.size(); ++i) {
        active[i] = state[i];
        old[i] = committed[i];
    }
    std::array<C3d8Kinematics, 8> quadrature;
    for (std::size_t q = 0; q < 8; ++q) {
        quadrature[q] = evaluate_cartesian_incremental_kinematics(geometry.points[q], active, committed, formulation);
        result.current_volume += quadrature[q].current_weighted_measure.value();
        result.committed_volume +=
            evaluate_cartesian_incremental_kinematics(geometry.points[q], old, committed, formulation)
                .current_weighted_measure.value();
    }
    if (formulation == StrainFormulation::small) {
        result.current_volume = geometry.reference_volume;
        result.committed_volume = geometry.reference_volume;
    }
    for (std::size_t q = 0; q < result.material_point_count; ++q) {
        const auto k =
            reduced ? evaluate_cartesian_incremental_kinematics(geometry.reduced_point, active, committed, formulation)
                    : quadrature[q];
        auto& point = result.points[q];
        point.strain_increment = {k.strain_increment.xx.value(),
            k.strain_increment.yy.value(),
            k.strain_increment.zz.value(),
            k.strain_increment.xy.value(),
            k.strain_increment.yz.value(),
            k.strain_increment.xz.value()};
        point.rotation = {k.rotation.xx.value(),
            k.rotation.xy.value(),
            k.rotation.xz.value(),
            k.rotation.yx.value(),
            k.rotation.yy.value(),
            k.rotation.yz.value(),
            k.rotation.zx.value(),
            k.rotation.zy.value(),
            k.rotation.zz.value()};
        point.current_weighted_measure = k.current_weighted_measure.value();
        if (!reduced) {
            point.temperature = state[hex8_node_gauss_permutation[q]];
            for (std::size_t n = 0; n < 8; ++n)
                for (std::size_t d = 0; d < 3; ++d)
                    point.thermal_gradient[n][d] = k.current_gradient[n][d].value();
        } else {
            for (std::size_t g = 0; g < 8; ++g) {
                const double measure = formulation == StrainFormulation::finite
                                           ? quadrature[g].current_weighted_measure.value()
                                           : geometry.points[g].weighted_measure;
                for (std::size_t n = 0; n < 8; ++n) {
                    point.temperature += measure * geometry.points[g].shape[n] * state[n] / result.current_volume;
                    for (std::size_t d = 0; d < 3; ++d)
                        point.thermal_gradient[n][d] +=
                            measure * quadrature[g].current_gradient[n][d].value() / result.current_volume;
                }
            }
        }
    }
    return result;
}
} // namespace fuelsim::c3d8_detail
