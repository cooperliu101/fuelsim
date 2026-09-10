#pragma once
#include "c3d8_types.hpp"

namespace fuelsim {
// Public geometry diagnostics used for host conservation accounting.
struct C3d8Kinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> current_gradient;
    adlite::Scalar current_weighted_measure{0.0};
};

C3d8Kinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);

} // namespace fuelsim
