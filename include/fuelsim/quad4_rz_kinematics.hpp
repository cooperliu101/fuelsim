#ifndef FUELSIM_QUAD4_RZ_KINEMATICS_HPP
#define FUELSIM_QUAD4_RZ_KINEMATICS_HPP

#include "fuelsim/material.hpp"
#include "fuelsim/quad4_rz.hpp"

#include <array>

namespace fuelsim {

enum class StrainFormulation {
    small,
    finite,
};

struct AxisymmetricKinematics final {
    std::array<adlite::Scalar, quad4_node_count> gradient_r;
    std::array<adlite::Scalar, quad4_node_count> gradient_z;
    adlite::Scalar radius;
    adlite::Scalar weighted_measure;
    adlite::Scalar strain_rr;
    adlite::Scalar strain_zz;
    adlite::Scalar strain_hoop;
    adlite::Scalar strain_rz;
    AxisymmetricRotation rotation;
};

AxisymmetricKinematics
evaluate_axisymmetric_kinematics(const RzQuadraturePoint& point,
                                 const LocalAdValues& state,
                                 StrainFormulation strain_formulation);

AxisymmetricKinematics evaluate_axisymmetric_incremental_kinematics(
    const RzQuadraturePoint& point, const LocalAdValues& current_state,
    const LocalValues& committed_state, StrainFormulation strain_formulation);

} // namespace fuelsim

#endif
