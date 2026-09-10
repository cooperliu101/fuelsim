#pragma once
#include "c3d8_types.hpp"
#include "detail/ad_local_system.hpp"
#include "detail/cartesian_kinematics.hpp"
#include "material_rotation.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim::element_detail {
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

using namespace cartesian_detail;
constexpr double gauss = 0.577350269189625764509148780501957456;
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {{{{-1.0, -1.0, -1.0}},
    {{1.0, -1.0, -1.0}},
    {{1.0, 1.0, -1.0}},
    {{-1.0, 1.0, -1.0}},
    {{-1.0, -1.0, 1.0}},
    {{1.0, -1.0, 1.0}},
    {{1.0, 1.0, 1.0}},
    {{-1.0, 1.0, 1.0}}}};
constexpr std::array<std::size_t, hex8_node_count> hex8_node_to_gauss = {0, 1, 3, 2, 4, 5, 7, 6};
cartesian_detail::ActiveMatrix3 displacement_gradient(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state);
cartesian_detail::Matrix3 deformation_gradient(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);

C3d8Kinematics evaluate_cartesian_kinematics_from_gradient(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);
cartesian_detail::Matrix3 multiply_matrices(const cartesian_detail::Matrix3& first,
    const cartesian_detail::Matrix3& second);
elements::C3d8Diagnostics
diagnose_hex8(const Hex8Geometry&, const Hex8LocalValues&, const Hex8LocalValues&, StrainFormulation, bool reduced);
} // namespace fuelsim::element_detail
