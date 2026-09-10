#pragma once
#include "cartesian_types.hpp"
#include "element_types.hpp"
#include "matrix3.hpp"
#include <adlite/adlite.hpp>
#include <array>

namespace fuelsim::cartesian_detail {
struct KinematicsCore final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    ActiveMatrix3 current_inverse{};
    ActiveMatrix3 midpoint_inverse{};
    adlite::Scalar current_determinant{1.0};
};

KinematicsCore evaluate_kinematics(const ActiveMatrix3& gradient,
    const Matrix3& committed_deformation,
    StrainFormulation strain_formulation);
// This increment-only helper populates strain_increment and rotation. It cannot reconstruct the current
// configuration, so current_inverse and current_determinant retain their default values and must not be read.
KinematicsCore evaluate_hughes_winget_increment(const ActiveMatrix3& central_displacement_gradient);

} // namespace fuelsim::cartesian_detail
