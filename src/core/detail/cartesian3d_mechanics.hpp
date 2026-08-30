#pragma once
#include "fuelsim/core/kinematics.hpp"
#include "fuelsim/core/material.hpp"
#include "fuelsim/core/mesh.hpp"
#include <adlite/adlite.hpp>
#include <array>

namespace fuelsim::cartesian_detail {
using Matrix3 = std::array<std::array<double, 3>, 3>;
using ActiveMatrix3 = std::array<std::array<adlite::Scalar, 3>, 3>;

double determinant(const Matrix3& matrix);
adlite::Scalar determinant(const ActiveMatrix3& matrix);
Matrix3 inverse(const Matrix3& matrix, double determinant_value);
ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value);
ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second);

struct KinematicsCore final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    ActiveMatrix3 current_inverse{};
    adlite::Scalar current_determinant{1.0};
};

KinematicsCore evaluate_kinematics(
    const ActiveMatrix3& gradient, const Matrix3& committed_deformation, StrainFormulation strain_formulation);
// This increment-only helper populates strain_increment and rotation. It cannot reconstruct the current
// configuration, so current_inverse and current_determinant retain their default values and must not be read.
KinematicsCore evaluate_hughes_winget_increment(const ActiveMatrix3& central_displacement_gradient);

MaterialFunctionContext material_context(double time, const CartesianPoint3& point);

struct CartesianStressTangent final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 6>, 6> tangent{};
    std::array<double, 6> thermal{};
};

CartesianStressTangent evaluate_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 6>& fed_strain, double temperature, double time_step,
    const CartesianMaterialPointState* committed_material, MaterialFunctionContext context);
} // namespace fuelsim::cartesian_detail
