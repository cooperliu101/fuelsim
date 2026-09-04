#include "cartesian3d_mechanics.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::cartesian_detail {
namespace {
ActiveMatrix3 identity_active_matrix() {
    ActiveMatrix3 result{};
    for (std::size_t index = 0; index < 3; ++index) result[index][index] = 1.0;
    return result;
}
} // namespace

double determinant(const Matrix3& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

adlite::Scalar determinant(const ActiveMatrix3& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

Matrix3 inverse(const Matrix3& matrix, double determinant_value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant_value,
                 (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant_value,
                 (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant_value}},
        {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant_value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant_value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant_value}},
        {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant_value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant_value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant_value}}}};
}

ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant_value,
                 (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant_value,
                 (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant_value}},
        {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant_value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant_value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant_value}},
        {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant_value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant_value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant_value}}}};
}

ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second) {
    ActiveMatrix3 result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k) result[i][j] += first[i][k] * second[k][j];
    return result;
}

KinematicsCore evaluate_hughes_winget_increment(const ActiveMatrix3& central_displacement_gradient) {
    KinematicsCore result{};
    const ActiveMatrix3& hughes_winget = central_displacement_gradient;
    ActiveMatrix3 spatial_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) spatial_strain[i][j] = 0.5 * (hughes_winget[i][j] + hughes_winget[j][i]);

    ActiveMatrix3 rotation_numerator = identity_active_matrix();
    ActiveMatrix3 rotation_denominator = identity_active_matrix();
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const adlite::Scalar half_spin = 0.25 * (hughes_winget[i][j] - hughes_winget[j][i]);
            rotation_numerator[i][j] += half_spin;
            rotation_denominator[i][j] -= half_spin;
        }
    const adlite::Scalar rotation_denominator_determinant = determinant(rotation_denominator);
    if (!std::isfinite(rotation_denominator_determinant.value()) || rotation_denominator_determinant.value() == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian rotation denominator is singular");
    const ActiveMatrix3 rotation_denominator_inverse = inverse(rotation_denominator, rotation_denominator_determinant);
    ActiveMatrix3 rotation{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                rotation[i][j] += rotation_numerator[i][k] * rotation_denominator_inverse[k][j];
    ActiveMatrix3 spatial_times_rotation{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k) spatial_times_rotation[i][j] += spatial_strain[i][k] * rotation[k][j];
    ActiveMatrix3 corotational_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                corotational_strain[i][j] += rotation[k][i] * spatial_times_rotation[k][j];
    result.strain_increment = {corotational_strain[0][0], corotational_strain[1][1], corotational_strain[2][2],
        corotational_strain[0][1], corotational_strain[1][2], corotational_strain[0][2]};
    result.rotation = {rotation[0][0], rotation[0][1], rotation[0][2], rotation[1][0], rotation[1][1], rotation[1][2],
        rotation[2][0], rotation[2][1], rotation[2][2]};
    return result;
}

KinematicsCore evaluate_kinematics(
    const ActiveMatrix3& gradient, const Matrix3& committed_deformation, StrainFormulation strain_formulation) {
    KinematicsCore result{};
    result.current_inverse = identity_active_matrix();
    if (strain_formulation == StrainFormulation::small) {
        result.strain_increment = {gradient[0][0], gradient[1][1], gradient[2][2],
            0.5 * (gradient[0][1] + gradient[1][0]), 0.5 * (gradient[1][2] + gradient[2][1]),
            0.5 * (gradient[0][2] + gradient[2][0])};
        return result;
    }
    ActiveMatrix3 current = gradient;
    for (std::size_t direction = 0; direction < 3; ++direction) current[direction][direction] += 1.0;
    result.current_determinant = determinant(current);
    if (!std::isfinite(result.current_determinant.value()) || !(result.current_determinant.value() > 0.0))
        throw std::domain_error("Finite-strain Cartesian deformation must preserve a positive Jacobian");
    result.current_inverse = inverse(current, result.current_determinant);
    const double old_determinant = determinant(committed_deformation);
    if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
        throw std::domain_error("Committed finite-strain Cartesian state requires a positive Jacobian");
    const double incremental_determinant = result.current_determinant.value() / old_determinant;
    if (!std::isfinite(incremental_determinant) || !(incremental_determinant > 0.0))
        throw std::domain_error("Incremental finite-strain Cartesian state requires a positive Jacobian");
    // (F_new F_old^-1 - I)(F_new F_old^-1 + I)^-1 is exactly
    // (F_new - F_old)(F_new + F_old)^-1.  The latter avoids an active
    // inverse and matrix product while retaining the same Hughes-Winget map.
    ActiveMatrix3 deformation_sum{}, deformation_difference{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            deformation_sum[i][j] = current[i][j] + committed_deformation[i][j];
            deformation_difference[i][j] = current[i][j] - committed_deformation[i][j];
        }
    const adlite::Scalar plus_determinant = determinant(deformation_sum);
    if (!std::isfinite(plus_determinant.value()) || plus_determinant.value() == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian increment has singular delta-F plus identity");
    const ActiveMatrix3 plus_inverse = inverse(deformation_sum, plus_determinant);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) result.midpoint_inverse[i][j] = 2.0 * plus_inverse[i][j];
    ActiveMatrix3 hughes_winget{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                hughes_winget[i][j] += 2.0 * deformation_difference[i][k] * plus_inverse[k][j];
    ActiveMatrix3 spatial_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) spatial_strain[i][j] = 0.5 * (hughes_winget[i][j] + hughes_winget[j][i]);

    ActiveMatrix3 rotation_numerator = identity_active_matrix();
    ActiveMatrix3 rotation_denominator = identity_active_matrix();
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const adlite::Scalar half_spin = 0.25 * (hughes_winget[i][j] - hughes_winget[j][i]);
            rotation_numerator[i][j] += half_spin;
            rotation_denominator[i][j] -= half_spin;
        }
    const adlite::Scalar rotation_denominator_determinant = determinant(rotation_denominator);
    if (!std::isfinite(rotation_denominator_determinant.value()) || rotation_denominator_determinant.value() == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian rotation denominator is singular");
    const ActiveMatrix3 rotation_denominator_inverse = inverse(rotation_denominator, rotation_denominator_determinant);
    ActiveMatrix3 rotation{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                rotation[i][j] += rotation_numerator[i][k] * rotation_denominator_inverse[k][j];
    ActiveMatrix3 spatial_times_rotation{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k) spatial_times_rotation[i][j] += spatial_strain[i][k] * rotation[k][j];
    ActiveMatrix3 corotational_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                corotational_strain[i][j] += rotation[k][i] * spatial_times_rotation[k][j];
    result.strain_increment = {corotational_strain[0][0], corotational_strain[1][1], corotational_strain[2][2],
        corotational_strain[0][1], corotational_strain[1][2], corotational_strain[0][2]};
    result.rotation = {rotation[0][0], rotation[0][1], rotation[0][2], rotation[1][0], rotation[1][1], rotation[1][2],
        rotation[2][0], rotation[2][1], rotation[2][2]};
    return result;
}

MaterialFunctionContext material_context(double time, const CartesianPoint3& point) {
    return {time, point.x, point.y, point.z};
}

CartesianStressTangent evaluate_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 6>& fed_strain, double temperature, double time_step,
    const CartesianMaterialPointState* committed_material, MaterialFunctionContext context) {
    std::array<double, 7> seeds{};
    for (std::size_t component = 0; component < 6; ++component) seeds[component] = fed_strain[component];
    seeds[6] = temperature;
    std::array<adlite::Scalar, 7> active{};
    adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    const SymmetricTensor3 strain{active[0], active[1], active[2], active[3], active[4], active[5]};
    const SymmetricTensor3 stress =
        committed_material == nullptr
            ? material.stress(strain, active[6], context)
            : material.response(strain, active[6], time_step, *committed_material, context).stress;
    const std::array<const adlite::Scalar*, 6> components = {
        &stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz};
    CartesianStressTangent result{};
    result.stress = {stress.xx.value(), stress.yy.value(), stress.zz.value(), stress.xy.value(), stress.yz.value(),
        stress.xz.value()};
    std::array<double, 7> derivatives{};
    for (std::size_t row = 0; row < 6; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 6; ++column) result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[6];
    }
    return result;
}
} // namespace fuelsim::cartesian_detail
