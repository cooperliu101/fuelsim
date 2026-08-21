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
    const ActiveMatrix3 incremental = multiply(current, inverse(committed_deformation, old_determinant));
    const adlite::Scalar incremental_determinant = determinant(incremental);
    if (!std::isfinite(incremental_determinant.value()) || !(incremental_determinant.value() > 0.0))
        throw std::domain_error("Incremental finite-strain Cartesian state requires a positive Jacobian");
    const ActiveMatrix3 incremental_inverse = inverse(incremental, incremental_determinant);
    ActiveMatrix3 cinv_minus_identity{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t k = 0; k < 3; ++k)
                cinv_minus_identity[i][j] += incremental_inverse[i][k] * incremental_inverse[j][k];
            if (i == j) cinv_minus_identity[i][j] -= 1.0;
        }
    ActiveMatrix3 strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            strain[i][j] = -0.5 * cinv_minus_identity[i][j];
            for (std::size_t k = 0; k < 3; ++k)
                strain[i][j] += 0.25 * cinv_minus_identity[i][k] * cinv_minus_identity[k][j];
        }
    result.strain_increment = {strain[0][0], strain[1][1], strain[2][2], strain[0][1], strain[1][2], strain[0][2]};
    const std::array<adlite::Scalar, 3> axial = {incremental_inverse[1][2] - incremental_inverse[2][1],
        incremental_inverse[2][0] - incremental_inverse[0][2], incremental_inverse[0][1] - incremental_inverse[1][0]};
    const adlite::Scalar q = 0.25 * (axial[0] * axial[0] + axial[1] * axial[1] + axial[2] * axial[2]);
    const adlite::Scalar trace_minus_one =
        incremental_inverse[0][0] + incremental_inverse[1][1] + incremental_inverse[2][2] - 1.0;
    const adlite::Scalar p = 0.25 * trace_minus_one * trace_minus_one, sum = p + q;
    if (!std::isfinite(sum.value()) || !(sum.value() > 0.0))
        throw std::domain_error("MOOSE Taylor finite-strain rotation has invalid Cartesian p+q");
    const adlite::Scalar p2 = p * p, p3 = p2 * p, p4 = p3 * p, sum2 = sum * sum, sum3 = sum2 * sum;
    const adlite::Scalar c1_squared = p + 3.0 * p2 * (1.0 - sum) / sum2 - 2.0 * p3 * (1.0 - sum) / sum3;
    if (!std::isfinite(c1_squared.value()) || !(c1_squared.value() > 0.0))
        throw std::domain_error("MOOSE Cartesian Rashid rotation has nonpositive C1 squared");
    const adlite::Scalar c1 = adlite::sqrt(c1_squared);
    adlite::Scalar c2;
    if (q.value() > 0.01)
        c2 = (1.0 - c1) / (4.0 * q);
    else {
        const adlite::Scalar q2 = q * q, q3 = q2 * q;
        c2 = 0.125 + q * 0.03125 * (p2 - 12.0 * (p - 1.0)) / p2 + q2 * (p - 2.0) * (p2 - 10.0 * p + 32.0) / p3 +
             q3 * (1104.0 - 992.0 * p + 376.0 * p2 - 72.0 * p3 + 5.0 * p4) / (512.0 * p4);
    }
    const adlite::Scalar c3_test = (p * q * (3.0 - q) + p3 + q * q) / sum3;
    if (!std::isfinite(c3_test.value()) || !(c3_test.value() > 0.0))
        throw std::domain_error("MOOSE Cartesian Rashid rotation has nonpositive C3 test");
    const adlite::Scalar c3 = 0.5 * adlite::sqrt(c3_test);
    ActiveMatrix3 rashid{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            rashid[i][j] = c2 * axial[i] * axial[j];
            if (i == j) rashid[i][j] += c1;
        }
    rashid[0][1] += c3 * axial[2];
    rashid[0][2] -= c3 * axial[1];
    rashid[1][0] -= c3 * axial[2];
    rashid[1][2] += c3 * axial[0];
    rashid[2][0] += c3 * axial[1];
    rashid[2][1] -= c3 * axial[0];
    result.rotation = {rashid[0][0], rashid[1][0], rashid[2][0], rashid[0][1], rashid[1][1], rashid[2][1], rashid[0][2],
        rashid[1][2], rashid[2][2]};
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
