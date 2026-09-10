#include "c3d20_kinematics.hpp"
#include "cartesian_kinematics.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::c3d20_detail {
using namespace cartesian_detail;

ActiveMatrix3 displacement_gradient(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state) {
    ActiveMatrix3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 20; ++node)
                result[component][direction] +=
                    point.displacement_gradient[node][direction] * state[8 + 20 * component + node];
    return result;
}

namespace {
Matrix3 deformation_gradient(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    Matrix3 result{};
    for (std::size_t component = 0; component < 3; ++component) {
        result[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 20; ++node)
                result[component][direction] +=
                    point.displacement_gradient[node][direction] * state[8 + 20 * component + node];
    }
    return result;
}

Matrix3 displacement_gradient_values(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    Matrix3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 20; ++node)
                result[component][direction] +=
                    point.displacement_gradient[node][direction] * state[8 + 20 * component + node];
    return result;
}

} // namespace

Hex20KinematicsValues evaluate_kinematics_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    Hex20KinematicsValues result{};
    const Matrix3 gradient = displacement_gradient_values(point, state);
    Matrix3 current_inverse{};
    for (std::size_t index = 0; index < 3; ++index)
        current_inverse[index][index] = 1.0;
    double current_determinant = 1.0;
    if (strain_formulation == StrainFormulation::small) {
        result.strain_increment = {gradient[0][0],
            gradient[1][1],
            gradient[2][2],
            0.5 * (gradient[0][1] + gradient[1][0]),
            0.5 * (gradient[1][2] + gradient[2][1]),
            0.5 * (gradient[0][2] + gradient[2][0])};
    } else {
        Matrix3 current = gradient;
        for (std::size_t index = 0; index < 3; ++index)
            current[index][index] += 1.0;
        current_determinant = determinant(current);
        if (!std::isfinite(current_determinant) || !(current_determinant > 0.0))
            throw std::domain_error("Finite-strain Cartesian deformation must preserve a positive Jacobian");
        current_inverse = inverse(current, current_determinant);
        const Matrix3 old = deformation_gradient(point, committed_state);
        const double old_determinant = determinant(old);
        if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
            throw std::domain_error("Committed finite-strain Cartesian state requires a positive Jacobian");
        const double incremental_determinant = current_determinant / old_determinant;
        if (!std::isfinite(incremental_determinant) || !(incremental_determinant > 0.0))
            throw std::domain_error("Incremental finite-strain Cartesian state requires a positive Jacobian");
        Matrix3 deformation_sum{}, deformation_difference{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                deformation_sum[i][j] = current[i][j] + old[i][j];
                deformation_difference[i][j] = current[i][j] - old[i][j];
            }
        const double plus_determinant = determinant(deformation_sum);
        if (!std::isfinite(plus_determinant) || plus_determinant == 0.0)
            throw std::domain_error("Abaqus Hughes-Winget Cartesian increment has singular delta-F plus identity");
        const Matrix3 plus_inverse = inverse(deformation_sum, plus_determinant);
        Matrix3 hughes_winget = cartesian_detail::multiply(deformation_difference, plus_inverse);
        for (auto& row : hughes_winget)
            for (double& value : row)
                value *= 2.0;
        Matrix3 spatial_strain{}, rotation_numerator{}, rotation_denominator{};
        for (std::size_t index = 0; index < 3; ++index) {
            rotation_numerator[index][index] = 1.0;
            rotation_denominator[index][index] = 1.0;
        }
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                spatial_strain[i][j] = 0.5 * (hughes_winget[i][j] + hughes_winget[j][i]);
                const double half_spin = 0.25 * (hughes_winget[i][j] - hughes_winget[j][i]);
                rotation_numerator[i][j] += half_spin;
                rotation_denominator[i][j] -= half_spin;
            }
        const double rotation_denominator_determinant = determinant(rotation_denominator);
        if (!std::isfinite(rotation_denominator_determinant) || rotation_denominator_determinant == 0.0)
            throw std::domain_error("Abaqus Hughes-Winget Cartesian rotation denominator is singular");
        const Matrix3 rotation = cartesian_detail::multiply(rotation_numerator,
            inverse(rotation_denominator, rotation_denominator_determinant));
        const Matrix3 spatial_times_rotation = cartesian_detail::multiply(spatial_strain, rotation);
        Matrix3 corotational_strain{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = i; j < 3; ++j)
                for (std::size_t k = 0; k < 3; ++k)
                    corotational_strain[i][j] += rotation[k][i] * spatial_times_rotation[k][j];
        result.strain_increment = {corotational_strain[0][0],
            corotational_strain[1][1],
            corotational_strain[2][2],
            corotational_strain[0][1],
            corotational_strain[1][2],
            corotational_strain[0][2]};
        result.rotation = {rotation[0][0],
            rotation[0][1],
            rotation[0][2],
            rotation[1][0],
            rotation[1][1],
            rotation[1][2],
            rotation[2][0],
            rotation[2][1],
            rotation[2][2]};
    }
    result.current_weighted_measure = point.weighted_measure * current_determinant;
    for (std::size_t node = 0; node < 20; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.current_gradient[node][direction] +=
                    point.displacement_gradient[node][reference] * current_inverse[reference][direction];
    return result;
}

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point,
    const ActiveMatrix3& gradient,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    Hex20Kinematics result{};
    const Matrix3 old = deformation_gradient(point, committed_state);
    const cartesian_detail::KinematicsCore core =
        cartesian_detail::evaluate_kinematics(gradient, old, strain_formulation);
    result.strain_increment = core.strain_increment;
    result.rotation = core.rotation;
    result.current_weighted_measure = point.weighted_measure * core.current_determinant;
    for (std::size_t node = 0; node < 20; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.current_gradient[node][direction] +=
                    point.displacement_gradient[node][reference] * core.current_inverse[reference][direction];
    if (strain_formulation == StrainFormulation::finite)
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                result.current_inverse_values[i][j] = core.current_inverse[i][j].value();
                result.midpoint_inverse_values[i][j] = core.midpoint_inverse[i][j].value();
            }
    return result;
}

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalAdValues& state,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    return evaluate_kinematics(point, displacement_gradient(point, state), committed_state, strain_formulation);
}

Hex20SourceMeasureValues evaluate_source_measure_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state) {
    Matrix3 current{};
    for (std::size_t direction = 0; direction < 3; ++direction)
        current[direction][direction] = 1.0;
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                current[component][direction] +=
                    point.source_displacement_gradient[node][direction] * state[8 + 20 * component + node];
    const double determinant_value = determinant(current);
    if (!std::isfinite(determinant_value) || !(determinant_value > 0.0))
        throw std::domain_error("Finite-strain HEX20 thermal source corner configuration requires a positive Jacobian");
    Matrix3 cofactor{};
    cofactor[0][0] = current[1][1] * current[2][2] - current[1][2] * current[2][1];
    cofactor[0][1] = current[1][2] * current[2][0] - current[1][0] * current[2][2];
    cofactor[0][2] = current[1][0] * current[2][1] - current[1][1] * current[2][0];
    cofactor[1][0] = current[0][2] * current[2][1] - current[0][1] * current[2][2];
    cofactor[1][1] = current[0][0] * current[2][2] - current[0][2] * current[2][0];
    cofactor[1][2] = current[0][1] * current[2][0] - current[0][0] * current[2][1];
    cofactor[2][0] = current[0][1] * current[1][2] - current[0][2] * current[1][1];
    cofactor[2][1] = current[0][2] * current[1][0] - current[0][0] * current[1][2];
    cofactor[2][2] = current[0][0] * current[1][1] - current[0][1] * current[1][0];
    for (auto& row : cofactor)
        for (double& value : row)
            value *= point.source_weighted_measure;
    return {cofactor, point.source_weighted_measure * determinant_value};
}

std::array<std::array<double, 3>, 8> temperature_shape_gradients(const Hex20MechanicalQuadraturePoint& point,
    const Matrix3& inverse_map) {
    std::array<std::array<double, 3>, 8> midpoint_gradient{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                midpoint_gradient[node][direction] +=
                    point.temperature_gradient[node][reference] * inverse_map[reference][direction];
    return midpoint_gradient;
}

Hex20FiniteThermalKinematicsValues evaluate_finite_thermal_kinematics_values(
    const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state) {
    Hex20FiniteThermalKinematicsValues result;
    const Matrix3 current = deformation_gradient(point, state), old = deformation_gradient(point, committed_state);
    const double current_determinant = determinant(current), old_determinant = determinant(old);
    if (!std::isfinite(current_determinant) || !(current_determinant > 0.0))
        throw std::domain_error("Finite-strain HEX20 thermal current configuration requires a positive Jacobian");
    if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
        throw std::domain_error("Finite-strain HEX20 thermal committed configuration requires a positive Jacobian");
    Matrix3 midpoint{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            midpoint[i][j] = 0.5 * (current[i][j] + old[i][j]);
    const double midpoint_determinant = determinant(midpoint);
    if (!std::isfinite(midpoint_determinant) || !(midpoint_determinant > 0.0))
        throw std::domain_error("Finite-strain HEX20 thermal midpoint configuration requires a positive Jacobian");
    const Matrix3 midpoint_inverse = inverse(midpoint, midpoint_determinant);
    result.midpoint_temperature_gradient = temperature_shape_gradients(point, midpoint_inverse);
    result.current_weighted_measure = point.weighted_measure * current_determinant;

    return result;
}

} // namespace fuelsim::c3d20_detail

namespace fuelsim {
using namespace c3d20_detail;
using cartesian_detail::determinant;

void validate_hex20_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    const double value = determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX20 deformation must preserve a positive Jacobian");
}
} // namespace fuelsim
