#include "fuelsim/core/cartesian3d_hex20.hpp"
#include "detail/ad_local_system.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
using ActiveMatrix3 = std::array<std::array<adlite::Scalar, 3>, 3>;
using Matrix3 = std::array<std::array<double, 3>, 3>;
constexpr double gauss3 = 0.774596669241483377035853079956479922;
constexpr double gauss2 = 0.577350269189625764509148780502;
constexpr std::array<double, 2> gauss2_points = {-gauss2, gauss2};
constexpr std::array<double, 2> gauss2_weights = {1.0, 1.0};
constexpr std::array<double, 3> gauss3_points = {-gauss3, 0.0, gauss3};
constexpr std::array<double, 3> gauss3_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
constexpr std::array<std::array<double, 3>, 8> corner_signs = {
    {{{-1.0, -1.0, -1.0}}, {{1.0, -1.0, -1.0}}, {{1.0, 1.0, -1.0}}, {{-1.0, 1.0, -1.0}}, {{-1.0, -1.0, 1.0}},
        {{1.0, -1.0, 1.0}}, {{1.0, 1.0, 1.0}}, {{-1.0, 1.0, 1.0}}}};

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

Matrix3 inverse(const Matrix3& matrix, double value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / value,
                 (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / value,
                 (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / value}},
        {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / value}},
        {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / value}}}};
}

ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / value,
                 (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / value,
                 (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / value}},
        {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / value}},
        {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / value}}}};
}

ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second) {
    ActiveMatrix3 result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k) result[i][j] += first[i][k] * second[k][j];
    return result;
}

void evaluate_hex20_shapes(double xi, double eta, double zeta, std::array<double, 20>& shape,
    std::array<std::array<double, 3>, 20>& derivative) {
    for (std::size_t node = 0; node < 8; ++node) {
        const double sx = corner_signs[node][0], sy = corner_signs[node][1], sz = corner_signs[node][2];
        const double ax = 1.0 + sx * xi, ay = 1.0 + sy * eta, az = 1.0 + sz * zeta;
        const double sum = sx * xi + sy * eta + sz * zeta - 2.0;
        shape[node] = 0.125 * ax * ay * az * sum;
        derivative[node][0] = 0.125 * sx * ay * az * (sum + ax);
        derivative[node][1] = 0.125 * sy * ax * az * (sum + ay);
        derivative[node][2] = 0.125 * sz * ax * ay * (sum + az);
    }
    const auto xi_edge = [&](std::size_t node, double sy, double sz) {
        shape[node] = 0.25 * (1.0 - xi * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta);
        derivative[node] = {{-0.5 * xi * (1.0 + sy * eta) * (1.0 + sz * zeta),
            0.25 * sy * (1.0 - xi * xi) * (1.0 + sz * zeta), 0.25 * sz * (1.0 - xi * xi) * (1.0 + sy * eta)}};
    };
    const auto eta_edge = [&](std::size_t node, double sx, double sz) {
        shape[node] = 0.25 * (1.0 - eta * eta) * (1.0 + sx * xi) * (1.0 + sz * zeta);
        derivative[node] = {{0.25 * sx * (1.0 - eta * eta) * (1.0 + sz * zeta),
            -0.5 * eta * (1.0 + sx * xi) * (1.0 + sz * zeta), 0.25 * sz * (1.0 - eta * eta) * (1.0 + sx * xi)}};
    };
    const auto zeta_edge = [&](std::size_t node, double sx, double sy) {
        shape[node] = 0.25 * (1.0 - zeta * zeta) * (1.0 + sx * xi) * (1.0 + sy * eta);
        derivative[node] = {{0.25 * sx * (1.0 - zeta * zeta) * (1.0 + sy * eta),
            0.25 * sy * (1.0 - zeta * zeta) * (1.0 + sx * xi), -0.5 * zeta * (1.0 + sx * xi) * (1.0 + sy * eta)}};
    };
    xi_edge(8, -1.0, -1.0);
    eta_edge(9, 1.0, -1.0);
    xi_edge(10, 1.0, -1.0);
    eta_edge(11, -1.0, -1.0);
    zeta_edge(12, -1.0, -1.0);
    zeta_edge(13, 1.0, -1.0);
    zeta_edge(14, 1.0, 1.0);
    zeta_edge(15, -1.0, 1.0);
    xi_edge(16, -1.0, 1.0);
    eta_edge(17, 1.0, 1.0);
    xi_edge(18, 1.0, 1.0);
    eta_edge(19, -1.0, 1.0);
}

void evaluate_hex8_temperature_shapes(double xi, double eta, double zeta, std::array<double, 8>& shape,
    std::array<std::array<double, 3>, 8>& derivative) {
    for (std::size_t node = 0; node < 8; ++node) {
        const double sx = corner_signs[node][0], sy = corner_signs[node][1], sz = corner_signs[node][2];
        shape[node] = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta);
        derivative[node] = {{0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
            0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta), 0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)}};
    }
}

ActiveMatrix3 displacement_gradient(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state) {
    ActiveMatrix3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 20; ++node)
                result[component][direction] +=
                    point.displacement_gradient[node][direction] * state[8 + 20 * component + node];
    return result;
}

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

struct Hex20Kinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, 20> current_gradient;
    adlite::Scalar current_weighted_measure{0.0};
};

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point, const ActiveMatrix3& gradient,
    const Hex20LocalValues& committed_state, StrainFormulation strain_formulation) {
    Hex20Kinematics result{};
    if (strain_formulation == StrainFormulation::small) {
        result.strain_increment = {gradient[0][0], gradient[1][1], gradient[2][2],
            0.5 * (gradient[0][1] + gradient[1][0]), 0.5 * (gradient[1][2] + gradient[2][1]),
            0.5 * (gradient[0][2] + gradient[2][0])};
        for (std::size_t node = 0; node < 20; ++node)
            for (std::size_t direction = 0; direction < 3; ++direction)
                result.current_gradient[node][direction] = point.displacement_gradient[node][direction];
        result.current_weighted_measure = point.weighted_measure;
        return result;
    }
    ActiveMatrix3 current = gradient;
    for (std::size_t direction = 0; direction < 3; ++direction) current[direction][direction] += 1.0;
    const adlite::Scalar current_determinant = determinant(current);
    if (!std::isfinite(current_determinant.value()) || !(current_determinant.value() > 0.0))
        throw std::domain_error("Finite-strain HEX20 deformation must preserve a positive Jacobian");
    const ActiveMatrix3 current_inverse = inverse(current, current_determinant);
    for (std::size_t node = 0; node < 20; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.current_gradient[node][direction] +=
                    point.displacement_gradient[node][reference] * current_inverse[reference][direction];
    result.current_weighted_measure = point.weighted_measure * current_determinant;
    const Matrix3 old = deformation_gradient(point, committed_state);
    const double old_determinant = determinant(old);
    if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
        throw std::domain_error("Committed finite-strain HEX20 state requires a positive Jacobian");
    const ActiveMatrix3 incremental = multiply(current, inverse(old, old_determinant));
    const adlite::Scalar incremental_determinant = determinant(incremental);
    if (!std::isfinite(incremental_determinant.value()) || !(incremental_determinant.value() > 0.0))
        throw std::domain_error("Incremental finite-strain HEX20 state requires a positive Jacobian");
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
        throw std::domain_error("MOOSE Taylor finite-strain rotation has invalid HEX20 p+q");
    const adlite::Scalar p2 = p * p, p3 = p2 * p, p4 = p3 * p, sum2 = sum * sum, sum3 = sum2 * sum;
    const adlite::Scalar c1_squared = p + 3.0 * p2 * (1.0 - sum) / sum2 - 2.0 * p3 * (1.0 - sum) / sum3;
    if (!std::isfinite(c1_squared.value()) || !(c1_squared.value() > 0.0))
        throw std::domain_error("MOOSE HEX20 Rashid rotation has nonpositive C1 squared");
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
        throw std::domain_error("MOOSE HEX20 Rashid rotation has nonpositive C3 test");
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

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state,
    const Hex20LocalValues& committed_state, StrainFormulation strain_formulation) {
    return evaluate_kinematics(point, displacement_gradient(point, state), committed_state, strain_formulation);
}

MaterialFunctionContext material_context(double time, const CartesianPoint3& point) {
    return {time, point.x, point.y, point.z};
}

struct CartesianStressTangent final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 6>, 6> tangent{};
    std::array<double, 6> thermal{};
};

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

adlite::Scalar interpolate_temperature(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node) result += point.temperature_shape[node] * state[node];
    return result;
}

adlite::Scalar interpolate_temperature(const Hex20ThermalQuadraturePoint& point, const Hex20LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node) result += point.temperature_shape[node] * state[node];
    return result;
}

void add_thermal_point_residual(const Hex20ThermalQuadraturePoint& point, const Hex20LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material, double time, double volumetric_heat_source,
    const Hex20LocalValues* committed_state, double time_step, Hex20LocalAdValues& residual,
    bool include_thermal_time_term) {
    const adlite::Scalar temperature = interpolate_temperature(point, state);
    std::array<adlite::Scalar, 3> temperature_gradient{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += point.temperature_gradient[node][direction] * state[node];
    const MaterialFunctionContext context = material_context(time, point.position);
    const adlite::Scalar conductivity = material.conductivity(temperature, context);
    adlite::Scalar temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr && include_thermal_time_term) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node)
            old_temperature += point.temperature_shape[node] * (*committed_state)[node];
        temperature_rate = (temperature - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(temperature, context);
    }
    for (std::size_t node = 0; node < 8; ++node) {
        adlite::Scalar conduction = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            conduction += point.temperature_gradient[node][direction] * temperature_gradient[direction];
        residual[node] +=
            point.weighted_measure *
            (conductivity * conduction + point.temperature_shape[node] * heat_capacity * temperature_rate -
                point.temperature_shape[node] * volumetric_heat_source);
    }
}

void add_mechanical_point_residual(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time,
    const Hex20LocalValues* committed_state, const CartesianMaterialPointState* committed_material, double time_step,
    Hex20LocalAdValues& residual) {
    const adlite::Scalar temperature = interpolate_temperature(point, state);
    const MaterialFunctionContext context = material_context(time, point.position);
    const Hex20LocalValues undeformed{};
    const Hex20LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const Hex20Kinematics kinematics = evaluate_kinematics(point, state, old_state, strain_formulation);
    SymmetricTensor3 stress;
    if (committed_material == nullptr) {
        stress = material.stress(kinematics.strain_increment, temperature, context);
        if (strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
    } else if (strain_formulation == StrainFormulation::finite) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node) old_temperature += point.temperature_shape[node] * old_state[node];
        stress = material
                     .incremental_response(kinematics.strain_increment, kinematics.rotation, temperature,
                         old_temperature, time_step, *committed_material, context)
                     .stress;
    } else {
        stress =
            material.response(kinematics.strain_increment, temperature, time_step, *committed_material, context).stress;
    }
    for (std::size_t node = 0; node < 20; ++node) {
        const adlite::Scalar gx = kinematics.current_gradient[node][0], gy = kinematics.current_gradient[node][1],
                             gz = kinematics.current_gradient[node][2];
        residual[8 + node] += kinematics.current_weighted_measure * (stress.xx * gx + stress.xy * gy + stress.xz * gz);
        residual[28 + node] += kinematics.current_weighted_measure * (stress.xy * gx + stress.yy * gy + stress.yz * gz);
        residual[48 + node] += kinematics.current_weighted_measure * (stress.xz * gx + stress.yz * gy + stress.zz * gz);
    }
}

void add_thermal_point_system(const Hex20ThermalQuadraturePoint& point, const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material, double time, double volumetric_heat_source,
    const Hex20LocalValues* committed_state, double time_step, Hex20LocalAdValues& residual,
    Hex20LocalJacobian& jacobian, bool include_thermal_time_term) {
    constexpr std::size_t point_width = 1, temperature_index = 0;
    std::array<adlite::Scalar, 1> active_temperature{};
    double temperature_value = 0.0;
    for (std::size_t node = 0; node < 8; ++node) temperature_value += point.temperature_shape[node] * state[node];
    active_temperature[0] = adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    std::array<adlite::Scalar, 3> temperature_gradient{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += point.temperature_gradient[node][direction] * state[node];
    const MaterialFunctionContext context = material_context(time, point.position);
    const adlite::Scalar conductivity = material.conductivity(active_temperature[0], context);
    adlite::Scalar temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr && include_thermal_time_term) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node)
            old_temperature += point.temperature_shape[node] * (*committed_state)[node];
        temperature_rate = (active_temperature[0] - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(active_temperature[0], context);
    }
    Hex20LocalAdValues point_residual{};
    point_residual.fill(adlite::Scalar(0.0));
    for (std::size_t node = 0; node < 8; ++node) {
        adlite::Scalar conduction = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            conduction += point.temperature_gradient[node][direction] * temperature_gradient[direction];
        point_residual[node] +=
            point.weighted_measure *
            (conductivity * conduction + point.temperature_shape[node] * heat_capacity * temperature_rate -
                point.temperature_shape[node] * volumetric_heat_source);
    }
    std::array<double, point_width> derivatives{};
    for (std::size_t node = 0; node < 8; ++node) {
        point_residual[node].copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t other = 0; other < 8; ++other) {
            double gradient_dot = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                gradient_dot +=
                    point.temperature_gradient[node][direction] * point.temperature_gradient[other][direction];
            jacobian[node * hex20_local_dof_count + other] +=
                point.weighted_measure * conductivity.value() * gradient_dot +
                point.temperature_shape[other] * derivatives[temperature_index];
        }
    }
    for (std::size_t row = 0; row < residual.size(); ++row) residual[row] += point_residual[row];
}

void add_mechanical_point_system(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time,
    const Hex20LocalValues* committed_state, const CartesianMaterialPointState* committed_material, double time_step,
    Hex20LocalAdValues& residual, Hex20LocalJacobian& jacobian) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    std::array<double, 9> gradient_values{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 20; ++node)
                gradient_values[component * 3 + direction] +=
                    point.displacement_gradient[node][direction] * state[8 + 20 * component + node];
    double temperature_value = 0.0;
    for (std::size_t node = 0; node < 8; ++node) temperature_value += point.temperature_shape[node] * state[node];
    ActiveMatrix3 active_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            active_gradient[component][direction] = adlite::Scalar::independent(
                gradient_values[component * 3 + direction], component * 3 + direction, point_width);
    const adlite::Scalar active_temperature =
        adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    const Hex20LocalValues undeformed{};
    const Hex20LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const Hex20Kinematics kinematics = evaluate_kinematics(point, active_gradient, old_state, strain_formulation);
    const MaterialFunctionContext context = material_context(time, point.position);
    const std::array<const adlite::Scalar*, 6> strain_components = {&kinematics.strain_increment.xx,
        &kinematics.strain_increment.yy, &kinematics.strain_increment.zz, &kinematics.strain_increment.xy,
        &kinematics.strain_increment.yz, &kinematics.strain_increment.xz};
    std::array<double, 6> fed_strain{};
    if (committed_material != nullptr && strain_formulation == StrainFormulation::finite) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node) old_temperature += point.temperature_shape[node] * old_state[node];
        if (!std::isfinite(old_temperature) || !(old_temperature > 0.0))
            throw std::domain_error("Incremental HEX20 material committed temperature must be finite and positive");
        MaterialFunctionContext old_context = context;
        old_context.time -= time_step;
        const SymmetricTensor3 old_imposed = material.eigenstrain(adlite::Scalar(old_temperature), old_context);
        const std::array<double, 6> imposed = {old_imposed.xx.value(), old_imposed.yy.value(), old_imposed.zz.value(),
            old_imposed.xy.value(), old_imposed.yz.value(), old_imposed.xz.value()};
        for (std::size_t component = 0; component < 6; ++component)
            fed_strain[component] = committed_material->elastic_strain[component] +
                                    strain_components[component]->value() + imposed[component] +
                                    committed_material->plastic_strain[component] +
                                    committed_material->creep_strain[component];
    } else {
        for (std::size_t component = 0; component < 6; ++component)
            fed_strain[component] = strain_components[component]->value();
    }
    const CartesianStressTangent tangent =
        evaluate_stress_tangent(material, fed_strain, temperature_value, time_step, committed_material, context);
    std::array<adlite::Scalar, 7> compose_inputs{};
    for (std::size_t component = 0; component < 6; ++component)
        compose_inputs[component] = *strain_components[component];
    compose_inputs[6] = active_temperature;
    const std::array<double, 6> stress_values = {tangent.stress.xx, tangent.stress.yy, tangent.stress.zz,
        tangent.stress.xy, tangent.stress.yz, tangent.stress.xz};
    std::array<double, 7> partials{};
    std::array<adlite::Scalar, 6> composed{};
    for (std::size_t component = 0; component < 6; ++component) {
        for (std::size_t column = 0; column < 6; ++column) partials[column] = tangent.tangent[component][column];
        partials[6] = tangent.thermal[component];
        composed[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    SymmetricTensor3 stress{composed[0], composed[1], composed[2], composed[3], composed[4], composed[5]};
    if (strain_formulation == StrainFormulation::finite) stress = rotate_cartesian_tensor(stress, kinematics.rotation);
    Hex20LocalAdValues point_residual{};
    point_residual.fill(adlite::Scalar(0.0));
    for (std::size_t node = 0; node < 20; ++node) {
        const adlite::Scalar gx = kinematics.current_gradient[node][0], gy = kinematics.current_gradient[node][1],
                             gz = kinematics.current_gradient[node][2];
        point_residual[8 + node] +=
            kinematics.current_weighted_measure * (stress.xx * gx + stress.xy * gy + stress.xz * gz);
        point_residual[28 + node] +=
            kinematics.current_weighted_measure * (stress.xy * gx + stress.yy * gy + stress.yz * gz);
        point_residual[48 + node] +=
            kinematics.current_weighted_measure * (stress.xz * gx + stress.yz * gy + stress.zz * gz);
    }
    std::array<double, point_width> derivatives{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 20; ++node) {
            const std::size_t row = 8 + 20 * component + node;
            point_residual[row].copy_derivatives(derivatives.data(), derivatives.size());
            for (std::size_t other = 0; other < 8; ++other)
                jacobian[row * hex20_local_dof_count + other] +=
                    derivatives[temperature_index] * point.temperature_shape[other];
            for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component)
                for (std::size_t other = 0; other < 20; ++other) {
                    double chained = 0.0;
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        chained += derivatives[displacement_component * 3 + direction] *
                                   point.displacement_gradient[other][direction];
                    jacobian[row * hex20_local_dof_count + 8 + 20 * displacement_component + other] += chained;
                }
        }
    for (std::size_t row = 0; row < residual.size(); ++row) residual[row] += point_residual[row];
}

Hex20LocalResidual compute_local(const CartesianThermoelasticData& data, const Hex20Geometry& geometry,
    const Hex20LocalValues& state, const Hex20LocalValues* committed_state, const CartesianMaterialHistory* history,
    double time_step, Hex20LocalJacobian* jacobian, bool include_thermal_time_term) {
    if (committed_state != nullptr && (!std::isfinite(time_step) || !(time_step > 0.0)))
        throw std::invalid_argument("HEX20 time step must be finite and positive");
    if (history != nullptr && history->size() != geometry.mechanical_points.size())
        throw std::invalid_argument("HEX20 material history must contain 27 integration points");
    Hex20LocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    if (jacobian == nullptr) {
        Hex20LocalAdValues active{};
        ad_local_system::make_passive(state.data(), state.size(), active.data());
        for (const Hex20ThermalQuadraturePoint& point : geometry.thermal_points)
            add_thermal_point_residual(point, active, data.material, data.time, data.volumetric_heat_source,
                committed_state, time_step, residual, include_thermal_time_term);
        for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q)
            add_mechanical_point_residual(geometry.mechanical_points[q], active, data.material, data.strain_formulation,
                data.time, committed_state, history == nullptr ? nullptr : &(*history)[q], time_step, residual);
    } else {
        jacobian->fill(0.0);
        for (const Hex20ThermalQuadraturePoint& point : geometry.thermal_points)
            add_thermal_point_system(point, state, data.material, data.time, data.volumetric_heat_source,
                committed_state, time_step, residual, *jacobian, include_thermal_time_term);
        for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q)
            add_mechanical_point_system(geometry.mechanical_points[q], state, data.material, data.strain_formulation,
                data.time, committed_state, history == nullptr ? nullptr : &(*history)[q], time_step, residual,
                *jacobian);
    }
    Hex20LocalResidual result{};
    ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    return result;
}

void evaluate_quad8_shapes(double xi, double eta, std::array<double, 8>& shape, std::array<double, 8>& derivative_xi,
    std::array<double, 8>& derivative_eta) {
    constexpr std::array<std::array<double, 2>, 4> signs = {
        {{{-1.0, -1.0}}, {{1.0, -1.0}}, {{1.0, 1.0}}, {{-1.0, 1.0}}}};
    for (std::size_t node = 0; node < 4; ++node) {
        const double sx = signs[node][0], sy = signs[node][1], ax = 1.0 + sx * xi, ay = 1.0 + sy * eta,
                     sum = sx * xi + sy * eta - 1.0;
        shape[node] = 0.25 * ax * ay * sum;
        derivative_xi[node] = 0.25 * sx * ay * (sum + ax);
        derivative_eta[node] = 0.25 * sy * ax * (sum + ay);
    }
    shape[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
    derivative_xi[4] = -xi * (1.0 - eta);
    derivative_eta[4] = -0.5 * (1.0 - xi * xi);
    shape[5] = 0.5 * (1.0 + xi) * (1.0 - eta * eta);
    derivative_xi[5] = 0.5 * (1.0 - eta * eta);
    derivative_eta[5] = -(1.0 + xi) * eta;
    shape[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
    derivative_xi[6] = -xi * (1.0 + eta);
    derivative_eta[6] = 0.5 * (1.0 - xi * xi);
    shape[7] = 0.5 * (1.0 - xi) * (1.0 - eta * eta);
    derivative_xi[7] = -0.5 * (1.0 - eta * eta);
    derivative_eta[7] = -(1.0 - xi) * eta;
}
} // namespace

Hex20Geometry make_hex20_geometry(const Hex20Coordinates& coordinates) {
    Hex20Geometry geometry{};
    std::size_t thermal_q = 0;
    for (std::size_t kz = 0; kz < 2; ++kz)
        for (std::size_t ky = 0; ky < 2; ++ky)
            for (std::size_t kx = 0; kx < 2; ++kx) {
                const double xi = gauss2_points[kx], eta = gauss2_points[ky], zeta = gauss2_points[kz];
                std::array<double, 20> displacement_shape{};
                std::array<std::array<double, 3>, 20> displacement_derivative{};
                evaluate_hex20_shapes(xi, eta, zeta, displacement_shape, displacement_derivative);
                std::array<double, 8> temperature_shape{};
                std::array<std::array<double, 3>, 8> temperature_derivative{};
                evaluate_hex8_temperature_shapes(xi, eta, zeta, temperature_shape, temperature_derivative);
                Matrix3 jacobian{};
                CartesianPoint3 position{0.0, 0.0, 0.0};
                for (std::size_t node = 0; node < 20; ++node) {
                    position.x += displacement_shape[node] * coordinates[node].x;
                    position.y += displacement_shape[node] * coordinates[node].y;
                    position.z += displacement_shape[node] * coordinates[node].z;
                    for (std::size_t physical = 0; physical < 3; ++physical) {
                        const double coordinate = physical == 0
                                                      ? coordinates[node].x
                                                      : (physical == 1 ? coordinates[node].y : coordinates[node].z);
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            jacobian[physical][natural] += coordinate * displacement_derivative[node][natural];
                    }
                }
                const double determinant_value = determinant(jacobian);
                if (!std::isfinite(determinant_value) || !(determinant_value > 0.0))
                    throw std::invalid_argument("Hex20Geometry requires a finite positive Jacobian determinant");
                const Matrix3 inverse_jacobian = inverse(jacobian, determinant_value);
                Hex20ThermalQuadraturePoint& point = geometry.thermal_points[thermal_q++];
                point.temperature_shape = temperature_shape;
                point.position = position;
                point.weighted_measure =
                    determinant_value * gauss2_weights[kx] * gauss2_weights[ky] * gauss2_weights[kz];
                for (std::size_t node = 0; node < 8; ++node)
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            point.temperature_gradient[node][physical] +=
                                temperature_derivative[node][natural] * inverse_jacobian[natural][physical];
            }
    std::size_t mechanical_q = 0;
    for (std::size_t kz = 0; kz < 3; ++kz)
        for (std::size_t ky = 0; ky < 3; ++ky)
            for (std::size_t kx = 0; kx < 3; ++kx) {
                const double xi = gauss3_points[kx], eta = gauss3_points[ky], zeta = gauss3_points[kz];
                std::array<double, 20> displacement_shape{};
                std::array<std::array<double, 3>, 20> displacement_derivative{};
                evaluate_hex20_shapes(xi, eta, zeta, displacement_shape, displacement_derivative);
                std::array<double, 8> temperature_shape{};
                std::array<std::array<double, 3>, 8> temperature_derivative{};
                evaluate_hex8_temperature_shapes(xi, eta, zeta, temperature_shape, temperature_derivative);
                Matrix3 jacobian{};
                CartesianPoint3 position{0.0, 0.0, 0.0};
                for (std::size_t node = 0; node < 20; ++node) {
                    position.x += displacement_shape[node] * coordinates[node].x;
                    position.y += displacement_shape[node] * coordinates[node].y;
                    position.z += displacement_shape[node] * coordinates[node].z;
                    for (std::size_t physical = 0; physical < 3; ++physical) {
                        const double coordinate = physical == 0
                                                      ? coordinates[node].x
                                                      : (physical == 1 ? coordinates[node].y : coordinates[node].z);
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            jacobian[physical][natural] += coordinate * displacement_derivative[node][natural];
                    }
                }
                const double determinant_value = determinant(jacobian);
                if (!std::isfinite(determinant_value) || !(determinant_value > 0.0))
                    throw std::invalid_argument("Hex20Geometry requires a finite positive Jacobian determinant");
                const Matrix3 inverse_jacobian = inverse(jacobian, determinant_value);
                Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[mechanical_q++];
                point.temperature_shape = temperature_shape;
                point.displacement_shape = displacement_shape;
                point.position = position;
                point.weighted_measure =
                    determinant_value * gauss3_weights[kx] * gauss3_weights[ky] * gauss3_weights[kz];
                for (std::size_t node = 0; node < 20; ++node)
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            point.displacement_gradient[node][physical] +=
                                displacement_derivative[node][natural] * inverse_jacobian[natural][physical];
            }
    return geometry;
}

Quad8FaceGeometry make_quad8_face_geometry(const Quad8FaceCoordinates& coordinates) {
    Quad8FaceGeometry geometry{};
    std::size_t q = 0;
    for (std::size_t ky = 0; ky < 2; ++ky)
        for (std::size_t kx = 0; kx < 2; ++kx) {
            const double xi = gauss2_points[kx], eta = gauss2_points[ky];
            Quad8FaceThermalQuadraturePoint& point = geometry.thermal_points[q++];
            std::array<double, quad8_face_displacement_node_count> displacement_shape{}, derivative_xi{},
                derivative_eta{};
            evaluate_quad8_shapes(xi, eta, displacement_shape, derivative_xi, derivative_eta);
            point.temperature_shape = {{0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 + eta), 0.25 * (1.0 - xi) * (1.0 + eta)}};
            CartesianPoint3 tangent_xi{}, tangent_eta{};
            for (std::size_t node = 0; node < 8; ++node) {
                tangent_xi.x += derivative_xi[node] * coordinates[node].x;
                tangent_xi.y += derivative_xi[node] * coordinates[node].y;
                tangent_xi.z += derivative_xi[node] * coordinates[node].z;
                tangent_eta.x += derivative_eta[node] * coordinates[node].x;
                tangent_eta.y += derivative_eta[node] * coordinates[node].y;
                tangent_eta.z += derivative_eta[node] * coordinates[node].z;
            }
            const CartesianPoint3 area{tangent_xi.y * tangent_eta.z - tangent_xi.z * tangent_eta.y,
                tangent_xi.z * tangent_eta.x - tangent_xi.x * tangent_eta.z,
                tangent_xi.x * tangent_eta.y - tangent_xi.y * tangent_eta.x};
            const double measure = std::sqrt(area.x * area.x + area.y * area.y + area.z * area.z);
            if (!std::isfinite(measure) || !(measure > 0.0))
                throw std::invalid_argument("Quad8FaceGeometry requires a finite positive area measure");
            point.quadrature_weight = gauss2_weights[kx] * gauss2_weights[ky];
            point.weighted_measure = measure * point.quadrature_weight;
        }
    q = 0;
    for (std::size_t ky = 0; ky < 3; ++ky)
        for (std::size_t kx = 0; kx < 3; ++kx) {
            const double xi = gauss3_points[kx], eta = gauss3_points[ky];
            Quad8FaceMechanicalQuadraturePoint& point = geometry.mechanical_points[q++];
            evaluate_quad8_shapes(xi, eta, point.displacement_shape, point.derivative_xi, point.derivative_eta);
            for (std::size_t node = 0; node < 8; ++node) {
                point.tangent_xi.x += point.derivative_xi[node] * coordinates[node].x;
                point.tangent_xi.y += point.derivative_xi[node] * coordinates[node].y;
                point.tangent_xi.z += point.derivative_xi[node] * coordinates[node].z;
                point.tangent_eta.x += point.derivative_eta[node] * coordinates[node].x;
                point.tangent_eta.y += point.derivative_eta[node] * coordinates[node].y;
                point.tangent_eta.z += point.derivative_eta[node] * coordinates[node].z;
            }
            const CartesianPoint3 area{
                point.tangent_xi.y * point.tangent_eta.z - point.tangent_xi.z * point.tangent_eta.y,
                point.tangent_xi.z * point.tangent_eta.x - point.tangent_xi.x * point.tangent_eta.z,
                point.tangent_xi.x * point.tangent_eta.y - point.tangent_xi.y * point.tangent_eta.x};
            const double measure = std::sqrt(area.x * area.x + area.y * area.y + area.z * area.z);
            if (!std::isfinite(measure) || !(measure > 0.0))
                throw std::invalid_argument("Quad8FaceGeometry requires a finite positive area measure");
            point.quadrature_weight = gauss3_weights[kx] * gauss3_weights[ky];
        }
    return geometry;
}

void validate_hex20_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    const double value = determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX20 deformation must preserve a positive Jacobian");
}

Hex20LocalResidual compute_hex20_thermoelastic(const CartesianThermoelasticData& data, const Hex20Geometry& geometry,
    const Hex20LocalValues& state, const Hex20LocalValues* committed_state, double time_step,
    Hex20LocalJacobian* jacobian) {
    return compute_local(data, geometry, state, committed_state, nullptr, time_step, jacobian, true);
}

Hex20LocalResidual compute_hex20_transient(const CartesianThermoelasticData& data, const Hex20Geometry& geometry,
    const Hex20LocalValues& state, const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material, double time_step, Hex20LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    return compute_local(
        data, geometry, state, &committed_state, &committed_material, time_step, jacobian, include_thermal_time_term);
}

CartesianMaterialHistory compute_hex20_transient_update(const CartesianThermoelasticData& data,
    const Hex20Geometry& geometry, const Hex20LocalValues& state, const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material, double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("HEX20 transient update time step must be finite and positive");
    if (committed_material.size() != geometry.mechanical_points.size())
        throw std::invalid_argument("HEX20 material history must contain 27 integration points");
    Hex20LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    CartesianMaterialHistory result(geometry.mechanical_points.size());
    for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q) {
        const Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[q];
        const adlite::Scalar temperature = interpolate_temperature(point, passive);
        const Hex20Kinematics kinematics =
            evaluate_kinematics(point, passive, committed_state, data.strain_formulation);
        CartesianInelasticStressResponse response;
        if (data.strain_formulation == StrainFormulation::finite) {
            double old_temperature = 0.0;
            for (std::size_t node = 0; node < 8; ++node)
                old_temperature += point.temperature_shape[node] * committed_state[node];
            response = data.material.incremental_response(kinematics.strain_increment, kinematics.rotation, temperature,
                old_temperature, time_step, committed_material[q], material_context(data.time, point.position));
        } else {
            response = data.material.response(kinematics.strain_increment, temperature, time_step,
                committed_material[q], material_context(data.time, point.position));
        }
        result[q] = response.trial_state;
    }
    return result;
}

std::array<SymmetricTensor3Values, hex20_mechanical_quadrature_point_count> compute_hex20_stress(
    const CartesianThermoelasticData& data, const Hex20Geometry& geometry, const Hex20LocalValues& state) {
    Hex20LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    std::array<SymmetricTensor3Values, 27> result{};
    for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q) {
        const Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[q];
        const adlite::Scalar temperature = interpolate_temperature(point, passive);
        const Hex20Kinematics kinematics =
            evaluate_kinematics(point, passive, Hex20LocalValues{}, data.strain_formulation);
        SymmetricTensor3 stress =
            data.material.stress(kinematics.strain_increment, temperature, material_context(data.time, point.position));
        if (data.strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
        result[q] = {stress.xx.value(), stress.yy.value(), stress.zz.value(), stress.xy.value(), stress.yz.value(),
            stress.xz.value()};
    }
    return result;
}

Quad8FaceLocalResidual compute_quad8_face_boundary(const Quad4FaceBoundaryData& data, const Quad8FaceGeometry& geometry,
    const Quad8FaceLocalValues& state, Quad8FaceLocalJacobian* jacobian) {
    Quad8FaceLocalAdValues ad_state{};
    if (jacobian == nullptr)
        ad_local_system::make_passive(state.data(), state.size(), ad_state.data());
    else
        ad_local_system::make_active(state.data(), state.size(), ad_state.data());
    Quad8FaceLocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    if (data.kind == Quad4FaceBoundaryKind::convection) {
        for (const Quad8FaceThermalQuadraturePoint& point : geometry.thermal_points) {
            adlite::Scalar temperature = 0.0;
            for (std::size_t node = 0; node < 4; ++node) temperature += point.temperature_shape[node] * ad_state[node];
            const adlite::Scalar heat_flux = data.load * (temperature - data.ambient_temperature);
            for (std::size_t node = 0; node < 4; ++node)
                residual[node] += point.weighted_measure * point.temperature_shape[node] * heat_flux;
        }
    } else {
        for (const Quad8FaceMechanicalQuadraturePoint& point : geometry.mechanical_points) {
            std::array<adlite::Scalar, 3> tangent_xi = {point.tangent_xi.x, point.tangent_xi.y, point.tangent_xi.z};
            std::array<adlite::Scalar, 3> tangent_eta = {point.tangent_eta.x, point.tangent_eta.y, point.tangent_eta.z};
            if (data.use_displaced_geometry)
                for (std::size_t node = 0; node < 8; ++node)
                    for (std::size_t component = 0; component < 3; ++component) {
                        tangent_xi[component] += point.derivative_xi[node] * ad_state[4 + 8 * component + node];
                        tangent_eta[component] += point.derivative_eta[node] * ad_state[4 + 8 * component + node];
                    }
            const std::array<adlite::Scalar, 3> area = {tangent_xi[1] * tangent_eta[2] - tangent_xi[2] * tangent_eta[1],
                tangent_xi[2] * tangent_eta[0] - tangent_xi[0] * tangent_eta[2],
                tangent_xi[0] * tangent_eta[1] - tangent_xi[1] * tangent_eta[0]};
            const adlite::Scalar measure = adlite::hypot(adlite::hypot(area[0], area[1]), area[2]);
            if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
                throw std::domain_error("Three-dimensional quadratic face requires a positive current measure");
            if (data.kind == Quad4FaceBoundaryKind::pressure) {
                for (std::size_t node = 0; node < 8; ++node) {
                    residual[4 + node] +=
                        data.load * point.displacement_shape[node] * area[0] * point.quadrature_weight;
                    residual[12 + node] +=
                        data.load * point.displacement_shape[node] * area[1] * point.quadrature_weight;
                    residual[20 + node] +=
                        data.load * point.displacement_shape[node] * area[2] * point.quadrature_weight;
                }
            } else if (data.kind == Quad4FaceBoundaryKind::traction) {
                const std::size_t offset = data.component == CartesianTractionComponent::x
                                               ? 4
                                               : (data.component == CartesianTractionComponent::y ? 12 : 20);
                for (std::size_t node = 0; node < 8; ++node)
                    residual[offset + node] -=
                        data.load * point.displacement_shape[node] * measure * point.quadrature_weight;
            }
        }
    }
    Quad8FaceLocalResidual values{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), values.data());
    else
        ad_local_system::extract_system(residual.data(), ad_state.size(), values.data(), jacobian->data());
    return values;
}
} // namespace fuelsim
