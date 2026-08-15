#include "fuelsim/hex8.hpp"
#include "ad_local_system.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
constexpr double gauss = 0.577350269189625764509148780501957456;
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {
    {{{-1.0, -1.0, -1.0}}, {{1.0, -1.0, -1.0}}, {{1.0, 1.0, -1.0}}, {{-1.0, 1.0, -1.0}}, {{-1.0, -1.0, 1.0}},
        {{1.0, -1.0, 1.0}}, {{1.0, 1.0, 1.0}}, {{-1.0, 1.0, 1.0}}}};

double determinant(const std::array<std::array<double, 3>, 3>& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

std::array<std::array<double, 3>, 3> inverse(
    const std::array<std::array<double, 3>, 3>& matrix, double determinant_value) {
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

using ActiveMatrix3 = std::array<std::array<adlite::Scalar, 3>, 3>;

adlite::Scalar determinant(const ActiveMatrix3& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
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

ActiveMatrix3 multiply(const ActiveMatrix3& first, const std::array<std::array<double, 3>, 3>& second) {
    ActiveMatrix3 result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k) result[i][j] += first[i][k] * second[k][j];
    return result;
}

ActiveMatrix3 displacement_gradient(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state) {
    ActiveMatrix3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                result[component][direction] += point.gradient[node][direction] * state[8 * (component + 1) + node];
    return result;
}

std::array<std::array<double, 3>, 3> deformation_gradient(
    const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    std::array<std::array<double, 3>, 3> result{};
    for (std::size_t component = 0; component < 3; ++component) {
        result[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                result[component][direction] += point.gradient[node][direction] * state[8 * (component + 1) + node];
    }
    return result;
}

adlite::Scalar interpolate_hex8(
    const std::array<double, 8>& coefficients, const Hex8LocalAdValues& state, std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node) result += coefficients[node] * state[offset + node];
    return result;
}

SymmetricTensor3 strain_at(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state) {
    adlite::Scalar ux_x = 0.0, ux_y = 0.0, ux_z = 0.0, uy_x = 0.0, uy_y = 0.0, uy_z = 0.0, uz_x = 0.0, uz_y = 0.0,
                   uz_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        ux_x += point.gradient[node][0] * state[8 + node];
        ux_y += point.gradient[node][1] * state[8 + node];
        ux_z += point.gradient[node][2] * state[8 + node];
        uy_x += point.gradient[node][0] * state[16 + node];
        uy_y += point.gradient[node][1] * state[16 + node];
        uy_z += point.gradient[node][2] * state[16 + node];
        uz_x += point.gradient[node][0] * state[24 + node];
        uz_y += point.gradient[node][1] * state[24 + node];
        uz_z += point.gradient[node][2] * state[24 + node];
    }
    return {ux_x, uy_y, uz_z, 0.5 * (ux_y + uy_x), 0.5 * (uy_z + uz_y), 0.5 * (ux_z + uz_x)};
}
} // namespace

SymmetricTensor3 rotate_cartesian_tensor(const SymmetricTensor3& tensor, const CartesianRotation& rotation) {
    const ActiveMatrix3 r = {{{rotation.xx, rotation.xy, rotation.xz}, {rotation.yx, rotation.yy, rotation.yz},
        {rotation.zx, rotation.zy, rotation.zz}}};
    const ActiveMatrix3 value = {
        {{tensor.xx, tensor.xy, tensor.xz}, {tensor.xy, tensor.yy, tensor.yz}, {tensor.xz, tensor.yz, tensor.zz}}};
    ActiveMatrix3 rotated{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                for (std::size_t l = 0; l < 3; ++l) rotated[i][j] += r[i][k] * value[k][l] * r[j][l];
    return {rotated[0][0], rotated[1][1], rotated[2][2], rotated[0][1], rotated[1][2], rotated[0][2]};
}

CartesianInelasticStressResponse IsotropicThermoelasticMaterial::incremental_response(
    const SymmetricTensor3& strain_increment, const CartesianRotation& rotation, const adlite::Scalar& temperature,
    double committed_temperature, double time_step, const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental Cartesian material committed temperature must be finite and positive");
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const SymmetricTensor3 old_imposed = eigenstrain(adlite::Scalar(committed_temperature), old_context);
    const SymmetricTensor3 synthetic_total{committed.elastic_strain[0] + strain_increment.xx + old_imposed.xx +
                                               committed.plastic_strain[0] + committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment.yy + old_imposed.yy + committed.plastic_strain[1] +
            committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment.zz + old_imposed.zz + committed.plastic_strain[2] +
            committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment.xy + old_imposed.xy + committed.plastic_strain[3] +
            committed.creep_strain[3],
        committed.elastic_strain[4] + strain_increment.yz + old_imposed.yz + committed.plastic_strain[4] +
            committed.creep_strain[4],
        committed.elastic_strain[5] + strain_increment.xz + old_imposed.xz + committed.plastic_strain[5] +
            committed.creep_strain[5]};
    CartesianInelasticStressResponse result = response(synthetic_total, temperature, time_step, committed, context);
    result.stress = rotate_cartesian_tensor(result.stress, rotation);
    std::array<double, 6>* histories[3] = {
        &result.trial_state.elastic_strain, &result.trial_state.plastic_strain, &result.trial_state.creep_strain};
    for (std::array<double, 6>* history : histories) {
        const SymmetricTensor3 rotated = rotate_cartesian_tensor(
            {(*history)[0], (*history)[1], (*history)[2], (*history)[3], (*history)[4], (*history)[5]}, rotation);
        *history = {rotated.xx.value(), rotated.yy.value(), rotated.zz.value(), rotated.xy.value(), rotated.yz.value(),
            rotated.xz.value()};
        for (double component : *history)
            if (!std::isfinite(component))
                throw std::domain_error("Rotated Cartesian material history must contain only finite values");
    }
    result.trial_state.stress = {result.stress.xx.value(), result.stress.yy.value(), result.stress.zz.value(),
        result.stress.xy.value(), result.stress.yz.value(), result.stress.xz.value()};
    return result;
}

void validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    const double value = determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX8 deformation must preserve a positive Jacobian");
}

CartesianKinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state, const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    CartesianKinematics result{};
    if (strain_formulation == StrainFormulation::small) {
        result.strain_increment = strain_at(point, current_state);
        for (std::size_t node = 0; node < 8; ++node)
            for (std::size_t direction = 0; direction < 3; ++direction)
                result.current_gradient[node][direction] = point.gradient[node][direction];
        result.current_weighted_measure = point.weighted_measure;
        return result;
    }
    ActiveMatrix3 current = displacement_gradient(point, current_state);
    for (std::size_t direction = 0; direction < 3; ++direction) current[direction][direction] += 1.0;
    const adlite::Scalar current_determinant = determinant(current);
    if (!std::isfinite(current_determinant.value()) || !(current_determinant.value() > 0.0))
        throw std::domain_error("Finite-strain HEX8 deformation must preserve a positive Jacobian");
    const ActiveMatrix3 current_inverse = inverse(current, current_determinant);
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.current_gradient[node][direction] +=
                    point.gradient[node][reference] * current_inverse[reference][direction];
    result.current_weighted_measure = point.weighted_measure * current_determinant;
    const std::array<std::array<double, 3>, 3> old = deformation_gradient(point, committed_state);
    const double old_determinant = determinant(old);
    if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
        throw std::domain_error("Committed finite-strain HEX8 state requires a positive Jacobian");
    const ActiveMatrix3 incremental = multiply(current, inverse(old, old_determinant));
    const adlite::Scalar incremental_determinant = determinant(incremental);
    if (!std::isfinite(incremental_determinant.value()) || !(incremental_determinant.value() > 0.0))
        throw std::domain_error("Incremental finite-strain HEX8 state requires a positive Jacobian");
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
        throw std::domain_error("MOOSE Taylor finite-strain rotation has invalid three-dimensional p+q");
    const adlite::Scalar p2 = p * p, p3 = p2 * p, p4 = p3 * p, sum2 = sum * sum, sum3 = sum2 * sum;
    const adlite::Scalar c1_squared = p + 3.0 * p2 * (1.0 - sum) / sum2 - 2.0 * p3 * (1.0 - sum) / sum3;
    if (!std::isfinite(c1_squared.value()) || !(c1_squared.value() > 0.0))
        throw std::domain_error("MOOSE three-dimensional Rashid rotation has nonpositive C1 squared");
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
        throw std::domain_error("MOOSE three-dimensional Rashid rotation has nonpositive C3 test");
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

namespace {
MaterialFunctionContext material_context(double time, const CartesianPoint3& point) {
    return {time, point.x, point.y, point.z};
}

void add_hex8_point_residual(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time,
    double volumetric_heat_source, const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material, double time_step, Hex8LocalAdValues& residual) {
    const adlite::Scalar temperature = interpolate_hex8(point.shape, state, 0);
    adlite::Scalar gradient_temperature_x = 0.0, gradient_temperature_y = 0.0, gradient_temperature_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        gradient_temperature_x += point.gradient[node][0] * state[node];
        gradient_temperature_y += point.gradient[node][1] * state[node];
        gradient_temperature_z += point.gradient[node][2] * state[node];
    }
    const MaterialFunctionContext context = material_context(time, point.position);
    const adlite::Scalar conductivity = material.conductivity(temperature, context);
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const CartesianKinematics kinematics =
        evaluate_cartesian_incremental_kinematics(point, state, old_state, strain_formulation);
    SymmetricTensor3 stress;
    if (committed_material == nullptr) {
        stress = material.stress(kinematics.strain_increment, temperature, context);
        if (strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
    } else if (strain_formulation == StrainFormulation::finite) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node) old_temperature += point.shape[node] * old_state[node];
        stress = material
                     .incremental_response(kinematics.strain_increment, kinematics.rotation, temperature,
                         old_temperature, time_step, *committed_material, context)
                     .stress;
    } else {
        stress =
            material.response(kinematics.strain_increment, temperature, time_step, *committed_material, context).stress;
    }
    adlite::Scalar temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node) old_temperature += point.shape[node] * (*committed_state)[node];
        temperature_rate = (temperature - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(temperature, context);
    }
    for (std::size_t node = 0; node < 8; ++node) {
        const double gradient_x = point.gradient[node][0], gradient_y = point.gradient[node][1],
                     gradient_z = point.gradient[node][2];
        residual[node] +=
            point.weighted_measure *
            (conductivity * (gradient_x * gradient_temperature_x + gradient_y * gradient_temperature_y +
                                gradient_z * gradient_temperature_z) +
                point.shape[node] * heat_capacity * temperature_rate - point.shape[node] * volumetric_heat_source);
        const adlite::Scalar current_gradient_x = kinematics.current_gradient[node][0];
        const adlite::Scalar current_gradient_y = kinematics.current_gradient[node][1];
        const adlite::Scalar current_gradient_z = kinematics.current_gradient[node][2];
        residual[8 + node] +=
            kinematics.current_weighted_measure *
            (stress.xx * current_gradient_x + stress.xy * current_gradient_y + stress.xz * current_gradient_z);
        residual[16 + node] +=
            kinematics.current_weighted_measure *
            (stress.xy * current_gradient_x + stress.yy * current_gradient_y + stress.yz * current_gradient_z);
        residual[24 + node] +=
            kinematics.current_weighted_measure *
            (stress.xz * current_gradient_x + stress.yz * current_gradient_y + stress.zz * current_gradient_z);
    }
}

std::array<SymmetricTensor3Values, 8> evaluate_hex8_stress(const Hex8Geometry& geometry, const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time) {
    Hex8LocalAdValues ad_state{};
    ad_local_system::make_passive(state.data(), state.size(), ad_state.data());
    std::array<SymmetricTensor3Values, 8> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = interpolate_hex8(point.shape, ad_state, 0);
        const CartesianKinematics kinematics =
            evaluate_cartesian_incremental_kinematics(point, ad_state, Hex8LocalValues{}, strain_formulation);
        SymmetricTensor3 stress =
            material.stress(kinematics.strain_increment, temperature, material_context(time, point.position));
        if (strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
        result[q] = {stress.xx.value(), stress.yy.value(), stress.zz.value(), stress.xy.value(), stress.yz.value(),
            stress.xz.value()};
    }
    return result;
}

Hex8LocalResidual compute_hex8_local(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state, const Hex8MaterialHistory* history,
    double time_step, Hex8LocalJacobian* jacobian) {
    if (committed_state != nullptr && (!std::isfinite(time_step) || time_step <= 0.0))
        throw std::invalid_argument("HEX8 time step must be finite and positive");
    Hex8LocalAdValues active{}, residual{};
    if (jacobian == nullptr)
        ad_local_system::make_passive(state.data(), state.size(), active.data());
    else
        ad_local_system::make_active(state.data(), state.size(), active.data());
    residual.fill(adlite::Scalar(0.0));
    for (std::size_t q = 0; q < geometry.points.size(); ++q)
        add_hex8_point_residual(geometry.points[q], active, data.material, data.strain_formulation, data.time,
            data.volumetric_heat_source, committed_state, history == nullptr ? nullptr : &(*history)[q], time_step,
            residual);
    Hex8LocalResidual result{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    else
        ad_local_system::extract_system(residual.data(), active.size(), result.data(), jacobian->data());
    return result;
}
} // namespace

Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates) {
    Hex8Geometry geometry{};
    std::size_t q = 0;
    for (double zeta : {-gauss, gauss}) {
        for (double eta : {-gauss, gauss}) {
            for (double xi : {-gauss, gauss}) {
                std::array<double, 8> shape{};
                std::array<std::array<double, 3>, 8> derivative{};
                for (std::size_t node = 0; node < 8; ++node) {
                    const double sx = hex8_signs[node][0], sy = hex8_signs[node][1], sz = hex8_signs[node][2];
                    shape[node] = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta);
                    derivative[node] = {{0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
                        0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta),
                        0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)}};
                }
                std::array<std::array<double, 3>, 3> jacobian{};
                CartesianPoint3 position{0.0, 0.0, 0.0};
                for (std::size_t node = 0; node < 8; ++node) {
                    const std::array<double, 3> coordinate = {
                        coordinates[node].x,
                        coordinates[node].y,
                        coordinates[node].z,
                    };
                    position.x += shape[node] * coordinates[node].x;
                    position.y += shape[node] * coordinates[node].y;
                    position.z += shape[node] * coordinates[node].z;
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            jacobian[physical][natural] += coordinate[physical] * derivative[node][natural];
                }
                const double determinant_value = determinant(jacobian);
                if (!std::isfinite(determinant_value) || !(determinant_value > 0.0))
                    throw std::invalid_argument("Hex8Geometry requires a finite positive Jacobian determinant");
                const std::array<std::array<double, 3>, 3> inverse_jacobian = inverse(jacobian, determinant_value);
                Hex8QuadraturePoint& point = geometry.points[q++];
                point.shape = shape;
                point.position = position;
                point.weighted_measure = determinant_value;
                for (std::size_t node = 0; node < 8; ++node) {
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            point.gradient[node][physical] +=
                                derivative[node][natural] * inverse_jacobian[natural][physical];
                }
            }
        }
    }
    return geometry;
}

Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates) {
    Quad4FaceGeometry geometry{};
    const std::array<std::array<double, 2>, 4> locations = {
        {{{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}}}};
    for (std::size_t q = 0; q < locations.size(); ++q) {
        const double xi = locations[q][0], eta = locations[q][1];
        const std::array<double, 4> shape = {{0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta), 0.25 * (1.0 - xi) * (1.0 + eta)}};
        const std::array<double, 4> derivative_xi = {
            {-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)}};
        const std::array<double, 4> derivative_eta = {
            {-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)}};
        CartesianPoint3 tangent_xi{0.0, 0.0, 0.0};
        CartesianPoint3 tangent_eta{0.0, 0.0, 0.0};
        for (std::size_t node = 0; node < 4; ++node) {
            tangent_xi.x += derivative_xi[node] * coordinates[node].x;
            tangent_xi.y += derivative_xi[node] * coordinates[node].y;
            tangent_xi.z += derivative_xi[node] * coordinates[node].z;
            tangent_eta.x += derivative_eta[node] * coordinates[node].x;
            tangent_eta.y += derivative_eta[node] * coordinates[node].y;
            tangent_eta.z += derivative_eta[node] * coordinates[node].z;
        }
        const CartesianPoint3 area_vector{tangent_xi.y * tangent_eta.z - tangent_xi.z * tangent_eta.y,
            tangent_xi.z * tangent_eta.x - tangent_xi.x * tangent_eta.z,
            tangent_xi.x * tangent_eta.y - tangent_xi.y * tangent_eta.x};
        const double measure =
            std::sqrt(area_vector.x * area_vector.x + area_vector.y * area_vector.y + area_vector.z * area_vector.z);
        if (!std::isfinite(measure) || !(measure > 0.0))
            throw std::invalid_argument("Quad4FaceGeometry requires a finite positive area measure");
        geometry.points[q] = {shape, derivative_xi, derivative_eta, tangent_xi, tangent_eta, measure};
    }
    return geometry;
}

Hex8LocalResidual compute_hex8_thermoelastic(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state, double time_step,
    Hex8LocalJacobian* jacobian) {
    return compute_hex8_local(data, geometry, state, committed_state, nullptr, time_step, jacobian);
}

Hex8LocalResidual compute_hex8_transient(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material,
    double time_step, Hex8LocalJacobian* jacobian) {
    return compute_hex8_local(data, geometry, state, &committed_state, &committed_material, time_step, jacobian);
}

Hex8MaterialHistory compute_hex8_transient_update(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material,
    double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("HEX8 transient update time step must be finite and positive");
    Hex8LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    Hex8MaterialHistory result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = interpolate_hex8(point.shape, passive, 0);
        const CartesianKinematics kinematics =
            evaluate_cartesian_incremental_kinematics(point, passive, committed_state, data.strain_formulation);
        CartesianInelasticStressResponse response;
        if (data.strain_formulation == StrainFormulation::finite) {
            double old_temperature = 0.0;
            for (std::size_t node = 0; node < 8; ++node) old_temperature += point.shape[node] * committed_state[node];
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

std::array<SymmetricTensor3Values, 8> compute_hex8_stress(
    const Hex8ThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    return evaluate_hex8_stress(geometry, state, data.material, data.strain_formulation, data.time);
}

Quad4FaceLocalResidual compute_quad4_face_boundary(const Quad4FaceBoundaryData& data, const Quad4FaceGeometry& geometry,
    const Quad4FaceLocalValues& state, Quad4FaceLocalJacobian* jacobian) {
    Quad4FaceLocalAdValues ad_state{};
    if (jacobian == nullptr)
        ad_local_system::make_passive(state.data(), state.size(), ad_state.data());
    else
        ad_local_system::make_active(state.data(), state.size(), ad_state.data());
    Quad4FaceLocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        std::array<adlite::Scalar, 3> tangent_xi = {point.tangent_xi.x, point.tangent_xi.y, point.tangent_xi.z};
        std::array<adlite::Scalar, 3> tangent_eta = {point.tangent_eta.x, point.tangent_eta.y, point.tangent_eta.z};
        if (data.use_displaced_geometry)
            for (std::size_t node = 0; node < 4; ++node)
                for (std::size_t component = 0; component < 3; ++component) {
                    tangent_xi[component] += point.derivative_xi[node] * ad_state[4 * (component + 1) + node];
                    tangent_eta[component] += point.derivative_eta[node] * ad_state[4 * (component + 1) + node];
                }
        const std::array<adlite::Scalar, 3> area = {tangent_xi[1] * tangent_eta[2] - tangent_xi[2] * tangent_eta[1],
            tangent_xi[2] * tangent_eta[0] - tangent_xi[0] * tangent_eta[2],
            tangent_xi[0] * tangent_eta[1] - tangent_xi[1] * tangent_eta[0]};
        const adlite::Scalar measure = adlite::hypot(adlite::hypot(area[0], area[1]), area[2]);
        if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
            throw std::domain_error("Three-dimensional mechanical face requires a positive current measure");
        if (data.kind == Quad4FaceBoundaryKind::pressure) {
            for (std::size_t node = 0; node < 4; ++node) {
                residual[4 + node] += data.load * point.shape[node] * area[0];
                residual[8 + node] += data.load * point.shape[node] * area[1];
                residual[12 + node] += data.load * point.shape[node] * area[2];
            }
        } else if (data.kind == Quad4FaceBoundaryKind::traction) {
            const std::size_t offset = data.component == CartesianTractionComponent::x
                                           ? 4
                                           : (data.component == CartesianTractionComponent::y ? 8 : 12);
            for (std::size_t node = 0; node < 4; ++node)
                residual[offset + node] -= data.load * measure * point.shape[node];
        } else {
            adlite::Scalar temperature = 0.0;
            for (std::size_t node = 0; node < 4; ++node) temperature += point.shape[node] * ad_state[node];
            const adlite::Scalar heat_flux = data.load * (temperature - data.ambient_temperature);
            for (std::size_t node = 0; node < 4; ++node)
                residual[node] += point.weighted_measure * point.shape[node] * heat_flux;
        }
    }
    Quad4FaceLocalResidual values{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), values.data());
    else
        ad_local_system::extract_system(residual.data(), ad_state.size(), values.data(), jacobian->data());
    return values;
}
} // namespace fuelsim
