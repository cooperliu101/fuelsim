#include "detail/ad_local_system.hpp"
#include "detail/cartesian3d_mechanics.hpp"
#include "fuelsim/core/cartesian3d_hex8.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {
constexpr double gauss = 0.577350269189625764509148780501957456;
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {
    {{{-1.0, -1.0, -1.0}}, {{1.0, -1.0, -1.0}}, {{1.0, 1.0, -1.0}}, {{-1.0, 1.0, -1.0}}, {{-1.0, -1.0, 1.0}},
        {{1.0, -1.0, 1.0}}, {{1.0, 1.0, 1.0}}, {{-1.0, 1.0, 1.0}}}};
constexpr std::array<std::size_t, hex8_node_count> hex8_node_to_gauss = {0, 1, 3, 2, 4, 5, 7, 6};
constexpr std::array<std::array<double, 4>, hex8_node_count> finite_reduced_hex8_raw_hourglass = {
    {{{1.0, -1.0, 1.0, -1.0}}, {{-1.0, -1.0, -1.0, 1.0}}, {{1.0, 1.0, -1.0, -1.0}}, {{-1.0, 1.0, 1.0, 1.0}},
        {{1.0, 1.0, -1.0, 1.0}}, {{-1.0, 1.0, 1.0, -1.0}}, {{1.0, -1.0, 1.0, 1.0}}, {{-1.0, -1.0, -1.0, -1.0}}}};

cartesian_detail::ActiveMatrix3 displacement_gradient(
    const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state) {
    cartesian_detail::ActiveMatrix3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                result[component][direction] += point.gradient[node][direction] * state[8 * (component + 1) + node];
    return result;
}

cartesian_detail::Matrix3 deformation_gradient(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    cartesian_detail::Matrix3 result{};
    for (std::size_t component = 0; component < 3; ++component) {
        result[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                result[component][direction] += point.gradient[node][direction] * state[8 * (component + 1) + node];
    }
    return result;
}

adlite::Scalar average_hex8_temperature(const Hex8LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node) result += state[node] / 8.0;
    return result;
}

double average_hex8_temperature(const Hex8LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < 8; ++node) result += state[node] / 8.0;
    return result;
}

adlite::Scalar reduced_hex8_temperature(const Hex8Geometry& geometry, const Hex8LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        result += geometry.reduced_capacity_points[node].weighted_measure / geometry.reference_volume * state[node];
    return result;
}

double reduced_hex8_temperature(const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        result += geometry.reduced_capacity_points[node].weighted_measure / geometry.reference_volume * state[node];
    return result;
}

double average_hex8_strain_trace(const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result += geometry.average_shape_gradient[node][component] * state[8 * (component + 1) + node];
    return result;
}

struct FiniteTracePoint final {
    adlite::Scalar strain_trace;
    adlite::Scalar midpoint_weighted_measure;
    adlite::Scalar current_weighted_measure;
};

FiniteTracePoint finite_trace_point(const Hex8QuadraturePoint& point, const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state) {
    const cartesian_detail::Matrix3 old = deformation_gradient(point, committed_state);
    const cartesian_detail::KinematicsCore kinematics =
        cartesian_detail::evaluate_kinematics(gradient, old, StrainFormulation::finite);
    cartesian_detail::ActiveMatrix3 midpoint{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            midpoint[i][j] = 0.5 * (gradient[i][j] + old[i][j]);
            if (i == j) midpoint[i][j] += 0.5;
        }
    const adlite::Scalar midpoint_determinant = cartesian_detail::determinant(midpoint);
    if (!std::isfinite(midpoint_determinant.value()) || !(midpoint_determinant.value() > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint configuration must preserve a positive Jacobian");
    return {kinematics.strain_increment.xx + kinematics.strain_increment.yy + kinematics.strain_increment.zz,
        point.weighted_measure * midpoint_determinant, point.weighted_measure * kinematics.current_determinant};
}

double finite_average_hex8_strain_trace(
    const Hex8Geometry& geometry, const Hex8LocalValues& state, const Hex8LocalValues& committed_state) {
    Hex8LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    double numerator = 0.0, volume = 0.0;
    for (const Hex8QuadraturePoint& point : geometry.points) {
        const FiniteTracePoint trace_point =
            finite_trace_point(point, displacement_gradient(point, passive), committed_state);
        numerator += trace_point.midpoint_weighted_measure.value() * trace_point.strain_trace.value();
        volume += trace_point.midpoint_weighted_measure.value();
    }
    if (!std::isfinite(volume) || !(volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint volume must be finite and positive");
    return numerator / volume;
}

struct FiniteAverageTraceSystem final {
    double value = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> displacement_derivatives{};
    double current_volume = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> current_volume_derivatives{};
};

FiniteAverageTraceSystem finite_average_hex8_strain_trace_system(
    const Hex8Geometry& geometry, const Hex8LocalValues& state, const Hex8LocalValues& committed_state) {
    constexpr std::size_t width = 9;
    FiniteAverageTraceSystem result;
    double numerator = 0.0, volume = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> numerator_derivatives{}, volume_derivatives{};
    for (const Hex8QuadraturePoint& point : geometry.points) {
        cartesian_detail::ActiveMatrix3 gradient{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                double value = 0.0;
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    value += point.gradient[node][direction] * state[8 * (component + 1) + node];
                gradient[component][direction] = adlite::Scalar::independent(value, component * 3 + direction, width);
            }
        const FiniteTracePoint trace_point = finite_trace_point(point, gradient, committed_state);
        const adlite::Scalar point_numerator = trace_point.midpoint_weighted_measure * trace_point.strain_trace;
        numerator += point_numerator.value();
        volume += trace_point.midpoint_weighted_measure.value();
        result.current_volume += trace_point.current_weighted_measure.value();
        std::array<double, width> point_numerator_derivatives{}, point_volume_derivatives{},
            point_current_volume_derivatives{};
        point_numerator.copy_derivatives(point_numerator_derivatives.data(), point_numerator_derivatives.size());
        trace_point.midpoint_weighted_measure.copy_derivatives(
            point_volume_derivatives.data(), point_volume_derivatives.size());
        trace_point.current_weighted_measure.copy_derivatives(
            point_current_volume_derivatives.data(), point_current_volume_derivatives.size());
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t direction = 0; direction < 3; ++direction) {
                    const double shape_chain = point.gradient[node][direction];
                    numerator_derivatives[node][component] +=
                        point_numerator_derivatives[component * 3 + direction] * shape_chain;
                    volume_derivatives[node][component] +=
                        point_volume_derivatives[component * 3 + direction] * shape_chain;
                    result.current_volume_derivatives[node][component] +=
                        point_current_volume_derivatives[component * 3 + direction] * shape_chain;
                }
    }
    if (!std::isfinite(volume) || !(volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint volume must be finite and positive");
    if (!std::isfinite(result.current_volume) || !(result.current_volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 current volume must be finite and positive");
    result.value = numerator / volume;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result.displacement_derivatives[node][component] =
                (numerator_derivatives[node][component] - result.value * volume_derivatives[node][component]) / volume;
    return result;
}

SymmetricTensor3 selectively_reduced_strain(
    const SymmetricTensor3& strain, double average_trace, StrainFormulation strain_formulation) {
    if (strain_formulation != StrainFormulation::small && strain_formulation != StrainFormulation::finite)
        return strain;
    const adlite::Scalar correction = (average_trace - strain.xx - strain.yy - strain.zz) / 3.0;
    return {strain.xx + correction, strain.yy + correction, strain.zz + correction, strain.xy, strain.yz, strain.xz};
}

SymmetricTensor3 expansion_adjusted_strain(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain, const adlite::Scalar& point_temperature, const adlite::Scalar& element_temperature,
    MaterialFunctionContext context) {
    const SymmetricTensor3 point_imposed = material.eigenstrain(point_temperature, context);
    const SymmetricTensor3 element_imposed = material.eigenstrain(element_temperature, context);
    return {strain.xx + point_imposed.xx - element_imposed.xx, strain.yy + point_imposed.yy - element_imposed.yy,
        strain.zz + point_imposed.zz - element_imposed.zz, strain.xy + point_imposed.xy - element_imposed.xy,
        strain.yz + point_imposed.yz - element_imposed.yz, strain.xz + point_imposed.xz - element_imposed.xz};
}

SymmetricTensor3 expansion_adjusted_increment(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment, const adlite::Scalar& point_temperature,
    const adlite::Scalar& element_temperature, double committed_point_temperature, double committed_element_temperature,
    double time_step, MaterialFunctionContext context) {
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const SymmetricTensor3 current_point = material.eigenstrain(point_temperature, context);
    const SymmetricTensor3 current_element = material.eigenstrain(element_temperature, context);
    const SymmetricTensor3 old_point = material.eigenstrain(adlite::Scalar(committed_point_temperature), old_context);
    const SymmetricTensor3 old_element =
        material.eigenstrain(adlite::Scalar(committed_element_temperature), old_context);
    return {strain_increment.xx + old_element.xx - old_point.xx + current_point.xx - current_element.xx,
        strain_increment.yy + old_element.yy - old_point.yy + current_point.yy - current_element.yy,
        strain_increment.zz + old_element.zz - old_point.zz + current_point.zz - current_element.zz,
        strain_increment.xy + old_element.xy - old_point.xy + current_point.xy - current_element.xy,
        strain_increment.yz + old_element.yz - old_point.yz + current_point.yz - current_element.yz,
        strain_increment.xz + old_element.xz - old_point.xz + current_point.xz - current_element.xz};
}

std::array<double, 6> eigenstrain_temperature_derivative(
    const IsotropicThermoelasticMaterial& material, double temperature, MaterialFunctionContext context) {
    const adlite::Scalar active_temperature = adlite::Scalar::independent(temperature, 0, 1);
    const SymmetricTensor3 imposed = material.eigenstrain(active_temperature, context);
    const std::array<const adlite::Scalar*, 6> components = {
        &imposed.xx, &imposed.yy, &imposed.zz, &imposed.xy, &imposed.yz, &imposed.xz};
    std::array<double, 6> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        components[component]->copy_derivatives(&result[component], 1);
    return result;
}

} // namespace

SymmetricTensor3 rotate_cartesian_tensor(const SymmetricTensor3& tensor, const CartesianRotation& rotation) {
    const cartesian_detail::ActiveMatrix3 r = {{{rotation.xx, rotation.xy, rotation.xz},
        {rotation.yx, rotation.yy, rotation.yz}, {rotation.zx, rotation.zy, rotation.zz}}};
    const cartesian_detail::ActiveMatrix3 value = {
        {{tensor.xx, tensor.xy, tensor.xz}, {tensor.xy, tensor.yy, tensor.yz}, {tensor.xz, tensor.yz, tensor.zz}}};
    cartesian_detail::ActiveMatrix3 rotated{};
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
    const double value = cartesian_detail::determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX8 deformation must preserve a positive Jacobian");
}

namespace {
CartesianKinematics evaluate_cartesian_kinematics_from_gradient(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient, const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    CartesianKinematics result{};
    const cartesian_detail::Matrix3 old = deformation_gradient(point, committed_state);
    const cartesian_detail::KinematicsCore core =
        cartesian_detail::evaluate_kinematics(gradient, old, strain_formulation);
    result.strain_increment = core.strain_increment;
    result.rotation = core.rotation;
    result.current_weighted_measure = point.weighted_measure * core.current_determinant;
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.current_gradient[node][direction] +=
                    point.gradient[node][reference] * core.current_inverse[reference][direction];
    return result;
}
} // namespace

CartesianKinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state, const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    return evaluate_cartesian_kinematics_from_gradient(
        point, displacement_gradient(point, current_state), committed_state, strain_formulation);
}

namespace {
using cartesian_detail::material_context;

adlite::Scalar small_strain_element_pressure(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry, const Hex8LocalAdValues& state, double average_strain_trace, double time) {
    adlite::Scalar average_bulk_modulus = 0.0;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = state[hex8_node_to_gauss[q]];
        const ActiveThermoelasticProperties properties =
            material.active_properties(temperature, material_context(time, point.position));
        average_bulk_modulus += point.weighted_measure / geometry.reference_volume *
                                (properties.lame_lambda + 2.0 * properties.shear_modulus / 3.0);
    }
    const adlite::Scalar average_temperature = average_hex8_temperature(state);
    const SymmetricTensor3 imposed =
        material.eigenstrain(average_temperature, material_context(time, geometry.selective_position));
    return average_bulk_modulus * (average_strain_trace - imposed.xx - imposed.yy - imposed.zz);
}

struct SmallStrainElementPressureSystem final {
    double value = 0.0, trace_derivative = 0.0;
    std::array<double, hex8_node_count> temperature_derivatives{};
};

SmallStrainElementPressureSystem small_strain_element_pressure_system(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry, const Hex8LocalValues& state, double average_strain_trace, double time) {
    SmallStrainElementPressureSystem result;
    double average_bulk_modulus = 0.0;
    std::array<double, hex8_node_count> bulk_modulus_temperature_derivatives{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const std::size_t material_node = hex8_node_to_gauss[q];
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = adlite::Scalar::independent(state[material_node], 0, 1);
        const ActiveThermoelasticProperties properties =
            material.active_properties(temperature, material_context(time, point.position));
        const adlite::Scalar bulk_modulus = properties.lame_lambda + 2.0 * properties.shear_modulus / 3.0;
        const double normalized_weight = point.weighted_measure / geometry.reference_volume;
        average_bulk_modulus += normalized_weight * bulk_modulus.value();
        bulk_modulus_temperature_derivatives[material_node] =
            normalized_weight * (bulk_modulus.is_active() ? bulk_modulus.derivative(0) : 0.0);
    }

    const double average_temperature = average_hex8_temperature(state);
    const adlite::Scalar active_average_temperature = adlite::Scalar::independent(average_temperature, 0, 1);
    const SymmetricTensor3 imposed =
        material.eigenstrain(active_average_temperature, material_context(time, geometry.selective_position));
    const adlite::Scalar imposed_trace = imposed.xx + imposed.yy + imposed.zz;
    const double elastic_trace = average_strain_trace - imposed_trace.value();
    result.value = average_bulk_modulus * elastic_trace;
    result.trace_derivative = average_bulk_modulus;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        result.temperature_derivatives[node] = bulk_modulus_temperature_derivatives[node] * elastic_trace -
                                               average_bulk_modulus *
                                                   (imposed_trace.is_active() ? imposed_trace.derivative(0) : 0.0) /
                                                   static_cast<double>(hex8_node_count);
    return result;
}

struct FiniteElementPressureSystem final {
    double value = 0.0;
    std::array<double, hex8_local_dof_count> derivatives{};
};

FiniteElementPressureSystem finite_element_pressure_system(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry, const Hex8LocalValues& state, double time, const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* committed_material, double time_step,
    const FiniteAverageTraceSystem& average_trace_system) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const double expansion_temperature = average_hex8_temperature(state);
    FiniteElementPressureSystem result;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const std::size_t material_node = hex8_node_to_gauss[q];
        cartesian_detail::ActiveMatrix3 gradient{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                double value = 0.0;
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    value += point.gradient[node][direction] * state[8 * (component + 1) + node];
                gradient[component][direction] =
                    adlite::Scalar::independent(value, component * 3 + direction, point_width);
            }
        const adlite::Scalar active_temperature =
            adlite::Scalar::independent(state[material_node], temperature_index, point_width);
        const CartesianKinematics kinematics =
            evaluate_cartesian_kinematics_from_gradient(point, gradient, old_state, StrainFormulation::finite);
        const SymmetricTensor3 constitutive_strain = selectively_reduced_strain(
            kinematics.strain_increment, average_trace_system.value, StrainFormulation::finite);
        const std::array<const adlite::Scalar*, 6> strain_components = {&constitutive_strain.xx,
            &constitutive_strain.yy, &constitutive_strain.zz, &constitutive_strain.xy, &constitutive_strain.yz,
            &constitutive_strain.xz};
        const MaterialFunctionContext context = material_context(time, point.position);
        const SymmetricTensor3 point_imposed = material.eigenstrain(adlite::Scalar(state[material_node]), context);
        const SymmetricTensor3 element_imposed = material.eigenstrain(adlite::Scalar(expansion_temperature), context);
        const std::array<double, 6> point_imposed_values = {point_imposed.xx.value(), point_imposed.yy.value(),
            point_imposed.zz.value(), point_imposed.xy.value(), point_imposed.yz.value(), point_imposed.xz.value()};
        const std::array<double, 6> element_imposed_values = {element_imposed.xx.value(), element_imposed.yy.value(),
            element_imposed.zz.value(), element_imposed.xy.value(), element_imposed.yz.value(),
            element_imposed.xz.value()};
        std::array<double, 6> fed_strain{};
        const CartesianMaterialPointState* point_committed_material =
            committed_material == nullptr ? nullptr : &(*committed_material)[q];
        if (point_committed_material != nullptr) {
            const double old_temperature = old_state[material_node];
            if (!std::isfinite(old_temperature) || !(old_temperature > 0.0))
                throw std::domain_error(
                    "Incremental Cartesian material committed temperature must be finite and positive");
            MaterialFunctionContext old_context = context;
            old_context.time -= time_step;
            const SymmetricTensor3 old_element_imposed =
                material.eigenstrain(adlite::Scalar(average_hex8_temperature(old_state)), old_context);
            const std::array<double, 6> old_element_values = {old_element_imposed.xx.value(),
                old_element_imposed.yy.value(), old_element_imposed.zz.value(), old_element_imposed.xy.value(),
                old_element_imposed.yz.value(), old_element_imposed.xz.value()};
            for (std::size_t component = 0; component < 6; ++component)
                fed_strain[component] = point_committed_material->elastic_strain[component] +
                                        strain_components[component]->value() + old_element_values[component] +
                                        point_committed_material->plastic_strain[component] +
                                        point_committed_material->creep_strain[component] +
                                        point_imposed_values[component] - element_imposed_values[component];
        } else {
            for (std::size_t component = 0; component < 6; ++component)
                fed_strain[component] = strain_components[component]->value() + point_imposed_values[component] -
                                        element_imposed_values[component];
        }
        const cartesian_detail::CartesianStressTangent tangent = cartesian_detail::evaluate_stress_tangent(
            material, fed_strain, state[material_node], time_step, point_committed_material, context);
        const std::array<double, 6> point_imposed_derivative =
                                        eigenstrain_temperature_derivative(material, state[material_node], context),
                                    element_imposed_derivative =
                                        eigenstrain_temperature_derivative(material, expansion_temperature, context);
        std::array<double, 6> pressure_strain_derivatives{};
        for (std::size_t component = 0; component < 6; ++component)
            pressure_strain_derivatives[component] =
                (tangent.tangent[0][component] + tangent.tangent[1][component] + tangent.tangent[2][component]) / 3.0;
        double point_temperature_derivative = (tangent.thermal[0] + tangent.thermal[1] + tangent.thermal[2]) / 3.0,
               expansion_temperature_derivative = 0.0;
        for (std::size_t component = 0; component < 6; ++component) {
            point_temperature_derivative +=
                pressure_strain_derivatives[component] * point_imposed_derivative[component];
            expansion_temperature_derivative -=
                pressure_strain_derivatives[component] * element_imposed_derivative[component];
        }
        const double average_trace_derivative =
            (pressure_strain_derivatives[0] + pressure_strain_derivatives[1] + pressure_strain_derivatives[2]) / 3.0;
        std::array<adlite::Scalar, 7> compose_inputs{};
        for (std::size_t component = 0; component < 6; ++component)
            compose_inputs[component] = *strain_components[component];
        compose_inputs[6] = active_temperature;
        std::array<double, 7> partials{};
        for (std::size_t component = 0; component < 6; ++component)
            partials[component] = pressure_strain_derivatives[component];
        partials[6] = point_temperature_derivative;
        const double pressure_value = (tangent.stress.xx + tangent.stress.yy + tangent.stress.zz) / 3.0;
        const adlite::Scalar pressure =
            adlite::compose(pressure_value, compose_inputs.data(), partials.data(), compose_inputs.size());
        const double normalized_weight = point.weighted_measure / geometry.reference_volume;
        result.value += normalized_weight * pressure.value();
        std::array<double, point_width> derivatives{};
        pressure.copy_derivatives(derivatives.data(), derivatives.size());
        result.derivatives[material_node] += normalized_weight * derivatives[temperature_index];
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            result.derivatives[node] += normalized_weight * expansion_temperature_derivative / 8.0;
            for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component) {
                double local_chain = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    local_chain +=
                        derivatives[displacement_component * 3 + direction] * point.gradient[node][direction];
                result.derivatives[8 * (displacement_component + 1) + node] +=
                    normalized_weight *
                    (local_chain + average_trace_derivative *
                                       average_trace_system.displacement_derivatives[node][displacement_component]);
            }
        }
    }
    return result;
}

void add_hex8_point_residual(const Hex8QuadraturePoint& point, std::size_t material_node,
    const Hex8LocalAdValues& state, const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation, double time, const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material, double time_step, double average_strain_trace,
    const adlite::Scalar& element_pressure, double reference_volume, double finite_current_volume,
    double finite_element_pressure, Hex8LocalAdValues& residual) {
    const adlite::Scalar temperature = state[material_node];
    const adlite::Scalar expansion_temperature = average_hex8_temperature(state);
    const MaterialFunctionContext context = material_context(time, point.position);
    const adlite::Scalar conductivity = material.conductivity(temperature, context);
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const CartesianKinematics kinematics =
        evaluate_cartesian_incremental_kinematics(point, state, old_state, strain_formulation);
    const SymmetricTensor3 constitutive_strain =
        selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, strain_formulation);
    adlite::Scalar gradient_temperature_x = 0.0, gradient_temperature_y = 0.0, gradient_temperature_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        gradient_temperature_x += kinematics.current_gradient[node][0] * state[node];
        gradient_temperature_y += kinematics.current_gradient[node][1] * state[node];
        gradient_temperature_z += kinematics.current_gradient[node][2] * state[node];
    }
    SymmetricTensor3 stress;
    if (committed_material == nullptr) {
        stress = material.stress(
            expansion_adjusted_strain(material, constitutive_strain, temperature, expansion_temperature, context),
            temperature, context);
        if (strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
    } else if (strain_formulation == StrainFormulation::finite) {
        const double old_temperature = old_state[material_node];
        const double old_expansion_temperature = average_hex8_temperature(old_state);
        stress = material
                     .incremental_response(
                         expansion_adjusted_increment(material, constitutive_strain, temperature, expansion_temperature,
                             old_temperature, old_expansion_temperature, time_step, context),
                         kinematics.rotation, temperature, old_temperature, time_step, *committed_material, context)
                     .stress;
    } else {
        stress = material
                     .response(expansion_adjusted_strain(
                                   material, constitutive_strain, temperature, expansion_temperature, context),
                         temperature, time_step, *committed_material, context)
                     .stress;
    }
    if (strain_formulation == StrainFormulation::small) {
        const adlite::Scalar point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        stress.xx += element_pressure - point_pressure;
        stress.yy += element_pressure - point_pressure;
        stress.zz += element_pressure - point_pressure;
    } else if (strain_formulation == StrainFormulation::finite) {
        const adlite::Scalar point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        const adlite::Scalar deviatoric_scale =
            point.weighted_measure * finite_current_volume / (reference_volume * kinematics.current_weighted_measure);
        stress.xx = deviatoric_scale * (stress.xx - point_pressure) + finite_element_pressure;
        stress.yy = deviatoric_scale * (stress.yy - point_pressure) + finite_element_pressure;
        stress.zz = deviatoric_scale * (stress.zz - point_pressure) + finite_element_pressure;
        stress.xy *= deviatoric_scale;
        stress.yz *= deviatoric_scale;
        stress.xz *= deviatoric_scale;
    }
    for (std::size_t node = 0; node < 8; ++node) {
        const adlite::Scalar current_gradient_x = kinematics.current_gradient[node][0];
        const adlite::Scalar current_gradient_y = kinematics.current_gradient[node][1];
        const adlite::Scalar current_gradient_z = kinematics.current_gradient[node][2];
        residual[node] += kinematics.current_weighted_measure * conductivity *
                          (current_gradient_x * gradient_temperature_x + current_gradient_y * gradient_temperature_y +
                              current_gradient_z * gradient_temperature_z);
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

// Assembles one quadrature point's residual and exact 32-by-32 Jacobian with narrow AD. The
// kinematics chain is seeded on the nine displacement-gradient components plus the point
// associated nodal material temperature (width 10), while the constitutive evaluation uses width 7 and is reattached
// with adlite::compose, so stress rotation and current-configuration geometry derivatives
// remain exact. The final chain from (H, T_point) to the 32 local DOFs is linear and applied
// in closed form: dH[i][j]/du[b][d] = delta(i, d) * dN_b/dX_j. Abaqus first-order heat elements evaluate
// temperature-dependent material properties at the node associated with each standard Gauss point, so the material-
// temperature chain is a Kronecker delta while the temperature gradient retains the standard shape-function chain.
// Abaqus-style small-strain selective integration replaces the local trace by its volume
// average; that separate average-gradient chain is added directly. First-order thermal
// expansion similarly uses an arithmetic nodal average and a separate 1/8 chain. Neither
// rule widens the kinematics or constitutive AD evaluations.
void add_hex8_point_system(const Hex8QuadraturePoint& point, std::size_t material_node, const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time,
    const Hex8LocalValues* committed_state, const CartesianMaterialPointState* committed_material, double time_step,
    Hex8LocalAdValues& residual, Hex8LocalJacobian& jacobian, double average_strain_trace,
    const std::array<std::array<double, 3>, hex8_node_count>& average_trace_displacement_derivatives,
    const SmallStrainElementPressureSystem& element_pressure, double reference_volume,
    const FiniteAverageTraceSystem& finite_average_trace_system,
    const FiniteElementPressureSystem& finite_element_pressure) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    std::array<double, 9> gradient_values{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                gradient_values[component * 3 + direction] +=
                    point.gradient[node][direction] * state[8 * (component + 1) + node];
    const double temperature_value = state[material_node];
    const double expansion_temperature_value = average_hex8_temperature(state);
    cartesian_detail::ActiveMatrix3 active_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            active_gradient[component][direction] = adlite::Scalar::independent(
                gradient_values[component * 3 + direction], component * 3 + direction, point_width);
    const adlite::Scalar active_temperature =
        adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const CartesianKinematics kinematics =
        evaluate_cartesian_kinematics_from_gradient(point, active_gradient, old_state, strain_formulation);
    const SymmetricTensor3 constitutive_strain =
        selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, strain_formulation);
    const MaterialFunctionContext context = material_context(time, point.position);
    const std::array<const adlite::Scalar*, 6> strain_components = {&constitutive_strain.xx, &constitutive_strain.yy,
        &constitutive_strain.zz, &constitutive_strain.xy, &constitutive_strain.yz, &constitutive_strain.xz};
    std::array<double, 6> fed_strain{};
    const SymmetricTensor3 point_imposed = material.eigenstrain(adlite::Scalar(temperature_value), context);
    const SymmetricTensor3 element_imposed = material.eigenstrain(adlite::Scalar(expansion_temperature_value), context);
    const std::array<double, 6> point_imposed_values = {point_imposed.xx.value(), point_imposed.yy.value(),
        point_imposed.zz.value(), point_imposed.xy.value(), point_imposed.yz.value(), point_imposed.xz.value()};
    const std::array<double, 6> element_imposed_values = {element_imposed.xx.value(), element_imposed.yy.value(),
        element_imposed.zz.value(), element_imposed.xy.value(), element_imposed.yz.value(), element_imposed.xz.value()};
    if (committed_material != nullptr && strain_formulation == StrainFormulation::finite) {
        const double old_temperature = old_state[material_node];
        if (!std::isfinite(old_temperature) || !(old_temperature > 0.0))
            throw std::domain_error("Incremental Cartesian material committed temperature must be finite and positive");
        MaterialFunctionContext old_context = context;
        old_context.time -= time_step;
        const SymmetricTensor3 old_element_imposed =
            material.eigenstrain(adlite::Scalar(average_hex8_temperature(old_state)), old_context);
        const std::array<double, 6> old_element_values = {old_element_imposed.xx.value(),
            old_element_imposed.yy.value(), old_element_imposed.zz.value(), old_element_imposed.xy.value(),
            old_element_imposed.yz.value(), old_element_imposed.xz.value()};
        for (std::size_t component = 0; component < 6; ++component)
            fed_strain[component] = committed_material->elastic_strain[component] +
                                    strain_components[component]->value() + old_element_values[component] +
                                    committed_material->plastic_strain[component] +
                                    committed_material->creep_strain[component] + point_imposed_values[component] -
                                    element_imposed_values[component];
    } else {
        for (std::size_t component = 0; component < 6; ++component)
            fed_strain[component] = strain_components[component]->value() + point_imposed_values[component] -
                                    element_imposed_values[component];
    }
    const cartesian_detail::CartesianStressTangent tangent = cartesian_detail::evaluate_stress_tangent(
        material, fed_strain, temperature_value, time_step, committed_material, context);
    const std::array<double, 6> point_imposed_derivative =
                                    eigenstrain_temperature_derivative(material, temperature_value, context),
                                element_imposed_derivative =
                                    eigenstrain_temperature_derivative(material, expansion_temperature_value, context);
    std::array<double, 6> point_thermal{}, expansion_thermal{};
    std::array<double, 6> average_trace_stress_derivative{};
    for (std::size_t row = 0; row < 6; ++row) {
        point_thermal[row] = tangent.thermal[row];
        for (std::size_t component = 0; component < 6; ++component) {
            point_thermal[row] += tangent.tangent[row][component] * point_imposed_derivative[component];
            expansion_thermal[row] -= tangent.tangent[row][component] * element_imposed_derivative[component];
        }
        if (strain_formulation == StrainFormulation::small || strain_formulation == StrainFormulation::finite)
            average_trace_stress_derivative[row] =
                (tangent.tangent[row][0] + tangent.tangent[row][1] + tangent.tangent[row][2]) / 3.0;
    }
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
        partials[6] = point_thermal[component];
        composed[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    SymmetricTensor3 stress{composed[0], composed[1], composed[2], composed[3], composed[4], composed[5]};
    SymmetricTensor3 expansion_stress_derivative{expansion_thermal[0], expansion_thermal[1], expansion_thermal[2],
        expansion_thermal[3], expansion_thermal[4], expansion_thermal[5]};
    if (strain_formulation == StrainFormulation::small) {
        const adlite::Scalar point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        stress.xx += element_pressure.value - point_pressure;
        stress.yy += element_pressure.value - point_pressure;
        stress.zz += element_pressure.value - point_pressure;
        average_trace_stress_derivative = {element_pressure.trace_derivative, element_pressure.trace_derivative,
            element_pressure.trace_derivative, 0.0, 0.0, 0.0};
        expansion_stress_derivative = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    }
    SymmetricTensor3 finite_point_deviatoric{};
    if (strain_formulation == StrainFormulation::finite) {
        stress = rotate_cartesian_tensor(stress, kinematics.rotation);
        const CartesianRotation passive_rotation{kinematics.rotation.xx.value(), kinematics.rotation.xy.value(),
            kinematics.rotation.xz.value(), kinematics.rotation.yx.value(), kinematics.rotation.yy.value(),
            kinematics.rotation.yz.value(), kinematics.rotation.zx.value(), kinematics.rotation.zy.value(),
            kinematics.rotation.zz.value()};
        expansion_stress_derivative = rotate_cartesian_tensor(expansion_stress_derivative, passive_rotation);
        const SymmetricTensor3 rotated_average_trace_derivative =
            rotate_cartesian_tensor({average_trace_stress_derivative[0], average_trace_stress_derivative[1],
                                        average_trace_stress_derivative[2], average_trace_stress_derivative[3],
                                        average_trace_stress_derivative[4], average_trace_stress_derivative[5]},
                passive_rotation);
        average_trace_stress_derivative = {rotated_average_trace_derivative.xx.value(),
            rotated_average_trace_derivative.yy.value(), rotated_average_trace_derivative.zz.value(),
            rotated_average_trace_derivative.xy.value(), rotated_average_trace_derivative.yz.value(),
            rotated_average_trace_derivative.xz.value()};
        const adlite::Scalar point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        finite_point_deviatoric = {stress.xx - point_pressure, stress.yy - point_pressure, stress.zz - point_pressure,
            stress.xy, stress.yz, stress.xz};
        const adlite::Scalar deviatoric_scale = point.weighted_measure * finite_average_trace_system.current_volume /
                                                (reference_volume * kinematics.current_weighted_measure);
        stress = {deviatoric_scale * finite_point_deviatoric.xx + finite_element_pressure.value,
            deviatoric_scale * finite_point_deviatoric.yy + finite_element_pressure.value,
            deviatoric_scale * finite_point_deviatoric.zz + finite_element_pressure.value,
            deviatoric_scale * finite_point_deviatoric.xy, deviatoric_scale * finite_point_deviatoric.yz,
            deviatoric_scale * finite_point_deviatoric.xz};
        const adlite::Scalar expansion_pressure =
            (expansion_stress_derivative.xx + expansion_stress_derivative.yy + expansion_stress_derivative.zz) / 3.0;
        expansion_stress_derivative = {deviatoric_scale.value() * (expansion_stress_derivative.xx - expansion_pressure),
            deviatoric_scale.value() * (expansion_stress_derivative.yy - expansion_pressure),
            deviatoric_scale.value() * (expansion_stress_derivative.zz - expansion_pressure),
            deviatoric_scale.value() * expansion_stress_derivative.xy,
            deviatoric_scale.value() * expansion_stress_derivative.yz,
            deviatoric_scale.value() * expansion_stress_derivative.xz};
        const double average_trace_pressure = (average_trace_stress_derivative[0] + average_trace_stress_derivative[1] +
                                                  average_trace_stress_derivative[2]) /
                                              3.0;
        average_trace_stress_derivative = {
            deviatoric_scale.value() * (average_trace_stress_derivative[0] - average_trace_pressure),
            deviatoric_scale.value() * (average_trace_stress_derivative[1] - average_trace_pressure),
            deviatoric_scale.value() * (average_trace_stress_derivative[2] - average_trace_pressure),
            deviatoric_scale.value() * average_trace_stress_derivative[3],
            deviatoric_scale.value() * average_trace_stress_derivative[4],
            deviatoric_scale.value() * average_trace_stress_derivative[5]};
    }
    const adlite::Scalar conductivity = material.conductivity(active_temperature, context);
    adlite::Scalar gradient_temperature_x = 0.0, gradient_temperature_y = 0.0, gradient_temperature_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        gradient_temperature_x += kinematics.current_gradient[node][0] * state[node];
        gradient_temperature_y += kinematics.current_gradient[node][1] * state[node];
        gradient_temperature_z += kinematics.current_gradient[node][2] * state[node];
    }
    Hex8LocalAdValues point_residual{};
    point_residual.fill(adlite::Scalar(0.0));
    for (std::size_t node = 0; node < 8; ++node) {
        const adlite::Scalar current_gradient_x = kinematics.current_gradient[node][0];
        const adlite::Scalar current_gradient_y = kinematics.current_gradient[node][1];
        const adlite::Scalar current_gradient_z = kinematics.current_gradient[node][2];
        point_residual[node] +=
            kinematics.current_weighted_measure * conductivity *
            (current_gradient_x * gradient_temperature_x + current_gradient_y * gradient_temperature_y +
                current_gradient_z * gradient_temperature_z);
        point_residual[8 + node] +=
            kinematics.current_weighted_measure *
            (stress.xx * current_gradient_x + stress.xy * current_gradient_y + stress.xz * current_gradient_z);
        point_residual[16 + node] +=
            kinematics.current_weighted_measure *
            (stress.xy * current_gradient_x + stress.yy * current_gradient_y + stress.yz * current_gradient_z);
        point_residual[24 + node] +=
            kinematics.current_weighted_measure *
            (stress.xz * current_gradient_x + stress.yz * current_gradient_y + stress.zz * current_gradient_z);
    }
    std::array<double, point_width> derivatives{};
    const double conductivity_value = conductivity.value();
    for (std::size_t node = 0; node < 8; ++node) {
        const adlite::Scalar current_gradient_x = kinematics.current_gradient[node][0];
        const adlite::Scalar current_gradient_y = kinematics.current_gradient[node][1];
        const adlite::Scalar current_gradient_z = kinematics.current_gradient[node][2];
        point_residual[node].copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t other = 0; other < 8; ++other) {
            double gradient_dot = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                gradient_dot += kinematics.current_gradient[node][direction].value() *
                                kinematics.current_gradient[other][direction].value();
            jacobian[node * 32 + other] +=
                kinematics.current_weighted_measure.value() * conductivity_value * gradient_dot +
                (other == material_node ? derivatives[temperature_index] : 0.0);
        }
        for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component) {
            for (std::size_t other = 0; other < 8; ++other) {
                double chained = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    chained += derivatives[displacement_component * 3 + direction] * point.gradient[other][direction];
                jacobian[node * 32 + 8 * (displacement_component + 1) + other] += chained;
            }
        }
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t row = 8 * (component + 1) + node;
            point_residual[row].copy_derivatives(derivatives.data(), derivatives.size());
            const double expansion_derivative =
                component == 0   ? kinematics.current_weighted_measure.value() *
                                       (expansion_stress_derivative.xx.value() * current_gradient_x.value() +
                                           expansion_stress_derivative.xy.value() * current_gradient_y.value() +
                                           expansion_stress_derivative.xz.value() * current_gradient_z.value())
                : component == 1 ? kinematics.current_weighted_measure.value() *
                                       (expansion_stress_derivative.xy.value() * current_gradient_x.value() +
                                           expansion_stress_derivative.yy.value() * current_gradient_y.value() +
                                           expansion_stress_derivative.yz.value() * current_gradient_z.value())
                                 : kinematics.current_weighted_measure.value() *
                                       (expansion_stress_derivative.xz.value() * current_gradient_x.value() +
                                           expansion_stress_derivative.yz.value() * current_gradient_y.value() +
                                           expansion_stress_derivative.zz.value() * current_gradient_z.value());
            const double average_trace_derivative =
                component == 0   ? kinematics.current_weighted_measure.value() *
                                       (average_trace_stress_derivative[0] * current_gradient_x.value() +
                                           average_trace_stress_derivative[3] * current_gradient_y.value() +
                                           average_trace_stress_derivative[5] * current_gradient_z.value())
                : component == 1 ? kinematics.current_weighted_measure.value() *
                                       (average_trace_stress_derivative[3] * current_gradient_x.value() +
                                           average_trace_stress_derivative[1] * current_gradient_y.value() +
                                           average_trace_stress_derivative[4] * current_gradient_z.value())
                                 : kinematics.current_weighted_measure.value() *
                                       (average_trace_stress_derivative[5] * current_gradient_x.value() +
                                           average_trace_stress_derivative[4] * current_gradient_y.value() +
                                           average_trace_stress_derivative[2] * current_gradient_z.value());
            const double pressure_temperature_gradient = component == 0   ? current_gradient_x.value()
                                                         : component == 1 ? current_gradient_y.value()
                                                                          : current_gradient_z.value();
            const double finite_current_volume_derivative =
                strain_formulation == StrainFormulation::finite
                    ? point.weighted_measure / reference_volume *
                          (component == 0      ? finite_point_deviatoric.xx.value() * current_gradient_x.value() +
                                                     finite_point_deviatoric.xy.value() * current_gradient_y.value() +
                                                     finite_point_deviatoric.xz.value() * current_gradient_z.value()
                              : component == 1 ? finite_point_deviatoric.xy.value() * current_gradient_x.value() +
                                                     finite_point_deviatoric.yy.value() * current_gradient_y.value() +
                                                     finite_point_deviatoric.yz.value() * current_gradient_z.value()
                                               : finite_point_deviatoric.xz.value() * current_gradient_x.value() +
                                                     finite_point_deviatoric.yz.value() * current_gradient_y.value() +
                                                     finite_point_deviatoric.zz.value() * current_gradient_z.value())
                    : 0.0;
            for (std::size_t other = 0; other < 8; ++other) {
                if (other == material_node) jacobian[row * 32 + other] += derivatives[temperature_index];
                jacobian[row * 32 + other] += expansion_derivative / 8.0;
                jacobian[row * 32 + other] += kinematics.current_weighted_measure.value() *
                                              element_pressure.temperature_derivatives[other] *
                                              pressure_temperature_gradient;
                if (strain_formulation == StrainFormulation::finite)
                    jacobian[row * 32 + other] += kinematics.current_weighted_measure.value() *
                                                  pressure_temperature_gradient *
                                                  finite_element_pressure.derivatives[other];
                for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component) {
                    double chained = 0.0;
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        chained +=
                            derivatives[displacement_component * 3 + direction] * point.gradient[other][direction];
                    jacobian[row * 32 + 8 * (displacement_component + 1) + other] += chained;
                    jacobian[row * 32 + 8 * (displacement_component + 1) + other] +=
                        average_trace_derivative *
                        average_trace_displacement_derivatives[other][displacement_component];
                    if (strain_formulation == StrainFormulation::finite) {
                        const std::size_t column = 8 * (displacement_component + 1) + other;
                        jacobian[row * 32 + column] +=
                            finite_current_volume_derivative *
                            finite_average_trace_system.current_volume_derivatives[other][displacement_component];
                        jacobian[row * 32 + column] += kinematics.current_weighted_measure.value() *
                                                       pressure_temperature_gradient *
                                                       finite_element_pressure.derivatives[column];
                    }
                }
            }
        }
    }
    for (std::size_t row = 0; row < residual.size(); ++row) residual[row] += point_residual[row];
}

adlite::Scalar hex8_nodal_volume_measure(const Hex8Geometry& geometry, std::size_t node, const Hex8LocalAdValues& state,
    StrainFormulation strain_formulation) {
    const Hex8QuadraturePoint& point = geometry.points[hex8_node_to_gauss[node]];
    if (strain_formulation == StrainFormulation::small) return point.weighted_measure;
    return evaluate_cartesian_incremental_kinematics(point, state, Hex8LocalValues{}, strain_formulation)
        .current_weighted_measure;
}

struct Hex8NodalVolumeSystem final {
    double value;
    std::array<double, 9> gradient_derivatives;
};

Hex8NodalVolumeSystem hex8_nodal_volume_system(const Hex8Geometry& geometry, std::size_t node,
    const Hex8LocalValues& state, StrainFormulation strain_formulation) {
    const Hex8QuadraturePoint& point = geometry.points[hex8_node_to_gauss[node]];
    Hex8NodalVolumeSystem result{point.weighted_measure, {}};
    if (strain_formulation == StrainFormulation::small) return result;
    constexpr std::size_t width = 9;
    cartesian_detail::ActiveMatrix3 active_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction) {
            double value = 0.0;
            for (std::size_t other = 0; other < hex8_node_count; ++other)
                value += point.gradient[other][direction] * state[8 * (component + 1) + other];
            active_gradient[component][direction] =
                adlite::Scalar::independent(value, component * 3 + direction, width);
        }
    const adlite::Scalar measure =
        evaluate_cartesian_kinematics_from_gradient(point, active_gradient, Hex8LocalValues{}, strain_formulation)
            .current_weighted_measure;
    result.value = measure.value();
    measure.copy_derivatives(result.gradient_derivatives.data(), result.gradient_derivatives.size());
    return result;
}

void add_hex8_nodal_body_source(const Hex8Geometry& geometry, const Hex8LocalAdValues& state,
    StrainFormulation strain_formulation, double volumetric_heat_source, Hex8LocalAdValues& residual) {
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        residual[node] -= hex8_nodal_volume_measure(geometry, node, state, strain_formulation) * volumetric_heat_source;
}

void add_hex8_nodal_body_source_system(const Hex8Geometry& geometry, const Hex8LocalValues& state,
    StrainFormulation strain_formulation, double volumetric_heat_source, Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian) {
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const Hex8QuadraturePoint& point = geometry.points[hex8_node_to_gauss[node]];
        const Hex8NodalVolumeSystem measure = hex8_nodal_volume_system(geometry, node, state, strain_formulation);
        residual[node] -= measure.value * volumetric_heat_source;
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t other = 0; other < hex8_node_count; ++other) {
                double derivative = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    derivative +=
                        measure.gradient_derivatives[component * 3 + direction] * point.gradient[other][direction];
                jacobian[node * hex8_local_dof_count + 8 * (component + 1) + other] -=
                    volumetric_heat_source * derivative;
            }
    }
}

void add_hex8_lumped_capacity(const Hex8Geometry& geometry, const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state, const IsotropicThermoelasticMaterial& material, double time,
    double time_step, StrainFormulation strain_formulation, Hex8LocalAdValues& residual) {
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const Hex8CapacityPoint& point = geometry.capacity_points[node];
        const adlite::Scalar capacity = material.heat_capacity(state[node], material_context(time, point.position));
        residual[node] += hex8_nodal_volume_measure(geometry, node, state, strain_formulation) * capacity *
                          (state[node] - committed_state[node]) / time_step;
    }
}

void add_hex8_lumped_capacity_system(const Hex8Geometry& geometry, const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state, const IsotropicThermoelasticMaterial& material, double time,
    double time_step, StrainFormulation strain_formulation, Hex8LocalAdValues& residual, Hex8LocalJacobian& jacobian) {
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const Hex8CapacityPoint& point = geometry.capacity_points[node];
        const Hex8QuadraturePoint& geometry_point = geometry.points[hex8_node_to_gauss[node]];
        const Hex8NodalVolumeSystem measure = hex8_nodal_volume_system(geometry, node, state, strain_formulation);
        const adlite::Scalar temperature = adlite::Scalar::independent(state[node], 0, 1);
        const adlite::Scalar capacity = material.heat_capacity(temperature, material_context(time, point.position));
        const adlite::Scalar rate = capacity * (temperature - committed_state[node]) / time_step;
        const adlite::Scalar term = measure.value * rate;
        residual[node] += term.value();
        jacobian[node * hex8_local_dof_count + node] += term.derivative(0);
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t other = 0; other < hex8_node_count; ++other) {
                double derivative = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    derivative += measure.gradient_derivatives[component * 3 + direction] *
                                  geometry_point.gradient[other][direction];
                jacobian[node * hex8_local_dof_count + 8 * (component + 1) + other] += rate.value() * derivative;
            }
    }
}

struct ActiveReducedHex8Geometry final {
    adlite::Scalar volume{0.0}, center_measure{0.0};
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> average_gradient{};
    std::array<adlite::Scalar, hex8_node_count> shape_measures{};
    std::array<std::array<adlite::Scalar, 4>, hex8_node_count> hourglass_shape{};
    std::array<adlite::Scalar, 4> thermal_hourglass_coefficients{};
};

ActiveReducedHex8Geometry active_reduced_hex8_geometry(const Hex8Geometry& reference,
    const std::array<std::array<adlite::Scalar, 3>, hex8_node_count>& displacement, const char* configuration_name) {
    ActiveReducedHex8Geometry result;
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> current_coordinates{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        current_coordinates[node][0] = reference.capacity_points[node].position.x + displacement[node][0];
        current_coordinates[node][1] = reference.capacity_points[node].position.y + displacement[node][1];
        current_coordinates[node][2] = reference.capacity_points[node].position.z + displacement[node][2];
    }

    for (const Hex8QuadraturePoint& point : reference.points) {
        cartesian_detail::ActiveMatrix3 deformation{};
        for (std::size_t component = 0; component < 3; ++component) {
            deformation[component][component] = 1.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    deformation[component][direction] +=
                        displacement[node][component] * point.gradient[node][direction];
        }
        const adlite::Scalar determinant = cartesian_detail::determinant(deformation);
        if (!std::isfinite(determinant.value()) || !(determinant.value() > 0.0))
            throw std::domain_error(std::string("C3D8RT ") + configuration_name +
                                    " configuration must preserve positive Jacobians at all integration points");
        const cartesian_detail::ActiveMatrix3 inverse = cartesian_detail::inverse(deformation, determinant);
        const adlite::Scalar measure = point.weighted_measure * determinant;
        result.volume += measure;
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            result.shape_measures[node] += measure * point.shape[node];
            for (std::size_t current_direction = 0; current_direction < 3; ++current_direction) {
                adlite::Scalar current_gradient = 0.0;
                for (std::size_t reference_direction = 0; reference_direction < 3; ++reference_direction)
                    current_gradient +=
                        point.gradient[node][reference_direction] * inverse[reference_direction][current_direction];
                result.average_gradient[node][current_direction] += measure * current_gradient;
            }
        }
    }
    if (!std::isfinite(result.volume.value()) || !(result.volume.value() > 0.0))
        throw std::domain_error(std::string("C3D8RT ") + configuration_name + " volume must be finite and positive");
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result.average_gradient[node][component] /= result.volume;

    cartesian_detail::ActiveMatrix3 center_jacobian{};
    for (std::size_t physical = 0; physical < 3; ++physical)
        for (std::size_t natural = 0; natural < 3; ++natural)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                center_jacobian[physical][natural] +=
                    current_coordinates[node][physical] * hex8_signs[node][natural] / 8.0;
    result.center_measure = 8.0 * cartesian_detail::determinant(center_jacobian);
    if (!std::isfinite(result.center_measure.value()) || !(result.center_measure.value() > 0.0))
        throw std::domain_error(
            std::string("C3D8RT ") + configuration_name + " center Jacobian must be finite and positive");

    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<adlite::Scalar, 3> projected_coordinate{};
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                projected_coordinate[component] +=
                    current_coordinates[node][component] * finite_reduced_hex8_raw_hourglass[node][mode];
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            result.hourglass_shape[node][mode] = finite_reduced_hex8_raw_hourglass[node][mode];
            for (std::size_t component = 0; component < 3; ++component)
                result.hourglass_shape[node][mode] -=
                    result.average_gradient[node][component] * projected_coordinate[component];
        }
    }

    cartesian_detail::ActiveMatrix3 inverse_effective_mapping{};
    for (std::size_t natural = 0; natural < 3; ++natural)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                inverse_effective_mapping[natural][physical] +=
                    hex8_signs[node][natural] * result.average_gradient[node][physical];
    const adlite::Scalar inverse_effective_determinant = cartesian_detail::determinant(inverse_effective_mapping);
    if (!std::isfinite(inverse_effective_determinant.value()) || inverse_effective_determinant.value() == 0.0)
        throw std::domain_error(std::string("C3D8RT ") + configuration_name + " effective mapping must be nonsingular");
    const cartesian_detail::ActiveMatrix3 effective_mapping =
        cartesian_detail::inverse(inverse_effective_mapping, inverse_effective_determinant);
    cartesian_detail::ActiveMatrix3 metric{};
    for (std::size_t first = 0; first < 3; ++first)
        for (std::size_t second = 0; second < 3; ++second)
            for (std::size_t physical = 0; physical < 3; ++physical)
                metric[first][second] += effective_mapping[physical][first] * effective_mapping[physical][second];
    const adlite::Scalar first_pivot = metric[0][0];
    const adlite::Scalar second_pivot = metric[1][1] - metric[0][1] * metric[0][1] / first_pivot;
    const adlite::Scalar leading_determinant = metric[0][0] * metric[1][1] - metric[0][1] * metric[0][1];
    const adlite::Scalar third_pivot =
        metric[2][2] - (metric[1][1] * metric[0][2] * metric[0][2] - 2.0 * metric[0][1] * metric[0][2] * metric[1][2] +
                           metric[0][0] * metric[1][2] * metric[1][2]) /
                           leading_determinant;
    if (!std::isfinite(first_pivot.value()) || !std::isfinite(second_pivot.value()) ||
        !std::isfinite(third_pivot.value()) || !(first_pivot.value() > 0.0) || !(second_pivot.value() > 0.0) ||
        !(third_pivot.value() > 0.0))
        throw std::domain_error(
            std::string("C3D8RT ") + configuration_name + " effective metric must be positive definite");
    const adlite::Scalar thermal_scale = result.volume / 192.0;
    const adlite::Scalar inverse_x = 1.0 / first_pivot, inverse_y = 1.0 / second_pivot, inverse_z = 1.0 / third_pivot;
    result.thermal_hourglass_coefficients = {thermal_scale * (inverse_x + inverse_y),
        thermal_scale * (inverse_x + inverse_z), thermal_scale * (inverse_x + inverse_z),
        thermal_scale * (inverse_x + inverse_y + inverse_z) / 3.0};
    return result;
}

std::array<std::array<adlite::Scalar, 3>, hex8_node_count> reduced_hex8_displacement(const Hex8LocalAdValues& state) {
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> result{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[node][component] = state[8 * (component + 1) + node];
    return result;
}

std::array<std::array<adlite::Scalar, 3>, hex8_node_count> reduced_hex8_midpoint_displacement(
    const Hex8LocalAdValues& state, const Hex8LocalValues& committed_state) {
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> result{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[node][component] =
                0.5 * (state[8 * (component + 1) + node] + committed_state[8 * (component + 1) + node]);
    return result;
}

cartesian_detail::ActiveMatrix3 reduced_hex8_central_gradient(
    const ActiveReducedHex8Geometry& midpoint, const Hex8LocalAdValues& state, const Hex8LocalValues& committed_state) {
    cartesian_detail::ActiveMatrix3 central_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                central_gradient[component][direction] +=
                    (state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node]) *
                    midpoint.average_gradient[node][direction];
    return central_gradient;
}

cartesian_detail::KinematicsCore reduced_hex8_finite_kinematics(
    const ActiveReducedHex8Geometry& midpoint, const Hex8LocalAdValues& state, const Hex8LocalValues& committed_state) {
    return cartesian_detail::evaluate_hughes_winget_increment(
        reduced_hex8_central_gradient(midpoint, state, committed_state));
}

struct ReducedFiniteMaterialLinearization final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 6>, 6> tangent{};
    std::array<double, 6> thermal{};
};

struct ReducedFiniteStressLinearization final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 10>, 6> tangent{};
};

ReducedFiniteMaterialLinearization reduced_finite_material_linearization(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3Values& strain_increment, double temperature, double committed_temperature, double time_step,
    const CartesianMaterialPointState* committed_material, MaterialFunctionContext context) {
    std::array<double, 7> seeds = {strain_increment.xx, strain_increment.yy, strain_increment.zz, strain_increment.xy,
        strain_increment.yz, strain_increment.xz, temperature};
    std::array<adlite::Scalar, 7> active{};
    adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    const SymmetricTensor3 strain{active[0], active[1], active[2], active[3], active[4], active[5]};
    SymmetricTensor3 stress;
    if (committed_material == nullptr)
        stress = material.stress(strain, active[6], context);
    else
        stress = material
                     .incremental_response(strain, CartesianRotation{}, active[6], committed_temperature, time_step,
                         *committed_material, context)
                     .stress;
    ReducedFiniteMaterialLinearization result;
    result.stress = {stress.xx.value(), stress.yy.value(), stress.zz.value(), stress.xy.value(), stress.yz.value(),
        stress.xz.value()};
    const std::array<const adlite::Scalar*, 6> components = {
        &stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz};
    std::array<double, 7> derivatives{};
    for (std::size_t row = 0; row < 6; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 6; ++column) result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[6];
    }
    return result;
}

ReducedFiniteStressLinearization reduced_finite_stress_linearization(
    const cartesian_detail::ActiveMatrix3& passive_gradient, double temperature,
    const ReducedFiniteMaterialLinearization& material_linearization) {
    std::array<double, 10> values{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            values[3 * component + direction] = passive_gradient[component][direction].value();
    values[9] = temperature;
    std::array<adlite::Scalar, 10> active{};
    adlite::seed_identity(values.data(), values.size(), active.data());
    cartesian_detail::ActiveMatrix3 gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            gradient[component][direction] = active[3 * component + direction];
    const cartesian_detail::KinematicsCore kinematics = cartesian_detail::evaluate_hughes_winget_increment(gradient);
    const std::array<adlite::Scalar, 7> material_inputs = {kinematics.strain_increment.xx,
        kinematics.strain_increment.yy, kinematics.strain_increment.zz, kinematics.strain_increment.xy,
        kinematics.strain_increment.yz, kinematics.strain_increment.xz, active[9]};
    const std::array<double, 6> stress_values = {material_linearization.stress.xx, material_linearization.stress.yy,
        material_linearization.stress.zz, material_linearization.stress.xy, material_linearization.stress.yz,
        material_linearization.stress.xz};
    std::array<adlite::Scalar, 6> material_stress{};
    std::array<double, 7> partials{};
    for (std::size_t row = 0; row < 6; ++row) {
        for (std::size_t column = 0; column < 6; ++column)
            partials[column] = material_linearization.tangent[row][column];
        partials[6] = material_linearization.thermal[row];
        material_stress[row] =
            adlite::compose(stress_values[row], material_inputs.data(), partials.data(), material_inputs.size());
    }
    const SymmetricTensor3 rotated =
        rotate_cartesian_tensor({material_stress[0], material_stress[1], material_stress[2], material_stress[3],
                                    material_stress[4], material_stress[5]},
            kinematics.rotation);
    const std::array<const adlite::Scalar*, 6> components = {
        &rotated.xx, &rotated.yy, &rotated.zz, &rotated.xy, &rotated.yz, &rotated.xz};
    ReducedFiniteStressLinearization result;
    result.stress = {rotated.xx.value(), rotated.yy.value(), rotated.zz.value(), rotated.xy.value(), rotated.yz.value(),
        rotated.xz.value()};
    for (std::size_t row = 0; row < 6; ++row)
        components[row]->copy_derivatives(result.tangent[row].data(), result.tangent[row].size());
    return result;
}

adlite::Scalar reduced_hex8_temperature(const ActiveReducedHex8Geometry& geometry, const Hex8LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        result += geometry.shape_measures[node] * state[node] / geometry.volume;
    return result;
}

Hex8LocalAdValues reduced_hex8_finite_residual(const CartesianThermoelasticData& data, const Hex8Geometry& reference,
    const Hex8LocalAdValues& state, const Hex8LocalValues& committed_state, double time_step,
    bool include_thermal_time_term, double initial_shear_modulus,
    const ReducedFiniteStressLinearization& stress_linearization) {
    const auto displacement = reduced_hex8_displacement(state);
    const auto midpoint_displacement = reduced_hex8_midpoint_displacement(state, committed_state);
    const ActiveReducedHex8Geometry current = active_reduced_hex8_geometry(reference, displacement, "current");
    const ActiveReducedHex8Geometry midpoint =
        active_reduced_hex8_geometry(reference, midpoint_displacement, "midpoint");
    const cartesian_detail::ActiveMatrix3 central_gradient =
        reduced_hex8_central_gradient(midpoint, state, committed_state);
    const adlite::Scalar temperature = reduced_hex8_temperature(current, state);
    const MaterialFunctionContext context = material_context(data.time, reference.reduced_point.position);

    std::array<adlite::Scalar, 10> compose_inputs{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            compose_inputs[3 * component + direction] = central_gradient[component][direction];
    compose_inputs[9] = temperature;
    const std::array<double, 6> stress_values = {stress_linearization.stress.xx, stress_linearization.stress.yy,
        stress_linearization.stress.zz, stress_linearization.stress.xy, stress_linearization.stress.yz,
        stress_linearization.stress.xz};
    std::array<adlite::Scalar, 6> composed_stress{};
    for (std::size_t row = 0; row < 6; ++row) {
        composed_stress[row] = adlite::compose(
            stress_values[row], compose_inputs.data(), stress_linearization.tangent[row].data(), compose_inputs.size());
    }
    const SymmetricTensor3 stress{composed_stress[0], composed_stress[1], composed_stress[2], composed_stress[3],
        composed_stress[4], composed_stress[5]};

    const adlite::Scalar conductivity = data.material.conductivity(temperature, context);
    std::array<adlite::Scalar, 3> temperature_gradient{};
    std::array<adlite::Scalar, 4> temperature_hourglass_amplitude{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += current.average_gradient[node][direction] * state[node];
        for (std::size_t mode = 0; mode < 4; ++mode)
            temperature_hourglass_amplitude[mode] += current.hourglass_shape[node][mode] * state[node];
    }

    Hex8LocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        adlite::Scalar uniform_thermal = 0.0, hourglass_thermal = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            uniform_thermal += current.average_gradient[node][direction] * temperature_gradient[direction];
        for (std::size_t mode = 0; mode < 4; ++mode)
            hourglass_thermal += current.hourglass_shape[node][mode] * current.thermal_hourglass_coefficients[mode] *
                                 temperature_hourglass_amplitude[mode];
        residual[node] += conductivity * (current.volume * uniform_thermal + hourglass_thermal) -
                          current.center_measure * data.volumetric_heat_source / 8.0;
        const adlite::Scalar gradient_x = current.average_gradient[node][0];
        const adlite::Scalar gradient_y = current.average_gradient[node][1];
        const adlite::Scalar gradient_z = current.average_gradient[node][2];
        residual[8 + node] +=
            current.volume * (stress.xx * gradient_x + stress.xy * gradient_y + stress.xz * gradient_z);
        residual[16 + node] +=
            current.volume * (stress.xy * gradient_x + stress.yy * gradient_y + stress.yz * gradient_z);
        residual[24 + node] +=
            current.volume * (stress.xz * gradient_x + stress.yz * gradient_y + stress.zz * gradient_z);
    }

    // Abaqus/Standard's finite-strain total-stiffness control transports each
    // reference modal displacement through the element-average deformation
    // before applying the reference hourglass metric.  Differentiate the
    // resulting quadratic modal energy exactly: the first force term is the
    // F D F^T push-forward, while the second is the derivative of F itself.
    // The fixed 0.005 factor is Abaqus's documented default and is not
    // calibrated from the verification probes.
    cartesian_detail::ActiveMatrix3 average_deformation{};
    for (std::size_t component = 0; component < 3; ++component) {
        average_deformation[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                average_deformation[component][direction] +=
                    state[8 * (component + 1) + node] * reference.average_shape_gradient[node][direction];
    }
    constexpr double abaqus_total_stiffness_factor = 0.005;
    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<adlite::Scalar, 3> reference_amplitude{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                reference_amplitude[component] +=
                    reference.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
        std::array<adlite::Scalar, 3> material_modal_force{};
        for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
            adlite::Scalar transported_amplitude = 0.0;
            for (std::size_t component = 0; component < 3; ++component)
                transported_amplitude +=
                    reference_amplitude[component] * average_deformation[component][material_direction];
            material_modal_force[material_direction] = abaqus_total_stiffness_factor * initial_shear_modulus *
                                                       reference.mechanical_hourglass_metrics[material_direction] *
                                                       transported_amplitude;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            adlite::Scalar spatial_modal_force = 0.0;
            for (std::size_t material_direction = 0; material_direction < 3; ++material_direction)
                spatial_modal_force +=
                    average_deformation[component][material_direction] * material_modal_force[material_direction];
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                adlite::Scalar deformation_gradient_term = 0.0;
                for (std::size_t material_direction = 0; material_direction < 3; ++material_direction)
                    deformation_gradient_term += reference.average_shape_gradient[node][material_direction] *
                                                 material_modal_force[material_direction];
                residual[8 * (component + 1) + node] += reference.hourglass_shape[node][mode] * spatial_modal_force +
                                                        reference_amplitude[component] * deformation_gradient_term;
            }
        }
    }

    if (include_thermal_time_term)
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            const adlite::Scalar capacity = data.material.heat_capacity(
                state[node], material_context(data.time, reference.capacity_points[node].position));
            residual[node] +=
                current.shape_measures[node] * capacity * (state[node] - committed_state[node]) / time_step;
        }
    return residual;
}

void add_reduced_hex8_finite_strain_system(const CartesianThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material, double time_step, bool include_thermal_time_term,
    Hex8LocalAdValues& residual, Hex8LocalJacobian* jacobian) {
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    Hex8LocalAdValues passive_state{};
    ad_local_system::make_passive(state.data(), state.size(), passive_state.data());
    const ActiveReducedHex8Geometry passive_current =
        active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive_state), "current");
    const ActiveReducedHex8Geometry passive_midpoint = active_reduced_hex8_geometry(
        geometry, reduced_hex8_midpoint_displacement(passive_state, old_state), "midpoint");
    const cartesian_detail::KinematicsCore passive_kinematics =
        reduced_hex8_finite_kinematics(passive_midpoint, passive_state, old_state);
    const adlite::Scalar passive_temperature = reduced_hex8_temperature(passive_current, passive_state);
    double old_temperature = data.initial_temperature;
    if (committed_material != nullptr) {
        Hex8LocalAdValues passive_old_state{};
        ad_local_system::make_passive(old_state.data(), old_state.size(), passive_old_state.data());
        const ActiveReducedHex8Geometry old_geometry =
            active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive_old_state), "committed");
        old_temperature = reduced_hex8_temperature(old_geometry, passive_old_state).value();
    }
    const SymmetricTensor3Values strain_values{passive_kinematics.strain_increment.xx.value(),
        passive_kinematics.strain_increment.yy.value(), passive_kinematics.strain_increment.zz.value(),
        passive_kinematics.strain_increment.xy.value(), passive_kinematics.strain_increment.yz.value(),
        passive_kinematics.strain_increment.xz.value()};
    const MaterialFunctionContext context = material_context(data.time, geometry.reduced_point.position);
    const ReducedFiniteMaterialLinearization material_linearization =
        reduced_finite_material_linearization(data.material, strain_values, passive_temperature.value(),
            old_temperature, time_step, committed_material, context);
    const ReducedFiniteStressLinearization stress_linearization =
        reduced_finite_stress_linearization(reduced_hex8_central_gradient(passive_midpoint, passive_state, old_state),
            passive_temperature.value(), material_linearization);
    const ActiveThermoelasticProperties initial_properties = data.material.active_properties(
        adlite::Scalar(data.initial_temperature), material_context(0.0, geometry.reduced_point.position));
    const double initial_shear_modulus = initial_properties.shear_modulus.value();
    if (!std::isfinite(initial_shear_modulus) || !(initial_shear_modulus > 0.0))
        throw std::invalid_argument("C3D8RT initial shear modulus must be finite and positive");

    const bool add_thermal_time_term = committed_state != nullptr && include_thermal_time_term;
    residual = reduced_hex8_finite_residual(data, geometry, passive_state, old_state, time_step, add_thermal_time_term,
        initial_shear_modulus, stress_linearization);
    if (jacobian == nullptr) return;
    jacobian->fill(0.0);
    // Current and midpoint geometry depend on all 24 displacement values, while material
    // properties and capacity depend on all eight temperatures.  Evaluate those two independent
    // blocks once each after composing the width-7 material tangent into the width-10 kinematics
    // chain above.  This preserves the prohibition on a complete 32-DOF identity seed and avoids
    // rebuilding the full reduced-integration residual once per Jacobian column.
    Hex8LocalAdValues displacement_state{};
    ad_local_system::make_passive(state.data(), state.size(), displacement_state.data());
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            const std::size_t column = 8 * (component + 1) + node;
            displacement_state[column] = adlite::Scalar::independent(state[column], 8 * component + node, 24);
        }
    const Hex8LocalAdValues displacement_residual = reduced_hex8_finite_residual(data, geometry, displacement_state,
        old_state, time_step, add_thermal_time_term, initial_shear_modulus, stress_linearization);
    std::array<double, 24> displacement_derivatives{};
    for (std::size_t row = 0; row < hex8_local_dof_count; ++row) {
        if (!displacement_residual[row].is_active()) continue;
        displacement_residual[row].copy_derivatives(displacement_derivatives.data(), displacement_derivatives.size());
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                (*jacobian)[row * hex8_local_dof_count + 8 * (component + 1) + node] =
                    displacement_derivatives[8 * component + node];
    }

    Hex8LocalAdValues temperature_state{};
    ad_local_system::make_passive(state.data(), state.size(), temperature_state.data());
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        temperature_state[node] = adlite::Scalar::independent(state[node], node, hex8_node_count);
    const Hex8LocalAdValues temperature_residual = reduced_hex8_finite_residual(data, geometry, temperature_state,
        old_state, time_step, add_thermal_time_term, initial_shear_modulus, stress_linearization);
    std::array<double, hex8_node_count> temperature_derivatives{};
    for (std::size_t row = 0; row < hex8_local_dof_count; ++row) {
        if (!temperature_residual[row].is_active()) continue;
        temperature_residual[row].copy_derivatives(temperature_derivatives.data(), temperature_derivatives.size());
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            (*jacobian)[row * hex8_local_dof_count + node] = temperature_derivatives[node];
    }
}

void add_reduced_hex8_small_strain_system(const CartesianThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material, double time_step, bool include_thermal_time_term,
    Hex8LocalAdValues& residual, Hex8LocalJacobian* jacobian) {
    if (data.strain_formulation != StrainFormulation::small)
        throw std::logic_error("C3D8RT small-strain integration received a non-small strain formulation");
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");

    constexpr std::size_t point_width = 10, temperature_index = 9;
    const Hex8QuadraturePoint& point = geometry.reduced_point;
    std::array<double, 9> gradient_values{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                gradient_values[component * 3 + direction] +=
                    point.gradient[node][direction] * state[8 * (component + 1) + node];
    const double temperature_value = reduced_hex8_temperature(geometry, state);
    cartesian_detail::ActiveMatrix3 active_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            active_gradient[component][direction] = adlite::Scalar::independent(
                gradient_values[component * 3 + direction], component * 3 + direction, point_width);
    const adlite::Scalar active_temperature =
        adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    const SymmetricTensor3 strain{active_gradient[0][0], active_gradient[1][1], active_gradient[2][2],
        0.5 * (active_gradient[0][1] + active_gradient[1][0]), 0.5 * (active_gradient[1][2] + active_gradient[2][1]),
        0.5 * (active_gradient[0][2] + active_gradient[2][0])};
    const std::array<const adlite::Scalar*, 6> strain_components = {
        &strain.xx, &strain.yy, &strain.zz, &strain.xy, &strain.yz, &strain.xz};
    std::array<double, 6> fed_strain{};
    for (std::size_t component = 0; component < fed_strain.size(); ++component)
        fed_strain[component] = strain_components[component]->value();
    const MaterialFunctionContext context = material_context(data.time, point.position);
    const ActiveThermoelasticProperties initial_properties = data.material.active_properties(
        adlite::Scalar(data.initial_temperature), material_context(0.0, point.position));
    const double reference_shear_modulus = initial_properties.shear_modulus.value();
    if (!std::isfinite(reference_shear_modulus) || !(reference_shear_modulus > 0.0))
        throw std::invalid_argument("C3D8RT initial shear modulus must be finite and positive");
    const cartesian_detail::CartesianStressTangent tangent = cartesian_detail::evaluate_stress_tangent(
        data.material, fed_strain, temperature_value, time_step, committed_material, context);
    std::array<adlite::Scalar, 7> compose_inputs{};
    for (std::size_t component = 0; component < 6; ++component)
        compose_inputs[component] = *strain_components[component];
    compose_inputs[6] = active_temperature;
    const std::array<double, 6> stress_values = {tangent.stress.xx, tangent.stress.yy, tangent.stress.zz,
        tangent.stress.xy, tangent.stress.yz, tangent.stress.xz};
    std::array<double, 7> partials{};
    std::array<adlite::Scalar, 6> composed_stress{};
    for (std::size_t component = 0; component < 6; ++component) {
        for (std::size_t column = 0; column < 6; ++column) partials[column] = tangent.tangent[component][column];
        partials[6] = tangent.thermal[component];
        composed_stress[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    const SymmetricTensor3 stress{composed_stress[0], composed_stress[1], composed_stress[2], composed_stress[3],
        composed_stress[4], composed_stress[5]};

    const adlite::Scalar active_conductivity_temperature = adlite::Scalar::independent(temperature_value, 0, 1);
    const adlite::Scalar conductivity = data.material.conductivity(active_conductivity_temperature, context);
    const double conductivity_derivative = conductivity.is_active() ? conductivity.derivative(0) : 0.0;
    std::array<double, 3> temperature_gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += point.gradient[node][direction] * state[node];
    std::array<double, 4> temperature_hourglass_amplitudes{};
    for (std::size_t mode = 0; mode < 4; ++mode)
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            temperature_hourglass_amplitudes[mode] += geometry.hourglass_shape[node][mode] * state[node];

    std::array<double, point_width> derivatives{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        double uniform_thermal = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            uniform_thermal += point.gradient[node][direction] * temperature_gradient[direction];
        double hourglass_thermal = 0.0;
        for (std::size_t mode = 0; mode < 4; ++mode)
            hourglass_thermal += geometry.hourglass_shape[node][mode] * geometry.thermal_hourglass_coefficients[mode] *
                                 temperature_hourglass_amplitudes[mode];
        const double thermal_operator = geometry.reference_volume * uniform_thermal + hourglass_thermal;
        residual[node] += conductivity.value() * thermal_operator -
                          geometry.reduced_body_source_measure * data.volumetric_heat_source / 8.0;

        const std::array<adlite::Scalar, 3> mechanical = {
            geometry.reference_volume * (stress.xx * point.gradient[node][0] + stress.xy * point.gradient[node][1] +
                                            stress.xz * point.gradient[node][2]),
            geometry.reference_volume * (stress.xy * point.gradient[node][0] + stress.yy * point.gradient[node][1] +
                                            stress.yz * point.gradient[node][2]),
            geometry.reference_volume * (stress.xz * point.gradient[node][0] + stress.yz * point.gradient[node][1] +
                                            stress.zz * point.gradient[node][2])};
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t row = 8 * (component + 1) + node;
            residual[row] += mechanical[component].value();
            if (jacobian == nullptr) continue;
            mechanical[component].copy_derivatives(derivatives.data(), derivatives.size());
            for (std::size_t other = 0; other < hex8_node_count; ++other) {
                (*jacobian)[row * hex8_local_dof_count + other] +=
                    derivatives[temperature_index] * geometry.reduced_capacity_points[other].weighted_measure /
                    geometry.reference_volume;
                for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component) {
                    double derivative = 0.0;
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        derivative +=
                            derivatives[displacement_component * 3 + direction] * point.gradient[other][direction];
                    (*jacobian)[row * hex8_local_dof_count + 8 * (displacement_component + 1) + other] += derivative;
                }
            }
        }

        if (jacobian != nullptr)
            for (std::size_t other = 0; other < hex8_node_count; ++other) {
                double uniform_tangent = 0.0, hourglass_tangent = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    uniform_tangent += point.gradient[node][direction] * point.gradient[other][direction];
                for (std::size_t mode = 0; mode < 4; ++mode)
                    hourglass_tangent += geometry.hourglass_shape[node][mode] *
                                         geometry.thermal_hourglass_coefficients[mode] *
                                         geometry.hourglass_shape[other][mode];
                (*jacobian)[node * hex8_local_dof_count + other] +=
                    conductivity.value() * (geometry.reference_volume * uniform_tangent + hourglass_tangent) +
                    conductivity_derivative * thermal_operator *
                        geometry.reduced_capacity_points[other].weighted_measure / geometry.reference_volume;
            }
    }

    constexpr double abaqus_total_stiffness_factor = 0.005;
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t mode = 0; mode < 4; ++mode) {
            double amplitude = 0.0;
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                amplitude += geometry.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
            const double stiffness = abaqus_total_stiffness_factor * reference_shear_modulus *
                                     geometry.mechanical_hourglass_metrics[component];
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                const std::size_t row = 8 * (component + 1) + node;
                residual[row] += stiffness * geometry.hourglass_shape[node][mode] * amplitude;
                if (jacobian != nullptr)
                    for (std::size_t other = 0; other < hex8_node_count; ++other)
                        (*jacobian)[row * hex8_local_dof_count + 8 * (component + 1) + other] +=
                            stiffness * geometry.hourglass_shape[node][mode] * geometry.hourglass_shape[other][mode];
            }
        }

    if (committed_state != nullptr && include_thermal_time_term)
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            const Hex8CapacityPoint& capacity_point = geometry.reduced_capacity_points[node];
            const MaterialFunctionContext capacity_context = material_context(data.time, capacity_point.position);
            const adlite::Scalar active_capacity_temperature = adlite::Scalar::independent(state[node], 0, 1);
            const adlite::Scalar capacity = data.material.heat_capacity(active_capacity_temperature, capacity_context);
            const adlite::Scalar capacity_term = capacity_point.weighted_measure * capacity *
                                                 (active_capacity_temperature - (*committed_state)[node]) / time_step;
            residual[node] += capacity_term.value();
            if (jacobian != nullptr) (*jacobian)[node * hex8_local_dof_count + node] += capacity_term.derivative(0);
        }
}

std::array<SymmetricTensor3Values, 8> evaluate_hex8_stress(const Hex8Geometry& geometry, const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time) {
    Hex8LocalAdValues ad_state{};
    ad_local_system::make_passive(state.data(), state.size(), ad_state.data());
    const adlite::Scalar expansion_temperature = average_hex8_temperature(ad_state);
    const Hex8LocalValues undeformed{};
    const double average_strain_trace = strain_formulation == StrainFormulation::finite
                                            ? finite_average_hex8_strain_trace(geometry, state, undeformed)
                                            : average_hex8_strain_trace(geometry, state);
    std::array<SymmetricTensor3Values, 8> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = ad_state[hex8_node_to_gauss[q]];
        const CartesianKinematics kinematics =
            evaluate_cartesian_incremental_kinematics(point, ad_state, Hex8LocalValues{}, strain_formulation);
        const SymmetricTensor3 constitutive_strain =
            selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, strain_formulation);
        const MaterialFunctionContext context = material_context(time, point.position);
        SymmetricTensor3 stress = material.stress(
            expansion_adjusted_strain(material, constitutive_strain, temperature, expansion_temperature, context),
            temperature, context);
        if (strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
        result[q] = {stress.xx.value(), stress.yy.value(), stress.zz.value(), stress.xy.value(), stress.yz.value(),
            stress.xz.value()};
    }
    return result;
}

Hex8LocalResidual compute_hex8_local(const CartesianThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state, const CartesianMaterialHistory* history,
    double time_step, Hex8LocalJacobian* jacobian, bool include_thermal_time_term) {
    if (committed_state != nullptr && (!std::isfinite(time_step) || time_step <= 0.0))
        throw std::invalid_argument("HEX8 time step must be finite and positive");
    const std::size_t expected_material_points =
        data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt ? 1U : geometry.points.size();
    if (history != nullptr && history->size() != expected_material_points)
        throw std::invalid_argument(data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt
                                        ? "C3D8RT material history must contain one integration point"
                                        : "C3D8T material history must contain eight integration points");
    Hex8LocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    if (data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt) {
        if (jacobian != nullptr) jacobian->fill(0.0);
        if (data.strain_formulation == StrainFormulation::finite)
            add_reduced_hex8_finite_strain_system(data, geometry, state, committed_state,
                history == nullptr ? nullptr : &history->front(), time_step, include_thermal_time_term, residual,
                jacobian);
        else
            add_reduced_hex8_small_strain_system(data, geometry, state, committed_state,
                history == nullptr ? nullptr : &history->front(), time_step, include_thermal_time_term, residual,
                jacobian);
        Hex8LocalResidual result{};
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
        return result;
    }
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    double average_strain_trace = average_hex8_strain_trace(geometry, state);
    std::array<std::array<double, 3>, hex8_node_count> average_trace_displacement_derivatives =
        geometry.average_shape_gradient;
    FiniteAverageTraceSystem finite_average_trace;
    if (data.strain_formulation == StrainFormulation::finite) {
        finite_average_trace = finite_average_hex8_strain_trace_system(geometry, state, old_state);
        average_strain_trace = finite_average_trace.value;
        average_trace_displacement_derivatives = finite_average_trace.displacement_derivatives;
    }
    const FiniteElementPressureSystem finite_element_pressure =
        data.strain_formulation == StrainFormulation::finite
            ? finite_element_pressure_system(
                  data.material, geometry, state, data.time, committed_state, history, time_step, finite_average_trace)
            : FiniteElementPressureSystem{};
    if (jacobian == nullptr) {
        Hex8LocalAdValues active{};
        ad_local_system::make_passive(state.data(), state.size(), active.data());
        const adlite::Scalar element_pressure =
            data.strain_formulation == StrainFormulation::small
                ? small_strain_element_pressure(data.material, geometry, active, average_strain_trace, data.time)
                : adlite::Scalar(0.0);
        for (std::size_t q = 0; q < geometry.points.size(); ++q)
            add_hex8_point_residual(geometry.points[q], hex8_node_to_gauss[q], active, data.material,
                data.strain_formulation, data.time, committed_state, history == nullptr ? nullptr : &(*history)[q],
                time_step, average_strain_trace, element_pressure, geometry.reference_volume,
                finite_average_trace.current_volume, finite_element_pressure.value, residual);
        add_hex8_nodal_body_source(geometry, active, data.strain_formulation, data.volumetric_heat_source, residual);
        if (committed_state != nullptr && include_thermal_time_term)
            add_hex8_lumped_capacity(geometry, active, *committed_state, data.material, data.time, time_step,
                data.strain_formulation, residual);
        Hex8LocalResidual result{};
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
        return result;
    }
    jacobian->fill(0.0);
    const SmallStrainElementPressureSystem element_pressure =
        data.strain_formulation == StrainFormulation::small
            ? small_strain_element_pressure_system(data.material, geometry, state, average_strain_trace, data.time)
            : SmallStrainElementPressureSystem{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q)
        add_hex8_point_system(geometry.points[q], hex8_node_to_gauss[q], state, data.material, data.strain_formulation,
            data.time, committed_state, history == nullptr ? nullptr : &(*history)[q], time_step, residual, *jacobian,
            average_strain_trace, average_trace_displacement_derivatives, element_pressure, geometry.reference_volume,
            finite_average_trace, finite_element_pressure);
    add_hex8_nodal_body_source_system(
        geometry, state, data.strain_formulation, data.volumetric_heat_source, residual, *jacobian);
    if (committed_state != nullptr && include_thermal_time_term)
        add_hex8_lumped_capacity_system(geometry, state, *committed_state, data.material, data.time, time_step,
            data.strain_formulation, residual, *jacobian);
    Hex8LocalResidual result{};
    ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    return result;
}
} // namespace

Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates) {
    Hex8Geometry geometry{};
    for (const CartesianPoint3& coordinate : coordinates) {
        geometry.selective_position.x += coordinate.x / 8.0;
        geometry.selective_position.y += coordinate.y / 8.0;
        geometry.selective_position.z += coordinate.z / 8.0;
    }
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
                const double determinant_value = cartesian_detail::determinant(jacobian);
                if (!std::isfinite(determinant_value) || !(determinant_value > 0.0))
                    throw std::invalid_argument("Hex8Geometry requires a finite positive Jacobian determinant");
                const cartesian_detail::Matrix3 inverse_jacobian =
                    cartesian_detail::inverse(jacobian, determinant_value);
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
    for (const Hex8QuadraturePoint& point : geometry.points) {
        geometry.reference_volume += point.weighted_measure;
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                geometry.average_shape_gradient[node][component] +=
                    point.weighted_measure * point.gradient[node][component];
    }
    if (!std::isfinite(geometry.reference_volume) || !(geometry.reference_volume > 0.0))
        throw std::invalid_argument("Hex8Geometry requires a finite positive reference volume");
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            geometry.average_shape_gradient[node][component] /= geometry.reference_volume;
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        geometry.capacity_points[node] = {
            coordinates[node], geometry.points[hex8_node_to_gauss[node]].weighted_measure};
        geometry.reduced_point.shape[node] = 1.0 / 8.0;
        geometry.reduced_point.gradient[node] = geometry.average_shape_gradient[node];
        geometry.reduced_capacity_points[node].position = coordinates[node];
        for (const Hex8QuadraturePoint& point : geometry.points)
            geometry.reduced_capacity_points[node].weighted_measure += point.weighted_measure * point.shape[node];
    }
    geometry.reduced_point.position = geometry.selective_position;
    geometry.reduced_point.weighted_measure = geometry.reference_volume;
    cartesian_detail::Matrix3 center_jacobian{};
    for (std::size_t physical = 0; physical < 3; ++physical)
        for (std::size_t natural = 0; natural < 3; ++natural)
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                const double coordinate = physical == 0   ? coordinates[node].x
                                          : physical == 1 ? coordinates[node].y
                                                          : coordinates[node].z;
                center_jacobian[physical][natural] += coordinate * hex8_signs[node][natural] / 8.0;
            }
    geometry.reduced_body_source_measure = 8.0 * cartesian_detail::determinant(center_jacobian);
    if (!std::isfinite(geometry.reduced_body_source_measure) || !(geometry.reduced_body_source_measure > 0.0))
        throw std::invalid_argument("Reduced HEX8 body-source integration requires a positive center Jacobian");

    constexpr std::array<std::array<double, 4>, hex8_node_count> raw_hourglass = {
        {{{1.0, -1.0, 1.0, -1.0}}, {{-1.0, -1.0, -1.0, 1.0}}, {{1.0, 1.0, -1.0, -1.0}}, {{-1.0, 1.0, 1.0, 1.0}},
            {{1.0, 1.0, -1.0, 1.0}}, {{-1.0, 1.0, 1.0, -1.0}}, {{1.0, -1.0, 1.0, 1.0}}, {{-1.0, -1.0, -1.0, -1.0}}}};
    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<double, 3> projected_coordinate{};
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            projected_coordinate[0] += coordinates[node].x * raw_hourglass[node][mode];
            projected_coordinate[1] += coordinates[node].y * raw_hourglass[node][mode];
            projected_coordinate[2] += coordinates[node].z * raw_hourglass[node][mode];
        }
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            geometry.hourglass_shape[node][mode] = raw_hourglass[node][mode];
            for (std::size_t component = 0; component < 3; ++component)
                geometry.hourglass_shape[node][mode] -=
                    geometry.average_shape_gradient[node][component] * projected_coordinate[component];
        }
    }

    cartesian_detail::Matrix3 inverse_effective_mapping{};
    for (std::size_t natural = 0; natural < 3; ++natural)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                inverse_effective_mapping[natural][physical] +=
                    hex8_signs[node][natural] * geometry.average_shape_gradient[node][physical];
    const double inverse_effective_determinant = cartesian_detail::determinant(inverse_effective_mapping);
    if (!std::isfinite(inverse_effective_determinant) || inverse_effective_determinant == 0.0)
        throw std::invalid_argument("Reduced HEX8 effective mapping must be nonsingular");
    const cartesian_detail::Matrix3 effective_mapping =
        cartesian_detail::inverse(inverse_effective_mapping, inverse_effective_determinant);
    cartesian_detail::Matrix3 metric{};
    for (std::size_t first = 0; first < 3; ++first)
        for (std::size_t second = 0; second < 3; ++second)
            for (std::size_t physical = 0; physical < 3; ++physical)
                metric[first][second] += effective_mapping[physical][first] * effective_mapping[physical][second];
    const double first_pivot = metric[0][0];
    const double second_pivot = metric[1][1] - metric[0][1] * metric[0][1] / first_pivot;
    const double leading_determinant = metric[0][0] * metric[1][1] - metric[0][1] * metric[0][1];
    const double third_pivot =
        metric[2][2] - (metric[1][1] * metric[0][2] * metric[0][2] - 2.0 * metric[0][1] * metric[0][2] * metric[1][2] +
                           metric[0][0] * metric[1][2] * metric[1][2]) /
                           leading_determinant;
    if (!std::isfinite(first_pivot) || !std::isfinite(second_pivot) || !std::isfinite(third_pivot) ||
        !(first_pivot > 0.0) || !(second_pivot > 0.0) || !(third_pivot > 0.0))
        throw std::invalid_argument("Reduced HEX8 effective metric must be positive definite");
    const double inverse_length_x_squared = 1.0 / first_pivot;
    const double inverse_length_y_squared = 1.0 / second_pivot;
    const double inverse_length_z_squared = 1.0 / third_pivot;
    const double thermal_scale = geometry.reference_volume / 192.0;
    geometry.thermal_hourglass_coefficients = {thermal_scale * (inverse_length_x_squared + inverse_length_y_squared),
        thermal_scale * (inverse_length_x_squared + inverse_length_z_squared),
        thermal_scale * (inverse_length_x_squared + inverse_length_z_squared),
        thermal_scale * (inverse_length_x_squared + inverse_length_y_squared + inverse_length_z_squared) / 3.0};
    for (std::size_t component = 0; component < 3; ++component) {
        double gradient_norm_squared = 0.0;
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            gradient_norm_squared +=
                geometry.average_shape_gradient[node][component] * geometry.average_shape_gradient[node][component];
        geometry.mechanical_hourglass_metrics[component] =
            std::sqrt(2.0) * geometry.reference_volume * gradient_norm_squared / 6.0;
    }
    return geometry;
}

Quad4FaceQuadraturePoint make_quad4_face_quadrature_point(
    const Quad4FaceCoordinates& coordinates, double xi, double eta, double quadrature_weight) {
    Quad4FaceQuadraturePoint point{};
    point.shape = {{0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)}};
    point.derivative_xi = {{-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)}};
    point.derivative_eta = {{-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)}};
    for (std::size_t node = 0; node < 4; ++node) {
        point.tangent_xi.x += point.derivative_xi[node] * coordinates[node].x;
        point.tangent_xi.y += point.derivative_xi[node] * coordinates[node].y;
        point.tangent_xi.z += point.derivative_xi[node] * coordinates[node].z;
        point.tangent_eta.x += point.derivative_eta[node] * coordinates[node].x;
        point.tangent_eta.y += point.derivative_eta[node] * coordinates[node].y;
        point.tangent_eta.z += point.derivative_eta[node] * coordinates[node].z;
    }
    point.normal_derivative_xi = point.derivative_xi;
    point.normal_derivative_eta = point.derivative_eta;
    point.normal_tangent_xi = point.tangent_xi;
    point.normal_tangent_eta = point.tangent_eta;
    const CartesianPoint3 area_vector{
        point.tangent_xi.y * point.tangent_eta.z - point.tangent_xi.z * point.tangent_eta.y,
        point.tangent_xi.z * point.tangent_eta.x - point.tangent_xi.x * point.tangent_eta.z,
        point.tangent_xi.x * point.tangent_eta.y - point.tangent_xi.y * point.tangent_eta.x};
    const double measure =
        std::sqrt(area_vector.x * area_vector.x + area_vector.y * area_vector.y + area_vector.z * area_vector.z);
    if (!std::isfinite(measure) || !(measure > 0.0) || !std::isfinite(quadrature_weight) || !(quadrature_weight > 0.0))
        throw std::invalid_argument("Quad4 face quadrature point requires finite positive measure and weight");
    point.weighted_measure = quadrature_weight * measure;
    return point;
}

Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates) {
    Quad4FaceGeometry geometry{};
    const std::array<std::array<double, 2>, 4> locations = {
        {{{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}}}};
    for (std::size_t q = 0; q < locations.size(); ++q)
        geometry.points[q] = make_quad4_face_quadrature_point(coordinates, locations[q][0], locations[q][1], 1.0);
    constexpr double nodal_area_location = 1.0 / 3.0;
    const std::array<std::array<double, 2>, 4> nodal_area_locations = {
        {{{-nodal_area_location, -nodal_area_location}}, {{nodal_area_location, -nodal_area_location}},
            {{nodal_area_location, nodal_area_location}}, {{-nodal_area_location, nodal_area_location}}}};
    for (std::size_t q = 0; q < nodal_area_locations.size(); ++q) {
        Quad4FaceQuadraturePoint& point = geometry.thermal_points[q];
        point =
            make_quad4_face_quadrature_point(coordinates, nodal_area_locations[q][0], nodal_area_locations[q][1], 1.0);
        point.shape.fill(0.0);
        point.shape[q] = 1.0;
        const Quad4FaceQuadraturePoint& normal_point = geometry.points[q];
        point.normal_derivative_xi = normal_point.derivative_xi;
        point.normal_derivative_eta = normal_point.derivative_eta;
        point.normal_tangent_xi = normal_point.tangent_xi;
        point.normal_tangent_eta = normal_point.tangent_eta;
        const CartesianPoint3 area{point.tangent_xi.y * point.tangent_eta.z - point.tangent_xi.z * point.tangent_eta.y,
            point.tangent_xi.z * point.tangent_eta.x - point.tangent_xi.x * point.tangent_eta.z,
            point.tangent_xi.x * point.tangent_eta.y - point.tangent_xi.y * point.tangent_eta.x};
        const CartesianPoint3 normal_area{normal_point.tangent_xi.y * normal_point.tangent_eta.z -
                                              normal_point.tangent_xi.z * normal_point.tangent_eta.y,
            normal_point.tangent_xi.z * normal_point.tangent_eta.x -
                normal_point.tangent_xi.x * normal_point.tangent_eta.z,
            normal_point.tangent_xi.x * normal_point.tangent_eta.y -
                normal_point.tangent_xi.y * normal_point.tangent_eta.x};
        const double normal_measure =
            std::sqrt(normal_area.x * normal_area.x + normal_area.y * normal_area.y + normal_area.z * normal_area.z);
        point.weighted_measure =
            (area.x * normal_area.x + area.y * normal_area.y + area.z * normal_area.z) / normal_measure;
        if (!std::isfinite(point.weighted_measure) || !(point.weighted_measure > 0.0))
            throw std::invalid_argument("Quad4 face thermal nodal integration requires a positive projected area");
    }
    return geometry;
}

Hex8LocalResidual compute_hex8_thermoelastic(const CartesianThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state, double time_step,
    Hex8LocalJacobian* jacobian) {
    return compute_hex8_local(data, geometry, state, committed_state, nullptr, time_step, jacobian, true);
}

Hex8LocalResidual compute_hex8_transient(const CartesianThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material, double time_step, Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    return compute_hex8_local(
        data, geometry, state, &committed_state, &committed_material, time_step, jacobian, include_thermal_time_term);
}

CartesianMaterialHistory compute_hex8_transient_update(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry, const Hex8LocalValues& state, const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material, double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("HEX8 transient update time step must be finite and positive");
    if (data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt) {
        if (committed_material.size() != 1)
            throw std::invalid_argument("C3D8RT material history must contain one integration point");
        Hex8LocalAdValues passive{};
        ad_local_system::make_passive(state.data(), state.size(), passive.data());
        const Hex8QuadraturePoint& point = geometry.reduced_point;
        const MaterialFunctionContext context = material_context(data.time, point.position);
        CartesianInelasticStressResponse response;
        if (data.strain_formulation == StrainFormulation::finite) {
            const ActiveReducedHex8Geometry current =
                active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive), "current");
            const ActiveReducedHex8Geometry midpoint = active_reduced_hex8_geometry(
                geometry, reduced_hex8_midpoint_displacement(passive, committed_state), "midpoint");
            const cartesian_detail::KinematicsCore kinematics =
                reduced_hex8_finite_kinematics(midpoint, passive, committed_state);
            const adlite::Scalar temperature = reduced_hex8_temperature(current, passive);
            Hex8LocalAdValues passive_old{};
            ad_local_system::make_passive(committed_state.data(), committed_state.size(), passive_old.data());
            const ActiveReducedHex8Geometry old_geometry =
                active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive_old), "committed");
            const double old_temperature = reduced_hex8_temperature(old_geometry, passive_old).value();
            response = data.material.incremental_response(kinematics.strain_increment, kinematics.rotation, temperature,
                old_temperature, time_step, committed_material.front(), context);
        } else {
            const adlite::Scalar temperature = reduced_hex8_temperature(geometry, passive);
            const CartesianKinematics kinematics =
                evaluate_cartesian_incremental_kinematics(point, passive, committed_state, StrainFormulation::small);
            response = data.material.response(
                kinematics.strain_increment, temperature, time_step, committed_material.front(), context);
        }
        return CartesianMaterialHistory{response.trial_state};
    }
    if (committed_material.size() != geometry.points.size())
        throw std::invalid_argument("C3D8T material history must contain eight integration points");
    Hex8LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    const adlite::Scalar expansion_temperature = average_hex8_temperature(passive);
    const double old_expansion_temperature = average_hex8_temperature(committed_state);
    const double average_strain_trace = data.strain_formulation == StrainFormulation::finite
                                            ? finite_average_hex8_strain_trace(geometry, state, committed_state)
                                            : average_hex8_strain_trace(geometry, state);
    CartesianMaterialHistory result(geometry.points.size());
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const std::size_t material_node = hex8_node_to_gauss[q];
        const adlite::Scalar temperature = passive[material_node];
        const CartesianKinematics kinematics =
            evaluate_cartesian_incremental_kinematics(point, passive, committed_state, data.strain_formulation);
        CartesianInelasticStressResponse response;
        if (data.strain_formulation == StrainFormulation::finite) {
            const double old_temperature = committed_state[material_node];
            const MaterialFunctionContext context = material_context(data.time, point.position);
            response = data.material.incremental_response(
                expansion_adjusted_increment(data.material,
                    selectively_reduced_strain(
                        kinematics.strain_increment, average_strain_trace, data.strain_formulation),
                    temperature, expansion_temperature, old_temperature, old_expansion_temperature, time_step, context),
                kinematics.rotation, temperature, old_temperature, time_step, committed_material[q], context);
        } else {
            const MaterialFunctionContext context = material_context(data.time, point.position);
            const SymmetricTensor3 constitutive_strain =
                selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, data.strain_formulation);
            response = data.material.response(expansion_adjusted_strain(data.material, constitutive_strain, temperature,
                                                  expansion_temperature, context),
                temperature, time_step, committed_material[q], context);
        }
        result[q] = response.trial_state;
    }
    return result;
}

double compute_hex8_mechanical_hourglass_energy(
    const CartesianThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    if (data.hex8_element_formulation != Hex8ElementFormulation::c3d8rt) return 0.0;
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    const ActiveThermoelasticProperties initial_properties = data.material.active_properties(
        adlite::Scalar(data.initial_temperature), material_context(0.0, geometry.reduced_point.position));
    const double initial_shear_modulus = initial_properties.shear_modulus.value();
    if (!std::isfinite(initial_shear_modulus) || !(initial_shear_modulus > 0.0))
        throw std::invalid_argument("C3D8RT initial shear modulus must be finite and positive");

    constexpr double abaqus_total_stiffness_factor = 0.005;
    double energy = 0.0;
    if (data.strain_formulation == StrainFormulation::small) {
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t mode = 0; mode < 4; ++mode) {
                double amplitude = 0.0;
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    amplitude += geometry.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
                const double stiffness = abaqus_total_stiffness_factor * initial_shear_modulus *
                                         geometry.mechanical_hourglass_metrics[component];
                energy += 0.5 * stiffness * amplitude * amplitude;
            }
        return energy;
    }
    if (data.strain_formulation != StrainFormulation::finite)
        throw std::invalid_argument("C3D8RT hourglass energy requires small or finite strain");

    std::array<std::array<double, 3>, 3> average_deformation{};
    for (std::size_t component = 0; component < 3; ++component) {
        average_deformation[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                average_deformation[component][direction] +=
                    state[8 * (component + 1) + node] * geometry.average_shape_gradient[node][direction];
    }
    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<double, 3> reference_amplitude{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                reference_amplitude[component] +=
                    geometry.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
        for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
            double transported_amplitude = 0.0;
            for (std::size_t component = 0; component < 3; ++component)
                transported_amplitude +=
                    reference_amplitude[component] * average_deformation[component][material_direction];
            const double stiffness = abaqus_total_stiffness_factor * initial_shear_modulus *
                                     geometry.mechanical_hourglass_metrics[material_direction];
            energy += 0.5 * stiffness * transported_amplitude * transported_amplitude;
        }
    }
    return energy;
}

std::array<SymmetricTensor3Values, 8> compute_hex8_stress(
    const CartesianThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    if (data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt) {
        Hex8LocalAdValues passive{};
        ad_local_system::make_passive(state.data(), state.size(), passive.data());
        SymmetricTensor3 stress;
        if (data.strain_formulation == StrainFormulation::finite) {
            const Hex8LocalValues undeformed{};
            const ActiveReducedHex8Geometry current =
                active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive), "current");
            const ActiveReducedHex8Geometry midpoint = active_reduced_hex8_geometry(
                geometry, reduced_hex8_midpoint_displacement(passive, undeformed), "midpoint");
            const cartesian_detail::KinematicsCore kinematics =
                reduced_hex8_finite_kinematics(midpoint, passive, undeformed);
            const adlite::Scalar temperature = reduced_hex8_temperature(current, passive);
            stress = rotate_cartesian_tensor(data.material.stress(kinematics.strain_increment, temperature,
                                                 material_context(data.time, geometry.reduced_point.position)),
                kinematics.rotation);
        } else {
            const adlite::Scalar temperature = reduced_hex8_temperature(geometry, passive);
            const CartesianKinematics kinematics = evaluate_cartesian_incremental_kinematics(
                geometry.reduced_point, passive, Hex8LocalValues{}, StrainFormulation::small);
            stress = data.material.stress(
                kinematics.strain_increment, temperature, material_context(data.time, geometry.reduced_point.position));
        }
        const SymmetricTensor3Values values{stress.xx.value(), stress.yy.value(), stress.zz.value(), stress.xy.value(),
            stress.yz.value(), stress.xz.value()};
        std::array<SymmetricTensor3Values, 8> result{};
        result.fill(values);
        return result;
    }
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
    const auto& points =
        data.kind == Quad4FaceBoundaryKind::surface_heat_flux || data.kind == Quad4FaceBoundaryKind::convection
            ? geometry.thermal_points
            : geometry.points;
    for (const Quad4FaceQuadraturePoint& point : points) {
        std::array<adlite::Scalar, 3> tangent_xi = {point.tangent_xi.x, point.tangent_xi.y, point.tangent_xi.z};
        std::array<adlite::Scalar, 3> tangent_eta = {point.tangent_eta.x, point.tangent_eta.y, point.tangent_eta.z};
        std::array<adlite::Scalar, 3> normal_tangent_xi = {
            point.normal_tangent_xi.x, point.normal_tangent_xi.y, point.normal_tangent_xi.z};
        std::array<adlite::Scalar, 3> normal_tangent_eta = {
            point.normal_tangent_eta.x, point.normal_tangent_eta.y, point.normal_tangent_eta.z};
        if (data.use_displaced_geometry)
            for (std::size_t node = 0; node < 4; ++node)
                for (std::size_t component = 0; component < 3; ++component) {
                    tangent_xi[component] += point.derivative_xi[node] * ad_state[4 * (component + 1) + node];
                    tangent_eta[component] += point.derivative_eta[node] * ad_state[4 * (component + 1) + node];
                    normal_tangent_xi[component] +=
                        point.normal_derivative_xi[node] * ad_state[4 * (component + 1) + node];
                    normal_tangent_eta[component] +=
                        point.normal_derivative_eta[node] * ad_state[4 * (component + 1) + node];
                }
        const std::array<adlite::Scalar, 3> area = {tangent_xi[1] * tangent_eta[2] - tangent_xi[2] * tangent_eta[1],
            tangent_xi[2] * tangent_eta[0] - tangent_xi[0] * tangent_eta[2],
            tangent_xi[0] * tangent_eta[1] - tangent_xi[1] * tangent_eta[0]};
        const adlite::Scalar raw_measure = adlite::hypot(adlite::hypot(area[0], area[1]), area[2]);
        if (!std::isfinite(raw_measure.value()) || !(raw_measure.value() > 0.0))
            throw std::domain_error("Three-dimensional face requires a positive current measure");
        const bool thermal =
            data.kind == Quad4FaceBoundaryKind::surface_heat_flux || data.kind == Quad4FaceBoundaryKind::convection;
        adlite::Scalar measure = raw_measure;
        if (thermal) {
            if (!data.use_displaced_geometry)
                measure = point.weighted_measure;
            else {
                const std::array<adlite::Scalar, 3> normal_area = {
                    normal_tangent_xi[1] * normal_tangent_eta[2] - normal_tangent_xi[2] * normal_tangent_eta[1],
                    normal_tangent_xi[2] * normal_tangent_eta[0] - normal_tangent_xi[0] * normal_tangent_eta[2],
                    normal_tangent_xi[0] * normal_tangent_eta[1] - normal_tangent_xi[1] * normal_tangent_eta[0]};
                const adlite::Scalar normal_measure =
                    adlite::hypot(adlite::hypot(normal_area[0], normal_area[1]), normal_area[2]);
                measure =
                    (area[0] * normal_area[0] + area[1] * normal_area[1] + area[2] * normal_area[2]) / normal_measure;
                if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
                    throw std::domain_error(
                        "Three-dimensional thermal face requires a positive current projected area");
            }
        }
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
        } else if (data.kind == Quad4FaceBoundaryKind::surface_heat_flux) {
            for (std::size_t node = 0; node < 4; ++node) residual[node] -= data.load * measure * point.shape[node];
        } else {
            adlite::Scalar temperature = 0.0;
            for (std::size_t node = 0; node < 4; ++node) temperature += point.shape[node] * ad_state[node];
            const adlite::Scalar heat_flux = data.load * (temperature - data.ambient_temperature);
            for (std::size_t node = 0; node < 4; ++node) residual[node] += measure * point.shape[node] * heat_flux;
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
