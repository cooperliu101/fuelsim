#include "c3d20_assembly.hpp"
#include "ad_local_system.hpp"
#include "c3d20_kinematics.hpp"
#include "c3d20_types.hpp"
#include "cartesian_material.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim::c3d20_detail {

namespace {
using cartesian_detail::ActiveMatrix3;
using cartesian_detail::determinant;
using cartesian_detail::inverse;
using cartesian_detail::material_context;
using cartesian_detail::Matrix3;

adlite::Scalar interpolate_temperature(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        result += point.temperature_shape[node] * state[node];
    return result;
}

double interpolate_temperature_values(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        result += point.temperature_shape[node] * state[node];
    return result;
}

double interpolate_temperature_values(const Hex20ThermalQuadraturePoint& point, const Hex20LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        result += point.temperature_shape[node] * state[node];
    return result;
}

void add_thermal_point_residual_values(const Hex20ThermalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double volumetric_heat_source,
    const Hex20LocalValues* committed_state,
    double time_step,
    Hex20LocalResidual& residual,
    bool include_thermal_time_term) {
    const double temperature = interpolate_temperature_values(point, state);
    std::array<double, 3> temperature_gradient{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += point.temperature_gradient[node][direction] * state[node];
    const MaterialFunctionContext context = material_context(time, point.position);
    const double conductivity = material.conductivity(adlite::Scalar(temperature), context).value();
    double temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr && include_thermal_time_term) {
        const double old_temperature = interpolate_temperature_values(point, *committed_state);
        temperature_rate = (temperature - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(adlite::Scalar(temperature), context).value();
    }
    for (std::size_t node = 0; node < 8; ++node) {
        double conduction = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            conduction += point.temperature_gradient[node][direction] * temperature_gradient[direction];
        residual[node] +=
            point.weighted_measure
            * (conductivity * conduction + point.temperature_shape[node] * heat_capacity * temperature_rate
                - point.temperature_shape[node] * volumetric_heat_source);
    }
}

void add_mechanical_point_residual_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    const Hex20LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    Hex20LocalResidual& residual) {
    const double temperature = interpolate_temperature_values(point, state);
    const MaterialFunctionContext context = material_context(time, point.position);
    const Hex20LocalValues undeformed{};
    const Hex20LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const Hex20KinematicsValues kinematics = evaluate_kinematics_values(point, state, old_state, strain_formulation);
    const SymmetricTensor3Values strain{kinematics.strain_increment.xx,
        kinematics.strain_increment.yy,
        kinematics.strain_increment.zz,
        kinematics.strain_increment.xy,
        kinematics.strain_increment.yz,
        kinematics.strain_increment.xz};
    SymmetricTensor3Values stress{};
    if (committed_material == nullptr) {
        stress = material.stress_values(strain, temperature, context);
        if (strain_formulation == StrainFormulation::finite)
            stress = cartesian_detail::rotate_cartesian_tensor_values(stress, kinematics.rotation);
    } else if (strain_formulation == StrainFormulation::finite) {
        const double old_temperature = interpolate_temperature_values(point, old_state);
        stress = material
                     .incremental_response_values(strain,
                         kinematics.rotation,
                         temperature,
                         old_temperature,
                         time_step,
                         *committed_material,
                         context)
                     .stress;
    } else {
        stress = material.response_values(strain, temperature, time_step, *committed_material, context).stress;
    }
    for (std::size_t node = 0; node < 20; ++node) {
        const double gx = kinematics.current_gradient[node][0], gy = kinematics.current_gradient[node][1],
                     gz = kinematics.current_gradient[node][2];
        residual[8 + node] += kinematics.current_weighted_measure * (stress.xx * gx + stress.xy * gy + stress.xz * gz);
        residual[28 + node] += kinematics.current_weighted_measure * (stress.xy * gx + stress.yy * gy + stress.yz * gz);
        residual[48 + node] += kinematics.current_weighted_measure * (stress.xz * gx + stress.yz * gy + stress.zz * gz);
    }
}

void add_thermal_point_system(const Hex20ThermalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double volumetric_heat_source,
    const Hex20LocalValues* committed_state,
    double time_step,
    Hex20LocalResidual& residual,
    Hex20LocalJacobian& jacobian,
    bool include_thermal_time_term) {
    constexpr std::size_t point_width = 1, temperature_index = 0;
    std::array<adlite::Scalar, 1> active_temperature{};
    double temperature_value = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        temperature_value += point.temperature_shape[node] * state[node];
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
    std::array<double, point_width> derivatives{};
    for (std::size_t node = 0; node < 8; ++node) {
        adlite::Scalar conduction = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            conduction += point.temperature_gradient[node][direction] * temperature_gradient[direction];
        const adlite::Scalar point_residual =
            point.weighted_measure
            * (conductivity * conduction + point.temperature_shape[node] * heat_capacity * temperature_rate
                - point.temperature_shape[node] * volumetric_heat_source);
        point_residual.copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t other = 0; other < 8; ++other) {
            double gradient_dot = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                gradient_dot +=
                    point.temperature_gradient[node][direction] * point.temperature_gradient[other][direction];
            jacobian[node * hex20_local_dof_count + other] +=
                point.weighted_measure * conductivity.value() * gradient_dot
                + point.temperature_shape[other] * derivatives[temperature_index];
        }
        residual[node] += point_residual.value();
    }
}

void add_finite_thermal_point_residual_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double volumetric_heat_source,
    const Hex20LocalValues* committed_state,
    double time_step,
    Hex20LocalResidual& residual,
    bool include_thermal_time_term) {
    const Hex20LocalValues undeformed{};
    const Hex20LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const Hex20FiniteThermalKinematicsValues kinematics =
        evaluate_finite_thermal_kinematics_values(point, state, old_state);
    const double temperature = interpolate_temperature_values(point, state);
    std::array<double, 3> temperature_gradient{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += kinematics.midpoint_temperature_gradient[node][direction] * state[node];
    const MaterialFunctionContext context = material_context(time, point.position);
    const double conductivity = material.conductivity(adlite::Scalar(temperature), context).value();
    double temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr && include_thermal_time_term) {
        const double old_temperature = interpolate_temperature_values(point, *committed_state);
        temperature_rate = (temperature - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(adlite::Scalar(temperature), context).value();
    }
    const double source_measure =
        volumetric_heat_source == 0.0 ? 0.0 : evaluate_source_measure_values(point, state).weighted_measure;
    for (std::size_t node = 0; node < 8; ++node) {
        double conduction = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            conduction += kinematics.midpoint_temperature_gradient[node][direction] * temperature_gradient[direction];
        residual[node] +=
            kinematics.current_weighted_measure
                * (conductivity * conduction + point.temperature_shape[node] * heat_capacity * temperature_rate)
            - source_measure * point.temperature_shape[node] * volumetric_heat_source;
    }
}

void add_finite_thermal_point_system(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double volumetric_heat_source,
    const Hex20LocalValues* committed_state,
    double time_step,
    const Hex20Kinematics& kinematics,
    const adlite::Scalar& active_temperature,
    Hex20LocalResidual& residual,
    Hex20LocalJacobian& jacobian,
    bool include_thermal_time_term) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    const MaterialFunctionContext context = material_context(time, point.position);
    const adlite::Scalar conductivity = material.conductivity(active_temperature, context);
    adlite::Scalar temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr && include_thermal_time_term) {
        const double old_temperature = interpolate_temperature_values(point, *committed_state);
        temperature_rate = (active_temperature - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(active_temperature, context);
    }
    const auto midpoint_gradient = temperature_shape_gradients(point, kinematics.midpoint_inverse_values);
    std::array<double, 3> temperature_gradient{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += midpoint_gradient[node][direction] * state[node];
    Hex20SourceMeasureValues source;
    if (volumetric_heat_source != 0.0)
        source = evaluate_source_measure_values(point, state);
    std::array<std::array<double, 8>, 3> source_displacement_derivative{};
    if (volumetric_heat_source != 0.0)
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t other = 0; other < 8; ++other)
                for (std::size_t direction = 0; direction < 3; ++direction)
                    source_displacement_derivative[component][other] +=
                        source.measure_derivative[component][direction]
                        * point.source_displacement_gradient[other][direction];
    std::array<double, 3> midpoint_projected_temperature{};
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t direction = 0; direction < 3; ++direction)
            midpoint_projected_temperature[row] +=
                kinematics.midpoint_inverse_values[row][direction] * temperature_gradient[direction];

    std::array<double, point_width> derivatives{};
    for (std::size_t node = 0; node < 8; ++node) {
        double conduction = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            conduction += midpoint_gradient[node][direction] * temperature_gradient[direction];
        const adlite::Scalar thermal_integrand =
            conductivity * conduction + point.temperature_shape[node] * heat_capacity * temperature_rate;
        thermal_integrand.copy_derivatives(derivatives.data(), derivatives.size());
        const double current_weight = kinematics.current_weighted_measure.value();
        std::array<double, 3> midpoint_projected_test{};
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t direction = 0; direction < 3; ++direction)
                midpoint_projected_test[row] +=
                    midpoint_gradient[node][direction] * kinematics.midpoint_inverse_values[row][direction];
        std::array<std::array<double, 3>, 3> full_gradient_derivative{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                const double conduction_derivative =
                    -0.5
                    * (midpoint_gradient[node][component] * midpoint_projected_temperature[direction]
                        + temperature_gradient[component] * midpoint_projected_test[direction]);
                full_gradient_derivative[component][direction] =
                    current_weight
                    * (kinematics.current_inverse_values[direction][component] * thermal_integrand.value()
                        + conductivity.value() * conduction_derivative);
            }
        for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component)
            for (std::size_t other = 0; other < 20; ++other) {
                double chained = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    chained += full_gradient_derivative[displacement_component][direction]
                               * point.displacement_gradient[other][direction];
                if (other < 8)
                    chained -= point.temperature_shape[node] * volumetric_heat_source
                               * source_displacement_derivative[displacement_component][other];
                jacobian[node * hex20_local_dof_count + 8 + 20 * displacement_component + other] += chained;
            }
        for (std::size_t other = 0; other < 8; ++other) {
            double gradient_dot = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                gradient_dot += midpoint_gradient[node][direction] * midpoint_gradient[other][direction];
            jacobian[node * hex20_local_dof_count + other] +=
                current_weight
                * (conductivity.value() * gradient_dot
                    + point.temperature_shape[other] * derivatives[temperature_index]);
        }
        residual[node] += current_weight * thermal_integrand.value()
                          - source.weighted_measure * point.temperature_shape[node] * volumetric_heat_source;
    }
}

void add_mechanical_point_system(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    double volumetric_heat_source,
    const Hex20LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    Hex20LocalResidual& residual,
    Hex20LocalJacobian& jacobian,
    bool include_thermal_time_term) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    std::array<double, 9> gradient_values{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 20; ++node)
                gradient_values[component * 3 + direction] +=
                    point.displacement_gradient[node][direction] * state[8 + 20 * component + node];
    double temperature_value = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        temperature_value += point.temperature_shape[node] * state[node];
    ActiveMatrix3 active_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            active_gradient[component][direction] =
                adlite::Scalar::independent(gradient_values[component * 3 + direction],
                    component * 3 + direction,
                    point_width);
    const adlite::Scalar active_temperature =
        adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    const Hex20LocalValues undeformed{};
    const Hex20LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const Hex20Kinematics kinematics = evaluate_kinematics(point, active_gradient, old_state, strain_formulation);
    if (strain_formulation == StrainFormulation::finite)
        add_finite_thermal_point_system(point,
            state,
            material,
            time,
            volumetric_heat_source,
            committed_state,
            time_step,
            kinematics,
            active_temperature,
            residual,
            jacobian,
            include_thermal_time_term);
    const MaterialFunctionContext context = material_context(time, point.position);
    const std::array<const adlite::Scalar*, 6> strain_components = {&kinematics.strain_increment.xx,
        &kinematics.strain_increment.yy,
        &kinematics.strain_increment.zz,
        &kinematics.strain_increment.xy,
        &kinematics.strain_increment.yz,
        &kinematics.strain_increment.xz};
    std::array<double, 6> fed_strain{};
    if (committed_material != nullptr && strain_formulation == StrainFormulation::finite) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node)
            old_temperature += point.temperature_shape[node] * old_state[node];
        if (!std::isfinite(old_temperature) || !(old_temperature > 0.0))
            throw std::domain_error("Incremental HEX20 material committed temperature must be finite and positive");
        MaterialFunctionContext old_context = context;
        old_context.time -= time_step;
        const SymmetricTensor3 old_imposed = material.eigenstrain(adlite::Scalar(old_temperature), old_context);
        const std::array<double, 6> imposed = {old_imposed.xx.value(),
            old_imposed.yy.value(),
            old_imposed.zz.value(),
            old_imposed.xy.value(),
            old_imposed.yz.value(),
            old_imposed.xz.value()};
        for (std::size_t component = 0; component < 6; ++component)
            fed_strain[component] = committed_material->elastic_strain[component]
                                    + strain_components[component]->value() + imposed[component]
                                    + committed_material->plastic_strain[component]
                                    + committed_material->creep_strain[component];
    } else {
        for (std::size_t component = 0; component < 6; ++component)
            fed_strain[component] = strain_components[component]->value();
    }
    const cartesian_detail::CartesianStressTangent tangent = cartesian_detail::evaluate_stress_tangent(material,
        fed_strain,
        temperature_value,
        time_step,
        committed_material,
        context);
    std::array<adlite::Scalar, 7> compose_inputs{};
    for (std::size_t component = 0; component < 6; ++component)
        compose_inputs[component] = *strain_components[component];
    compose_inputs[6] = active_temperature;
    const std::array<double, 6> stress_values = {tangent.stress.xx,
        tangent.stress.yy,
        tangent.stress.zz,
        tangent.stress.xy,
        tangent.stress.yz,
        tangent.stress.xz};
    std::array<double, 7> partials{};
    std::array<adlite::Scalar, 6> composed{};
    for (std::size_t component = 0; component < 6; ++component) {
        for (std::size_t column = 0; column < 6; ++column)
            partials[column] = tangent.tangent[component][column];
        partials[6] = tangent.thermal[component];
        composed[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    SymmetricTensor3 stress{composed[0], composed[1], composed[2], composed[3], composed[4], composed[5]};
    if (strain_formulation == StrainFormulation::finite)
        stress = rotate_cartesian_tensor(stress, kinematics.rotation);
    std::array<adlite::Scalar, 60> point_residual{};
    point_residual.fill(adlite::Scalar(0.0));
    for (std::size_t node = 0; node < 20; ++node) {
        const adlite::Scalar gx = kinematics.current_gradient[node][0], gy = kinematics.current_gradient[node][1],
                             gz = kinematics.current_gradient[node][2];
        point_residual[node] +=
            kinematics.current_weighted_measure * (stress.xx * gx + stress.xy * gy + stress.xz * gz);
        point_residual[20 + node] +=
            kinematics.current_weighted_measure * (stress.xy * gx + stress.yy * gy + stress.yz * gz);
        point_residual[40 + node] +=
            kinematics.current_weighted_measure * (stress.xz * gx + stress.yz * gy + stress.zz * gz);
    }
    std::array<double, point_width> derivatives{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 20; ++node) {
            const std::size_t row = 8 + 20 * component + node;
            point_residual[20 * component + node].copy_derivatives(derivatives.data(), derivatives.size());
            for (std::size_t other = 0; other < 8; ++other)
                jacobian[row * hex20_local_dof_count + other] +=
                    derivatives[temperature_index] * point.temperature_shape[other];
            for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component)
                for (std::size_t other = 0; other < 20; ++other) {
                    double chained = 0.0;
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        chained += derivatives[displacement_component * 3 + direction]
                                   * point.displacement_gradient[other][direction];
                    jacobian[row * hex20_local_dof_count + 8 + 20 * displacement_component + other] += chained;
                }
            residual[row] += point_residual[20 * component + node].value();
        }
}

Hex20LocalResidual compute_local(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    Hex20LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    if (committed_state != nullptr && (!std::isfinite(time_step) || !(time_step > 0.0)))
        throw std::invalid_argument("HEX20 time step must be finite and positive");
    if (history != nullptr && history->size() != geometry.mechanical_points.size())
        throw std::invalid_argument("HEX20 material history must match the active integration point count");
    if (jacobian == nullptr) {
        Hex20LocalResidual result{};
        if (data.strain_formulation == StrainFormulation::finite)
            for (const Hex20MechanicalQuadraturePoint& point : geometry.mechanical_points)
                add_finite_thermal_point_residual_values(point,
                    state,
                    data.material,
                    data.time,
                    data.volumetric_heat_source,
                    committed_state,
                    time_step,
                    result,
                    include_thermal_time_term);
        else
            for (const Hex20ThermalQuadraturePoint& point : geometry.thermal_points)
                add_thermal_point_residual_values(point,
                    state,
                    data.material,
                    data.time,
                    data.volumetric_heat_source,
                    committed_state,
                    time_step,
                    result,
                    include_thermal_time_term);
        for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q)
            add_mechanical_point_residual_values(geometry.mechanical_points[q],
                state,
                data.material,
                data.strain_formulation,
                data.time,
                committed_state,
                history == nullptr ? nullptr : &(*history)[q],
                time_step,
                result);
        return result;
    }
    Hex20LocalResidual residual{};
    jacobian->fill(0.0);
    if (data.strain_formulation != StrainFormulation::finite)
        for (const Hex20ThermalQuadraturePoint& point : geometry.thermal_points)
            add_thermal_point_system(point,
                state,
                data.material,
                data.time,
                data.volumetric_heat_source,
                committed_state,
                time_step,
                residual,
                *jacobian,
                include_thermal_time_term);
    for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q)
        add_mechanical_point_system(geometry.mechanical_points[q],
            state,
            data.material,
            data.strain_formulation,
            data.time,
            data.volumetric_heat_source,
            committed_state,
            history == nullptr ? nullptr : &(*history)[q],
            time_step,
            residual,
            *jacobian,
            include_thermal_time_term);
    return residual;
}

Hex20LocalResidual compute_hex20_thermoelastic(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues* committed_state,
    double time_step,
    Hex20LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    return compute_local(data,
        geometry,
        state,
        committed_state,
        nullptr,
        time_step,
        jacobian,
        include_thermal_time_term);
}

Hex20LocalResidual compute_hex20_transient(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex20LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    return compute_local(data,
        geometry,
        state,
        &committed_state,
        &committed_material,
        time_step,
        jacobian,
        include_thermal_time_term);
}

CartesianMaterialHistory compute_hex20_transient_update(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("HEX20 transient update time step must be finite and positive");
    if (committed_material.size() != geometry.mechanical_points.size())
        throw std::invalid_argument("HEX20 material history must match the active integration point count");
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
            response = data.material.incremental_response(kinematics.strain_increment,
                kinematics.rotation,
                temperature,
                old_temperature,
                time_step,
                committed_material[q],
                material_context(data.time, point.position));
        } else {
            response = data.material.response(kinematics.strain_increment,
                temperature,
                time_step,
                committed_material[q],
                material_context(data.time, point.position));
        }
        result[q] = response.trial_state;
    }
    return result;
}

std::vector<SymmetricTensor3Values>
compute_hex20_stress(const elements::C3d20Input& data, const Hex20Geometry& geometry, const Hex20LocalValues& state) {
    Hex20LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    std::vector<SymmetricTensor3Values> result(geometry.mechanical_points.size());
    for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q) {
        const Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[q];
        const adlite::Scalar temperature = interpolate_temperature(point, passive);
        const Hex20Kinematics kinematics =
            evaluate_kinematics(point, passive, Hex20LocalValues{}, data.strain_formulation);
        SymmetricTensor3 stress =
            data.material.stress(kinematics.strain_increment, temperature, material_context(data.time, point.position));
        if (data.strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
        result[q] = {stress.xx.value(),
            stress.yy.value(),
            stress.zz.value(),
            stress.xy.value(),
            stress.yz.value(),
            stress.xz.value()};
    }
    return result;
}

} // namespace

elements::C3d20Result evaluate(const elements::C3d20Input& input, elements::ElementRequest request) {
    const auto& data = input;
    elements::C3d20Result result;
    if (request.residual || request.jacobian) {
        auto* tangent = request.jacobian ? &result.jacobian : nullptr;
        if (input.committed_history)
            result.residual = compute_hex20_transient(data,
                input.geometry,
                input.state,
                input.committed_state,
                *input.committed_history,
                input.time_step,
                tangent,
                input.include_thermal_time_term);
        else
            result.residual = compute_hex20_thermoelastic(data,
                input.geometry,
                input.state,
                input.time_step > 0 ? &input.committed_state : nullptr,
                input.time_step,
                tangent,
                input.include_thermal_time_term);
    }
    if (request.history && input.committed_history)
        result.history = compute_hex20_transient_update(data,
            input.geometry,
            input.state,
            input.committed_state,
            *input.committed_history,
            input.time_step);
    if (request.stress)
        result.stress = compute_hex20_stress(data, input.geometry, input.state);
    return result;
}
} // namespace fuelsim::c3d20_detail
