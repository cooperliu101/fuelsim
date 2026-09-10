#include "c3d8t.hpp"
#include "ad_local_system.hpp"
#include "c3d_common.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {
using namespace c3d8_detail;
using namespace cartesian_detail;

struct FiniteTracePoint final {
    adlite::Scalar strain_trace;
    adlite::Scalar midpoint_weighted_measure;
    adlite::Scalar current_weighted_measure;
};

struct FiniteAverageTraceValues final {
    double value = 0.0;
    double current_volume = 0.0;
};

struct FiniteAverageTraceSystem final {
    double value = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> displacement_derivatives{};
    double current_volume = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> current_volume_derivatives{};
};

struct SmallStrainElementPressureSystem final {
    double value = 0.0, trace_derivative = 0.0;
    std::array<double, hex8_node_count> temperature_derivatives{};
};

struct FiniteElementPressureSystem final {
    double value = 0.0;
    std::array<double, hex8_local_dof_count> derivatives{};
};

struct FinitePointSystemCache final {
    adlite::Scalar active_temperature;
    C3d8Kinematics kinematics;
    cartesian_detail::CartesianStressTangent tangent;
};

struct FinitePointResidualCache final {
    cartesian_detail::Matrix3 displacement_gradient{};
    cartesian_detail::Matrix3 current_deformation{};
    cartesian_detail::Matrix3 rotation{};
    std::array<std::array<double, 3>, hex8_node_count> current_gradient{};
    SymmetricTensor3Values strain_increment{};
    SymmetricTensor3Values stress{};
    double current_weighted_measure = 0.0;
};

struct Hex8NodalVolumeSystem final {
    double value;
    std::array<double, 9> gradient_derivatives;
};

adlite::Scalar average_hex8_temperature(const Hex8LocalAdValues& state);
double average_hex8_temperature(const Hex8LocalValues& state);
double average_hex8_strain_trace(const Hex8Geometry& geometry, const Hex8LocalValues& state);
FiniteTracePoint finite_trace_point(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state);
FiniteAverageTraceValues finite_average_hex8_strain_trace_values(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state);
double finite_average_hex8_strain_trace(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state);
SymmetricTensor3
selectively_reduced_strain(const SymmetricTensor3& strain, double average_trace, StrainFormulation strain_formulation);
SymmetricTensor3 expansion_adjusted_strain(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain,
    const adlite::Scalar& point_temperature,
    const adlite::Scalar& element_temperature,
    MaterialFunctionContext context);
SymmetricTensor3 expansion_adjusted_increment(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment,
    const adlite::Scalar& point_temperature,
    const adlite::Scalar& element_temperature,
    double committed_point_temperature,
    double committed_element_temperature,
    double time_step,
    MaterialFunctionContext context);
std::array<double, 6> eigenstrain_temperature_derivative(const IsotropicThermoelasticMaterial& material,
    double temperature,
    MaterialFunctionContext context);
adlite::Scalar small_strain_element_pressure(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalAdValues& state,
    double average_strain_trace,
    double time);
SmallStrainElementPressureSystem small_strain_element_pressure_system(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    double average_strain_trace,
    double time);
FinitePointResidualCache finite_point_residual_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state);
FiniteAverageTraceValues prepare_finite_point_residuals(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    std::array<FinitePointResidualCache, hex8_node_count>& point_residuals);
CartesianRotation active_rotation(const cartesian_detail::Matrix3& rotation);
double prepare_finite_point_stresses(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* committed_material,
    double time_step,
    double average_strain_trace,
    std::array<FinitePointResidualCache, hex8_node_count>& point_residuals);
void add_finite_hex8_point_residual(const Hex8QuadraturePoint& point,
    std::size_t material_node,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    const FinitePointResidualCache& point_residual,
    double reference_volume,
    double current_volume,
    double element_pressure,
    Hex8LocalResidual& residual);
Hex8LocalResidual c3d8t_finite_residual_values(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    bool include_thermal_time_term);
FiniteElementPressureSystem finite_element_pressure_system(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* committed_material,
    double time_step,
    FiniteAverageTraceSystem& average_trace_system,
    std::array<FinitePointSystemCache, hex8_node_count>& point_systems);
void add_hex8_point_residual(const Hex8QuadraturePoint& point,
    std::size_t material_node,
    const Hex8LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    double average_strain_trace,
    const adlite::Scalar& element_pressure,
    double reference_volume,
    double finite_current_volume,
    double finite_element_pressure,
    Hex8LocalAdValues& residual);
void add_hex8_point_system(const Hex8QuadraturePoint& point,
    std::size_t material_node,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian,
    double average_strain_trace,
    const std::array<std::array<double, 3>, hex8_node_count>& average_trace_displacement_derivatives,
    const SmallStrainElementPressureSystem& element_pressure,
    double reference_volume,
    const FiniteAverageTraceSystem& finite_average_trace_system,
    const FiniteElementPressureSystem& finite_element_pressure,
    const FinitePointSystemCache* finite_point_system);
adlite::Scalar hex8_nodal_volume_measure(const Hex8Geometry& geometry,
    std::size_t node,
    const Hex8LocalAdValues& state,
    StrainFormulation strain_formulation);
Hex8NodalVolumeSystem hex8_nodal_volume_system(const Hex8Geometry& geometry,
    std::size_t node,
    const Hex8LocalValues& state,
    StrainFormulation strain_formulation);
void add_hex8_nodal_body_source(const Hex8Geometry& geometry,
    const Hex8LocalAdValues& state,
    StrainFormulation strain_formulation,
    double volumetric_heat_source,
    Hex8LocalAdValues& residual);
void add_hex8_nodal_body_source_system(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    StrainFormulation strain_formulation,
    double volumetric_heat_source,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian);
void add_hex8_lumped_capacity(const Hex8Geometry& geometry,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double time_step,
    StrainFormulation strain_formulation,
    Hex8LocalAdValues& residual);
void add_hex8_lumped_capacity_system(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double time_step,
    StrainFormulation strain_formulation,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian);
void assemble_c3d8t_finite_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian);
void assemble_c3d8t_small_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian);
std::array<SymmetricTensor3Values, 8> evaluate_hex8_stress(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time);
Hex8LocalResidual compute_hex8_local(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term);
Hex8LocalResidual compute_hex8_thermoelastic(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term);
Hex8LocalResidual compute_hex8_transient(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term);
CartesianMaterialHistory compute_hex8_transient_update(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step);
std::array<SymmetricTensor3Values, 8>
compute_hex8_stress(const elements::C3d8Input& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);

adlite::Scalar average_hex8_temperature(const Hex8LocalAdValues& state) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        result += state[node] / 8.0;
    return result;
}

double average_hex8_temperature(const Hex8LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        result += state[node] / 8.0;
    return result;
}

double average_hex8_strain_trace(const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result += geometry.average_shape_gradient[node][component] * state[8 * (component + 1) + node];
    return result;
}

FiniteTracePoint finite_trace_point(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state) {
    const cartesian_detail::Matrix3 old = deformation_gradient(point, committed_state);
    const cartesian_detail::KinematicsCore kinematics =
        cartesian_detail::evaluate_kinematics(gradient, old, StrainFormulation::finite);
    cartesian_detail::ActiveMatrix3 midpoint{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            midpoint[i][j] = 0.5 * (gradient[i][j] + old[i][j]);
            if (i == j)
                midpoint[i][j] += 0.5;
        }
    const adlite::Scalar midpoint_determinant = cartesian_detail::determinant(midpoint);
    if (!std::isfinite(midpoint_determinant.value()) || !(midpoint_determinant.value() > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint configuration must preserve a positive Jacobian");
    return {kinematics.strain_increment.xx + kinematics.strain_increment.yy + kinematics.strain_increment.zz,
        point.weighted_measure * midpoint_determinant,
        point.weighted_measure * kinematics.current_determinant};
}

FiniteAverageTraceValues finite_average_hex8_strain_trace_values(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state) {
    Hex8LocalAdValues passive{};
    ad_local_system::make_passive(state.data(), state.size(), passive.data());
    double numerator = 0.0, volume = 0.0;
    FiniteAverageTraceValues result;
    for (const Hex8QuadraturePoint& point : geometry.points) {
        const FiniteTracePoint trace_point =
            finite_trace_point(point, displacement_gradient(point, passive), committed_state);
        numerator += trace_point.midpoint_weighted_measure.value() * trace_point.strain_trace.value();
        volume += trace_point.midpoint_weighted_measure.value();
        result.current_volume += trace_point.current_weighted_measure.value();
    }
    if (!std::isfinite(volume) || !(volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint volume must be finite and positive");
    if (!std::isfinite(result.current_volume) || !(result.current_volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 current volume must be finite and positive");
    result.value = numerator / volume;
    return result;
}

double finite_average_hex8_strain_trace(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state) {
    return finite_average_hex8_strain_trace_values(geometry, state, committed_state).value;
}

SymmetricTensor3
selectively_reduced_strain(const SymmetricTensor3& strain, double average_trace, StrainFormulation strain_formulation) {
    if (strain_formulation != StrainFormulation::small && strain_formulation != StrainFormulation::finite)
        return strain;
    const adlite::Scalar correction = (average_trace - strain.xx - strain.yy - strain.zz) / 3.0;
    return {strain.xx + correction, strain.yy + correction, strain.zz + correction, strain.xy, strain.yz, strain.xz};
}

SymmetricTensor3 expansion_adjusted_strain(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain,
    const adlite::Scalar& point_temperature,
    const adlite::Scalar& element_temperature,
    MaterialFunctionContext context) {
    const SymmetricTensor3 point_imposed = material.eigenstrain(point_temperature, context);
    const SymmetricTensor3 element_imposed = material.eigenstrain(element_temperature, context);
    return {strain.xx + point_imposed.xx - element_imposed.xx,
        strain.yy + point_imposed.yy - element_imposed.yy,
        strain.zz + point_imposed.zz - element_imposed.zz,
        strain.xy + point_imposed.xy - element_imposed.xy,
        strain.yz + point_imposed.yz - element_imposed.yz,
        strain.xz + point_imposed.xz - element_imposed.xz};
}

SymmetricTensor3 expansion_adjusted_increment(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment,
    const adlite::Scalar& point_temperature,
    const adlite::Scalar& element_temperature,
    double committed_point_temperature,
    double committed_element_temperature,
    double time_step,
    MaterialFunctionContext context) {
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

std::array<double, 6> eigenstrain_temperature_derivative(const IsotropicThermoelasticMaterial& material,
    double temperature,
    MaterialFunctionContext context) {
    const adlite::Scalar active_temperature = adlite::Scalar::independent(temperature, 0, 1);
    const SymmetricTensor3 imposed = material.eigenstrain(active_temperature, context);
    const std::array<const adlite::Scalar*, 6> components =
        {&imposed.xx, &imposed.yy, &imposed.zz, &imposed.xy, &imposed.yz, &imposed.xz};
    std::array<double, 6> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        components[component]->copy_derivatives(&result[component], 1);
    return result;
}

adlite::Scalar small_strain_element_pressure(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalAdValues& state,
    double average_strain_trace,
    double time) {
    adlite::Scalar average_bulk_modulus = 0.0;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = state[hex8_node_gauss_permutation[q]];
        const ActiveThermoelasticProperties properties =
            material.active_properties(temperature, material_context(time, point.position));
        average_bulk_modulus += point.weighted_measure / geometry.reference_volume
                                * (properties.lame_lambda + 2.0 * properties.shear_modulus / 3.0);
    }
    const adlite::Scalar average_temperature = average_hex8_temperature(state);
    const SymmetricTensor3 imposed =
        material.eigenstrain(average_temperature, material_context(time, geometry.selective_position));
    return average_bulk_modulus * (average_strain_trace - imposed.xx - imposed.yy - imposed.zz);
}

SmallStrainElementPressureSystem small_strain_element_pressure_system(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    double average_strain_trace,
    double time) {
    SmallStrainElementPressureSystem result;
    double average_bulk_modulus = 0.0;
    std::array<double, hex8_node_count> bulk_modulus_temperature_derivatives{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const std::size_t material_node = hex8_node_gauss_permutation[q];
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
        result.temperature_derivatives[node] = bulk_modulus_temperature_derivatives[node] * elastic_trace
                                               - average_bulk_modulus
                                                     * (imposed_trace.is_active() ? imposed_trace.derivative(0) : 0.0)
                                                     / static_cast<double>(hex8_node_count);
    return result;
}

FinitePointResidualCache finite_point_residual_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state) {
    FinitePointResidualCache result;
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                result.displacement_gradient[component][direction] +=
                    point.gradient[node][direction] * state[8 * (component + 1) + node];
    result.current_deformation = result.displacement_gradient;
    for (std::size_t direction = 0; direction < 3; ++direction)
        result.current_deformation[direction][direction] += 1.0;
    const double current_determinant = cartesian_detail::determinant(result.current_deformation);
    if (!std::isfinite(current_determinant) || !(current_determinant > 0.0))
        throw std::domain_error("Finite-strain Cartesian deformation must preserve a positive Jacobian");
    const cartesian_detail::Matrix3 current_inverse =
        cartesian_detail::inverse(result.current_deformation, current_determinant);
    const cartesian_detail::Matrix3 old_deformation = deformation_gradient(point, committed_state);
    const double old_determinant = cartesian_detail::determinant(old_deformation);
    if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
        throw std::domain_error("Committed finite-strain Cartesian state requires a positive Jacobian");
    const double incremental_determinant = current_determinant / old_determinant;
    if (!std::isfinite(incremental_determinant) || !(incremental_determinant > 0.0))
        throw std::domain_error("Incremental finite-strain Cartesian state requires a positive Jacobian");
    // Match the active Hughes-Winget arithmetic order so residual-only and Jacobian calls return the same primal
    // values while this path keeps every derivative array out of the residual callback.
    cartesian_detail::Matrix3 deformation_sum{}, deformation_difference{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            deformation_sum[i][j] = result.current_deformation[i][j] + old_deformation[i][j];
            deformation_difference[i][j] = result.current_deformation[i][j] - old_deformation[i][j];
        }
    const double plus_determinant = cartesian_detail::determinant(deformation_sum);
    if (!std::isfinite(plus_determinant) || plus_determinant == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian increment has singular delta-F plus identity");
    const cartesian_detail::Matrix3 plus_inverse = cartesian_detail::inverse(deformation_sum, plus_determinant);
    cartesian_detail::Matrix3 hughes_winget{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                hughes_winget[i][j] += 2.0 * deformation_difference[i][k] * plus_inverse[k][j];
    cartesian_detail::Matrix3 spatial_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            spatial_strain[i][j] = 0.5 * (hughes_winget[i][j] + hughes_winget[j][i]);

    cartesian_detail::Matrix3 rotation_numerator{}, rotation_denominator{};
    for (std::size_t i = 0; i < 3; ++i) {
        rotation_numerator[i][i] = 1.0;
        rotation_denominator[i][i] = 1.0;
        for (std::size_t j = 0; j < 3; ++j) {
            const double half_spin = 0.25 * (hughes_winget[i][j] - hughes_winget[j][i]);
            rotation_numerator[i][j] += half_spin;
            rotation_denominator[i][j] -= half_spin;
        }
    }
    const double rotation_denominator_determinant = cartesian_detail::determinant(rotation_denominator);
    if (!std::isfinite(rotation_denominator_determinant) || rotation_denominator_determinant == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian rotation denominator is singular");
    result.rotation = cartesian_detail::multiply(rotation_numerator,
        cartesian_detail::inverse(rotation_denominator, rotation_denominator_determinant));
    cartesian_detail::Matrix3 spatial_times_rotation{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                spatial_times_rotation[i][j] += spatial_strain[i][k] * result.rotation[k][j];
    cartesian_detail::Matrix3 corotational_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                corotational_strain[i][j] += result.rotation[k][i] * spatial_times_rotation[k][j];
    result.strain_increment = {corotational_strain[0][0],
        corotational_strain[1][1],
        corotational_strain[2][2],
        corotational_strain[0][1],
        corotational_strain[1][2],
        corotational_strain[0][2]};
    result.current_weighted_measure = point.weighted_measure * current_determinant;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.current_gradient[node][direction] +=
                    point.gradient[node][reference] * current_inverse[reference][direction];
    return result;
}

FiniteAverageTraceValues prepare_finite_point_residuals(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    std::array<FinitePointResidualCache, hex8_node_count>& point_residuals) {
    FiniteAverageTraceValues result;
    double numerator = 0.0, midpoint_volume = 0.0;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        FinitePointResidualCache& point_residual = point_residuals[q];
        point_residual = finite_point_residual_kinematics(point, state, committed_state);
        const cartesian_detail::Matrix3 old_deformation = deformation_gradient(point, committed_state);
        cartesian_detail::Matrix3 midpoint{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                midpoint[i][j] = 0.5 * (point_residual.displacement_gradient[i][j] + old_deformation[i][j]);
                if (i == j)
                    midpoint[i][j] += 0.5;
            }
        const double point_midpoint_volume = point.weighted_measure * cartesian_detail::determinant(midpoint);
        if (!std::isfinite(point_midpoint_volume) || !(point_midpoint_volume > 0.0))
            throw std::domain_error(
                "Abaqus finite-strain HEX8 midpoint configuration must preserve a positive Jacobian");
        numerator += point_midpoint_volume
                     * (point_residual.strain_increment.xx + point_residual.strain_increment.yy
                         + point_residual.strain_increment.zz);
        midpoint_volume += point_midpoint_volume;
        result.current_volume += point_residual.current_weighted_measure;
    }
    if (!std::isfinite(midpoint_volume) || !(midpoint_volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint volume must be finite and positive");
    if (!std::isfinite(result.current_volume) || !(result.current_volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 current volume must be finite and positive");
    result.value = numerator / midpoint_volume;
    return result;
}

CartesianRotation active_rotation(const cartesian_detail::Matrix3& rotation) {
    return {rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2]};
}

double prepare_finite_point_stresses(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* committed_material,
    double time_step,
    double average_strain_trace,
    std::array<FinitePointResidualCache, hex8_node_count>& point_residuals) {
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const adlite::Scalar expansion_temperature(average_hex8_temperature(state));
    const double old_expansion_temperature = average_hex8_temperature(old_state);
    double element_pressure = 0.0;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const std::size_t material_node = hex8_node_gauss_permutation[q];
        FinitePointResidualCache& point_residual = point_residuals[q];
        const adlite::Scalar temperature(state[material_node]);
        const SymmetricTensor3 strain{point_residual.strain_increment.xx,
            point_residual.strain_increment.yy,
            point_residual.strain_increment.zz,
            point_residual.strain_increment.xy,
            point_residual.strain_increment.yz,
            point_residual.strain_increment.xz};
        const SymmetricTensor3 constitutive_strain =
            selectively_reduced_strain(strain, average_strain_trace, StrainFormulation::finite);
        const MaterialFunctionContext context = material_context(time, point.position);
        SymmetricTensor3 stress;
        if (committed_material == nullptr) {
            stress = material.stress(
                expansion_adjusted_strain(material, constitutive_strain, temperature, expansion_temperature, context),
                temperature,
                context);
            stress = rotate_cartesian_tensor(stress, active_rotation(point_residual.rotation));
        } else {
            stress = evaluate_incremental_cartesian_response(material,
                expansion_adjusted_increment(material,
                    constitutive_strain,
                    temperature,
                    expansion_temperature,
                    old_state[material_node],
                    old_expansion_temperature,
                    time_step,
                    context),
                temperature,
                old_state[material_node],
                time_step,
                (*committed_material)[q],
                context)
                         .stress;
            stress = rotate_cartesian_tensor(stress, active_rotation(point_residual.rotation));
        }
        point_residual.stress = {stress.xx.value(),
            stress.yy.value(),
            stress.zz.value(),
            stress.xy.value(),
            stress.yz.value(),
            stress.xz.value()};
        element_pressure += point.weighted_measure / geometry.reference_volume
                            * (point_residual.stress.xx + point_residual.stress.yy + point_residual.stress.zz) / 3.0;
    }
    return element_pressure;
}

void add_finite_hex8_point_residual(const Hex8QuadraturePoint& point,
    std::size_t material_node,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    const FinitePointResidualCache& point_residual,
    double reference_volume,
    double current_volume,
    double element_pressure,
    Hex8LocalResidual& residual) {
    const double conductivity =
        material.conductivity(adlite::Scalar(state[material_node]), material_context(time, point.position)).value();
    std::array<double, 3> temperature_gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += point_residual.current_gradient[node][direction] * state[node];
    const double point_pressure =
        (point_residual.stress.xx + point_residual.stress.yy + point_residual.stress.zz) / 3.0;
    const double deviatoric_scale =
        point.weighted_measure * current_volume / (reference_volume * point_residual.current_weighted_measure);
    const SymmetricTensor3Values stress{deviatoric_scale * (point_residual.stress.xx - point_pressure)
                                            + element_pressure,
        deviatoric_scale * (point_residual.stress.yy - point_pressure) + element_pressure,
        deviatoric_scale * (point_residual.stress.zz - point_pressure) + element_pressure,
        deviatoric_scale * point_residual.stress.xy,
        deviatoric_scale * point_residual.stress.yz,
        deviatoric_scale * point_residual.stress.xz};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const double gradient_x = point_residual.current_gradient[node][0];
        const double gradient_y = point_residual.current_gradient[node][1];
        const double gradient_z = point_residual.current_gradient[node][2];
        residual[node] += point_residual.current_weighted_measure * conductivity
                          * (gradient_x * temperature_gradient[0] + gradient_y * temperature_gradient[1]
                              + gradient_z * temperature_gradient[2]);
        residual[8 + node] += point_residual.current_weighted_measure
                              * (stress.xx * gradient_x + stress.xy * gradient_y + stress.xz * gradient_z);
        residual[16 + node] += point_residual.current_weighted_measure
                               * (stress.xy * gradient_x + stress.yy * gradient_y + stress.yz * gradient_z);
        residual[24 + node] += point_residual.current_weighted_measure
                               * (stress.xz * gradient_x + stress.yz * gradient_y + stress.zz * gradient_z);
    }
}

Hex8LocalResidual c3d8t_finite_residual_values(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    bool include_thermal_time_term) {
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    std::array<FinitePointResidualCache, hex8_node_count> point_residuals{};
    const FiniteAverageTraceValues average_trace =
        prepare_finite_point_residuals(geometry, state, old_state, point_residuals);
    const double element_pressure = prepare_finite_point_stresses(data.material,
        geometry,
        state,
        data.time,
        committed_state,
        history,
        time_step,
        average_trace.value,
        point_residuals);

    Hex8LocalResidual residual{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q)
        add_finite_hex8_point_residual(geometry.points[q],
            hex8_node_gauss_permutation[q],
            state,
            data.material,
            data.time,
            point_residuals[q],
            geometry.reference_volume,
            average_trace.current_volume,
            element_pressure,
            residual);
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const double nodal_measure = point_residuals[hex8_node_gauss_permutation[node]].current_weighted_measure;
        residual[node] -= nodal_measure * data.volumetric_heat_source;
        if (committed_state != nullptr && include_thermal_time_term) {
            const double capacity = data.material
                                        .heat_capacity(adlite::Scalar(state[node]),
                                            material_context(data.time, geometry.capacity_points[node].position))
                                        .value();
            residual[node] += nodal_measure * capacity * (state[node] - (*committed_state)[node]) / time_step;
        }
    }
    return residual;
}

FiniteElementPressureSystem finite_element_pressure_system(const IsotropicThermoelasticMaterial& material,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* committed_material,
    double time_step,
    FiniteAverageTraceSystem& average_trace_system,
    std::array<FinitePointSystemCache, hex8_node_count>& point_systems) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const double expansion_temperature = average_hex8_temperature(state);
    FiniteElementPressureSystem result;
    double trace_numerator = 0.0, midpoint_volume = 0.0, average_trace_pressure_derivative = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> trace_numerator_derivatives{}, midpoint_volume_derivatives{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const std::size_t material_node = hex8_node_gauss_permutation[q];
        FinitePointSystemCache& point_system = point_systems[q];
        cartesian_detail::ActiveMatrix3 gradient{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                double value = 0.0;
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    value += point.gradient[node][direction] * state[8 * (component + 1) + node];
                gradient[component][direction] =
                    adlite::Scalar::independent(value, component * 3 + direction, point_width);
            }
        point_system.active_temperature =
            adlite::Scalar::independent(state[material_node], temperature_index, point_width);
        point_system.kinematics =
            evaluate_cartesian_kinematics_from_gradient(point, gradient, old_state, StrainFormulation::finite);
        const C3d8Kinematics& kinematics = point_system.kinematics;
        const cartesian_detail::Matrix3 old_gradient = deformation_gradient(point, old_state);
        cartesian_detail::ActiveMatrix3 midpoint{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                midpoint[i][j] = 0.5 * (gradient[i][j] + old_gradient[i][j]);
                if (i == j)
                    midpoint[i][j] += 0.5;
            }
        const adlite::Scalar midpoint_determinant = cartesian_detail::determinant(midpoint);
        if (!std::isfinite(midpoint_determinant.value()) || !(midpoint_determinant.value() > 0.0))
            throw std::domain_error(
                "Abaqus finite-strain HEX8 midpoint configuration must preserve a positive Jacobian");
        const adlite::Scalar point_midpoint_volume = point.weighted_measure * midpoint_determinant;
        const adlite::Scalar point_trace =
            kinematics.strain_increment.xx + kinematics.strain_increment.yy + kinematics.strain_increment.zz;
        const adlite::Scalar point_trace_numerator = point_midpoint_volume * point_trace;
        trace_numerator += point_trace_numerator.value();
        midpoint_volume += point_midpoint_volume.value();
        std::array<double, point_width> point_trace_numerator_derivatives{}, point_midpoint_volume_derivatives{},
            point_current_volume_derivatives{};
        point_trace_numerator.copy_derivatives(point_trace_numerator_derivatives.data(),
            point_trace_numerator_derivatives.size());
        point_midpoint_volume.copy_derivatives(point_midpoint_volume_derivatives.data(),
            point_midpoint_volume_derivatives.size());
        kinematics.current_weighted_measure.copy_derivatives(point_current_volume_derivatives.data(),
            point_current_volume_derivatives.size());
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t direction = 0; direction < 3; ++direction) {
                    const double shape_chain = point.gradient[node][direction];
                    trace_numerator_derivatives[node][component] +=
                        point_trace_numerator_derivatives[component * 3 + direction] * shape_chain;
                    midpoint_volume_derivatives[node][component] +=
                        point_midpoint_volume_derivatives[component * 3 + direction] * shape_chain;
                    average_trace_system.current_volume_derivatives[node][component] +=
                        point_current_volume_derivatives[component * 3 + direction] * shape_chain;
                }
        const SymmetricTensor3 constitutive_strain = selectively_reduced_strain(kinematics.strain_increment,
            average_trace_system.value,
            StrainFormulation::finite);
        const std::array<const adlite::Scalar*, 6> strain_components = {&constitutive_strain.xx,
            &constitutive_strain.yy,
            &constitutive_strain.zz,
            &constitutive_strain.xy,
            &constitutive_strain.yz,
            &constitutive_strain.xz};
        const MaterialFunctionContext context = material_context(time, point.position);
        const SymmetricTensor3 point_imposed = material.eigenstrain(adlite::Scalar(state[material_node]), context);
        const SymmetricTensor3 element_imposed = material.eigenstrain(adlite::Scalar(expansion_temperature), context);
        const std::array<double, 6> point_imposed_values = {point_imposed.xx.value(),
            point_imposed.yy.value(),
            point_imposed.zz.value(),
            point_imposed.xy.value(),
            point_imposed.yz.value(),
            point_imposed.xz.value()};
        const std::array<double, 6> element_imposed_values = {element_imposed.xx.value(),
            element_imposed.yy.value(),
            element_imposed.zz.value(),
            element_imposed.xy.value(),
            element_imposed.yz.value(),
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
                old_element_imposed.yy.value(),
                old_element_imposed.zz.value(),
                old_element_imposed.xy.value(),
                old_element_imposed.yz.value(),
                old_element_imposed.xz.value()};
            for (std::size_t component = 0; component < 6; ++component)
                fed_strain[component] = point_committed_material->elastic_strain[component]
                                        + strain_components[component]->value() + old_element_values[component]
                                        + point_committed_material->plastic_strain[component]
                                        + point_committed_material->creep_strain[component]
                                        + point_imposed_values[component] - element_imposed_values[component];
        } else {
            for (std::size_t component = 0; component < 6; ++component)
                fed_strain[component] = strain_components[component]->value() + point_imposed_values[component]
                                        - element_imposed_values[component];
        }
        point_system.tangent = cartesian_detail::evaluate_stress_tangent(material,
            fed_strain,
            state[material_node],
            time_step,
            point_committed_material,
            context);
        const cartesian_detail::CartesianStressTangent& tangent = point_system.tangent;
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
        compose_inputs[6] = point_system.active_temperature;
        std::array<double, 7> partials{};
        for (std::size_t component = 0; component < 6; ++component)
            partials[component] = pressure_strain_derivatives[component];
        partials[6] = point_temperature_derivative;
        const double pressure_value = (tangent.stress.xx + tangent.stress.yy + tangent.stress.zz) / 3.0;
        const adlite::Scalar pressure =
            adlite::compose(pressure_value, compose_inputs.data(), partials.data(), compose_inputs.size());
        const double normalized_weight = point.weighted_measure / geometry.reference_volume;
        result.value += normalized_weight * pressure.value();
        average_trace_pressure_derivative += normalized_weight * average_trace_derivative;
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
                result.derivatives[8 * (displacement_component + 1) + node] += normalized_weight * local_chain;
            }
        }
    }
    if (!std::isfinite(midpoint_volume) || !(midpoint_volume > 0.0))
        throw std::domain_error("Abaqus finite-strain HEX8 midpoint volume must be finite and positive");
    const double active_average_trace = trace_numerator / midpoint_volume;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component) {
            average_trace_system.displacement_derivatives[node][component] =
                (trace_numerator_derivatives[node][component]
                    - active_average_trace * midpoint_volume_derivatives[node][component])
                / midpoint_volume;
            result.derivatives[8 * (component + 1) + node] +=
                average_trace_pressure_derivative * average_trace_system.displacement_derivatives[node][component];
        }
    return result;
}

void add_hex8_point_residual(const Hex8QuadraturePoint& point,
    std::size_t material_node,
    const Hex8LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    double average_strain_trace,
    const adlite::Scalar& element_pressure,
    double reference_volume,
    double finite_current_volume,
    double finite_element_pressure,
    Hex8LocalAdValues& residual) {
    const adlite::Scalar temperature = state[material_node];
    const adlite::Scalar expansion_temperature = average_hex8_temperature(state);
    const MaterialFunctionContext context = material_context(time, point.position);
    const adlite::Scalar conductivity = material.conductivity(temperature, context);
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const C3d8Kinematics kinematics =
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
            temperature,
            context);
        if (strain_formulation == StrainFormulation::finite)
            stress = rotate_cartesian_tensor(stress, kinematics.rotation);
    } else if (strain_formulation == StrainFormulation::finite) {
        const double old_temperature = old_state[material_node];
        const double old_expansion_temperature = average_hex8_temperature(old_state);
        stress = material
                     .incremental_response(expansion_adjusted_increment(material,
                                               constitutive_strain,
                                               temperature,
                                               expansion_temperature,
                                               old_temperature,
                                               old_expansion_temperature,
                                               time_step,
                                               context),
                         kinematics.rotation,
                         temperature,
                         old_temperature,
                         time_step,
                         *committed_material,
                         context)
                     .stress;
    } else {
        stress = material
                     .response(expansion_adjusted_strain(material,
                                   constitutive_strain,
                                   temperature,
                                   expansion_temperature,
                                   context),
                         temperature,
                         time_step,
                         *committed_material,
                         context)
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
        residual[node] += kinematics.current_weighted_measure * conductivity
                          * (current_gradient_x * gradient_temperature_x + current_gradient_y * gradient_temperature_y
                              + current_gradient_z * gradient_temperature_z);
        residual[8 + node] +=
            kinematics.current_weighted_measure
            * (stress.xx * current_gradient_x + stress.xy * current_gradient_y + stress.xz * current_gradient_z);
        residual[16 + node] +=
            kinematics.current_weighted_measure
            * (stress.xy * current_gradient_x + stress.yy * current_gradient_y + stress.yz * current_gradient_z);
        residual[24 + node] +=
            kinematics.current_weighted_measure
            * (stress.xz * current_gradient_x + stress.yz * current_gradient_y + stress.zz * current_gradient_z);
    }
}

void add_hex8_point_system(const Hex8QuadraturePoint& point,
    std::size_t material_node,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian,
    double average_strain_trace,
    const std::array<std::array<double, 3>, hex8_node_count>& average_trace_displacement_derivatives,
    const SmallStrainElementPressureSystem& element_pressure,
    double reference_volume,
    const FiniteAverageTraceSystem& finite_average_trace_system,
    const FiniteElementPressureSystem& finite_element_pressure,
    const FinitePointSystemCache* finite_point_system) {
    constexpr std::size_t point_width = 10, temperature_index = 9;
    const double temperature_value = state[material_node];
    const double expansion_temperature_value = average_hex8_temperature(state);
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    adlite::Scalar local_active_temperature;
    C3d8Kinematics local_kinematics;
    if (finite_point_system == nullptr) {
        cartesian_detail::ActiveMatrix3 active_gradient{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                double value = 0.0;
                for (std::size_t node = 0; node < 8; ++node)
                    value += point.gradient[node][direction] * state[8 * (component + 1) + node];
                active_gradient[component][direction] =
                    adlite::Scalar::independent(value, component * 3 + direction, point_width);
            }
        local_active_temperature = adlite::Scalar::independent(temperature_value, temperature_index, point_width);
        local_kinematics =
            evaluate_cartesian_kinematics_from_gradient(point, active_gradient, old_state, strain_formulation);
    }
    const adlite::Scalar& active_temperature =
        finite_point_system == nullptr ? local_active_temperature : finite_point_system->active_temperature;
    const C3d8Kinematics& kinematics =
        finite_point_system == nullptr ? local_kinematics : finite_point_system->kinematics;
    const SymmetricTensor3 constitutive_strain =
        selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, strain_formulation);
    const MaterialFunctionContext context = material_context(time, point.position);
    const std::array<const adlite::Scalar*, 6> strain_components = {&constitutive_strain.xx,
        &constitutive_strain.yy,
        &constitutive_strain.zz,
        &constitutive_strain.xy,
        &constitutive_strain.yz,
        &constitutive_strain.xz};
    cartesian_detail::CartesianStressTangent local_tangent;
    if (finite_point_system == nullptr) {
        std::array<double, 6> fed_strain{};
        const SymmetricTensor3 point_imposed = material.eigenstrain(adlite::Scalar(temperature_value), context);
        const SymmetricTensor3 element_imposed =
            material.eigenstrain(adlite::Scalar(expansion_temperature_value), context);
        const std::array<double, 6> point_imposed_values = {point_imposed.xx.value(),
            point_imposed.yy.value(),
            point_imposed.zz.value(),
            point_imposed.xy.value(),
            point_imposed.yz.value(),
            point_imposed.xz.value()};
        const std::array<double, 6> element_imposed_values = {element_imposed.xx.value(),
            element_imposed.yy.value(),
            element_imposed.zz.value(),
            element_imposed.xy.value(),
            element_imposed.yz.value(),
            element_imposed.xz.value()};
        if (committed_material != nullptr && strain_formulation == StrainFormulation::finite) {
            const double old_temperature = old_state[material_node];
            if (!std::isfinite(old_temperature) || !(old_temperature > 0.0))
                throw std::domain_error(
                    "Incremental Cartesian material committed temperature must be finite and positive");
            MaterialFunctionContext old_context = context;
            old_context.time -= time_step;
            const SymmetricTensor3 old_element_imposed =
                material.eigenstrain(adlite::Scalar(average_hex8_temperature(old_state)), old_context);
            const std::array<double, 6> old_element_values = {old_element_imposed.xx.value(),
                old_element_imposed.yy.value(),
                old_element_imposed.zz.value(),
                old_element_imposed.xy.value(),
                old_element_imposed.yz.value(),
                old_element_imposed.xz.value()};
            for (std::size_t component = 0; component < 6; ++component)
                fed_strain[component] = committed_material->elastic_strain[component]
                                        + strain_components[component]->value() + old_element_values[component]
                                        + committed_material->plastic_strain[component]
                                        + committed_material->creep_strain[component] + point_imposed_values[component]
                                        - element_imposed_values[component];
        } else {
            for (std::size_t component = 0; component < 6; ++component)
                fed_strain[component] = strain_components[component]->value() + point_imposed_values[component]
                                        - element_imposed_values[component];
        }
        local_tangent = cartesian_detail::evaluate_stress_tangent(material,
            fed_strain,
            temperature_value,
            time_step,
            committed_material,
            context);
    }
    const cartesian_detail::CartesianStressTangent& tangent =
        finite_point_system == nullptr ? local_tangent : finite_point_system->tangent;
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
        partials[6] = point_thermal[component];
        composed[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    SymmetricTensor3 stress{composed[0], composed[1], composed[2], composed[3], composed[4], composed[5]};
    SymmetricTensor3Values expansion_stress_derivative{expansion_thermal[0],
        expansion_thermal[1],
        expansion_thermal[2],
        expansion_thermal[3],
        expansion_thermal[4],
        expansion_thermal[5]};
    if (strain_formulation == StrainFormulation::small) {
        const adlite::Scalar point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        stress.xx += element_pressure.value - point_pressure;
        stress.yy += element_pressure.value - point_pressure;
        stress.zz += element_pressure.value - point_pressure;
        average_trace_stress_derivative = {element_pressure.trace_derivative,
            element_pressure.trace_derivative,
            element_pressure.trace_derivative,
            0.0,
            0.0,
            0.0};
        expansion_stress_derivative = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    }
    SymmetricTensor3 finite_point_deviatoric{};
    if (strain_formulation == StrainFormulation::finite) {
        stress = rotate_cartesian_tensor(stress, kinematics.rotation);
        const cartesian_detail::Matrix3 passive_rotation = {
            {{kinematics.rotation.xx.value(), kinematics.rotation.xy.value(), kinematics.rotation.xz.value()},
                {kinematics.rotation.yx.value(), kinematics.rotation.yy.value(), kinematics.rotation.yz.value()},
                {kinematics.rotation.zx.value(), kinematics.rotation.zy.value(), kinematics.rotation.zz.value()}}};
        expansion_stress_derivative = rotate_cartesian_tensor_values(expansion_stress_derivative, passive_rotation);
        const SymmetricTensor3Values rotated_average_trace_derivative =
            rotate_cartesian_tensor_values({average_trace_stress_derivative[0],
                                               average_trace_stress_derivative[1],
                                               average_trace_stress_derivative[2],
                                               average_trace_stress_derivative[3],
                                               average_trace_stress_derivative[4],
                                               average_trace_stress_derivative[5]},
                passive_rotation);
        average_trace_stress_derivative = {rotated_average_trace_derivative.xx,
            rotated_average_trace_derivative.yy,
            rotated_average_trace_derivative.zz,
            rotated_average_trace_derivative.xy,
            rotated_average_trace_derivative.yz,
            rotated_average_trace_derivative.xz};
        const adlite::Scalar point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        finite_point_deviatoric = {stress.xx - point_pressure,
            stress.yy - point_pressure,
            stress.zz - point_pressure,
            stress.xy,
            stress.yz,
            stress.xz};
        const adlite::Scalar deviatoric_scale = point.weighted_measure * finite_average_trace_system.current_volume
                                                / (reference_volume * kinematics.current_weighted_measure);
        stress = {deviatoric_scale * finite_point_deviatoric.xx + finite_element_pressure.value,
            deviatoric_scale * finite_point_deviatoric.yy + finite_element_pressure.value,
            deviatoric_scale * finite_point_deviatoric.zz + finite_element_pressure.value,
            deviatoric_scale * finite_point_deviatoric.xy,
            deviatoric_scale * finite_point_deviatoric.yz,
            deviatoric_scale * finite_point_deviatoric.xz};
        const double expansion_pressure =
            (expansion_stress_derivative.xx + expansion_stress_derivative.yy + expansion_stress_derivative.zz) / 3.0;
        expansion_stress_derivative = {deviatoric_scale.value() * (expansion_stress_derivative.xx - expansion_pressure),
            deviatoric_scale.value() * (expansion_stress_derivative.yy - expansion_pressure),
            deviatoric_scale.value() * (expansion_stress_derivative.zz - expansion_pressure),
            deviatoric_scale.value() * expansion_stress_derivative.xy,
            deviatoric_scale.value() * expansion_stress_derivative.yz,
            deviatoric_scale.value() * expansion_stress_derivative.xz};
        const double average_trace_pressure = (average_trace_stress_derivative[0] + average_trace_stress_derivative[1]
                                                  + average_trace_stress_derivative[2])
                                              / 3.0;
        average_trace_stress_derivative = {deviatoric_scale.value()
                                               * (average_trace_stress_derivative[0] - average_trace_pressure),
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
            kinematics.current_weighted_measure * conductivity
            * (current_gradient_x * gradient_temperature_x + current_gradient_y * gradient_temperature_y
                + current_gradient_z * gradient_temperature_z);
        point_residual[8 + node] +=
            kinematics.current_weighted_measure
            * (stress.xx * current_gradient_x + stress.xy * current_gradient_y + stress.xz * current_gradient_z);
        point_residual[16 + node] +=
            kinematics.current_weighted_measure
            * (stress.xy * current_gradient_x + stress.yy * current_gradient_y + stress.yz * current_gradient_z);
        point_residual[24 + node] +=
            kinematics.current_weighted_measure
            * (stress.xz * current_gradient_x + stress.yz * current_gradient_y + stress.zz * current_gradient_z);
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
                gradient_dot += kinematics.current_gradient[node][direction].value()
                                * kinematics.current_gradient[other][direction].value();
            jacobian[node * 32 + other] +=
                kinematics.current_weighted_measure.value() * conductivity_value * gradient_dot
                + (other == material_node ? derivatives[temperature_index] : 0.0);
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
                component == 0   ? kinematics.current_weighted_measure.value()
                                       * (expansion_stress_derivative.xx * current_gradient_x.value()
                                           + expansion_stress_derivative.xy * current_gradient_y.value()
                                           + expansion_stress_derivative.xz * current_gradient_z.value())
                : component == 1 ? kinematics.current_weighted_measure.value()
                                       * (expansion_stress_derivative.xy * current_gradient_x.value()
                                           + expansion_stress_derivative.yy * current_gradient_y.value()
                                           + expansion_stress_derivative.yz * current_gradient_z.value())
                                 : kinematics.current_weighted_measure.value()
                                       * (expansion_stress_derivative.xz * current_gradient_x.value()
                                           + expansion_stress_derivative.yz * current_gradient_y.value()
                                           + expansion_stress_derivative.zz * current_gradient_z.value());
            const double average_trace_derivative =
                component == 0   ? kinematics.current_weighted_measure.value()
                                       * (average_trace_stress_derivative[0] * current_gradient_x.value()
                                           + average_trace_stress_derivative[3] * current_gradient_y.value()
                                           + average_trace_stress_derivative[5] * current_gradient_z.value())
                : component == 1 ? kinematics.current_weighted_measure.value()
                                       * (average_trace_stress_derivative[3] * current_gradient_x.value()
                                           + average_trace_stress_derivative[1] * current_gradient_y.value()
                                           + average_trace_stress_derivative[4] * current_gradient_z.value())
                                 : kinematics.current_weighted_measure.value()
                                       * (average_trace_stress_derivative[5] * current_gradient_x.value()
                                           + average_trace_stress_derivative[4] * current_gradient_y.value()
                                           + average_trace_stress_derivative[2] * current_gradient_z.value());
            const double pressure_temperature_gradient = component == 0   ? current_gradient_x.value()
                                                         : component == 1 ? current_gradient_y.value()
                                                                          : current_gradient_z.value();
            const double finite_current_volume_derivative =
                strain_formulation == StrainFormulation::finite
                    ? point.weighted_measure / reference_volume
                          * (component == 0    ? finite_point_deviatoric.xx.value() * current_gradient_x.value()
                                                     + finite_point_deviatoric.xy.value() * current_gradient_y.value()
                                                     + finite_point_deviatoric.xz.value() * current_gradient_z.value()
                              : component == 1 ? finite_point_deviatoric.xy.value() * current_gradient_x.value()
                                                     + finite_point_deviatoric.yy.value() * current_gradient_y.value()
                                                     + finite_point_deviatoric.yz.value() * current_gradient_z.value()
                                               : finite_point_deviatoric.xz.value() * current_gradient_x.value()
                                                     + finite_point_deviatoric.yz.value() * current_gradient_y.value()
                                                     + finite_point_deviatoric.zz.value() * current_gradient_z.value())
                    : 0.0;
            for (std::size_t other = 0; other < 8; ++other) {
                if (other == material_node)
                    jacobian[row * 32 + other] += derivatives[temperature_index];
                jacobian[row * 32 + other] += expansion_derivative / 8.0;
                jacobian[row * 32 + other] += kinematics.current_weighted_measure.value()
                                              * element_pressure.temperature_derivatives[other]
                                              * pressure_temperature_gradient;
                if (strain_formulation == StrainFormulation::finite)
                    jacobian[row * 32 + other] += kinematics.current_weighted_measure.value()
                                                  * pressure_temperature_gradient
                                                  * finite_element_pressure.derivatives[other];
                for (std::size_t displacement_component = 0; displacement_component < 3; ++displacement_component) {
                    double chained = 0.0;
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        chained +=
                            derivatives[displacement_component * 3 + direction] * point.gradient[other][direction];
                    jacobian[row * 32 + 8 * (displacement_component + 1) + other] += chained;
                    jacobian[row * 32 + 8 * (displacement_component + 1) + other] +=
                        average_trace_derivative
                        * average_trace_displacement_derivatives[other][displacement_component];
                    if (strain_formulation == StrainFormulation::finite) {
                        const std::size_t column = 8 * (displacement_component + 1) + other;
                        jacobian[row * 32 + column] +=
                            finite_current_volume_derivative
                            * finite_average_trace_system.current_volume_derivatives[other][displacement_component];
                        jacobian[row * 32 + column] += kinematics.current_weighted_measure.value()
                                                       * pressure_temperature_gradient
                                                       * finite_element_pressure.derivatives[column];
                    }
                }
            }
        }
    }
    for (std::size_t row = 0; row < residual.size(); ++row)
        residual[row] += point_residual[row];
}

adlite::Scalar hex8_nodal_volume_measure(const Hex8Geometry& geometry,
    std::size_t node,
    const Hex8LocalAdValues& state,
    StrainFormulation strain_formulation) {
    const Hex8QuadraturePoint& point = geometry.points[hex8_node_gauss_permutation[node]];
    if (strain_formulation == StrainFormulation::small)
        return point.weighted_measure;
    cartesian_detail::ActiveMatrix3 current = displacement_gradient(point, state);
    for (std::size_t direction = 0; direction < 3; ++direction)
        current[direction][direction] += 1.0;
    const adlite::Scalar determinant = cartesian_detail::determinant(current);
    if (!std::isfinite(determinant.value()) || !(determinant.value() > 0.0))
        throw std::domain_error("Finite-strain Cartesian deformation must preserve a positive Jacobian");
    return point.weighted_measure * determinant;
}

Hex8NodalVolumeSystem hex8_nodal_volume_system(const Hex8Geometry& geometry,
    std::size_t node,
    const Hex8LocalValues& state,
    StrainFormulation strain_formulation) {
    const Hex8QuadraturePoint& point = geometry.points[hex8_node_gauss_permutation[node]];
    Hex8NodalVolumeSystem result{point.weighted_measure, {}};
    if (strain_formulation == StrainFormulation::small)
        return result;
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
    for (std::size_t direction = 0; direction < 3; ++direction)
        active_gradient[direction][direction] += 1.0;
    const adlite::Scalar determinant = cartesian_detail::determinant(active_gradient);
    if (!std::isfinite(determinant.value()) || !(determinant.value() > 0.0))
        throw std::domain_error("Finite-strain Cartesian deformation must preserve a positive Jacobian");
    const adlite::Scalar measure = point.weighted_measure * determinant;
    result.value = measure.value();
    measure.copy_derivatives(result.gradient_derivatives.data(), result.gradient_derivatives.size());
    return result;
}

void add_hex8_nodal_body_source(const Hex8Geometry& geometry,
    const Hex8LocalAdValues& state,
    StrainFormulation strain_formulation,
    double volumetric_heat_source,
    Hex8LocalAdValues& residual) {
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        residual[node] -= hex8_nodal_volume_measure(geometry, node, state, strain_formulation) * volumetric_heat_source;
}

void add_hex8_nodal_body_source_system(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    StrainFormulation strain_formulation,
    double volumetric_heat_source,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian) {
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const Hex8QuadraturePoint& point = geometry.points[hex8_node_gauss_permutation[node]];
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

void add_hex8_lumped_capacity(const Hex8Geometry& geometry,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double time_step,
    StrainFormulation strain_formulation,
    Hex8LocalAdValues& residual) {
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const Hex8CapacityPoint& point = geometry.capacity_points[node];
        const adlite::Scalar capacity = material.heat_capacity(state[node], material_context(time, point.position));
        residual[node] += hex8_nodal_volume_measure(geometry, node, state, strain_formulation) * capacity
                          * (state[node] - committed_state[node]) / time_step;
    }
}

void add_hex8_lumped_capacity_system(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const IsotropicThermoelasticMaterial& material,
    double time,
    double time_step,
    StrainFormulation strain_formulation,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian& jacobian) {
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        const Hex8CapacityPoint& point = geometry.capacity_points[node];
        const Hex8QuadraturePoint& geometry_point = geometry.points[hex8_node_gauss_permutation[node]];
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
                    derivative += measure.gradient_derivatives[component * 3 + direction]
                                  * geometry_point.gradient[other][direction];
                jacobian[node * hex8_local_dof_count + 8 * (component + 1) + other] += rate.value() * derivative;
            }
    }
}

void assemble_c3d8t_finite_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian) {
    if (jacobian == nullptr) {
        const Hex8LocalResidual passive_residual = c3d8t_finite_residual_values(data,
            geometry,
            state,
            committed_state,
            history,
            time_step,
            include_thermal_time_term);
        for (std::size_t row = 0; row < hex8_local_dof_count; ++row)
            residual[row] = passive_residual[row];
        return;
    }

    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    FiniteAverageTraceSystem average_trace;
    const FiniteAverageTraceValues values = finite_average_hex8_strain_trace_values(geometry, state, old_state);
    average_trace.value = values.value;
    average_trace.current_volume = values.current_volume;

    FiniteElementPressureSystem element_pressure;
    std::array<FinitePointSystemCache, hex8_node_count> point_systems{};
    element_pressure = finite_element_pressure_system(data.material,
        geometry,
        state,
        data.time,
        committed_state,
        history,
        time_step,
        average_trace,
        point_systems);

    jacobian->fill(0.0);
    const SmallStrainElementPressureSystem small_pressure{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q)
        add_hex8_point_system(geometry.points[q],
            hex8_node_gauss_permutation[q],
            state,
            data.material,
            StrainFormulation::finite,
            data.time,
            committed_state,
            history == nullptr ? nullptr : &(*history)[q],
            time_step,
            residual,
            *jacobian,
            average_trace.value,
            average_trace.displacement_derivatives,
            small_pressure,
            geometry.reference_volume,
            average_trace,
            element_pressure,
            &point_systems[q]);
    add_hex8_nodal_body_source_system(geometry,
        state,
        StrainFormulation::finite,
        data.volumetric_heat_source,
        residual,
        *jacobian);
    if (committed_state != nullptr && include_thermal_time_term)
        add_hex8_lumped_capacity_system(geometry,
            state,
            *committed_state,
            data.material,
            data.time,
            time_step,
            StrainFormulation::finite,
            residual,
            *jacobian);
}

void assemble_c3d8t_small_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian) {
    const double average_trace = average_hex8_strain_trace(geometry, state);
    if (jacobian == nullptr) {
        Hex8LocalAdValues passive{};
        ad_local_system::make_passive(state.data(), state.size(), passive.data());
        const adlite::Scalar element_pressure =
            small_strain_element_pressure(data.material, geometry, passive, average_trace, data.time);
        for (std::size_t q = 0; q < geometry.points.size(); ++q)
            add_hex8_point_residual(geometry.points[q],
                hex8_node_gauss_permutation[q],
                passive,
                data.material,
                StrainFormulation::small,
                data.time,
                committed_state,
                history == nullptr ? nullptr : &(*history)[q],
                time_step,
                average_trace,
                element_pressure,
                geometry.reference_volume,
                0.0,
                0.0,
                residual);
        add_hex8_nodal_body_source(geometry, passive, StrainFormulation::small, data.volumetric_heat_source, residual);
        if (committed_state != nullptr && include_thermal_time_term)
            add_hex8_lumped_capacity(geometry,
                passive,
                *committed_state,
                data.material,
                data.time,
                time_step,
                StrainFormulation::small,
                residual);
        return;
    }

    jacobian->fill(0.0);
    const SmallStrainElementPressureSystem element_pressure =
        small_strain_element_pressure_system(data.material, geometry, state, average_trace, data.time);
    const FiniteAverageTraceSystem finite_average_trace{};
    const FiniteElementPressureSystem finite_element_pressure{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q)
        add_hex8_point_system(geometry.points[q],
            hex8_node_gauss_permutation[q],
            state,
            data.material,
            StrainFormulation::small,
            data.time,
            committed_state,
            history == nullptr ? nullptr : &(*history)[q],
            time_step,
            residual,
            *jacobian,
            average_trace,
            geometry.average_shape_gradient,
            element_pressure,
            geometry.reference_volume,
            finite_average_trace,
            finite_element_pressure,
            nullptr);
    add_hex8_nodal_body_source_system(geometry,
        state,
        StrainFormulation::small,
        data.volumetric_heat_source,
        residual,
        *jacobian);
    if (committed_state != nullptr && include_thermal_time_term)
        add_hex8_lumped_capacity_system(geometry,
            state,
            *committed_state,
            data.material,
            data.time,
            time_step,
            StrainFormulation::small,
            residual,
            *jacobian);
}

std::array<SymmetricTensor3Values, 8> evaluate_hex8_stress(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time) {
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
        const adlite::Scalar temperature = ad_state[hex8_node_gauss_permutation[q]];
        const C3d8Kinematics kinematics =
            evaluate_cartesian_incremental_kinematics(point, ad_state, Hex8LocalValues{}, strain_formulation);
        const SymmetricTensor3 constitutive_strain =
            selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, strain_formulation);
        const MaterialFunctionContext context = material_context(time, point.position);
        SymmetricTensor3 stress = material.stress(
            expansion_adjusted_strain(material, constitutive_strain, temperature, expansion_temperature, context),
            temperature,
            context);
        if (strain_formulation == StrainFormulation::finite)
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

Hex8LocalResidual compute_hex8_local(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    if (committed_state != nullptr && (!std::isfinite(time_step) || time_step <= 0.0))
        throw std::invalid_argument("HEX8 time step must be finite and positive");
    if (history != nullptr && history->size() != 8)
        throw std::invalid_argument("C3D8T material history has the wrong integration point count");
    Hex8LocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    if (data.strain_formulation == StrainFormulation::finite)
        assemble_c3d8t_finite_strain_system(data,
            geometry,
            state,
            committed_state,
            history,
            time_step,
            include_thermal_time_term,
            residual,
            jacobian);
    else
        assemble_c3d8t_small_strain_system(data,
            geometry,
            state,
            committed_state,
            history,
            time_step,
            include_thermal_time_term,
            residual,
            jacobian);
    Hex8LocalResidual result{};
    ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    return result;
}

Hex8LocalResidual compute_hex8_thermoelastic(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    return compute_hex8_local(data,
        geometry,
        state,
        committed_state,
        nullptr,
        time_step,
        jacobian,
        include_thermal_time_term);
}

Hex8LocalResidual compute_hex8_transient(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    return compute_hex8_local(data,
        geometry,
        state,
        &committed_state,
        &committed_material,
        time_step,
        jacobian,
        include_thermal_time_term);
}

CartesianMaterialHistory compute_hex8_transient_update(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("HEX8 transient update time step must be finite and positive");

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
        const std::size_t material_node = hex8_node_gauss_permutation[q];
        const adlite::Scalar temperature = passive[material_node];
        const C3d8Kinematics kinematics =
            evaluate_cartesian_incremental_kinematics(point, passive, committed_state, data.strain_formulation);
        CartesianInelasticStressResponse response;
        if (data.strain_formulation == StrainFormulation::finite) {
            const double old_temperature = committed_state[material_node];
            const MaterialFunctionContext context = material_context(data.time, point.position);
            response = data.material.incremental_response(expansion_adjusted_increment(data.material,
                                                              selectively_reduced_strain(kinematics.strain_increment,
                                                                  average_strain_trace,
                                                                  data.strain_formulation),
                                                              temperature,
                                                              expansion_temperature,
                                                              old_temperature,
                                                              old_expansion_temperature,
                                                              time_step,
                                                              context),
                kinematics.rotation,
                temperature,
                old_temperature,
                time_step,
                committed_material[q],
                context);
        } else {
            const MaterialFunctionContext context = material_context(data.time, point.position);
            const SymmetricTensor3 constitutive_strain =
                selectively_reduced_strain(kinematics.strain_increment, average_strain_trace, data.strain_formulation);
            response = data.material.response(expansion_adjusted_strain(data.material,
                                                  constitutive_strain,
                                                  temperature,
                                                  expansion_temperature,
                                                  context),
                temperature,
                time_step,
                committed_material[q],
                context);
        }
        result[q] = response.trial_state;
    }
    return result;
}

std::array<SymmetricTensor3Values, 8>
compute_hex8_stress(const elements::C3d8Input& data, const Hex8Geometry& geometry, const Hex8LocalValues& state) {

    return evaluate_hex8_stress(geometry, state, data.material, data.strain_formulation, data.time);
}
} // namespace
} // namespace fuelsim

namespace fuelsim::elements {
C3d8Result evaluate_c3d8t(const C3d8Input& input, ElementRequest request) {
    const auto& data = input;
    C3d8Result result;
    if (request.residual || request.jacobian) {
        auto* tangent = request.jacobian ? &result.jacobian : nullptr;
        if (input.committed_history)
            result.residual = compute_hex8_transient(data,
                input.geometry,
                input.state,
                input.committed_state,
                *input.committed_history,
                input.time_step,
                tangent,
                input.include_thermal_time_term);
        else
            result.residual = compute_hex8_thermoelastic(data,
                input.geometry,
                input.state,
                input.time_step > 0 ? &input.committed_state : nullptr,
                input.time_step,
                tangent,
                input.include_thermal_time_term);
    }
    if (request.history && input.committed_history)
        result.history = compute_hex8_transient_update(data,
            input.geometry,
            input.state,
            input.committed_state,
            *input.committed_history,
            input.time_step);
    if (request.history && input.committed_history)
        c3d8_detail::set_history_geometry(input, false, result);
    if (request.stress)
        result.stress = compute_hex8_stress(data, input.geometry, input.state);
    return result;
}
} // namespace fuelsim::elements
