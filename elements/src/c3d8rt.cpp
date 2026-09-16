#include "c3d8rt.hpp"
#include "ad_local_system.hpp"
#include "c3d_common.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {
using namespace c3d8_detail;
using namespace cartesian_detail;
constexpr double abaqus_total_stiffness_factor = 0.005;

struct ActiveReducedHex8Geometry final {
    adlite::Scalar volume{0.0}, center_measure{0.0};
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> average_gradient{};
    std::array<adlite::Scalar, hex8_node_count> shape_measures{};
    std::array<std::array<adlite::Scalar, 4>, hex8_node_count> hourglass_shape{};
    std::array<adlite::Scalar, 4> thermal_hourglass_coefficients{};
};

struct ReducedHex8GeometryValues final {
    ReducedHex8Metric thermal_metric;
    double volume = 0.0, center_measure = 0.0;
    std::array<std::array<double, 3>, hex8_node_count> average_gradient{};
    std::array<double, hex8_node_count> shape_measures{};
    std::array<std::array<double, 4>, hex8_node_count> hourglass_shape{};
    std::array<double, 4> thermal_hourglass_coefficients{};
};

constexpr std::size_t reduced_displacement_dof_count = 24;
using ReducedDisplacementDerivatives = std::array<double, reduced_displacement_dof_count>;

struct ReducedHex8GeometryDerivatives final {
    ReducedDisplacementDerivatives volume{}, center_measure{};
    std::array<std::array<ReducedDisplacementDerivatives, 3>, hex8_node_count> average_gradient{};
    std::array<ReducedDisplacementDerivatives, hex8_node_count> shape_measures{};
    std::array<std::array<ReducedDisplacementDerivatives, 4>, hex8_node_count> hourglass_shape{};
    std::array<ReducedDisplacementDerivatives, 4> thermal_hourglass_coefficients{};
};

struct ReducedFiniteKinematicsValues final {
    SymmetricTensor3Values strain_increment;
    cartesian_detail::Matrix3 rotation{};
};

struct ReducedFiniteStressLinearization final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 10>, 6> tangent{};
    std::array<double, 6> committed_temperature_tangent{};
};

adlite::Scalar reduced_hex8_temperature(const Hex8Geometry& geometry, const Hex8LocalAdValues& state);
double reduced_hex8_temperature(const Hex8Geometry& geometry, const Hex8LocalValues& state);
ReducedHex8GeometryValues reduced_hex8_geometry_values(const Hex8Geometry& reference,
    const std::array<std::array<double, 3>, hex8_node_count>& displacement,
    const char* configuration_name);
ReducedHex8GeometryDerivatives reduced_hex8_geometry_derivatives(const Hex8Geometry& reference,
    const std::array<std::array<double, 3>, hex8_node_count>& displacement,
    const ReducedHex8GeometryValues& values,
    double displacement_derivative_scale);
ActiveReducedHex8Geometry active_reduced_hex8_geometry(const Hex8Geometry& reference,
    const std::array<std::array<adlite::Scalar, 3>, hex8_node_count>& displacement,
    const char* configuration_name);
std::array<std::array<adlite::Scalar, 3>, hex8_node_count> reduced_hex8_displacement(const Hex8LocalAdValues& state);
std::array<std::array<adlite::Scalar, 3>, hex8_node_count>
reduced_hex8_midpoint_displacement(const Hex8LocalAdValues& state, const Hex8LocalValues& committed_state);
std::array<std::array<double, 3>, hex8_node_count> reduced_hex8_displacement_values(const Hex8LocalValues& state);
std::array<std::array<double, 3>, hex8_node_count>
reduced_hex8_midpoint_displacement_values(const Hex8LocalValues& state, const Hex8LocalValues& committed_state);
std::array<std::array<double, 3>, hex8_node_count> reduced_hex8_pushed_reference_gradient(const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    cartesian_detail::Matrix3& average_deformation);
cartesian_detail::Matrix3 reduced_hex8_average_deformation(const Hex8Geometry& reference, const Hex8LocalValues& state);
cartesian_detail::Matrix3 reduced_hex8_central_gradient_values(const Hex8Geometry& reference,
    const ReducedHex8GeometryValues& midpoint,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state);
ReducedFiniteKinematicsValues reduced_hex8_finite_kinematics_values(const cartesian_detail::Matrix3& central_gradient);
cartesian_detail::ActiveMatrix3 reduced_hex8_central_gradient(const Hex8Geometry& reference,
    const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state);
cartesian_detail::KinematicsCore reduced_hex8_finite_kinematics(const Hex8Geometry& reference,
    const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state);
CartesianStressTangent reduced_finite_material_linearization(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3Values& strain_increment,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context);
ReducedFiniteStressLinearization reduced_finite_stress_linearization(const cartesian_detail::Matrix3& passive_gradient,
    double temperature,
    const CartesianStressTangent& material_linearization);
adlite::Scalar reduced_hex8_temperature(const ActiveReducedHex8Geometry& geometry, const Hex8LocalAdValues& state);
double reduced_hex8_temperature_value(const ReducedHex8GeometryValues& geometry, const Hex8LocalValues& state);
SymmetricTensor3Values reduced_finite_stress_values(const IsotropicThermoelasticMaterial& material,
    const ReducedFiniteKinematicsValues& kinematics,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context);
Hex8LocalResidual reduced_hex8_finite_residual_values(const elements::C3d8Input& data,
    const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    double time_step,
    bool include_thermal_time_term,
    double initial_shear_modulus,
    const ReducedHex8GeometryValues& current,
    const SymmetricTensor3Values& stress);
void add_reduced_hex8_finite_jacobian(const elements::C3d8Input& data,
    const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    double time_step,
    bool include_thermal_time_term,
    double initial_shear_modulus,
    const ReducedHex8GeometryValues& current,
    const ReducedHex8GeometryDerivatives& current_derivatives,
    const ReducedHex8GeometryValues& midpoint,
    const ReducedHex8GeometryDerivatives& midpoint_derivatives,
    const ReducedFiniteStressLinearization& stress_linearization,
    Hex8LocalJacobian& jacobian);
void assemble_c3d8rt_finite_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalResidual& residual,
    Hex8LocalJacobian* jacobian);
void assemble_c3d8rt_small_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalResidual& residual,
    Hex8LocalJacobian* jacobian);
Hex8LocalResidual compute_hex8_local(const elements::C3d8Input& data, Hex8LocalJacobian* jacobian);
CartesianMaterialHistory compute_hex8_transient_update(const elements::C3d8Input& data);
std::array<SymmetricTensor3Values, 8> compute_hex8_stress(const elements::C3d8Input& data);

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

ReducedHex8GeometryValues reduced_hex8_geometry_values(const Hex8Geometry& reference,
    const std::array<std::array<double, 3>, hex8_node_count>& displacement,
    const char* configuration_name) {
    ReducedHex8GeometryValues result;
    std::array<std::array<double, 3>, hex8_node_count> current_coordinates{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        current_coordinates[node][0] = reference.capacity_points[node].position.x + displacement[node][0];
        current_coordinates[node][1] = reference.capacity_points[node].position.y + displacement[node][1];
        current_coordinates[node][2] = reference.capacity_points[node].position.z + displacement[node][2];
    }

    for (const Hex8QuadraturePoint& point : reference.points) {
        cartesian_detail::Matrix3 deformation{};
        for (std::size_t component = 0; component < 3; ++component) {
            deformation[component][component] = 1.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    deformation[component][direction] +=
                        displacement[node][component] * point.gradient[node][direction];
        }
        const double determinant = cartesian_detail::determinant(deformation);
        if (!std::isfinite(determinant) || !(determinant > 0.0))
            throw std::domain_error(std::string("C3D8RT ") + configuration_name
                                    + " configuration must preserve positive Jacobians at all integration points");
        const cartesian_detail::Matrix3 inverse = cartesian_detail::inverse(deformation, determinant);
        const double measure = point.weighted_measure * determinant;
        result.volume += measure;
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            result.shape_measures[node] += measure * point.shape[node];
            for (std::size_t current_direction = 0; current_direction < 3; ++current_direction) {
                double current_gradient = 0.0;
                for (std::size_t reference_direction = 0; reference_direction < 3; ++reference_direction)
                    current_gradient +=
                        point.gradient[node][reference_direction] * inverse[reference_direction][current_direction];
                result.average_gradient[node][current_direction] += measure * current_gradient;
            }
        }
    }
    if (!std::isfinite(result.volume) || !(result.volume > 0.0))
        throw std::domain_error(std::string("C3D8RT ") + configuration_name + " volume must be finite and positive");
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result.average_gradient[node][component] /= result.volume;

    cartesian_detail::Matrix3 center_jacobian{};
    for (std::size_t physical = 0; physical < 3; ++physical)
        for (std::size_t natural = 0; natural < 3; ++natural)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                center_jacobian[physical][natural] +=
                    current_coordinates[node][physical] * hex8_signs[node][natural] / 8.0;
    result.center_measure = 8.0 * cartesian_detail::determinant(center_jacobian);
    if (!std::isfinite(result.center_measure) || !(result.center_measure > 0.0))
        throw std::domain_error(
            std::string("C3D8RT ") + configuration_name + " center Jacobian must be finite and positive");

    result.hourglass_shape = hex8_hourglass_shape(current_coordinates, result.average_gradient);

    const cartesian_detail::Matrix3 center_inverse =
        cartesian_detail::inverse(center_jacobian, result.center_measure / 8.0);
    std::array<std::array<double, 3>, hex8_node_count> center_gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t natural = 0; natural < 3; ++natural)
                center_gradient[node][physical] += hex8_signs[node][natural] * center_inverse[natural][physical] / 8.0;
    try {
        // B528/B532 and the B544 capacity-isolated replay identify the center
        // metric and center measure for these coefficients. The uniform term
        // and hourglass projection still use the volume-averaged gradient.
        result.thermal_metric = reduced_hex8_metric(center_gradient);
        result.thermal_hourglass_coefficients =
            reduced_hex8_thermal_hourglass_coefficients(result.thermal_metric, result.center_measure);
    } catch (const std::domain_error& error) {
        throw std::domain_error(std::string("C3D8RT ") + configuration_name + " " + error.what());
    }
    return result;
}

ReducedHex8GeometryDerivatives reduced_hex8_geometry_derivatives(const Hex8Geometry& reference,
    const std::array<std::array<double, 3>, hex8_node_count>& displacement,
    const ReducedHex8GeometryValues& values,
    double displacement_derivative_scale) {
    ReducedHex8GeometryDerivatives result;
    std::array<std::array<ReducedDisplacementDerivatives, 3>, hex8_node_count> gradient_numerator_derivatives{};
    for (const Hex8QuadraturePoint& point : reference.points) {
        cartesian_detail::Matrix3 deformation{};
        for (std::size_t component = 0; component < 3; ++component) {
            deformation[component][component] = 1.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    deformation[component][direction] +=
                        displacement[node][component] * point.gradient[node][direction];
        }
        const double determinant = cartesian_detail::determinant(deformation);
        const cartesian_detail::Matrix3 inverse = cartesian_detail::inverse(deformation, determinant);
        const double measure = point.weighted_measure * determinant;
        std::array<std::array<double, 3>, hex8_node_count> current_gradient{};
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            for (std::size_t current_direction = 0; current_direction < 3; ++current_direction)
                for (std::size_t reference_direction = 0; reference_direction < 3; ++reference_direction)
                    current_gradient[node][current_direction] +=
                        point.gradient[node][reference_direction] * inverse[reference_direction][current_direction];
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t active_node = 0; active_node < hex8_node_count; ++active_node) {
                const std::size_t column = 8 * component + active_node;
                double determinant_derivative = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    determinant_derivative += inverse[direction][component] * point.gradient[active_node][direction];
                determinant_derivative *= determinant * displacement_derivative_scale;
                const double measure_derivative = point.weighted_measure * determinant_derivative;
                result.volume[column] += measure_derivative;
                for (std::size_t node = 0; node < hex8_node_count; ++node) {
                    result.shape_measures[node][column] += measure_derivative * point.shape[node];
                    for (std::size_t direction = 0; direction < 3; ++direction) {
                        const double gradient_derivative = -displacement_derivative_scale
                                                           * current_gradient[node][component]
                                                           * current_gradient[active_node][direction];
                        gradient_numerator_derivatives[node][direction][column] +=
                            measure_derivative * current_gradient[node][direction] + measure * gradient_derivative;
                    }
                }
            }
    }
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column)
                result.average_gradient[node][direction][column] =
                    (gradient_numerator_derivatives[node][direction][column]
                        - values.average_gradient[node][direction] * result.volume[column])
                    / values.volume;

    std::array<std::array<double, 3>, hex8_node_count> current_coordinates{};
    cartesian_detail::Matrix3 center_jacobian{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        current_coordinates[node] = {reference.capacity_points[node].position.x + displacement[node][0],
            reference.capacity_points[node].position.y + displacement[node][1],
            reference.capacity_points[node].position.z + displacement[node][2]};
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t natural = 0; natural < 3; ++natural)
                center_jacobian[physical][natural] +=
                    current_coordinates[node][physical] * hex8_signs[node][natural] / 8.0;
    }
    const double center_determinant = values.center_measure / 8.0;
    const cartesian_detail::Matrix3 center_inverse = cartesian_detail::inverse(center_jacobian, center_determinant);
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t active_node = 0; active_node < hex8_node_count; ++active_node) {
            const std::size_t column = 8 * component + active_node;
            for (std::size_t natural = 0; natural < 3; ++natural)
                result.center_measure[column] += displacement_derivative_scale * center_determinant
                                                 * center_inverse[natural][component]
                                                 * hex8_signs[active_node][natural];
        }

    std::array<std::array<double, 3>, 4> projected_coordinates{};
    for (std::size_t mode = 0; mode < 4; ++mode)
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                projected_coordinates[mode][component] +=
                    current_coordinates[node][component] * hex8_raw_hourglass[node][mode];
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t mode = 0; mode < 4; ++mode)
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t active_node = 0; active_node < hex8_node_count; ++active_node) {
                    const std::size_t column = 8 * component + active_node;
                    double derivative = displacement_derivative_scale * values.average_gradient[node][component]
                                        * hex8_raw_hourglass[active_node][mode];
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        derivative +=
                            result.average_gradient[node][direction][column] * projected_coordinates[mode][direction];
                    result.hourglass_shape[node][mode][column] = -derivative;
                }

    const auto& effective_mapping = values.thermal_metric.mapping;
    const auto& metric = values.thermal_metric.metric;
    const double first_pivot = values.thermal_metric.pivots[0], second_pivot = values.thermal_metric.pivots[1],
                 third_pivot = values.thermal_metric.pivots[2], leading = values.thermal_metric.leading,
                 numerator = values.thermal_metric.numerator;
    const double inverse_x = 1.0 / first_pivot, inverse_y = 1.0 / second_pivot, inverse_z = 1.0 / third_pivot;
    const double thermal_scale = values.center_measure / 192.0;
    const std::array<double, 4> sums = {inverse_x + inverse_y,
        inverse_x + inverse_z,
        inverse_x + inverse_z,
        (inverse_x + inverse_y + inverse_z) / 3.0};
    for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
        // With center gradients the effective mapping is the center Jacobian:
        // sum(sign_i * grad_center_i) = inverse(J_center). Differentiate that
        // mapping directly instead of applying a volume-average gradient chain.
        const std::size_t component = column / hex8_node_count, active_node = column % hex8_node_count;
        cartesian_detail::Matrix3 effective_derivative{};
        for (std::size_t natural = 0; natural < 3; ++natural)
            effective_derivative[component][natural] =
                displacement_derivative_scale * hex8_signs[active_node][natural] / 8.0;
        cartesian_detail::Matrix3 metric_derivative{};
        for (std::size_t first = 0; first < 3; ++first)
            for (std::size_t second = 0; second < 3; ++second)
                for (std::size_t physical = 0; physical < 3; ++physical) {
                    metric_derivative[first][second] +=
                        effective_derivative[physical][first] * effective_mapping[physical][second]
                        + effective_mapping[physical][first] * effective_derivative[physical][second];
                }
        const double first_pivot_derivative = metric_derivative[0][0];
        const double second_pivot_derivative =
            metric_derivative[1][1] - 2.0 * metric[0][1] * metric_derivative[0][1] / first_pivot
            + metric[0][1] * metric[0][1] * first_pivot_derivative / (first_pivot * first_pivot);
        const double leading_derivative = metric_derivative[0][0] * metric[1][1]
                                          + metric[0][0] * metric_derivative[1][1]
                                          - 2.0 * metric[0][1] * metric_derivative[0][1];
        const double numerator_derivative = metric_derivative[1][1] * metric[0][2] * metric[0][2]
                                            + 2.0 * metric[1][1] * metric[0][2] * metric_derivative[0][2]
                                            - 2.0
                                                  * (metric_derivative[0][1] * metric[0][2] * metric[1][2]
                                                      + metric[0][1] * metric_derivative[0][2] * metric[1][2]
                                                      + metric[0][1] * metric[0][2] * metric_derivative[1][2])
                                            + metric_derivative[0][0] * metric[1][2] * metric[1][2]
                                            + 2.0 * metric[0][0] * metric[1][2] * metric_derivative[1][2];
        const double third_pivot_derivative = metric_derivative[2][2] - numerator_derivative / leading
                                              + numerator * leading_derivative / (leading * leading);
        const double inverse_x_derivative = -first_pivot_derivative * inverse_x * inverse_x;
        const double inverse_y_derivative = -second_pivot_derivative * inverse_y * inverse_y;
        const double inverse_z_derivative = -third_pivot_derivative * inverse_z * inverse_z;
        const double thermal_scale_derivative = result.center_measure[column] / 192.0;
        const std::array<double, 4> sum_derivatives = {inverse_x_derivative + inverse_y_derivative,
            inverse_x_derivative + inverse_z_derivative,
            inverse_x_derivative + inverse_z_derivative,
            (inverse_x_derivative + inverse_y_derivative + inverse_z_derivative) / 3.0};
        for (std::size_t mode = 0; mode < 4; ++mode)
            result.thermal_hourglass_coefficients[mode][column] =
                thermal_scale_derivative * sums[mode] + thermal_scale * sum_derivatives[mode];
    }
    return result;
}

ActiveReducedHex8Geometry active_reduced_hex8_geometry(const Hex8Geometry& reference,
    const std::array<std::array<adlite::Scalar, 3>, hex8_node_count>& displacement,
    const char* configuration_name) {
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
            throw std::domain_error(std::string("C3D8RT ") + configuration_name
                                    + " configuration must preserve positive Jacobians at all integration points");
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

    result.hourglass_shape = hex8_hourglass_shape(current_coordinates, result.average_gradient);

    const cartesian_detail::ActiveMatrix3 center_inverse =
        cartesian_detail::inverse(center_jacobian, result.center_measure / 8.0);
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> center_gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t natural = 0; natural < 3; ++natural)
                center_gradient[node][physical] += hex8_signs[node][natural] * center_inverse[natural][physical] / 8.0;
    try {
        result.thermal_hourglass_coefficients =
            reduced_hex8_thermal_hourglass_coefficients(center_gradient, result.center_measure);
    } catch (const std::domain_error& error) {
        throw std::domain_error(std::string("C3D8RT ") + configuration_name + " " + error.what());
    }
    return result;
}

std::array<std::array<adlite::Scalar, 3>, hex8_node_count> reduced_hex8_displacement(const Hex8LocalAdValues& state) {
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> result{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[node][component] = state[8 * (component + 1) + node];
    return result;
}

std::array<std::array<adlite::Scalar, 3>, hex8_node_count>
reduced_hex8_midpoint_displacement(const Hex8LocalAdValues& state, const Hex8LocalValues& committed_state) {
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> result{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[node][component] =
                0.5 * (state[8 * (component + 1) + node] + committed_state[8 * (component + 1) + node]);
    return result;
}

std::array<std::array<double, 3>, hex8_node_count> reduced_hex8_displacement_values(const Hex8LocalValues& state) {
    std::array<std::array<double, 3>, hex8_node_count> result{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[node][component] = state[8 * (component + 1) + node];
    return result;
}

std::array<std::array<double, 3>, hex8_node_count>
reduced_hex8_midpoint_displacement_values(const Hex8LocalValues& state, const Hex8LocalValues& committed_state) {
    std::array<std::array<double, 3>, hex8_node_count> result{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[node][component] =
                0.5 * (state[8 * (component + 1) + node] + committed_state[8 * (component + 1) + node]);
    return result;
}

cartesian_detail::Matrix3 reduced_hex8_central_gradient_values(const Hex8Geometry& reference,
    const ReducedHex8GeometryValues& midpoint,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state) {
    (void)reduced_hex8_average_deformation(reference, state);
    (void)reduced_hex8_average_deformation(reference, committed_state);
    Hex8LocalValues midpoint_state{};
    for (std::size_t dof = 8; dof < hex8_local_dof_count; ++dof)
        midpoint_state[dof] = 0.5 * (state[dof] + committed_state[dof]);
    cartesian_detail::Matrix3 midpoint_deformation{};
    const auto pushed_gradient =
        reduced_hex8_pushed_reference_gradient(reference, midpoint_state, midpoint_deformation);
    cartesian_detail::Matrix3 central_gradient{};
    double midpoint_trace = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                const double increment =
                    state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node];
                central_gradient[component][direction] += increment * pushed_gradient[node][direction];
                if (component == direction)
                    midpoint_trace += increment * midpoint.average_gradient[node][direction];
            }
    const double trace_correction =
        (midpoint_trace - central_gradient[0][0] - central_gradient[1][1] - central_gradient[2][2]) / 3.0;
    for (std::size_t direction = 0; direction < 3; ++direction)
        central_gradient[direction][direction] += trace_correction;
    return central_gradient;
}

ReducedFiniteKinematicsValues reduced_hex8_finite_kinematics_values(const cartesian_detail::Matrix3& central_gradient) {
    ReducedFiniteKinematicsValues result;
    std::array<double, 6> strain{};
    cartesian_detail::hughes_winget_rotation(central_gradient, result.rotation, strain);
    result.strain_increment = {strain[0], strain[1], strain[2], strain[3], strain[4], strain[5]};
    return result;
}

cartesian_detail::ActiveMatrix3 reduced_hex8_central_gradient(const Hex8Geometry& reference,
    const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state) {
    cartesian_detail::ActiveMatrix3 midpoint_deformation{}, current_deformation{};
    cartesian_detail::Matrix3 committed_deformation{};
    for (std::size_t component = 0; component < 3; ++component) {
        midpoint_deformation[component][component] = 1.0;
        current_deformation[component][component] = 1.0;
        committed_deformation[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                midpoint_deformation[component][direction] +=
                    0.5 * (state[8 * (component + 1) + node] + committed_state[8 * (component + 1) + node])
                    * reference.average_shape_gradient[node][direction];
                current_deformation[component][direction] +=
                    state[8 * (component + 1) + node] * reference.average_shape_gradient[node][direction];
                committed_deformation[component][direction] +=
                    committed_state[8 * (component + 1) + node] * reference.average_shape_gradient[node][direction];
            }
    }
    const double current_determinant = cartesian_detail::determinant(current_deformation).value();
    const double committed_determinant = cartesian_detail::determinant(committed_deformation);
    if (!std::isfinite(current_determinant) || !(current_determinant > 0.0) || !std::isfinite(committed_determinant)
        || !(committed_determinant > 0.0))
        throw std::domain_error("C3D8RT current and committed average deformations require positive Jacobians");
    const adlite::Scalar determinant = cartesian_detail::determinant(midpoint_deformation);
    if (!std::isfinite(determinant.value()) || !(determinant.value() > 0.0))
        throw std::domain_error("C3D8RT midpoint average deformation must preserve a finite positive Jacobian");
    const auto inverse = cartesian_detail::inverse(midpoint_deformation, determinant);
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> pushed_gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t current_direction = 0; current_direction < 3; ++current_direction)
            for (std::size_t reference_direction = 0; reference_direction < 3; ++reference_direction)
                pushed_gradient[node][current_direction] += reference.average_shape_gradient[node][reference_direction]
                                                            * inverse[reference_direction][current_direction];
    cartesian_detail::ActiveMatrix3 central_gradient{};
    adlite::Scalar midpoint_trace = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                const adlite::Scalar increment =
                    state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node];
                central_gradient[component][direction] += increment * pushed_gradient[node][direction];
                if (component == direction)
                    midpoint_trace += increment * midpoint.average_gradient[node][direction];
            }
    const adlite::Scalar trace_correction =
        (midpoint_trace - central_gradient[0][0] - central_gradient[1][1] - central_gradient[2][2]) / 3.0;
    for (std::size_t direction = 0; direction < 3; ++direction)
        central_gradient[direction][direction] += trace_correction;
    return central_gradient;
}

cartesian_detail::KinematicsCore reduced_hex8_finite_kinematics(const Hex8Geometry& reference,
    const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state) {
    return cartesian_detail::evaluate_hughes_winget_increment(
        reduced_hex8_central_gradient(reference, midpoint, state, committed_state));
}

CartesianStressTangent reduced_finite_material_linearization(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3Values& strain_increment,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context) {
    std::array<double, 7> seeds = {strain_increment.xx,
        strain_increment.yy,
        strain_increment.zz,
        strain_increment.xy,
        strain_increment.yz,
        strain_increment.xz,
        temperature};
    std::array<adlite::Scalar, 7> active{};
    adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    const SymmetricTensor3 strain{active[0], active[1], active[2], active[3], active[4], active[5]};
    SymmetricTensor3 stress;
    if (committed_material == nullptr)
        stress = material.stress(strain, active[6], context);
    else
        stress = material
                     .incremental_response(strain,
                         CartesianRotation{},
                         active[6],
                         committed_temperature,
                         time_step,
                         *committed_material,
                         context)
                     .stress;
    CartesianStressTangent result;
    result.stress = {stress.xx.value(),
        stress.yy.value(),
        stress.zz.value(),
        stress.xy.value(),
        stress.yz.value(),
        stress.xz.value()};
    const std::array<const adlite::Scalar*, 6> components =
        {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz};
    std::array<double, 7> derivatives{};
    for (std::size_t row = 0; row < 6; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 6; ++column)
            result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[6];
    }
    return result;
}

ReducedFiniteStressLinearization reduced_finite_stress_linearization(const cartesian_detail::Matrix3& passive_gradient,
    double temperature,
    const CartesianStressTangent& material_linearization) {
    std::array<double, 10> values{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            values[3 * component + direction] = passive_gradient[component][direction];
    values[9] = temperature;
    std::array<adlite::Scalar, 10> active{};
    adlite::seed_identity(values.data(), values.size(), active.data());
    cartesian_detail::ActiveMatrix3 gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            gradient[component][direction] = active[3 * component + direction];
    const cartesian_detail::KinematicsCore kinematics = cartesian_detail::evaluate_hughes_winget_increment(gradient);
    const SymmetricTensor3 rotated = rotate_cartesian_tensor(compose_cartesian_stress(material_linearization,
                                                                 kinematics.strain_increment,
                                                                 active[9],
                                                                 material_linearization.thermal),
        kinematics.rotation);
    const std::array<const adlite::Scalar*, 6> components =
        {&rotated.xx, &rotated.yy, &rotated.zz, &rotated.xy, &rotated.yz, &rotated.xz};
    ReducedFiniteStressLinearization result;
    result.stress = {rotated.xx.value(),
        rotated.yy.value(),
        rotated.zz.value(),
        rotated.xy.value(),
        rotated.yz.value(),
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

double reduced_hex8_temperature_value(const ReducedHex8GeometryValues& geometry, const Hex8LocalValues& state) {
    double result = 0.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        result += geometry.shape_measures[node] * state[node] / geometry.volume;
    return result;
}

SymmetricTensor3Values reduced_finite_stress_values(const IsotropicThermoelasticMaterial& material,
    const ReducedFiniteKinematicsValues& kinematics,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context) {
    const SymmetricTensor3 strain{kinematics.strain_increment.xx,
        kinematics.strain_increment.yy,
        kinematics.strain_increment.zz,
        kinematics.strain_increment.xy,
        kinematics.strain_increment.yz,
        kinematics.strain_increment.xz};
    const adlite::Scalar active_temperature(temperature);
    SymmetricTensor3 stress;
    if (committed_material == nullptr)
        stress = material.stress(strain, active_temperature, context);
    else
        stress = evaluate_incremental_cartesian_response(material,
            strain,
            active_temperature,
            committed_temperature,
            time_step,
            *committed_material,
            context)
                     .stress;
    return rotate_cartesian_tensor_values({stress.xx.value(),
                                              stress.yy.value(),
                                              stress.zz.value(),
                                              stress.xy.value(),
                                              stress.yz.value(),
                                              stress.xz.value()},
        kinematics.rotation);
}

cartesian_detail::Matrix3 reduced_hex8_average_deformation(const Hex8Geometry& reference,
    const Hex8LocalValues& state) {
    cartesian_detail::Matrix3 average_deformation{};
    for (std::size_t component = 0; component < 3; ++component) {
        average_deformation[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                average_deformation[component][direction] +=
                    state[8 * (component + 1) + node] * reference.average_shape_gradient[node][direction];
    }
    const double determinant = cartesian_detail::determinant(average_deformation);
    if (!std::isfinite(determinant) || !(determinant > 0.0))
        throw std::domain_error("C3D8RT average deformation must preserve a finite positive Jacobian");
    return average_deformation;
}

std::array<std::array<double, 3>, hex8_node_count> reduced_hex8_pushed_reference_gradient(const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    cartesian_detail::Matrix3& average_deformation) {
    average_deformation = reduced_hex8_average_deformation(reference, state);
    const double determinant = cartesian_detail::determinant(average_deformation);
    const cartesian_detail::Matrix3 inverse = cartesian_detail::inverse(average_deformation, determinant);
    std::array<std::array<double, 3>, hex8_node_count> gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t current_direction = 0; current_direction < 3; ++current_direction)
            for (std::size_t reference_direction = 0; reference_direction < 3; ++reference_direction)
                gradient[node][current_direction] += reference.average_shape_gradient[node][reference_direction]
                                                     * inverse[reference_direction][current_direction];
    return gradient;
}

Hex8LocalResidual reduced_hex8_finite_residual_values(const elements::C3d8Input& data,
    const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    double time_step,
    bool include_thermal_time_term,
    double initial_shear_modulus,
    const ReducedHex8GeometryValues& current,
    const SymmetricTensor3Values& stress) {
    const double temperature = reduced_hex8_temperature_value(current, state);
    const MaterialFunctionContext context = material_context(data.time, reference.reduced_point.position);
    const double conductivity = data.material.conductivity(adlite::Scalar(temperature), context).value();
    std::array<double, 3> temperature_gradient{};
    std::array<double, 4> temperature_hourglass_amplitude{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        for (std::size_t direction = 0; direction < 3; ++direction)
            temperature_gradient[direction] += current.average_gradient[node][direction] * state[node];
        for (std::size_t mode = 0; mode < 4; ++mode)
            temperature_hourglass_amplitude[mode] += current.hourglass_shape[node][mode] * state[node];
    }

    cartesian_detail::Matrix3 average_deformation{};
    const auto pushed_gradient = reduced_hex8_pushed_reference_gradient(reference, state, average_deformation);
    const double pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
    Hex8LocalResidual residual{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        double uniform_thermal = 0.0, hourglass_thermal = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            uniform_thermal += current.average_gradient[node][direction] * temperature_gradient[direction];
        for (std::size_t mode = 0; mode < 4; ++mode)
            hourglass_thermal += current.hourglass_shape[node][mode] * current.thermal_hourglass_coefficients[mode]
                                 * temperature_hourglass_amplitude[mode];
        residual[node] += conductivity * (current.volume * uniform_thermal + hourglass_thermal)
                          - current.center_measure * data.volumetric_heat_source / 8.0;
        // Native B544 replay distinguishes the deviatoric and volumetric virtual gradients.
        const double gradient_x = pushed_gradient[node][0];
        const double gradient_y = pushed_gradient[node][1];
        const double gradient_z = pushed_gradient[node][2];
        residual[8 + node] += current.volume
                              * ((stress.xx - pressure) * gradient_x + stress.xy * gradient_y + stress.xz * gradient_z
                                  + pressure * current.average_gradient[node][0]);
        residual[16 + node] += current.volume
                               * (stress.xy * gradient_x + (stress.yy - pressure) * gradient_y + stress.yz * gradient_z
                                   + pressure * current.average_gradient[node][1]);
        residual[24 + node] += current.volume
                               * (stress.xz * gradient_x + stress.yz * gradient_y + (stress.zz - pressure) * gradient_z
                                   + pressure * current.average_gradient[node][2]);
    }

    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<double, 3> reference_amplitude{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                reference_amplitude[component] +=
                    reference.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
        std::array<double, 3> material_modal_force{};
        for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
            double transported_amplitude = 0.0;
            for (std::size_t component = 0; component < 3; ++component)
                transported_amplitude +=
                    reference_amplitude[component] * average_deformation[component][material_direction];
            material_modal_force[material_direction] = abaqus_total_stiffness_factor * initial_shear_modulus
                                                       * reference.mechanical_hourglass_metrics[material_direction]
                                                       * transported_amplitude;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            double spatial_modal_force = 0.0;
            for (std::size_t material_direction = 0; material_direction < 3; ++material_direction)
                spatial_modal_force +=
                    average_deformation[component][material_direction] * material_modal_force[material_direction];
            for (std::size_t node = 0; node < hex8_node_count; ++node) {
                double deformation_gradient_term = 0.0;
                for (std::size_t material_direction = 0; material_direction < 3; ++material_direction)
                    deformation_gradient_term += reference.average_shape_gradient[node][material_direction]
                                                 * material_modal_force[material_direction];
                residual[8 * (component + 1) + node] += reference.hourglass_shape[node][mode] * spatial_modal_force
                                                        + reference_amplitude[component] * deformation_gradient_term;
            }
        }
    }

    if (include_thermal_time_term)
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            const Hex8CapacityPoint& capacity_point = reference.reduced_capacity_points[node];
            const double capacity = data.material
                                        .reference_heat_capacity(adlite::Scalar(state[node]),
                                            data.initial_temperature,
                                            material_context(data.time, capacity_point.position))
                                        .value();
            residual[node] +=
                capacity_point.weighted_measure * capacity * (state[node] - committed_state[node]) / time_step;
        }
    return residual;
}

void add_reduced_hex8_finite_jacobian(const elements::C3d8Input& data,
    const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    double time_step,
    bool include_thermal_time_term,
    double initial_shear_modulus,
    const ReducedHex8GeometryValues& current,
    const ReducedHex8GeometryDerivatives& current_derivatives,
    const ReducedHex8GeometryValues& midpoint,
    const ReducedHex8GeometryDerivatives& midpoint_derivatives,
    const ReducedFiniteStressLinearization& stress_linearization,
    Hex8LocalJacobian& jacobian) {
    jacobian.fill(0.0);
    std::array<ReducedDisplacementDerivatives, 9> central_gradient_derivatives{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
                double derivative = 0.0;
                const std::size_t column_component = column / 8, column_node = column % 8;
                if (component == column_component)
                    derivative += midpoint.average_gradient[column_node][direction];
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    derivative += (state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node])
                                  * midpoint_derivatives.average_gradient[node][direction][column];
                central_gradient_derivatives[3 * component + direction][column] = derivative;
            }

    Hex8LocalValues midpoint_state{};
    for (std::size_t dof = 8; dof < hex8_local_dof_count; ++dof)
        midpoint_state[dof] = 0.5 * (state[dof] + committed_state[dof]);
    cartesian_detail::Matrix3 midpoint_deformation{};
    const auto midpoint_pushed_gradient =
        reduced_hex8_pushed_reference_gradient(reference, midpoint_state, midpoint_deformation);
    cartesian_detail::Matrix3 pushed_increment{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                pushed_increment[component][direction] +=
                    (state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node])
                    * midpoint_pushed_gradient[node][direction];
    for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
        const std::size_t column_component = column / 8, column_node = column % 8;
        const double midpoint_trace_derivative = central_gradient_derivatives[0][column]
                                                 + central_gradient_derivatives[4][column]
                                                 + central_gradient_derivatives[8][column];
        double pushed_trace_derivative = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t direction = 0; direction < 3; ++direction) {
                const double derivative =
                    ((component == column_component ? 1.0 : 0.0) - 0.5 * pushed_increment[component][column_component])
                    * midpoint_pushed_gradient[column_node][direction];
                central_gradient_derivatives[3 * component + direction][column] = derivative;
                if (component == direction)
                    pushed_trace_derivative += derivative;
            }
        const double trace_correction = (midpoint_trace_derivative - pushed_trace_derivative) / 3.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            central_gradient_derivatives[4 * direction][column] += trace_correction;
    }

    const double temperature = reduced_hex8_temperature_value(current, state);
    const double committed_temperature = reduced_hex8_temperature_value(current, committed_state);
    ReducedDisplacementDerivatives temperature_displacement_derivatives{};
    ReducedDisplacementDerivatives committed_temperature_displacement_derivatives{};
    std::array<double, hex8_node_count> temperature_temperature_derivatives{};
    for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
        double numerator_derivative = 0.0, committed_numerator_derivative = 0.0;
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            numerator_derivative += current_derivatives.shape_measures[node][column] * state[node];
            committed_numerator_derivative += current_derivatives.shape_measures[node][column] * committed_state[node];
        }
        temperature_displacement_derivatives[column] =
            (numerator_derivative - temperature * current_derivatives.volume[column]) / current.volume;
        committed_temperature_displacement_derivatives[column] =
            (committed_numerator_derivative - committed_temperature * current_derivatives.volume[column])
            / current.volume;
    }
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        temperature_temperature_derivatives[node] = current.shape_measures[node] / current.volume;

    const std::array<double, 6> stress_values = {stress_linearization.stress.xx,
        stress_linearization.stress.yy,
        stress_linearization.stress.zz,
        stress_linearization.stress.xy,
        stress_linearization.stress.yz,
        stress_linearization.stress.xz};
    std::array<ReducedDisplacementDerivatives, 6> stress_displacement_derivatives{};
    std::array<std::array<double, hex8_node_count>, 6> stress_temperature_derivatives{};
    for (std::size_t stress_component = 0; stress_component < 6; ++stress_component) {
        for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
            double derivative =
                stress_linearization.tangent[stress_component][9] * temperature_displacement_derivatives[column]
                + stress_linearization.committed_temperature_tangent[stress_component]
                      * committed_temperature_displacement_derivatives[column];
            for (std::size_t gradient = 0; gradient < 9; ++gradient)
                derivative += stress_linearization.tangent[stress_component][gradient]
                              * central_gradient_derivatives[gradient][column];
            stress_displacement_derivatives[stress_component][column] = derivative;
        }
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            stress_temperature_derivatives[stress_component][node] =
                stress_linearization.tangent[stress_component][9] * temperature_temperature_derivatives[node];
    }

    cartesian_detail::Matrix3 average_deformation{};
    const auto pushed_gradient = reduced_hex8_pushed_reference_gradient(reference, state, average_deformation);
    const double pressure = (stress_values[0] + stress_values[1] + stress_values[2]) / 3.0;
    ReducedDisplacementDerivatives pressure_displacement_derivatives{};
    std::array<double, hex8_node_count> pressure_temperature_derivatives{};
    for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column)
        pressure_displacement_derivatives[column] =
            (stress_displacement_derivatives[0][column] + stress_displacement_derivatives[1][column]
                + stress_displacement_derivatives[2][column])
            / 3.0;
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        pressure_temperature_derivatives[node] =
            (stress_temperature_derivatives[0][node] + stress_temperature_derivatives[1][node]
                + stress_temperature_derivatives[2][node])
            / 3.0;
    const std::array<std::array<double, 3>, 3> deviatoric_stress = {
        {{stress_values[0] - pressure, stress_values[3], stress_values[5]},
            {stress_values[3], stress_values[1] - pressure, stress_values[4]},
            {stress_values[5], stress_values[4], stress_values[2] - pressure}}};
    constexpr std::array<std::array<std::size_t, 3>, 3> stress_indices = {{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};

    const MaterialFunctionContext context = material_context(data.time, reference.reduced_point.position);
    const adlite::Scalar active_temperature = adlite::Scalar::independent(temperature, 0, 1);
    const adlite::Scalar active_conductivity = data.material.conductivity(active_temperature, context);
    const double conductivity = active_conductivity.value();
    const double conductivity_temperature_derivative =
        active_conductivity.is_active() ? active_conductivity.derivative(0) : 0.0;
    std::array<double, 3> temperature_gradient{};
    std::array<ReducedDisplacementDerivatives, 3> temperature_gradient_displacement_derivatives{};
    std::array<double, 4> temperature_hourglass_amplitude{};
    std::array<ReducedDisplacementDerivatives, 4> amplitude_displacement_derivatives{};
    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        for (std::size_t direction = 0; direction < 3; ++direction) {
            temperature_gradient[direction] += current.average_gradient[node][direction] * state[node];
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column)
                temperature_gradient_displacement_derivatives[direction][column] +=
                    current_derivatives.average_gradient[node][direction][column] * state[node];
        }
        for (std::size_t mode = 0; mode < 4; ++mode) {
            temperature_hourglass_amplitude[mode] += current.hourglass_shape[node][mode] * state[node];
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column)
                amplitude_displacement_derivatives[mode][column] +=
                    current_derivatives.hourglass_shape[node][mode][column] * state[node];
        }
    }

    for (std::size_t node = 0; node < hex8_node_count; ++node) {
        double uniform_thermal = 0.0, hourglass_thermal = 0.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            uniform_thermal += current.average_gradient[node][direction] * temperature_gradient[direction];
        for (std::size_t mode = 0; mode < 4; ++mode)
            hourglass_thermal += current.hourglass_shape[node][mode] * current.thermal_hourglass_coefficients[mode]
                                 * temperature_hourglass_amplitude[mode];
        const double conduction_measure = current.volume * uniform_thermal + hourglass_thermal;
        for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
            double uniform_derivative = 0.0, hourglass_derivative = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                uniform_derivative +=
                    current_derivatives.average_gradient[node][direction][column] * temperature_gradient[direction]
                    + current.average_gradient[node][direction]
                          * temperature_gradient_displacement_derivatives[direction][column];
            for (std::size_t mode = 0; mode < 4; ++mode)
                hourglass_derivative += (current_derivatives.hourglass_shape[node][mode][column]
                                                * current.thermal_hourglass_coefficients[mode]
                                            + current.hourglass_shape[node][mode]
                                                  * current_derivatives.thermal_hourglass_coefficients[mode][column])
                                            * temperature_hourglass_amplitude[mode]
                                        + current.hourglass_shape[node][mode]
                                              * current.thermal_hourglass_coefficients[mode]
                                              * amplitude_displacement_derivatives[mode][column];
            const double measure_derivative = current_derivatives.volume[column] * uniform_thermal
                                              + current.volume * uniform_derivative + hourglass_derivative;
            jacobian[node * hex8_local_dof_count + 8 + column] +=
                conductivity_temperature_derivative * temperature_displacement_derivatives[column] * conduction_measure
                + conductivity * measure_derivative
                - current_derivatives.center_measure[column] * data.volumetric_heat_source / 8.0;
        }
        for (std::size_t temperature_node = 0; temperature_node < hex8_node_count; ++temperature_node) {
            double uniform_derivative = 0.0, hourglass_derivative = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                uniform_derivative +=
                    current.average_gradient[node][direction] * current.average_gradient[temperature_node][direction];
            for (std::size_t mode = 0; mode < 4; ++mode)
                hourglass_derivative += current.hourglass_shape[node][mode]
                                        * current.thermal_hourglass_coefficients[mode]
                                        * current.hourglass_shape[temperature_node][mode];
            jacobian[node * hex8_local_dof_count + temperature_node] +=
                conductivity_temperature_derivative * temperature_temperature_derivatives[temperature_node]
                    * conduction_measure
                + conductivity * (current.volume * uniform_derivative + hourglass_derivative);
        }

        for (std::size_t row_component = 0; row_component < 3; ++row_component) {
            const std::size_t row = 8 * (row_component + 1) + node;
            double force_density = pressure * current.average_gradient[node][row_component];
            for (std::size_t direction = 0; direction < 3; ++direction)
                force_density += deviatoric_stress[row_component][direction] * pushed_gradient[node][direction];
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
                const std::size_t column_component = column / 8, column_node = column % 8;
                const double pressure_derivative = pressure_displacement_derivatives[column];
                double force_density_derivative =
                    pressure_derivative * current.average_gradient[node][row_component]
                    + pressure * current_derivatives.average_gradient[node][row_component][column];
                for (std::size_t direction = 0; direction < 3; ++direction) {
                    const std::size_t stress_index = stress_indices[row_component][direction];
                    double stress_derivative = stress_displacement_derivatives[stress_index][column];
                    if (row_component == direction)
                        stress_derivative -= pressure_derivative;
                    const double gradient_derivative =
                        -pushed_gradient[node][column_component] * pushed_gradient[column_node][direction];
                    force_density_derivative += stress_derivative * pushed_gradient[node][direction]
                                                + deviatoric_stress[row_component][direction] * gradient_derivative;
                }
                jacobian[row * hex8_local_dof_count + 8 + column] +=
                    current_derivatives.volume[column] * force_density + current.volume * force_density_derivative;
            }
            for (std::size_t temperature_node = 0; temperature_node < hex8_node_count; ++temperature_node) {
                const double pressure_derivative = pressure_temperature_derivatives[temperature_node];
                double derivative = pressure_derivative * current.average_gradient[node][row_component];
                for (std::size_t direction = 0; direction < 3; ++direction) {
                    const std::size_t stress_index = stress_indices[row_component][direction];
                    double stress_derivative = stress_temperature_derivatives[stress_index][temperature_node];
                    if (row_component == direction)
                        stress_derivative -= pressure_derivative;
                    derivative += stress_derivative * pushed_gradient[node][direction];
                }
                jacobian[row * hex8_local_dof_count + temperature_node] += current.volume * derivative;
            }
        }
    }

    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<double, 3> reference_amplitude{}, material_modal_force{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                reference_amplitude[component] +=
                    reference.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
        for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
            double transported_amplitude = 0.0;
            for (std::size_t component = 0; component < 3; ++component)
                transported_amplitude +=
                    reference_amplitude[component] * average_deformation[component][material_direction];
            material_modal_force[material_direction] = abaqus_total_stiffness_factor * initial_shear_modulus
                                                       * reference.mechanical_hourglass_metrics[material_direction]
                                                       * transported_amplitude;
        }
        for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
            const std::size_t column_component = column / 8, column_node = column % 8;
            std::array<double, 3> material_modal_force_derivative{};
            for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
                double transported_derivative = reference.hourglass_shape[column_node][mode]
                                                * average_deformation[column_component][material_direction];
                transported_derivative += reference_amplitude[column_component]
                                          * reference.average_shape_gradient[column_node][material_direction];
                material_modal_force_derivative[material_direction] =
                    abaqus_total_stiffness_factor * initial_shear_modulus
                    * reference.mechanical_hourglass_metrics[material_direction] * transported_derivative;
            }
            for (std::size_t row_component = 0; row_component < 3; ++row_component) {
                double spatial_force = 0.0, spatial_force_derivative = 0.0;
                for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
                    spatial_force += average_deformation[row_component][material_direction]
                                     * material_modal_force[material_direction];
                    spatial_force_derivative += average_deformation[row_component][material_direction]
                                                * material_modal_force_derivative[material_direction];
                    if (row_component == column_component)
                        spatial_force_derivative += reference.average_shape_gradient[column_node][material_direction]
                                                    * material_modal_force[material_direction];
                }
                for (std::size_t node = 0; node < hex8_node_count; ++node) {
                    double deformation_term = 0.0, deformation_term_derivative = 0.0;
                    for (std::size_t material_direction = 0; material_direction < 3; ++material_direction) {
                        deformation_term += reference.average_shape_gradient[node][material_direction]
                                            * material_modal_force[material_direction];
                        deformation_term_derivative += reference.average_shape_gradient[node][material_direction]
                                                       * material_modal_force_derivative[material_direction];
                    }
                    double derivative = reference.hourglass_shape[node][mode] * spatial_force_derivative
                                        + reference_amplitude[row_component] * deformation_term_derivative;
                    if (row_component == column_component)
                        derivative += reference.hourglass_shape[column_node][mode] * deformation_term;
                    jacobian[(8 * (row_component + 1) + node) * hex8_local_dof_count + 8 + column] += derivative;
                }
            }
        }
    }

    if (include_thermal_time_term)
        for (std::size_t node = 0; node < hex8_node_count; ++node) {
            const adlite::Scalar nodal_temperature = adlite::Scalar::independent(state[node], 0, 1);
            const adlite::Scalar capacity = data.material.reference_heat_capacity(nodal_temperature,
                data.initial_temperature,
                material_context(data.time, reference.reduced_capacity_points[node].position));
            const double capacity_derivative = capacity.is_active() ? capacity.derivative(0) : 0.0;
            jacobian[node * hex8_local_dof_count + node] +=
                reference.reduced_capacity_points[node].weighted_measure
                * (capacity.value() + capacity_derivative * (state[node] - committed_state[node])) / time_step;
        }
}

void assemble_c3d8rt_finite_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalResidual& residual,
    Hex8LocalJacobian* jacobian) {
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const ReducedHex8GeometryValues current =
        reduced_hex8_geometry_values(geometry, reduced_hex8_displacement_values(state), "current");
    const ReducedHex8GeometryValues midpoint =
        reduced_hex8_geometry_values(geometry, reduced_hex8_midpoint_displacement_values(state, old_state), "midpoint");
    const cartesian_detail::Matrix3 central_gradient =
        reduced_hex8_central_gradient_values(geometry, midpoint, state, old_state);
    const ReducedFiniteKinematicsValues kinematics = reduced_hex8_finite_kinematics_values(central_gradient);
    const double temperature = reduced_hex8_temperature_value(current, state);
    double old_temperature = data.initial_temperature;
    if (committed_material != nullptr) {
        (void)reduced_hex8_geometry_values(geometry, reduced_hex8_displacement_values(old_state), "committed");
        old_temperature = reduced_hex8_temperature_value(current, old_state);
    }
    const MaterialFunctionContext context = material_context(data.time, geometry.reduced_point.position);
    const SymmetricTensor3Values stress = reduced_finite_stress_values(data.material,
        kinematics,
        temperature,
        old_temperature,
        time_step,
        committed_material,
        context);
    const ActiveThermoelasticProperties initial_properties =
        data.material.active_properties(adlite::Scalar(data.initial_temperature),
            material_context(0.0, geometry.reduced_point.position));
    const double initial_shear_modulus = initial_properties.shear_modulus.value();
    if (!std::isfinite(initial_shear_modulus) || !(initial_shear_modulus > 0.0))
        throw std::invalid_argument("C3D8RT initial shear modulus must be finite and positive");
    const bool add_thermal_time_term = committed_state != nullptr && include_thermal_time_term;
    residual = reduced_hex8_finite_residual_values(data,
        geometry,
        state,
        old_state,
        time_step,
        add_thermal_time_term,
        initial_shear_modulus,
        current,
        stress);
    if (jacobian == nullptr)
        return;

    jacobian->fill(0.0);
    const CartesianStressTangent material_linearization = reduced_finite_material_linearization(data.material,
        kinematics.strain_increment,
        temperature,
        old_temperature,
        time_step,
        committed_material,
        context);
    ReducedFiniteStressLinearization stress_linearization =
        reduced_finite_stress_linearization(central_gradient, temperature, material_linearization);
    if (committed_material != nullptr) {
        MaterialFunctionContext old_context = context;
        old_context.time -= time_step;
        const SymmetricTensor3 old_eigenstrain =
            data.material.eigenstrain(adlite::Scalar::independent(old_temperature, 0, 1), old_context);
        const std::array<const adlite::Scalar*, 6> components = {&old_eigenstrain.xx,
            &old_eigenstrain.yy,
            &old_eigenstrain.zz,
            &old_eigenstrain.xy,
            &old_eigenstrain.yz,
            &old_eigenstrain.xz};
        std::array<double, 6> eigenstrain_derivative{}, material_stress_derivative{};
        for (std::size_t component = 0; component < 6; ++component)
            eigenstrain_derivative[component] =
                components[component]->is_active() ? components[component]->derivative(0) : 0.0;
        for (std::size_t row = 0; row < 6; ++row)
            for (std::size_t column = 0; column < 6; ++column)
                material_stress_derivative[row] +=
                    material_linearization.tangent[row][column] * eigenstrain_derivative[column];
        const SymmetricTensor3Values derivative = rotate_cartesian_tensor_values({material_stress_derivative[0],
                                                                                     material_stress_derivative[1],
                                                                                     material_stress_derivative[2],
                                                                                     material_stress_derivative[3],
                                                                                     material_stress_derivative[4],
                                                                                     material_stress_derivative[5]},
            kinematics.rotation);
        stress_linearization.committed_temperature_tangent =
            {derivative.xx, derivative.yy, derivative.zz, derivative.xy, derivative.yz, derivative.xz};
    }
    const ReducedHex8GeometryDerivatives current_derivatives =
        reduced_hex8_geometry_derivatives(geometry, reduced_hex8_displacement_values(state), current, 1.0);
    const ReducedHex8GeometryDerivatives midpoint_derivatives = reduced_hex8_geometry_derivatives(geometry,
        reduced_hex8_midpoint_displacement_values(state, old_state),
        midpoint,
        0.5);
    add_reduced_hex8_finite_jacobian(data,
        geometry,
        state,
        old_state,
        time_step,
        add_thermal_time_term,
        initial_shear_modulus,
        current,
        current_derivatives,
        midpoint,
        midpoint_derivatives,
        stress_linearization,
        *jacobian);
}

void assemble_c3d8rt_small_strain_system(const elements::C3d8Input& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalResidual& residual,
    Hex8LocalJacobian* jacobian) {
    if (data.strain_formulation != StrainFormulation::small)
        throw std::logic_error("C3D8RT small-strain integration received a non-small strain formulation");
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    if (jacobian != nullptr)
        jacobian->fill(0.0);

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
            active_gradient[component][direction] =
                adlite::Scalar::independent(gradient_values[component * 3 + direction],
                    component * 3 + direction,
                    point_width);
    const adlite::Scalar active_temperature =
        adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    const SymmetricTensor3 strain{active_gradient[0][0],
        active_gradient[1][1],
        active_gradient[2][2],
        0.5 * (active_gradient[0][1] + active_gradient[1][0]),
        0.5 * (active_gradient[1][2] + active_gradient[2][1]),
        0.5 * (active_gradient[0][2] + active_gradient[2][0])};
    const std::array<const adlite::Scalar*, 6> strain_components =
        {&strain.xx, &strain.yy, &strain.zz, &strain.xy, &strain.yz, &strain.xz};
    std::array<double, 6> fed_strain{};
    for (std::size_t component = 0; component < fed_strain.size(); ++component)
        fed_strain[component] = strain_components[component]->value();
    const MaterialFunctionContext context = material_context(data.time, point.position);
    const ActiveThermoelasticProperties initial_properties =
        data.material.active_properties(adlite::Scalar(data.initial_temperature),
            material_context(0.0, point.position));
    const double reference_shear_modulus = initial_properties.shear_modulus.value();
    if (!std::isfinite(reference_shear_modulus) || !(reference_shear_modulus > 0.0))
        throw std::invalid_argument("C3D8RT initial shear modulus must be finite and positive");
    const fuelsim::CartesianStressTangent tangent = fuelsim::evaluate_stress_tangent(data.material,
        fed_strain,
        temperature_value,
        time_step,
        committed_material,
        context);
    const SymmetricTensor3 stress = compose_cartesian_stress(tangent, strain, active_temperature, tangent.thermal);

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
            hourglass_thermal += geometry.hourglass_shape[node][mode] * geometry.thermal_hourglass_coefficients[mode]
                                 * temperature_hourglass_amplitudes[mode];
        const double thermal_operator = geometry.reference_volume * uniform_thermal + hourglass_thermal;
        residual[node] += conductivity.value() * thermal_operator
                          - geometry.reduced_body_source_measure * data.volumetric_heat_source / 8.0;

        const std::array<adlite::Scalar, 3> mechanical = {geometry.reference_volume
                                                              * (stress.xx * point.gradient[node][0]
                                                                  + stress.xy * point.gradient[node][1]
                                                                  + stress.xz * point.gradient[node][2]),
            geometry.reference_volume
                * (stress.xy * point.gradient[node][0] + stress.yy * point.gradient[node][1]
                    + stress.yz * point.gradient[node][2]),
            geometry.reference_volume
                * (stress.xz * point.gradient[node][0] + stress.yz * point.gradient[node][1]
                    + stress.zz * point.gradient[node][2])};
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t row = 8 * (component + 1) + node;
            residual[row] += mechanical[component].value();
            if (jacobian == nullptr)
                continue;
            mechanical[component].copy_derivatives(derivatives.data(), derivatives.size());
            for (std::size_t other = 0; other < hex8_node_count; ++other) {
                (*jacobian)[row * hex8_local_dof_count + other] +=
                    derivatives[temperature_index] * geometry.reduced_capacity_points[other].weighted_measure
                    / geometry.reference_volume;
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
                    hourglass_tangent += geometry.hourglass_shape[node][mode]
                                         * geometry.thermal_hourglass_coefficients[mode]
                                         * geometry.hourglass_shape[other][mode];
                (*jacobian)[node * hex8_local_dof_count + other] +=
                    conductivity.value() * (geometry.reference_volume * uniform_tangent + hourglass_tangent)
                    + conductivity_derivative * thermal_operator
                          * geometry.reduced_capacity_points[other].weighted_measure / geometry.reference_volume;
            }
    }

    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t mode = 0; mode < 4; ++mode) {
            double amplitude = 0.0;
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                amplitude += geometry.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
            const double stiffness = abaqus_total_stiffness_factor * reference_shear_modulus
                                     * geometry.mechanical_hourglass_metrics[component];
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
            const adlite::Scalar capacity = data.material.reference_heat_capacity(active_capacity_temperature,
                data.initial_temperature,
                capacity_context);
            const adlite::Scalar capacity_term = capacity_point.weighted_measure * capacity
                                                 * (active_capacity_temperature - (*committed_state)[node]) / time_step;
            residual[node] += capacity_term.value();
            if (jacobian != nullptr)
                (*jacobian)[node * hex8_local_dof_count + node] += capacity_term.derivative(0);
        }
}

Hex8LocalResidual compute_hex8_local(const elements::C3d8Input& data, Hex8LocalJacobian* jacobian) {
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    const auto* history = data.committed_history;
    const auto* committed_state = history || data.time_step > 0.0 ? &data.committed_state : nullptr;
    const double time_step = data.time_step;
    const bool include_thermal_time_term = data.include_thermal_time_term;
    if (committed_state != nullptr && (!std::isfinite(time_step) || time_step <= 0.0))
        throw std::invalid_argument("HEX8 time step must be finite and positive");
    if (history != nullptr && history->size() != 1)
        throw std::invalid_argument("C3D8RT material history has the wrong integration point count");
    Hex8LocalResidual residual{};
    if (data.strain_formulation == StrainFormulation::finite)
        assemble_c3d8rt_finite_strain_system(data,
            geometry,
            state,
            committed_state,
            history == nullptr ? nullptr : &history->front(),
            time_step,
            include_thermal_time_term,
            residual,
            jacobian);
    else
        assemble_c3d8rt_small_strain_system(data,
            geometry,
            state,
            committed_state,
            history == nullptr ? nullptr : &history->front(),
            time_step,
            include_thermal_time_term,
            residual,
            jacobian);
    return residual;
}

CartesianMaterialHistory compute_hex8_transient_update(const elements::C3d8Input& data) {
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    const auto& committed_state = data.committed_state;
    const auto& committed_material = *data.committed_history;
    const double time_step = data.time_step;
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("HEX8 transient update time step must be finite and positive");
    {
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
            const ActiveReducedHex8Geometry midpoint = active_reduced_hex8_geometry(geometry,
                reduced_hex8_midpoint_displacement(passive, committed_state),
                "midpoint");
            const cartesian_detail::KinematicsCore kinematics =
                reduced_hex8_finite_kinematics(geometry, midpoint, passive, committed_state);
            const adlite::Scalar temperature = reduced_hex8_temperature(current, passive);
            Hex8LocalAdValues passive_old{};
            ad_local_system::make_passive(committed_state.data(), committed_state.size(), passive_old.data());
            (void)active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive_old), "committed");
            const double old_temperature = reduced_hex8_temperature(current, passive_old).value();
            response = data.material.incremental_response(kinematics.strain_increment,
                kinematics.rotation,
                temperature,
                old_temperature,
                time_step,
                committed_material.front(),
                context);
        } else {
            const adlite::Scalar temperature = reduced_hex8_temperature(geometry, passive);
            const C3d8Kinematics kinematics =
                evaluate_cartesian_incremental_kinematics(point, passive, committed_state, StrainFormulation::small);
            response = data.material.response(kinematics.strain_increment,
                temperature,
                time_step,
                committed_material.front(),
                context);
        }
        return CartesianMaterialHistory{response.trial_state};
    }
}

std::array<SymmetricTensor3Values, 8> compute_hex8_stress(const elements::C3d8Input& data) {
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    {
        Hex8LocalAdValues passive{};
        ad_local_system::make_passive(state.data(), state.size(), passive.data());
        SymmetricTensor3 stress;
        if (data.strain_formulation == StrainFormulation::finite) {
            const Hex8LocalValues undeformed{};
            const ActiveReducedHex8Geometry current =
                active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive), "current");
            const ActiveReducedHex8Geometry midpoint = active_reduced_hex8_geometry(geometry,
                reduced_hex8_midpoint_displacement(passive, undeformed),
                "midpoint");
            const cartesian_detail::KinematicsCore kinematics =
                reduced_hex8_finite_kinematics(geometry, midpoint, passive, undeformed);
            const adlite::Scalar temperature = reduced_hex8_temperature(current, passive);
            stress = rotate_cartesian_tensor(data.material.stress(kinematics.strain_increment,
                                                 temperature,
                                                 material_context(data.time, geometry.reduced_point.position)),
                kinematics.rotation);
        } else {
            const adlite::Scalar temperature = reduced_hex8_temperature(geometry, passive);
            const C3d8Kinematics kinematics = evaluate_cartesian_incremental_kinematics(geometry.reduced_point,
                passive,
                Hex8LocalValues{},
                StrainFormulation::small);
            stress = data.material.stress(kinematics.strain_increment,
                temperature,
                material_context(data.time, geometry.reduced_point.position));
        }
        const SymmetricTensor3Values values{stress.xx.value(),
            stress.yy.value(),
            stress.zz.value(),
            stress.xy.value(),
            stress.yz.value(),
            stress.xz.value()};
        std::array<SymmetricTensor3Values, 8> result{};
        result.fill(values);
        return result;
    }
}
} // namespace

double elements::c3d8rt_hourglass_energy(const C3d8Input& data) {
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    const ActiveThermoelasticProperties initial_properties =
        data.material.active_properties(adlite::Scalar(data.initial_temperature),
            material_context(0.0, geometry.reduced_point.position));
    const double initial_shear_modulus = initial_properties.shear_modulus.value();
    if (!std::isfinite(initial_shear_modulus) || !(initial_shear_modulus > 0.0))
        throw std::invalid_argument("C3D8RT initial shear modulus must be finite and positive");

    double energy = 0.0;
    if (data.strain_formulation == StrainFormulation::small) {
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t mode = 0; mode < 4; ++mode) {
                double amplitude = 0.0;
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    amplitude += geometry.hourglass_shape[node][mode] * state[8 * (component + 1) + node];
                const double stiffness = abaqus_total_stiffness_factor * initial_shear_modulus
                                         * geometry.mechanical_hourglass_metrics[component];
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
            const double stiffness = abaqus_total_stiffness_factor * initial_shear_modulus
                                     * geometry.mechanical_hourglass_metrics[material_direction];
            energy += 0.5 * stiffness * transported_amplitude * transported_amplitude;
        }
    }
    return energy;
}
} // namespace fuelsim

namespace fuelsim::elements {
C3d8Result evaluate_c3d8rt(const C3d8Input& input, ElementRequest request) {
    C3d8Result result;
    if (request.residual || request.jacobian)
        result.residual = compute_hex8_local(input, request.jacobian ? &result.jacobian : nullptr);
    if (request.history && input.committed_history) {
        result.history = compute_hex8_transient_update(input);
        c3d8_detail::set_history_geometry(input, true, result);
    }
    if (request.stress)
        result.stress = compute_hex8_stress(input);
    if (request.residual || request.jacobian)
        add_hex8_body_acceleration(input, result.residual);
    return result;
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Hex8Geometry make_c3d8rt_geometry(const Hex8Coordinates& coordinates) {
    return c3d8_detail::make_hex8_geometry(coordinates);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
void validate_c3d8rt_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    c3d8_detail::validate_cartesian_deformation(point, state);
}
} // namespace fuelsim::elements
