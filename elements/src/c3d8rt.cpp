#include "c3d8rt.hpp"
#include "detail/hex8_geometry.hpp"

namespace fuelsim {
namespace {
using namespace element_detail;
constexpr std::array<std::array<double, 4>, hex8_node_count> finite_reduced_hex8_raw_hourglass = {
    {{{1.0, -1.0, 1.0, -1.0}},
        {{-1.0, -1.0, -1.0, 1.0}},
        {{1.0, 1.0, -1.0, -1.0}},
        {{-1.0, 1.0, 1.0, 1.0}},
        {{1.0, 1.0, -1.0, 1.0}},
        {{-1.0, 1.0, 1.0, -1.0}},
        {{1.0, -1.0, 1.0, 1.0}},
        {{-1.0, -1.0, -1.0, -1.0}}}};

struct ActiveReducedHex8Geometry final {
    adlite::Scalar volume{0.0}, center_measure{0.0};
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> average_gradient{};
    std::array<adlite::Scalar, hex8_node_count> shape_measures{};
    std::array<std::array<adlite::Scalar, 4>, hex8_node_count> hourglass_shape{};
    std::array<adlite::Scalar, 4> thermal_hourglass_coefficients{};
};

struct ReducedHex8GeometryValues final {
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

struct ReducedFiniteMaterialLinearization final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 6>, 6> tangent{};
    std::array<double, 6> thermal{};
};

struct ReducedFiniteStressLinearization final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 10>, 6> tangent{};
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
cartesian_detail::Matrix3 reduced_hex8_central_gradient_values(const ReducedHex8GeometryValues& midpoint,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state);
ReducedFiniteKinematicsValues reduced_hex8_finite_kinematics_values(const cartesian_detail::Matrix3& central_gradient);
cartesian_detail::ActiveMatrix3 reduced_hex8_central_gradient(const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state);
cartesian_detail::KinematicsCore reduced_hex8_finite_kinematics(const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state);
ReducedFiniteMaterialLinearization reduced_finite_material_linearization(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3Values& strain_increment,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context);
ReducedFiniteStressLinearization reduced_finite_stress_linearization(const cartesian_detail::Matrix3& passive_gradient,
    double temperature,
    const ReducedFiniteMaterialLinearization& material_linearization);
adlite::Scalar reduced_hex8_temperature(const ActiveReducedHex8Geometry& geometry, const Hex8LocalAdValues& state);
double reduced_hex8_temperature_value(const ReducedHex8GeometryValues& geometry, const Hex8LocalValues& state);
SymmetricTensor3Values reduced_finite_stress_values(const IsotropicThermoelasticMaterial& material,
    const ReducedFiniteKinematicsValues& kinematics,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context);
Hex8LocalResidual reduced_hex8_finite_residual_values(const CartesianThermoelasticData& data,
    const Hex8Geometry& reference,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    double time_step,
    bool include_thermal_time_term,
    double initial_shear_modulus,
    const ReducedHex8GeometryValues& current,
    const SymmetricTensor3Values& stress);
void add_reduced_hex8_finite_jacobian(const CartesianThermoelasticData& data,
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
void assemble_c3d8rt_finite_strain_system(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian);
void assemble_c3d8rt_small_strain_system(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian);
Hex8LocalResidual compute_hex8_local(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term);
Hex8LocalResidual compute_hex8_thermoelastic(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term);
Hex8LocalResidual compute_hex8_transient(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term);
CartesianMaterialHistory compute_hex8_transient_update(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step);
std::array<SymmetricTensor3Values, 8>
compute_hex8_stress(const CartesianThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);

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

    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<double, 3> projected_coordinate{};
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

    cartesian_detail::Matrix3 inverse_effective_mapping{};
    for (std::size_t natural = 0; natural < 3; ++natural)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                inverse_effective_mapping[natural][physical] +=
                    hex8_signs[node][natural] * result.average_gradient[node][physical];
    const double inverse_effective_determinant = cartesian_detail::determinant(inverse_effective_mapping);
    if (!std::isfinite(inverse_effective_determinant) || inverse_effective_determinant == 0.0)
        throw std::domain_error(std::string("C3D8RT ") + configuration_name + " effective mapping must be nonsingular");
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
        metric[2][2]
        - (metric[1][1] * metric[0][2] * metric[0][2] - 2.0 * metric[0][1] * metric[0][2] * metric[1][2]
              + metric[0][0] * metric[1][2] * metric[1][2])
              / leading_determinant;
    if (!std::isfinite(first_pivot) || !std::isfinite(second_pivot) || !std::isfinite(third_pivot)
        || !(first_pivot > 0.0) || !(second_pivot > 0.0) || !(third_pivot > 0.0))
        throw std::domain_error(
            std::string("C3D8RT ") + configuration_name + " effective metric must be positive definite");
    const double thermal_scale = result.volume / 192.0;
    const double inverse_x = 1.0 / first_pivot, inverse_y = 1.0 / second_pivot, inverse_z = 1.0 / third_pivot;
    result.thermal_hourglass_coefficients = {thermal_scale * (inverse_x + inverse_y),
        thermal_scale * (inverse_x + inverse_z),
        thermal_scale * (inverse_x + inverse_z),
        thermal_scale * (inverse_x + inverse_y + inverse_z) / 3.0};
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
                    current_coordinates[node][component] * finite_reduced_hex8_raw_hourglass[node][mode];
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t mode = 0; mode < 4; ++mode)
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t active_node = 0; active_node < hex8_node_count; ++active_node) {
                    const std::size_t column = 8 * component + active_node;
                    double derivative = displacement_derivative_scale * values.average_gradient[node][component]
                                        * finite_reduced_hex8_raw_hourglass[active_node][mode];
                    for (std::size_t direction = 0; direction < 3; ++direction)
                        derivative +=
                            result.average_gradient[node][direction][column] * projected_coordinates[mode][direction];
                    result.hourglass_shape[node][mode][column] = -derivative;
                }

    cartesian_detail::Matrix3 inverse_effective_mapping{};
    for (std::size_t natural = 0; natural < 3; ++natural)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                inverse_effective_mapping[natural][physical] +=
                    hex8_signs[node][natural] * values.average_gradient[node][physical];
    const cartesian_detail::Matrix3 effective_mapping =
        cartesian_detail::inverse(inverse_effective_mapping, cartesian_detail::determinant(inverse_effective_mapping));
    cartesian_detail::Matrix3 metric{};
    for (std::size_t first = 0; first < 3; ++first)
        for (std::size_t second = 0; second < 3; ++second)
            for (std::size_t physical = 0; physical < 3; ++physical)
                metric[first][second] += effective_mapping[physical][first] * effective_mapping[physical][second];
    const double first_pivot = metric[0][0];
    const double second_pivot = metric[1][1] - metric[0][1] * metric[0][1] / first_pivot;
    const double leading = metric[0][0] * metric[1][1] - metric[0][1] * metric[0][1];
    const double numerator = metric[1][1] * metric[0][2] * metric[0][2]
                             - 2.0 * metric[0][1] * metric[0][2] * metric[1][2]
                             + metric[0][0] * metric[1][2] * metric[1][2];
    const double third_pivot = metric[2][2] - numerator / leading;
    const double inverse_x = 1.0 / first_pivot, inverse_y = 1.0 / second_pivot, inverse_z = 1.0 / third_pivot;
    const double thermal_scale = values.volume / 192.0;
    const std::array<double, 4> sums = {inverse_x + inverse_y,
        inverse_x + inverse_z,
        inverse_x + inverse_z,
        (inverse_x + inverse_y + inverse_z) / 3.0};
    for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
        cartesian_detail::Matrix3 mapping_derivative{}, effective_derivative{};
        for (std::size_t natural = 0; natural < 3; ++natural)
            for (std::size_t physical = 0; physical < 3; ++physical)
                for (std::size_t node = 0; node < hex8_node_count; ++node)
                    mapping_derivative[natural][physical] +=
                        hex8_signs[node][natural] * result.average_gradient[node][physical][column];
        for (std::size_t first = 0; first < 3; ++first)
            for (std::size_t second = 0; second < 3; ++second)
                for (std::size_t i = 0; i < 3; ++i)
                    for (std::size_t j = 0; j < 3; ++j)
                        effective_derivative[first][second] -=
                            effective_mapping[first][i] * mapping_derivative[i][j] * effective_mapping[j][second];
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
        const double thermal_scale_derivative = result.volume[column] / 192.0;
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
        metric[2][2]
        - (metric[1][1] * metric[0][2] * metric[0][2] - 2.0 * metric[0][1] * metric[0][2] * metric[1][2]
              + metric[0][0] * metric[1][2] * metric[1][2])
              / leading_determinant;
    if (!std::isfinite(first_pivot.value()) || !std::isfinite(second_pivot.value())
        || !std::isfinite(third_pivot.value()) || !(first_pivot.value() > 0.0) || !(second_pivot.value() > 0.0)
        || !(third_pivot.value() > 0.0))
        throw std::domain_error(
            std::string("C3D8RT ") + configuration_name + " effective metric must be positive definite");
    const adlite::Scalar thermal_scale = result.volume / 192.0;
    const adlite::Scalar inverse_x = 1.0 / first_pivot, inverse_y = 1.0 / second_pivot, inverse_z = 1.0 / third_pivot;
    result.thermal_hourglass_coefficients = {thermal_scale * (inverse_x + inverse_y),
        thermal_scale * (inverse_x + inverse_z),
        thermal_scale * (inverse_x + inverse_z),
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

cartesian_detail::Matrix3 reduced_hex8_central_gradient_values(const ReducedHex8GeometryValues& midpoint,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state) {
    cartesian_detail::Matrix3 central_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                central_gradient[component][direction] +=
                    (state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node])
                    * midpoint.average_gradient[node][direction];
    return central_gradient;
}

ReducedFiniteKinematicsValues reduced_hex8_finite_kinematics_values(const cartesian_detail::Matrix3& central_gradient) {
    ReducedFiniteKinematicsValues result;
    cartesian_detail::Matrix3 spatial_strain{}, rotation_numerator{}, rotation_denominator{};
    for (std::size_t i = 0; i < 3; ++i) {
        rotation_numerator[i][i] = 1.0;
        rotation_denominator[i][i] = 1.0;
        for (std::size_t j = 0; j < 3; ++j) {
            spatial_strain[i][j] = 0.5 * (central_gradient[i][j] + central_gradient[j][i]);
            const double half_spin = 0.25 * (central_gradient[i][j] - central_gradient[j][i]);
            rotation_numerator[i][j] += half_spin;
            rotation_denominator[i][j] -= half_spin;
        }
    }
    const double denominator_determinant = cartesian_detail::determinant(rotation_denominator);
    if (!std::isfinite(denominator_determinant) || denominator_determinant == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian rotation denominator is singular");
    result.rotation =
        multiply_matrices(rotation_numerator, cartesian_detail::inverse(rotation_denominator, denominator_determinant));
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
    return result;
}

cartesian_detail::ActiveMatrix3 reduced_hex8_central_gradient(const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state) {
    cartesian_detail::ActiveMatrix3 central_gradient{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                central_gradient[component][direction] +=
                    (state[8 * (component + 1) + node] - committed_state[8 * (component + 1) + node])
                    * midpoint.average_gradient[node][direction];
    return central_gradient;
}

cartesian_detail::KinematicsCore reduced_hex8_finite_kinematics(const ActiveReducedHex8Geometry& midpoint,
    const Hex8LocalAdValues& state,
    const Hex8LocalValues& committed_state) {
    return cartesian_detail::evaluate_hughes_winget_increment(
        reduced_hex8_central_gradient(midpoint, state, committed_state));
}

ReducedFiniteMaterialLinearization reduced_finite_material_linearization(const IsotropicThermoelasticMaterial& material,
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
    ReducedFiniteMaterialLinearization result;
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
    const ReducedFiniteMaterialLinearization& material_linearization) {
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
    const std::array<adlite::Scalar, 7> material_inputs = {kinematics.strain_increment.xx,
        kinematics.strain_increment.yy,
        kinematics.strain_increment.zz,
        kinematics.strain_increment.xy,
        kinematics.strain_increment.yz,
        kinematics.strain_increment.xz,
        active[9]};
    const std::array<double, 6> stress_values = {material_linearization.stress.xx,
        material_linearization.stress.yy,
        material_linearization.stress.zz,
        material_linearization.stress.xy,
        material_linearization.stress.yz,
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
    const SymmetricTensor3 rotated = rotate_cartesian_tensor({material_stress[0],
                                                                 material_stress[1],
                                                                 material_stress[2],
                                                                 material_stress[3],
                                                                 material_stress[4],
                                                                 material_stress[5]},
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

Hex8LocalResidual reduced_hex8_finite_residual_values(const CartesianThermoelasticData& data,
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
        const double gradient_x = current.average_gradient[node][0];
        const double gradient_y = current.average_gradient[node][1];
        const double gradient_z = current.average_gradient[node][2];
        residual[8 + node] +=
            current.volume * (stress.xx * gradient_x + stress.xy * gradient_y + stress.xz * gradient_z);
        residual[16 + node] +=
            current.volume * (stress.xy * gradient_x + stress.yy * gradient_y + stress.yz * gradient_z);
        residual[24 + node] +=
            current.volume * (stress.xz * gradient_x + stress.yz * gradient_y + stress.zz * gradient_z);
    }

    cartesian_detail::Matrix3 average_deformation{};
    for (std::size_t component = 0; component < 3; ++component) {
        average_deformation[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                average_deformation[component][direction] +=
                    state[8 * (component + 1) + node] * reference.average_shape_gradient[node][direction];
    }
    constexpr double abaqus_total_stiffness_factor = 0.005;
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
            const double capacity = data.material
                                        .heat_capacity(adlite::Scalar(state[node]),
                                            material_context(data.time, reference.capacity_points[node].position))
                                        .value();
            residual[node] +=
                current.shape_measures[node] * capacity * (state[node] - committed_state[node]) / time_step;
        }
    return residual;
}

void add_reduced_hex8_finite_jacobian(const CartesianThermoelasticData& data,
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

    const double temperature = reduced_hex8_temperature_value(current, state);
    ReducedDisplacementDerivatives temperature_displacement_derivatives{};
    std::array<double, hex8_node_count> temperature_temperature_derivatives{};
    for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
        double numerator_derivative = 0.0;
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            numerator_derivative += current_derivatives.shape_measures[node][column] * state[node];
        temperature_displacement_derivatives[column] =
            (numerator_derivative - temperature * current_derivatives.volume[column]) / current.volume;
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
                stress_linearization.tangent[stress_component][9] * temperature_displacement_derivatives[column];
            for (std::size_t gradient = 0; gradient < 9; ++gradient)
                derivative += stress_linearization.tangent[stress_component][gradient]
                              * central_gradient_derivatives[gradient][column];
            stress_displacement_derivatives[stress_component][column] = derivative;
        }
        for (std::size_t node = 0; node < hex8_node_count; ++node)
            stress_temperature_derivatives[stress_component][node] =
                stress_linearization.tangent[stress_component][9] * temperature_temperature_derivatives[node];
    }

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

        const std::array<double, 3> gradient = {current.average_gradient[node][0],
            current.average_gradient[node][1],
            current.average_gradient[node][2]};
        const std::array<std::array<double, 3>, 3> stress_matrix = {
            {{stress_values[0], stress_values[3], stress_values[5]},
                {stress_values[3], stress_values[1], stress_values[4]},
                {stress_values[5], stress_values[4], stress_values[2]}}};
        constexpr std::array<std::array<std::size_t, 3>, 3> stress_indices = {{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
        for (std::size_t row_component = 0; row_component < 3; ++row_component) {
            const std::size_t row = 8 * (row_component + 1) + node;
            double force_density = 0.0;
            for (std::size_t direction = 0; direction < 3; ++direction)
                force_density += stress_matrix[row_component][direction] * gradient[direction];
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column) {
                double force_density_derivative = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction) {
                    const std::size_t stress_index = stress_indices[row_component][direction];
                    force_density_derivative +=
                        stress_displacement_derivatives[stress_index][column] * gradient[direction]
                        + stress_matrix[row_component][direction]
                              * current_derivatives.average_gradient[node][direction][column];
                }
                jacobian[row * hex8_local_dof_count + 8 + column] +=
                    current_derivatives.volume[column] * force_density + current.volume * force_density_derivative;
            }
            for (std::size_t temperature_node = 0; temperature_node < hex8_node_count; ++temperature_node) {
                double derivative = 0.0;
                for (std::size_t direction = 0; direction < 3; ++direction)
                    derivative +=
                        stress_temperature_derivatives[stress_indices[row_component][direction]][temperature_node]
                        * gradient[direction];
                jacobian[row * hex8_local_dof_count + temperature_node] += current.volume * derivative;
            }
        }
    }

    cartesian_detail::Matrix3 average_deformation{};
    for (std::size_t component = 0; component < 3; ++component) {
        average_deformation[component][component] = 1.0;
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                average_deformation[component][direction] +=
                    state[8 * (component + 1) + node] * reference.average_shape_gradient[node][direction];
    }
    constexpr double abaqus_total_stiffness_factor = 0.005;
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
            const adlite::Scalar capacity = data.material.heat_capacity(nodal_temperature,
                material_context(data.time, reference.capacity_points[node].position));
            const double rate = capacity.value() * (state[node] - committed_state[node]) / time_step;
            for (std::size_t column = 0; column < reduced_displacement_dof_count; ++column)
                jacobian[node * hex8_local_dof_count + 8 + column] +=
                    current_derivatives.shape_measures[node][column] * rate;
            const double capacity_derivative = capacity.is_active() ? capacity.derivative(0) : 0.0;
            jacobian[node * hex8_local_dof_count + node] +=
                current.shape_measures[node]
                * (capacity.value() + capacity_derivative * (state[node] - committed_state[node])) / time_step;
        }
}

void assemble_c3d8rt_finite_strain_system(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
    Hex8LocalJacobian* jacobian) {
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    const Hex8LocalValues undeformed{};
    const Hex8LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const ReducedHex8GeometryValues current =
        reduced_hex8_geometry_values(geometry, reduced_hex8_displacement_values(state), "current");
    const ReducedHex8GeometryValues midpoint =
        reduced_hex8_geometry_values(geometry, reduced_hex8_midpoint_displacement_values(state, old_state), "midpoint");
    const cartesian_detail::Matrix3 central_gradient = reduced_hex8_central_gradient_values(midpoint, state, old_state);
    const ReducedFiniteKinematicsValues kinematics = reduced_hex8_finite_kinematics_values(central_gradient);
    const double temperature = reduced_hex8_temperature_value(current, state);
    double old_temperature = data.initial_temperature;
    if (committed_material != nullptr) {
        const ReducedHex8GeometryValues old_geometry =
            reduced_hex8_geometry_values(geometry, reduced_hex8_displacement_values(old_state), "committed");
        old_temperature = reduced_hex8_temperature_value(old_geometry, old_state);
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
    const Hex8LocalResidual passive_residual = reduced_hex8_finite_residual_values(data,
        geometry,
        state,
        old_state,
        time_step,
        add_thermal_time_term,
        initial_shear_modulus,
        current,
        stress);
    for (std::size_t row = 0; row < hex8_local_dof_count; ++row)
        residual[row] = passive_residual[row];
    if (jacobian == nullptr)
        return;

    jacobian->fill(0.0);
    const ReducedFiniteMaterialLinearization material_linearization =
        reduced_finite_material_linearization(data.material,
            kinematics.strain_increment,
            temperature,
            old_temperature,
            time_step,
            committed_material,
            context);
    const ReducedFiniteStressLinearization stress_linearization =
        reduced_finite_stress_linearization(central_gradient, temperature, material_linearization);
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

void assemble_c3d8rt_small_strain_system(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialPointState* committed_material,
    double time_step,
    bool include_thermal_time_term,
    Hex8LocalAdValues& residual,
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
    const cartesian_detail::CartesianStressTangent tangent = cartesian_detail::evaluate_stress_tangent(data.material,
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
    std::array<adlite::Scalar, 6> composed_stress{};
    for (std::size_t component = 0; component < 6; ++component) {
        for (std::size_t column = 0; column < 6; ++column)
            partials[column] = tangent.tangent[component][column];
        partials[6] = tangent.thermal[component];
        composed_stress[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    const SymmetricTensor3 stress{composed_stress[0],
        composed_stress[1],
        composed_stress[2],
        composed_stress[3],
        composed_stress[4],
        composed_stress[5]};

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

    constexpr double abaqus_total_stiffness_factor = 0.005;
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
            const adlite::Scalar capacity = data.material.heat_capacity(active_capacity_temperature, capacity_context);
            const adlite::Scalar capacity_term = capacity_point.weighted_measure * capacity
                                                 * (active_capacity_temperature - (*committed_state)[node]) / time_step;
            residual[node] += capacity_term.value();
            if (jacobian != nullptr)
                (*jacobian)[node * hex8_local_dof_count + node] += capacity_term.derivative(0);
        }
}

Hex8LocalResidual compute_hex8_local(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state,
    const CartesianMaterialHistory* history,
    double time_step,
    Hex8LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    if (committed_state != nullptr && (!std::isfinite(time_step) || time_step <= 0.0))
        throw std::invalid_argument("HEX8 time step must be finite and positive");
    if (history != nullptr && history->size() != 1)
        throw std::invalid_argument("C3D8RT material history has the wrong integration point count");
    Hex8LocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
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
    Hex8LocalResidual result{};
    ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    return result;
}

Hex8LocalResidual compute_hex8_thermoelastic(const CartesianThermoelasticData& data,
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

Hex8LocalResidual compute_hex8_transient(const CartesianThermoelasticData& data,
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

CartesianMaterialHistory compute_hex8_transient_update(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step) {
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
                reduced_hex8_finite_kinematics(midpoint, passive, committed_state);
            const adlite::Scalar temperature = reduced_hex8_temperature(current, passive);
            Hex8LocalAdValues passive_old{};
            ad_local_system::make_passive(committed_state.data(), committed_state.size(), passive_old.data());
            const ActiveReducedHex8Geometry old_geometry =
                active_reduced_hex8_geometry(geometry, reduced_hex8_displacement(passive_old), "committed");
            const double old_temperature = reduced_hex8_temperature(old_geometry, passive_old).value();
            response = data.material.incremental_response(kinematics.strain_increment,
                kinematics.rotation,
                temperature,
                old_temperature,
                time_step,
                committed_material.front(),
                context);
        } else {
            const adlite::Scalar temperature = reduced_hex8_temperature(geometry, passive);
            const CartesianKinematics kinematics =
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

std::array<SymmetricTensor3Values, 8> compute_hex8_stress(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state) {
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
                reduced_hex8_finite_kinematics(midpoint, passive, undeformed);
            const adlite::Scalar temperature = reduced_hex8_temperature(current, passive);
            stress = rotate_cartesian_tensor(data.material.stress(kinematics.strain_increment,
                                                 temperature,
                                                 material_context(data.time, geometry.reduced_point.position)),
                kinematics.rotation);
        } else {
            const adlite::Scalar temperature = reduced_hex8_temperature(geometry, passive);
            const CartesianKinematics kinematics = evaluate_cartesian_incremental_kinematics(geometry.reduced_point,
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

double compute_c3d8_mechanical_hourglass_energy(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state) {
    if (data.hex8_element_formulation != Hex8ElementFormulation::c3d8rt)
        return 0.0;
    if (!std::isfinite(data.initial_temperature) || !(data.initial_temperature > 0.0))
        throw std::invalid_argument("C3D8RT requires a finite positive initial temperature");
    const ActiveThermoelasticProperties initial_properties =
        data.material.active_properties(adlite::Scalar(data.initial_temperature),
            material_context(0.0, geometry.reduced_point.position));
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

#include "c3d8rt.hpp"
#include <stdexcept>

namespace fuelsim::elements {
C3d8Result evaluate_c3d8rt(const C3d8Input& input, ElementRequest request) {
    const CartesianThermoelasticData data{input.material,
        input.volumetric_heat_source,
        input.time,
        input.strain_formulation,
        Hex8ElementFormulation::c3d8rt,
        input.initial_temperature};
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
    if (request.stress)
        result.stress = compute_hex8_stress(data, input.geometry, input.state);
    return result;
}
} // namespace fuelsim::elements
