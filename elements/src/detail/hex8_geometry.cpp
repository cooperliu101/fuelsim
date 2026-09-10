#include "hex8_geometry.hpp"

namespace fuelsim::element_detail {
cartesian_detail::ActiveMatrix3 displacement_gradient(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& state) {
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

CartesianKinematics evaluate_cartesian_kinematics_from_gradient(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state,
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

cartesian_detail::Matrix3 multiply_matrices(const cartesian_detail::Matrix3& first,
    const cartesian_detail::Matrix3& second) {
    cartesian_detail::Matrix3 result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                result[i][j] += first[i][k] * second[k][j];
    return result;
}
} // namespace fuelsim::element_detail

namespace fuelsim {
using namespace element_detail;

void validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    const double value = cartesian_detail::determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX8 deformation must preserve a positive Jacobian");
}

CartesianKinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    return evaluate_cartesian_kinematics_from_gradient(point,
        displacement_gradient(point, current_state),
        committed_state,
        strain_formulation);
}

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
        geometry.capacity_points[node] = {coordinates[node],
            geometry.points[hex8_node_to_gauss[node]].weighted_measure};
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

    constexpr std::array<std::array<double, 4>, hex8_node_count> raw_hourglass = {{{{1.0, -1.0, 1.0, -1.0}},
        {{-1.0, -1.0, -1.0, 1.0}},
        {{1.0, 1.0, -1.0, -1.0}},
        {{-1.0, 1.0, 1.0, 1.0}},
        {{1.0, 1.0, -1.0, 1.0}},
        {{-1.0, 1.0, 1.0, -1.0}},
        {{1.0, -1.0, 1.0, 1.0}},
        {{-1.0, -1.0, -1.0, -1.0}}}};
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
        metric[2][2]
        - (metric[1][1] * metric[0][2] * metric[0][2] - 2.0 * metric[0][1] * metric[0][2] * metric[1][2]
              + metric[0][0] * metric[1][2] * metric[1][2])
              / leading_determinant;
    if (!std::isfinite(first_pivot) || !std::isfinite(second_pivot) || !std::isfinite(third_pivot)
        || !(first_pivot > 0.0) || !(second_pivot > 0.0) || !(third_pivot > 0.0))
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
} // namespace fuelsim
