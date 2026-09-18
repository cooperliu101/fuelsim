#include "c3d_common.hpp"
#include "ad_local_system.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim::cartesian_detail {
namespace {
template <typename Scalar> Scalar determinant_impl(const std::array<std::array<Scalar, 3>, 3>& matrix) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>,
        "Local matrix operations only support double and adlite::Scalar");
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
           - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
           + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

template <typename Scalar>
std::array<std::array<Scalar, 3>, 3> inverse_impl(const std::array<std::array<Scalar, 3>, 3>& matrix,
    const Scalar& determinant_value) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>,
        "Local matrix operations only support double and adlite::Scalar");
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

template <typename Scalar>
std::array<std::array<Scalar, 3>, 3> multiply_impl(const std::array<std::array<Scalar, 3>, 3>& first,
    const Matrix3& second) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>,
        "Local matrix operations only support double and adlite::Scalar");
    std::array<std::array<Scalar, 3>, 3> result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                result[i][j] += first[i][k] * second[k][j];
    return result;
}

template <typename Scalar>
std::array<std::array<Scalar, 3>, 3> central_increment_gradient_impl(
    const std::array<std::array<Scalar, 3>, 3>& current,
    const Matrix3& committed,
    std::array<std::array<Scalar, 3>, 3>* midpoint_inverse) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>);
    std::array<std::array<Scalar, 3>, 3> sum{}, difference{}, gradient{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            sum[i][j] = current[i][j] + committed[i][j];
            difference[i][j] = current[i][j] - committed[i][j];
        }
    const Scalar det = determinant(sum);
    using std::isfinite;
    if (!isfinite(det) || det == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian increment has singular delta-F plus identity");
    const auto plus_inverse = inverse(sum, det);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            if (midpoint_inverse)
                (*midpoint_inverse)[i][j] = 2.0 * plus_inverse[i][j];
            for (std::size_t k = 0; k < 3; ++k)
                gradient[i][j] += 2.0 * difference[i][k] * plus_inverse[k][j];
        }
    return gradient;
}
} // namespace

double determinant(const Matrix3& matrix) {
    return determinant_impl(matrix);
}

adlite::Scalar determinant(const ActiveMatrix3& matrix) {
    return determinant_impl(matrix);
}

Matrix3 inverse(const Matrix3& matrix, double determinant_value) {
    return inverse_impl(matrix, determinant_value);
}

ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value) {
    return inverse_impl(matrix, determinant_value);
}

ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second) {
    return multiply_impl(first, second);
}

Matrix3 multiply(const Matrix3& first, const Matrix3& second) {
    return multiply_impl(first, second);
}

Matrix3 central_increment_gradient(const Matrix3& current, const Matrix3& committed, Matrix3* midpoint_inverse) {
    return central_increment_gradient_impl(current, committed, midpoint_inverse);
}

ActiveMatrix3
central_increment_gradient(const ActiveMatrix3& current, const Matrix3& committed, ActiveMatrix3* midpoint_inverse) {
    return central_increment_gradient_impl(current, committed, midpoint_inverse);
}
} // namespace fuelsim::cartesian_detail

namespace fuelsim::cartesian_detail {
namespace {
ActiveMatrix3 identity_active_matrix() {
    ActiveMatrix3 result{};
    for (std::size_t index = 0; index < 3; ++index)
        result[index][index] = 1.0;
    return result;
}
} // namespace

KinematicsCore evaluate_hughes_winget_increment(const ActiveMatrix3& central_displacement_gradient) {
    KinematicsCore result{};
    const ActiveMatrix3& hughes_winget = central_displacement_gradient;
    ActiveMatrix3 rotation{};
    std::array<adlite::Scalar, 6> strain{};
    hughes_winget_rotation(hughes_winget, rotation, strain);
    result.strain_increment = {strain[0], strain[1], strain[2], strain[3], strain[4], strain[5]};
    result.rotation = {rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2]};
    return result;
}

KinematicsCore evaluate_kinematics(const ActiveMatrix3& gradient,
    const Matrix3& committed_deformation,
    StrainFormulation strain_formulation) {
    KinematicsCore result{};
    result.current_inverse = identity_active_matrix();
    if (strain_formulation == StrainFormulation::small) {
        result.strain_increment = {gradient[0][0],
            gradient[1][1],
            gradient[2][2],
            0.5 * (gradient[0][1] + gradient[1][0]),
            0.5 * (gradient[1][2] + gradient[2][1]),
            0.5 * (gradient[0][2] + gradient[2][0])};
        return result;
    }
    ActiveMatrix3 current = gradient;
    for (std::size_t direction = 0; direction < 3; ++direction)
        current[direction][direction] += 1.0;
    result.current_determinant = determinant(current);
    if (!std::isfinite(result.current_determinant.value()) || !(result.current_determinant.value() > 0.0))
        throw std::domain_error("Finite-strain Cartesian deformation must preserve a positive Jacobian");
    result.current_inverse = inverse(current, result.current_determinant);
    const double old_determinant = determinant(committed_deformation);
    if (!std::isfinite(old_determinant) || !(old_determinant > 0.0))
        throw std::domain_error("Committed finite-strain Cartesian state requires a positive Jacobian");
    const double incremental_determinant = result.current_determinant.value() / old_determinant;
    if (!std::isfinite(incremental_determinant) || !(incremental_determinant > 0.0))
        throw std::domain_error("Incremental finite-strain Cartesian state requires a positive Jacobian");
    const ActiveMatrix3 hughes_winget =
        central_increment_gradient(current, committed_deformation, &result.midpoint_inverse);
    ActiveMatrix3 rotation{};
    std::array<adlite::Scalar, 6> strain{};
    hughes_winget_rotation(hughes_winget, rotation, strain);
    result.strain_increment = {strain[0], strain[1], strain[2], strain[3], strain[4], strain[5]};
    result.rotation = {rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2]};
    return result;
}

} // namespace fuelsim::cartesian_detail

namespace fuelsim {
using namespace c3d8_detail;
using namespace cartesian_detail;

namespace {
constexpr double gauss = 0.577350269189625764509148780501957456;
} // namespace

Hex8Geometry c3d8_detail::make_hex8_geometry(const Hex8Coordinates& coordinates) {
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
                hex8_shape_values(xi, eta, zeta, shape, derivative);
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
            geometry.points[hex8_node_gauss_permutation[node]].weighted_measure};
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

    std::array<std::array<double, 3>, 8> nodal_coordinates{};
    for (std::size_t node = 0; node < 8; ++node)
        nodal_coordinates[node] = {coordinates[node].x, coordinates[node].y, coordinates[node].z};
    geometry.hourglass_shape = hex8_hourglass_shape(nodal_coordinates, geometry.average_shape_gradient);

    const cartesian_detail::Matrix3 center_inverse =
        cartesian_detail::inverse(center_jacobian, geometry.reduced_body_source_measure / 8.0);
    std::array<std::array<double, 3>, hex8_node_count> center_gradient{};
    for (std::size_t node = 0; node < hex8_node_count; ++node)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t natural = 0; natural < 3; ++natural)
                center_gradient[node][physical] += hex8_signs[node][natural] * center_inverse[natural][physical] / 8.0;
    try {
        // Thermal stabilization uses the center metric and center measure;
        // its projection and the mechanical metrics retain volume averages.
        geometry.thermal_hourglass_coefficients =
            reduced_hex8_thermal_hourglass_coefficients(center_gradient, geometry.reduced_body_source_measure);
    } catch (const std::domain_error& error) {
        throw std::invalid_argument(std::string("Reduced HEX8 ") + error.what());
    }
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

namespace fuelsim::c3d8_detail {
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

C3d8Kinematics evaluate_cartesian_kinematics_from_gradient(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    C3d8Kinematics result{};
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

C3d8Kinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    return evaluate_cartesian_kinematics_from_gradient(point,
        displacement_gradient(point, current_state),
        committed_state,
        strain_formulation);
}
} // namespace fuelsim::c3d8_detail

namespace fuelsim {
using namespace c3d8_detail;

void c3d8_detail::validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    const double value = cartesian_detail::determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX8 deformation must preserve a positive Jacobian");
}
} // namespace fuelsim

namespace fuelsim::c3d8_detail {
void set_history_geometry(const elements::C3d8Input& input, bool reduced, elements::C3d8Result& result) {
    const auto& geometry = input.geometry;
    Hex8LocalAdValues current{}, old{};
    for (std::size_t i = 0; i < current.size(); ++i) {
        current[i] = input.state[i];
        old[i] = input.committed_state[i];
    }
    std::array<C3d8Kinematics, 8> points;
    for (std::size_t q = 0; q < points.size(); ++q) {
        points[q] = evaluate_cartesian_incremental_kinematics(geometry.points[q],
            current,
            input.committed_state,
            input.strain_formulation);
        result.current_volume += points[q].current_weighted_measure.value();
        result.committed_volume += evaluate_cartesian_incremental_kinematics(geometry.points[q],
            old,
            input.committed_state,
            input.strain_formulation)
                                       .current_weighted_measure.value();
    }
    if (input.strain_formulation == StrainFormulation::small) {
        result.current_volume = geometry.reference_volume;
        result.committed_volume = geometry.reference_volume;
    }
    for (std::size_t q = 0; q < result.history.size(); ++q) {
        const auto k = reduced ? evaluate_cartesian_incremental_kinematics(geometry.reduced_point,
                                     current,
                                     input.committed_state,
                                     input.strain_formulation)
                               : points[q];
        result.incremental_rotations[q] = {k.rotation.xx.value(),
            k.rotation.xy.value(),
            k.rotation.xz.value(),
            k.rotation.yx.value(),
            k.rotation.yy.value(),
            k.rotation.yz.value(),
            k.rotation.zx.value(),
            k.rotation.zy.value(),
            k.rotation.zz.value()};
    }
}
} // namespace fuelsim::c3d8_detail

namespace fuelsim::c3d8_detail {
namespace {
template <typename Scalar>
void reduced_hex8_metric_impl(const std::array<std::array<Scalar, 3>, 8>& average_gradient,
    std::array<std::array<Scalar, 3>, 3>& effective_mapping,
    std::array<std::array<Scalar, 3>, 3>& metric,
    std::array<Scalar, 3>& pivots,
    Scalar& leading_determinant,
    Scalar& numerator) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>);
    using std::isfinite;
    std::array<std::array<Scalar, 3>, 3> inverse_effective_mapping{};
    for (std::size_t natural = 0; natural < 3; ++natural)
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t node = 0; node < hex8_node_count; ++node)
                inverse_effective_mapping[natural][physical] +=
                    hex8_signs[node][natural] * average_gradient[node][physical];
    const Scalar inverse_effective_determinant = cartesian_detail::determinant(inverse_effective_mapping);
    if (!isfinite(inverse_effective_determinant) || inverse_effective_determinant == 0.0)
        throw std::domain_error("effective mapping must be nonsingular");
    effective_mapping = cartesian_detail::inverse(inverse_effective_mapping, inverse_effective_determinant);
    metric = {};
    for (std::size_t first = 0; first < 3; ++first)
        for (std::size_t second = 0; second < 3; ++second)
            for (std::size_t physical = 0; physical < 3; ++physical)
                metric[first][second] += effective_mapping[physical][first] * effective_mapping[physical][second];
    const Scalar first_pivot = metric[0][0];
    const Scalar second_pivot = metric[1][1] - metric[0][1] * metric[0][1] / first_pivot;
    leading_determinant = metric[0][0] * metric[1][1] - metric[0][1] * metric[0][1];
    numerator = metric[1][1] * metric[0][2] * metric[0][2] - 2.0 * metric[0][1] * metric[0][2] * metric[1][2]
                + metric[0][0] * metric[1][2] * metric[1][2];
    const Scalar third_pivot = metric[2][2] - numerator / leading_determinant;
    if (!isfinite(first_pivot) || !isfinite(second_pivot) || !isfinite(third_pivot) || !(first_pivot > 0.0)
        || !(second_pivot > 0.0) || !(third_pivot > 0.0))
        throw std::domain_error("effective metric must be positive definite");
    pivots = {first_pivot, second_pivot, third_pivot};
}

template <typename Scalar>
std::array<Scalar, 4> reduced_hex8_thermal_hourglass_coefficients_impl(const std::array<Scalar, 3>& pivots,
    const Scalar& volume) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>);
    const auto& first_pivot = pivots[0];
    const auto& second_pivot = pivots[1];
    const auto& third_pivot = pivots[2];
    const Scalar thermal_scale = volume / 192.0;
    const Scalar inverse_x = 1.0 / first_pivot, inverse_y = 1.0 / second_pivot, inverse_z = 1.0 / third_pivot;
    return {thermal_scale * (inverse_x + inverse_y),
        thermal_scale * (inverse_x + inverse_z),
        thermal_scale * (inverse_x + inverse_z),
        thermal_scale * (inverse_x + inverse_y + inverse_z) / 3.0};
}

template <typename Scalar>
std::array<std::array<Scalar, 4>, 8> hex8_hourglass_shape_impl(const std::array<std::array<Scalar, 3>, 8>& coordinates,
    const std::array<std::array<Scalar, 3>, 8>& average_gradient) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>);
    std::array<std::array<Scalar, 4>, 8> result{};
    for (std::size_t mode = 0; mode < 4; ++mode) {
        std::array<Scalar, 3> projected{};
        for (std::size_t node = 0; node < 8; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                projected[component] += coordinates[node][component] * hex8_raw_hourglass[node][mode];
        for (std::size_t node = 0; node < 8; ++node) {
            result[node][mode] = hex8_raw_hourglass[node][mode];
            for (std::size_t component = 0; component < 3; ++component)
                result[node][mode] -= average_gradient[node][component] * projected[component];
        }
    }
    return result;
}
} // namespace

std::array<double, 4> reduced_hex8_thermal_hourglass_coefficients(
    const std::array<std::array<double, 3>, hex8_node_count>& average_gradient,
    double volume) {
    return reduced_hex8_thermal_hourglass_coefficients(reduced_hex8_metric(average_gradient), volume);
}

std::array<adlite::Scalar, 4> reduced_hex8_thermal_hourglass_coefficients(
    const std::array<std::array<adlite::Scalar, 3>, hex8_node_count>& average_gradient,
    const adlite::Scalar& volume) {
    cartesian_detail::ActiveMatrix3 mapping{}, metric{};
    std::array<adlite::Scalar, 3> pivots{};
    adlite::Scalar leading, numerator;
    reduced_hex8_metric_impl(average_gradient, mapping, metric, pivots, leading, numerator);
    return reduced_hex8_thermal_hourglass_coefficients_impl(pivots, volume);
}

void hex8_shape_values(double xi,
    double eta,
    double zeta,
    std::array<double, 8>& shape,
    std::array<std::array<double, 3>, 8>& derivative) {
    for (std::size_t node = 0; node < 8; ++node) {
        const double sx = hex8_signs[node][0], sy = hex8_signs[node][1], sz = hex8_signs[node][2];
        shape[node] = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta);
        derivative[node] = {{0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
            0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta),
            0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)}};
    }
}

ReducedHex8Metric reduced_hex8_metric(const std::array<std::array<double, 3>, 8>& average_gradient) {
    ReducedHex8Metric result;
    reduced_hex8_metric_impl(average_gradient,
        result.mapping,
        result.metric,
        result.pivots,
        result.leading,
        result.numerator);
    return result;
}

std::array<double, 4> reduced_hex8_thermal_hourglass_coefficients(const ReducedHex8Metric& metric, double volume) {
    return reduced_hex8_thermal_hourglass_coefficients_impl(metric.pivots, volume);
}

std::array<std::array<double, 4>, 8> hex8_hourglass_shape(const std::array<std::array<double, 3>, 8>& coordinates,
    const std::array<std::array<double, 3>, 8>& average_gradient) {
    return hex8_hourglass_shape_impl(coordinates, average_gradient);
}

std::array<std::array<adlite::Scalar, 4>, 8> hex8_hourglass_shape(
    const std::array<std::array<adlite::Scalar, 3>, 8>& coordinates,
    const std::array<std::array<adlite::Scalar, 3>, 8>& average_gradient) {
    return hex8_hourglass_shape_impl(coordinates, average_gradient);
}
} // namespace fuelsim::c3d8_detail

namespace fuelsim::elements {
void add_hex8_body_acceleration(const C3d8Input& input, Hex8LocalResidual& residual) {
    for (double value : input.body_acceleration)
        if (!std::isfinite(value))
            throw std::invalid_argument("Body acceleration must be finite");
    if (input.body_acceleration == std::array<double, 3>{})
        return;
    // Use the exact cancellation of current density and current volume.
    for (const auto& point : input.geometry.points) {
        const double mass = point.weighted_measure
                            * input.material.initial_density(input.initial_temperature,
                                {0.0, point.position.x, point.position.y, point.position.z});
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < 8; ++node)
                residual[8 + component * 8 + node] -= mass * point.shape[node] * input.body_acceleration[component];
    }
}
} // namespace fuelsim::elements

namespace fuelsim::cartesian_detail {
void evaluate_hex20_shapes(double xi,
    double eta,
    double zeta,
    std::array<double, 20>& shape,
    std::array<std::array<double, 3>, 20>& derivative) {
    for (std::size_t node = 0; node < 8; ++node) {
        const double sx = c3d8_detail::hex8_signs[node][0], sy = c3d8_detail::hex8_signs[node][1],
                     sz = c3d8_detail::hex8_signs[node][2];
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
            0.25 * sy * (1.0 - xi * xi) * (1.0 + sz * zeta),
            0.25 * sz * (1.0 - xi * xi) * (1.0 + sy * eta)}};
    };
    const auto eta_edge = [&](std::size_t node, double sx, double sz) {
        shape[node] = 0.25 * (1.0 - eta * eta) * (1.0 + sx * xi) * (1.0 + sz * zeta);
        derivative[node] = {{0.25 * sx * (1.0 - eta * eta) * (1.0 + sz * zeta),
            -0.5 * eta * (1.0 + sx * xi) * (1.0 + sz * zeta),
            0.25 * sz * (1.0 - eta * eta) * (1.0 + sx * xi)}};
    };
    const auto zeta_edge = [&](std::size_t node, double sx, double sy) {
        shape[node] = 0.25 * (1.0 - zeta * zeta) * (1.0 + sx * xi) * (1.0 + sy * eta);
        derivative[node] = {{0.25 * sx * (1.0 - zeta * zeta) * (1.0 + sy * eta),
            0.25 * sy * (1.0 - zeta * zeta) * (1.0 + sx * xi),
            -0.5 * zeta * (1.0 + sx * xi) * (1.0 + sy * eta)}};
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

} // namespace fuelsim::cartesian_detail
