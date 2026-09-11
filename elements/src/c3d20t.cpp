#include "c3d20t.hpp"
#include "ad_local_system.hpp"
#include "c3d_common.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
using namespace cartesian_detail;
constexpr double gauss3 = 0.774596669241483377035853079956479922;
constexpr double gauss2 = 0.577350269189625764509148780502;
constexpr std::array<double, 2> gauss2_points = {-gauss2, gauss2};
constexpr std::array<double, 2> gauss2_weights = {1.0, 1.0};
constexpr std::array<double, 3> gauss3_points = {-gauss3, 0.0, gauss3};
constexpr std::array<double, 3> gauss3_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
constexpr std::array<std::array<double, 3>, 8> corner_signs = {{{{-1.0, -1.0, -1.0}},
    {{1.0, -1.0, -1.0}},
    {{1.0, 1.0, -1.0}},
    {{-1.0, 1.0, -1.0}},
    {{-1.0, -1.0, 1.0}},
    {{1.0, -1.0, 1.0}},
    {{1.0, 1.0, 1.0}},
    {{-1.0, 1.0, 1.0}}}};

void evaluate_hex20_shapes(double xi,
    double eta,
    double zeta,
    std::array<double, 20>& shape,
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

void evaluate_hex8_temperature_shapes(double xi,
    double eta,
    double zeta,
    std::array<double, 8>& shape,
    std::array<std::array<double, 3>, 8>& derivative) {
    for (std::size_t node = 0; node < 8; ++node) {
        const double sx = corner_signs[node][0], sy = corner_signs[node][1], sz = corner_signs[node][2];
        shape[node] = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta);
        derivative[node] = {{0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
            0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta),
            0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)}};
    }
}

struct Hex20ReferenceMapping final {
    std::array<double, 20> displacement_shape;
    std::array<std::array<double, 3>, 20> displacement_derivative;
    std::array<double, 8> temperature_shape;
    std::array<std::array<double, 3>, 8> temperature_derivative;
    CartesianPoint3 position;
    Matrix3 inverse_jacobian;
    double determinant;
};

Hex20ReferenceMapping evaluate_hex20_mapping(const Hex20Coordinates& coordinates, double xi, double eta, double zeta) {
    Hex20ReferenceMapping result{};
    evaluate_hex20_shapes(xi, eta, zeta, result.displacement_shape, result.displacement_derivative);
    evaluate_hex8_temperature_shapes(xi, eta, zeta, result.temperature_shape, result.temperature_derivative);
    Matrix3 jacobian{};
    for (std::size_t node = 0; node < 20; ++node) {
        result.position.x += result.displacement_shape[node] * coordinates[node].x;
        result.position.y += result.displacement_shape[node] * coordinates[node].y;
        result.position.z += result.displacement_shape[node] * coordinates[node].z;
        const std::array<double, 3> coordinate = {coordinates[node].x, coordinates[node].y, coordinates[node].z};
        for (std::size_t physical = 0; physical < 3; ++physical)
            for (std::size_t natural = 0; natural < 3; ++natural)
                jacobian[physical][natural] += coordinate[physical] * result.displacement_derivative[node][natural];
    }
    result.determinant = determinant(jacobian);
    if (!std::isfinite(result.determinant) || !(result.determinant > 0.0))
        throw std::invalid_argument("Hex20Geometry requires a finite positive Jacobian determinant");
    result.inverse_jacobian = inverse(jacobian, result.determinant);
    return result;
}

} // namespace

Hex20Geometry elements::make_c3d20t_geometry(const Hex20Coordinates& coordinates, C3d20Quadrature quadrature) {
    if (quadrature != C3d20Quadrature::full && quadrature != C3d20Quadrature::reduced)
        throw std::invalid_argument("Invalid C3D20 quadrature");
    const std::size_t order = quadrature == C3d20Quadrature::full ? 3 : 2;
    if (order != 2 && order != 3)
        throw std::invalid_argument("Quadratic hexahedron requires Gauss order two or three");
    const bool reduced = order == 2;
    const double* points = reduced ? gauss2_points.data() : gauss3_points.data();
    const double* weights = reduced ? gauss2_weights.data() : gauss3_weights.data();
    Hex20Geometry geometry{};
    std::size_t thermal_q = 0;
    for (std::size_t kz = 0; kz < 2; ++kz)
        for (std::size_t ky = 0; ky < 2; ++ky)
            for (std::size_t kx = 0; kx < 2; ++kx) {
                const double xi = gauss2_points[kx], eta = gauss2_points[ky], zeta = gauss2_points[kz];
                const Hex20ReferenceMapping mapping = evaluate_hex20_mapping(coordinates, xi, eta, zeta);
                Hex20ThermalQuadraturePoint& point = geometry.thermal_points[thermal_q++];
                point.temperature_shape = mapping.temperature_shape;
                point.position = mapping.position;
                point.weighted_measure =
                    mapping.determinant * gauss2_weights[kx] * gauss2_weights[ky] * gauss2_weights[kz];
                for (std::size_t node = 0; node < 8; ++node)
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            point.temperature_gradient[node][physical] += mapping.temperature_derivative[node][natural]
                                                                          * mapping.inverse_jacobian[natural][physical];
            }
    geometry.mechanical_points.resize(order * order * order);
    std::size_t mechanical_q = 0;
    for (std::size_t kz = 0; kz < order; ++kz)
        for (std::size_t ky = 0; ky < order; ++ky)
            for (std::size_t kx = 0; kx < order; ++kx) {
                const double xi = points[kx], eta = points[ky], zeta = points[kz];
                const Hex20ReferenceMapping mapping = evaluate_hex20_mapping(coordinates, xi, eta, zeta);
                Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[mechanical_q++];
                point.temperature_shape = mapping.temperature_shape;
                point.displacement_shape = mapping.displacement_shape;
                point.position = mapping.position;
                point.weighted_measure = mapping.determinant * weights[kx] * weights[ky] * weights[kz];
                Matrix3 source_jacobian{};
                for (std::size_t node = 0; node < 8; ++node) {
                    const std::array<double, 3> coordinate = {coordinates[node].x,
                        coordinates[node].y,
                        coordinates[node].z};
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            source_jacobian[physical][natural] +=
                                coordinate[physical] * mapping.temperature_derivative[node][natural];
                }
                const double source_determinant = determinant(source_jacobian);
                if (!std::isfinite(source_determinant) || !(source_determinant > 0.0))
                    throw std::invalid_argument(
                        "Hex20Geometry thermal source corners require a finite positive Jacobian determinant");
                const Matrix3 source_inverse = inverse(source_jacobian, source_determinant);
                point.source_weighted_measure = source_determinant * weights[kx] * weights[ky] * weights[kz];
                for (std::size_t node = 0; node < 8; ++node)
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural) {
                            point.temperature_gradient[node][physical] += mapping.temperature_derivative[node][natural]
                                                                          * mapping.inverse_jacobian[natural][physical];
                            point.source_displacement_gradient[node][physical] +=
                                mapping.temperature_derivative[node][natural] * source_inverse[natural][physical];
                        }
                for (std::size_t node = 0; node < 20; ++node)
                    for (std::size_t physical = 0; physical < 3; ++physical)
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            point.displacement_gradient[node][physical] +=
                                mapping.displacement_derivative[node][natural]
                                * mapping.inverse_jacobian[natural][physical];
            }
    return geometry;
}
} // namespace fuelsim

namespace fuelsim::c3d20_detail {
namespace {
using cartesian_detail::ActiveMatrix3;
using cartesian_detail::Matrix3;

struct Hex20KinematicsValues final {
    SymmetricTensor3Values strain_increment{};
    CartesianRotation rotation{};
    std::array<std::array<double, 3>, 20> current_gradient{};
    double current_weighted_measure = 0.0;
};

struct Hex20Kinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, 20> current_gradient;
    Matrix3 current_inverse_values{}, midpoint_inverse_values{};
    adlite::Scalar current_weighted_measure{0.0};
};

struct Hex20SourceMeasureValues final {
    Matrix3 measure_derivative{};
    double weighted_measure = 0.0;
};

struct Hex20FiniteThermalKinematicsValues final {
    std::array<std::array<double, 3>, 8> midpoint_temperature_gradient{};
    double current_weighted_measure = 0.0;
};

ActiveMatrix3 displacement_gradient(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalAdValues& state);

Hex20KinematicsValues evaluate_kinematics_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation);

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point,
    const ActiveMatrix3& gradient,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation);

Hex20Kinematics evaluate_kinematics(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalAdValues& state,
    const Hex20LocalValues& committed_state,
    StrainFormulation strain_formulation);

Hex20SourceMeasureValues evaluate_source_measure_values(const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state);

std::array<std::array<double, 3>, 8> temperature_shape_gradients(const Hex20MechanicalQuadraturePoint& point,
    const Matrix3& inverse_map);

Hex20FiniteThermalKinematicsValues evaluate_finite_thermal_kinematics_values(
    const Hex20MechanicalQuadraturePoint& point,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state);
} // namespace
} // namespace fuelsim::c3d20_detail

namespace fuelsim::c3d20_detail {
namespace {
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
        Matrix3 rotation{};
        std::array<double, 6> strain{};
        cartesian_detail::hughes_winget_rotation(hughes_winget, rotation, strain);
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

} // namespace
} // namespace fuelsim::c3d20_detail

namespace fuelsim {
using namespace c3d20_detail;
using cartesian_detail::determinant;

void elements::validate_c3d20t_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state) {
    const double value = determinant(deformation_gradient(point, state));
    if (!std::isfinite(value) || !(value > 0.0))
        throw std::domain_error("Finite-strain HEX20 deformation must preserve a positive Jacobian");
}
} // namespace fuelsim

namespace fuelsim::c3d20_detail {

namespace {
using cartesian_detail::ActiveMatrix3;
using cartesian_detail::determinant;
using cartesian_detail::inverse;
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
            stress = fuelsim::rotate_cartesian_tensor_values(stress, kinematics.rotation);
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
    const fuelsim::CartesianStressTangent tangent = fuelsim::evaluate_stress_tangent(material,
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

} // namespace fuelsim::c3d20_detail

namespace fuelsim::elements {
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request, C3d20Quadrature quadrature) {
    using namespace c3d20_detail;
    if (quadrature != C3d20Quadrature::full && quadrature != C3d20Quadrature::reduced)
        throw std::invalid_argument("Invalid C3D20 quadrature");
    const std::size_t expected = quadrature == C3d20Quadrature::full ? 27 : 8;
    if (input.geometry.mechanical_points.size() != expected)
        throw std::invalid_argument("C3D20 geometry does not match the selected quadrature");
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
} // namespace fuelsim::elements

namespace fuelsim::elements {
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request) {
    return evaluate_c3d20t(input, request, C3d20Quadrature::full);
}

Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates) {
    return make_c3d20t_geometry(coordinates, C3d20Quadrature::full);
}
} // namespace fuelsim::elements
