#include "quad8_face.hpp"
#include "ad_local_system.hpp"
#include "contact_common.hpp"
#include "contact_types.hpp"
#include "quad4_face.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace fuelsim {
namespace {
constexpr double gauss3 = 0.774596669241483377035853079956479922;
constexpr double gauss2 = 0.577350269189625764509148780502;
constexpr std::array<double, 2> gauss2_points = {-gauss2, gauss2};
constexpr std::array<double, 2> gauss2_weights = {1.0, 1.0};
constexpr std::array<double, 3> gauss3_points = {-gauss3, 0.0, gauss3};
constexpr std::array<double, 3> gauss3_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

void evaluate_quad8_shapes(double xi,
    double eta,
    std::array<double, 8>& shape,
    std::array<double, 8>& derivative_xi,
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

Quad8FaceGeometry make_quad8_face_geometry(const Quad8FaceCoordinates& coordinates) {
    Quad8FaceGeometry geometry{};
    std::size_t q = 0;
    for (std::size_t ky = 0; ky < 2; ++ky)
        for (std::size_t kx = 0; kx < 2; ++kx) {
            const double xi = gauss2_points[kx], eta = gauss2_points[ky];
            Quad8FaceThermalQuadraturePoint& point = geometry.thermal_points[q++];
            evaluate_quad8_shapes(xi, eta, point.displacement_shape, point.derivative_xi, point.derivative_eta);
            point.temperature_shape = {{0.25 * (1.0 - xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 + eta),
                0.25 * (1.0 - xi) * (1.0 + eta)}};
            CartesianPoint3 tangent_xi{}, tangent_eta{};
            for (std::size_t node = 0; node < 8; ++node) {
                tangent_xi.x += point.derivative_xi[node] * coordinates[node].x;
                tangent_xi.y += point.derivative_xi[node] * coordinates[node].y;
                tangent_xi.z += point.derivative_xi[node] * coordinates[node].z;
                tangent_eta.x += point.derivative_eta[node] * coordinates[node].x;
                tangent_eta.y += point.derivative_eta[node] * coordinates[node].y;
                tangent_eta.z += point.derivative_eta[node] * coordinates[node].z;
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
        for (std::size_t kx = 0; kx < 3; ++kx)
            geometry.mechanical_points[q++] = make_quad8_face_mechanical_point(coordinates,
                gauss3_points[kx],
                gauss3_points[ky],
                gauss3_weights[kx] * gauss3_weights[ky]);
    return geometry;
}

Quad8FaceMechanicalQuadraturePoint make_quad8_face_mechanical_point(const Quad8FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight) {
    Quad8FaceMechanicalQuadraturePoint point{};
    evaluate_quad8_shapes(xi, eta, point.displacement_shape, point.derivative_xi, point.derivative_eta);
    for (std::size_t node = 0; node < 8; ++node) {
        point.tangent_xi.x += point.derivative_xi[node] * coordinates[node].x;
        point.tangent_xi.y += point.derivative_xi[node] * coordinates[node].y;
        point.tangent_xi.z += point.derivative_xi[node] * coordinates[node].z;
        point.tangent_eta.x += point.derivative_eta[node] * coordinates[node].x;
        point.tangent_eta.y += point.derivative_eta[node] * coordinates[node].y;
        point.tangent_eta.z += point.derivative_eta[node] * coordinates[node].z;
    }
    const CartesianPoint3 area{point.tangent_xi.y * point.tangent_eta.z - point.tangent_xi.z * point.tangent_eta.y,
        point.tangent_xi.z * point.tangent_eta.x - point.tangent_xi.x * point.tangent_eta.z,
        point.tangent_xi.x * point.tangent_eta.y - point.tangent_xi.y * point.tangent_eta.x};
    const double measure = std::sqrt(area.x * area.x + area.y * area.y + area.z * area.z);
    if (!std::isfinite(measure) || !(measure > 0.0) || !std::isfinite(quadrature_weight) || !(quadrature_weight > 0.0))
        throw std::invalid_argument("Quad8 face mechanical point requires a finite positive weighted measure");
    point.quadrature_weight = quadrature_weight;
    return point;
}

Quad8FaceLocalResidual compute_quad8_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad8FaceGeometry& geometry,
    const Quad8FaceLocalValues& state,
    Quad8FaceLocalJacobian* jacobian) {
    Quad8FaceLocalAdValues ad_state{};
    if (jacobian == nullptr)
        ad_local_system::make_passive(state.data(), state.size(), ad_state.data());
    else
        ad_local_system::make_active(state.data(), state.size(), ad_state.data());
    Quad8FaceLocalAdValues residual{};
    residual.fill(adlite::Scalar(0.0));
    if (data.kind == Quad4FaceBoundaryKind::surface_heat_flux)
        throw std::invalid_argument("HEX20 surface heat flux is not implemented");
    if (data.kind == Quad4FaceBoundaryKind::convection) {
        for (const Quad8FaceThermalQuadraturePoint& point : geometry.thermal_points) {
            adlite::Scalar temperature = 0.0;
            for (std::size_t node = 0; node < 4; ++node)
                temperature += point.temperature_shape[node] * ad_state[node];
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

namespace fuelsim {
namespace {
using ActivePoint3 = std::array<adlite::Scalar, 3>;
adlite::Scalar quad8_disk_fraction_ad(const Quad8ToQuad8HeatGeometry& geometry,
    const std::array<ActivePoint3, 16>& nodes);

struct Quad8ShapeValues final {
    std::array<adlite::Scalar, 8> shape{}, derivative_xi{}, derivative_eta{};
    std::array<adlite::Scalar, 8> second_xi{}, second_xi_eta{}, second_eta{};
};

struct SurfaceProjection8 final {
    bool projected = false;
    bool xi_clamped = false, eta_clamped = false;
    std::array<adlite::Scalar, 8> primary_shape{};
    ActivePoint3 primary_point{}, normal{}, tangent_xi{};
    adlite::Scalar gap{0.0}, xi{0.0}, eta{0.0};
};

struct CartesianContactAdValue8 final {
    bool projected = false, sliding = false;
    std::array<adlite::Scalar, 8> primary_shape{};
    ActivePoint3 normal{}, tangent_first{}, tangential_slip{}, elastic_tangential_slip{}, tangential_traction_vector{};
    adlite::Scalar gap{0.0}, pressure{0.0}, tributary_area{0.0}, contact_force{0.0}, tangential_traction{0.0},
        tangential_force{0.0}, friction_dissipation{0.0};
};

struct SurfaceBasis final {
    ActivePoint3 first{}, second{};
};

struct DoubleQuad8ShapeValues final {
    std::array<double, 8> shape{}, derivative_xi{}, derivative_eta{};
    std::array<double, 8> second_xi{}, second_xi_eta{}, second_eta{};
};

Quad8SurfaceContactLocalAdValues make_ad_state(const Quad8SurfaceContactLocalValues& state, bool derivatives);

ActivePoint3 subtract(const ActivePoint3& first, const ActivePoint3& second) {
    return {first[0] - second[0], first[1] - second[1], first[2] - second[2]};
}

adlite::Scalar dot(const ActivePoint3& first, const ActivePoint3& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

ActivePoint3 cross(const ActivePoint3& first, const ActivePoint3& second) {
    return {first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

adlite::Scalar norm(const ActivePoint3& value) {
    return adlite::hypot(adlite::hypot(value[0], value[1]), value[2]);
}

void quad8_shape(const adlite::Scalar& xi, const adlite::Scalar& eta, Quad8ShapeValues& result) {
    const adlite::Scalar xm = 1.0 - xi, xp = 1.0 + xi, ym = 1.0 - eta, yp = 1.0 + eta;
    result.shape = {0.25 * xm * ym * (-xi - eta - 1.0),
        0.25 * xp * ym * (xi - eta - 1.0),
        0.25 * xp * yp * (xi + eta - 1.0),
        0.25 * xm * yp * (-xi + eta - 1.0),
        0.5 * (1.0 - xi * xi) * ym,
        0.5 * xp * (1.0 - eta * eta),
        0.5 * (1.0 - xi * xi) * yp,
        0.5 * xm * (1.0 - eta * eta)};
    result.derivative_xi = {0.25 * ym * (2.0 * xi + eta),
        0.25 * ym * (2.0 * xi - eta),
        0.25 * yp * (2.0 * xi + eta),
        0.25 * yp * (2.0 * xi - eta),
        -xi * ym,
        0.5 * (1.0 - eta * eta),
        -xi * yp,
        -0.5 * (1.0 - eta * eta)};
    result.derivative_eta = {0.25 * xm * (xi + 2.0 * eta),
        0.25 * xp * (-xi + 2.0 * eta),
        0.25 * xp * (xi + 2.0 * eta),
        0.25 * xm * (-xi + 2.0 * eta),
        -0.5 * (1.0 - xi * xi),
        -xp * eta,
        0.5 * (1.0 - xi * xi),
        -xm * eta};
    result.second_xi = {0.5 * ym, 0.5 * ym, 0.5 * yp, 0.5 * yp, -ym, 0.0, -yp, 0.0};
    result.second_xi_eta = {0.25 * (1.0 - 2.0 * xi - 2.0 * eta),
        -0.25 * (1.0 + 2.0 * xi - 2.0 * eta),
        0.25 * (1.0 + 2.0 * xi + 2.0 * eta),
        0.25 * (-1.0 + 2.0 * xi - 2.0 * eta),
        xi,
        -eta,
        -xi,
        eta};
    result.second_eta = {0.5 * xm, 0.5 * xp, 0.5 * xp, 0.5 * xm, 0.0, -xp, 0.0, -xm};
}

void double_quad8_shape(double xi, double eta, DoubleQuad8ShapeValues& result) {
    const double xm = 1.0 - xi, xp = 1.0 + xi, ym = 1.0 - eta, yp = 1.0 + eta;
    result.shape = {0.25 * xm * ym * (-xi - eta - 1.0),
        0.25 * xp * ym * (xi - eta - 1.0),
        0.25 * xp * yp * (xi + eta - 1.0),
        0.25 * xm * yp * (-xi + eta - 1.0),
        0.5 * (1.0 - xi * xi) * ym,
        0.5 * xp * (1.0 - eta * eta),
        0.5 * (1.0 - xi * xi) * yp,
        0.5 * xm * (1.0 - eta * eta)};
    result.derivative_xi = {0.25 * ym * (2.0 * xi + eta),
        0.25 * ym * (2.0 * xi - eta),
        0.25 * yp * (2.0 * xi + eta),
        0.25 * yp * (2.0 * xi - eta),
        -xi * ym,
        0.5 * (1.0 - eta * eta),
        -xi * yp,
        -0.5 * (1.0 - eta * eta)};
    result.derivative_eta = {0.25 * xm * (xi + 2.0 * eta),
        0.25 * xp * (-xi + 2.0 * eta),
        0.25 * xp * (xi + 2.0 * eta),
        0.25 * xm * (-xi + 2.0 * eta),
        -0.5 * (1.0 - xi * xi),
        -xp * eta,
        0.5 * (1.0 - xi * xi),
        -xm * eta};
    result.second_xi = {0.5 * ym, 0.5 * ym, 0.5 * yp, 0.5 * yp, -ym, 0.0, -yp, 0.0};
    result.second_xi_eta = {0.25 * (1.0 - 2.0 * xi - 2.0 * eta),
        -0.25 * (1.0 + 2.0 * xi - 2.0 * eta),
        0.25 * (1.0 + 2.0 * xi + 2.0 * eta),
        0.25 * (-1.0 + 2.0 * xi - 2.0 * eta),
        xi,
        -eta,
        -xi,
        eta};
    result.second_eta = {0.5 * xm, 0.5 * xp, 0.5 * xp, 0.5 * xm, 0.0, -xp, 0.0, -xm};
}

CartesianPoint3 double_interpolate(const std::array<CartesianPoint3, 8>& coordinates,
    const std::array<double, 8>& coefficients) {
    CartesianPoint3 result{};
    for (std::size_t node = 0; node < coordinates.size(); ++node) {
        result.x += coefficients[node] * coordinates[node].x;
        result.y += coefficients[node] * coordinates[node].y;
        result.z += coefficients[node] * coordinates[node].z;
    }
    return result;
}

CartesianPoint3 double_subtract(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

double double_dot(const CartesianPoint3& first, const CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

CartesianPoint3 double_cross(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

void quad4_temperature_shape(const adlite::Scalar& xi,
    const adlite::Scalar& eta,
    std::array<adlite::Scalar, 4>& shape) {
    shape = {0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)};
}

std::array<ActivePoint3, 16> current_nodes(const std::array<CartesianPoint3, 8>& secondary,
    const std::array<CartesianPoint3, 8>& primary,
    const Quad8SurfaceContactLocalAdValues& state) {
    std::array<ActivePoint3, 16> result{};
    for (std::size_t node = 0; node < 16; ++node) {
        const CartesianPoint3& reference = node < 8 ? secondary[node] : primary[node - 8];
        result[node] = {reference.x + state[8 + node], reference.y + state[24 + node], reference.z + state[40 + node]};
    }
    return result;
}

ActivePoint3 interpolate_point(const std::array<ActivePoint3, 16>& nodes,
    std::size_t offset,
    const std::array<adlite::Scalar, 8>& shape) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += shape[node] * nodes[offset + node][component];
    return result;
}

ActivePoint3 interpolate_primary(const std::array<ActivePoint3, 16>& nodes,
    const std::array<adlite::Scalar, 8>& coefficients) {
    return interpolate_point(nodes, 8, coefficients);
}

SurfaceProjection8 project_to_primary(const ActivePoint3& secondary_point,
    const std::array<ActivePoint3, 16>& nodes,
    double normal_orientation,
    bool allow_extrapolation = false) {
    adlite::Scalar xi = 0.0, eta = 0.0;
    bool converged = false;
    // Translate the Newton geometry to one primary node so global translations
    // do not set the attainable projection increment through cancellation.
    std::array<ActivePoint3, 16> local_nodes{};
    for (std::size_t node = 8; node < 16; ++node)
        local_nodes[node] = subtract(nodes[node], nodes[8]);
    const ActivePoint3 local_secondary = subtract(secondary_point, nodes[8]);
    for (std::size_t iteration = 0; iteration < 16; ++iteration) {
        Quad8ShapeValues values;
        quad8_shape(xi, eta, values);
        const ActivePoint3 point = interpolate_primary(local_nodes, values.shape);
        const ActivePoint3 tangent_xi = interpolate_primary(local_nodes, values.derivative_xi);
        const ActivePoint3 tangent_eta = interpolate_primary(local_nodes, values.derivative_eta);
        const ActivePoint3 tangent_xi_xi = interpolate_primary(local_nodes, values.second_xi);
        const ActivePoint3 tangent_xi_eta = interpolate_primary(local_nodes, values.second_xi_eta);
        const ActivePoint3 tangent_eta_eta = interpolate_primary(local_nodes, values.second_eta);
        const ActivePoint3 difference = subtract(local_secondary, point);
        const adlite::Scalar residual_xi = dot(difference, tangent_xi), residual_eta = dot(difference, tangent_eta);
        const adlite::Scalar jacobian_xi_xi = -dot(tangent_xi, tangent_xi) + dot(difference, tangent_xi_xi);
        const adlite::Scalar jacobian_xi_eta = -dot(tangent_eta, tangent_xi) + dot(difference, tangent_xi_eta);
        const adlite::Scalar jacobian_eta_xi = -dot(tangent_xi, tangent_eta) + dot(difference, tangent_xi_eta);
        const adlite::Scalar jacobian_eta_eta = -dot(tangent_eta, tangent_eta) + dot(difference, tangent_eta_eta);
        const adlite::Scalar determinant = jacobian_xi_xi * jacobian_eta_eta - jacobian_xi_eta * jacobian_eta_xi;
        if (!std::isfinite(determinant.value()) || std::abs(determinant.value()) <= std::numeric_limits<double>::min())
            throw std::domain_error("HEX20 contact projection has a singular Q8 surface Jacobian");
        const adlite::Scalar delta_xi =
                                 (-residual_xi * jacobian_eta_eta + jacobian_xi_eta * residual_eta) / determinant,
                             delta_eta = (-jacobian_xi_xi * residual_eta + jacobian_eta_xi * residual_xi) / determinant;
        xi += delta_xi;
        eta += delta_eta;
        if (contact_common::projection_increment_converged(delta_xi,
                delta_eta,
                xi,
                eta,
                quad8_surface_contact_local_dof_count)) {
            converged = true;
            break;
        }
    }
    if (!converged)
        throw std::domain_error("Contact projection Newton iteration did not converge");
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(xi.value()) || !std::isfinite(eta.value())
        || (!allow_extrapolation
            && (xi.value() < -1.0 - tolerance || xi.value() > 1.0 + tolerance || eta.value() < -1.0 - tolerance
                || eta.value() > 1.0 + tolerance)))
        return {};
    const bool xi_clamped = !allow_extrapolation && (xi.value() < -1.0 || xi.value() > 1.0);
    const bool eta_clamped = !allow_extrapolation && (eta.value() < -1.0 || eta.value() > 1.0);
    if (xi_clamped)
        xi = xi.value() < 0.0 ? -1.0 : 1.0;
    if (eta_clamped)
        eta = eta.value() < 0.0 ? -1.0 : 1.0;
    Quad8ShapeValues values;
    quad8_shape(xi, eta, values);
    const ActivePoint3 primary_point = interpolate_primary(nodes, values.shape);
    const ActivePoint3 tangent_xi = interpolate_primary(nodes, values.derivative_xi);
    const ActivePoint3 tangent_eta = interpolate_primary(nodes, values.derivative_eta);
    const ActivePoint3 area_vector = cross(tangent_xi, tangent_eta);
    const adlite::Scalar measure = norm(area_vector);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("HEX20 contact primary Q8 face has a nonpositive current measure");
    SurfaceProjection8 result;
    result.projected = true;
    result.xi_clamped = xi_clamped;
    result.eta_clamped = eta_clamped;
    result.primary_shape = values.shape;
    result.primary_point = primary_point;
    result.tangent_xi = tangent_xi;
    result.xi = xi;
    result.eta = eta;
    for (std::size_t component = 0; component < 3; ++component)
        result.normal[component] = normal_orientation * area_vector[component] / measure;
    result.gap = dot(subtract(primary_point, secondary_point), result.normal);
    return result;
}

adlite::Scalar current_surface_measure(const std::array<ActivePoint3, 16>& nodes,
    const std::array<adlite::Scalar, 8>& derivative_xi,
    const std::array<adlite::Scalar, 8>& derivative_eta,
    std::size_t offset) {
    const ActivePoint3 tangent_xi = interpolate_point(nodes, offset, derivative_xi);
    const ActivePoint3 tangent_eta = interpolate_point(nodes, offset, derivative_eta);
    const adlite::Scalar measure = norm(cross(tangent_xi, tangent_eta));
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("HEX20 contact secondary Q8 face has a nonpositive current measure");
    return measure;
}

adlite::Scalar temperature(const Quad8SurfaceContactLocalAdValues& state,
    std::size_t offset,
    const std::array<adlite::Scalar, 4>& shape) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        result += shape[node] * state[offset + node];
    return result;
}

std::array<adlite::Scalar, 8> active_values(const std::array<double, 8>& values) {
    std::array<adlite::Scalar, 8> result{};
    for (std::size_t node = 0; node < 8; ++node)
        result[node] = values[node];
    return result;
}

std::array<adlite::Scalar, 4> active_temperature_values(const std::array<double, 4>& values) {
    std::array<adlite::Scalar, 4> result{};
    for (std::size_t node = 0; node < 4; ++node)
        result[node] = values[node];
    return result;
}

struct HeatAdValue8 final {
    bool projected = false;
    std::array<adlite::Scalar, 4> primary_temperature_shape{};
    adlite::Scalar gap{0.0}, heat_flux{0.0}, weighted_measure{0.0};
    adlite::Scalar secondary_temperature{0.0}, primary_temperature{0.0};
    adlite::Scalar transfer_fraction{1.0};
};

HeatAdValue8 evaluate_heat_geometry(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalAdValues& state,
    bool disk_transfer = false) {
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    Quad8ShapeValues secondary_shape;
    secondary_shape.shape = {};
    for (std::size_t node = 0; node < 8; ++node)
        secondary_shape.shape[node] = geometry.secondary_displacement_shape[node];
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape.shape);
    const adlite::Scalar fraction = disk_transfer ? quad8_disk_fraction_ad(geometry, nodes) : adlite::Scalar(1.0);
    if (fraction.value() == 0.0) {
        HeatAdValue8 empty{};
        empty.projected = true;
        empty.transfer_fraction = fraction;
        return empty;
    }
    const SurfaceProjection8 projection =
        project_to_primary(secondary_point, nodes, geometry.normal_orientation, disk_transfer);
    if (!projection.projected)
        return {};
    std::array<adlite::Scalar, 4> primary_temperature_shape{};
    quad4_temperature_shape(projection.xi, projection.eta, primary_temperature_shape);
    adlite::Scalar primary_temperature = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        primary_temperature += primary_temperature_shape[node] * state[4 + node];
    const adlite::Scalar secondary_temperature =
        temperature(state, 0, active_temperature_values(geometry.secondary_temperature_shape));
    const adlite::Scalar measure = current_surface_measure(nodes,
        active_values(geometry.secondary_derivative_xi),
        active_values(geometry.secondary_derivative_eta),
        0);
    return {true,
        primary_temperature_shape,
        projection.gap,
        0.0,
        fraction * measure * geometry.quadrature_weight,
        secondary_temperature,
        primary_temperature,
        fraction};
}

HeatAdValue8 evaluate_heat(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalAdValues& state) {
    HeatAdValue8 result = evaluate_heat_geometry(geometry, state);
    if (result.projected)
        result.heat_flux = contact_common::gap_conductance(properties,
                               result.gap,
                               result.secondary_temperature,
                               result.primary_temperature)
                           * (result.secondary_temperature - result.primary_temperature);
    return result;
}

CartesianContactAdValue8 evaluate_mechanical(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history);

void apply_friction(const NormalContactProperties& properties,
    const ContactPointHistory& history,
    const ActivePoint3& transported_history,
    const ActivePoint3& transported_total_history,
    const ActivePoint3& relative_increment,
    CartesianContactAdValue8& result) {
    const auto friction = contact_common::friction_return(properties,
        history,
        transported_history,
        transported_total_history,
        relative_increment,
        result.normal,
        result.pressure,
        result.tributary_area);
    result.sliding = friction.sliding;
    result.tangential_slip = friction.tangential_slip;
    result.elastic_tangential_slip = friction.elastic_tangential_slip;
    result.tangential_traction_vector = friction.tangential_traction_vector;
    result.tangential_traction = friction.tangential_traction;
    result.tangential_force = friction.tangential_force;
    result.friction_dissipation = friction.friction_dissipation;
}

ActivePoint3 stored_history(const ContactPointHistory& history) {
    return {history.cartesian_elastic_tangential_slip[0],
        history.cartesian_elastic_tangential_slip[1],
        history.cartesian_elastic_tangential_slip[2]};
}

ActivePoint3 stored_total_history(const ContactPointHistory& history) {
    return {history.cartesian_total_tangential_slip[0],
        history.cartesian_total_tangential_slip[1],
        history.cartesian_total_tangential_slip[2]};
}

SurfaceBasis surface_basis(const SurfaceProjection8& projection) {
    const adlite::Scalar measure = norm(projection.tangent_xi);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("HEX20 primary contact surface has an undefined convected tangent");
    SurfaceBasis result;
    for (std::size_t component = 0; component < 3; ++component)
        result.first[component] = projection.tangent_xi[component] / measure;
    result.second = cross(projection.normal, result.first);
    return result;
}

SurfaceBasis stored_surface_basis(const ContactPointHistory& history, const SurfaceBasis& fallback) {
    if (!history.cartesian_tangent_basis_initialized)
        return fallback;
    SurfaceBasis result;
    for (std::size_t component = 0; component < 3; ++component) {
        result.first[component] = history.cartesian_contact_tangent_first[component];
        result.second[component] = history.cartesian_contact_normal[component];
    }
    result.second = cross(result.second, result.first);
    return result;
}

ActivePoint3
transport_surface_vector(const ActivePoint3& vector, const SurfaceBasis& current, const SurfaceBasis& committed) {
    const adlite::Scalar first_component = dot(vector, committed.first),
                         second_component = dot(vector, committed.second);
    ActivePoint3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        result[component] = first_component * current.first[component] + second_component * current.second[component];
    return result;
}

ActivePoint3 relative_displacement(const Quad8SurfaceContactLocalAdValues& state,
    const std::array<adlite::Scalar, 8>& secondary_shape,
    const std::array<adlite::Scalar, 8>& primary_shape) {
    ActivePoint3 result{};
    for (std::size_t component = 0; component < 3; ++component) {
        const std::size_t offset = 8 + 16 * component;
        for (std::size_t node = 0; node < 8; ++node) {
            result[component] += secondary_shape[node] * state[offset + node];
            result[component] -= primary_shape[node] * state[offset + 8 + node];
        }
    }
    return result;
}

ActivePoint3 incremental_relative_displacement(const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalAdValues& committed_state,
    const std::array<adlite::Scalar, 8>& secondary_shape,
    const std::array<adlite::Scalar, 8>& primary_shape) {
    const ActivePoint3 current = relative_displacement(state, secondary_shape, primary_shape),
                       committed = relative_displacement(committed_state, secondary_shape, primary_shape);
    return subtract(current, committed);
}

ActivePoint3 relative_position(const std::array<ActivePoint3, 16>& nodes,
    const std::array<adlite::Scalar, 8>& secondary_shape,
    const std::array<adlite::Scalar, 8>& primary_shape) {
    return subtract(interpolate_point(nodes, 0, secondary_shape), interpolate_point(nodes, 8, primary_shape));
}

ActivePoint3 objective_surface_increment(const std::array<ActivePoint3, 16>& nodes,
    const std::array<ActivePoint3, 16>& committed_nodes,
    const std::array<adlite::Scalar, 8>& secondary_shape,
    const std::array<adlite::Scalar, 8>& primary_shape,
    const SurfaceBasis& current_basis,
    const SurfaceBasis& committed_basis) {
    const ActivePoint3 current = relative_position(nodes, secondary_shape, primary_shape),
                       committed = relative_position(committed_nodes, secondary_shape, primary_shape);
    const adlite::Scalar first_increment = dot(current, current_basis.first) - dot(committed, committed_basis.first),
                         second_increment = dot(current, current_basis.second) - dot(committed, committed_basis.second);
    ActivePoint3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        result[component] =
            first_increment * current_basis.first[component] + second_increment * current_basis.second[component];
    return result;
}

SurfaceProjection8 small_sliding_projection(const std::array<ActivePoint3, 16>& nodes,
    const Quad8ToQuad8MechanicalGeometry& geometry) {
    const ActivePoint3 secondary_point =
                           interpolate_point(nodes, 0, active_values(geometry.secondary_displacement_shape)),
                       primary_point = interpolate_point(nodes, 8, active_values(geometry.primary_displacement_shape)),
                       primary_tangent_xi = interpolate_point(nodes, 8, active_values(geometry.primary_derivative_xi)),
                       primary_tangent_eta =
                           interpolate_point(nodes, 8, active_values(geometry.primary_derivative_eta)),
                       primary_area = cross(primary_tangent_xi, primary_tangent_eta);
    const adlite::Scalar primary_measure = norm(primary_area);
    if (!std::isfinite(primary_measure.value()) || !(primary_measure.value() > 0.0))
        throw std::domain_error("HEX20 small-sliding primary tangent plane has a nonpositive current measure");
    SurfaceProjection8 result;
    result.projected = true;
    result.primary_shape = active_values(geometry.primary_displacement_shape);
    result.primary_point = primary_point;
    result.tangent_xi = primary_tangent_xi;
    for (std::size_t component = 0; component < result.normal.size(); ++component)
        result.normal[component] = geometry.normal_orientation * primary_area[component] / primary_measure;
    result.gap = dot(subtract(primary_point, secondary_point), result.normal);
    return result;
}

SurfaceProjection8 finite_sliding_committed_projection(const std::array<ActivePoint3, 16>& nodes,
    const std::array<adlite::Scalar, 8>& secondary_shape,
    const adlite::Scalar& xi,
    const adlite::Scalar& eta,
    double normal_orientation) {
    Quad8ShapeValues values;
    quad8_shape(xi, eta, values);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape),
                       primary_point = interpolate_point(nodes, 8, values.shape),
                       tangent_xi = interpolate_point(nodes, 8, values.derivative_xi),
                       tangent_eta = interpolate_point(nodes, 8, values.derivative_eta),
                       area = cross(tangent_xi, tangent_eta);
    const adlite::Scalar measure = norm(area);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("HEX20 finite-sliding primary surface has a nonpositive committed measure");
    SurfaceProjection8 result;
    result.projected = true;
    result.primary_shape = values.shape;
    result.primary_point = primary_point;
    result.tangent_xi = tangent_xi;
    result.xi = xi;
    result.eta = eta;
    for (std::size_t component = 0; component < 3; ++component)
        result.normal[component] = normal_orientation * area[component] / measure;
    result.gap = dot(subtract(primary_point, secondary_point), result.normal);
    return result;
}

CartesianContactAdValue8 evaluate_surface_mechanical(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const std::array<adlite::Scalar, 8> secondary_shape = active_values(geometry.secondary_displacement_shape);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape);
    const SurfaceProjection8 projection = geometry.finite_sliding
                                              ? project_to_primary(secondary_point, nodes, geometry.normal_orientation)
                                              : small_sliding_projection(nodes, geometry);
    if (!projection.projected)
        return {};
    CartesianContactAdValue8 result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    result.pressure = adlite::max(-properties.penalty * result.gap, adlite::Scalar(0.0));
    result.tributary_area = geometry.quadrature_weight
                            * current_surface_measure(nodes,
                                active_values(geometry.secondary_derivative_xi),
                                active_values(geometry.secondary_derivative_eta),
                                0);
    result.contact_force = result.pressure * result.tributary_area;
    if (properties.friction_coefficient != 0.0) {
        const Quad8SurfaceContactLocalAdValues committed_ad_state = make_ad_state(committed_state, false);
        const std::array<ActivePoint3, 16> committed_nodes =
            current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, committed_ad_state);
        const SurfaceProjection8 committed_projection = geometry.finite_sliding
                                                            ? finite_sliding_committed_projection(committed_nodes,
                                                                  secondary_shape,
                                                                  projection.xi,
                                                                  projection.eta,
                                                                  geometry.normal_orientation)
                                                            : small_sliding_projection(committed_nodes, geometry);
        const SurfaceBasis current_basis = surface_basis(projection),
                           committed_coordinate_basis = surface_basis(committed_projection),
                           committed_contact_basis = stored_surface_basis(history, committed_coordinate_basis);
        result.tangent_first = current_basis.first;
        apply_friction(properties,
            history,
            transport_surface_vector(stored_history(history), current_basis, committed_contact_basis),
            transport_surface_vector(stored_total_history(history), current_basis, committed_contact_basis),
            objective_surface_increment(nodes,
                committed_nodes,
                secondary_shape,
                result.primary_shape,
                current_basis,
                committed_contact_basis),
            result);
    }
    return result;
}

CartesianContactAdValue8 evaluate_mechanical(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const ActivePoint3 secondary_point = nodes[geometry.secondary_local_node];
    const SurfaceProjection8 projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected)
        return {};
    CartesianContactAdValue8 result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    result.pressure = adlite::max(-properties.penalty * result.gap, adlite::Scalar(0.0));
    adlite::Scalar face_measure = 0.0, raw_sum = 0.0;
    std::array<adlite::Scalar, 8> raw{};
    for (std::size_t q = 0; q < quad8_surface_contact_quadrature_point_count; ++q) {
        const adlite::Scalar measure = current_surface_measure(nodes,
            active_values(geometry.secondary_derivatives_xi[q]),
            active_values(geometry.secondary_derivatives_eta[q]),
            0);
        const std::array<adlite::Scalar, 8> shape = active_values(geometry.secondary_shapes[q]);
        face_measure += geometry.secondary_quadrature_weights[q] * measure;
        for (std::size_t node = 0; node < 8; ++node)
            raw[node] += geometry.secondary_quadrature_weights[q] * measure
                         * (geometry.nodal_area_rule == Quad8NodalAreaRule::positive_lumped ? shape[node] * shape[node]
                                                                                            : shape[node]);
    }
    for (const adlite::Scalar& value : raw)
        raw_sum += value;
    if (!std::isfinite(face_measure.value()) || !(face_measure.value() > 0.0) || !std::isfinite(raw_sum.value())
        || !(raw_sum.value() > 0.0))
        throw std::domain_error("HEX20 contact Q8 face measure is nonpositive");
    if (geometry.nodal_area_rule == Quad8NodalAreaRule::positive_lumped
        && (!std::isfinite(raw[geometry.secondary_local_node].value())
            || !(raw[geometry.secondary_local_node].value() > 0.0)))
        throw std::domain_error("HEX20 contact positive-lumped Q8 nodal area is nonpositive");
    result.tributary_area = face_measure * raw[geometry.secondary_local_node] / raw_sum;
    result.contact_force = result.pressure * result.tributary_area;
    if (geometry.nodal_area_rule == Quad8NodalAreaRule::consistent_shape && !(result.tributary_area.value() > 0.0))
        return result;
    std::array<adlite::Scalar, 8> secondary_shape{};
    secondary_shape[geometry.secondary_local_node] = 1.0;
    const Quad8SurfaceContactLocalAdValues committed_ad_state = make_ad_state(committed_state, false);
    apply_friction(properties,
        history,
        stored_history(history),
        stored_total_history(history),
        incremental_relative_displacement(state, committed_ad_state, secondary_shape, result.primary_shape),
        result);
    return result;
}

Quad8SurfaceContactLocalAdValues make_ad_state(const Quad8SurfaceContactLocalValues& state, bool derivatives) {
    Quad8SurfaceContactLocalAdValues result{};
    if (derivatives)
        ad_local_system::make_active(state.data(), state.size(), result.data());
    else
        ad_local_system::make_passive(state.data(), state.size(), result.data());
    return result;
}

Quad8SurfaceContactLocalResidual extract(const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalAdValues& residual,
    Quad8SurfaceContactLocalJacobian* jacobian) {
    Quad8SurfaceContactLocalResidual values{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), values.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), values.data(), jacobian->data());
    return values;
}

// Circular averaging of extrapolated primary shapes. Coordinates are scaled
// by the current secondary-face radius before circle/polygon integration.
double quad8_disk_fraction_double(const Quad8ToQuad8HeatGeometry& geometry,
    const std::array<std::array<double, 3>, 16>& nodes) {
    std::array<double, 3> center{}, first{}, second{};
    for (std::size_t i = 0; i < 8; ++i)
        for (std::size_t c = 0; c < 3; ++c) {
            center[c] += geometry.secondary_displacement_shape[i] * nodes[i][c];
            first[c] += geometry.secondary_derivative_xi[i] * nodes[i][c];
            second[c] += geometry.secondary_derivative_eta[i] * nodes[i][c];
        }
    double first_squared = 0.0;
    for (std::size_t c = 0; c < 3; ++c)
        first_squared += first[c] * first[c];
    const double first_length = std::sqrt(first_squared);
    if (!(first_length > 0.0))
        throw std::domain_error("Disk transfer has a singular secondary tangent");
    for (auto& value : first)
        value /= first_length;
    double mixed = 0.0;
    for (std::size_t c = 0; c < 3; ++c)
        mixed += first[c] * second[c];
    double second_squared = 0.0;
    for (std::size_t c = 0; c < 3; ++c) {
        second[c] -= mixed * first[c];
        second_squared += second[c] * second[c];
    }
    const double second_length = std::sqrt(second_squared);
    if (!(second_length > 0.0))
        throw std::domain_error("Disk transfer has a singular secondary tangent plane");
    for (auto& value : second)
        value /= second_length;
    double minimum_edge_squared = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        double squared = 0.0;
        for (std::size_t c = 0; c < 3; ++c) {
            const auto difference = nodes[(i + 1) % 4][c] - nodes[i][c];
            squared += difference * difference;
        }
        if (i == 0 || squared < minimum_edge_squared)
            minimum_edge_squared = squared;
    }
    // H20 disconnected-face identification and independent nonmatching matrix.
    const double radius = 0.028276366456418445 * std::sqrt(minimum_edge_squared);
    if (!(radius > 0.0))
        throw std::domain_error("Disk transfer has a nonpositive secondary radius");
    constexpr std::array<std::size_t, 8> boundary{0, 4, 1, 5, 2, 6, 3, 7};
    std::array<std::array<double, 2>, 8> polygon{};
    bool vertex_inside = false, origin_inside = false, intersects = false;
    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            const auto delta = nodes[8 + boundary[i]][c] - center[c];
            polygon[i][0] += delta * first[c] / radius;
            polygon[i][1] += delta * second[c] / radius;
        }
        const auto squared = polygon[i][0] * polygon[i][0] + polygon[i][1] * polygon[i][1];
        vertex_inside = vertex_inside || squared < 1.0;
    }
    double area = 0.0;
    for (std::size_t i = 0; i < 8; ++i) {
        const auto& p = polygon[i];
        const auto& q = polygon[(i + 1) % 8];
        if ((p[1] > 0.0) != (q[1] > 0.0)) {
            const double crossing = p[0] - p[1] * (q[0] - p[0]) / (q[1] - p[1]);
            if (crossing > 0.0)
                origin_inside = !origin_inside;
        }
        const std::array<double, 2> d{q[0] - p[0], q[1] - p[1]};
        const double a = d[0] * d[0] + d[1] * d[1];
        if (a == 0.0)
            continue;
        const double b = 2.0 * (p[0] * d[0] + p[1] * d[1]), c = p[0] * p[0] + p[1] * p[1] - 1.0;
        const double discriminant = b * b - 4.0 * a * c;
        std::array<double, 4> cuts{0.0, 1.0, 0.0, 0.0};
        std::size_t count = 2;
        if (discriminant > 0.0) {
            const double root = std::sqrt(discriminant);
            for (const double sign : {-1.0, 1.0}) {
                const double t = (-b + sign * root) / (2.0 * a);
                if (t > 0.0 && t < 1.0) {
                    cuts[count++] = t;
                    intersects = true;
                }
            }
        }
        // At most four cuts; bounded insertion avoids GCC's small-array
        // std::sort bounds warning in independent non-LTO Release builds.
        for (std::size_t i = 1; i < cuts.size() && i < count; ++i)
            for (std::size_t j = i; j > 0 && cuts[j] < cuts[j - 1]; --j)
                std::swap(cuts[j], cuts[j - 1]);
        for (std::size_t part = 0; part + 1 < count; ++part) {
            const double middle = 0.5 * (cuts[part] + cuts[part + 1]);
            const std::array<double, 2> u{p[0] + cuts[part] * d[0], p[1] + cuts[part] * d[1]},
                v{p[0] + cuts[part + 1] * d[0], p[1] + cuts[part + 1] * d[1]},
                m{p[0] + middle * d[0], p[1] + middle * d[1]};
            const double product = u[0] * v[1] - u[1] * v[0];
            if ((m[0] * m[0] + m[1] * m[1]) < 1.0)
                area += 0.5 * product;
            else
                area += 0.5 * std::atan2(product, u[0] * v[0] + u[1] * v[1]);
        }
    }
    if (!intersects && !vertex_inside)
        return origin_inside ? 1.0 : 0.0;
    constexpr double pi = 3.141592653589793238462643383279502884;
    if (!std::isfinite(area))
        throw std::domain_error("Disk transfer has nonfinite intersection area");
    return (area < 0.0 ? -area : area) / pi;
}

// Circular averaging of extrapolated primary shapes. Coordinates are scaled
// by the current secondary-face radius before circle/polygon integration.
adlite::Scalar quad8_disk_fraction_ad(const Quad8ToQuad8HeatGeometry& geometry,
    const std::array<std::array<adlite::Scalar, 3>, 16>& nodes) {
    std::array<adlite::Scalar, 3> center{}, first{}, second{};
    for (std::size_t i = 0; i < 8; ++i)
        for (std::size_t c = 0; c < 3; ++c) {
            center[c] += geometry.secondary_displacement_shape[i] * nodes[i][c];
            first[c] += geometry.secondary_derivative_xi[i] * nodes[i][c];
            second[c] += geometry.secondary_derivative_eta[i] * nodes[i][c];
        }
    adlite::Scalar first_squared = 0.0;
    for (std::size_t c = 0; c < 3; ++c)
        first_squared += first[c] * first[c];
    const adlite::Scalar first_length = adlite::sqrt(first_squared);
    if (!(first_length.value() > 0.0))
        throw std::domain_error("Disk transfer has a singular secondary tangent");
    for (auto& value : first)
        value /= first_length;
    adlite::Scalar mixed = 0.0;
    for (std::size_t c = 0; c < 3; ++c)
        mixed += first[c] * second[c];
    adlite::Scalar second_squared = 0.0;
    for (std::size_t c = 0; c < 3; ++c) {
        second[c] -= mixed * first[c];
        second_squared += second[c] * second[c];
    }
    const adlite::Scalar second_length = adlite::sqrt(second_squared);
    if (!(second_length.value() > 0.0))
        throw std::domain_error("Disk transfer has a singular secondary tangent plane");
    for (auto& value : second)
        value /= second_length;
    adlite::Scalar minimum_edge_squared = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        adlite::Scalar squared = 0.0;
        for (std::size_t c = 0; c < 3; ++c) {
            const auto difference = nodes[(i + 1) % 4][c] - nodes[i][c];
            squared += difference * difference;
        }
        if (i == 0 || squared.value() < minimum_edge_squared.value())
            minimum_edge_squared = squared;
    }
    // H20 disconnected-face identification and independent nonmatching matrix.
    const adlite::Scalar radius = 0.028276366456418445 * adlite::sqrt(minimum_edge_squared);
    if (!(radius.value() > 0.0))
        throw std::domain_error("Disk transfer has a nonpositive secondary radius");
    constexpr std::array<std::size_t, 8> boundary{0, 4, 1, 5, 2, 6, 3, 7};
    std::array<std::array<adlite::Scalar, 2>, 8> polygon{};
    bool vertex_inside = false, origin_inside = false, intersects = false;
    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            const auto delta = nodes[8 + boundary[i]][c] - center[c];
            polygon[i][0] += delta * first[c] / radius;
            polygon[i][1] += delta * second[c] / radius;
        }
        const auto squared = polygon[i][0] * polygon[i][0] + polygon[i][1] * polygon[i][1];
        vertex_inside = vertex_inside || squared.value() < 1.0;
    }
    adlite::Scalar area = 0.0;
    for (std::size_t i = 0; i < 8; ++i) {
        const auto& p = polygon[i];
        const auto& q = polygon[(i + 1) % 8];
        if ((p[1].value() > 0.0) != (q[1].value() > 0.0)) {
            const double crossing =
                p[0].value() - p[1].value() * (q[0].value() - p[0].value()) / (q[1].value() - p[1].value());
            if (crossing > 0.0)
                origin_inside = !origin_inside;
        }
        const std::array<adlite::Scalar, 2> d{q[0] - p[0], q[1] - p[1]};
        const adlite::Scalar a = d[0] * d[0] + d[1] * d[1];
        if (a.value() == 0.0)
            continue;
        const adlite::Scalar b = 2.0 * (p[0] * d[0] + p[1] * d[1]), c = p[0] * p[0] + p[1] * p[1] - 1.0;
        const adlite::Scalar discriminant = b * b - 4.0 * a * c;
        std::array<adlite::Scalar, 4> cuts{0.0, 1.0, 0.0, 0.0};
        std::size_t count = 2;
        if (discriminant.value() > 0.0) {
            const adlite::Scalar root = adlite::sqrt(discriminant);
            for (const double sign : {-1.0, 1.0}) {
                const adlite::Scalar t = (-b + sign * root) / (2.0 * a);
                if (t.value() > 0.0 && t.value() < 1.0) {
                    cuts[count++] = t;
                    intersects = true;
                }
            }
        }
        for (std::size_t i = 1; i < cuts.size() && i < count; ++i)
            for (std::size_t j = i; j > 0 && cuts[j].value() < cuts[j - 1].value(); --j)
                std::swap(cuts[j], cuts[j - 1]);
        for (std::size_t part = 0; part + 1 < count; ++part) {
            const adlite::Scalar middle = 0.5 * (cuts[part] + cuts[part + 1]);
            const std::array<adlite::Scalar, 2> u{p[0] + cuts[part] * d[0], p[1] + cuts[part] * d[1]},
                v{p[0] + cuts[part + 1] * d[0], p[1] + cuts[part + 1] * d[1]},
                m{p[0] + middle * d[0], p[1] + middle * d[1]};
            const adlite::Scalar product = u[0] * v[1] - u[1] * v[0];
            if ((m[0] * m[0] + m[1] * m[1]).value() < 1.0)
                area += 0.5 * product;
            else
                area += 0.5 * adlite::atan2(product, u[0] * v[0] + u[1] * v[1]);
        }
    }
    if (!intersects && !vertex_inside)
        return origin_inside ? 1.0 : 0.0;
    constexpr double pi = 3.141592653589793238462643383279502884;
    if (!std::isfinite(area.value()))
        throw std::domain_error("Disk transfer has nonfinite intersection area");
    return (area.value() < 0.0 ? -area : area) / pi;
}

} // namespace

std::vector<double> compute_quad8_gap_heat_patch(const GapHeatProperties& properties,
    const std::vector<Quad8HeatPatchSample>& samples,
    const std::vector<double>& state,
    std::vector<double>* jacobian,
    CartesianHeatQuadratureValue* value) {
    if (samples.empty() || state.empty())
        throw std::invalid_argument("HEX20 thermal patch requires samples and a local state");
    const std::size_t size = state.size();
    // Moments are area, area*gap, area*T_secondary, area*T_primary, then
    // the area-weighted signed temperature test functions for each patch DOF.
    std::vector<double> moments(size + 4, 0.0);
    std::vector<double> derivatives(jacobian == nullptr ? 0 : (size + 4) * size, 0.0);
    std::vector<std::vector<std::size_t>> groups;
    std::map<std::size_t, std::size_t> disk_groups;
    std::size_t plain_group = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (!samples[i].disk_transfer) {
            if (plain_group == std::numeric_limits<std::size_t>::max()) {
                plain_group = groups.size();
                groups.push_back({});
            }
            groups[plain_group].push_back(i);
            continue;
        }
        const auto found = disk_groups.emplace(samples[i].transfer_group, groups.size());
        if (found.second)
            groups.push_back({});
        groups[found.first->second].push_back(i);
    }
    for (const auto& group : groups) {
        std::vector<double> group_moments(moments.size(), 0.0), group_derivatives(derivatives.size(), 0.0);
        std::vector<double> fraction_derivatives(jacobian == nullptr ? 0 : size, 0.0);
        const bool disk = samples[group.front()].disk_transfer;
        double fraction = disk ? 0.0 : 1.0;
        for (const std::size_t sample_index : group) {
            const Quad8HeatPatchSample& sample = samples[sample_index];
            Quad8SurfaceContactLocalValues local{};
            for (std::size_t column = 0; column < local.size(); ++column) {
                if (sample.local_dofs[column] >= size)
                    throw std::invalid_argument("HEX20 thermal patch sample has an invalid local DOF");
                local[column] = state[sample.local_dofs[column]];
            }
            const HeatAdValue8 point = evaluate_heat_geometry(sample.geometry,
                make_ad_state(local, jacobian != nullptr),
                sample.disk_transfer);
            if (!point.projected)
                throw std::domain_error("HEX20 thermal patch lost a required primary projection");
            if (sample.disk_transfer && point.transfer_fraction.value() == 0.0)
                continue;
            if (!std::isfinite(point.weighted_measure.value()) || !(point.weighted_measure.value() > 0.0))
                throw std::domain_error("HEX20 thermal patch requires positive finite sample measures");
            if (sample.disk_transfer) {
                fraction += point.transfer_fraction.value();
                if (jacobian != nullptr) {
                    std::array<double, quad8_surface_contact_local_dof_count> partials{};
                    point.transfer_fraction.copy_derivatives(partials.data(), partials.size());
                    for (std::size_t column = 0; column < local.size(); ++column)
                        fraction_derivatives[sample.local_dofs[column]] += partials[column];
                }
            }
            std::array<adlite::Scalar, 12> local_moments{};
            local_moments[0] = point.weighted_measure;
            local_moments[1] = point.weighted_measure * point.gap;
            local_moments[2] = point.weighted_measure * point.secondary_temperature;
            local_moments[3] = point.weighted_measure * point.primary_temperature;
            for (std::size_t node = 0; node < 4; ++node) {
                local_moments[4 + node] = point.weighted_measure * sample.geometry.secondary_temperature_shape[node];
                local_moments[8 + node] = -point.weighted_measure * point.primary_temperature_shape[node];
            }
            for (std::size_t row = 0; row < local_moments.size(); ++row) {
                const std::size_t target = row < 4 ? row : 4 + sample.local_dofs[row - 4];
                group_moments[target] += local_moments[row].value();
                if (jacobian != nullptr) {
                    std::array<double, quad8_surface_contact_local_dof_count> local_derivatives{};
                    local_moments[row].copy_derivatives(local_derivatives.data(), local_derivatives.size());
                    for (std::size_t column = 0; column < local.size(); ++column)
                        group_derivatives[target * size + sample.local_dofs[column]] += local_derivatives[column];
                }
            }
        }
        if (!(fraction > 0.0) || !std::isfinite(fraction))
            throw std::domain_error("HEX20 thermal disk sample has no primary support");
        for (std::size_t row = 0; row < moments.size(); ++row) {
            const double normalized = group_moments[row] / fraction;
            moments[row] += normalized;
            if (jacobian != nullptr)
                for (std::size_t column = 0; column < size; ++column)
                    derivatives[row * size + column] +=
                        (group_derivatives[row * size + column] - normalized * fraction_derivatives[column]) / fraction;
        }
    }
    const double area = moments[0];
    if (!std::isfinite(area) || !(area > 0.0))
        throw std::domain_error("HEX20 thermal patch has invalid total area");
    const std::array<double, 3> averages{moments[1] / area, moments[2] / area, moments[3] / area};
    std::array<adlite::Scalar, 3> active{};
    if (jacobian == nullptr)
        ad_local_system::make_passive(averages.data(), averages.size(), active.data());
    else
        ad_local_system::make_active(averages.data(), averages.size(), active.data());
    const adlite::Scalar flux =
        contact_common::gap_conductance(properties, active[0], active[1], active[2]) * (active[1] - active[2]);
    std::vector<double> result(size, 0.0);
    for (std::size_t row = 0; row < size; ++row)
        result[row] = moments[4 + row] * flux.value();
    if (value != nullptr)
        *value = {true, averages[0], flux.value(), area};
    if (jacobian != nullptr) {
        std::array<double, 3> flux_derivatives{};
        flux.copy_derivatives(flux_derivatives.data(), flux_derivatives.size());
        jacobian->assign(size * size, 0.0);
        for (std::size_t column = 0; column < size; ++column) {
            double flux_derivative = 0.0;
            for (std::size_t moment = 0; moment < 3; ++moment)
                flux_derivative +=
                    flux_derivatives[moment]
                    * (derivatives[(moment + 1) * size + column] - averages[moment] * derivatives[column]) / area;
            for (std::size_t row = 0; row < size; ++row)
                (*jacobian)[row * size + column] =
                    derivatives[(4 + row) * size + column] * flux.value() + moments[4 + row] * flux_derivative;
        }
    }
    return result;
}

Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_gap_heat(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    Quad8SurfaceContactLocalJacobian* jacobian) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad8SurfaceContactLocalAdValues residual{};
    const HeatAdValue8 value = evaluate_heat(properties, geometry, ad_state);
    if (value.projected) {
        for (std::size_t node = 0; node < 4; ++node) {
            residual[node] += value.weighted_measure * geometry.secondary_temperature_shape[node] * value.heat_flux;
            residual[4 + node] -= value.weighted_measure * value.primary_temperature_shape[node] * value.heat_flux;
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianHeatQuadratureValue compute_quad8_to_quad8_gap_heat_value(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state) {
    const HeatAdValue8 value = evaluate_heat(properties, geometry, make_ad_state(state, false));
    return {value.projected, value.gap.value(), value.heat_flux.value(), value.weighted_measure.value()};
}

ContactProjectionValue compute_quad8_to_quad8_heat_projection(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    Quad8ShapeValues shape;
    for (std::size_t node = 0; node < 8; ++node)
        shape.shape[node] = geometry.secondary_displacement_shape[node];
    const SurfaceProjection8 projection =
        project_to_primary(interpolate_point(nodes, 0, shape.shape), nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

std::array<adlite::Scalar, 8> compute_quad8_primary_shape_derivatives(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    bool allow_extrapolation) {
    const auto active = make_ad_state(state, false);
    const auto nodes = current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, active);
    const auto shape = active_values(geometry.secondary_displacement_shape);
    const auto projection =
        project_to_primary(interpolate_point(nodes, 0, shape), nodes, geometry.normal_orientation, allow_extrapolation);
    if (!projection.projected)
        throw std::domain_error("HEX20 finite-sliding shape derivative lost its primary projection");
    // Differentiate the converged two-coordinate projection equations, rather
    // than carrying all interface directions through every Newton iteration.
    std::array<CartesianPoint3, 8> primary{};
    for (std::size_t node = 0; node < 8; ++node)
        primary[node] = {nodes[8 + node][0].value(), nodes[8 + node][1].value(), nodes[8 + node][2].value()};
    const auto secondary = interpolate_point(nodes, 0, shape);
    DoubleQuad8ShapeValues values;
    double_quad8_shape(projection.xi.value(), projection.eta.value(), values);
    const CartesianPoint3 point = double_interpolate(primary, values.shape),
                          first = double_interpolate(primary, values.derivative_xi),
                          second = double_interpolate(primary, values.derivative_eta),
                          first_first = double_interpolate(primary, values.second_xi),
                          first_second = double_interpolate(primary, values.second_xi_eta),
                          second_second = double_interpolate(primary, values.second_eta),
                          difference{secondary[0].value() - point.x,
                              secondary[1].value() - point.y,
                              secondary[2].value() - point.z};
    const double a = -double_dot(first, first) + double_dot(difference, first_first),
                 b = -double_dot(first, second) + double_dot(difference, first_second),
                 c = -double_dot(second, second) + double_dot(difference, second_second), determinant = a * c - b * b;
    if (!std::isfinite(determinant) || std::abs(determinant) <= std::numeric_limits<double>::min())
        throw std::domain_error("HEX20 contact projection has a singular Q8 surface Jacobian");
    const std::array<double, 3> first_components{first.x, first.y, first.z},
        second_components{second.x, second.y, second.z},
        difference_components{difference.x, difference.y, difference.z};
    std::array<std::array<double, quad8_surface_contact_local_dof_count>, 8> derivatives{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 16; ++node) {
            const bool is_secondary = node < 8;
            const std::size_t local = is_secondary ? node : node - 8;
            const double coefficient =
                is_secondary ? geometry.secondary_displacement_shape[local] : -values.shape[local];
            const double rhs_first =
                coefficient * first_components[component]
                + (is_secondary ? 0.0 : values.derivative_xi[local] * difference_components[component]);
            const double rhs_second =
                coefficient * second_components[component]
                + (is_secondary ? 0.0 : values.derivative_eta[local] * difference_components[component]);
            const double dxi = projection.xi_clamped ? 0.0 : (-c * rhs_first + b * rhs_second) / determinant;
            const double deta = projection.eta_clamped ? 0.0 : (b * rhs_first - a * rhs_second) / determinant;
            for (std::size_t output = 0; output < 8; ++output)
                derivatives[output][8 + 16 * component + node] =
                    values.derivative_xi[output] * dxi + values.derivative_eta[output] * deta;
        }
    std::array<adlite::Scalar, 8> result{};
    for (std::size_t node = 0; node < 8; ++node)
        result[node] = adlite::Scalar::seeded(projection.primary_shape[node].value(),
            derivatives[node].data(),
            quad8_surface_contact_local_dof_count);
    return result;
}

std::array<double, 8> compute_quad8_disk_transfer(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    std::array<adlite::Scalar, 8>* derivatives) {
    std::array<std::array<double, 3>, 16> nodes{};
    std::array<CartesianPoint3, 8> secondary{}, primary{};
    for (std::size_t i = 0; i < 8; ++i) {
        const auto& s = geometry.secondary_coordinates[i];
        const auto& p = geometry.primary_coordinates[i];
        nodes[i] = {s.x + state[8 + i], s.y + state[24 + i], s.z + state[40 + i]};
        nodes[8 + i] = {p.x + state[16 + i], p.y + state[32 + i], p.z + state[48 + i]};
        secondary[i] = {nodes[i][0], nodes[i][1], nodes[i][2]};
        primary[i] = {nodes[8 + i][0], nodes[8 + i][1], nodes[8 + i][2]};
    }
    const double fraction = quad8_disk_fraction_double(geometry, nodes);
    std::array<double, 8> result{};
    if (fraction == 0.0) {
        if (derivatives)
            derivatives->fill(adlite::Scalar(0.0));
        return result;
    }
    const auto projection = compute_quad8_reference_projection(secondary,
        primary,
        geometry.secondary_displacement_shape,
        geometry.normal_orientation,
        true);
    if (!projection.projected)
        throw std::domain_error("Disk transfer lost its extrapolated primary surface");
    for (std::size_t i = 0; i < 8; ++i)
        result[i] = fraction * projection.primary_shape[i];
    if (derivatives) {
        const auto active = make_ad_state(state, true);
        const auto active_nodes = current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, active);
        const auto active_fraction = quad8_disk_fraction_ad(geometry, active_nodes);
        const auto shapes = compute_quad8_primary_shape_derivatives(geometry, state, true);
        for (std::size_t i = 0; i < 8; ++i)
            (*derivatives)[i] = adlite::compose(result[i], active_fraction * shapes[i], 1.0);
    }
    return result;
}

Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_contact(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad8SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue8 value =
        evaluate_surface_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        const std::array<adlite::Scalar, 8> secondary_shape = active_values(geometry.secondary_displacement_shape);
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 + 16 * component;
            const adlite::Scalar force = value.contact_force * value.normal[component]
                                         + value.tributary_area * value.tangential_traction_vector[component];
            for (std::size_t node = 0; node < 8; ++node) {
                residual[offset + node] += secondary_shape[node] * force;
                residual[offset + 8 + node] -= value.primary_shape[node] * force;
            }
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianContactPointValue compute_quad8_to_quad8_contact_value(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const CartesianContactAdValue8 value =
        evaluate_surface_mechanical(properties, geometry, make_ad_state(state, false), committed_state, history);
    return {value.projected,
        value.gap.value(),
        value.pressure.value(),
        value.tributary_area.value(),
        value.contact_force.value(),
        value.tangential_traction.value(),
        value.tangential_force.value(),
        value.friction_dissipation.value(),
        {value.normal[0].value(), value.normal[1].value(), value.normal[2].value()},
        {value.tangent_first[0].value(), value.tangent_first[1].value(), value.tangent_first[2].value()},
        {value.tangential_traction_vector[0].value(),
            value.tangential_traction_vector[1].value(),
            value.tangential_traction_vector[2].value()},
        {value.tangential_slip[0].value(), value.tangential_slip[1].value(), value.tangential_slip[2].value()},
        {value.elastic_tangential_slip[0].value(),
            value.elastic_tangential_slip[1].value(),
            value.elastic_tangential_slip[2].value()},
        value.sliding};
}

ContactProjectionValue compute_quad8_to_quad8_contact_projection(const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    if (!geometry.finite_sliding)
        return {true, small_sliding_projection(nodes, geometry).gap.value()};
    const ActivePoint3 secondary_point =
        interpolate_point(nodes, 0, active_values(geometry.secondary_displacement_shape));
    const SurfaceProjection8 projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

Quad8ReferenceProjectionValue compute_quad8_reference_projection(
    const std::array<CartesianPoint3, 8>& secondary_coordinates,
    const std::array<CartesianPoint3, 8>& primary_coordinates,
    const std::array<double, 8>& secondary_shape,
    double normal_orientation,
    bool allow_extrapolation) {
    const CartesianPoint3 secondary_point = double_interpolate(secondary_coordinates, secondary_shape);
    Quad8ReferenceProjectionValue result{};
    double xi = 0.0, eta = 0.0;
    for (std::size_t iteration = 0; iteration < 16; ++iteration) {
        DoubleQuad8ShapeValues values;
        double_quad8_shape(xi, eta, values);
        const CartesianPoint3 point = double_interpolate(primary_coordinates, values.shape),
                              tangent_xi = double_interpolate(primary_coordinates, values.derivative_xi),
                              tangent_eta = double_interpolate(primary_coordinates, values.derivative_eta),
                              tangent_xi_xi = double_interpolate(primary_coordinates, values.second_xi),
                              tangent_xi_eta = double_interpolate(primary_coordinates, values.second_xi_eta),
                              tangent_eta_eta = double_interpolate(primary_coordinates, values.second_eta),
                              difference = double_subtract(secondary_point, point);
        const double residual_xi = double_dot(difference, tangent_xi),
                     residual_eta = double_dot(difference, tangent_eta),
                     jacobian_xi_xi = -double_dot(tangent_xi, tangent_xi) + double_dot(difference, tangent_xi_xi),
                     jacobian_xi_eta = -double_dot(tangent_eta, tangent_xi) + double_dot(difference, tangent_xi_eta),
                     jacobian_eta_xi = -double_dot(tangent_xi, tangent_eta) + double_dot(difference, tangent_xi_eta),
                     jacobian_eta_eta = -double_dot(tangent_eta, tangent_eta) + double_dot(difference, tangent_eta_eta),
                     determinant = jacobian_xi_xi * jacobian_eta_eta - jacobian_xi_eta * jacobian_eta_xi;
        if (!std::isfinite(determinant) || std::abs(determinant) <= std::numeric_limits<double>::min())
            throw std::domain_error("HEX20 contact projection has a singular Q8 surface Jacobian");
        xi += (-residual_xi * jacobian_eta_eta + jacobian_xi_eta * residual_eta) / determinant;
        eta += (-jacobian_xi_xi * residual_eta + jacobian_eta_xi * residual_xi) / determinant;
    }
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(xi) || !std::isfinite(eta)
        || (!allow_extrapolation
            && (xi < -1.0 - tolerance || xi > 1.0 + tolerance || eta < -1.0 - tolerance || eta > 1.0 + tolerance)))
        return result;
    if (!allow_extrapolation) {
        xi = std::max(-1.0, std::min(1.0, xi));
        eta = std::max(-1.0, std::min(1.0, eta));
    }
    DoubleQuad8ShapeValues values;
    double_quad8_shape(xi, eta, values);
    const CartesianPoint3 primary_point = double_interpolate(primary_coordinates, values.shape),
                          tangent_xi = double_interpolate(primary_coordinates, values.derivative_xi),
                          tangent_eta = double_interpolate(primary_coordinates, values.derivative_eta),
                          area_vector = double_cross(tangent_xi, tangent_eta);
    const double measure = std::sqrt(double_dot(area_vector, area_vector));
    if (!std::isfinite(measure) || !(measure > 0.0))
        throw std::domain_error("HEX20 contact primary Q8 face has a nonpositive current measure");
    result.projected = true;
    result.primary_shape = values.shape;
    result.primary_derivative_xi = values.derivative_xi;
    result.primary_derivative_eta = values.derivative_eta;
    result.normal = {normal_orientation * area_vector.x / measure,
        normal_orientation * area_vector.y / measure,
        normal_orientation * area_vector.z / measure};
    result.gap = double_dot(double_subtract(primary_point, secondary_point), result.normal);
    return result;
}

Quad8SurfaceContactLocalResidual compute_node_to_quad8_contact(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad8SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue8 value =
        evaluate_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 + 16 * component;
            const adlite::Scalar force = value.contact_force * value.normal[component]
                                         + value.tributary_area * value.tangential_traction_vector[component];
            residual[offset + geometry.secondary_local_node] += force;
            for (std::size_t node = 0; node < 8; ++node)
                residual[offset + 8 + node] -= value.primary_shape[node] * force;
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianContactPointValue compute_node_to_quad8_contact_value(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const CartesianContactAdValue8 value =
        evaluate_mechanical(properties, geometry, make_ad_state(state, false), committed_state, history);
    return {value.projected,
        value.gap.value(),
        value.pressure.value(),
        value.tributary_area.value(),
        value.contact_force.value(),
        value.tangential_traction.value(),
        value.tangential_force.value(),
        value.friction_dissipation.value(),
        {value.normal[0].value(), value.normal[1].value(), value.normal[2].value()},
        {},
        {value.tangential_traction_vector[0].value(),
            value.tangential_traction_vector[1].value(),
            value.tangential_traction_vector[2].value()},
        {value.tangential_slip[0].value(), value.tangential_slip[1].value(), value.tangential_slip[2].value()},
        {value.elastic_tangential_slip[0].value(),
            value.elastic_tangential_slip[1].value(),
            value.elastic_tangential_slip[2].value()},
        value.sliding};
}

ContactProjectionValue compute_node_to_quad8_contact_projection(const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection8 projection =
        project_to_primary(nodes[geometry.secondary_local_node], nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

} // namespace fuelsim
