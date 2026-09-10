#include "quad8_face_boundary.hpp"
#include "boundary_types.hpp"
#include "detail/ad_local_system.hpp"
#include <cmath>
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
