#include "quad4_face_boundary.hpp"
#include "boundary_types.hpp"
#include "detail/ad_local_system.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
constexpr double gauss = 0.577350269189625764509148780501957456;
}

Quad4FaceQuadraturePoint make_quad4_face_quadrature_point(const Quad4FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight) {
    Quad4FaceQuadraturePoint point{};
    point.shape = {{0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
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
    const CartesianPoint3 area_vector{point.tangent_xi.y * point.tangent_eta.z
                                          - point.tangent_xi.z * point.tangent_eta.y,
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
    const std::array<std::array<double, 2>, 4> nodal_area_locations = {{{{-nodal_area_location, -nodal_area_location}},
        {{nodal_area_location, -nodal_area_location}},
        {{nodal_area_location, nodal_area_location}},
        {{-nodal_area_location, nodal_area_location}}}};
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
        const CartesianPoint3 normal_area{normal_point.tangent_xi.y * normal_point.tangent_eta.z
                                              - normal_point.tangent_xi.z * normal_point.tangent_eta.y,
            normal_point.tangent_xi.z * normal_point.tangent_eta.x
                - normal_point.tangent_xi.x * normal_point.tangent_eta.z,
            normal_point.tangent_xi.x * normal_point.tangent_eta.y
                - normal_point.tangent_xi.y * normal_point.tangent_eta.x};
        const double normal_measure =
            std::sqrt(normal_area.x * normal_area.x + normal_area.y * normal_area.y + normal_area.z * normal_area.z);
        point.weighted_measure =
            (area.x * normal_area.x + area.y * normal_area.y + area.z * normal_area.z) / normal_measure;
        if (!std::isfinite(point.weighted_measure) || !(point.weighted_measure > 0.0))
            throw std::invalid_argument("Quad4 face thermal nodal integration requires a positive projected area");
    }
    return geometry;
}

Quad4FaceLocalResidual compute_quad4_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad4FaceGeometry& geometry,
    const Quad4FaceLocalValues& state,
    Quad4FaceLocalJacobian* jacobian) {
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
        std::array<adlite::Scalar, 3> normal_tangent_xi = {point.normal_tangent_xi.x,
            point.normal_tangent_xi.y,
            point.normal_tangent_xi.z};
        std::array<adlite::Scalar, 3> normal_tangent_eta = {point.normal_tangent_eta.x,
            point.normal_tangent_eta.y,
            point.normal_tangent_eta.z};
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
                const std::array<adlite::Scalar, 3> normal_area = {normal_tangent_xi[1] * normal_tangent_eta[2]
                                                                       - normal_tangent_xi[2] * normal_tangent_eta[1],
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
            for (std::size_t node = 0; node < 4; ++node)
                residual[node] -= data.load * measure * point.shape[node];
        } else {
            adlite::Scalar temperature = 0.0;
            for (std::size_t node = 0; node < 4; ++node)
                temperature += point.shape[node] * ad_state[node];
            const adlite::Scalar heat_flux = data.load * (temperature - data.ambient_temperature);
            for (std::size_t node = 0; node < 4; ++node)
                residual[node] += measure * point.shape[node] * heat_flux;
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
