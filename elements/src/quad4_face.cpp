#include "quad4_face.hpp"
#include "ad_local_system.hpp"
#include "contact_common.hpp"
#include "contact_types.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

namespace fuelsim {
namespace {
using ActivePoint3 = std::array<adlite::Scalar, 3>;

struct SurfaceProjection final {
    bool projected = false;
    std::array<adlite::Scalar, 4> primary_shape{};
    ActivePoint3 primary_point{}, normal{}, tangent_xi{};
    adlite::Scalar gap{0.0}, distance{0.0}, xi{0.0}, eta{0.0};
};

struct CartesianContactAdValue final {
    bool projected = false, sliding = false;
    std::array<adlite::Scalar, 4> primary_shape{};
    ActivePoint3 normal{}, tangent_first{}, tangential_slip{}, elastic_tangential_slip{}, tangential_traction_vector{};
    adlite::Scalar gap{0.0}, pressure{0.0}, tributary_area{0.0}, contact_force{0.0}, tangential_traction{0.0},
        tangential_force{0.0}, friction_dissipation{0.0};
};

struct SurfaceBasis final {
    ActivePoint3 first{}, second{};
};

struct DoubleSurfaceProjection final {
    bool projected = false;
    bool xi_constrained = false, eta_constrained = false;
    std::array<double, 4> primary_shape{};
    CartesianPoint3 primary_point{}, normal{}, tangent_xi{};
    double gap = 0.0, distance = 0.0, xi = 0.0, eta = 0.0;
};

CartesianPoint3 subtract_points(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

double dot_points(const CartesianPoint3& first, const CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

CartesianPoint3 cross_points(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

double norm_point(const CartesianPoint3& value) {
    return std::hypot(std::hypot(value.x, value.y), value.z);
}

CartesianPoint3 interpolate_points(const std::array<CartesianPoint3, 4>& nodes, const std::array<double, 4>& shape) {
    CartesianPoint3 result{};
    for (std::size_t node = 0; node < 4; ++node) {
        result.x += shape[node] * nodes[node].x;
        result.y += shape[node] * nodes[node].y;
        result.z += shape[node] * nodes[node].z;
    }
    return result;
}

void quad4_shape_values(double xi,
    double eta,
    std::array<double, 4>& shape,
    std::array<double, 4>& derivative_xi,
    std::array<double, 4>& derivative_eta) {
    shape = {0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)};
    derivative_xi = {-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)};
    derivative_eta = {-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)};
}

DoubleSurfaceProjection project_to_primary_double(const CartesianPoint3& secondary_point,
    const std::array<CartesianPoint3, 4>& primary_nodes,
    double normal_orientation,
    bool bounded_closest) {
    double xi = 0.0, eta = 0.0;
    constexpr std::array<double, 4> mixed_coefficients = {0.25, -0.25, 0.25, -0.25};
    const CartesianPoint3 mixed = interpolate_points(primary_nodes, mixed_coefficients);
    for (std::size_t iteration = 0; iteration < 12; ++iteration) {
        std::array<double, 4> shape{}, derivative_xi{}, derivative_eta{};
        quad4_shape_values(xi, eta, shape, derivative_xi, derivative_eta);
        const CartesianPoint3 point = interpolate_points(primary_nodes, shape),
                              tangent_xi = interpolate_points(primary_nodes, derivative_xi),
                              tangent_eta = interpolate_points(primary_nodes, derivative_eta),
                              difference = subtract_points(secondary_point, point);
        const double residual_xi = dot_points(difference, tangent_xi),
                     residual_eta = dot_points(difference, tangent_eta),
                     jacobian_xi_xi = -dot_points(tangent_xi, tangent_xi),
                     jacobian_eta_eta = -dot_points(tangent_eta, tangent_eta),
                     jacobian_xi_eta = -dot_points(tangent_eta, tangent_xi) + dot_points(difference, mixed),
                     jacobian_eta_xi = -dot_points(tangent_xi, tangent_eta) + dot_points(difference, mixed),
                     determinant = jacobian_xi_xi * jacobian_eta_eta - jacobian_xi_eta * jacobian_eta_xi;
        if (!std::isfinite(determinant) || std::abs(determinant) <= std::numeric_limits<double>::min())
            throw std::domain_error("Three-dimensional contact projection has a singular surface Jacobian");
        const double delta_xi = (-residual_xi * jacobian_eta_eta + jacobian_xi_eta * residual_eta) / determinant,
                     delta_eta = (-jacobian_xi_xi * residual_eta + jacobian_eta_xi * residual_xi) / determinant;
        xi += delta_xi;
        eta += delta_eta;
        const double coordinate_scale = std::max({1.0, std::abs(xi), std::abs(eta)});
        if (std::max(std::abs(delta_xi), std::abs(delta_eta))
            <= 64.0 * std::numeric_limits<double>::epsilon() * coordinate_scale)
            break;
    }
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(xi) || !std::isfinite(eta))
        return {};
    const auto projection_at =
        [&](double coordinate_xi, double coordinate_eta, bool xi_constrained = false, bool eta_constrained = false) {
            std::array<double, 4> shape{}, derivative_xi{}, derivative_eta{};
            quad4_shape_values(coordinate_xi, coordinate_eta, shape, derivative_xi, derivative_eta);
            const CartesianPoint3 primary_point = interpolate_points(primary_nodes, shape),
                                  tangent_xi = interpolate_points(primary_nodes, derivative_xi),
                                  tangent_eta = interpolate_points(primary_nodes, derivative_eta),
                                  area_vector = cross_points(tangent_xi, tangent_eta);
            const double measure = norm_point(area_vector);
            if (!std::isfinite(measure) || !(measure > 0.0))
                throw std::domain_error("Three-dimensional contact primary face has a nonpositive current measure");
            DoubleSurfaceProjection result;
            result.projected = true;
            result.xi_constrained = xi_constrained;
            result.eta_constrained = eta_constrained;
            result.primary_shape = shape;
            result.primary_point = primary_point;
            result.tangent_xi = tangent_xi;
            result.xi = coordinate_xi;
            result.eta = coordinate_eta;
            result.normal = {normal_orientation * area_vector.x / measure,
                normal_orientation * area_vector.y / measure,
                normal_orientation * area_vector.z / measure};
            const CartesianPoint3 separation = subtract_points(primary_point, secondary_point);
            result.gap = dot_points(separation, result.normal);
            result.distance = norm_point(separation);
            return result;
        };
    const bool interior =
        xi >= -1.0 - tolerance && xi <= 1.0 + tolerance && eta >= -1.0 - tolerance && eta <= 1.0 + tolerance;
    if (!bounded_closest) {
        if (!interior) {
            DoubleSurfaceProjection result;
            result.xi = xi;
            result.eta = eta;
            return result;
        }
        return projection_at(std::max(-1.0, std::min(1.0, xi)), std::max(-1.0, std::min(1.0, eta)));
    }
    DoubleSurfaceProjection closest;
    if (interior)
        closest = projection_at(xi, eta);
    const auto consider_edge =
        [&](std::size_t first_node, std::size_t second_node, bool varying_xi, double fixed_coordinate) {
            const CartesianPoint3 center = {0.5 * (primary_nodes[first_node].x + primary_nodes[second_node].x),
                                      0.5 * (primary_nodes[first_node].y + primary_nodes[second_node].y),
                                      0.5 * (primary_nodes[first_node].z + primary_nodes[second_node].z)},
                                  half_tangent = {0.5 * (primary_nodes[second_node].x - primary_nodes[first_node].x),
                                      0.5 * (primary_nodes[second_node].y - primary_nodes[first_node].y),
                                      0.5 * (primary_nodes[second_node].z - primary_nodes[first_node].z)};
            const double metric = dot_points(half_tangent, half_tangent);
            if (!std::isfinite(metric) || !(metric > 0.0))
                throw std::domain_error("Three-dimensional contact primary face has a degenerate current edge");
            double edge_coordinate = dot_points(subtract_points(secondary_point, center), half_tangent) / metric;
            const bool edge_constrained = edge_coordinate < -1.0 || edge_coordinate > 1.0;
            edge_coordinate = std::max(-1.0, std::min(1.0, edge_coordinate));
            const DoubleSurfaceProjection candidate =
                varying_xi ? projection_at(edge_coordinate, fixed_coordinate, edge_constrained, true)
                           : projection_at(fixed_coordinate, edge_coordinate, true, edge_constrained);
            if (!closest.projected || candidate.distance < closest.distance)
                closest = candidate;
        };
    consider_edge(0, 1, true, -1.0);
    consider_edge(1, 2, false, 1.0);
    consider_edge(3, 2, true, 1.0);
    consider_edge(0, 3, false, -1.0);
    return closest;
}

Quad4SurfaceContactLocalAdValues make_ad_state(const Quad4SurfaceContactLocalValues& state, bool derivatives);

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

std::array<ActivePoint3, 8> current_nodes(const std::array<CartesianPoint3, 4>& secondary,
    const std::array<CartesianPoint3, 4>& primary,
    const Quad4SurfaceContactLocalAdValues& state) {
    std::array<ActivePoint3, 8> result{};
    for (std::size_t node = 0; node < 8; ++node) {
        const CartesianPoint3& reference = node < 4 ? secondary[node] : primary[node - 4];
        result[node] = {reference.x + state[8 + node], reference.y + state[16 + node], reference.z + state[24 + node]};
    }
    return result;
}

ActivePoint3
interpolate_point(const std::array<ActivePoint3, 8>& nodes, std::size_t offset, const std::array<double, 4>& shape) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 4; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += shape[node] * nodes[offset + node][component];
    return result;
}

void quad4_shape(const adlite::Scalar& xi,
    const adlite::Scalar& eta,
    std::array<adlite::Scalar, 4>& shape,
    std::array<adlite::Scalar, 4>& derivative_xi,
    std::array<adlite::Scalar, 4>& derivative_eta) {
    shape = {0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)};
    derivative_xi = {-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)};
    derivative_eta = {-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)};
}

ActivePoint3 interpolate_primary(const std::array<ActivePoint3, 8>& nodes,
    const std::array<adlite::Scalar, 4>& coefficients) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 4; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += coefficients[node] * nodes[4 + node][component];
    return result;
}

SurfaceProjection project_to_primary(const ActivePoint3& secondary_point,
    const std::array<ActivePoint3, 8>& nodes,
    double normal_orientation,
    bool bounded_closest = false) {
    adlite::Scalar xi = 0.0, eta = 0.0;
    bool converged = false;
    // Translate the Newton geometry to one primary node so global translations
    // do not set the attainable projection increment through cancellation.
    std::array<ActivePoint3, 8> local_nodes{};
    for (std::size_t node = 4; node < 8; ++node)
        local_nodes[node] = subtract(nodes[node], nodes[4]);
    const ActivePoint3 local_secondary = subtract(secondary_point, nodes[4]);
    constexpr std::array<double, 4> mixed_coefficients = {0.25, -0.25, 0.25, -0.25};
    for (std::size_t iteration = 0; iteration < 12; ++iteration) {
        std::array<adlite::Scalar, 4> shape{}, derivative_xi{}, derivative_eta{};
        quad4_shape(xi, eta, shape, derivative_xi, derivative_eta);
        const ActivePoint3 point = interpolate_primary(local_nodes, shape),
                           tangent_xi = interpolate_primary(local_nodes, derivative_xi),
                           tangent_eta = interpolate_primary(local_nodes, derivative_eta),
                           mixed = interpolate_point(local_nodes, 4, mixed_coefficients),
                           difference = subtract(local_secondary, point);
        const adlite::Scalar residual_xi = dot(difference, tangent_xi), residual_eta = dot(difference, tangent_eta),
                             jacobian_xi_xi = -dot(tangent_xi, tangent_xi),
                             jacobian_eta_eta = -dot(tangent_eta, tangent_eta),
                             jacobian_xi_eta = -dot(tangent_eta, tangent_xi) + dot(difference, mixed),
                             jacobian_eta_xi = -dot(tangent_xi, tangent_eta) + dot(difference, mixed),
                             determinant = jacobian_xi_xi * jacobian_eta_eta - jacobian_xi_eta * jacobian_eta_xi;
        if (!std::isfinite(determinant.value()) || std::abs(determinant.value()) <= std::numeric_limits<double>::min())
            throw std::domain_error("Three-dimensional contact projection has a singular surface Jacobian");
        const adlite::Scalar delta_xi =
                                 (-residual_xi * jacobian_eta_eta + jacobian_xi_eta * residual_eta) / determinant,
                             delta_eta = (-jacobian_xi_xi * residual_eta + jacobian_eta_xi * residual_xi) / determinant;
        xi += delta_xi;
        eta += delta_eta;
        if (contact_common::projection_increment_converged(delta_xi,
                delta_eta,
                xi,
                eta,
                quad4_surface_contact_local_dof_count)) {
            converged = true;
            break;
        }
    }
    if (!converged)
        throw std::domain_error("Contact projection Newton iteration did not converge");
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(xi.value()) || !std::isfinite(eta.value()))
        return {};
    const auto projection_at = [&](const adlite::Scalar& coordinate_xi, const adlite::Scalar& coordinate_eta) {
        std::array<adlite::Scalar, 4> shape{}, derivative_xi{}, derivative_eta{};
        quad4_shape(coordinate_xi, coordinate_eta, shape, derivative_xi, derivative_eta);
        const ActivePoint3 primary_point = interpolate_primary(nodes, shape),
                           tangent_xi = interpolate_primary(nodes, derivative_xi),
                           tangent_eta = interpolate_primary(nodes, derivative_eta),
                           area_vector = cross(tangent_xi, tangent_eta);
        const adlite::Scalar measure = norm(area_vector);
        if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
            throw std::domain_error("Three-dimensional contact primary face has a nonpositive current measure");
        SurfaceProjection result;
        result.projected = true;
        result.primary_shape = shape;
        result.primary_point = primary_point;
        result.tangent_xi = tangent_xi;
        result.xi = coordinate_xi;
        result.eta = coordinate_eta;
        for (std::size_t component = 0; component < 3; ++component)
            result.normal[component] = normal_orientation * area_vector[component] / measure;
        result.gap = dot(subtract(primary_point, secondary_point), result.normal);
        result.distance = norm(subtract(primary_point, secondary_point));
        return result;
    };
    const bool interior = xi.value() >= -1.0 - tolerance && xi.value() <= 1.0 + tolerance
                          && eta.value() >= -1.0 - tolerance && eta.value() <= 1.0 + tolerance;
    if (!bounded_closest) {
        if (!interior) {
            SurfaceProjection result;
            result.xi = xi;
            result.eta = eta;
            return result;
        }
        if (xi.value() < -1.0)
            xi = -1.0;
        if (xi.value() > 1.0)
            xi = 1.0;
        if (eta.value() < -1.0)
            eta = -1.0;
        if (eta.value() > 1.0)
            eta = 1.0;
        return projection_at(xi, eta);
    }

    SurfaceProjection closest;
    if (interior)
        closest = projection_at(xi, eta);
    const auto consider_edge =
        [&](std::size_t first_node, std::size_t second_node, bool varying_xi, double fixed_coordinate) {
            ActivePoint3 center{}, half_tangent{};
            for (std::size_t component = 0; component < 3; ++component) {
                center[component] = 0.5 * (nodes[4 + first_node][component] + nodes[4 + second_node][component]);
                half_tangent[component] = 0.5 * (nodes[4 + second_node][component] - nodes[4 + first_node][component]);
            }
            const adlite::Scalar metric = dot(half_tangent, half_tangent);
            if (!std::isfinite(metric.value()) || !(metric.value() > 0.0))
                throw std::domain_error("Three-dimensional contact primary face has a degenerate current edge");
            adlite::Scalar edge_coordinate = dot(subtract(secondary_point, center), half_tangent) / metric;
            if (edge_coordinate.value() < -1.0)
                edge_coordinate = -1.0;
            if (edge_coordinate.value() > 1.0)
                edge_coordinate = 1.0;
            const SurfaceProjection candidate = varying_xi ? projection_at(edge_coordinate, fixed_coordinate)
                                                           : projection_at(fixed_coordinate, edge_coordinate);
            if (!closest.projected || candidate.distance.value() < closest.distance.value())
                closest = candidate;
        };
    consider_edge(0, 1, true, -1.0);
    consider_edge(1, 2, false, 1.0);
    consider_edge(3, 2, true, 1.0);
    consider_edge(0, 3, false, -1.0);
    return closest;
}

SurfaceProjection project_to_primary_from_double(const ActivePoint3& secondary_point,
    const std::array<ActivePoint3, 8>& nodes,
    double normal_orientation,
    const DoubleSurfaceProjection& double_projection) {
    if (!double_projection.projected)
        return {};
    adlite::Scalar xi = double_projection.xi, eta = double_projection.eta;
    if (!double_projection.xi_constrained && !double_projection.eta_constrained) {
        std::array<adlite::Scalar, 4> shape{}, derivative_xi{}, derivative_eta{};
        quad4_shape(xi, eta, shape, derivative_xi, derivative_eta);
        constexpr std::array<double, 4> mixed_coefficients = {0.25, -0.25, 0.25, -0.25};
        const ActivePoint3 point = interpolate_primary(nodes, shape),
                           tangent_xi = interpolate_primary(nodes, derivative_xi),
                           tangent_eta = interpolate_primary(nodes, derivative_eta),
                           mixed = interpolate_point(nodes, 4, mixed_coefficients),
                           difference = subtract(secondary_point, point);
        const adlite::Scalar residual_xi = dot(difference, tangent_xi), residual_eta = dot(difference, tangent_eta),
                             jacobian_xi_xi = -dot(tangent_xi, tangent_xi),
                             jacobian_eta_eta = -dot(tangent_eta, tangent_eta),
                             jacobian_xi_eta = -dot(tangent_eta, tangent_xi) + dot(difference, mixed),
                             jacobian_eta_xi = -dot(tangent_xi, tangent_eta) + dot(difference, mixed),
                             determinant = jacobian_xi_xi * jacobian_eta_eta - jacobian_xi_eta * jacobian_eta_xi;
        if (!std::isfinite(determinant.value()) || std::abs(determinant.value()) <= std::numeric_limits<double>::min())
            throw std::domain_error("Three-dimensional contact projection has a singular surface Jacobian");
        xi += (-residual_xi * jacobian_eta_eta + jacobian_xi_eta * residual_eta) / determinant;
        eta += (-jacobian_xi_xi * residual_eta + jacobian_eta_xi * residual_xi) / determinant;
    } else if (!double_projection.xi_constrained || !double_projection.eta_constrained) {
        const bool varying_xi = !double_projection.xi_constrained;
        const bool negative_fixed = varying_xi ? double_projection.eta < 0.0 : double_projection.xi < 0.0;
        const std::size_t first_node = varying_xi ? (negative_fixed ? 0 : 3) : (negative_fixed ? 0 : 1),
                          second_node = varying_xi ? (negative_fixed ? 1 : 2) : (negative_fixed ? 3 : 2);
        ActivePoint3 center{}, half_tangent{};
        for (std::size_t component = 0; component < 3; ++component) {
            center[component] = 0.5 * (nodes[4 + first_node][component] + nodes[4 + second_node][component]);
            half_tangent[component] = 0.5 * (nodes[4 + second_node][component] - nodes[4 + first_node][component]);
        }
        const adlite::Scalar metric = dot(half_tangent, half_tangent);
        if (!std::isfinite(metric.value()) || !(metric.value() > 0.0))
            throw std::domain_error("Three-dimensional contact primary face has a degenerate current edge");
        const adlite::Scalar coordinate = dot(subtract(secondary_point, center), half_tangent) / metric;
        if (varying_xi)
            xi = coordinate;
        else
            eta = coordinate;
    }

    std::array<adlite::Scalar, 4> shape{}, derivative_xi{}, derivative_eta{};
    quad4_shape(xi, eta, shape, derivative_xi, derivative_eta);
    const ActivePoint3 primary_point = interpolate_primary(nodes, shape),
                       tangent_xi = interpolate_primary(nodes, derivative_xi),
                       tangent_eta = interpolate_primary(nodes, derivative_eta),
                       area_vector = cross(tangent_xi, tangent_eta);
    const adlite::Scalar measure = norm(area_vector);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("Three-dimensional contact primary face has a nonpositive current measure");
    SurfaceProjection result;
    result.projected = true;
    result.primary_shape = shape;
    result.primary_point = primary_point;
    result.tangent_xi = tangent_xi;
    result.xi = xi;
    result.eta = eta;
    for (std::size_t component = 0; component < 3; ++component)
        result.normal[component] = normal_orientation * area_vector[component] / measure;
    result.gap = dot(subtract(primary_point, secondary_point), result.normal);
    result.distance = norm(subtract(primary_point, secondary_point));
    return result;
}

adlite::Scalar current_surface_measure(const std::array<ActivePoint3, 8>& nodes,
    std::size_t offset,
    const std::array<double, 4>& derivative_xi,
    const std::array<double, 4>& derivative_eta) {
    const ActivePoint3 tangent_xi = interpolate_point(nodes, offset, derivative_xi),
                       tangent_eta = interpolate_point(nodes, offset, derivative_eta);
    const adlite::Scalar measure = norm(cross(tangent_xi, tangent_eta));
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("Three-dimensional contact secondary face has a nonpositive current measure");
    return measure;
}

std::array<adlite::Scalar, 4> active_values(const std::array<double, 4>& values) {
    std::array<adlite::Scalar, 4> result{};
    for (std::size_t node = 0; node < 4; ++node)
        result[node] = values[node];
    return result;
}

ActivePoint3 interpolate_point(const std::array<ActivePoint3, 8>& nodes,
    std::size_t offset,
    const std::array<adlite::Scalar, 4>& shape) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 4; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += shape[node] * nodes[offset + node][component];
    return result;
}

SurfaceBasis surface_basis(const SurfaceProjection& projection, const ActivePoint3& normal) {
    ActivePoint3 tangent = projection.tangent_xi;
    const adlite::Scalar normal_component = dot(tangent, normal);
    for (std::size_t component = 0; component < 3; ++component)
        tangent[component] -= normal_component * normal[component];
    const adlite::Scalar measure = norm(tangent);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("Three-dimensional primary contact surface has an undefined convected tangent");
    SurfaceBasis result;
    for (std::size_t component = 0; component < 3; ++component)
        result.first[component] = tangent[component] / measure;
    result.second = cross(normal, result.first);
    return result;
}

ActivePoint3 secondary_average_normal(const std::array<ActivePoint3, 8>& nodes,
    const Quad4ToQuad4MechanicalGeometry& geometry) {
    const ActivePoint3 tangent_xi = interpolate_point(nodes, 0, geometry.secondary_normal_derivative_xi),
                       tangent_eta = interpolate_point(nodes, 0, geometry.secondary_normal_derivative_eta),
                       area = cross(tangent_xi, tangent_eta);
    const adlite::Scalar measure = norm(area);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("Three-dimensional secondary contact surface has an undefined averaged normal");
    ActivePoint3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        result[component] = geometry.secondary_normal_orientation * area[component] / measure;
    return result;
}

struct SurfacePullbackCovectors final {
    ActivePoint3 first, second;
};

SurfacePullbackCovectors surface_pullback_covectors(const std::array<ActivePoint3, 8>& current_nodes_value,
    const std::array<ActivePoint3, 8>& reference_nodes,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    double tangent_orientation,
    const ActivePoint3& current_normal,
    const ActivePoint3& current_first) {
    const ActivePoint3 current_xi = interpolate_point(current_nodes_value, 0, geometry.secondary_derivative_xi),
                       current_eta = interpolate_point(current_nodes_value, 0, geometry.secondary_derivative_eta),
                       reference_xi = interpolate_point(reference_nodes, 0, geometry.secondary_derivative_xi),
                       reference_eta = interpolate_point(reference_nodes, 0, geometry.secondary_derivative_eta),
                       current_second = cross(current_normal, current_first),
                       reference_area = cross(reference_xi, reference_eta);
    const adlite::Scalar reference_area_measure = norm(reference_area), reference_xi_measure = norm(reference_xi);
    if (!std::isfinite(reference_area_measure.value()) || !(reference_area_measure.value() > 0.0)
        || !std::isfinite(reference_xi_measure.value()) || !(reference_xi_measure.value() > 0.0))
        throw std::domain_error("Three-dimensional averaged friction has a degenerate reference tangent metric");
    ActivePoint3 reference_normal{}, reference_first{};
    for (std::size_t component = 0; component < 3; ++component) {
        reference_normal[component] =
            geometry.secondary_normal_orientation * reference_area[component] / reference_area_measure;
        reference_first[component] = tangent_orientation * reference_xi[component] / reference_xi_measure;
    }
    const ActivePoint3 reference_second = cross(reference_normal, reference_first);
    const adlite::Scalar current_11 = dot(current_xi, current_first), current_12 = dot(current_eta, current_first),
                         current_21 = dot(current_xi, current_second), current_22 = dot(current_eta, current_second),
                         reference_11 = dot(reference_xi, reference_first),
                         reference_12 = dot(reference_eta, reference_first),
                         reference_21 = dot(reference_xi, reference_second),
                         reference_22 = dot(reference_eta, reference_second),
                         current_determinant = current_11 * current_22 - current_12 * current_21;
    if (!std::isfinite(current_determinant.value()) || !(std::abs(current_determinant.value()) > 0.0))
        throw std::domain_error("Three-dimensional averaged friction has a singular current tangent metric");
    const adlite::Scalar first_first = (reference_11 * current_22 - reference_12 * current_21) / current_determinant,
                         first_second = (-reference_11 * current_12 + reference_12 * current_11) / current_determinant,
                         second_first = (reference_21 * current_22 - reference_22 * current_21) / current_determinant,
                         second_second = (-reference_21 * current_12 + reference_22 * current_11) / current_determinant;
    SurfacePullbackCovectors result{};
    for (std::size_t component = 0; component < 3; ++component) {
        result.first[component] = first_first * current_first[component] + first_second * current_second[component];
        result.second[component] = second_first * current_first[component] + second_second * current_second[component];
    }
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

ActivePoint3 relative_position(const std::array<ActivePoint3, 8>& nodes,
    const std::array<adlite::Scalar, 4>& secondary_shape,
    const std::array<adlite::Scalar, 4>& primary_shape) {
    return subtract(interpolate_point(nodes, 0, secondary_shape), interpolate_point(nodes, 4, primary_shape));
}

ActivePoint3 objective_surface_increment(const std::array<ActivePoint3, 8>& nodes,
    const std::array<ActivePoint3, 8>& committed_nodes,
    const std::array<adlite::Scalar, 4>& secondary_shape,
    const std::array<adlite::Scalar, 4>& primary_shape,
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

SurfaceProjection finite_sliding_committed_projection(const std::array<ActivePoint3, 8>& nodes,
    const std::array<adlite::Scalar, 4>& secondary_shape,
    const adlite::Scalar& xi,
    const adlite::Scalar& eta,
    double normal_orientation) {
    std::array<adlite::Scalar, 4> shape{}, derivative_xi{}, derivative_eta{};
    quad4_shape(xi, eta, shape, derivative_xi, derivative_eta);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape),
                       primary_point = interpolate_point(nodes, 4, shape),
                       tangent_xi = interpolate_point(nodes, 4, derivative_xi),
                       tangent_eta = interpolate_point(nodes, 4, derivative_eta), area = cross(tangent_xi, tangent_eta);
    const adlite::Scalar measure = norm(area);
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("HEX8 finite-sliding primary surface has a nonpositive committed measure");
    SurfaceProjection result;
    result.projected = true;
    result.primary_shape = shape;
    result.primary_point = primary_point;
    result.tangent_xi = tangent_xi;
    result.xi = xi;
    result.eta = eta;
    for (std::size_t component = 0; component < 3; ++component)
        result.normal[component] = normal_orientation * area[component] / measure;
    result.gap = dot(subtract(primary_point, secondary_point), result.normal);
    return result;
}

void apply_surface_friction(const NormalContactProperties& properties,
    const ContactPointHistory& history,
    const ActivePoint3& transported_history,
    const ActivePoint3& transported_total_history,
    const ActivePoint3& relative_increment,
    CartesianContactAdValue& result) {
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

adlite::Scalar
temperature(const Quad4SurfaceContactLocalAdValues& state, std::size_t offset, const std::array<double, 4>& shape) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        result += shape[node] * state[offset + node];
    return result;
}

struct HeatAdValue final {
    bool projected = false;
    std::array<adlite::Scalar, 4> primary_shape{};
    adlite::Scalar gap{0.0}, heat_flux{0.0}, weighted_measure{0.0};
};

HeatAdValue evaluate_heat(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalAdValues& state) {
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, geometry.secondary_shape);
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected)
        return {};
    const adlite::Scalar secondary_temperature = temperature(state, 0, geometry.secondary_shape);
    adlite::Scalar primary_temperature = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        primary_temperature += projection.primary_shape[node] * state[4 + node];
    const adlite::Scalar conductance = contact_common::gap_conductance(properties,
                             projection.gap,
                             secondary_temperature,
                             primary_temperature),
                         heat_flux = conductance * (secondary_temperature - primary_temperature),
                         measure = current_surface_measure(nodes,
                             0,
                             geometry.secondary_derivative_xi,
                             geometry.secondary_derivative_eta);
    return {true, projection.primary_shape, projection.gap, heat_flux, measure};
}

CartesianContactAdValue evaluate_surface_mechanical(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalAdValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const std::array<adlite::Scalar, 4> secondary_shape = active_values(geometry.secondary_shape);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape);
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected)
        return {};
    const ActivePoint3 normal = secondary_average_normal(nodes, geometry);
    CartesianContactAdValue result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = normal;
    result.gap = dot(subtract(projection.primary_point, secondary_point), result.normal);
    result.pressure = adlite::max(-properties.penalty * result.gap, adlite::Scalar(0.0));
    result.tributary_area =
        geometry.quadrature_weight
        * current_surface_measure(nodes, 0, geometry.secondary_derivative_xi, geometry.secondary_derivative_eta);
    result.contact_force = result.pressure * result.tributary_area;
    if (properties.friction_coefficient != 0.0) {
        const Quad4SurfaceContactLocalAdValues committed_ad_state = make_ad_state(committed_state, false);
        const std::array<ActivePoint3, 8> committed_nodes =
            current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, committed_ad_state);
        const SurfaceProjection committed_projection = finite_sliding_committed_projection(committed_nodes,
            secondary_shape,
            projection.xi,
            projection.eta,
            geometry.normal_orientation);
        const ActivePoint3 committed_normal = secondary_average_normal(committed_nodes, geometry);
        const SurfaceBasis current_basis = surface_basis(projection, normal),
                           committed_coordinate_basis = surface_basis(committed_projection, committed_normal),
                           committed_contact_basis = stored_surface_basis(history, committed_coordinate_basis);
        result.tangent_first = current_basis.first;
        apply_surface_friction(properties,
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

CartesianContactAdValue evaluate_mechanical(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalAdValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const ActivePoint3 secondary_point = nodes[geometry.secondary_local_node];
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected)
        return {};
    CartesianContactAdValue result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    const adlite::Scalar multiplier =
        properties.augmented_lagrangian ? adlite::Scalar(history.normal_multiplier) : adlite::Scalar(0.0);
    result.pressure = adlite::max(multiplier - properties.penalty * result.gap, adlite::Scalar(0.0));
    for (std::size_t q = 0; q < 4; ++q)
        result.tributary_area += geometry.secondary_shapes[q][geometry.secondary_local_node]
                                 * current_surface_measure(nodes,
                                     0,
                                     geometry.secondary_derivatives_xi[q],
                                     geometry.secondary_derivatives_eta[q]);
    result.contact_force = result.pressure * result.tributary_area;
    if (properties.friction_coefficient == 0.0)
        return result;
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] = history.cartesian_total_tangential_slip[component];
    const adlite::Scalar stored_total_normal = dot(result.tangential_slip, result.normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] -= stored_total_normal * result.normal[component];
    if (!(result.pressure.value() > 0.0))
        return result;
    ActivePoint3 relative_increment{};
    for (std::size_t component = 0; component < 3; ++component) {
        const std::size_t offset = 8 * (component + 1);
        relative_increment[component] =
            state[offset + geometry.secondary_local_node] - committed_state[offset + geometry.secondary_local_node];
        for (std::size_t node = 0; node < 4; ++node)
            relative_increment[component] -=
                projection.primary_shape[node] * (state[offset + 4 + node] - committed_state[offset + 4 + node]);
    }
    const adlite::Scalar normal_increment = dot(relative_increment, result.normal);
    for (std::size_t component = 0; component < 3; ++component) {
        const adlite::Scalar total_history_component = history.cartesian_total_tangential_slip[component];
        const adlite::Scalar history_component = history.cartesian_elastic_tangential_slip[component];
        result.tangential_slip[component] =
            total_history_component + relative_increment[component] - normal_increment * result.normal[component];
        result.elastic_tangential_slip[component] =
            history_component + relative_increment[component] - normal_increment * result.normal[component];
    }
    const adlite::Scalar total_history_normal = dot(result.tangential_slip, result.normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] -= total_history_normal * result.normal[component];
    const adlite::Scalar history_normal = dot(result.elastic_tangential_slip, result.normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= history_normal * result.normal[component];
    const ActivePoint3 trial_elastic_tangential_slip = result.elastic_tangential_slip;
    ActivePoint3 trial_traction{};
    for (std::size_t component = 0; component < 3; ++component)
        trial_traction[component] = properties.penalty * result.elastic_tangential_slip[component];
    const adlite::Scalar trial_magnitude = norm(trial_traction),
                         sliding_limit = properties.friction_coefficient * result.pressure;
    if (trial_magnitude.value() < sliding_limit.value()
        || (trial_magnitude.value() == sliding_limit.value() && !history.sliding)) {
        result.tangential_traction_vector = trial_traction;
        result.tangential_traction = trial_magnitude;
    } else {
        if (!(trial_magnitude.value() > 0.0))
            throw std::domain_error("Three-dimensional sliding contact has an undefined tangential direction");
        for (std::size_t component = 0; component < 3; ++component) {
            result.tangential_traction_vector[component] = sliding_limit * trial_traction[component] / trial_magnitude;
            result.elastic_tangential_slip[component] =
                result.tangential_traction_vector[component] / properties.penalty;
        }
        result.tangential_traction = sliding_limit;
        result.sliding = true;
        for (std::size_t component = 0; component < 3; ++component)
            result.friction_dissipation +=
                result.tangential_traction_vector[component]
                * (trial_elastic_tangential_slip[component] - result.elastic_tangential_slip[component]);
        result.friction_dissipation *= result.tributary_area;
    }
    result.tangential_force = result.tangential_traction * result.tributary_area;
    return result;
}

Quad4SurfaceContactLocalAdValues make_ad_state(const Quad4SurfaceContactLocalValues& state, bool derivatives) {
    Quad4SurfaceContactLocalAdValues result{};
    if (derivatives)
        ad_local_system::make_active(state.data(), state.size(), result.data());
    else
        ad_local_system::make_passive(state.data(), state.size(), result.data());
    return result;
}

Quad4SurfaceContactLocalResidual extract(const Quad4SurfaceContactLocalAdValues& state,
    const Quad4SurfaceContactLocalAdValues& residual,
    Quad4SurfaceContactLocalJacobian* jacobian) {
    Quad4SurfaceContactLocalResidual values{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), values.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), values.data(), jacobian->data());
    return values;
}
} // namespace

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_gap_heat(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    Quad4SurfaceContactLocalJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad4SurfaceContactLocalAdValues residual{};
    const HeatAdValue value = evaluate_heat(properties, geometry, ad_state);
    if (value.projected)
        for (std::size_t node = 0; node < 4; ++node) {
            residual[node] += value.weighted_measure * geometry.secondary_shape[node] * value.heat_flux;
            residual[4 + node] -= value.weighted_measure * value.primary_shape[node] * value.heat_flux;
        }
    return extract(ad_state, residual, jacobian);
}

CartesianHeatQuadratureValue compute_quad4_to_quad4_gap_heat_value(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state) {
    const HeatAdValue value = evaluate_heat(properties, geometry, make_ad_state(state, false));
    return {value.projected, value.gap.value(), value.heat_flux.value(), value.weighted_measure.value()};
}

ContactProjectionValue compute_quad4_to_quad4_heat_projection(const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection projection =
        project_to_primary(interpolate_point(nodes, 0, geometry.secondary_shape), nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

Quad4ReferenceProjectionValue compute_quad4_reference_projection(
    const std::array<CartesianPoint3, 4>& secondary_coordinates,
    const std::array<CartesianPoint3, 4>& primary_coordinates,
    const std::array<double, 4>& secondary_shape,
    double normal_orientation) {
    const DoubleSurfaceProjection projection =
        project_to_primary_double(interpolate_points(secondary_coordinates, secondary_shape),
            primary_coordinates,
            normal_orientation,
            false);
    Quad4ReferenceProjectionValue result{};
    result.projected = projection.projected;
    result.xi = projection.xi;
    result.eta = projection.eta;
    if (!projection.projected)
        return result;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    result.distance = projection.distance;
    return result;
}

Quad4ReferenceProjectionValue compute_quad4_reference_closest_projection(
    const std::array<CartesianPoint3, 4>& secondary_coordinates,
    const std::array<CartesianPoint3, 4>& primary_coordinates,
    const std::array<double, 4>& secondary_shape,
    double normal_orientation) {
    const DoubleSurfaceProjection projection =
        project_to_primary_double(interpolate_points(secondary_coordinates, secondary_shape),
            primary_coordinates,
            normal_orientation,
            true);
    Quad4ReferenceProjectionValue result{};
    result.projected = projection.projected;
    result.xi = projection.xi;
    result.eta = projection.eta;
    if (!projection.projected)
        return result;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    result.distance = projection.distance;
    return result;
}

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_contact(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad4SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue value =
        evaluate_surface_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        const std::array<adlite::Scalar, 4> secondary_shape = active_values(geometry.secondary_shape);
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 * (component + 1);
            const adlite::Scalar force = value.contact_force * value.normal[component]
                                         + value.tributary_area * value.tangential_traction_vector[component];
            for (std::size_t node = 0; node < 4; ++node) {
                residual[offset + node] += secondary_shape[node] * force;
                residual[offset + 4 + node] -= value.primary_shape[node] * force;
            }
        }
    }
    return extract(ad_state, residual, jacobian);
}

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_tangential_force_geometry(
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const std::array<double, 4>& secondary_distribution,
    double tangent_orientation,
    double first_traction,
    double second_traction,
    const Quad4SurfaceContactLocalValues& state,
    Quad4SurfaceContactLocalJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad4SurfaceContactLocalAdValues residual{};
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, geometry.secondary_shape);
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected)
        return extract(ad_state, residual, jacobian);
    const ActivePoint3 normal = secondary_average_normal(nodes, geometry),
                       tangent = interpolate_point(nodes, 0, geometry.secondary_derivative_xi);
    const adlite::Scalar tangent_measure = norm(tangent);
    if (!std::isfinite(tangent_measure.value()) || !(tangent_measure.value() > 0.0))
        throw std::domain_error("Three-dimensional averaged friction has an undefined current tangent");
    ActivePoint3 first{};
    for (std::size_t component = 0; component < 3; ++component)
        first[component] = tangent_orientation * tangent[component] / tangent_measure;
    const ActivePoint3 second = cross(normal, first);
    const adlite::Scalar area =
        geometry.quadrature_weight
        * current_surface_measure(nodes, 0, geometry.secondary_derivative_xi, geometry.secondary_derivative_eta);
    for (std::size_t component = 0; component < 3; ++component) {
        const std::size_t offset = 8 * (component + 1);
        const adlite::Scalar force = area * (first_traction * first[component] + second_traction * second[component]);
        for (std::size_t node = 0; node < 4; ++node) {
            residual[offset + node] += secondary_distribution[node] * force;
            residual[offset + 4 + node] -= projection.primary_shape[node] * force;
        }
    }
    return extract(ad_state, residual, jacobian);
}

Quad4AveragedFrictionGeometryValue compute_quad4_averaged_friction_geometry_value(
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const std::array<double, 4>& secondary_distribution,
    double tangent_orientation,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    Quad4AveragedFrictionGeometryJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr),
                                           committed_ad_state = make_ad_state(committed_state, false);
    const Quad4SurfaceContactLocalValues reference_state{};
    const std::array<ActivePoint3, 8> nodes = current_nodes(geometry.secondary_coordinates,
                                          geometry.primary_coordinates,
                                          ad_state),
                                      committed_nodes = current_nodes(geometry.secondary_coordinates,
                                          geometry.primary_coordinates,
                                          committed_ad_state),
                                      reference_nodes = current_nodes(geometry.secondary_coordinates,
                                          geometry.primary_coordinates,
                                          make_ad_state(reference_state, false));
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, geometry.secondary_shape);
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected)
        return {};
    const ActivePoint3 normal = secondary_average_normal(nodes, geometry),
                       tangent = interpolate_point(nodes, 0, geometry.secondary_derivative_xi);
    const adlite::Scalar tangent_measure = norm(tangent);
    if (!std::isfinite(tangent_measure.value()) || !(tangent_measure.value() > 0.0))
        throw std::domain_error("Three-dimensional averaged friction has an undefined current tangent");
    ActivePoint3 first{};
    for (std::size_t component = 0; component < 3; ++component)
        first[component] = tangent_orientation * tangent[component] / tangent_measure;
    const SurfacePullbackCovectors pullback =
        surface_pullback_covectors(nodes, reference_nodes, geometry, tangent_orientation, normal, first);
    const adlite::Scalar area =
        geometry.quadrature_weight
        * current_surface_measure(nodes, 0, geometry.secondary_derivative_xi, geometry.secondary_derivative_eta);
    const ActivePoint3 distributed_secondary = interpolate_point(nodes, 0, secondary_distribution),
                       committed_secondary = interpolate_point(committed_nodes, 0, secondary_distribution),
                       committed_primary = interpolate_point(committed_nodes, 4, projection.primary_shape),
                       separation = subtract(projection.primary_point, distributed_secondary),
                       increment = subtract(separation, subtract(committed_primary, committed_secondary));
    Quad4SurfaceContactLocalAdValues output{};
    output[0] = area;
    for (std::size_t component = 0; component < 3; ++component) {
        output[1 + component] = area * normal[component];
        output[4 + component] = area * first[component];
        output[7 + component] = area * separation[component];
    }
    output[10] = area * dot(increment, pullback.first);
    output[11] = area * dot(increment, pullback.second);
    Quad4SurfaceContactLocalJacobian local_jacobian{};
    const Quad4SurfaceContactLocalResidual values =
        extract(ad_state, output, jacobian == nullptr ? nullptr : &local_jacobian);
    if (jacobian != nullptr)
        for (std::size_t row = 0; row < 12; ++row)
            for (std::size_t column = 0; column < quad4_surface_contact_local_dof_count; ++column)
                (*jacobian)[row * quad4_surface_contact_local_dof_count + column] =
                    local_jacobian[row * quad4_surface_contact_local_dof_count + column];
    return {true,
        values[0],
        {values[1], values[2], values[3]},
        {values[4], values[5], values[6]},
        {values[7], values[8], values[9]},
        {values[10], values[11]}};
}

Quad4NormalForceAreaValue compute_quad4_to_quad4_normal_force_area(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    Quad4NormalForceAreaJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    const CartesianContactAdValue value = evaluate_surface_mechanical(properties, geometry, ad_state, {}, {});
    Quad4SurfaceContactLocalAdValues output{};
    output[0] = value.contact_force;
    output[1] = value.tributary_area;
    Quad4SurfaceContactLocalJacobian full_jacobian{};
    const Quad4SurfaceContactLocalResidual values =
        extract(ad_state, output, jacobian == nullptr ? nullptr : &full_jacobian);
    if (jacobian != nullptr)
        for (std::size_t row = 0; row < 2; ++row)
            for (std::size_t column = 0; column < quad4_surface_contact_local_dof_count; ++column)
                (*jacobian)[row * quad4_surface_contact_local_dof_count + column] =
                    full_jacobian[row * quad4_surface_contact_local_dof_count + column];
    return {value.projected, values[0], values[1]};
}

Quad4FiniteRegionNormalGeometryValue compute_quad4_finite_region_normal_geometry(
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    Quad4FiniteRegionNormalGeometryJacobian* jacobian) {
    std::array<CartesianPoint3, 4> secondary_nodes = geometry.secondary_coordinates,
                                   primary_nodes = geometry.primary_coordinates;
    for (std::size_t node = 0; node < 4; ++node) {
        secondary_nodes[node].x += state[8 + node];
        secondary_nodes[node].y += state[16 + node];
        secondary_nodes[node].z += state[24 + node];
        primary_nodes[node].x += state[12 + node];
        primary_nodes[node].y += state[20 + node];
        primary_nodes[node].z += state[28 + node];
    }
    const CartesianPoint3 secondary_point_value = interpolate_points(secondary_nodes, geometry.secondary_shape);
    const DoubleSurfaceProjection double_projection =
        project_to_primary_double(secondary_point_value, primary_nodes, geometry.normal_orientation, true);
    if (jacobian == nullptr) {
        const DoubleSurfaceProjection& projection = double_projection;
        if (!projection.projected)
            return {};
        const CartesianPoint3 normal_tangent_xi =
                                  interpolate_points(secondary_nodes, geometry.secondary_normal_derivative_xi),
                              normal_tangent_eta =
                                  interpolate_points(secondary_nodes, geometry.secondary_normal_derivative_eta),
                              normal_area_vector = cross_points(normal_tangent_xi, normal_tangent_eta),
                              area_tangent_xi = interpolate_points(secondary_nodes, geometry.secondary_derivative_xi),
                              area_tangent_eta = interpolate_points(secondary_nodes, geometry.secondary_derivative_eta),
                              area_vector = cross_points(area_tangent_xi, area_tangent_eta);
        const double normal_measure = norm_point(normal_area_vector), area_measure = norm_point(area_vector);
        if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
            throw std::domain_error("Three-dimensional secondary contact surface has an undefined averaged normal");
        if (!std::isfinite(area_measure) || !(area_measure > 0.0))
            throw std::domain_error("Three-dimensional contact secondary face has a nonpositive current measure");
        const CartesianPoint3 normal = {geometry.secondary_normal_orientation * normal_area_vector.x / normal_measure,
            geometry.secondary_normal_orientation * normal_area_vector.y / normal_measure,
            geometry.secondary_normal_orientation * normal_area_vector.z / normal_measure};
        Quad4FiniteRegionNormalGeometryValue result{};
        result.projected = true;
        result.area = geometry.quadrature_weight * area_measure;
        result.gap_integral =
            result.area * dot_points(subtract_points(projection.primary_point, secondary_point_value), normal);
        const std::array<double, 3> normal_components = {normal.x, normal.y, normal.z};
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 * (component + 1);
            const double force = result.area * normal_components[component];
            for (std::size_t node = 0; node < 4; ++node) {
                result.unit_pressure_residual[offset + node] += geometry.secondary_shape[node] * force;
                result.unit_pressure_residual[offset + 4 + node] -= projection.primary_shape[node] * force;
            }
        }
        return result;
    }
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const std::array<adlite::Scalar, 4> secondary_shape = active_values(geometry.secondary_shape);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape);
    const SurfaceProjection projection =
        project_to_primary_from_double(secondary_point, nodes, geometry.normal_orientation, double_projection);
    if (!projection.projected)
        return {};

    const ActivePoint3 normal = secondary_average_normal(nodes, geometry);
    const adlite::Scalar gap = dot(subtract(projection.primary_point, secondary_point), normal),
                         area = geometry.quadrature_weight
                                * current_surface_measure(nodes,
                                    0,
                                    geometry.secondary_derivative_xi,
                                    geometry.secondary_derivative_eta);
    std::array<adlite::Scalar, quad4_finite_region_normal_geometry_output_count> output{};
    output[0] = area;
    output[1] = area * gap;
    for (std::size_t component = 0; component < 3; ++component) {
        const std::size_t offset = 8 * (component + 1);
        const adlite::Scalar force = area * normal[component];
        for (std::size_t node = 0; node < 4; ++node) {
            output[2 + offset + node] += secondary_shape[node] * force;
            output[2 + offset + 4 + node] -= projection.primary_shape[node] * force;
        }
    }

    Quad4FiniteRegionNormalGeometryValue result{};
    result.projected = true;
    result.area = output[0].value();
    result.gap_integral = output[1].value();
    for (std::size_t row = 0; row < quad4_surface_contact_local_dof_count; ++row)
        result.unit_pressure_residual[row] = output[2 + row].value();
    if (jacobian != nullptr) {
        std::array<double, quad4_surface_contact_local_dof_count> derivatives{};
        for (std::size_t row = 0; row < output.size(); ++row) {
            output[row].copy_derivatives(derivatives.data(), derivatives.size());
            for (std::size_t column = 0; column < derivatives.size(); ++column)
                (*jacobian)[row * quad4_surface_contact_local_dof_count + column] = derivatives[column];
        }
    }
    return result;
}

CartesianContactPointValue compute_quad4_to_quad4_contact_value(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const CartesianContactAdValue value =
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

ContactProjectionValue compute_quad4_to_quad4_contact_projection(const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection projection =
        project_to_primary(interpolate_point(nodes, 0, active_values(geometry.secondary_shape)),
            nodes,
            geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

Quad4SurfaceContactLocalResidual compute_node_to_quad4_contact(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad4SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue value = evaluate_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 * (component + 1);
            const adlite::Scalar force = value.contact_force * value.normal[component]
                                         + value.tributary_area * value.tangential_traction_vector[component];
            residual[offset + geometry.secondary_local_node] += force;
            for (std::size_t node = 0; node < 4; ++node)
                residual[offset + 4 + node] -= value.primary_shape[node] * force;
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianContactPointValue compute_node_to_quad4_contact_value(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history) {
    const CartesianContactAdValue value =
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

ContactProjectionValue compute_node_to_quad4_contact_projection(const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection projection =
        project_to_primary(nodes[geometry.secondary_local_node], nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}
} // namespace fuelsim
