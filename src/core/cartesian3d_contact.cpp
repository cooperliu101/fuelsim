#include "detail/ad_local_system.hpp"
#include "fuelsim/core/contact.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim {
namespace {
using ActivePoint3 = std::array<adlite::Scalar, 3>;

struct SurfaceProjection final {
    bool projected = false;
    std::array<adlite::Scalar, 4> primary_shape{};
    ActivePoint3 primary_point{}, normal{};
    adlite::Scalar gap{0.0};
};

struct CartesianContactAdValue final {
    bool projected = false, sliding = false;
    std::array<adlite::Scalar, 4> primary_shape{};
    ActivePoint3 normal{}, elastic_tangential_slip{}, tangential_traction_vector{};
    adlite::Scalar gap{0.0}, pressure{0.0}, tributary_area{0.0}, contact_force{0.0}, tangential_traction{0.0},
        tangential_force{0.0};
};

ActivePoint3 subtract(const ActivePoint3& first, const ActivePoint3& second) {
    return {first[0] - second[0], first[1] - second[1], first[2] - second[2]};
}

adlite::Scalar dot(const ActivePoint3& first, const ActivePoint3& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

ActivePoint3 cross(const ActivePoint3& first, const ActivePoint3& second) {
    return {first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

adlite::Scalar norm(const ActivePoint3& value) { return adlite::hypot(adlite::hypot(value[0], value[1]), value[2]); }

std::array<ActivePoint3, 8> current_nodes(const std::array<CartesianPoint3, 4>& secondary,
    const std::array<CartesianPoint3, 4>& primary, const Quad4SurfaceContactLocalAdValues& state) {
    std::array<ActivePoint3, 8> result{};
    for (std::size_t node = 0; node < 8; ++node) {
        const CartesianPoint3& reference = node < 4 ? secondary[node] : primary[node - 4];
        result[node] = {reference.x + state[8 + node], reference.y + state[16 + node], reference.z + state[24 + node]};
    }
    return result;
}

ActivePoint3 interpolate_point(
    const std::array<ActivePoint3, 8>& nodes, std::size_t offset, const std::array<double, 4>& shape) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 4; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += shape[node] * nodes[offset + node][component];
    return result;
}

void quad4_shape(const adlite::Scalar& xi, const adlite::Scalar& eta, std::array<adlite::Scalar, 4>& shape,
    std::array<adlite::Scalar, 4>& derivative_xi, std::array<adlite::Scalar, 4>& derivative_eta) {
    shape = {0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)};
    derivative_xi = {-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)};
    derivative_eta = {-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)};
}

ActivePoint3 interpolate_primary(
    const std::array<ActivePoint3, 8>& nodes, const std::array<adlite::Scalar, 4>& coefficients) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 4; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += coefficients[node] * nodes[4 + node][component];
    return result;
}

SurfaceProjection project_to_primary(
    const ActivePoint3& secondary_point, const std::array<ActivePoint3, 8>& nodes, double normal_orientation) {
    adlite::Scalar xi = 0.0, eta = 0.0;
    constexpr std::array<double, 4> mixed_coefficients = {0.25, -0.25, 0.25, -0.25};
    for (std::size_t iteration = 0; iteration < 12; ++iteration) {
        std::array<adlite::Scalar, 4> shape{}, derivative_xi{}, derivative_eta{};
        quad4_shape(xi, eta, shape, derivative_xi, derivative_eta);
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
        const adlite::Scalar delta_xi =
                                 (-residual_xi * jacobian_eta_eta + jacobian_xi_eta * residual_eta) / determinant,
                             delta_eta = (-jacobian_xi_xi * residual_eta + jacobian_eta_xi * residual_xi) / determinant;
        xi += delta_xi;
        eta += delta_eta;
    }
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(xi.value()) || !std::isfinite(eta.value()) || xi.value() < -1.0 - tolerance ||
        xi.value() > 1.0 + tolerance || eta.value() < -1.0 - tolerance || eta.value() > 1.0 + tolerance)
        return {};
    if (xi.value() < -1.0) xi = -1.0;
    if (xi.value() > 1.0) xi = 1.0;
    if (eta.value() < -1.0) eta = -1.0;
    if (eta.value() > 1.0) eta = 1.0;
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
    for (std::size_t component = 0; component < 3; ++component)
        result.normal[component] = normal_orientation * area_vector[component] / measure;
    result.gap = dot(subtract(primary_point, secondary_point), result.normal);
    return result;
}

adlite::Scalar current_surface_measure(const std::array<ActivePoint3, 8>& nodes, std::size_t offset,
    const std::array<double, 4>& derivative_xi, const std::array<double, 4>& derivative_eta) {
    const ActivePoint3 tangent_xi = interpolate_point(nodes, offset, derivative_xi),
                       tangent_eta = interpolate_point(nodes, offset, derivative_eta);
    const adlite::Scalar measure = norm(cross(tangent_xi, tangent_eta));
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("Three-dimensional contact secondary face has a nonpositive current measure");
    return measure;
}

adlite::Scalar temperature(
    const Quad4SurfaceContactLocalAdValues& state, std::size_t offset, const std::array<double, 4>& shape) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 4; ++node) result += shape[node] * state[offset + node];
    return result;
}

struct HeatAdValue final {
    bool projected = false;
    std::array<adlite::Scalar, 4> primary_shape{};
    adlite::Scalar gap{0.0}, heat_flux{0.0}, weighted_measure{0.0};
};

HeatAdValue evaluate_heat(const GapHeatProperties& properties, const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalAdValues& state) {
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, geometry.secondary_shape);
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected) return {};
    const adlite::Scalar secondary_temperature = temperature(state, 0, geometry.secondary_shape);
    adlite::Scalar primary_temperature = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        primary_temperature += projection.primary_shape[node] * state[4 + node];
    const adlite::Scalar thermal_gap = adlite::max(projection.gap, adlite::Scalar(properties.minimum_gap)),
                         conductance = properties.gap_conductivity / thermal_gap,
                         heat_flux = conductance * (secondary_temperature - primary_temperature),
                         measure = current_surface_measure(
                             nodes, 0, geometry.secondary_derivative_xi, geometry.secondary_derivative_eta);
    return {true, projection.primary_shape, projection.gap, heat_flux, measure};
}

CartesianContactAdValue evaluate_mechanical(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalAdValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history) {
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const ActivePoint3 secondary_point = nodes[geometry.secondary_local_node];
    const SurfaceProjection projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected) return {};
    CartesianContactAdValue result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    const adlite::Scalar multiplier =
        properties.augmented_lagrangian ? adlite::Scalar(history.normal_multiplier) : adlite::Scalar(0.0);
    result.pressure = adlite::max(multiplier - properties.penalty * result.gap, adlite::Scalar(0.0));
    for (std::size_t q = 0; q < 4; ++q)
        result.tributary_area += geometry.secondary_shapes[q][geometry.secondary_local_node] *
                                 current_surface_measure(nodes, 0, geometry.secondary_derivatives_xi[q],
                                     geometry.secondary_derivatives_eta[q]);
    result.contact_force = result.pressure * result.tributary_area;
    if (properties.friction_coefficient == 0.0 || !(result.pressure.value() > 0.0)) return result;
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
        const adlite::Scalar history_component = history.cartesian_elastic_tangential_slip[component];
        result.elastic_tangential_slip[component] =
            history_component + relative_increment[component] - normal_increment * result.normal[component];
    }
    const adlite::Scalar history_normal = dot(result.elastic_tangential_slip, result.normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= history_normal * result.normal[component];
    ActivePoint3 trial_traction{};
    for (std::size_t component = 0; component < 3; ++component)
        trial_traction[component] = properties.penalty * result.elastic_tangential_slip[component];
    const adlite::Scalar trial_magnitude = norm(trial_traction),
                         sliding_limit = properties.friction_coefficient * result.pressure;
    if (trial_magnitude.value() < sliding_limit.value() ||
        (trial_magnitude.value() == sliding_limit.value() && !history.sliding)) {
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
    const Quad4SurfaceContactLocalAdValues& residual, Quad4SurfaceContactLocalJacobian* jacobian) {
    Quad4SurfaceContactLocalResidual values{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), values.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), values.data(), jacobian->data());
    return values;
}
} // namespace

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_gap_heat(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
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
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    const HeatAdValue value = evaluate_heat(properties, geometry, make_ad_state(state, false));
    return {value.projected, value.gap.value(), value.heat_flux.value(), value.weighted_measure.value()};
}

ContactProjectionValue compute_quad4_to_quad4_heat_projection(
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection projection =
        project_to_primary(interpolate_point(nodes, 0, geometry.secondary_shape), nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

Quad4SurfaceContactLocalResidual compute_node_to_quad4_contact(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad4SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue value = evaluate_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 * (component + 1);
            const adlite::Scalar force = value.contact_force * value.normal[component] +
                                         value.tributary_area * value.tangential_traction_vector[component];
            residual[offset + geometry.secondary_local_node] += force;
            for (std::size_t node = 0; node < 4; ++node)
                residual[offset + 4 + node] -= value.primary_shape[node] * force;
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianContactPointValue compute_node_to_quad4_contact_value(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history) {
    const CartesianContactAdValue value =
        evaluate_mechanical(properties, geometry, make_ad_state(state, false), committed_state, history);
    return {value.projected, value.gap.value(), value.pressure.value(), value.tributary_area.value(),
        value.contact_force.value(), value.tangential_traction.value(), value.tangential_force.value(),
        {value.normal[0].value(), value.normal[1].value(), value.normal[2].value()},
        {value.elastic_tangential_slip[0].value(), value.elastic_tangential_slip[1].value(),
            value.elastic_tangential_slip[2].value()},
        value.sliding};
}

ContactProjectionValue compute_node_to_quad4_contact_projection(
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    const Quad4SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 8> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection projection =
        project_to_primary(nodes[geometry.secondary_local_node], nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}
} // namespace fuelsim
