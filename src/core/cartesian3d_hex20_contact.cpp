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

struct Quad8ShapeValues final {
    std::array<adlite::Scalar, 8> shape{}, derivative_xi{}, derivative_eta{};
    std::array<adlite::Scalar, 8> second_xi{}, second_xi_eta{}, second_eta{};
};

struct SurfaceProjection8 final {
    bool projected = false;
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

Quad8SurfaceContactLocalAdValues make_ad_state(const Quad8SurfaceContactLocalValues& state, bool derivatives);

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

void quad8_shape(const adlite::Scalar& xi, const adlite::Scalar& eta, Quad8ShapeValues& result) {
    const adlite::Scalar xm = 1.0 - xi, xp = 1.0 + xi, ym = 1.0 - eta, yp = 1.0 + eta;
    result.shape = {0.25 * xm * ym * (-xi - eta - 1.0), 0.25 * xp * ym * (xi - eta - 1.0),
        0.25 * xp * yp * (xi + eta - 1.0), 0.25 * xm * yp * (-xi + eta - 1.0), 0.5 * (1.0 - xi * xi) * ym,
        0.5 * xp * (1.0 - eta * eta), 0.5 * (1.0 - xi * xi) * yp, 0.5 * xm * (1.0 - eta * eta)};
    result.derivative_xi = {0.25 * ym * (2.0 * xi + eta), 0.25 * ym * (2.0 * xi - eta), 0.25 * yp * (2.0 * xi + eta),
        0.25 * yp * (2.0 * xi - eta), -xi * ym, 0.5 * (1.0 - eta * eta), -xi * yp, -0.5 * (1.0 - eta * eta)};
    result.derivative_eta = {0.25 * xm * (xi + 2.0 * eta), 0.25 * xp * (-xi + 2.0 * eta), 0.25 * xp * (xi + 2.0 * eta),
        0.25 * xm * (-xi + 2.0 * eta), -0.5 * (1.0 - xi * xi), -xp * eta, 0.5 * (1.0 - xi * xi), -xm * eta};
    result.second_xi = {0.5 * ym, 0.5 * ym, 0.5 * yp, 0.5 * yp, -ym, 0.0, -yp, 0.0};
    result.second_xi_eta = {0.25 * (1.0 - 2.0 * xi - 2.0 * eta), -0.25 * (1.0 + 2.0 * xi - 2.0 * eta),
        0.25 * (1.0 + 2.0 * xi + 2.0 * eta), 0.25 * (-1.0 + 2.0 * xi - 2.0 * eta), xi, -eta, -xi, eta};
    result.second_eta = {0.5 * xm, 0.5 * xp, 0.5 * xp, 0.5 * xm, 0.0, -xp, 0.0, -xm};
}

void quad4_temperature_shape(
    const adlite::Scalar& xi, const adlite::Scalar& eta, std::array<adlite::Scalar, 4>& shape) {
    shape = {0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta)};
}

std::array<ActivePoint3, 16> current_nodes(const std::array<CartesianPoint3, 8>& secondary,
    const std::array<CartesianPoint3, 8>& primary, const Quad8SurfaceContactLocalAdValues& state) {
    std::array<ActivePoint3, 16> result{};
    for (std::size_t node = 0; node < 16; ++node) {
        const CartesianPoint3& reference = node < 8 ? secondary[node] : primary[node - 8];
        result[node] = {reference.x + state[8 + node], reference.y + state[24 + node], reference.z + state[40 + node]};
    }
    return result;
}

ActivePoint3 interpolate_point(
    const std::array<ActivePoint3, 16>& nodes, std::size_t offset, const std::array<adlite::Scalar, 8>& shape) {
    ActivePoint3 result{};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t component = 0; component < 3; ++component)
            result[component] += shape[node] * nodes[offset + node][component];
    return result;
}

ActivePoint3 interpolate_primary(
    const std::array<ActivePoint3, 16>& nodes, const std::array<adlite::Scalar, 8>& coefficients) {
    return interpolate_point(nodes, 8, coefficients);
}

SurfaceProjection8 project_to_primary(
    const ActivePoint3& secondary_point, const std::array<ActivePoint3, 16>& nodes, double normal_orientation) {
    adlite::Scalar xi = 0.0, eta = 0.0;
    for (std::size_t iteration = 0; iteration < 16; ++iteration) {
        Quad8ShapeValues values;
        quad8_shape(xi, eta, values);
        const ActivePoint3 point = interpolate_primary(nodes, values.shape);
        const ActivePoint3 tangent_xi = interpolate_primary(nodes, values.derivative_xi);
        const ActivePoint3 tangent_eta = interpolate_primary(nodes, values.derivative_eta);
        const ActivePoint3 tangent_xi_xi = interpolate_primary(nodes, values.second_xi);
        const ActivePoint3 tangent_xi_eta = interpolate_primary(nodes, values.second_xi_eta);
        const ActivePoint3 tangent_eta_eta = interpolate_primary(nodes, values.second_eta);
        const ActivePoint3 difference = subtract(secondary_point, point);
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
    }
    constexpr double tolerance = 1.0e-10;
    if (!std::isfinite(xi.value()) || !std::isfinite(eta.value()) || xi.value() < -1.0 - tolerance ||
        xi.value() > 1.0 + tolerance || eta.value() < -1.0 - tolerance || eta.value() > 1.0 + tolerance)
        return {};
    if (xi.value() < -1.0) xi = -1.0;
    if (xi.value() > 1.0) xi = 1.0;
    if (eta.value() < -1.0) eta = -1.0;
    if (eta.value() > 1.0) eta = 1.0;
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
    const std::array<adlite::Scalar, 8>& derivative_xi, const std::array<adlite::Scalar, 8>& derivative_eta,
    std::size_t offset) {
    const ActivePoint3 tangent_xi = interpolate_point(nodes, offset, derivative_xi);
    const ActivePoint3 tangent_eta = interpolate_point(nodes, offset, derivative_eta);
    const adlite::Scalar measure = norm(cross(tangent_xi, tangent_eta));
    if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
        throw std::domain_error("HEX20 contact secondary Q8 face has a nonpositive current measure");
    return measure;
}

adlite::Scalar temperature(
    const Quad8SurfaceContactLocalAdValues& state, std::size_t offset, const std::array<adlite::Scalar, 4>& shape) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 4; ++node) result += shape[node] * state[offset + node];
    return result;
}

std::array<adlite::Scalar, 8> active_values(const std::array<double, 8>& values) {
    std::array<adlite::Scalar, 8> result{};
    for (std::size_t node = 0; node < 8; ++node) result[node] = values[node];
    return result;
}

std::array<adlite::Scalar, 4> active_temperature_values(const std::array<double, 4>& values) {
    std::array<adlite::Scalar, 4> result{};
    for (std::size_t node = 0; node < 4; ++node) result[node] = values[node];
    return result;
}

adlite::Scalar gap_conductance(const GapHeatProperties& properties, const adlite::Scalar& gap,
    const adlite::Scalar& secondary_temperature, const adlite::Scalar& primary_temperature) {
    if (properties.law == GapHeatConductanceLaw::gas_gap) {
        const adlite::Scalar thermal_gap = adlite::max(gap, adlite::Scalar(properties.minimum_gap));
        return properties.gap_conductivity / thermal_gap;
    }
    const adlite::Scalar pressure = adlite::max(-properties.contact_penalty * gap, adlite::Scalar(0.0));
    const adlite::Scalar average_temperature = 0.5 * (secondary_temperature + primary_temperature);
    const adlite::Scalar result =
        properties.conductance + properties.clearance_derivative * gap + properties.pressure_derivative * pressure +
        properties.temperature_derivative * (average_temperature - properties.reference_temperature);
    if (!std::isfinite(result.value()) || result.value() < 0.0)
        throw std::domain_error("HEX20 affine gap conductance must be finite and nonnegative");
    return result;
}

struct HeatAdValue8 final {
    bool projected = false;
    std::array<adlite::Scalar, 4> primary_temperature_shape{};
    adlite::Scalar gap{0.0}, heat_flux{0.0}, weighted_measure{0.0};
};

HeatAdValue8 evaluate_heat(const GapHeatProperties& properties, const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalAdValues& state) {
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    Quad8ShapeValues secondary_shape;
    secondary_shape.shape = {};
    for (std::size_t node = 0; node < 8; ++node)
        secondary_shape.shape[node] = geometry.secondary_displacement_shape[node];
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape.shape);
    const SurfaceProjection8 projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected) return {};
    std::array<adlite::Scalar, 4> primary_temperature_shape{};
    quad4_temperature_shape(projection.xi, projection.eta, primary_temperature_shape);
    adlite::Scalar primary_temperature = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        primary_temperature += primary_temperature_shape[node] * state[4 + node];
    const adlite::Scalar secondary_temperature =
        temperature(state, 0, active_temperature_values(geometry.secondary_temperature_shape));
    const adlite::Scalar heat_flux =
        gap_conductance(properties, projection.gap, secondary_temperature, primary_temperature) *
        (secondary_temperature - primary_temperature);
    const adlite::Scalar measure = current_surface_measure(
        nodes, active_values(geometry.secondary_derivative_xi), active_values(geometry.secondary_derivative_eta), 0);
    return {true, primary_temperature_shape, projection.gap, heat_flux, measure * geometry.quadrature_weight};
}

CartesianContactAdValue8 evaluate_mechanical(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history);

void apply_friction(const NormalContactProperties& properties, const ContactPointHistory& history,
    const ActivePoint3& transported_history, const ActivePoint3& transported_total_history,
    const ActivePoint3& relative_increment, CartesianContactAdValue8& result) {
    if (properties.friction_coefficient == 0.0 || !(result.pressure.value() > 0.0)) return;
    const adlite::Scalar normal_increment = dot(relative_increment, result.normal);
    for (std::size_t component = 0; component < 3; ++component) {
        result.tangential_slip[component] = transported_total_history[component] + relative_increment[component] -
                                            normal_increment * result.normal[component];
        result.elastic_tangential_slip[component] = transported_history[component] + relative_increment[component] -
                                                    normal_increment * result.normal[component];
    }
    const adlite::Scalar total_history_normal = dot(result.tangential_slip, result.normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] -= total_history_normal * result.normal[component];
    const adlite::Scalar history_normal = dot(result.elastic_tangential_slip, result.normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= history_normal * result.normal[component];
    const ActivePoint3 trial_elastic_tangential_slip = result.elastic_tangential_slip;
    const adlite::Scalar sliding_limit = properties.friction_coefficient * result.pressure,
                         stick_stiffness = properties.maximum_elastic_slip > 0.0
                                               ? sliding_limit / properties.maximum_elastic_slip
                                               : adlite::Scalar(properties.penalty);
    ActivePoint3 trial_traction{};
    for (std::size_t component = 0; component < 3; ++component)
        trial_traction[component] = stick_stiffness * result.elastic_tangential_slip[component];
    const adlite::Scalar trial_magnitude = norm(trial_traction);
    if (trial_magnitude.value() < sliding_limit.value() ||
        (trial_magnitude.value() == sliding_limit.value() && !history.sliding)) {
        result.tangential_traction_vector = trial_traction;
        result.tangential_traction = trial_magnitude;
    } else {
        if (!(trial_magnitude.value() > 0.0))
            throw std::domain_error("HEX20 sliding contact has an undefined tangential direction");
        for (std::size_t component = 0; component < 3; ++component) {
            result.tangential_traction_vector[component] = sliding_limit * trial_traction[component] / trial_magnitude;
            result.elastic_tangential_slip[component] = result.tangential_traction_vector[component] / stick_stiffness;
        }
        result.tangential_traction = sliding_limit;
        result.sliding = true;
        for (std::size_t component = 0; component < 3; ++component)
            result.friction_dissipation +=
                result.tangential_traction_vector[component] *
                (trial_elastic_tangential_slip[component] - result.elastic_tangential_slip[component]);
        result.friction_dissipation *= result.tributary_area;
    }
    result.tangential_force = result.tangential_traction * result.tributary_area;
}

ActivePoint3 stored_history(const ContactPointHistory& history) {
    return {history.cartesian_elastic_tangential_slip[0], history.cartesian_elastic_tangential_slip[1],
        history.cartesian_elastic_tangential_slip[2]};
}

ActivePoint3 stored_total_history(const ContactPointHistory& history) {
    return {history.cartesian_total_tangential_slip[0], history.cartesian_total_tangential_slip[1],
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
    if (!history.cartesian_tangent_basis_initialized) return fallback;
    SurfaceBasis result;
    for (std::size_t component = 0; component < 3; ++component) {
        result.first[component] = history.cartesian_contact_tangent_first[component];
        result.second[component] = history.cartesian_contact_normal[component];
    }
    result.second = cross(result.second, result.first);
    return result;
}

ActivePoint3 transport_surface_vector(
    const ActivePoint3& vector, const SurfaceBasis& current, const SurfaceBasis& committed) {
    const adlite::Scalar first_component = dot(vector, committed.first),
                         second_component = dot(vector, committed.second);
    ActivePoint3 result{};
    for (std::size_t component = 0; component < 3; ++component)
        result[component] = first_component * current.first[component] + second_component * current.second[component];
    return result;
}

ActivePoint3 relative_displacement(const Quad8SurfaceContactLocalAdValues& state,
    const std::array<adlite::Scalar, 8>& secondary_shape, const std::array<adlite::Scalar, 8>& primary_shape) {
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
    const Quad8SurfaceContactLocalAdValues& committed_state, const std::array<adlite::Scalar, 8>& secondary_shape,
    const std::array<adlite::Scalar, 8>& primary_shape) {
    const ActivePoint3 current = relative_displacement(state, secondary_shape, primary_shape),
                       committed = relative_displacement(committed_state, secondary_shape, primary_shape);
    return subtract(current, committed);
}

ActivePoint3 relative_position(const std::array<ActivePoint3, 16>& nodes,
    const std::array<adlite::Scalar, 8>& secondary_shape, const std::array<adlite::Scalar, 8>& primary_shape) {
    return subtract(interpolate_point(nodes, 0, secondary_shape), interpolate_point(nodes, 8, primary_shape));
}

ActivePoint3 objective_surface_increment(const std::array<ActivePoint3, 16>& nodes,
    const std::array<ActivePoint3, 16>& committed_nodes, const std::array<adlite::Scalar, 8>& secondary_shape,
    const std::array<adlite::Scalar, 8>& primary_shape, const SurfaceBasis& current_basis,
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

SurfaceProjection8 small_sliding_projection(
    const std::array<ActivePoint3, 16>& nodes, const Quad8ToQuad8MechanicalGeometry& geometry) {
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
    const std::array<adlite::Scalar, 8>& secondary_shape, const adlite::Scalar& xi, const adlite::Scalar& eta,
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
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history) {
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const std::array<adlite::Scalar, 8> secondary_shape = active_values(geometry.secondary_displacement_shape);
    const ActivePoint3 secondary_point = interpolate_point(nodes, 0, secondary_shape);
    const SurfaceProjection8 projection = geometry.finite_sliding
                                              ? project_to_primary(secondary_point, nodes, geometry.normal_orientation)
                                              : small_sliding_projection(nodes, geometry);
    if (!projection.projected) return {};
    CartesianContactAdValue8 result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    result.pressure = adlite::max(-properties.penalty * result.gap, adlite::Scalar(0.0));
    result.tributary_area =
        geometry.quadrature_weight * current_surface_measure(nodes, active_values(geometry.secondary_derivative_xi),
                                         active_values(geometry.secondary_derivative_eta), 0);
    result.contact_force = result.pressure * result.tributary_area;
    if (properties.friction_coefficient != 0.0 && result.pressure.value() > 0.0) {
        const Quad8SurfaceContactLocalAdValues committed_ad_state = make_ad_state(committed_state, false);
        const std::array<ActivePoint3, 16> committed_nodes =
            current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, committed_ad_state);
        const SurfaceProjection8 committed_projection =
            geometry.finite_sliding ? finite_sliding_committed_projection(committed_nodes, secondary_shape,
                                          projection.xi, projection.eta, geometry.normal_orientation)
                                    : small_sliding_projection(committed_nodes, geometry);
        const SurfaceBasis current_basis = surface_basis(projection),
                           committed_coordinate_basis = surface_basis(committed_projection),
                           committed_contact_basis = stored_surface_basis(history, committed_coordinate_basis);
        result.tangent_first = current_basis.first;
        apply_friction(properties, history,
            transport_surface_vector(stored_history(history), current_basis, committed_contact_basis),
            transport_surface_vector(stored_total_history(history), current_basis, committed_contact_basis),
            objective_surface_increment(
                nodes, committed_nodes, secondary_shape, result.primary_shape, current_basis, committed_contact_basis),
            result);
    }
    return result;
}

CartesianContactAdValue8 evaluate_mechanical(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalAdValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history) {
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, state);
    const ActivePoint3 secondary_point = nodes[geometry.secondary_local_node];
    const SurfaceProjection8 projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    if (!projection.projected) return {};
    CartesianContactAdValue8 result;
    result.projected = true;
    result.primary_shape = projection.primary_shape;
    result.normal = projection.normal;
    result.gap = projection.gap;
    result.pressure = adlite::max(-properties.penalty * result.gap, adlite::Scalar(0.0));
    adlite::Scalar face_measure = 0.0, raw_sum = 0.0;
    std::array<adlite::Scalar, 8> raw{};
    for (std::size_t q = 0; q < quad8_surface_contact_quadrature_point_count; ++q) {
        const adlite::Scalar measure =
            current_surface_measure(nodes, active_values(geometry.secondary_derivatives_xi[q]),
                active_values(geometry.secondary_derivatives_eta[q]), 0);
        const std::array<adlite::Scalar, 8> shape = active_values(geometry.secondary_shapes[q]);
        face_measure += geometry.secondary_quadrature_weights[q] * measure;
        for (std::size_t node = 0; node < 8; ++node)
            raw[node] += geometry.secondary_quadrature_weights[q] * measure *
                         (geometry.nodal_area_rule == Quad8NodalAreaRule::positive_lumped ? shape[node] * shape[node]
                                                                                          : shape[node]);
    }
    for (const adlite::Scalar& value : raw) raw_sum += value;
    if (!std::isfinite(face_measure.value()) || !(face_measure.value() > 0.0) || !std::isfinite(raw_sum.value()) ||
        !(raw_sum.value() > 0.0))
        throw std::domain_error("HEX20 contact Q8 face measure is nonpositive");
    if (geometry.nodal_area_rule == Quad8NodalAreaRule::positive_lumped &&
        (!std::isfinite(raw[geometry.secondary_local_node].value()) ||
            !(raw[geometry.secondary_local_node].value() > 0.0)))
        throw std::domain_error("HEX20 contact positive-lumped Q8 nodal area is nonpositive");
    result.tributary_area = face_measure * raw[geometry.secondary_local_node] / raw_sum;
    result.contact_force = result.pressure * result.tributary_area;
    if (geometry.nodal_area_rule == Quad8NodalAreaRule::consistent_shape && !(result.tributary_area.value() > 0.0))
        return result;
    std::array<adlite::Scalar, 8> secondary_shape{};
    secondary_shape[geometry.secondary_local_node] = 1.0;
    const Quad8SurfaceContactLocalAdValues committed_ad_state = make_ad_state(committed_state, false);
    apply_friction(properties, history, stored_history(history), stored_total_history(history),
        incremental_relative_displacement(state, committed_ad_state, secondary_shape, result.primary_shape), result);
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
    const Quad8SurfaceContactLocalAdValues& residual, Quad8SurfaceContactLocalJacobian* jacobian) {
    Quad8SurfaceContactLocalResidual values{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), values.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), values.data(), jacobian->data());
    return values;
}
} // namespace

Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_gap_heat(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
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
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    const HeatAdValue8 value = evaluate_heat(properties, geometry, make_ad_state(state, false));
    return {value.projected, value.gap.value(), value.heat_flux.value(), value.weighted_measure.value()};
}

ContactProjectionValue compute_quad8_to_quad8_heat_projection(
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    Quad8ShapeValues shape;
    for (std::size_t node = 0; node < 8; ++node) shape.shape[node] = geometry.secondary_displacement_shape[node];
    const SurfaceProjection8 projection =
        project_to_primary(interpolate_point(nodes, 0, shape.shape), nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_contact(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad8SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue8 value =
        evaluate_surface_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        const std::array<adlite::Scalar, 8> secondary_shape = active_values(geometry.secondary_displacement_shape);
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 + 16 * component;
            const adlite::Scalar force = value.contact_force * value.normal[component] +
                                         value.tributary_area * value.tangential_traction_vector[component];
            for (std::size_t node = 0; node < 8; ++node) {
                residual[offset + node] += secondary_shape[node] * force;
                residual[offset + 8 + node] -= value.primary_shape[node] * force;
            }
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianContactPointValue compute_quad8_to_quad8_contact_value(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history) {
    const CartesianContactAdValue8 value =
        evaluate_surface_mechanical(properties, geometry, make_ad_state(state, false), committed_state, history);
    return {value.projected, value.gap.value(), value.pressure.value(), value.tributary_area.value(),
        value.contact_force.value(), value.tangential_traction.value(), value.tangential_force.value(),
        value.friction_dissipation.value(), {value.normal[0].value(), value.normal[1].value(), value.normal[2].value()},
        {value.tangent_first[0].value(), value.tangent_first[1].value(), value.tangent_first[2].value()},
        {value.tangential_traction_vector[0].value(), value.tangential_traction_vector[1].value(),
            value.tangential_traction_vector[2].value()},
        {value.tangential_slip[0].value(), value.tangential_slip[1].value(), value.tangential_slip[2].value()},
        {value.elastic_tangential_slip[0].value(), value.elastic_tangential_slip[1].value(),
            value.elastic_tangential_slip[2].value()},
        value.sliding};
}

ContactProjectionValue compute_quad8_to_quad8_contact_projection(
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    if (!geometry.finite_sliding) return {true, small_sliding_projection(nodes, geometry).gap.value()};
    const ActivePoint3 secondary_point =
        interpolate_point(nodes, 0, active_values(geometry.secondary_displacement_shape));
    const SurfaceProjection8 projection = project_to_primary(secondary_point, nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

Quad8ReferenceProjectionValue compute_quad8_reference_projection(
    const std::array<CartesianPoint3, 8>& secondary_coordinates,
    const std::array<CartesianPoint3, 8>& primary_coordinates, const std::array<double, 8>& secondary_shape,
    double normal_orientation) {
    const Quad8SurfaceContactLocalValues state{};
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(secondary_coordinates, primary_coordinates, make_ad_state(state, false));
    const SurfaceProjection8 projection =
        project_to_primary(interpolate_point(nodes, 0, active_values(secondary_shape)), nodes, normal_orientation);
    Quad8ReferenceProjectionValue result{};
    result.projected = projection.projected;
    if (!projection.projected) return result;
    Quad8ShapeValues primary_values;
    quad8_shape(projection.xi, projection.eta, primary_values);
    for (std::size_t node = 0; node < result.primary_shape.size(); ++node)
        result.primary_shape[node] = projection.primary_shape[node].value();
    for (std::size_t node = 0; node < result.primary_derivative_xi.size(); ++node) {
        result.primary_derivative_xi[node] = primary_values.derivative_xi[node].value();
        result.primary_derivative_eta[node] = primary_values.derivative_eta[node].value();
    }
    result.normal = {projection.normal[0].value(), projection.normal[1].value(), projection.normal[2].value()};
    result.gap = projection.gap.value();
    return result;
}

Quad8SurfaceContactLocalResidual compute_node_to_quad8_contact(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, jacobian != nullptr);
    Quad8SurfaceContactLocalAdValues residual{};
    const CartesianContactAdValue8 value =
        evaluate_mechanical(properties, geometry, ad_state, committed_state, history);
    if (value.projected) {
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = 8 + 16 * component;
            const adlite::Scalar force = value.contact_force * value.normal[component] +
                                         value.tributary_area * value.tangential_traction_vector[component];
            residual[offset + geometry.secondary_local_node] += force;
            for (std::size_t node = 0; node < 8; ++node)
                residual[offset + 8 + node] -= value.primary_shape[node] * force;
        }
    }
    return extract(ad_state, residual, jacobian);
}

CartesianContactPointValue compute_node_to_quad8_contact_value(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history) {
    const CartesianContactAdValue8 value =
        evaluate_mechanical(properties, geometry, make_ad_state(state, false), committed_state, history);
    return {value.projected, value.gap.value(), value.pressure.value(), value.tributary_area.value(),
        value.contact_force.value(), value.tangential_traction.value(), value.tangential_force.value(),
        value.friction_dissipation.value(), {value.normal[0].value(), value.normal[1].value(), value.normal[2].value()},
        {},
        {value.tangential_traction_vector[0].value(), value.tangential_traction_vector[1].value(),
            value.tangential_traction_vector[2].value()},
        {value.tangential_slip[0].value(), value.tangential_slip[1].value(), value.tangential_slip[2].value()},
        {value.elastic_tangential_slip[0].value(), value.elastic_tangential_slip[1].value(),
            value.elastic_tangential_slip[2].value()},
        value.sliding};
}

ContactProjectionValue compute_node_to_quad8_contact_projection(
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    const Quad8SurfaceContactLocalAdValues ad_state = make_ad_state(state, false);
    const std::array<ActivePoint3, 16> nodes =
        current_nodes(geometry.secondary_coordinates, geometry.primary_coordinates, ad_state);
    const SurfaceProjection8 projection =
        project_to_primary(nodes[geometry.secondary_local_node], nodes, geometry.normal_orientation);
    return {projection.projected, projection.projected ? projection.gap.value() : 0.0};
}

} // namespace fuelsim
