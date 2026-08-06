#include "fuelsim/interface.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

struct HeatAdQuadratureValue final {
    adlite::Scalar gap;
    adlite::Scalar heat_flux;
    adlite::Scalar weighted_measure;
    adlite::Scalar primary_shape_0;
    adlite::Scalar primary_shape_1;
};

struct ContactAdValue final {
    bool projected = false;
    adlite::Scalar gap = 0.0;
    adlite::Scalar pressure = 0.0;
    adlite::Scalar tributary_area = 0.0;
    adlite::Scalar tributary_length = 0.0;
    adlite::Scalar contact_force = 0.0;
    adlite::Scalar primary_shape_0 = 0.0;
    adlite::Scalar primary_shape_1 = 0.0;
    adlite::Scalar normal_r = 0.0;
    adlite::Scalar normal_z = 0.0;
    adlite::Scalar tangent_r = 0.0;
    adlite::Scalar tangent_z = 0.0;
    adlite::Scalar tangential_traction = 0.0;
    adlite::Scalar tangential_force = 0.0;
    adlite::Scalar elastic_tangential_slip = 0.0;
    bool sliding = false;
};

bool finite_point(const RzPoint& point) {
    return std::isfinite(point.r) && std::isfinite(point.z);
}

void validate_line(const Line2InterfaceSideCoordinates& coordinates,
                   const char* name) {
    for (const RzPoint& point : coordinates) {
        if (!finite_point(point))
            throw std::invalid_argument(std::string(name) +
                                        " requires finite coordinates");
        if (!(point.r >= 0.0))
            throw std::invalid_argument(std::string(name) +
                                        " requires nonnegative radii");
    }
    const double dr = coordinates[1].r - coordinates[0].r;
    const double dz = coordinates[1].z - coordinates[0].z;
    if (!(std::hypot(dr, dz) > 0.0))
        throw std::invalid_argument(std::string(name) +
                                    " requires a nonzero line length");
}

double reference_projection_fraction(
    const RzPoint& secondary,
    const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double tangent_r =
        primary_coordinates[1].r - primary_coordinates[0].r;
    const double tangent_z =
        primary_coordinates[1].z - primary_coordinates[0].z;
    const double length_squared =
        tangent_r * tangent_r + tangent_z * tangent_z;
    return ((secondary.r - primary_coordinates[0].r) * tangent_r +
            (secondary.z - primary_coordinates[0].z) * tangent_z) /
           length_squared;
}

double reference_normal_orientation(
    const RzPoint& secondary,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    double primary_fraction) {
    const double tangent_r =
        primary_coordinates[1].r - primary_coordinates[0].r;
    const double tangent_z =
        primary_coordinates[1].z - primary_coordinates[0].z;
    const double length = std::hypot(tangent_r, tangent_z);
    const double delta_r =
        primary_coordinates[0].r + primary_fraction * tangent_r - secondary.r;
    const double delta_z =
        primary_coordinates[0].z + primary_fraction * tangent_z - secondary.z;
    const double raw_gap =
        delta_r * tangent_z / length - delta_z * tangent_r / length;
    if (raw_gap == 0.0)
        throw std::invalid_argument(
            "Contact surfaces require a positive reference normal gap");
    return raw_gap > 0.0 ? 1.0 : -1.0;
}

adlite::Scalar
interpolate(const std::array<double, line2_interface_side_node_count>& shape,
            const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar value = 0.0;
    for (std::size_t node = 0; node < line2_interface_side_node_count; ++node)
        value += shape[node] * state[offset + node];
    return value;
}

HeatAdQuadratureValue
evaluate_heat_quadrature(const Line2RzHeatGeometry& geometry,
                         const Line2RzHeatQuadraturePoint& point,
                         const LocalAdValues& state,
                         const GapHeatProperties& properties) {
    if (geometry.radial_reference_geometry) {
        const adlite::Scalar secondary_temperature =
            interpolate(point.secondary_shape, state, 0);
        const adlite::Scalar primary_temperature =
            interpolate(point.primary_shape, state, 2);
        const adlite::Scalar secondary_radius =
            point.secondary_shape[0] *
                (geometry.secondary_coordinates[0].r + state[4]) +
            point.secondary_shape[1] *
                (geometry.secondary_coordinates[1].r + state[5]);
        const adlite::Scalar primary_radius =
            point.primary_shape[0] *
                (geometry.primary_coordinates[0].r + state[6]) +
            point.primary_shape[1] *
                (geometry.primary_coordinates[1].r + state[7]);
        const adlite::Scalar gap = primary_radius - secondary_radius;
        const adlite::Scalar thermal_gap =
            adlite::max(gap, adlite::Scalar(properties.minimum_gap));
        const adlite::Scalar heat_flux =
            properties.gap_conductivity / thermal_gap *
            (secondary_temperature - primary_temperature);
        const adlite::Scalar secondary_radius_0 =
            geometry.secondary_coordinates[0].r + state[4];
        const adlite::Scalar secondary_radius_1 =
            geometry.secondary_coordinates[1].r + state[5];
        const adlite::Scalar secondary_axial_0 =
            geometry.secondary_coordinates[0].z + state[8];
        const adlite::Scalar secondary_axial_1 =
            geometry.secondary_coordinates[1].z + state[9];
        const adlite::Scalar surface_jacobian =
            0.5 * adlite::hypot(secondary_radius_1 - secondary_radius_0,
                                secondary_axial_1 - secondary_axial_0);
        return {gap,
                heat_flux,
                2.0 * pi * secondary_radius * surface_jacobian *
                    point.integration_weight,
                adlite::Scalar(point.primary_shape[0]),
                adlite::Scalar(point.primary_shape[1])};
    }

    const adlite::Scalar secondary_radius =
        point.secondary_shape[0] *
            (geometry.secondary_coordinates[0].r + state[4]) +
        point.secondary_shape[1] *
            (geometry.secondary_coordinates[1].r + state[5]);
    const adlite::Scalar secondary_axial =
        point.secondary_shape[0] *
            (geometry.secondary_coordinates[0].z + state[8]) +
        point.secondary_shape[1] *
            (geometry.secondary_coordinates[1].z + state[9]);
    const adlite::Scalar primary_radius_0 =
        geometry.primary_coordinates[0].r + state[6];
    const adlite::Scalar primary_radius_1 =
        geometry.primary_coordinates[1].r + state[7];
    const adlite::Scalar primary_axial_0 =
        geometry.primary_coordinates[0].z + state[10];
    const adlite::Scalar primary_axial_1 =
        geometry.primary_coordinates[1].z + state[11];
    const adlite::Scalar tangent_r = primary_radius_1 - primary_radius_0;
    const adlite::Scalar tangent_z = primary_axial_1 - primary_axial_0;
    const adlite::Scalar tangent_length = adlite::hypot(tangent_r, tangent_z);
    adlite::Scalar primary_fraction =
        ((secondary_radius - primary_radius_0) * tangent_r +
         (secondary_axial - primary_axial_0) * tangent_z) /
        (tangent_length * tangent_length);
    constexpr double projection_tolerance = 1.0e-12;
    if (primary_fraction.value() < -projection_tolerance ||
        primary_fraction.value() > 1.0 + projection_tolerance)
        throw std::domain_error(
            "Thermal contact quadrature point left its primary segment; the "
            "current small-sliding projection is no longer valid");
    if (primary_fraction.value() < 0.0)
        primary_fraction = 0.0;
    else if (primary_fraction.value() > 1.0)
        primary_fraction = 1.0;
    const adlite::Scalar primary_shape_0 = 1.0 - primary_fraction;
    const adlite::Scalar primary_shape_1 = primary_fraction;
    const adlite::Scalar primary_radius =
        primary_shape_0 * primary_radius_0 +
        primary_shape_1 * primary_radius_1;
    const adlite::Scalar primary_axial =
        primary_shape_0 * primary_axial_0 +
        primary_shape_1 * primary_axial_1;
    const adlite::Scalar normal_r =
        point.normal_orientation * tangent_z / tangent_length;
    const adlite::Scalar normal_z =
        -point.normal_orientation * tangent_r / tangent_length;
    const adlite::Scalar gap =
        (primary_radius - secondary_radius) * normal_r +
        (primary_axial - secondary_axial) * normal_z;

    const adlite::Scalar secondary_temperature =
        interpolate(point.secondary_shape, state, 0);
    const adlite::Scalar primary_temperature =
        primary_shape_0 * state[2] + primary_shape_1 * state[3];
    const adlite::Scalar thermal_gap =
        adlite::max(gap, adlite::Scalar(properties.minimum_gap));
    const adlite::Scalar conductance =
        properties.gap_conductivity / thermal_gap;
    const adlite::Scalar heat_flux =
        conductance * (secondary_temperature - primary_temperature);

    const adlite::Scalar secondary_radius_0 =
        geometry.secondary_coordinates[0].r + state[4];
    const adlite::Scalar secondary_radius_1 =
        geometry.secondary_coordinates[1].r + state[5];
    const adlite::Scalar secondary_axial_0 =
        geometry.secondary_coordinates[0].z + state[8];
    const adlite::Scalar secondary_axial_1 =
        geometry.secondary_coordinates[1].z + state[9];
    const adlite::Scalar dr_dxi =
        0.5 * (secondary_radius_1 - secondary_radius_0);
    const adlite::Scalar dz_dxi = 0.5 * (secondary_axial_1 - secondary_axial_0);
    const adlite::Scalar surface_jacobian =
        adlite::sqrt(dr_dxi * dr_dxi + dz_dxi * dz_dxi);
    const adlite::Scalar weighted_measure =
        2.0 * pi * secondary_radius * surface_jacobian *
        point.integration_weight;

    return {gap, heat_flux, weighted_measure, primary_shape_0,
            primary_shape_1};
}

bool projection_is_inside(double fraction, bool includes_second_endpoint) {
    if (fraction < 0.0)
        return false;
    if (includes_second_endpoint)
        return fraction <= 1.0;
    return fraction < 1.0;
}

bool clamp_owned_chain_endpoint(
    adlite::Scalar& fraction, double reference_fraction,
    bool primary_segment_is_first, bool includes_second_endpoint) {
    constexpr double endpoint_tolerance = 1.0e-12;
    const double value = fraction.value();
    if (primary_segment_is_first && value < 0.0 &&
        std::abs(reference_fraction) <= endpoint_tolerance) {
        fraction = 0.0;
        return true;
    }
    if (includes_second_endpoint && value > 1.0 &&
        std::abs(reference_fraction - 1.0) <= endpoint_tolerance) {
        fraction = 1.0;
        return true;
    }
    return false;
}

void evaluate_friction(ContactAdValue& value,
                       const NodeToLineRzContactGeometry& geometry,
                       const LocalAdValues& state,
                       const LocalValues& committed_state,
                       const ContactPointHistory& history,
                       const NormalContactProperties& properties) {
    if (properties.friction_coefficient == 0.0 ||
        !(value.pressure.value() > 0.0))
        return;

    const std::size_t secondary = geometry.secondary_local_node;
    const adlite::Scalar secondary_increment_r =
        state[4 + secondary] - committed_state[4 + secondary];
    const adlite::Scalar secondary_increment_z =
        state[8 + secondary] - committed_state[8 + secondary];
    const adlite::Scalar primary_increment_r =
        value.primary_shape_0 * (state[6] - committed_state[6]) +
        value.primary_shape_1 * (state[7] - committed_state[7]);
    const adlite::Scalar primary_increment_z =
        value.primary_shape_0 * (state[10] - committed_state[10]) +
        value.primary_shape_1 * (state[11] - committed_state[11]);
    const adlite::Scalar tangential_increment =
        (secondary_increment_r - primary_increment_r) * value.tangent_r +
        (secondary_increment_z - primary_increment_z) * value.tangent_z;
    const adlite::Scalar trial_slip =
        history.elastic_tangential_slip + tangential_increment;
    const adlite::Scalar trial_traction = properties.penalty * trial_slip;
    const adlite::Scalar sliding_limit =
        properties.friction_coefficient * value.pressure;

    const double trial_magnitude = std::abs(trial_traction.value());
    if (trial_magnitude < sliding_limit.value() ||
        (trial_magnitude == sliding_limit.value() && !history.sliding)) {
        value.tangential_traction = trial_traction;
        value.elastic_tangential_slip = trial_slip;
    } else {
        const double direction = trial_traction.value() < 0.0 ? -1.0 : 1.0;
        value.tangential_traction = direction * sliding_limit;
        value.elastic_tangential_slip =
            value.tangential_traction / properties.penalty;
        value.sliding = true;
    }
    value.tangential_force =
        value.tangential_traction * value.tributary_area;
}

ContactAdValue evaluate_contact(const NodeToLineRzContactGeometry& geometry,
                                const LocalAdValues& state,
                                const LocalValues& committed_state,
                                const ContactPointHistory& history,
                                const NormalContactProperties& properties) {
    const std::size_t secondary = geometry.secondary_local_node;
    const std::size_t other = secondary == 0 ? 1 : 0;

    if (geometry.radial_reference_geometry) {
        const adlite::Scalar secondary_z =
            geometry.secondary_edge_coordinates[secondary].z +
            state[8 + secondary];
        const adlite::Scalar primary_z_0 =
            geometry.primary_segment_coordinates[0].z + state[10];
        const adlite::Scalar primary_z_1 =
            geometry.primary_segment_coordinates[1].z + state[11];
        adlite::Scalar fraction =
            (secondary_z - primary_z_0) / (primary_z_1 - primary_z_0);
        bool projected = projection_is_inside(
            fraction.value(),
            geometry.primary_segment_includes_second_endpoint);
        if (!projected)
            projected = clamp_owned_chain_endpoint(
                fraction, geometry.reference_primary_fraction,
                geometry.primary_segment_is_first,
                geometry.primary_segment_includes_second_endpoint);
        if (!projected)
            return {};
        const adlite::Scalar shape_0 = 1.0 - fraction;
        const adlite::Scalar shape_1 = fraction;
        const adlite::Scalar secondary_radius =
            geometry.secondary_edge_coordinates[secondary].r +
            state[4 + secondary];
        const adlite::Scalar primary_radius =
            shape_0 *
                (geometry.primary_segment_coordinates[0].r + state[6]) +
            shape_1 *
                (geometry.primary_segment_coordinates[1].r + state[7]);
        const adlite::Scalar gap = primary_radius - secondary_radius;
        const adlite::Scalar pressure =
            properties.penalty * adlite::max(-gap, adlite::Scalar(0.0));
        const adlite::Scalar other_radius =
            geometry.secondary_edge_coordinates[other].r + state[4 + other];
        const adlite::Scalar other_z =
            geometry.secondary_edge_coordinates[other].z + state[8 + other];
        const adlite::Scalar edge_length =
            adlite::hypot(other_radius - secondary_radius,
                          other_z - secondary_z);
        const adlite::Scalar tributary_length = 0.5 * edge_length;
        const adlite::Scalar tributary_area =
            2.0 * pi * secondary_radius * tributary_length;
        const adlite::Scalar force = pressure * tributary_area;
        ContactAdValue result;
        result.projected = true;
        result.gap = gap;
        result.pressure = pressure;
        result.tributary_area = tributary_area;
        result.tributary_length = tributary_length;
        result.contact_force = force;
        result.primary_shape_0 = shape_0;
        result.primary_shape_1 = shape_1;
        result.normal_r = 1.0;
        result.tangent_z = 1.0;
        evaluate_friction(result, geometry, state, committed_state, history,
                          properties);
        return result;
    }

    const adlite::Scalar secondary_radius =
        geometry.secondary_edge_coordinates[secondary].r + state[4 + secondary];
    const adlite::Scalar secondary_z =
        geometry.secondary_edge_coordinates[secondary].z + state[8 + secondary];
    const adlite::Scalar primary_radius_0 =
        geometry.primary_segment_coordinates[0].r + state[6];
    const adlite::Scalar primary_radius_1 =
        geometry.primary_segment_coordinates[1].r + state[7];
    const adlite::Scalar primary_z_0 =
        geometry.primary_segment_coordinates[0].z + state[10];
    const adlite::Scalar primary_z_1 =
        geometry.primary_segment_coordinates[1].z + state[11];
    const adlite::Scalar tangent_r = primary_radius_1 - primary_radius_0;
    const adlite::Scalar tangent_z = primary_z_1 - primary_z_0;
    const adlite::Scalar tangent_length = adlite::hypot(tangent_r, tangent_z);
    adlite::Scalar primary_fraction =
        ((secondary_radius - primary_radius_0) * tangent_r +
         (secondary_z - primary_z_0) * tangent_z) /
        (tangent_length * tangent_length);
    bool projected = projection_is_inside(
        primary_fraction.value(),
        geometry.primary_segment_includes_second_endpoint);
    if (!projected)
        projected = clamp_owned_chain_endpoint(
            primary_fraction, geometry.reference_primary_fraction,
            geometry.primary_segment_is_first,
            geometry.primary_segment_includes_second_endpoint);

    if (!projected) {
        return {};
    }

    const adlite::Scalar primary_shape_0 = 1.0 - primary_fraction;
    const adlite::Scalar primary_shape_1 = primary_fraction;
    const adlite::Scalar primary_radius =
        primary_shape_0 * primary_radius_0 + primary_shape_1 * primary_radius_1;
    const adlite::Scalar primary_z =
        primary_shape_0 * primary_z_0 + primary_shape_1 * primary_z_1;
    const adlite::Scalar normal_r =
        geometry.normal_orientation * tangent_z / tangent_length;
    const adlite::Scalar normal_z =
        -geometry.normal_orientation * tangent_r / tangent_length;
    const adlite::Scalar gap =
        (primary_radius - secondary_radius) * normal_r +
        (primary_z - secondary_z) * normal_z;
    const adlite::Scalar penetration = adlite::max(-gap, adlite::Scalar(0.0));
    const adlite::Scalar pressure = properties.penalty * penetration;

    const adlite::Scalar other_radius =
        geometry.secondary_edge_coordinates[other].r + state[4 + other];
    const adlite::Scalar other_z =
        geometry.secondary_edge_coordinates[other].z + state[8 + other];
    const adlite::Scalar dr = other_radius - secondary_radius;
    const adlite::Scalar dz = other_z - secondary_z;
    const adlite::Scalar edge_length = adlite::sqrt(dr * dr + dz * dz);
    const adlite::Scalar tributary_length = 0.5 * edge_length;
    const adlite::Scalar tributary_area =
        2.0 * pi * 0.5 * edge_length *
        (2.0 * secondary_radius + other_radius) / 3.0;
    const adlite::Scalar contact_force = pressure * tributary_area;

    ContactAdValue result;
    result.projected = true;
    result.gap = gap;
    result.pressure = pressure;
    result.tributary_area = tributary_area;
    result.tributary_length = tributary_length;
    result.contact_force = contact_force;
    result.primary_shape_0 = primary_shape_0;
    result.primary_shape_1 = primary_shape_1;
    result.normal_r = normal_r;
    result.normal_z = normal_z;
    result.tangent_r = tangent_r / tangent_length;
    result.tangent_z = tangent_z / tangent_length;
    evaluate_friction(result, geometry, state, committed_state, history,
                      properties);
    return result;
}

} // namespace

Line2RzHeatGeometry make_line2_rz_heat_geometry(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates) {
    return make_line2_rz_heat_geometry(secondary_coordinates,
                                       primary_coordinates, -1.0, 1.0);
}

Line2RzHeatGeometry make_line2_rz_heat_geometry(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    double secondary_coordinate_lower, double secondary_coordinate_upper) {
    validate_line(secondary_coordinates, "Line2RzHeatGeometry secondary");
    validate_line(primary_coordinates, "Line2RzHeatGeometry primary");
    if (!std::isfinite(secondary_coordinate_lower) ||
        !std::isfinite(secondary_coordinate_upper) ||
        secondary_coordinate_lower < -1.0 ||
        secondary_coordinate_upper > 1.0 ||
        !(secondary_coordinate_upper > secondary_coordinate_lower))
        throw std::invalid_argument(
            "Line2RzHeatGeometry requires a nonempty secondary interval in "
            "[-1,1]");

    constexpr double gauss = 0.577350269189625764509148780501957456;
    const std::array<double, line2_interface_quadrature_point_count> locations =
        {-gauss, gauss};

    Line2RzHeatGeometry geometry{};
    geometry.secondary_coordinates = secondary_coordinates;
    geometry.primary_coordinates = primary_coordinates;
    geometry.radial_reference_geometry =
        secondary_coordinates[0].r == secondary_coordinates[1].r &&
        primary_coordinates[0].r == primary_coordinates[1].r &&
        primary_coordinates[0].r > secondary_coordinates[0].r;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const double secondary_xi =
            0.5 * ((1.0 - locations[q]) * secondary_coordinate_lower +
                   (1.0 + locations[q]) * secondary_coordinate_upper);
        const std::array<double, 2> secondary_shape = {
            0.5 * (1.0 - secondary_xi),
            0.5 * (1.0 + secondary_xi),
        };
        const RzPoint secondary_point = {
            secondary_shape[0] * secondary_coordinates[0].r +
                secondary_shape[1] * secondary_coordinates[1].r,
            secondary_shape[0] * secondary_coordinates[0].z +
                secondary_shape[1] * secondary_coordinates[1].z,
        };
        const double primary_fraction =
            reference_projection_fraction(secondary_point,
                                          primary_coordinates);
        if (primary_fraction < 0.0 || primary_fraction > 1.0)
            throw std::invalid_argument(
                "Line2RzHeatGeometry secondary Gauss point does not project "
                "inside the primary segment");

        Line2RzHeatQuadraturePoint& point = geometry.points[q];
        point.secondary_shape = secondary_shape;
        point.primary_shape = {
            1.0 - primary_fraction,
            primary_fraction,
        };
        point.integration_weight =
            0.5 * (secondary_coordinate_upper - secondary_coordinate_lower);
        point.normal_orientation = reference_normal_orientation(
            secondary_point, primary_coordinates, primary_fraction);
    }
    return geometry;
}

Line2RzGapHeatKernel::Line2RzGapHeatKernel(GapHeatProperties properties)
    : _properties(properties) {
    if (!std::isfinite(_properties.gap_conductivity) ||
        !(_properties.gap_conductivity >= 0.0))
        throw std::invalid_argument(
            "GapHeatProperties gap_conductivity must be finite and "
            "nonnegative");
    if (!std::isfinite(_properties.minimum_gap) ||
        !(_properties.minimum_gap > 0.0))
        throw std::invalid_argument(
            "GapHeatProperties minimum_gap must be finite and positive");
}

const GapHeatProperties& Line2RzGapHeatKernel::properties() const noexcept {
    return _properties;
}

LocalResidual
Line2RzGapHeatKernel::residual(const Line2RzHeatGeometry& geometry,
                               const LocalValues& state) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];
    LocalAdValues passive_residual{};
    residual_ad(geometry, passive_state, passive_residual);

    LocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = passive_residual[row].value();
    return result;
}

LocalSystem Line2RzGapHeatKernel::linearize(const Line2RzHeatGeometry& geometry,
                                            const LocalValues& state) const {
    LocalAdValues active_state{};
    adlite::seed_identity(state.data(), state.size(), active_state.data());
    LocalAdValues active_residual{};
    residual_ad(geometry, active_state, active_residual);

    LocalSystem result{};
    adlite::extract_jacobian(active_residual.data(), active_residual.size(),
                             active_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

HeatQuadratureValues
Line2RzGapHeatKernel::quadrature_values(const Line2RzHeatGeometry& geometry,
                                        const LocalValues& state) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];

    HeatQuadratureValues result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const HeatAdQuadratureValue value = evaluate_heat_quadrature(
            geometry, geometry.points[q], passive_state, _properties);
        result[q] = {
            value.gap.value(),
            value.heat_flux.value(),
            value.weighted_measure.value(),
        };
    }
    return result;
}

void Line2RzGapHeatKernel::residual_ad(const Line2RzHeatGeometry& geometry,
                                       const LocalAdValues& state,
                                       LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));

    for (const Line2RzHeatQuadraturePoint& point : geometry.points) {
        const HeatAdQuadratureValue value =
            evaluate_heat_quadrature(geometry, point, state, _properties);
        for (std::size_t node = 0; node < line2_interface_side_node_count;
             ++node) {
            residual[node] += value.weighted_measure *
                              point.secondary_shape[node] * value.heat_flux;
            const adlite::Scalar primary_shape =
                node == 0 ? value.primary_shape_0 : value.primary_shape_1;
            residual[2 + node] -=
                value.weighted_measure * primary_shape * value.heat_flux;
        }
    }
}

NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates,
    std::size_t secondary_local_node,
    bool primary_segment_is_first,
    bool primary_segment_includes_upper_endpoint) {
    validate_line(secondary_edge_coordinates,
                  "NodeToLineRzContactGeometry secondary");
    validate_line(primary_segment_coordinates,
                  "NodeToLineRzContactGeometry primary");
    if (secondary_local_node >= line2_interface_side_node_count)
        throw std::invalid_argument(
            "NodeToLineRzContactGeometry secondary node is out of range");
    const RzPoint secondary_point =
        secondary_edge_coordinates[secondary_local_node];
    const double primary_fraction = reference_projection_fraction(
        secondary_point, primary_segment_coordinates);
    const double closest_fraction =
        std::max(0.0, std::min(1.0, primary_fraction));
    const double orientation = reference_normal_orientation(
        secondary_point, primary_segment_coordinates, closest_fraction);

    return {
        secondary_edge_coordinates,
        primary_segment_coordinates,
        secondary_local_node,
        primary_segment_is_first,
        primary_segment_includes_upper_endpoint,
        orientation,
        primary_fraction,
        secondary_edge_coordinates[0].r ==
                secondary_edge_coordinates[1].r &&
            primary_segment_coordinates[0].r ==
                primary_segment_coordinates[1].r &&
            primary_segment_coordinates[0].r >
                secondary_edge_coordinates[secondary_local_node].r,
    };
}

NodeToLineRzContactKernel::NodeToLineRzContactKernel(
    NormalContactProperties properties)
    : _properties(properties) {
    if (!std::isfinite(_properties.penalty) || !(_properties.penalty >= 0.0))
        throw std::invalid_argument(
            "NormalContactProperties penalty must be finite and nonnegative");
    if (!std::isfinite(_properties.friction_coefficient) ||
        !(_properties.friction_coefficient >= 0.0))
        throw std::invalid_argument(
            "NormalContactProperties friction_coefficient must be finite and "
            "nonnegative");
    if (_properties.friction_coefficient > 0.0 &&
        !(_properties.penalty > 0.0))
        throw std::invalid_argument(
            "Frictional contact requires a positive penalty");
}

const NormalContactProperties&
NodeToLineRzContactKernel::properties() const noexcept {
    return _properties;
}

LocalResidual
NodeToLineRzContactKernel::residual(const NodeToLineRzContactGeometry& geometry,
                                    const LocalValues& state,
                                    const LocalValues& committed_state,
                                    const ContactPointHistory& history) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];
    LocalAdValues passive_residual{};
    residual_ad(geometry, passive_state, committed_state, history,
                passive_residual);

    LocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = passive_residual[row].value();
    return result;
}

LocalSystem NodeToLineRzContactKernel::linearize(
    const NodeToLineRzContactGeometry& geometry,
    const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history) const {
    LocalAdValues active_state{};
    adlite::seed_identity(state.data(), state.size(), active_state.data());
    LocalAdValues active_residual{};
    residual_ad(geometry, active_state, committed_state, history,
                active_residual);

    LocalSystem result{};
    adlite::extract_jacobian(active_residual.data(), active_residual.size(),
                             active_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

ContactPointValue
NodeToLineRzContactKernel::value(const NodeToLineRzContactGeometry& geometry,
                                 const LocalValues& state,
                                 const LocalValues& committed_state,
                                 const ContactPointHistory& history) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];
    const ContactAdValue result =
        evaluate_contact(geometry, passive_state, committed_state, history,
                         _properties);
    return {
        result.projected,
        result.gap.value(),
        result.pressure.value(),
        result.tributary_area.value(),
        result.tributary_length.value(),
        result.contact_force.value(),
        result.tangential_traction.value(),
        result.tangential_force.value(),
        result.elastic_tangential_slip.value(),
        result.sliding,
    };
}

ContactPointHistory NodeToLineRzContactKernel::trial_history(
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state,
    const LocalValues& committed_state,
    const ContactPointHistory& history) const {
    const ContactPointValue trial =
        value(geometry, state, committed_state, history);
    if (!trial.projected)
        throw std::domain_error(
            "Cannot update friction history for an unprojected contact node");
    return {trial.elastic_tangential_slip, trial.sliding};
}

void NodeToLineRzContactKernel::residual_ad(
    const NodeToLineRzContactGeometry& geometry, const LocalAdValues& state,
    const LocalValues& committed_state, const ContactPointHistory& history,
    LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));
    const ContactAdValue value = evaluate_contact(
        geometry, state, committed_state, history, _properties);
    if (!value.projected)
        return;

    const std::size_t secondary = geometry.secondary_local_node;
    residual[4 + secondary] += value.contact_force * value.normal_r;
    residual[6] -= value.primary_shape_0 * value.contact_force * value.normal_r;
    residual[7] -= value.primary_shape_1 * value.contact_force * value.normal_r;
    residual[8 + secondary] += value.contact_force * value.normal_z;
    residual[10] -= value.primary_shape_0 * value.contact_force * value.normal_z;
    residual[11] -= value.primary_shape_1 * value.contact_force * value.normal_z;
    if (_properties.friction_coefficient == 0.0)
        return;
    residual[4 + secondary] += value.tangential_force * value.tangent_r;
    residual[6] -=
        value.primary_shape_0 * value.tangential_force * value.tangent_r;
    residual[7] -=
        value.primary_shape_1 * value.tangential_force * value.tangent_r;
    residual[8 + secondary] += value.tangential_force * value.tangent_z;
    residual[10] -=
        value.primary_shape_0 * value.tangential_force * value.tangent_z;
    residual[11] -=
        value.primary_shape_1 * value.tangential_force * value.tangent_z;
}

} // namespace fuelsim
