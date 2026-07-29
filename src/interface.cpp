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
};

struct ContactAdValue final {
    bool projected;
    adlite::Scalar gap;
    adlite::Scalar pressure;
    adlite::Scalar tributary_area;
    adlite::Scalar tributary_length;
    adlite::Scalar contact_force;
    adlite::Scalar primary_shape_0;
    adlite::Scalar primary_shape_1;
};

bool finite_point(const RzPoint& point) {
    return std::isfinite(point.r) && std::isfinite(point.z);
}

void validate_cylindrical_line(const Line2InterfaceSideCoordinates& coordinates,
                               const char* name) {
    for (const RzPoint& point : coordinates) {
        if (!finite_point(point))
            throw std::invalid_argument(std::string(name) +
                                        " requires finite coordinates");
        if (!(point.r > 0.0))
            throw std::invalid_argument(std::string(name) +
                                        " requires positive radii");
    }
    if (!(coordinates[1].z > coordinates[0].z))
        throw std::invalid_argument(std::string(name) +
                                    " requires bottom-to-top nodes");
    if (coordinates[0].r != coordinates[1].r)
        throw std::invalid_argument(std::string(name) +
                                    " requires a cylindrical axial line");
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
    const adlite::Scalar fuel_temperature =
        interpolate(point.fuel_shape, state, 0);
    const adlite::Scalar cladding_temperature =
        interpolate(point.cladding_shape, state, 2);
    const adlite::Scalar fuel_radial_displacement =
        interpolate(point.fuel_shape, state, 4);
    const adlite::Scalar cladding_radial_displacement =
        interpolate(point.cladding_shape, state, 6);

    const adlite::Scalar fuel_radius =
        point.fuel_reference_radius + fuel_radial_displacement;
    const adlite::Scalar cladding_radius =
        point.cladding_reference_radius + cladding_radial_displacement;
    const adlite::Scalar gap = cladding_radius - fuel_radius;
    const adlite::Scalar thermal_gap =
        adlite::max(gap, adlite::Scalar(properties.minimum_gap));
    const adlite::Scalar conductance =
        properties.gap_conductivity / thermal_gap;
    const adlite::Scalar heat_flux =
        conductance * (fuel_temperature - cladding_temperature);

    const adlite::Scalar fuel_radius_0 =
        geometry.fuel_coordinates[0].r + state[4];
    const adlite::Scalar fuel_radius_1 =
        geometry.fuel_coordinates[1].r + state[5];
    const adlite::Scalar fuel_axial_0 =
        geometry.fuel_coordinates[0].z + state[8];
    const adlite::Scalar fuel_axial_1 =
        geometry.fuel_coordinates[1].z + state[9];
    const adlite::Scalar dr_dxi = 0.5 * (fuel_radius_1 - fuel_radius_0);
    const adlite::Scalar dz_dxi = 0.5 * (fuel_axial_1 - fuel_axial_0);
    const adlite::Scalar surface_jacobian =
        adlite::sqrt(dr_dxi * dr_dxi + dz_dxi * dz_dxi);
    const adlite::Scalar weighted_measure =
        2.0 * pi * fuel_radius * surface_jacobian;

    return {gap, heat_flux, weighted_measure};
}

bool projection_is_inside(double secondary_z, double primary_z_0,
                          double primary_z_1, bool includes_upper_endpoint) {
    if (secondary_z < primary_z_0)
        return false;
    if (includes_upper_endpoint)
        return secondary_z <= primary_z_1;
    return secondary_z < primary_z_1;
}

ContactAdValue evaluate_contact(const NodeToLineRzContactGeometry& geometry,
                                const LocalAdValues& state,
                                const NormalContactProperties& properties) {
    const std::size_t secondary = geometry.secondary_local_node;
    const std::size_t other = secondary == 0 ? 1 : 0;

    const adlite::Scalar fuel_z =
        geometry.fuel_edge_coordinates[secondary].z + state[8 + secondary];
    const adlite::Scalar primary_z_0 =
        geometry.cladding_segment_coordinates[0].z + state[10];
    const adlite::Scalar primary_z_1 =
        geometry.cladding_segment_coordinates[1].z + state[11];
    const bool projected = projection_is_inside(
        fuel_z.value(), primary_z_0.value(), primary_z_1.value(),
        geometry.cladding_segment_includes_upper_endpoint);

    if (!projected) {
        return {
            false,
            adlite::Scalar(0.0),
            adlite::Scalar(0.0),
            adlite::Scalar(0.0),
            adlite::Scalar(0.0),
            adlite::Scalar(0.0),
            adlite::Scalar(0.0),
            adlite::Scalar(0.0),
        };
    }

    const adlite::Scalar axial_fraction =
        (fuel_z - primary_z_0) / (primary_z_1 - primary_z_0);
    const adlite::Scalar primary_shape_0 = 1.0 - axial_fraction;
    const adlite::Scalar primary_shape_1 = axial_fraction;

    const adlite::Scalar fuel_radius =
        geometry.fuel_edge_coordinates[secondary].r + state[4 + secondary];
    const adlite::Scalar primary_radius_0 =
        geometry.cladding_segment_coordinates[0].r + state[6];
    const adlite::Scalar primary_radius_1 =
        geometry.cladding_segment_coordinates[1].r + state[7];
    const adlite::Scalar primary_radius =
        primary_shape_0 * primary_radius_0 + primary_shape_1 * primary_radius_1;
    const adlite::Scalar gap = primary_radius - fuel_radius;
    const adlite::Scalar penetration = adlite::max(-gap, adlite::Scalar(0.0));
    const adlite::Scalar pressure = properties.penalty * penetration;

    const adlite::Scalar other_radius =
        geometry.fuel_edge_coordinates[other].r + state[4 + other];
    const adlite::Scalar other_z =
        geometry.fuel_edge_coordinates[other].z + state[8 + other];
    const adlite::Scalar dr = other_radius - fuel_radius;
    const adlite::Scalar dz = other_z - fuel_z;
    const adlite::Scalar edge_length = adlite::sqrt(dr * dr + dz * dz);
    const adlite::Scalar tributary_length = 0.5 * edge_length;
    const adlite::Scalar tributary_area =
        2.0 * pi * fuel_radius * tributary_length;
    const adlite::Scalar contact_force = pressure * tributary_area;

    return {
        true,
        gap,
        pressure,
        tributary_area,
        tributary_length,
        contact_force,
        primary_shape_0,
        primary_shape_1,
    };
}

} // namespace

Line2RzHeatGeometry make_line2_rz_heat_geometry(
    const Line2InterfaceSideCoordinates& fuel_coordinates,
    const Line2InterfaceSideCoordinates& cladding_coordinates) {
    validate_cylindrical_line(fuel_coordinates, "Line2RzHeatGeometry fuel");
    validate_cylindrical_line(cladding_coordinates,
                              "Line2RzHeatGeometry cladding");
    if (!(cladding_coordinates[0].r > fuel_coordinates[0].r))
        throw std::invalid_argument(
            "Line2RzHeatGeometry requires a positive reference gap");

    constexpr double gauss = 0.577350269189625764509148780501957456;
    const std::array<double, line2_interface_quadrature_point_count> locations =
        {-gauss, gauss};

    Line2RzHeatGeometry geometry{};
    geometry.fuel_coordinates = fuel_coordinates;
    geometry.cladding_coordinates = cladding_coordinates;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const double fuel_xi = locations[q];
        const std::array<double, 2> fuel_shape = {
            0.5 * (1.0 - fuel_xi),
            0.5 * (1.0 + fuel_xi),
        };
        const double fuel_z = fuel_shape[0] * fuel_coordinates[0].z +
                              fuel_shape[1] * fuel_coordinates[1].z;
        const double cladding_fraction =
            (fuel_z - cladding_coordinates[0].z) /
            (cladding_coordinates[1].z - cladding_coordinates[0].z);
        if (cladding_fraction < 0.0 || cladding_fraction > 1.0)
            throw std::invalid_argument(
                "Line2RzHeatGeometry fuel Gauss point does not project "
                "inside the cladding segment");

        Line2RzHeatQuadraturePoint& point = geometry.points[q];
        point.fuel_shape = fuel_shape;
        point.cladding_shape = {
            1.0 - cladding_fraction,
            cladding_fraction,
        };
        point.fuel_reference_radius = fuel_shape[0] * fuel_coordinates[0].r +
                                      fuel_shape[1] * fuel_coordinates[1].r;
        point.cladding_reference_radius =
            point.cladding_shape[0] * cladding_coordinates[0].r +
            point.cladding_shape[1] * cladding_coordinates[1].r;
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
            residual[node] += value.weighted_measure * point.fuel_shape[node] *
                              value.heat_flux;
            residual[2 + node] -= value.weighted_measure *
                                  point.cladding_shape[node] * value.heat_flux;
        }
    }
}

NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& fuel_edge_coordinates,
    const Line2InterfaceSideCoordinates& cladding_segment_coordinates,
    std::size_t secondary_local_node,
    bool cladding_segment_includes_upper_endpoint) {
    validate_cylindrical_line(fuel_edge_coordinates,
                              "NodeToLineRzContactGeometry fuel");
    validate_cylindrical_line(cladding_segment_coordinates,
                              "NodeToLineRzContactGeometry cladding");
    if (secondary_local_node >= line2_interface_side_node_count)
        throw std::invalid_argument(
            "NodeToLineRzContactGeometry secondary node is out of range");
    if (!(cladding_segment_coordinates[0].r >
          fuel_edge_coordinates[secondary_local_node].r))
        throw std::invalid_argument(
            "NodeToLineRzContactGeometry requires a positive reference gap");

    return {
        fuel_edge_coordinates,
        cladding_segment_coordinates,
        secondary_local_node,
        cladding_segment_includes_upper_endpoint,
    };
}

NodeToLineRzContactKernel::NodeToLineRzContactKernel(
    NormalContactProperties properties)
    : _properties(properties) {
    if (!std::isfinite(_properties.penalty) || !(_properties.penalty >= 0.0))
        throw std::invalid_argument(
            "NormalContactProperties penalty must be finite and nonnegative");
}

const NormalContactProperties&
NodeToLineRzContactKernel::properties() const noexcept {
    return _properties;
}

LocalResidual
NodeToLineRzContactKernel::residual(const NodeToLineRzContactGeometry& geometry,
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

LocalSystem NodeToLineRzContactKernel::linearize(
    const NodeToLineRzContactGeometry& geometry,
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

ContactPointValue
NodeToLineRzContactKernel::value(const NodeToLineRzContactGeometry& geometry,
                                 const LocalValues& state) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];
    const ContactAdValue result =
        evaluate_contact(geometry, passive_state, _properties);
    return {
        result.projected,
        result.gap.value(),
        result.pressure.value(),
        result.tributary_area.value(),
        result.tributary_length.value(),
        result.contact_force.value(),
    };
}

void NodeToLineRzContactKernel::residual_ad(
    const NodeToLineRzContactGeometry& geometry, const LocalAdValues& state,
    LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));
    const ContactAdValue value = evaluate_contact(geometry, state, _properties);
    if (!value.projected)
        return;

    const std::size_t secondary = geometry.secondary_local_node;
    residual[4 + secondary] += value.contact_force;
    residual[6] -= value.primary_shape_0 * value.contact_force;
    residual[7] -= value.primary_shape_1 * value.contact_force;
}

} // namespace fuelsim
