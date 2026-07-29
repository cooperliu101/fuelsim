#include "fuelsim/interface.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

struct InterfaceAdQuadratureValue final {
    adlite::Scalar gap;
    adlite::Scalar heat_flux;
    adlite::Scalar pressure;
};

bool finite_point(const RzPoint& point) {
    return std::isfinite(point.r) && std::isfinite(point.z);
}

adlite::Scalar
interpolate(const std::array<double, line2_interface_side_node_count>& shape,
            const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar value = 0.0;
    for (std::size_t node = 0; node < line2_interface_side_node_count; ++node)
        value += shape[node] * state[offset + node];
    return value;
}

InterfaceAdQuadratureValue
evaluate_quadrature_value(const Line2RzInterfaceQuadraturePoint& point,
                          const LocalAdValues& state,
                          const GapContactProperties& properties) {
    const adlite::Scalar fuel_temperature = interpolate(point.shape, state, 0);
    const adlite::Scalar clad_temperature = interpolate(point.shape, state, 2);
    const adlite::Scalar fuel_radial_displacement =
        interpolate(point.shape, state, 4);
    const adlite::Scalar clad_radial_displacement =
        interpolate(point.shape, state, 6);

    const adlite::Scalar fuel_radius =
        point.fuel_radius + fuel_radial_displacement;
    const adlite::Scalar clad_radius =
        point.clad_radius + clad_radial_displacement;
    const adlite::Scalar gap = clad_radius - fuel_radius;
    const adlite::Scalar thermal_gap =
        adlite::max(gap, adlite::Scalar(properties.minimum_gap));
    const adlite::Scalar conductance =
        properties.gap_conductivity / thermal_gap;
    const adlite::Scalar heat_flux =
        conductance * (fuel_temperature - clad_temperature);
    const adlite::Scalar penetration = adlite::max(-gap, adlite::Scalar(0.0));
    const adlite::Scalar pressure = properties.penalty * penetration;

    return {gap, heat_flux, pressure};
}

} // namespace

Line2RzInterfaceGeometry make_line2_rz_interface_geometry(
    const Line2InterfaceSideCoordinates& fuel_coordinates,
    const Line2InterfaceSideCoordinates& clad_coordinates) {
    for (const RzPoint& point : fuel_coordinates) {
        if (!finite_point(point))
            throw std::invalid_argument(
                "Line2RzInterfaceGeometry requires finite fuel coordinates");
        if (!(point.r > 0.0))
            throw std::invalid_argument(
                "Line2RzInterfaceGeometry requires positive fuel radii");
    }
    for (const RzPoint& point : clad_coordinates) {
        if (!finite_point(point))
            throw std::invalid_argument(
                "Line2RzInterfaceGeometry requires finite clad coordinates");
        if (!(point.r > 0.0))
            throw std::invalid_argument(
                "Line2RzInterfaceGeometry requires positive clad radii");
    }

    if (!(fuel_coordinates[1].z > fuel_coordinates[0].z))
        throw std::invalid_argument(
            "Line2RzInterfaceGeometry requires bottom-to-top fuel nodes");
    if (clad_coordinates[0].z != fuel_coordinates[0].z ||
        clad_coordinates[1].z != fuel_coordinates[1].z)
        throw std::invalid_argument(
            "Line2RzInterfaceGeometry requires matching axial coordinates");
    if (fuel_coordinates[0].r != fuel_coordinates[1].r ||
        clad_coordinates[0].r != clad_coordinates[1].r)
        throw std::invalid_argument(
            "Line2RzInterfaceGeometry requires axial cylindrical sides");
    if (!(clad_coordinates[0].r > fuel_coordinates[0].r) ||
        !(clad_coordinates[1].r > fuel_coordinates[1].r))
        throw std::invalid_argument(
            "Line2RzInterfaceGeometry requires a positive reference gap");

    constexpr double gauss = 0.577350269189625764509148780501957456;
    const std::array<double, line2_interface_quadrature_point_count> locations =
        {-gauss, gauss};
    const double line_jacobian =
        0.5 * (fuel_coordinates[1].z - fuel_coordinates[0].z);

    Line2RzInterfaceGeometry geometry{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const double xi = locations[q];
        Line2RzInterfaceQuadraturePoint& point = geometry.points[q];
        point.shape = {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
        for (std::size_t node = 0; node < line2_interface_side_node_count;
             ++node) {
            point.fuel_radius += point.shape[node] * fuel_coordinates[node].r;
            point.clad_radius += point.shape[node] * clad_coordinates[node].r;
        }
        point.weighted_measure = 2.0 * pi * point.fuel_radius * line_jacobian;
    }

    return geometry;
}

Line2RzGapContactKernel::Line2RzGapContactKernel(
    GapContactProperties properties)
    : properties_(properties) {
    if (!std::isfinite(properties_.gap_conductivity) ||
        !(properties_.gap_conductivity >= 0.0))
        throw std::invalid_argument(
            "GapContactProperties gap_conductivity must be finite and "
            "nonnegative");
    if (!std::isfinite(properties_.minimum_gap) ||
        !(properties_.minimum_gap > 0.0))
        throw std::invalid_argument(
            "GapContactProperties minimum_gap must be finite and positive");
    if (!std::isfinite(properties_.penalty) || !(properties_.penalty >= 0.0))
        throw std::invalid_argument(
            "GapContactProperties penalty must be finite and nonnegative");
}

const GapContactProperties&
Line2RzGapContactKernel::properties() const noexcept {
    return properties_;
}

LocalResidual
Line2RzGapContactKernel::residual(const Line2RzInterfaceGeometry& geometry,
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

LocalSystem
Line2RzGapContactKernel::linearize(const Line2RzInterfaceGeometry& geometry,
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

InterfaceQuadratureValues Line2RzGapContactKernel::quadrature_values(
    const Line2RzInterfaceGeometry& geometry, const LocalValues& state) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];

    InterfaceQuadratureValues result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const InterfaceAdQuadratureValue value = evaluate_quadrature_value(
            geometry.points[q], passive_state, properties_);
        result[q] = {
            value.gap.value(),
            value.heat_flux.value(),
            value.pressure.value(),
        };
    }
    return result;
}

void Line2RzGapContactKernel::residual_ad(
    const Line2RzInterfaceGeometry& geometry, const LocalAdValues& state,
    LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));

    for (const Line2RzInterfaceQuadraturePoint& point : geometry.points) {
        const InterfaceAdQuadratureValue value =
            evaluate_quadrature_value(point, state, properties_);

        for (std::size_t node = 0; node < line2_interface_side_node_count;
             ++node) {
            const double weighted_shape =
                point.weighted_measure * point.shape[node];

            residual[node] += weighted_shape * value.heat_flux;
            residual[2 + node] -= weighted_shape * value.heat_flux;

            residual[4 + node] += weighted_shape * value.pressure;
            residual[6 + node] -= weighted_shape * value.pressure;
        }
    }
}

} // namespace fuelsim
