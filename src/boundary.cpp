#include "fuelsim/boundary.hpp"

#include <adlite/adlite.hpp>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss = 0.577350269189625764509148780501957456;

void validate_properties(const ConvectionProperties& properties) {
    if (!std::isfinite(properties.heat_transfer_coefficient) ||
        properties.heat_transfer_coefficient < 0.0)
        throw std::invalid_argument(
            "Convection heat-transfer coefficient must be nonnegative");
    if (!std::isfinite(properties.ambient_temperature) ||
        !(properties.ambient_temperature > 0.0))
        throw std::invalid_argument(
            "Convection ambient temperature must be positive");
}

void validate_properties(const PressureProperties& properties) {
    if (!std::isfinite(properties.pressure) || properties.pressure < 0.0)
        throw std::invalid_argument(
            "Pressure magnitude must be finite and nonnegative");
}

void validate_edge(const std::array<RzPoint, 2>& coordinates,
                   const std::array<std::size_t, 2>& local_nodes,
                   const char* name) {
    if (local_nodes[0] >= 4 || local_nodes[1] >= 4 ||
        local_nodes[0] == local_nodes[1])
        throw std::invalid_argument(std::string(name) +
                                    " requires two distinct Quad4 local nodes");
    for (const RzPoint& point : coordinates) {
        if (!std::isfinite(point.r) || !std::isfinite(point.z) || point.r < 0.0)
            throw std::invalid_argument(
                std::string(name) +
                " coordinates must be finite with nonnegative radius");
    }
    const double dr = coordinates[1].r - coordinates[0].r;
    const double dz = coordinates[1].z - coordinates[0].z;
    if (!(std::hypot(dr, dz) > 0.0))
        throw std::invalid_argument(std::string(name) +
                                    " length must be positive");
}

} // namespace

Line2RzConvectionGeometry make_line2_rz_convection_geometry(
    const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes) {
    validate_edge(coordinates, local_nodes, "Convection edge");
    return {coordinates, local_nodes};
}

Line2RzPressureGeometry make_line2_rz_pressure_geometry(
    const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes) {
    validate_edge(coordinates, local_nodes, "Pressure edge");
    return {coordinates, local_nodes};
}

Line2RzPressureKernel::Line2RzPressureKernel(
    PressureProperties properties)
    : _properties(properties) {
    validate_properties(_properties);
}

const PressureProperties& Line2RzPressureKernel::properties() const noexcept {
    return _properties;
}

void Line2RzPressureKernel::set_properties(PressureProperties properties) {
    validate_properties(properties);
    _properties = properties;
}

void Line2RzPressureKernel::residual_ad(
    const Line2RzPressureGeometry& geometry, const LocalAdValues& state,
    LocalAdValues& residual) const {
    residual.fill(adlite::Scalar(0.0));
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const double xi : locations) {
        const std::array<double, 2> shape = {0.5 * (1.0 - xi),
                                             0.5 * (1.0 + xi)};
        std::array<adlite::Scalar, 2> radius{};
        std::array<adlite::Scalar, 2> axial{};
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
            const std::size_t local = geometry.local_nodes[edge_node];
            radius[edge_node] = geometry.coordinates[edge_node].r;
            axial[edge_node] = geometry.coordinates[edge_node].z;
            if (_properties.use_displaced_geometry) {
                radius[edge_node] += state[4 + local];
                axial[edge_node] += state[8 + local];
            }
        }
        const adlite::Scalar current_radius =
            shape[0] * radius[0] + shape[1] * radius[1];
        if (!std::isfinite(current_radius.value()) ||
            !(current_radius.value() > 0.0))
            throw std::domain_error(
                "Pressure edge current radius must be finite and positive");

        // Region boundary edges retain the parent Quad4 counter-clockwise
        // ordering.  (dz/dxi, -dr/dxi) is therefore outward normal times the
        // line Jacobian, so no normalization is needed.
        const adlite::Scalar outward_r = 0.5 * (axial[1] - axial[0]);
        const adlite::Scalar outward_z = -0.5 * (radius[1] - radius[0]);
        const adlite::Scalar measure = 2.0 * pi * current_radius;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
            const std::size_t local = geometry.local_nodes[edge_node];
            residual[4 + local] += measure * _properties.pressure *
                                   shape[edge_node] * outward_r;
            residual[8 + local] += measure * _properties.pressure *
                                   shape[edge_node] * outward_z;
        }
    }
}

LocalResidual Line2RzPressureKernel::residual(
    const Line2RzPressureGeometry& geometry, const LocalValues& state) const {
    LocalAdValues ad_state{};
    for (std::size_t dof = 0; dof < local_dof_count; ++dof)
        ad_state[dof] = state[dof];
    LocalAdValues ad_residual{};
    residual_ad(geometry, ad_state, ad_residual);
    LocalResidual result{};
    for (std::size_t dof = 0; dof < local_dof_count; ++dof)
        result[dof] = ad_residual[dof].value();
    return result;
}

LocalSystem Line2RzPressureKernel::linearize(
    const Line2RzPressureGeometry& geometry, const LocalValues& state) const {
    LocalAdValues ad_state{};
    adlite::seed_identity(state.data(), state.size(), ad_state.data());
    LocalAdValues ad_residual{};
    residual_ad(geometry, ad_state, ad_residual);
    LocalSystem result{};
    adlite::extract_jacobian(ad_residual.data(), ad_residual.size(),
                             ad_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

Line2RzConvectionKernel::Line2RzConvectionKernel(
    ConvectionProperties properties)
    : _properties(properties) {
    validate_properties(_properties);
}

const ConvectionProperties&
Line2RzConvectionKernel::properties() const noexcept {
    return _properties;
}

void Line2RzConvectionKernel::set_properties(ConvectionProperties properties) {
    validate_properties(properties);
    _properties = properties;
}

void Line2RzConvectionKernel::residual_ad(
    const Line2RzConvectionGeometry& geometry, const LocalAdValues& state,
    LocalAdValues& residual) const {
    residual.fill(adlite::Scalar(0.0));
    const double dr = geometry.coordinates[1].r - geometry.coordinates[0].r;
    const double dz = geometry.coordinates[1].z - geometry.coordinates[0].z;
    const double line_jacobian = 0.5 * std::hypot(dr, dz);
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const double xi : locations) {
        const std::array<double, 2> shape = {0.5 * (1.0 - xi),
                                             0.5 * (1.0 + xi)};
        const double radius = shape[0] * geometry.coordinates[0].r +
                              shape[1] * geometry.coordinates[1].r;
        adlite::Scalar temperature = 0.0;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node)
            temperature +=
                shape[edge_node] * state[geometry.local_nodes[edge_node]];
        const adlite::Scalar heat_flux =
            _properties.heat_transfer_coefficient *
            (temperature - _properties.ambient_temperature);
        const double measure = 2.0 * pi * radius * line_jacobian;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node)
            residual[geometry.local_nodes[edge_node]] +=
                measure * shape[edge_node] * heat_flux;
    }
}

LocalResidual
Line2RzConvectionKernel::residual(const Line2RzConvectionGeometry& geometry,
                                  const LocalValues& state) const {
    LocalAdValues ad_state{};
    for (std::size_t dof = 0; dof < local_dof_count; ++dof)
        ad_state[dof] = state[dof];
    LocalAdValues ad_residual{};
    residual_ad(geometry, ad_state, ad_residual);
    LocalResidual result{};
    for (std::size_t dof = 0; dof < local_dof_count; ++dof)
        result[dof] = ad_residual[dof].value();
    return result;
}

LocalSystem
Line2RzConvectionKernel::linearize(const Line2RzConvectionGeometry& geometry,
                                   const LocalValues& state) const {
    LocalAdValues ad_state{};
    adlite::seed_identity(state.data(), state.size(), ad_state.data());
    LocalAdValues ad_residual{};
    residual_ad(geometry, ad_state, ad_residual);
    LocalSystem result{};
    adlite::extract_jacobian(ad_residual.data(), ad_residual.size(),
                             ad_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

} // namespace fuelsim
