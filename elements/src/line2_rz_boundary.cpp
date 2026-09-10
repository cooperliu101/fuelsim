#include "line2_rz_boundary.hpp"
#include "boundary_types.hpp"
#include "detail/line2_rz_geometry.hpp"
#include "detail/line2_rz_local_system.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss = 0.577350269189625764509148780501957456;
using line2_rz_detail::validate_line;

void validate_edge(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes,
    const char* name) {
    if (local_nodes[0] >= 4 || local_nodes[1] >= 4 || local_nodes[0] == local_nodes[1])
        throw std::invalid_argument(std::string(name) + " requires two distinct Quad4 local nodes");
    validate_line(coordinates, name);
}

void displaced_edge_coordinates(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes,
    const Cax4LocalAdValues& state,
    bool use_displaced_geometry,
    std::array<adlite::Scalar, 2>& radius,
    std::array<adlite::Scalar, 2>& axial) {
    for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
        const std::size_t local = local_nodes[edge_node];
        radius[edge_node] = coordinates[edge_node].r;
        axial[edge_node] = coordinates[edge_node].z;
        if (use_displaced_geometry) {
            radius[edge_node] += state[4 + local];
            axial[edge_node] += state[8 + local];
        }
    }
}
} // namespace

Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes) {
    validate_edge(coordinates, local_nodes, "Boundary edge");
    return {coordinates, local_nodes};
}

namespace {
void compute_line2_rz_mechanical_residual_ad(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const Cax4LocalAdValues& state,
    Cax4LocalAdValues& residual) {
    residual.fill(adlite::Scalar(0.0));
    std::array<adlite::Scalar, 2> radius{};
    std::array<adlite::Scalar, 2> axial{};
    displaced_edge_coordinates(geometry.coordinates,
        geometry.local_nodes,
        state,
        data.use_displaced_geometry,
        radius,
        axial);
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const double xi : locations) {
        const std::array<double, 2> shape = {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
        const adlite::Scalar current_radius = shape[0] * radius[0] + shape[1] * radius[1];
        const adlite::Scalar dr_dxi = 0.5 * (radius[1] - radius[0]), dz_dxi = 0.5 * (axial[1] - axial[0]);
        const adlite::Scalar measure = data.kind == Line2RzBoundaryKind::pressure
                                           ? 2.0 * pi * current_radius
                                           : 2.0 * pi * current_radius * adlite::hypot(dr_dxi, dz_dxi);
        if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
            throw std::domain_error("Mechanical boundary current measure must be finite and positive");
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
            const std::size_t local = geometry.local_nodes[edge_node];
            if (data.kind == Line2RzBoundaryKind::pressure) {
                residual[4 + local] += measure * data.load * shape[edge_node] * dz_dxi;
                residual[8 + local] -= measure * data.load * shape[edge_node] * dr_dxi;
            } else {
                const std::size_t offset = data.component == TractionComponent::radial ? 4 : 8;
                residual[offset + local] -= measure * data.load * shape[edge_node];
            }
        }
    }
}

void compute_line2_rz_convection_residual_ad(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const Cax4LocalAdValues& state,
    Cax4LocalAdValues& residual) {
    residual.fill(adlite::Scalar(0.0));
    const double dr = geometry.coordinates[1].r - geometry.coordinates[0].r,
                 dz = geometry.coordinates[1].z - geometry.coordinates[0].z, line_jacobian = 0.5 * std::hypot(dr, dz);
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const double xi : locations) {
        const std::array<double, 2> shape = {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
        const double radius = shape[0] * geometry.coordinates[0].r + shape[1] * geometry.coordinates[1].r;
        adlite::Scalar temperature = 0.0;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node)
            temperature += shape[edge_node] * state[geometry.local_nodes[edge_node]];
        const adlite::Scalar heat_flux = data.load * (temperature - data.ambient);
        const double measure = 2.0 * pi * radius * line_jacobian;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node)
            residual[geometry.local_nodes[edge_node]] += measure * shape[edge_node] * heat_flux;
    }
}
} // namespace

Cax4LocalResidual compute_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const Cax4LocalValues& state,
    Cax4LocalJacobian* jacobian) {
    const Cax4LocalAdValues ad_state = line2_rz_detail::ad_state(state, jacobian != nullptr);
    Cax4LocalAdValues ad_residual{};
    if (data.kind == Line2RzBoundaryKind::convection)
        compute_line2_rz_convection_residual_ad(data, geometry, ad_state, ad_residual);
    else
        compute_line2_rz_mechanical_residual_ad(data, geometry, ad_state, ad_residual);
    return line2_rz_detail::values(ad_state, ad_residual, jacobian);
}

} // namespace fuelsim
