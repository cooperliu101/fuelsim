#pragma once
#include "thermal_types.hpp"

namespace fuelsim::elements {
struct ThermalSurfacePoint final {
    ThermalPoint point;
    std::array<double, 3> tangent_xi{}, tangent_eta{};
};

ThermalSurfacePoint
thermal_surface_point(const std::vector<CartesianPoint3>& coordinates, bool axisymmetric, double xi, double eta);
// Line nodes are endpoints followed by the midpoint. Face nodes are corners followed by edge midpoints.
ThermalGeometry make_thermal_boundary_geometry(const std::vector<CartesianPoint3>& coordinates, bool axisymmetric);
ThermalResult evaluate_thermal_boundary(const ThermalGeometry& geometry,
    const std::vector<double>& temperature,
    double inward_flux,
    double coefficient,
    double ambient,
    bool jacobian);
} // namespace fuelsim::elements
