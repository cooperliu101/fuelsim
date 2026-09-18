#pragma once
#include "contact_types.hpp"
#include "thermal_boundary.hpp"

namespace fuelsim::elements {
struct ThermalInterfacePoint final {
    // Each primary temperature entry identifies its candidate face and local node.
    std::vector<std::array<std::size_t, 2>> primary_nodes;
    std::vector<double> secondary_shape, primary_shape;
    double gap = 0.0, measure = 0.0;
};

std::vector<ThermalInterfacePoint> make_thermal_interface_points(const std::vector<CartesianPoint3>& secondary,
    const std::vector<std::vector<CartesianPoint3>>& primary,
    bool axisymmetric);
ThermalResult evaluate_thermal_interface(const ThermalInterfacePoint& point,
    const GapHeatProperties& material,
    const std::vector<double>& temperatures,
    bool jacobian);
} // namespace fuelsim::elements
