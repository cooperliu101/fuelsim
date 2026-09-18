#pragma once
#include "contact_types.hpp"
#include <array>
#include <vector>

namespace fuelsim::elements {
struct Line3PlaneContactResult final {
    bool projected = false;
    double gap = 0.0, pressure = 0.0, area = 0.0, heat_rate = 0.0, force = 0.0;
};

struct PlaneAveragingEdge final {
    std::array<std::size_t, 3> nodes;
    std::size_t local_node;
    double thickness;
};

struct PlaneAveragedContactGeometry final {
    std::vector<std::array<double, 2>> coordinates;
    std::vector<PlaneAveragingEdge> secondary;
    std::vector<std::array<std::size_t, 3>> primary;
    double penalty = 0.0;
    // Nonempty for a thermal corner-neighborhood constraint. State temperatures
    // follow the displacement pairs, in this list's order; midsides have no T.
    std::vector<std::size_t> temperature_nodes;
    GapHeatProperties heat{};
};

struct PlaneSurfaceContactResult final {
    std::vector<double> residual, jacobian;
    Line3PlaneContactResult point;
};

// Interface-local displacement values are [ux0,uy0,ux1,uy1,...].
PlaneSurfaceContactResult evaluate_plane_averaged_contact(const PlaneAveragedContactGeometry& geometry,
    const std::vector<double>& state,
    bool jacobian);
} // namespace fuelsim::elements
