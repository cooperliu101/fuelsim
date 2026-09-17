#pragma once
#include "contact_types.hpp"
#include <array>
#include <vector>

namespace fuelsim::elements {
// [Ts0,Ts1,Tp0,Tp1,usx(3),usy(3),upx(3),upy(3),qs(3),qp(3)].
using Line3PlaneValues = std::array<double, 22>;

struct Line3PlaneContactInput final {
    std::array<std::array<double, 2>, 3> secondary, primary;
    std::array<double, 2> secondary_reference, primary_reference;
    double secondary_thickness, primary_thickness;
    bool secondary_finite, primary_finite;
    const Line3PlaneValues& state;
    double coordinate, weight;
    GapHeatProperties heat;
    double penalty;
    bool thermal, mechanical;
    // Two quadratic endpoint projection equations have at most four roots.
    // Five fixed interval slots keep assembly and output sizes state independent.
    std::size_t segment = 0;
};

struct Line3PlaneContactResult final {
    Line3PlaneValues residual{};
    std::array<double, 484> jacobian{};
    bool projected = false;
    double gap = 0.0, pressure = 0.0, area = 0.0, heat_rate = 0.0, force = 0.0;
    double distance_squared = 0.0, primary_coordinate = 0.0;
    double secondary_coordinate = 0.0, interval_begin = 0.0, interval_end = 0.0;
};

Line3PlaneContactResult evaluate_line3_plane_contact(const Line3PlaneContactInput& input, bool jacobian);

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
