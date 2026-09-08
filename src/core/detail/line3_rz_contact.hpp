#pragma once
#include "fuelsim/core/contact.hpp"
#include <vector>

namespace fuelsim::rz8 {
std::pair<bool, double> project_line3(const std::array<RzPoint, 3>& primary, RzPoint point);

struct Line3ContactGeometry final {
    std::array<RzPoint, 3> secondary, primary;
    bool primary_first = false, primary_last = false;
    double coordinate = 0, weight = 0;
    std::size_t secondary_node = 0;
    bool mechanical = false;
    bool nodal_heat = false;
};

struct Line3ContactResult final {
    std::array<double, 16> residual{};
    std::array<double, 256> jacobian{};
    ContactPointValue mechanical{};
    HeatQuadratureValue thermal{};
    double primary_coordinate = 0;
};

Line3ContactResult compute_line3_contact(const Line3ContactGeometry& geometry,
    const GapHeatProperties& heat,
    const NormalContactProperties& mechanical,
    const std::vector<double>& state,
    const std::vector<double>& old,
    const ContactPointHistory& history,
    bool jacobian);
} // namespace fuelsim::rz8
