#pragma once
#include "contact_types.hpp"
#include "element_types.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <vector>

namespace fuelsim::elements {
enum class Line3RzBoundaryKind { pressure, traction, heat_flux, convection };

struct Line3RzBoundaryData final {
    Line3RzBoundaryKind kind;
    TractionComponent component;
    double load, ambient;
    bool use_displaced_geometry;
};

struct Line3RzBoundaryResult final {
    std::array<double, 8> residual{};
    std::array<double, 64> jacobian{};
};

// Local order: [T0,T1,ur0,ur1,ur_mid,uz0,uz1,uz_mid]. The first two coordinates
// are endpoints, the third is the midside node. load is the convection
// coefficient for convection; all controls are evaluated by the caller.
Line3RzBoundaryResult compute_line3_rz_boundary(const Line3RzBoundaryData& data,
    const std::array<RzPoint, 3>& coordinates,
    const std::vector<double>& state,
    bool jacobian = false);
} // namespace fuelsim::elements

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
