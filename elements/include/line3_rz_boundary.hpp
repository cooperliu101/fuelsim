#pragma once
#include "rz_quad4.hpp"
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
