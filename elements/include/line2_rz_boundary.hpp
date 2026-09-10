#pragma once
#include "boundary_types.hpp"
#include "cax4_types.hpp"

namespace fuelsim {
struct Line2RzBoundaryGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};

Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes);

enum class Line2RzBoundaryKind { pressure, traction, convection };

struct Line2RzBoundaryData final {
    Line2RzBoundaryKind kind;
    TractionComponent component;
    double load, ambient;
    bool use_displaced_geometry;
};

Cax4LocalResidual compute_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const Cax4LocalValues& state,
    Cax4LocalJacobian* jacobian = nullptr);
} // namespace fuelsim
