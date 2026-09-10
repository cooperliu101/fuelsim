#pragma once

namespace fuelsim {
enum class CartesianTractionComponent { x, y, z };

enum class Quad4FaceBoundaryKind { pressure, traction, surface_heat_flux, convection };

struct Quad4FaceBoundaryData final {
    Quad4FaceBoundaryKind kind;
    CartesianTractionComponent component;
    double load, ambient_temperature;
    bool use_displaced_geometry = false;
};

enum class TractionComponent { radial, axial };
} // namespace fuelsim
