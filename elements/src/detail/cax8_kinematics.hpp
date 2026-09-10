#pragma once
#include "cax8_types.hpp"

namespace fuelsim::cax8_detail {
struct PointKinematics final {
    std::array<adlite::Scalar, 6> active;
    std::array<Quad8RzValues, 6> chain{};
    std::array<double, 6> old{};
    std::array<adlite::Scalar, 4> strain;
    AxisymmetricRotation rotation;
    adlite::Scalar frr, frz, fzr, fzz, det, radius, measure;
};

struct MechanicalGradient final {
    adlite::Scalar radial, axial, hoop;
};

struct ThermalGeometry final {
    std::array<adlite::Scalar, 4> radial, axial;
};

struct SourceGeometry final {
    double measure = 0.0;
    Quad8RzValues derivative{};
};

PointKinematics evaluate_kinematics(const Quad8RzPoint& point,
    const Quad8RzValues& state,
    const Quad8RzValues& committed,
    bool finite,
    bool jacobian);
MechanicalGradient
mechanical_gradient(const Quad8RzPoint& point, const PointKinematics& kinematics, std::size_t node, bool finite);
ThermalGeometry thermal_geometry(const Quad8RzPoint& point, const PointKinematics& kinematics, bool finite);
SourceGeometry source_geometry(const Quad8RzPoint& point, const Quad8RzValues& state, bool finite);
} // namespace fuelsim::cax8_detail
