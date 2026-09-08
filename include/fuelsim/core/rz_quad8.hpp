#pragma once
#include "fuelsim/core/rz_quad4.hpp"

namespace fuelsim {
using Quad8RzCoordinates = std::array<RzPoint, 8>;
using Quad8RzValues = std::array<double, 20>;
using Quad8RzJacobian = std::array<double, 400>;
using Quad8MaterialHistory = std::array<MaterialPointState, 9>;

struct Quad8RzPoint final {
    std::array<double, 8> shape{}, gradient_r{}, gradient_z{};
    std::array<double, 4> temperature_shape{}, temperature_gradient_r{}, temperature_gradient_z{};
    std::array<double, 4> source_gradient_r{}, source_gradient_z{};
    double source_radius = 0.0, source_measure = 0.0;
    double radius = 0.0, axial_coordinate = 0.0, weighted_measure = 0.0;
};

struct Quad8RzGeometry final {
    Quad8RzCoordinates coordinates;
    std::array<Quad8RzPoint, 9> points;
    std::size_t point_count = 9;
};

struct Quad8RzResult final {
    Quad8RzValues residual{};
    Quad8RzJacobian jacobian{};
    Quad8MaterialHistory history{};
    double generated_heat_rate = 0.0, stored_heat_rate = 0.0;
};

Quad8RzPoint evaluate_quad8_rz_point(const Quad8RzCoordinates& coordinates, double xi, double eta, double weight);
Quad8RzGeometry make_quad8_rz_geometry(const Quad8RzCoordinates& coordinates,
    RzElementFormulation formulation = RzElementFormulation::cax8t);
Quad8RzResult compute_quad8_rz(const Quad4RzData& data,
    const Quad8RzGeometry& geometry,
    const Quad8RzValues& state,
    const Quad8RzValues& committed,
    const Quad8MaterialHistory* history,
    double time_step,
    bool jacobian,
    bool thermal_time = true);
} // namespace fuelsim
