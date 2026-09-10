#pragma once
#include "element_types.hpp"
#include "material.hpp"
#include "material_types.hpp"

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

Quad8RzPoint evaluate_quad8_rz_point(const Quad8RzCoordinates& coordinates, double xi, double eta, double weight);

} // namespace fuelsim

namespace fuelsim::elements {
struct Cax8Input final {
    const IsotropicThermoelasticMaterial& material;
    const Quad8RzGeometry& geometry;
    const Quad8RzValues& state;
    const Quad8RzValues& committed_state;
    const Quad8MaterialHistory* committed_history = nullptr;
    double time_step = 0.0;
    double time = 0.0;
    double volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 600.0;
};

struct Cax8Result final {
    std::array<AxisymmetricStressValues, 9> stress{};
    Quad8RzValues residual{};
    Quad8RzJacobian jacobian{};
    Quad8MaterialHistory history{};
    double generated_heat_rate = 0.0, stored_heat_rate = 0.0;
};
} // namespace fuelsim::elements
