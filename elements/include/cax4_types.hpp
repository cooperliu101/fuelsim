#pragma once
#include "axisymmetric_geometry.hpp"

namespace fuelsim::elements {
// Field-major local order: [T0..T3, ur0..ur3, uz0..uz3]. All inputs are borrowed
// for this call only. Null history selects a steady thermoelastic evaluation.
struct Cax4Input final {
    const IsotropicThermoelasticMaterial& material;
    const Quad4RzGeometry& geometry;
    const LocalValues& state;
    const LocalValues& committed_state;
    const Quad4MaterialHistory* committed_history = nullptr;
    double time_step = 0.0;
    double time = 0.0;
    double volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 600.0;
};

struct Cax4Result final {
    std::array<AxisymmetricStressValues, 4> stress{};
    LocalResidual residual{};
    LocalJacobian jacobian{};       // Row-major d(residual)/d(state); zero when not requested.
    Quad4MaterialHistory history{}; // Trial only; the caller owns acceptance/rollback.
    double stored_heat_rate = 0.0, generated_heat_rate = 0.0;
};

} // namespace fuelsim::elements
