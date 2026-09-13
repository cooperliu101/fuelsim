#pragma once
#include "element_types.hpp"
#include "material.hpp"
#include <array>
#include <cstddef>

namespace fuelsim {
using Cax4LocalDofs = std::array<std::size_t, 12>;
constexpr std::size_t cax4_local_dof_count = 12;
constexpr std::size_t cax4_local_jacobian_size = cax4_local_dof_count * cax4_local_dof_count;
using Cax4LocalValues = std::array<double, cax4_local_dof_count>;
using Cax4LocalResidual = std::array<double, cax4_local_dof_count>;
using Cax4LocalJacobian = std::array<double, cax4_local_jacobian_size>;
using Cax4LocalAdValues = std::array<adlite::Scalar, cax4_local_dof_count>;
constexpr std::size_t quad4_node_count = 4;
constexpr std::size_t quad4_local_dof_count = cax4_local_dof_count;
using Quad4Coordinates = std::array<RzPoint, quad4_node_count>;

struct RzQuadraturePoint final {
    std::array<double, quad4_node_count> shape, gradient_r, gradient_z;
    double radius, axial_coordinate, weighted_measure;
};

struct Quad4RzGeometry final {
    std::array<RzQuadraturePoint, 4> points;
    Quad4Coordinates coordinates{};
};

using Quad4MaterialHistory = std::array<MaterialPointState, 4>;

} // namespace fuelsim

namespace fuelsim::elements {
// Field-major local order: [T0..T3, ur0..ur3, uz0..uz3]. All inputs are borrowed
// for this call only. Null history selects a steady thermoelastic evaluation.
struct Cax4Input final {
    const IsotropicThermoelasticMaterial& material;
    const Quad4RzGeometry& geometry;
    const Cax4LocalValues& state;
    const Cax4LocalValues& committed_state;
    const Quad4MaterialHistory* committed_history = nullptr;
    double time_step = 0.0;
    double time = 0.0;
    double volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 600.0;
    std::array<double, 3> body_acceleration{};
};

struct Cax4Result final {
    std::array<AxisymmetricStressValues, 4> stress{};
    Cax4LocalResidual residual{};
    Cax4LocalJacobian jacobian{};   // Row-major d(residual)/d(state); zero when not requested.
    Quad4MaterialHistory history{}; // Trial only; the caller owns acceptance/rollback.
    double stored_heat_rate = 0.0, generated_heat_rate = 0.0;
};

} // namespace fuelsim::elements
