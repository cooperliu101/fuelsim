#pragma once
#include "core/mesh.hpp"
#include "material.hpp"
#include "rz_geometry.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>

namespace fuelsim {
using LocalDofs = std::array<std::size_t, local_dof_count>;

struct Quad4RzData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source = 0.0, time = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    RzElementFormulation element_formulation = RzElementFormulation::quad4;
    double initial_temperature = 600.0;
};

LocalResidual compute_quad4_rz_thermoelastic(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian = nullptr);
std::array<AxisymmetricStressValues, 4> compute_quad4_rz_thermoelastic_stress(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state);
LocalResidual compute_quad4_rz_transient(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& current_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material,
    double time_step,
    LocalJacobian* jacobian = nullptr,
    bool include_thermal_time_term = true);
Quad4MaterialHistory compute_quad4_rz_transient_update(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& converged_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material,
    double time_step);

struct Line2RzBoundaryGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};
enum class TractionComponent { radial, axial };
Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes);
enum class Line2RzBoundaryKind { pressure, traction, convection };

struct Line2RzBoundaryData final {
    Line2RzBoundaryKind kind;
    TractionComponent component;
    double load, ambient;
    bool use_displaced_geometry;
};

LocalResidual compute_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian = nullptr);
} // namespace fuelsim
