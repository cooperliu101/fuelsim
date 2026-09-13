#pragma once
#include "element_types.hpp"
#include "material.hpp"
#include <array>
#include <cstddef>

namespace fuelsim {
constexpr std::size_t cax2t_gps_local_dof_count = 6;
constexpr std::size_t cax2t_gps_material_point_count = 2;
using Cax2tGpsLocalValues = std::array<double, cax2t_gps_local_dof_count>;
using Cax2tGpsLocalResidual = std::array<double, cax2t_gps_local_dof_count>;
using Cax2tGpsLocalJacobian = std::array<double, cax2t_gps_local_dof_count * cax2t_gps_local_dof_count>;
using Cax2tGpsMaterialHistory = std::array<MaterialPointState, cax2t_gps_material_point_count>;

struct Cax2tGpsGeometry final {
    std::array<double, 2> radii;
    double z_lower, z_upper;
};
} // namespace fuelsim

namespace fuelsim::elements {
// This is a fuelsim radial generalized-plane-strain model, not an Abaqus element alias.
Cax2tGpsGeometry make_cax2t_gps_geometry(const std::array<double, 2>& radii, double z_lower, double z_upper);

// Local order: [T0, T1, ur0, ur1, w_lower, w_upper]. The last two values belong
// to the axial end sections of a physical body, not to the two radial nodes.
// All radial elements in the same body and axial segment must share them.
// The caller imposes ur=0 at an axis node and owns all acceptance/rollback.
struct Cax2tGpsInput final {
    const IsotropicThermoelasticMaterial& material;
    const Cax2tGpsGeometry& geometry;
    const Cax2tGpsLocalValues& state;
    const Cax2tGpsLocalValues& committed_state;
    const Cax2tGpsMaterialHistory* committed_history = nullptr;
    double time_step = 0.0;
    double time = 0.0;
    double volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 600.0;
};

struct Cax2tGpsResult final {
    std::array<AxisymmetricStressValues, cax2t_gps_material_point_count> stress{};
    Cax2tGpsLocalResidual residual{};
    Cax2tGpsLocalJacobian jacobian{};  // Row-major d(residual)/d(state).
    Cax2tGpsMaterialHistory history{}; // Trial only; input histories are never modified.
    double stored_heat_rate = 0.0, generated_heat_rate = 0.0;
};

// Two radial Gauss points use paired radial-endpoint temperatures for elastic and
// inelastic properties. Eigenstrain uses the arithmetic mean of the two radial
// temperatures. Conductivity uses paired endpoint temperatures; heat capacity is
// row-sum lumped with nodal properties and rates. Mechanical hoop strain uses a reference-volume mean;
// finite strain averages the hoop stretch before its midpoint increment. Mechanical
// point weights and conduction weights scale with the whole-element volume ratio.
// Finite conduction gradients use the incremental midpoint configuration; heat
// capacity uses fixed initial mass and current specific heat; source measures
// retain the actual pointwise current geometry.
Cax2tGpsResult evaluate_cax2t_gps(const Cax2tGpsInput& input, ElementRequest request = {});
} // namespace fuelsim::elements
