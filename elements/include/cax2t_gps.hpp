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
};

struct Cax2tGpsResult final {
    std::array<AxisymmetricStressValues, cax2t_gps_material_point_count> stress{};
    Cax2tGpsLocalResidual residual{};
    Cax2tGpsLocalJacobian jacobian{};  // Row-major d(residual)/d(state).
    Cax2tGpsMaterialHistory history{}; // Trial only; input histories are never modified.
    double stored_heat_rate = 0.0, generated_heat_rate = 0.0;
};

// Two radial Gauss points, consistent heat capacity, and integration-point
// material temperature, except eigenstrain temperature uses the arithmetic mean
// of the two radial nodes. Finite strain uses diagonal radial/axial/hoop stretches,
// midpoint strain increments, and current-configuration mechanical and thermal measures.
Cax2tGpsResult evaluate_cax2t_gps(const Cax2tGpsInput& input, ElementRequest request = {});
} // namespace fuelsim::elements
