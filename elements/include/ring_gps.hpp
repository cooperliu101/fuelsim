#pragma once
#include "contact_types.hpp"
#include <array>
#include <cstddef>

namespace fuelsim {
inline constexpr std::size_t ring_gps_local_dof_count = 8;
inline constexpr std::size_t ring_gps_quadrature_point_count = 2;
using RingGpsLocalValues = std::array<double, ring_gps_local_dof_count>;
using RingGpsLocalResidual = RingGpsLocalValues;
using RingGpsLocalJacobian = std::array<double, ring_gps_local_dof_count * ring_gps_local_dof_count>;
using RingGpsHistory = std::array<ContactPointHistory, ring_gps_quadrature_point_count>;

struct RingGpsGeometry final {
    double primary_reference_radius, secondary_reference_radius, z_lower, z_upper;
};

struct RingGpsPointValue final {
    double gap = 0.0, pressure = 0.0, heat_flux = 0.0, weighted_measure = 0.0;
    double signed_tangential_traction = 0.0, elastic_tangential_slip = 0.0, total_tangential_slip = 0.0;
    double friction_dissipation = 0.0;
    bool sliding = false;
};
} // namespace fuelsim

namespace fuelsim::elements {
// Local order: [Tp, Ts, urp, urs, wp_lower, wp_upper, ws_lower, ws_upper].
// Primary is the inner surface of the outer cylinder. The caller supplies the
// common axial interval and owns all history acceptance and failure recovery.
// Axial contact pairing remains fixed in the reference interval (small sliding).
// Finite strain uses the current secondary radius and axial length for the common
// surface measure; it does not enable finite sliding between different slices.
struct RingGpsInput final {
    const RingGpsGeometry& geometry;
    const RingGpsLocalValues& state;
    const RingGpsLocalValues& committed_state;
    const RingGpsHistory& committed_history;
    bool thermal = false;
    bool mechanical = false;
    GapHeatProperties heat{0.0, 1.0};
    NormalContactProperties normal{0.0};
    StrainFormulation strain_formulation = StrainFormulation::small;
};

struct RingGpsResult final {
    RingGpsLocalResidual residual{};
    RingGpsLocalJacobian jacobian{}; // Row-major; zero when not requested.
    RingGpsHistory history{};        // Trial only; zero when not requested.
    std::array<RingGpsPointValue, ring_gps_quadrature_point_count> points{};
    // Heat rate is positive from secondary to primary. Tangential force is the
    // signed internal resistance conjugate to secondary-minus-primary motion.
    double heat_rate = 0.0, normal_force = 0.0, tangential_force = 0.0, friction_dissipation = 0.0;
};

RingGpsResult evaluate_ring_gps(const RingGpsInput& input, ElementRequest request = {});
} // namespace fuelsim::elements
