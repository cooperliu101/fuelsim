#pragma once
#include "element_types.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {
enum class GapHeatConductanceLaw { gas_gap, affine };

struct GapHeatProperties final {
    double gap_conductivity, minimum_gap;
    GapHeatConductanceLaw law = GapHeatConductanceLaw::gas_gap;
    double conductance = 0.0, clearance_derivative = 0.0, pressure_derivative = 0.0, temperature_derivative = 0.0,
           reference_temperature = 0.0, contact_penalty = 0.0;
};

struct HeatQuadratureValue final {
    bool projected;
    double gap;
    double heat_flux, weighted_measure;
};

struct ContactProjectionValue final {
    bool projected;
    double gap;
};

struct NormalContactProperties final {
    double penalty;
    double friction_coefficient = 0.0;
    bool augmented_lagrangian = false;
    double maximum_elastic_slip = 0.0;
};

struct ContactPointHistory final {
    double elastic_tangential_slip = 0.0;
    bool sliding = false;
    double normal_multiplier = 0.0;
    std::array<double, 3> cartesian_elastic_tangential_slip{};
    // Abaqus CSLIP-compatible accumulated relative motion: add only while contact is active, transport it with the
    // current tangent basis, and add no relative motion while the projected contact is open.
    std::array<double, 3> cartesian_total_tangential_slip{};
    bool cartesian_tangent_basis_initialized = false;
    std::array<double, 3> cartesian_contact_normal{}, cartesian_contact_tangent_first{};
    // Signed accumulated relative motion in the current RZ tangent direction; frozen while open.
    double total_tangential_slip = 0.0;
};

struct ContactPointValue final {
    bool projected;
    double gap, pressure, tributary_area, tributary_length, contact_force, tangential_traction, tangential_force,
        elastic_tangential_slip;
    bool sliding;
    double total_tangential_slip = 0.0;
};

struct CartesianHeatQuadratureValue final {
    bool projected;
    double gap, heat_flux, weighted_measure;
};

struct CartesianContactPointValue final {
    bool projected;
    double gap, pressure, tributary_area, contact_force, tangential_traction, tangential_force, friction_dissipation;
    std::array<double, 3> normal, tangent_first, tangential_traction_vector, tangential_slip, elastic_tangential_slip;
    bool sliding;
};
} // namespace fuelsim

namespace fuelsim {
enum class Quad8NodalAreaRule {
    positive_lumped,
    consistent_shape,
};
} // namespace fuelsim
