#pragma once
#include "fuelsim/core/contact_area_rule.hpp"
#include "fuelsim/core/kinematics.hpp"
#include "fuelsim/core/material.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace fuelsim {
enum class Field {
    temperature,
    radial_displacement,
    axial_displacement,
    displacement_x,
    displacement_y,
    displacement_z,
};

class PiecewiseLinearTimeTable final {
  public:
    PiecewiseLinearTimeTable(std::string name, std::vector<double> times, std::vector<double> values);

    const std::string& name() const noexcept { return _name; }

    const std::vector<double>& times() const noexcept { return _times; }

    const std::vector<double>& values() const noexcept { return _values; }

    double value(double time) const;

  private:
    std::string _name;
    std::vector<double> _times, _values;
};

struct RegionDefinition final {
    std::string name, block;
    ThermoelasticProperties material;
    double volumetric_heat_source, initial_temperature;
    std::int64_t block_id = -1;
    std::string heat_source_function{};
    StrainFormulation strain_formulation = StrainFormulation::small;
};
enum class MechanicalContactFormulation {
    penalty,
    augmented_lagrangian,
};

struct ContactDefinition final {
    std::string name, primary, secondary;
    bool thermal, mechanical;
    double gap_conductivity, minimum_gap, penalty;
    double friction_coefficient = 0.0;
    bool automatic_penalty = false;
    double penalty_factor = 1.0;
    MechanicalContactFormulation mechanical_formulation = MechanicalContactFormulation::penalty;
    Quad8NodalAreaRule quad8_nodal_area_rule = Quad8NodalAreaRule::positive_lumped;
    double penetration_tolerance = 1.0e-8;
    std::size_t maximum_augmented_iterations = 20;
};

struct AugmentedContactUpdate final {
    bool converged = true, update_allowed = true;
    double maximum_penetration = 0.0, maximum_constraint_violation = 0.0, penetration_tolerance = 0.0;
};
enum class BoundaryConditionType {
    dirichlet,
    pressure,
    traction,
    convection,
};

struct BoundaryConditionDefinition final {
    std::string name;
    BoundaryConditionType type;
    std::string boundary;
    Field field;
    double value;
    bool scale_with_load = false;
    std::string function{};
    double heat_transfer_coefficient = 0.0, ambient_temperature = 0.0;
    std::string coefficient_function{};
    std::string ambient_temperature_function{};
    bool use_displaced_geometry = false;
    bool configuration_explicit = false;
};

struct SpatialDefinition final {
    std::vector<RegionDefinition> regions;
    std::vector<ContactDefinition> contacts;
    std::vector<BoundaryConditionDefinition> boundary_conditions;
    std::vector<PiecewiseLinearTimeTable> time_tables{};
};

struct ContactNodeSummary final {
    double r, z;
    bool projected;
    std::size_t primary_segment;
    double gap, pressure, tributary_area, tributary_length, contact_force, tangential_traction, tangential_force,
        elastic_tangential_slip;
    bool sliding;
};

struct CartesianContactNodeSummary final {
    double x, y, z;
    bool projected;
    std::size_t primary_face;
    double gap, pressure, tributary_area, contact_force, tangential_traction, tangential_force;
    std::array<double, 3> elastic_tangential_slip;
    bool sliding;
};

struct InterfaceSummary final {
    double minimum_gap = std::numeric_limits<double>::infinity();
    double minimum_contact_gap = std::numeric_limits<double>::infinity();
    double maximum_contact_pressure = 0.0, total_heat_rate = 0.0, total_contact_force = 0.0,
           total_tangential_force = 0.0;
    std::size_t projected_contact_nodes = 0, unprojected_contact_nodes = 0, active_contact_nodes = 0;
};
enum class SpatialContributionType {
    volume,
    thermal_contact,
    mechanical_contact,
    pressure,
    traction,
    convection,
};
} // namespace fuelsim
