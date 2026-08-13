#ifndef FUELSIM_SPATIAL_DEFINITION_HPP
#define FUELSIM_SPATIAL_DEFINITION_HPP
#include "fuelsim/boundary.hpp"
#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/quad4_rz_kinematics.hpp"
#include "fuelsim/time_table.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
namespace fuelsim {
struct RegionDefinition final {
    std::string name;
    std::string block;
    ThermoelasticProperties material;
    double volumetric_heat_source;
    double initial_temperature;
    std::int64_t block_id = -1;
    std::string heat_source_function{};
    StrainFormulation strain_formulation = StrainFormulation::small;
};
enum class MechanicalContactFormulation {
    penalty,
    augmented_lagrangian,
};
struct ContactDefinition final {
    std::string name;
    std::string primary;
    std::string secondary;
    bool thermal;
    bool mechanical;
    double gap_conductivity;
    double minimum_gap;
    double penalty;
    double friction_coefficient = 0.0;
    bool automatic_penalty = false;
    double penalty_factor = 1.0;
    MechanicalContactFormulation mechanical_formulation = MechanicalContactFormulation::penalty;
    double penetration_tolerance = 1.0e-8;
    std::size_t maximum_augmented_iterations = 20;
};
struct AugmentedContactUpdate final {
    bool converged = true;
    bool update_allowed = true;
    double maximum_penetration = 0.0;
    double maximum_constraint_violation = 0.0;
    double penetration_tolerance = 0.0;
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
    double heat_transfer_coefficient = 0.0;
    double ambient_temperature = 0.0;
    std::string coefficient_function{};
    std::string ambient_temperature_function{};
    bool use_displaced_geometry = false;
};
struct SpatialDefinition final {
    std::vector<RegionDefinition> regions;
    std::vector<ContactDefinition> contacts;
    std::vector<BoundaryConditionDefinition> boundary_conditions;
    std::vector<PiecewiseLinearTimeTable> time_tables{};
};
struct ContactNodeSummary final {
    double r;
    double z;
    bool projected;
    std::size_t primary_segment;
    double gap;
    double pressure;
    double tributary_area;
    double tributary_length;
    double contact_force;
    double tangential_traction;
    double tangential_force;
    double elastic_tangential_slip;
    bool sliding;
};
struct InterfaceSummary final {
    double minimum_gap;
    double minimum_contact_gap;
    double maximum_contact_pressure;
    double total_heat_rate;
    double total_contact_force;
    double total_tangential_force;
    std::size_t projected_contact_nodes;
    std::size_t unprojected_contact_nodes;
    std::size_t active_contact_nodes;
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
#endif
