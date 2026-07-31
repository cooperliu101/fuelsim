#ifndef FUELSIM_CASE_INPUT_HPP
#define FUELSIM_CASE_INPUT_HPP

#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/input_file.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/steady_fuel_cladding_problem.hpp"
#include "fuelsim/transient_fuel_cladding_problem.hpp"

#include <cstddef>
#include <string>

namespace fuelsim {

enum class CaseProblem {
    steady_fuel_cladding,
    transient_fuel_cladding,
};

struct MeshRegionInput final {
    std::string block;
    RzBoundaryNames boundaries;
};

struct FuelCladdingMeshInput final {
    std::string file;
    MeshRegionInput fuel;
    MeshRegionInput cladding;
};

struct FuelCladdingPhysicsInput final {
    double initial_temperature;
    double outer_temperature;
    double final_heat_source;
    double heat_source_ramp_time;
    double gap_conductivity;
    double minimum_gap;
    double contact_penalty;
};

struct SteadyExecutionInput final {
    std::size_t load_steps;
};

struct TransientExecutionInput final {
    double end_time;
    double initial_time_step;
    double minimum_time_step;
    double maximum_time_step;
    double growth_factor;
    double cutback_factor;
    std::size_t maximum_cutbacks;
};

struct NonlinearSolverInput final {
    double absolute_tolerance;
    double relative_tolerance;
    double step_tolerance;
    int maximum_iterations;
};

struct CaseOutputInput final {
    bool console;
    std::string csv_file;
};

struct FuelSimCaseDefinition final {
    int version;
    CaseProblem problem;
    FuelCladdingMeshInput mesh;
    ThermoelasticProperties fuel_thermoelastic;
    ThermoelasticProperties cladding_thermoelastic;
    TransientInelasticProperties fuel_transient;
    TransientInelasticProperties cladding_transient;
    FuelCladdingPhysicsInput physics;
    SteadyExecutionInput steady_execution;
    TransientExecutionInput transient_execution;
    NonlinearSolverInput solver;
    CaseOutputInput outputs;

    SteadyFuelCladdingParameters
    steady_parameters(const StructuredRzMesh& fuel_mesh,
                      const StructuredRzMesh& cladding_mesh) const;
    TransientFuelCladdingParameters
    transient_parameters(const StructuredRzMesh& fuel_mesh,
                         const StructuredRzMesh& cladding_mesh) const;
};

class CaseInputReader final {
  public:
    static FuelSimCaseDefinition read(const std::string& path);
};

} // namespace fuelsim

#endif
