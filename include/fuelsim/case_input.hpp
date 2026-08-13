#pragma once
#include "fuelsim/material_functions.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "fuelsim/transient_problem.hpp"
#include <cstddef>
#include <string>
#include <vector>
namespace fuelsim {
class InputDocument;
InputDocument parse_input_file(const std::string& path);
struct InputEntry final {
    std::string key, value;
    std::size_t line;
};
class InputSection final {
  public:
    const std::string& path() const noexcept { return _path; }
    std::size_t line() const noexcept { return _line; }
    const std::vector<InputEntry>& entries() const noexcept { return _entries; }
    const InputEntry& entry(const std::string& key) const;

  private:
    friend InputDocument parse_input_file(const std::string& path);
    std::string _path;
    std::size_t _line = 0;
    std::vector<InputEntry> _entries;
};
class InputDocument final {
  public:
    const std::string& source_path() const noexcept { return _source_path; }
    const std::vector<InputSection>& sections() const noexcept { return _sections; }
    const InputSection& section(const std::string& path) const;

  private:
    friend InputDocument parse_input_file(const std::string& path);
    std::string _source_path;
    std::vector<InputSection> _sections;
};
enum class CaseProblem {
    steady,
    transient,
};
enum class CaseGeometry {
    axisymmetric_rz,
    cartesian_3d,
};
struct CaseRegionDefinition final {
    RegionDefinition spatial;
    TransientInelasticProperties transient_material;
};
struct SteadyExecutionInput final {
    std::size_t load_steps;
    double cutback_factor;
    std::size_t maximum_cutbacks;
    double minimum_load_increment;
};
struct TransientExecutionInput final {
    double end_time, initial_time_step, minimum_time_step, maximum_time_step, growth_factor, cutback_factor;
    std::size_t maximum_cutbacks;
    double load_ramp_time;
    std::string restart_file;
    std::size_t target_nonlinear_iterations, iteration_window;
    double time_error_relative_tolerance, temperature_time_absolute_tolerance, displacement_time_absolute_tolerance,
        time_error_safety_factor, strain_history_time_absolute_tolerance, stress_history_time_absolute_tolerance;
};
struct NonlinearSolverInput final {
    double absolute_tolerance, relative_tolerance, step_tolerance;
    int maximum_iterations;
    std::string linear_solver, preconditioner;
    double linear_relative_tolerance;
    int maximum_linear_iterations;
    bool backtracking_fallback, field_residual_scaling;
    double residual_reduction_tolerance, temperature_residual_absolute_tolerance,
        mechanical_residual_absolute_tolerance, temperature_residual_scale, mechanical_residual_scale;
};
struct CaseOutputInput final {
    bool console;
    std::string csv_file, exodus_file;
    std::size_t exodus_interval;
    std::string history_file;
    std::size_t history_interval, progress_interval;
    std::string checkpoint_file;
    std::size_t checkpoint_interval;
};
struct FuelSimCaseDefinition final {
    int version;
    CaseProblem problem;
    CaseGeometry geometry;
    std::string mesh_file;
    std::vector<CaseRegionDefinition> regions;
    std::vector<ContactDefinition> contacts;
    std::vector<BoundaryConditionDefinition> boundary_conditions;
    std::vector<PiecewiseLinearTimeTable> time_tables;
    SteadyExecutionInput steady_execution;
    TransientExecutionInput transient_execution;
    NonlinearSolverInput solver;
    CaseOutputInput outputs;
    SpatialDefinition spatial_definition() const;
    TransientProblemDefinition transient_definition() const;
};
FuelSimCaseDefinition read_case_input(const std::string& path);
FuelSimCaseDefinition read_case_input(const std::string& path, const MaterialFunctionRegistry& registry);
} // namespace fuelsim
