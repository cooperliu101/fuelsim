#ifndef FUELSIM_CASE_INPUT_HPP
#define FUELSIM_CASE_INPUT_HPP
#include "fuelsim/input_file.hpp"
#include "fuelsim/material_functions.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "fuelsim/transient_problem.hpp"
#include <cstddef>
#include <string>
#include <vector>
namespace fuelsim {
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
    double end_time;
    double initial_time_step;
    double minimum_time_step;
    double maximum_time_step;
    double growth_factor;
    double cutback_factor;
    std::size_t maximum_cutbacks;
    double load_ramp_time;
    std::string restart_file;
    std::size_t target_nonlinear_iterations;
    std::size_t iteration_window;
    double time_error_relative_tolerance;
    double temperature_time_absolute_tolerance;
    double displacement_time_absolute_tolerance;
    double time_error_safety_factor;
    double strain_history_time_absolute_tolerance;
    double stress_history_time_absolute_tolerance;
};
struct NonlinearSolverInput final {
    double absolute_tolerance;
    double relative_tolerance;
    double step_tolerance;
    int maximum_iterations;
    std::string linear_solver;
    std::string preconditioner;
    double linear_relative_tolerance;
    int maximum_linear_iterations;
    bool backtracking_fallback;
    bool field_residual_scaling;
    double residual_reduction_tolerance;
    double temperature_residual_absolute_tolerance;
    double mechanical_residual_absolute_tolerance;
    double temperature_residual_scale;
    double mechanical_residual_scale;
};
struct CaseOutputInput final {
    bool console;
    std::string csv_file;
    std::string exodus_file;
    std::size_t exodus_interval;
    std::string history_file;
    std::size_t history_interval;
    std::size_t progress_interval;
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
class CaseInputReader final {
  public:
    static FuelSimCaseDefinition read(const std::string& path);
    static FuelSimCaseDefinition read(const std::string& path, const MaterialFunctionRegistry& registry);
};
} // namespace fuelsim
#endif
