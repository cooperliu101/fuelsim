#ifndef FUELSIM_CASE_INPUT_HPP
#define FUELSIM_CASE_INPUT_HPP

#include "fuelsim/input_file.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim {

enum class CaseProblem {
    steady,
    transient,
};

struct CaseRegionDefinition final {
    RegionDefinition spatial;
    TransientInelasticProperties transient_material;
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
    double load_ramp_time;
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
    std::string mesh_file;
    std::vector<CaseRegionDefinition> regions;
    std::vector<ContactDefinition> contacts;
    std::vector<BoundaryConditionDefinition> boundary_conditions;
    SteadyExecutionInput steady_execution;
    TransientExecutionInput transient_execution;
    NonlinearSolverInput solver;
    CaseOutputInput outputs;

    SteadyProblemDefinition steady_definition() const;
    TransientProblemDefinition transient_definition() const;
};

class CaseInputReader final {
  public:
    static FuelSimCaseDefinition read(const std::string& path);
};

} // namespace fuelsim

#endif
