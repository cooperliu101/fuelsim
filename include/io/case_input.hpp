#pragma once
#include "core/modal_beam.hpp"
#include "core/spatial_definition.hpp"
#include "material_functions.hpp"
#include "solver/solve_workflows.hpp"
#include <cstddef>
#include <string>

namespace fuelsim {
enum class CaseProblem {
    steady,
    transient,
};
enum class CaseGeometry {
    axisymmetric_rz,
    cartesian_3d,
    axisymmetric_1d,
};

struct CaseOutputInput final {
    bool console = true;
    std::string csv_file, exodus_file;
    std::size_t exodus_interval = 1;
    std::string history_file;
    std::size_t history_interval = 1, progress_interval = 1;
    std::string checkpoint_file;
    std::size_t checkpoint_interval = 1;
};

struct FuelSimCaseDefinition final {
    int version;
    std::size_t section_modes = 0;
    std::vector<ModalEndRegion> section_end_regions;
    CaseProblem problem;
    CaseGeometry geometry;
    std::string mesh_file;
    SpatialDefinition spatial;
    SteadyLoadOptions steady_execution;
    TransientTimeOptions transient_execution;
    std::string restart_file;
    SolverOptions solver;
    CaseOutputInput outputs;
};

FuelSimCaseDefinition read_case_input(const std::string& path);
FuelSimCaseDefinition read_case_input(const std::string& path, const MaterialFunctionRegistry& registry);
} // namespace fuelsim
