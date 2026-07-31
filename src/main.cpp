#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/steady_fuel_cladding_solver.hpp"
#include "fuelsim/transient_fuel_cladding_solver.hpp"

#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

class CaseOutput final {
  public:
    explicit CaseOutput(const fuelsim::CaseOutputInput& options)
        : _console(options.console) {
        if (!options.csv_file.empty()) {
            _csv.open(options.csv_file, std::ios::out | std::ios::trunc);
            if (!_csv)
                throw std::runtime_error("Could not open CSV output file '" +
                                         options.csv_file + "'");
            _csv << "metric,value\n";
            _csv << std::scientific << std::setprecision(12);
        }
        if (_console)
            std::cout << std::boolalpha << std::scientific
                      << std::setprecision(12);
    }

    void value(const std::string& key, const std::string& data) {
        if (_console)
            std::cout << key << '=' << data << '\n';
        if (_csv)
            _csv << key << ',' << data << '\n';
    }

    void value(const std::string& key, const char* data) {
        value(key, std::string(data));
    }

    void value(const std::string& key, double data) {
        if (_console)
            std::cout << key << '=' << data << '\n';
        if (_csv)
            _csv << key << ',' << data << '\n';
    }

    void value(const std::string& key, std::size_t data) {
        if (_console)
            std::cout << key << '=' << data << '\n';
        if (_csv)
            _csv << key << ',' << data << '\n';
    }

    void value(const std::string& key, int data) {
        if (_console)
            std::cout << key << '=' << data << '\n';
        if (_csv)
            _csv << key << ',' << data << '\n';
    }

    void value(const std::string& key, bool data) {
        if (_console)
            std::cout << key << '=' << std::boolalpha << data << '\n';
        if (_csv)
            _csv << key << ',' << std::boolalpha << data << '\n';
    }

  private:
    bool _console;
    std::ofstream _csv;
};

std::string extract_input_path(int& argc, char** argv) {
    std::string input_path;
    int output = 1;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string value = argv[argument];
        if (value != "-i") {
            argv[output++] = argv[argument];
            continue;
        }
        if (!input_path.empty())
            throw std::invalid_argument("fuelsim accepts exactly one -i file");
        if (argument + 1 >= argc)
            throw std::invalid_argument("fuelsim -i requires an input file");
        input_path = argv[++argument];
    }
    if (input_path.empty())
        throw std::invalid_argument(
            "Usage: fuelsim -i <case.fsi> [PETSc options]");
    argc = output;
    argv[argc] = nullptr;
    return input_path;
}

fuelsim::SolverOptions
solver_options(const fuelsim::NonlinearSolverInput& input) {
    return {
        input.absolute_tolerance,
        input.relative_tolerance,
        input.step_tolerance,
        input.maximum_iterations,
    };
}

struct ImportedFuelCladdingMeshes final {
    fuelsim::StructuredRzMesh fuel;
    fuelsim::StructuredRzMesh cladding;
};

ImportedFuelCladdingMeshes
read_meshes(const fuelsim::FuelSimCaseDefinition& definition) {
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh.file);
    return {
        fuelsim::StructuredRzMesh::from_unstructured_block(
            source, definition.mesh.fuel.block,
            definition.mesh.fuel.boundaries),
        fuelsim::StructuredRzMesh::from_unstructured_block(
            source, definition.mesh.cladding.block,
            definition.mesh.cladding.boundaries),
    };
}

bool run_steady(const fuelsim::FuelSimCaseDefinition& definition,
                const ImportedFuelCladdingMeshes& meshes, CaseOutput& output) {
    const fuelsim::SteadyFuelCladdingParameters parameters =
        definition.steady_parameters(meshes.fuel, meshes.cladding);
    const fuelsim::SteadyFuelCladdingLoadStepper load_stepper;
    const fuelsim::SteadyFuelCladdingLoadResult result =
        load_stepper.solve(parameters, meshes.fuel, meshes.cladding,
                           definition.steady_execution.load_steps,
                           solver_options(definition.solver));
    const fuelsim::SteadyFuelCladdingProblem problem(parameters, meshes.fuel,
                                                     meshes.cladding);
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.solve.state);

    output.value("problem", "steady_fuel_cladding");
    output.value("completed", result.completed && result.solve.converged);
    output.value("convergence_reason", fuelsim::petsc_convergence_reason_name(
                                           result.solve.convergence_reason));
    output.value("load_steps_completed", result.completed_steps);
    output.value("nonlinear_iterations_total",
                 result.total_nonlinear_iterations);
    output.value("residual_norm", result.solve.residual_norm);
    output.value("minimum_gap", interface.minimum_gap);
    output.value("maximum_contact_pressure",
                 interface.maximum_contact_pressure);
    output.value("total_contact_force", interface.total_contact_force);
    output.value("projected_contact_nodes", interface.projected_contact_nodes);
    output.value("active_contact_nodes", interface.active_contact_nodes);
    output.value("petsc_workspace_setups",
                 result.aggregate_timing.workspace_setups);
    output.value("total_seconds", result.total_seconds);
    return result.completed && result.solve.converged;
}

bool run_transient(const fuelsim::FuelSimCaseDefinition& definition,
                   const ImportedFuelCladdingMeshes& meshes,
                   CaseOutput& output) {
    fuelsim::TransientFuelCladdingProblem problem(
        definition.transient_parameters(meshes.fuel, meshes.cladding),
        meshes.fuel, meshes.cladding);
    const fuelsim::TransientTimeOptions time_options = {
        definition.transient_execution.end_time,
        definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step,
        definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor,
        definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks,
        definition.physics.heat_source_ramp_time,
    };
    const fuelsim::TransientFuelCladdingTimeStepper time_stepper;
    const fuelsim::TransientFuelCladdingResult result = time_stepper.solve(
        problem, time_options, solver_options(definition.solver));
    const fuelsim::RegionInelasticSummary fuel_history =
        problem.summarize_fuel_history();
    const fuelsim::RegionInelasticSummary cladding_history =
        problem.summarize_cladding_history();
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.committed_state);

    output.value("problem", "transient_fuel_cladding");
    output.value("completed", result.completed);
    output.value("committed_time", result.committed_time);
    output.value("accepted_steps", result.accepted_steps.size());
    output.value("total_cutbacks", result.total_cutbacks);
    output.value("nonlinear_iterations_total",
                 result.total_nonlinear_iterations);
    output.value("fuel_maximum_equivalent_plastic_strain",
                 fuel_history.maximum_equivalent_plastic_strain);
    output.value("fuel_maximum_equivalent_creep_strain",
                 fuel_history.maximum_equivalent_creep_strain);
    output.value("cladding_maximum_equivalent_plastic_strain",
                 cladding_history.maximum_equivalent_plastic_strain);
    output.value("cladding_maximum_equivalent_creep_strain",
                 cladding_history.maximum_equivalent_creep_strain);
    output.value("minimum_gap", interface.minimum_gap);
    output.value("maximum_contact_pressure",
                 interface.maximum_contact_pressure);
    output.value("total_contact_force", interface.total_contact_force);
    output.value("projected_contact_nodes", interface.projected_contact_nodes);
    output.value("active_contact_nodes", interface.active_contact_nodes);
    output.value("petsc_workspace_setups",
                 result.aggregate_timing.workspace_setups);
    output.value("total_seconds", result.total_seconds);
    return result.completed;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const std::string input_path = extract_input_path(argc, argv);
        const fuelsim::FuelSimCaseDefinition definition =
            fuelsim::CaseInputReader::read(input_path);
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim input-driven axisymmetric fuel-cladding solver\n");
        const ImportedFuelCladdingMeshes meshes = read_meshes(definition);
        CaseOutput output(definition.outputs);
        output.value("input_file", input_path);
        output.value("mesh_file", definition.mesh.file);

        const bool completed =
            definition.problem == fuelsim::CaseProblem::steady_fuel_cladding
                ? run_steady(definition, meshes, output)
                : run_transient(definition, meshes, output);
        return completed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim failed: " << error.what() << '\n';
        return 1;
    }
}
