#include "fuelsim/case_input.hpp"
#include "fuelsim/checkpoint_io.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"

#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

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
            _csv << "metric,value\n"
                 << std::scientific << std::setprecision(12);
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
    return {input.absolute_tolerance, input.relative_tolerance,
            input.step_tolerance, input.maximum_iterations};
}

void write_interface_summary(const std::string& name,
                             const fuelsim::InterfaceSummary& summary,
                             CaseOutput& output) {
    const std::string prefix = "contact." + name + ".";
    output.value(prefix + "minimum_gap", summary.minimum_gap);
    output.value(prefix + "maximum_contact_pressure",
                 summary.maximum_contact_pressure);
    output.value(prefix + "total_heat_rate", summary.total_heat_rate);
    output.value(prefix + "total_contact_force", summary.total_contact_force);
    output.value(prefix + "projected_contact_nodes",
                 summary.projected_contact_nodes);
    output.value(prefix + "active_contact_nodes", summary.active_contact_nodes);
}

bool run_steady(const fuelsim::FuelSimCaseDefinition& definition,
                const fuelsim::UnstructuredQuad4Mesh& source,
                CaseOutput& output) {
    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    const fuelsim::SteadyResult result =
        fuelsim::solve_steady(problem, definition.steady_execution.load_steps,
                              solver_options(definition.solver));
    output.value("problem", "steady");
    output.value("completed", result.completed && result.solve.converged);
    output.value("convergence_reason", fuelsim::petsc_convergence_reason_name(
                                           result.solve.convergence_reason));
    output.value("regions", problem.region_count());
    output.value("contacts", problem.contact_count());
    output.value("load_steps_completed", result.completed_steps);
    output.value("nonlinear_iterations_total",
                 result.total_nonlinear_iterations);
    output.value("residual_norm", result.solve.residual_norm);
    output.value("petsc_workspace_setups",
                 result.aggregate_timing.workspace_setups);
    output.value("total_seconds", result.total_seconds);
    for (std::size_t contact = 0; contact < problem.contact_count(); ++contact)
        write_interface_summary(
            problem.contact(contact).name,
            problem.summarize_interface(contact, result.solve.state), output);
    if (result.completed && result.solve.converged &&
        !definition.outputs.exodus_file.empty())
        fuelsim::ExodusResultsIo::write_steady(definition.outputs.exodus_file,
                                               source, problem,
                                               result.solve.state);
    return result.completed && result.solve.converged;
}

class TransientOutputObserver final : public fuelsim::TransientStepObserver {
  public:
    TransientOutputObserver(fuelsim::ExodusTransientResultsWriter* results,
                            std::string checkpoint_file,
                            std::size_t checkpoint_interval)
        : _results(results), _checkpoint_file(std::move(checkpoint_file)),
          _checkpoint_interval(checkpoint_interval), _accepted_steps(0) {}

    void accepted_step(const fuelsim::TransientProblem& problem,
                       const fuelsim::TransientAcceptedStep&) override {
        ++_accepted_steps;
        if (_results != nullptr)
            _results->append(problem);
        if (!_checkpoint_file.empty() &&
            _accepted_steps % _checkpoint_interval == 0)
            fuelsim::TransientCheckpointIo::write(_checkpoint_file, problem);
    }

  private:
    fuelsim::ExodusTransientResultsWriter* _results;
    std::string _checkpoint_file;
    std::size_t _checkpoint_interval;
    std::size_t _accepted_steps;
};

bool run_transient(const fuelsim::FuelSimCaseDefinition& definition,
                   const fuelsim::UnstructuredQuad4Mesh& source,
                   CaseOutput& output) {
    fuelsim::TransientProblem problem(definition.transient_definition(),
                                      source);
    if (!definition.transient_execution.restart_file.empty())
        fuelsim::TransientCheckpointIo::restore(
            definition.transient_execution.restart_file, problem);
    std::unique_ptr<fuelsim::ExodusTransientResultsWriter> results;
    if (!definition.outputs.exodus_file.empty()) {
        results = std::make_unique<fuelsim::ExodusTransientResultsWriter>(
            definition.outputs.exodus_file, source, problem);
        results->append(problem);
    }
    TransientOutputObserver observer(results.get(),
                                     definition.outputs.checkpoint_file,
                                     definition.outputs.checkpoint_interval);
    const fuelsim::TransientTimeOptions time_options = {
        definition.transient_execution.end_time,
        definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step,
        definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor,
        definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks,
        definition.transient_execution.load_ramp_time};
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options, solver_options(definition.solver), &observer);
    if (!definition.outputs.checkpoint_file.empty())
        fuelsim::TransientCheckpointIo::write(
            definition.outputs.checkpoint_file, problem);
    output.value("problem", "transient");
    output.value("completed", result.completed);
    output.value("regions", problem.region_count());
    output.value("contacts", definition.contacts.size());
    output.value("committed_time", result.committed_time);
    output.value("accepted_steps", result.accepted_steps.size());
    output.value("total_cutbacks", result.total_cutbacks);
    output.value("nonlinear_iterations_total",
                 result.total_nonlinear_iterations);
    output.value("petsc_workspace_setups",
                 result.aggregate_timing.workspace_setups);
    output.value("total_seconds", result.total_seconds);
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const fuelsim::RegionInelasticSummary summary =
            problem.summarize_region_history(region);
        const std::string prefix =
            "region." + problem.region(region).name + ".";
        output.value(prefix + "maximum_equivalent_plastic_strain",
                     summary.maximum_equivalent_plastic_strain);
        output.value(prefix + "maximum_equivalent_creep_strain",
                     summary.maximum_equivalent_creep_strain);
    }
    for (std::size_t contact = 0; contact < definition.contacts.size();
         ++contact)
        write_interface_summary(
            definition.contacts[contact].name,
            problem.summarize_interface(contact, result.committed_state),
            output);
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
            "fuelsim input-driven axisymmetric multi-region solver\n");
        const fuelsim::UnstructuredQuad4Mesh source =
            fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
        CaseOutput output(definition.outputs);
        output.value("input_file", input_path);
        output.value("mesh_file", definition.mesh_file);
        const bool completed =
            definition.problem == fuelsim::CaseProblem::steady
                ? run_steady(definition, source, output)
                : run_transient(definition, source, output);
        return completed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim failed: " << error.what() << '\n';
        return 1;
    }
}
