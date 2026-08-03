#include "fuelsim/case_input.hpp"
#include "fuelsim/checkpoint_io.hpp"
#include "fuelsim/diagnostics.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

class CaseOutput final {
  public:
    explicit CaseOutput(const fuelsim::CaseOutputInput& options,
                        bool force_console = false, bool active = true)
        : _console(active && (options.console || force_console)) {
        if (active && !options.csv_file.empty()) {
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

struct CommandLine final {
    std::string input_path;
    bool check_jacobian = false;
};

CommandLine extract_command_line(int& argc, char** argv) {
    CommandLine result;
    int output = 1;
    for (int argument = 1; argument < argc; ++argument) {
        const std::string value = argv[argument];
        if (value == "--check-jacobian") {
            if (result.check_jacobian)
                throw std::invalid_argument(
                    "fuelsim accepts --check-jacobian only once");
            result.check_jacobian = true;
            continue;
        }
        if (value != "-i") {
            argv[output++] = argv[argument];
            continue;
        }
        if (!result.input_path.empty())
            throw std::invalid_argument("fuelsim accepts exactly one -i file");
        if (argument + 1 >= argc)
            throw std::invalid_argument("fuelsim -i requires an input file");
        result.input_path = argv[++argument];
    }
    if (result.input_path.empty())
        throw std::invalid_argument(
            "Usage: fuelsim -i <case.fsi> [--check-jacobian] [PETSc options]");
    argc = output;
    argv[argc] = nullptr;
    return result;
}

fuelsim::SolverOptions
solver_options(const fuelsim::NonlinearSolverInput& input) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = input.absolute_tolerance;
    result.relative_tolerance = input.relative_tolerance;
    result.step_tolerance = input.step_tolerance;
    result.maximum_iterations = input.maximum_iterations;
    if (input.linear_solver == "direct")
        result.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    else if (input.linear_solver == "gmres")
        result.linear_solver = fuelsim::SolverOptions::LinearSolver::gmres;
    if (input.preconditioner == "lu")
        result.preconditioner =
            fuelsim::SolverOptions::Preconditioner::lu;
    else if (input.preconditioner == "block_jacobi")
        result.preconditioner =
            fuelsim::SolverOptions::Preconditioner::block_jacobi;
    else if (input.preconditioner == "field_split")
        result.preconditioner =
            fuelsim::SolverOptions::Preconditioner::field_split;
    else if (input.preconditioner == "hypre")
        result.preconditioner =
            fuelsim::SolverOptions::Preconditioner::hypre;
    result.linear_relative_tolerance = input.linear_relative_tolerance;
    result.maximum_linear_iterations = input.maximum_linear_iterations;
    result.backtracking_fallback = input.backtracking_fallback;
    result.field_residual_scaling = input.field_residual_scaling;
    result.residual_reduction_tolerance =
        input.residual_reduction_tolerance;
    result.temperature_residual_absolute_tolerance =
        input.temperature_residual_absolute_tolerance;
    result.mechanical_residual_absolute_tolerance =
        input.mechanical_residual_absolute_tolerance;
    return result;
}

void write_solver_diagnostics(const fuelsim::SolveResult& solve,
                              CaseOutput& output) {
    output.value("nonlinear_attempts", solve.nonlinear_attempts);
    output.value("used_backtracking_fallback",
                 solve.used_backtracking_fallback);
    if (solve.used_backtracking_fallback) {
        output.value("basic_failure_category",
                     fuelsim::solve_failure_category_name(
                         solve.basic_failure_category));
        if (!solve.basic_failure_message.empty())
            output.value("basic_failure_message",
                         solve.basic_failure_message);
    }
    const std::array<const char*, 3> fields = {"temperature", "radial",
                                               "axial"};
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const std::string prefix =
            "residual." + std::string(fields[field]) + ".";
        output.value(prefix + "initial_l2",
                     solve.initial_field_residual_norms[field]);
        output.value(prefix + "reference_l2",
                     solve.field_residual_reference_norms[field]);
        output.value(prefix + "final_l2",
                     solve.final_field_residual_norms[field]);
        output.value(prefix + "scaling",
                     solve.field_residual_scalings[field]);
        output.value(prefix + "final_scaled_l2",
                     solve.final_scaled_field_residual_norms[field]);
    }
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
    output.value(prefix + "unprojected_contact_nodes",
                 summary.unprojected_contact_nodes);
    output.value(prefix + "active_contact_nodes", summary.active_contact_nodes);
}

std::vector<double> diagnostic_direction(const fuelsim::DofMap& dof_map) {
    std::vector<double> result(dof_map.dof_count(), 0.0);
    for (std::size_t node = 0; node < dof_map.node_count(); ++node) {
        const double index = static_cast<double>(node % 7);
        result[dof_map.temperature(node)] = 0.25 + 0.05 * index;
        result[dof_map.radial_displacement(node)] =
            1.0e-6 * (0.4 + 0.1 * index);
        result[dof_map.axial_displacement(node)] =
            -1.0e-6 * (0.3 + 0.07 * index);
    }
    return result;
}

bool write_jacobian_check(const fuelsim::NonlinearProblem& problem,
                          const fuelsim::DofMap& dof_map,
                          const std::vector<double>& state,
                          CaseOutput& output) {
    const fuelsim::DirectionalJacobianCheck check =
        fuelsim::check_directional_jacobian(
            problem, dof_map, state, diagnostic_direction(dof_map), 1.0e-4);
    const std::array<const char*, 3> fields = {"temperature", "radial",
                                               "axial"};
    bool passed = true;
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const std::string prefix =
            "jacobian." + std::string(fields[field]) + ".";
        const double reference =
            check.finite_difference_directional_derivative.l2[field];
        const double difference = check.difference.l2[field];
        const double relative =
            reference > 0.0
                ? difference / reference
                : (difference == 0.0 ? 0.0
                                     : std::numeric_limits<double>::infinity());
        output.value(prefix + "residual_l2", check.residual.l2[field]);
        output.value(prefix + "analytic_l2",
                     check.analytic_directional_derivative.l2[field]);
        output.value(prefix + "finite_difference_l2", reference);
        output.value(prefix + "difference_l2", difference);
        output.value(prefix + "relative_l2", relative);
        output.value(prefix + "maximum_absolute_difference",
                     check.difference.maximum_absolute[field]);
        passed =
            difference <= 1.0e-6 * (1.0 + reference) &&
            check.difference.maximum_absolute[field] <=
                1.0e-6 * (1.0 + check.finite_difference_directional_derivative
                                    .maximum_absolute[field]) &&
            passed;
    }
    output.value("jacobian.check_passed", passed);
    return passed;
}

bool run_steady(const fuelsim::FuelSimCaseDefinition& definition,
                const fuelsim::UnstructuredQuad4Mesh& source,
                CaseOutput& output, bool check_jacobian,
                const fuelsim::PetscSession& session) {
    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    if (check_jacobian) {
        problem.set_load_factor(1.0);
        return write_jacobian_check(problem, problem.dof_map(),
                                    problem.initial_state(), output);
    }
    const fuelsim::SteadyResult result =
        fuelsim::solve_steady(problem,
                              {definition.steady_execution.load_steps,
                               definition.steady_execution.cutback_factor,
                               definition.steady_execution.maximum_cutbacks,
                               definition.steady_execution.minimum_load_increment},
                              solver_options(definition.solver));
    output.value("problem", "steady");
    output.value("completed", result.completed && result.solve.converged);
    output.value("convergence_reason", fuelsim::petsc_convergence_reason_name(
                                           result.solve.convergence_reason));
    output.value("regions", problem.region_count());
    output.value("contacts", problem.contact_count());
    output.value("load_steps_completed", result.completed_steps);
    output.value("rejected_load_steps", result.rejected_steps.size());
    output.value("load_cutbacks", result.total_cutbacks);
    output.value("nonlinear_iterations_total",
                 result.total_nonlinear_iterations);
    output.value("residual_norm", result.solve.residual_norm);
    write_solver_diagnostics(result.solve, output);
    output.value("failure_category", fuelsim::solve_failure_category_name(
                                         result.solve.failure_category));
    if (!result.solve.failure_message.empty())
        output.value("failure_message", result.solve.failure_message);
    output.value("petsc_workspace_setups",
                 result.aggregate_timing.workspace_setups);
    output.value("total_seconds", result.total_seconds);
    if (result.completed && result.solve.converged) {
        for (std::size_t contact = 0; contact < problem.contact_count();
             ++contact)
            write_interface_summary(
                problem.contact(contact).name,
                problem.summarize_interface(contact, result.solve.state),
                output);
    }
    if (result.completed && result.solve.converged &&
        !definition.outputs.exodus_file.empty())
        session.collective_root_action([&]() {
            fuelsim::ExodusResultsIo::write_steady(
                definition.outputs.exodus_file, source, problem,
                result.solve.state);
        });
    return result.completed && result.solve.converged;
}

class TransientOutputObserver final : public fuelsim::TransientStepObserver {
  public:
    TransientOutputObserver(fuelsim::ExodusTransientResultsWriter* results,
                            fuelsim::EngineeringHistoryWriter* history,
                            std::string checkpoint_file,
                            std::size_t exodus_interval,
                            std::size_t history_interval,
                            std::size_t progress_interval,
                            std::size_t checkpoint_interval,
                            const fuelsim::PetscSession& session,
                            CaseOutput& output)
        : _results(results), _history(history),
          _checkpoint_file(std::move(checkpoint_file)),
          _exodus_interval(exodus_interval),
          _history_interval(history_interval),
          _progress_interval(progress_interval),
          _checkpoint_interval(checkpoint_interval), _accepted_steps(0),
          _exodus_at_latest(true), _history_at_latest(true),
          _last_time_step(0.0), _last_next_time_step(0.0),
          _last_nonlinear_iterations(0),
          _session(session), _output(output) {}

    void accepted_step(const fuelsim::TransientProblem& problem,
                       const fuelsim::TransientAcceptedStep& step) override {
        ++_accepted_steps;
        _last_time_step = step.time_step;
        _last_next_time_step = step.next_time_step;
        _last_nonlinear_iterations = step.nonlinear_iterations;
        _exodus_at_latest = false;
        _history_at_latest = false;
        _session.collective_root_action([&]() {
            if (_results != nullptr &&
                _accepted_steps % _exodus_interval == 0) {
                _results->append(problem);
                _exodus_at_latest = true;
            }
            if (_history != nullptr &&
                _accepted_steps % _history_interval == 0) {
                _history->append(problem, step.time_step, step.next_time_step,
                                 step.nonlinear_iterations);
                _history_at_latest = true;
            }
            if (_accepted_steps % _progress_interval == 0) {
                _output.value("progress.accepted_steps", _accepted_steps);
                _output.value("progress.time", step.time);
                _output.value("progress.time_step", step.time_step);
                _output.value("progress.next_time_step", step.next_time_step);
                _output.value("progress.nonlinear_iterations",
                              step.nonlinear_iterations);
                _output.value("progress.cutbacks", step.cutbacks);
                _output.value("progress.time_error_estimate",
                              step.time_error_estimate);
            }
            if (!_checkpoint_file.empty() &&
                _accepted_steps % _checkpoint_interval == 0)
                fuelsim::TransientCheckpointIo::write(
                    _checkpoint_file, problem, step.next_time_step);
        });
    }

    void finalize(const fuelsim::TransientProblem& problem,
                  double next_time_step) {
        _session.collective_root_action([&]() {
            if (_results != nullptr && !_exodus_at_latest)
                _results->append(problem);
            if (_history != nullptr && !_history_at_latest)
                _history->append(problem, _last_time_step,
                                 _last_next_time_step,
                                 _last_nonlinear_iterations);
            if (!_checkpoint_file.empty())
                fuelsim::TransientCheckpointIo::write(
                    _checkpoint_file, problem, next_time_step);
        });
        _exodus_at_latest = true;
        _history_at_latest = true;
    }

  private:
    fuelsim::ExodusTransientResultsWriter* _results;
    fuelsim::EngineeringHistoryWriter* _history;
    std::string _checkpoint_file;
    std::size_t _exodus_interval;
    std::size_t _history_interval;
    std::size_t _progress_interval;
    std::size_t _checkpoint_interval;
    std::size_t _accepted_steps;
    bool _exodus_at_latest;
    bool _history_at_latest;
    double _last_time_step;
    double _last_next_time_step;
    int _last_nonlinear_iterations;
    const fuelsim::PetscSession& _session;
    CaseOutput& _output;
};

bool run_transient(const fuelsim::FuelSimCaseDefinition& definition,
                   const fuelsim::UnstructuredQuad4Mesh& source,
                   CaseOutput& output, bool check_jacobian,
                   const fuelsim::PetscSession& session) {
    fuelsim::TransientProblem problem(definition.transient_definition(),
                                      source);
    double restart_time_step = 0.0;
    if (!definition.transient_execution.restart_file.empty())
        restart_time_step = fuelsim::TransientCheckpointIo::restore(
            definition.transient_execution.restart_file, problem);
    if (check_jacobian) {
        if (!(problem.committed_time() <
              definition.transient_execution.end_time))
            throw std::invalid_argument(
                "Jacobian check requires a remaining transient time step");
        double end_time = std::min(
            definition.transient_execution.end_time,
            problem.committed_time() +
                (restart_time_step > 0.0
                     ? restart_time_step
                     : definition.transient_execution.initial_time_step));
        for (const double event : problem.time_events()) {
            if (event > problem.committed_time() && event < end_time) {
                end_time = event;
                break;
            }
        }
        const double load_factor =
            definition.transient_execution.load_ramp_time == 0.0
                ? 1.0
                : std::min(end_time /
                               definition.transient_execution.load_ramp_time,
                           1.0);
        problem.begin_time_step({end_time, load_factor});
        std::vector<double> state = problem.committed_solution();
        for (const fuelsim::DirichletCondition& condition :
             problem.dirichlet_conditions())
            state.at(condition.dof) = condition.value;
        const bool passed =
            write_jacobian_check(problem, problem.dof_map(), state, output);
        problem.rollback_time_step();
        return passed;
    }
    std::unique_ptr<fuelsim::ExodusTransientResultsWriter> results;
    std::unique_ptr<fuelsim::EngineeringHistoryWriter> history;
    std::string results_path;
    if (!definition.outputs.exodus_file.empty())
        session.collective_root_action([&]() {
            results_path =
                definition.transient_execution.restart_file.empty()
                    ? definition.outputs.exodus_file
                    : fuelsim::next_results_segment_path(
                          definition.outputs.exodus_file);
            results =
                std::make_unique<fuelsim::ExodusTransientResultsWriter>(
                    results_path, source, problem);
            results->append(problem);
        });
    if (!results_path.empty())
        output.value("results_file", results_path);
    std::string history_path;
    if (!definition.outputs.history_file.empty())
        session.collective_root_action([&]() {
            history_path =
                definition.transient_execution.restart_file.empty()
                    ? definition.outputs.history_file
                    : fuelsim::next_results_segment_path(
                          definition.outputs.history_file);
            history = std::make_unique<fuelsim::EngineeringHistoryWriter>(
                history_path, problem);
            history->append(problem, 0.0,
                            restart_time_step > 0.0
                                ? restart_time_step
                                : definition.transient_execution
                                      .initial_time_step,
                            0);
        });
    if (!history_path.empty())
        output.value("history_file", history_path);
    TransientOutputObserver observer(results.get(), history.get(),
                                     definition.outputs.checkpoint_file,
                                     definition.outputs.exodus_interval,
                                     definition.outputs.history_interval,
                                     definition.outputs.progress_interval,
                                     definition.outputs.checkpoint_interval,
                                     session, output);
    const fuelsim::TransientTimeOptions time_options = {
        definition.transient_execution.end_time,
        restart_time_step > 0.0
            ? restart_time_step
            : definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step,
        definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor,
        definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks,
        definition.transient_execution.load_ramp_time,
        definition.transient_execution.target_nonlinear_iterations,
        definition.transient_execution.iteration_window,
        definition.transient_execution.time_error_relative_tolerance,
        definition.transient_execution.temperature_time_absolute_tolerance,
        definition.transient_execution.displacement_time_absolute_tolerance,
        definition.transient_execution.time_error_safety_factor};
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options, solver_options(definition.solver), &observer);
    observer.finalize(problem, result.next_time_step);
    output.value("problem", "transient");
    output.value("completed", result.completed);
    output.value(
        "termination_reason",
        fuelsim::transient_termination_reason_name(result.termination_reason));
    output.value("regions", problem.region_count());
    output.value("contacts", definition.contacts.size());
    output.value("committed_time", result.committed_time);
    output.value("next_time_step", result.next_time_step);
    output.value("accepted_steps", result.accepted_steps.size());
    output.value("time_error_rejections", result.time_error_rejections);
    output.value("rejected_steps", result.rejected_steps.size());
    if (!result.rejected_steps.empty()) {
        const fuelsim::TransientRejectedStep& rejected =
            result.rejected_steps.back();
        output.value("last_rejected.attempted_end_time",
                     rejected.attempted_end_time);
        output.value("last_rejected.time_step", rejected.time_step);
        output.value("last_rejected.cutback_index", rejected.cutback_index);
        output.value("last_rejected.nonlinear_iterations",
                     rejected.nonlinear_iterations);
        output.value("last_rejected.convergence_reason",
                     fuelsim::petsc_convergence_reason_name(
                         rejected.convergence_reason));
        output.value("last_rejected.residual_norm", rejected.residual_norm);
        output.value("last_rejected.time_error_estimate",
                     rejected.time_error_estimate);
        output.value("last_rejected.failure_category",
                     fuelsim::solve_failure_category_name(
                         rejected.failure_category));
        if (!rejected.failure_message.empty())
            output.value("last_rejected.failure_message",
                         rejected.failure_message);
    }
    output.value("total_cutbacks", result.total_cutbacks);
    output.value("nonlinear_iterations_total",
                 result.total_nonlinear_iterations);
    output.value("petsc_workspace_setups",
                 result.aggregate_timing.workspace_setups);
    write_solver_diagnostics(result.last_attempt, output);
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
        const CommandLine command = extract_command_line(argc, argv);
        const fuelsim::FuelSimCaseDefinition definition =
            fuelsim::CaseInputReader::read(command.input_path);
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim input-driven axisymmetric multi-region solver\n");
        const bool root_rank = session.rank() == 0;
        const fuelsim::UnstructuredQuad4Mesh source =
            fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
        std::unique_ptr<CaseOutput> output;
        if (!root_rank)
            output = std::make_unique<CaseOutput>(
                definition.outputs, command.check_jacobian, false);
        session.collective_root_action([&]() {
            output = std::make_unique<CaseOutput>(
                definition.outputs, command.check_jacobian, true);
        });
        output->value("input_file", command.input_path);
        output->value("mesh_file", definition.mesh_file);
        output->value("mpi_ranks", session.size());
        const bool completed =
            definition.problem == fuelsim::CaseProblem::steady
                ? run_steady(definition, source, *output,
                             command.check_jacobian, session)
                : run_transient(definition, source, *output,
                                command.check_jacobian, session);
        return completed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim failed: " << error.what() << '\n';
        return 1;
    }
}
