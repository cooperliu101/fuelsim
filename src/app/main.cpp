#include "core/nonlinear_problem.hpp"
#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "solver/petsc_solver.hpp"
#include "solver/solve_workflows.hpp"
#include <algorithm>
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
#include <vector>

namespace fuelsim {
namespace {
class CaseOutput final {
  public:
    explicit CaseOutput(bool console) : _console(console) {
        if (_console)
            std::cout << std::scientific << std::setprecision(12);
    }

    CaseOutput(const CaseOutputInput& options, bool force_console, bool active)
        : CaseOutput(active && (options.console || force_console)) {
        if (!active || options.csv_file.empty())
            return;
        _csv.open(options.csv_file, std::ios::out | std::ios::trunc);
        if (!_csv)
            throw std::runtime_error("Could not open CSV output file '" + options.csv_file + "'");
        _csv << "metric,value\n" << std::scientific << std::setprecision(12);
    }

    void value(const std::string& key, const std::string& data) {
        if (_console)
            std::cout << key << '=' << data << '\n';
        if (_csv)
            _csv << key << ',' << data << '\n';
    }

    void value(const std::string& key, const char* data) { value(key, std::string(data)); }

    void value(const std::string& key, double data) {
        if (_console)
            std::cout << key << '=' << data << '\n';
        if (_csv)
            _csv << key << ',' << data << '\n';
    }

    void value(const std::string& key, std::size_t data) { value(key, std::to_string(data)); }

    void value(const std::string& key, int data) { value(key, std::to_string(data)); }

    void value(const std::string& key, bool data) { value(key, data ? "true" : "false"); }

  private:
    bool _console;
    std::ofstream _csv;
};

void write_conservation_summary(const std::string& prefix,
    const TransientConservationSummary& summary,
    CaseOutput& output) {
    for (const TransientConservationField& field : transient_conservation_fields)
        output.value(prefix + field.name, summary.*field.member);
}

void write_time_error_components(const std::string& prefix,
    const TransientTimeErrorEstimate& estimate,
    CaseOutput& output) {
    for (const TransientFieldTimeError& field : estimate.nodal_fields)
        output.value(prefix + field.name, field.value);
    output.value(prefix + "elastic_strain", estimate.elastic_strain);
    output.value(prefix + "plastic_strain", estimate.plastic_strain);
    output.value(prefix + "creep_strain", estimate.creep_strain);
    output.value(prefix + "equivalent_plastic_strain", estimate.equivalent_plastic_strain);
    output.value(prefix + "equivalent_creep_strain", estimate.equivalent_creep_strain);
    output.value(prefix + "stress", estimate.stress);
    output.value(prefix + "contact_friction", estimate.contact_friction);
    output.value(prefix + "contact_normal_multiplier", estimate.contact_normal_multiplier);
}

class TransientOutputObserver final : public TransientStepObserver {
  public:
    TransientOutputObserver(ExodusTransientResultsWriter* results,
        EngineeringHistoryWriter* history,
        const CaseOutputInput& options,
        const PetscSession& session,
        CaseOutput& progress_output)
        : _results(results), _history(history), _options(options), _session(session),
          _progress_output(progress_output) {}

    void accepted_step(const TransientProblem& problem, const TransientAcceptedStep& step) override;
    void finalize(const TransientProblem& problem, double next_time_step);

  private:
    ExodusTransientResultsWriter* _results;
    EngineeringHistoryWriter* _history;
    const CaseOutputInput& _options;
    std::size_t _accepted_steps = 0;
    bool _exodus_at_latest = true, _history_at_latest = true;
    double _last_time_step = 0.0, _last_next_time_step = 0.0;
    int _last_nonlinear_iterations = 0;
    bool _checkpoint_at_latest = false;
    const PetscSession& _session;
    CaseOutput& _progress_output;
};

void TransientOutputObserver::accepted_step(const TransientProblem& problem, const TransientAcceptedStep& step) {
    ++_accepted_steps;
    _last_time_step = step.time_step;
    _last_next_time_step = step.next_time_step;
    _last_nonlinear_iterations = step.nonlinear_iterations;
    _exodus_at_latest = false;
    _history_at_latest = false;
    _checkpoint_at_latest = false;
    _session.collective_root_action([&]() {
        if (_results != nullptr && _accepted_steps % _options.exodus_interval == 0) {
            _results->append(problem);
            _exodus_at_latest = true;
        }
        if (_history != nullptr && _accepted_steps % _options.history_interval == 0) {
            _history->append(problem, step.time_step, step.next_time_step, step.nonlinear_iterations);
            _history_at_latest = true;
        }
        if (_accepted_steps % _options.progress_interval == 0) {
            _progress_output.value("progress.accepted_steps", _accepted_steps);
            _progress_output.value("progress.time", step.time);
            _progress_output.value("progress.time_step", step.time_step);
            _progress_output.value("progress.next_time_step", step.next_time_step);
            _progress_output.value("progress.nonlinear_iterations", step.nonlinear_iterations);
            _progress_output.value("progress.linear_iterations", step.linear_iterations);
            _progress_output.value("progress.cutbacks", step.cutbacks);
            _progress_output.value("progress.time_error_estimate", step.time_error_estimate);
            write_time_error_components("progress.time_error.", step.time_error_components, _progress_output);
            write_conservation_summary("progress.conservation.", step.conservation, _progress_output);
        }
        if (!_options.checkpoint_file.empty() && _accepted_steps % _options.checkpoint_interval == 0) {
            write_transient_checkpoint(_options.checkpoint_file, problem, step.next_time_step);
            _checkpoint_at_latest = true;
        }
    });
}

void TransientOutputObserver::finalize(const TransientProblem& problem, double next_time_step) {
    _session.collective_root_action([&]() {
        if (_results != nullptr && !_exodus_at_latest)
            _results->append(problem);
        if (_history != nullptr && !_history_at_latest)
            _history->append(problem, _last_time_step, _last_next_time_step, _last_nonlinear_iterations);
        if (!_options.checkpoint_file.empty() && !_checkpoint_at_latest)
            write_transient_checkpoint(_options.checkpoint_file, problem, next_time_step);
    });
    _exodus_at_latest = true;
    _history_at_latest = true;
    _checkpoint_at_latest = true;
}

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
                throw std::invalid_argument("fuelsim accepts --check-jacobian only once");
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
        throw std::invalid_argument("Usage: fuelsim -i <case.fsi> [--check-jacobian] [PETSc options]");
    argc = output;
    argv[argc] = nullptr;
    return result;
}

void write_memory_diagnostics(const std::string& prefix, const SolveTiming& timing, CaseOutput& output) {
    output.value(prefix + "initial_resident_bytes", timing.initial_resident_bytes);
    output.value(prefix + "setup_resident_bytes", timing.setup_resident_bytes);
    output.value(prefix + "solve_resident_bytes", timing.solve_resident_bytes);
    output.value(prefix + "final_resident_bytes", timing.final_resident_bytes);
    output.value(prefix + "minimum_peak_resident_bytes", timing.minimum_peak_resident_bytes);
    output.value(prefix + "maximum_peak_resident_bytes", timing.maximum_peak_resident_bytes);
    output.value(prefix + "total_peak_resident_bytes", timing.total_peak_resident_bytes);
}

void write_solver_diagnostics(const SolveResult& solve, bool augmented_contact, CaseOutput& output) {
    output.value("nonlinear_attempts", solve.nonlinear_attempts);
    output.value("linear_iterations", solve.linear_iterations);
    output.value("used_backtracking_fallback", solve.used_backtracking_fallback);
    output.value("augmented_lagrangian_iterations", solve.augmented_lagrangian_iterations);
    if (augmented_contact)
        output.value("maximum_contact_penetration", solve.maximum_contact_penetration);
    output.value("global_state_dofs", solve.global_state_dofs);
    output.value("maximum_shadow_state_dofs", solve.maximum_shadow_state_dofs);
    output.value("total_shadow_state_dofs", solve.total_shadow_state_dofs);
    output.value("total_remote_shadow_state_dofs", solve.total_remote_shadow_state_dofs);
    write_memory_diagnostics("memory.", solve.timing, output);
    if (solve.used_backtracking_fallback) {
        output.value("initial_failure_category", solve_failure_category_name(solve.initial_failure_category));
        if (!solve.initial_failure_message.empty())
            output.value("initial_failure_message", solve.initial_failure_message);
    }
    for (std::size_t field = 0; field < solve.field_names.size(); ++field) {
        const std::string prefix = "residual." + solve.field_names[field] + ".";
        output.value(prefix + "initial_l2", solve.initial_field_residual_norms[field]);
        output.value(prefix + "reference_l2", solve.field_residual_reference_norms[field]);
        output.value(prefix + "final_l2", solve.final_field_residual_norms[field]);
        output.value(prefix + "scaling", solve.field_residual_scalings[field]);
        output.value(prefix + "final_scaled_l2", solve.final_scaled_field_residual_norms[field]);
    }
}

void write_interface_summary(const std::string& name, const InterfaceSummary& summary, CaseOutput& output) {
    const std::string prefix = "contact." + name + ".";
    output.value(prefix + "minimum_gap", summary.minimum_gap);
    output.value(prefix + "maximum_contact_pressure", summary.maximum_contact_pressure);
    output.value(prefix + "total_heat_rate", summary.total_heat_rate);
    output.value(prefix + "total_contact_force", summary.total_contact_force);
    output.value(prefix + "total_tangential_force", summary.total_tangential_force);
    output.value(prefix + "projected_contact_nodes", summary.projected_contact_nodes);
    output.value(prefix + "unprojected_contact_nodes", summary.unprojected_contact_nodes);
    output.value(prefix + "active_contact_nodes", summary.active_contact_nodes);
}

std::string output_segment_path(const std::string& restart_file, const std::string& output_file) {
    return restart_file.empty() ? output_file : next_results_segment_path(output_file);
}

bool write_jacobian_check(const NonlinearProblem& problem, const std::vector<double>& state, CaseOutput& output) {
    constexpr double step = 1.0e-4;
    problem.validate_discretization();
    problem.validate_state(state);
    std::vector<double> direction(problem.dof_count()), analytic(problem.dof_count()), plus = state, minus = state;
    for (const FieldDescriptor& field : problem.field_layout())
        for (std::size_t dof = field.begin; dof < field.end; ++dof) {
            const double index = static_cast<double>((dof - field.begin) % 7);
            direction[dof] =
                field.category == FieldCategory::thermal ? 0.25 + 0.05 * index : 1.0e-6 * (0.4 + 0.1 * index);
            plus[dof] += step * direction[dof];
            minus[dof] -= step * direction[dof];
        }
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
        problem.evaluate_contribution(entry, state, workspace, true);
        const std::size_t local_count = workspace.dofs.size();
        for (std::size_t row = 0; row < local_count; ++row)
            for (std::size_t column = 0; column < local_count; ++column)
                analytic[workspace.dofs[row]] +=
                    workspace.jacobian[row * local_count + column] * direction[workspace.dofs[column]];
    }
    const auto residual = [&](const std::vector<double>& values) {
        std::vector<double> result(problem.dof_count());
        ContributionWorkspace local_workspace;
        problem.validate_state(values);
        for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
            problem.evaluate_contribution(entry, values, local_workspace, false);
            for (std::size_t local = 0; local < local_workspace.dofs.size(); ++local)
                result[local_workspace.dofs[local]] += local_workspace.residual[local];
        }
        for (const DirichletCondition& condition : problem.dirichlet_conditions())
            result[condition.dof] = values[condition.dof] - condition.value;
        return result;
    };
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        analytic[condition.dof] = direction[condition.dof];
    const std::vector<double> current_residual = residual(state), plus_residual = residual(plus),
                              minus_residual = residual(minus);
    bool passed = true;
    for (const FieldDescriptor& field : problem.field_layout()) {
        double residual_l2 = 0.0, analytic_l2 = 0.0, reference = 0.0, difference = 0.0, reference_maximum = 0.0,
               difference_maximum = 0.0;
        for (std::size_t dof = field.begin; dof < field.end; ++dof) {
            const double finite_difference = (plus_residual[dof] - minus_residual[dof]) / (2.0 * step),
                         error = analytic[dof] - finite_difference;
            residual_l2 = std::hypot(residual_l2, current_residual[dof]);
            analytic_l2 = std::hypot(analytic_l2, analytic[dof]);
            reference = std::hypot(reference, finite_difference);
            difference = std::hypot(difference, error);
            reference_maximum = std::max(reference_maximum, std::abs(finite_difference));
            difference_maximum = std::max(difference_maximum, std::abs(error));
        }
        const std::string prefix = "jacobian." + field.name + ".";
        const double relative = reference > 0.0 ? difference / reference
                                                : (difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
        output.value(prefix + "residual_l2", residual_l2);
        output.value(prefix + "analytic_l2", analytic_l2);
        output.value(prefix + "finite_difference_l2", reference);
        output.value(prefix + "difference_l2", difference);
        output.value(prefix + "relative_l2", relative);
        output.value(prefix + "maximum_absolute_difference", difference_maximum);
        passed = difference <= 1.0e-6 * (1.0 + reference) && difference_maximum <= 1.0e-6 * (1.0 + reference_maximum)
                 && passed;
    }
    output.value("jacobian.check_passed", passed);
    return passed;
}

void write_configuration_warnings(const spatial_detail::SpatialLayout& spatial, const PetscSession& session) {
    if (spatial.configuration_warnings().empty())
        return;
    session.collective_root_action([&]() {
        for (const std::string& warning : spatial.configuration_warnings())
            std::cerr << "fuelsim warning: " << warning << '\n';
    });
}

bool run_steady(const FuelSimCaseDefinition& definition,
    const UnstructuredBar2Mesh* bar2_source,
    const UnstructuredQuad4Mesh* rz_source,
    const UnstructuredHex8Mesh* hex_source,
    const UnstructuredHex20Mesh* hex20_source,
    const UnstructuredQuad8Mesh* quad8_source,
    CaseOutput& output,
    bool check_jacobian,
    const PetscSession& session) {
    std::unique_ptr<SteadyProblem> problem_storage;
    if (bar2_source)
        problem_storage = std::make_unique<SteadyProblem>(definition.spatial, *bar2_source);
    else if (hex20_source != nullptr)
        problem_storage = std::make_unique<SteadyProblem>(definition.spatial, *hex20_source);
    else if (hex_source != nullptr)
        problem_storage = std::make_unique<SteadyProblem>(definition.spatial, *hex_source);
    else if (quad8_source)
        problem_storage = std::make_unique<SteadyProblem>(definition.spatial, *quad8_source);
    else
        problem_storage = std::make_unique<SteadyProblem>(definition.spatial, *rz_source);
    SteadyProblem& problem = *problem_storage;
    if (bar2_source)
        write_configuration_warnings(BackendAccess::radial_spatial(problem), session);
    else if (hex_source != nullptr || hex20_source != nullptr)
        write_configuration_warnings(BackendAccess::cartesian_spatial(problem), session);
    else if (quad8_source)
        write_configuration_warnings(BackendAccess::quad8_spatial(problem), session);
    else
        write_configuration_warnings(BackendAccess::steady(problem).spatial, session);
    if (check_jacobian) {
        problem.set_load_factor(1.0);
        return write_jacobian_check(problem, problem.initial_state(), output);
    }
    const SteadyResult result = solve_steady(problem, definition.steady_execution, definition.solver);
    output.value("problem", "steady");
    output.value("completed", result.completed && result.solve.converged);
    output.value("convergence_reason", petsc_convergence_reason_name(result.solve.convergence_reason));
    output.value("regions", definition.spatial.regions.size());
    output.value("contacts", definition.spatial.contacts.size());
    output.value("load_steps_completed", result.completed_steps);
    output.value("rejected_load_steps", result.rejected_steps.size());
    output.value("load_cutbacks", result.total_cutbacks);
    output.value("nonlinear_iterations_total", result.total_nonlinear_iterations);
    output.value("linear_iterations_total", result.total_linear_iterations);
    output.value("residual_evaluations_total", result.aggregate_timing.residual_evaluations);
    output.value("jacobian_evaluations_total", result.aggregate_timing.jacobian_evaluations);
    output.value("used_small_strain_predictor", result.used_small_strain_predictor);
    output.value("load_predictor_attempts", result.load_predictor_attempts);
    output.value("load_predictor_fallbacks", result.load_predictor_fallbacks);
    if (result.used_small_strain_predictor) {
        output.value("predictor_nonlinear_iterations", result.predictor_nonlinear_iterations);
        output.value("predictor_linear_iterations", result.predictor_linear_iterations);
        output.value("predictor_jacobian_evaluations", result.predictor_timing.jacobian_evaluations);
    }
    output.value("residual_norm", result.solve.residual_norm);
    write_solver_diagnostics(result.solve, problem.uses_augmented_contact(), output);
    write_memory_diagnostics("aggregate_memory.", result.aggregate_timing, output);
    output.value("failure_category", solve_failure_category_name(result.solve.failure_category));
    if (!result.solve.failure_message.empty())
        output.value("failure_message", result.solve.failure_message);
    output.value("petsc_workspace_setups", result.aggregate_timing.workspace_setups);
    output.value("total_seconds", result.total_seconds);
    if (result.completed && result.solve.converged) {
        for (std::size_t contact = 0; contact < definition.spatial.contacts.size(); ++contact) {
            InterfaceSummary summary;
            if (bar2_source)
                summary = BackendAccess::radial_spatial(problem).summarize_interface(contact, result.solve.state);
            else if (hex_source != nullptr || hex20_source != nullptr)
                summary = BackendAccess::cartesian_spatial(problem).summarize_interface(contact, result.solve.state);
            else if (quad8_source)
                summary = BackendAccess::quad8_spatial(problem).summarize_interface(contact, result.solve.state);
            else
                summary = BackendAccess::steady(problem).spatial.summarize_interface(contact, result.solve.state);
            write_interface_summary(definition.spatial.contacts[contact].name, summary, output);
        }
    }
    if (result.completed && result.solve.converged && !definition.outputs.exodus_file.empty())
        session.collective_root_action([&]() {
            if (bar2_source)
                write_steady_results(definition.outputs.exodus_file, *bar2_source, problem, result.solve.state);
            else if (hex20_source != nullptr)
                write_steady_results(definition.outputs.exodus_file, *hex20_source, problem, result.solve.state);
            else if (hex_source != nullptr)
                write_steady_results(definition.outputs.exodus_file, *hex_source, problem, result.solve.state);
            else if (quad8_source)
                write_steady_results(definition.outputs.exodus_file, *quad8_source, problem, result.solve.state);
            else
                write_steady_results(definition.outputs.exodus_file, *rz_source, problem, result.solve.state);
        });
    return result.completed && result.solve.converged;
}

bool run_transient(const FuelSimCaseDefinition& definition,
    const UnstructuredBar2Mesh* bar2_source,
    const UnstructuredQuad4Mesh* rz_source,
    const UnstructuredHex8Mesh* hex_source,
    const UnstructuredHex20Mesh* hex20_source,
    const UnstructuredQuad8Mesh* quad8_source,
    CaseOutput& output,
    bool check_jacobian,
    const PetscSession& session) {
    std::unique_ptr<TransientProblem> problem_storage;
    if (bar2_source)
        problem_storage = std::make_unique<TransientProblem>(definition.spatial, *bar2_source);
    else if (hex20_source != nullptr)
        problem_storage = std::make_unique<TransientProblem>(definition.spatial, *hex20_source);
    else if (hex_source != nullptr)
        problem_storage = std::make_unique<TransientProblem>(definition.spatial, *hex_source);
    else if (quad8_source)
        problem_storage = std::make_unique<TransientProblem>(definition.spatial, *quad8_source);
    else
        problem_storage = std::make_unique<TransientProblem>(definition.spatial, *rz_source);
    TransientProblem& problem = *problem_storage;
    if (bar2_source)
        write_configuration_warnings(BackendAccess::radial_spatial(problem), session);
    else if (hex_source != nullptr || hex20_source != nullptr)
        write_configuration_warnings(BackendAccess::cartesian_spatial(problem), session);
    else if (quad8_source)
        write_configuration_warnings(BackendAccess::quad8_spatial(problem), session);
    else
        write_configuration_warnings(BackendAccess::transient(problem).spatial, session);
    double restart_time_step = 0.0;
    if (!definition.restart_file.empty())
        restart_time_step = restore_transient_checkpoint(definition.restart_file, problem);
    const double first_time_step =
        restart_time_step > 0.0 ? restart_time_step : definition.transient_execution.initial_time_step;
    if (check_jacobian) {
        if (!(problem.committed_time() < definition.transient_execution.end_time))
            throw std::invalid_argument("Jacobian check requires a remaining transient time step");
        double end_time = std::min(definition.transient_execution.end_time, problem.committed_time() + first_time_step);
        for (const double event : problem.time_events()) {
            if (event > problem.committed_time() && event < end_time) {
                end_time = event;
                break;
            }
        }
        const double load_factor = definition.transient_execution.load_ramp_time == 0.0
                                       ? 1.0
                                       : std::min(end_time / definition.transient_execution.load_ramp_time, 1.0);
        problem.begin_time_step({end_time, load_factor, definition.transient_execution.include_thermal_time_term});
        std::vector<double> state = problem.committed_solution();
        for (const DirichletCondition& condition : problem.dirichlet_conditions())
            state.at(condition.dof) = condition.value;
        const bool passed = write_jacobian_check(problem, state, output);
        problem.rollback_time_step();
        return passed;
    }
    std::unique_ptr<ExodusTransientResultsWriter> results;
    std::unique_ptr<EngineeringHistoryWriter> history;
    std::string results_path;
    if (!definition.outputs.exodus_file.empty())
        session.collective_root_action([&]() {
            results_path = output_segment_path(definition.restart_file, definition.outputs.exodus_file);
            if (bar2_source)
                results = std::make_unique<ExodusTransientResultsWriter>(results_path, *bar2_source, problem);
            else if (hex20_source != nullptr)
                results = std::make_unique<ExodusTransientResultsWriter>(results_path, *hex20_source, problem);
            else if (hex_source != nullptr)
                results = std::make_unique<ExodusTransientResultsWriter>(results_path, *hex_source, problem);
            else if (quad8_source)
                results = std::make_unique<ExodusTransientResultsWriter>(results_path, *quad8_source, problem);
            else
                results = std::make_unique<ExodusTransientResultsWriter>(results_path, *rz_source, problem);
            results->append(problem);
        });
    if (!results_path.empty())
        output.value("results_file", results_path);
    std::string history_path;
    if (!definition.outputs.history_file.empty())
        session.collective_root_action([&]() {
            history_path = output_segment_path(definition.restart_file, definition.outputs.history_file);
            history = std::make_unique<EngineeringHistoryWriter>(history_path, problem);
            history->append(problem, 0.0, first_time_step, 0);
        });
    if (!history_path.empty())
        output.value("history_file", history_path);
    CaseOutput progress_output(definition.outputs.console && session.rank() == 0);
    TransientOutputObserver observer(results.get(), history.get(), definition.outputs, session, progress_output);
    TransientTimeOptions time_options = definition.transient_execution;
    time_options.initial_time_step = first_time_step;
    const TransientResult result = solve_transient(problem, time_options, definition.solver, &observer);
    observer.finalize(problem, result.next_time_step);
    output.value("problem", "transient");
    output.value("completed", result.completed);
    output.value("termination_reason", transient_termination_reason_name(result.termination_reason));
    output.value("regions", definition.spatial.regions.size());
    output.value("contacts", definition.spatial.contacts.size());
    output.value("committed_time", result.committed_time);
    output.value("next_time_step", result.next_time_step);
    output.value("accepted_steps", result.accepted_steps.size());
    output.value("contributions", problem.contribution_count());
    if (!result.accepted_steps.empty()) {
        double minimum_step = result.accepted_steps.front().time_step;
        double maximum_step = minimum_step, maximum_time_error = 0.0;
        for (const auto& step : result.accepted_steps) {
            minimum_step = std::min(minimum_step, step.time_step);
            maximum_step = std::max(maximum_step, step.time_step);
            maximum_time_error = std::max(maximum_time_error, step.time_error_estimate);
        }
        output.value("minimum_accepted_time_step", minimum_step);
        output.value("maximum_accepted_time_step", maximum_step);
        output.value("maximum_accepted_time_error_estimate", maximum_time_error);
    }
    output.value("time_error_rejections", result.time_error_rejections);
    output.value("rejected_steps", result.rejected_steps.size());
    if (!result.rejected_steps.empty()) {
        const TransientRejectedStep& rejected = result.rejected_steps.back();
        output.value("last_rejected.attempted_end_time", rejected.attempted_end_time);
        output.value("last_rejected.time_step", rejected.time_step);
        output.value("last_rejected.cutback_index", rejected.cutback_index);
        output.value("last_rejected.nonlinear_iterations", rejected.nonlinear_iterations);
        output.value("last_rejected.linear_iterations", rejected.linear_iterations);
        output.value("last_rejected.convergence_reason", petsc_convergence_reason_name(rejected.convergence_reason));
        output.value("last_rejected.residual_norm", rejected.residual_norm);
        output.value("last_rejected.time_error_estimate", rejected.time_error_estimate);
        write_time_error_components("last_rejected.time_error.", rejected.time_error_components, output);
        output.value("last_rejected.failure_category", solve_failure_category_name(rejected.failure_category));
        if (!rejected.failure_message.empty())
            output.value("last_rejected.failure_message", rejected.failure_message);
    }
    output.value("total_cutbacks", result.total_cutbacks);
    output.value("nonlinear_iterations_total", result.total_nonlinear_iterations);
    output.value("linear_iterations_total", result.total_linear_iterations);
    output.value("residual_evaluations_total", result.aggregate_timing.residual_evaluations);
    output.value("jacobian_evaluations_total", result.aggregate_timing.jacobian_evaluations);
    output.value("petsc_workspace_setups", result.aggregate_timing.workspace_setups);
    write_solver_diagnostics(result.last_attempt, problem.uses_augmented_contact(), output);
    write_memory_diagnostics("aggregate_memory.", result.aggregate_timing, output);
    output.value("total_seconds", result.total_seconds);
    write_conservation_summary("conservation.", problem.last_conservation_summary(), output);
    for (std::size_t region = 0; region < definition.spatial.regions.size(); ++region) {
        const RegionStateSummary summary = problem.summarize_region(region);
        const std::string prefix = "region." + definition.spatial.regions[region].name + ".";
        output.value(prefix + "maximum_equivalent_plastic_strain", summary.maximum_equivalent_plastic_strain);
        output.value(prefix + "maximum_equivalent_creep_strain", summary.maximum_equivalent_creep_strain);
    }
    for (std::size_t contact = 0; contact < definition.spatial.contacts.size(); ++contact) {
        InterfaceSummary summary;
        if (bar2_source)
            summary = BackendAccess::radial_spatial(problem).summarize_interface(contact, result.committed_state);
        else if (hex_source != nullptr || hex20_source != nullptr)
            summary = BackendAccess::cartesian_spatial(problem).summarize_interface(contact, result.committed_state);
        else if (quad8_source)
            summary = BackendAccess::quad8_spatial(problem).summarize_interface(contact, result.committed_state);
        else
            summary = BackendAccess::transient(problem).spatial.summarize_interface(contact, result.committed_state);
        write_interface_summary(definition.spatial.contacts[contact].name, summary, output);
    }
    return result.completed;
}

int run_application(int argc, char** argv) {
    try {
        const CommandLine command = extract_command_line(argc, argv);
        const FuelSimCaseDefinition definition = read_case_input(command.input_path);
        PetscSession session(argc, argv, "fuelsim input-driven multi-region thermo-mechanics solver\n");
        const bool root_rank = session.rank() == 0;
        std::unique_ptr<UnstructuredBar2Mesh> bar2_source;
        std::unique_ptr<UnstructuredQuad4Mesh> rz_source;
        std::unique_ptr<UnstructuredQuad8Mesh> quad8_source;
        std::unique_ptr<UnstructuredHex8Mesh> hex_source;
        std::unique_ptr<UnstructuredHex20Mesh> hex20_source;
        if (definition.geometry == CaseGeometry::axisymmetric_1d)
            bar2_source = std::make_unique<UnstructuredBar2Mesh>(read_exodus_bar2(definition.mesh_file));
        else if (definition.geometry == CaseGeometry::cartesian_3d)
            if (exodus_uses_hex20(definition.mesh_file))
                hex20_source = std::make_unique<UnstructuredHex20Mesh>(read_exodus_hex20(definition.mesh_file));
            else
                hex_source = std::make_unique<UnstructuredHex8Mesh>(read_exodus_hex8(definition.mesh_file));
        else if (exodus_uses_quad8(definition.mesh_file))
            quad8_source = std::make_unique<UnstructuredQuad8Mesh>(read_exodus_quad8(definition.mesh_file));
        else
            rz_source = std::make_unique<UnstructuredQuad4Mesh>(read_exodus_quad4(definition.mesh_file));
        std::unique_ptr<CaseOutput> output;
        if (!root_rank)
            output = std::make_unique<CaseOutput>(definition.outputs, command.check_jacobian, false);
        session.collective_root_action(
            [&]() { output = std::make_unique<CaseOutput>(definition.outputs, command.check_jacobian, true); });
        output->value("input_file", command.input_path);
        output->value("mesh_file", definition.mesh_file);
        output->value("mpi_ranks", session.size());
        const bool completed = definition.problem == CaseProblem::steady ? run_steady(definition,
                                                                               bar2_source.get(),
                                                                               rz_source.get(),
                                                                               hex_source.get(),
                                                                               hex20_source.get(),
                                                                               quad8_source.get(),
                                                                               *output,
                                                                               command.check_jacobian,
                                                                               session)
                                                                         : run_transient(definition,
                                                                               bar2_source.get(),
                                                                               rz_source.get(),
                                                                               hex_source.get(),
                                                                               hex20_source.get(),
                                                                               quad8_source.get(),
                                                                               *output,
                                                                               command.check_jacobian,
                                                                               session);
        return completed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim failed: " << error.what() << '\n';
        return 1;
    }
}
} // namespace
} // namespace fuelsim

int main(int argc, char** argv) {
    return fuelsim::run_application(argc, argv);
}
