#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void write_reference(const std::string& path, const std::vector<double>& state) {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output)
        throw std::runtime_error("Could not write MPI reference: " + path);
    output << state.size() << '\n' << std::setprecision(17);
    for (double value : state)
        output << value << '\n';
    if (!output)
        throw std::runtime_error("Could not complete MPI reference: " + path);
}

std::vector<double> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MPI reference: " + path);
    std::size_t count = 0;
    input >> count;
    std::vector<double> result(count, 0.0);
    for (double& value : result)
        input >> value;
    if (!input)
        throw std::runtime_error("MPI reference is incomplete: " + path);
    return result;
}

std::vector<double> flatten_transient_state(const fuelsim::TransientProblem& problem) {
    std::vector<double> result = {problem.committed_time(), problem.committed_load_factor()};
    result.insert(result.end(), problem.committed_solution().begin(), problem.committed_solution().end());
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        for (std::size_t element = 0; element < problem.region_mesh(region).elements().size(); ++element) {
            const fuelsim::Quad4MaterialHistory& history = problem.material_history(region, element);
            const auto& stresses = problem.material_stress(region, element);
            for (std::size_t q = 0; q < history.size(); ++q) {
                result.insert(result.end(), history[q].elastic_strain.begin(), history[q].elastic_strain.end());
                result.insert(result.end(), history[q].plastic_strain.begin(), history[q].plastic_strain.end());
                result.insert(result.end(), history[q].creep_strain.begin(), history[q].creep_strain.end());
                result.push_back(history[q].equivalent_plastic_strain);
                result.push_back(history[q].equivalent_creep_strain);
                result.push_back(stresses[q].rr);
                result.push_back(stresses[q].zz);
                result.push_back(stresses[q].hoop);
                result.push_back(stresses[q].rz);
            }
        }
    }
    for (const auto& contact : problem.committed_state().contact_histories) {
        for (const fuelsim::ContactPointHistory& history : contact) {
            result.push_back(history.elastic_tangential_slip);
            result.push_back(history.normal_multiplier);
            result.push_back(history.sliding ? 1.0 : 0.0);
        }
    }
    return result;
}

void compare_reference(const std::string& path, const std::vector<double>& state, double tolerance) {
    const std::vector<double> reference = read_reference(path);
    if (reference.size() != state.size())
        throw std::runtime_error("MPI reference state size differs");
    double maximum_absolute = 0.0;
    double maximum_scaled = 0.0;
    std::size_t maximum_scaled_index = 0;
    for (std::size_t value = 0; value < reference.size(); ++value) {
        const double difference = std::abs(state[value] - reference[value]);
        maximum_absolute = std::max(maximum_absolute, difference);
        const double scaled = difference / (1.0 + std::abs(reference[value]));
        if (scaled > maximum_scaled) {
            maximum_scaled = scaled;
            maximum_scaled_index = value;
        }
    }
    if (!(maximum_scaled < tolerance)) {
        std::ostringstream message;
        message << std::scientific << std::setprecision(12) << "one/multi-rank state difference exceeds tolerance: "
                << "maximum absolute=" << maximum_absolute << ", maximum scaled=" << maximum_scaled
                << ", index=" << maximum_scaled_index << ", one-rank=" << reference[maximum_scaled_index]
                << ", multi-rank=" << state[maximum_scaled_index];
        throw std::runtime_error(message.str());
    }
    std::cout << std::scientific << std::setprecision(12) << "mpi_equivalence_maximum_absolute=" << maximum_absolute
              << '\n'
              << "mpi_equivalence_maximum_scaled=" << maximum_scaled << '\n';
}

struct TransientStateLayout final {
    std::size_t _node_count = 0;
    std::size_t _quadrature_point_count = 0;
    std::size_t _contact_point_count = 0;
    std::size_t _flattened_size = 0;
};

TransientStateLayout transient_state_layout(const fuelsim::TransientProblem& problem) {
    if (problem.dof_count() % 3 != 0)
        throw std::runtime_error("Transient MPI state does not contain three complete nodal fields");

    const fuelsim::TransientCommittedState committed = problem.committed_state();
    if (committed.solution.size() != problem.dof_count())
        throw std::runtime_error("Transient MPI committed solution size differs from the problem degree-of-freedom "
                                 "count");
    if (committed.material_histories.size() != problem.region_count() ||
        committed.material_stresses.size() != problem.region_count())
        throw std::runtime_error("Transient MPI material-state region layout differs from the problem");
    if (committed.contact_histories.size() != problem.definition().spatial.contacts.size())
        throw std::runtime_error("Transient MPI contact-history pair count differs from the problem");

    TransientStateLayout layout;
    layout._node_count = problem.dof_count() / 3;
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const std::size_t element_count = problem.region_mesh(region).elements().size();
        if (committed.material_histories[region].size() != element_count ||
            committed.material_stresses[region].size() != element_count)
            throw std::runtime_error("Transient MPI material-state element layout differs from the problem");
        for (std::size_t element = 0; element < element_count; ++element) {
            const fuelsim::Quad4MaterialHistory& history = committed.material_histories[region][element];
            const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses =
                committed.material_stresses[region][element];
            if (history.size() != stresses.size())
                throw std::runtime_error("Transient MPI material history and stress quadrature layouts differ");
            layout._quadrature_point_count += history.size();
        }
    }
    for (std::size_t contact = 0; contact < committed.contact_histories.size(); ++contact) {
        if (committed.contact_histories[contact].size() != problem.contact_secondary_source_nodes(contact).size())
            throw std::runtime_error("Transient MPI contact-history point count differs from the secondary boundary");
        layout._contact_point_count += committed.contact_histories[contact].size();
    }

    constexpr std::size_t values_per_quadrature_point = 18;
    constexpr std::size_t values_per_contact_point = 3;
    layout._flattened_size = 2 + problem.dof_count() + values_per_quadrature_point * layout._quadrature_point_count +
                             values_per_contact_point * layout._contact_point_count;
    return layout;
}

class DifferenceSummary final {
  public:
    DifferenceSummary(std::string name, std::string unit, double absolute_tolerance, double relative_tolerance,
                      bool exact)
        : _name(std::move(name)), _unit(std::move(unit)), _absolute_tolerance(absolute_tolerance),
          _relative_tolerance(relative_tolerance), _exact(exact) {}

    void add(double reference, double actual, std::size_t flattened_index) {
        ++_count;
        double difference = 0.0;
        double tolerance_ratio = 0.0;
        if (!std::isfinite(reference) || !std::isfinite(actual)) {
            difference = std::numeric_limits<double>::infinity();
            tolerance_ratio = std::numeric_limits<double>::infinity();
        } else {
            difference = std::abs(actual - reference);
            if (_exact)
                tolerance_ratio = difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
            else
                tolerance_ratio = difference / (_absolute_tolerance + _relative_tolerance * std::abs(reference));
        }
        if (_count == 1 || difference > _maximum_absolute_difference) {
            _maximum_absolute_difference = difference;
            _maximum_absolute_index = flattened_index;
        }
        if (_count == 1 || tolerance_ratio > _maximum_tolerance_ratio) {
            _maximum_tolerance_ratio = tolerance_ratio;
            _maximum_tolerance_index = flattened_index;
            _maximum_tolerance_reference = reference;
            _maximum_tolerance_actual = actual;
        }
    }

    bool passed() const noexcept {
        return _maximum_tolerance_ratio <= 1.0;
    }

    double maximum_tolerance_ratio() const noexcept {
        return _maximum_tolerance_ratio;
    }

    void print() const {
        std::cout << std::scientific << std::setprecision(12) << "transient_mpi_" << _name
                  << "_maximum_absolute=" << _maximum_absolute_difference << ", unit=" << _unit
                  << ", maximum_absolute_index=" << _maximum_absolute_index
                  << ", maximum_tolerance_ratio=" << _maximum_tolerance_ratio << ", count=" << _count << '\n';
    }

    std::string failure_message() const {
        std::ostringstream message;
        message << std::scientific << std::setprecision(12) << "transient one/multi-rank " << _name
                << " difference exceeds " << (_exact ? "exact equality" : "the absolute/relative tolerance")
                << ": flattened index=" << _maximum_tolerance_index << ", one-rank=" << _maximum_tolerance_reference
                << ", multi-rank=" << _maximum_tolerance_actual
                << ", absolute difference=" << std::abs(_maximum_tolerance_actual - _maximum_tolerance_reference)
                << ", maximum tolerance ratio=" << _maximum_tolerance_ratio << ", unit=" << _unit;
        if (!_exact)
            message << ", absolute tolerance=" << _absolute_tolerance << ", relative tolerance=" << _relative_tolerance;
        return message.str();
    }

  private:
    std::string _name;
    std::string _unit;
    double _absolute_tolerance;
    double _relative_tolerance;
    bool _exact;
    std::size_t _count = 0;
    double _maximum_absolute_difference = 0.0;
    std::size_t _maximum_absolute_index = 0;
    double _maximum_tolerance_ratio = 0.0;
    std::size_t _maximum_tolerance_index = 0;
    double _maximum_tolerance_reference = 0.0;
    double _maximum_tolerance_actual = 0.0;
};

void add_difference(DifferenceSummary& summary, const std::vector<double>& reference, const std::vector<double>& state,
                    std::size_t flattened_index) {
    summary.add(reference[flattened_index], state[flattened_index], flattened_index);
}

void compare_transient_reference(const std::string& path, const std::vector<double>& state,
                                 const fuelsim::TransientProblem& problem, bool rank_sensitive_adaptive_path) {
    const std::vector<double> reference = read_reference(path);
    const TransientStateLayout layout = transient_state_layout(problem);
    if (state.size() != layout._flattened_size || reference.size() != layout._flattened_size) {
        std::ostringstream message;
        message << "Transient MPI flattened-state layout differs: expected=" << layout._flattened_size
                << ", one-rank=" << reference.size() << ", multi-rank=" << state.size()
                << ", nodes=" << layout._node_count << ", quadrature points=" << layout._quadrature_point_count
                << ", contact points=" << layout._contact_point_count;
        throw std::runtime_error(message.str());
    }

    DifferenceSummary time("time", "s", 1.0e-12, 1.0e-12, false);
    DifferenceSummary load_factor("load_factor", "dimensionless", 1.0e-12, 1.0e-12, false);
    // The integrated adaptive path can differ by one accepted step when a
    // rank-dependent nonlinear iteration count crosses its iteration-based
    // step-size threshold. Its field gates remain at least one order of
    // magnitude tighter than the case's absolute time-error scales. Fixed-step
    // transient regressions retain or strengthen their prior field gates.
    const double temperature_absolute = rank_sensitive_adaptive_path ? 1.0e-3 : 1.0e-10;
    const double displacement_absolute = 1.0e-10;
    const double strain_absolute = rank_sensitive_adaptive_path ? 1.0e-8 : 1.0e-10;
    // The fixed-step case contains stress components near zero; distributed
    // summation changes them by about 1e-7 Pa while nonzero stresses remain
    // constrained by the 1e-10 relative gate.
    const double stress_absolute = rank_sensitive_adaptive_path ? 1.0e3 : 1.0e-7;
    const double contact_slip_absolute = rank_sensitive_adaptive_path ? 1.0e-11 : 1.0e-10;
    const double contact_multiplier_absolute = rank_sensitive_adaptive_path ? 1.0e-11 : 1.0e-10;
    const double field_relative = rank_sensitive_adaptive_path ? 1.0e-7 : 1.0e-10;
    const double temperature_relative = 1.0e-10;
    const double displacement_relative = rank_sensitive_adaptive_path ? 1.0e-8 : 1.0e-10;
    const double contact_relative = rank_sensitive_adaptive_path ? 1.0e-8 : 1.0e-10;
    DifferenceSummary temperature("temperature", "K", temperature_absolute, temperature_relative, false);
    DifferenceSummary radial_displacement("radial_displacement", "m", displacement_absolute, displacement_relative,
                                          false);
    DifferenceSummary axial_displacement("axial_displacement", "m", displacement_absolute, displacement_relative,
                                         false);
    DifferenceSummary elastic_strain("elastic_strain", "dimensionless", strain_absolute, field_relative, false);
    DifferenceSummary plastic_strain("plastic_strain", "dimensionless", strain_absolute, field_relative, false);
    DifferenceSummary creep_strain("creep_strain", "dimensionless", strain_absolute, field_relative, false);
    DifferenceSummary equivalent_plastic_strain("equivalent_plastic_strain", "dimensionless", strain_absolute,
                                                field_relative, false);
    DifferenceSummary equivalent_creep_strain("equivalent_creep_strain", "dimensionless", strain_absolute,
                                              field_relative, false);
    DifferenceSummary stress("stress", "Pa", stress_absolute, field_relative, false);
    DifferenceSummary contact_slip("contact_elastic_tangential_slip", "m", contact_slip_absolute, contact_relative,
                                   false);
    DifferenceSummary contact_multiplier("contact_normal_multiplier", "Pa", contact_multiplier_absolute,
                                         contact_relative, false);
    DifferenceSummary contact_sliding("contact_sliding", "boolean", 0.0, 0.0, true);

    std::size_t index = 0;
    add_difference(time, reference, state, index++);
    add_difference(load_factor, reference, state, index++);
    for (std::size_t node = 0; node < layout._node_count; ++node)
        add_difference(temperature, reference, state, index++);
    for (std::size_t node = 0; node < layout._node_count; ++node)
        add_difference(radial_displacement, reference, state, index++);
    for (std::size_t node = 0; node < layout._node_count; ++node)
        add_difference(axial_displacement, reference, state, index++);
    for (std::size_t point = 0; point < layout._quadrature_point_count; ++point) {
        for (std::size_t component = 0; component < 4; ++component)
            add_difference(elastic_strain, reference, state, index++);
        for (std::size_t component = 0; component < 4; ++component)
            add_difference(plastic_strain, reference, state, index++);
        for (std::size_t component = 0; component < 4; ++component)
            add_difference(creep_strain, reference, state, index++);
        add_difference(equivalent_plastic_strain, reference, state, index++);
        add_difference(equivalent_creep_strain, reference, state, index++);
        for (std::size_t component = 0; component < 4; ++component)
            add_difference(stress, reference, state, index++);
    }
    for (std::size_t point = 0; point < layout._contact_point_count; ++point) {
        add_difference(contact_slip, reference, state, index++);
        add_difference(contact_multiplier, reference, state, index++);
        add_difference(contact_sliding, reference, state, index++);
    }
    if (index != layout._flattened_size)
        throw std::runtime_error("Transient MPI field comparison did not consume the complete flattened state");

    const std::array<DifferenceSummary*, 14> summaries = {
        &time,           &load_factor,    &temperature,        &radial_displacement,       &axial_displacement,
        &elastic_strain, &plastic_strain, &creep_strain,       &equivalent_plastic_strain, &equivalent_creep_strain,
        &stress,         &contact_slip,   &contact_multiplier, &contact_sliding,
    };
    const DifferenceSummary* worst_failure = nullptr;
    for (const DifferenceSummary* summary : summaries) {
        summary->print();
        if (!summary->passed() &&
            (worst_failure == nullptr || summary->maximum_tolerance_ratio() > worst_failure->maximum_tolerance_ratio()))
            worst_failure = summary;
    }
    if (worst_failure != nullptr)
        throw std::runtime_error(worst_failure->failure_message());
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_mpi_equivalence_benchmark "
                     "<write|compare|compare_field_split|compare_block_jacobi|"
                     "compare_hypre|write_transient|compare_transient|"
                     "write_transient_integrated|compare_transient_integrated|"
                     "test_io_failure> "
                     "<reference> "
                     "<case.fsi> [PETSc options]\n";
        return 2;
    }
    try {
        const std::string mode = argv[1];
        const std::string reference_path = argv[2];
        const std::string input_path = argv[3];
        for (int index = 4; index < argc; ++index)
            argv[index - 3] = argv[index];
        argc -= 3;
        argv[argc] = nullptr;

        fuelsim::PetscSession session(argc, argv, "fuelsim one/multi-rank equivalence benchmark\n");
        if (mode == "test_io_failure") {
            if (session.size() != 2)
                throw std::invalid_argument("Collective I/O failure test requires two ranks");
            bool caught = false;
            try {
                session.collective_root_action([]() { throw std::runtime_error("intentional root I/O failure"); });
            } catch (const std::runtime_error& error) {
                caught = std::string(error.what()).find("collective root-rank I/O failed") != std::string::npos;
            }
            if (!caught)
                throw std::runtime_error("Collective root I/O failure did not reach every rank");
            session.collective_root_action([]() {});
            if (session.rank() == 0)
                std::cout << "[PASS] root I/O failure reached every rank\n";
            return 0;
        }
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::CaseInputReader::read(input_path);
        const fuelsim::UnstructuredQuad4Mesh source = fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
        fuelsim::SolverOptions options;
        options.absolute_tolerance = definition.solver.absolute_tolerance;
        options.relative_tolerance = definition.solver.relative_tolerance;
        options.step_tolerance = definition.solver.step_tolerance;
        options.maximum_iterations = definition.solver.maximum_iterations;
        options.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
        options.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
        options.backtracking_fallback = definition.solver.backtracking_fallback;
        options.field_residual_scaling = definition.solver.field_residual_scaling;
        options.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
        options.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
        options.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
        options.temperature_residual_scale = definition.solver.temperature_residual_scale;
        options.mechanical_residual_scale = definition.solver.mechanical_residual_scale;
        if (definition.problem == fuelsim::CaseProblem::transient) {
            const bool write_transient = mode == "write_transient" || mode == "write_transient_integrated";
            const bool compare_transient = mode == "compare_transient" || mode == "compare_transient_integrated";
            const bool integrated_path = mode == "write_transient_integrated" || mode == "compare_transient_integrated";
            if (!write_transient && !compare_transient)
                throw std::invalid_argument("Transient MPI case requires write_transient or "
                                            "compare_transient mode");
            fuelsim::TransientProblem problem(definition.transient_definition(), source);
            const fuelsim::TransientExecutionInput& execution = definition.transient_execution;
            const fuelsim::TransientResult result = fuelsim::solve_transient(
                problem,
                {execution.end_time, execution.initial_time_step, execution.minimum_time_step,
                 execution.maximum_time_step, execution.growth_factor, execution.cutback_factor,
                 execution.maximum_cutbacks, execution.load_ramp_time, execution.target_nonlinear_iterations,
                 execution.iteration_window, execution.time_error_relative_tolerance,
                 execution.temperature_time_absolute_tolerance, execution.displacement_time_absolute_tolerance,
                 execution.time_error_safety_factor, execution.strain_history_time_absolute_tolerance,
                 execution.stress_history_time_absolute_tolerance},
                options);
            if (!result.completed || result.aggregate_timing.workspace_setups != 1) {
                std::ostringstream message;
                message << "Transient MPI equivalence solve failed or rebuilt "
                           "its workspace: completed="
                        << result.completed << ", setups=" << result.aggregate_timing.workspace_setups
                        << ", category=" << fuelsim::solve_failure_category_name(result.last_attempt.failure_category)
                        << ", message=" << result.last_attempt.failure_message;
                throw std::runtime_error(message.str());
            }
            if (integrated_path) {
                const std::size_t accepted_steps = result.accepted_steps.size();
                const bool expected_steps = session.size() <= 2
                                                ? accepted_steps == 97
                                                : session.size() == 4 && (accepted_steps == 97 || accepted_steps == 98);
                if (!expected_steps || result.time_error_rejections != 0) {
                    std::ostringstream message;
                    message << "Integrated transient MPI adaptive path differs from its qualified envelope: ranks="
                            << session.size() << ", accepted steps=" << accepted_steps
                            << ", time-error rejections=" << result.time_error_rejections;
                    throw std::runtime_error(message.str());
                }
            }
            const fuelsim::SolveResult& shadow = result.last_attempt;
            if (shadow.global_state_dofs != problem.dof_count())
                throw std::runtime_error("Transient shadow-state global size is incorrect");
            if (session.size() == 1 && shadow.maximum_shadow_state_dofs != problem.dof_count())
                throw std::runtime_error("One-rank transient solve does not cover the full state");
            if (session.size() > 1) {
                const std::size_t rank = static_cast<std::size_t>(session.rank());
                const std::size_t ranks = static_cast<std::size_t>(session.size());
                const std::size_t expected_begin = problem.contribution_count() * rank / ranks;
                const std::size_t expected_end = problem.contribution_count() * (rank + 1) / ranks;
                if (shadow.local_contribution_begin != expected_begin || shadow.local_contribution_end != expected_end)
                    throw std::runtime_error("Transient MPI contribution partition differs from ownership contract");
                if (shadow.maximum_shadow_state_dofs > problem.dof_count() ||
                    shadow.total_shadow_state_dofs > ranks * problem.dof_count() ||
                    shadow.total_remote_shadow_state_dofs == 0)
                    throw std::runtime_error("Multi-rank transient shadow-state bounds failed: global=" +
                                             std::to_string(problem.dof_count()) +
                                             ", maximum=" + std::to_string(shadow.maximum_shadow_state_dofs) +
                                             ", total=" + std::to_string(shadow.total_shadow_state_dofs) +
                                             ", remote=" + std::to_string(shadow.total_remote_shadow_state_dofs));
            }
            if (session.rank() == 0) {
                double minimum_time_step = std::numeric_limits<double>::infinity();
                double maximum_time_step = 0.0;
                double maximum_time_error = 0.0;
                for (const fuelsim::TransientAcceptedStep& step : result.accepted_steps) {
                    minimum_time_step = std::min(minimum_time_step, step.time_step);
                    maximum_time_step = std::max(maximum_time_step, step.time_step);
                    maximum_time_error = std::max(maximum_time_error, step.time_error_estimate);
                }
                std::cout << "transient_accepted_steps=" << result.accepted_steps.size() << '\n'
                          << "transient_time_error_rejections=" << result.time_error_rejections << '\n'
                          << "transient_minimum_time_step=" << minimum_time_step << '\n'
                          << "transient_maximum_time_step=" << maximum_time_step << '\n'
                          << "transient_maximum_time_error=" << maximum_time_error << '\n'
                          << "transient_total_nonlinear_iterations=" << result.total_nonlinear_iterations << '\n';
            }
            const std::vector<double> state = flatten_transient_state(problem);
            if (write_transient) {
                if (session.size() != 1)
                    throw std::invalid_argument("Transient MPI reference requires one rank");
                write_reference(reference_path, state);
                std::cout << "[PASS] wrote transient one-rank MPI reference\n";
                return 0;
            }
            if (session.size() < 2)
                throw std::invalid_argument("Transient MPI comparison requires at least two ranks");
            if (session.rank() == 0) {
                compare_transient_reference(reference_path, state, problem, integrated_path);
                std::cout << "transient_maximum_shadow_state_dofs=" << shadow.maximum_shadow_state_dofs << '\n'
                          << "transient_total_remote_shadow_state_dofs=" << shadow.total_remote_shadow_state_dofs
                          << '\n';
                std::cout << "[PASS] transient nodal, integration-point, and contact-history multi-rank equivalence\n";
            }
            return 0;
        }

        fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
        const bool field_split = mode == "compare_field_split";
        const bool block_jacobi = mode == "compare_block_jacobi";
        const bool hypre = mode == "compare_hypre";
        options.linear_solver = field_split || block_jacobi || hypre ? fuelsim::SolverOptions::LinearSolver::gmres
                                                                     : fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner = field_split    ? fuelsim::SolverOptions::Preconditioner::field_split
                                 : block_jacobi ? fuelsim::SolverOptions::Preconditioner::block_jacobi
                                 : hypre        ? fuelsim::SolverOptions::Preconditioner::hypre
                                                : fuelsim::SolverOptions::Preconditioner::lu;
        const fuelsim::SteadyResult result = fuelsim::solve_steady(
            problem,
            {field_split ? 2U : definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
             definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
            options);
        if (!result.completed || !result.solve.converged)
            throw std::runtime_error("MPI equivalence solve did not converge");
        if (result.aggregate_timing.workspace_setups != 1)
            throw std::runtime_error("MPI equivalence solve did not reuse one workspace");
        if (result.solve.global_state_dofs != problem.dof_count())
            throw std::runtime_error("Steady shadow-state global size is incorrect");
        if (session.size() == 1 && result.solve.maximum_shadow_state_dofs != problem.dof_count())
            throw std::runtime_error("One-rank steady solve does not cover the full state");

        if (mode == "write") {
            if (session.size() != 1)
                throw std::invalid_argument("MPI reference must be written with one rank");
            write_reference(reference_path, result.solve.state);
            std::cout << "[PASS] wrote one-rank MPI reference\n";
            return 0;
        }
        if (mode != "compare" && !field_split && !block_jacobi && !hypre)
            throw std::invalid_argument("Unknown MPI equivalence mode: " + mode);
        if (session.size() != 2)
            throw std::invalid_argument("MPI comparison must run with exactly two ranks");
        const std::size_t expected_begin =
            problem.contribution_count() * static_cast<std::size_t>(result.solve.mpi_rank) / 2U;
        const std::size_t expected_end =
            problem.contribution_count() * static_cast<std::size_t>(result.solve.mpi_rank + 1) / 2U;
        if (result.solve.local_contribution_begin != expected_begin ||
            result.solve.local_contribution_end != expected_end)
            throw std::runtime_error("MPI contribution partition differs from ownership contract");
        if (!(result.solve.total_shadow_state_dofs < 2 * problem.dof_count()) ||
            result.solve.total_remote_shadow_state_dofs == 0)
            throw std::runtime_error(
                "Two-rank steady shadow-state bounds failed: global=" + std::to_string(problem.dof_count()) +
                ", maximum=" + std::to_string(result.solve.maximum_shadow_state_dofs) +
                ", total=" + std::to_string(result.solve.total_shadow_state_dofs) +
                ", remote=" + std::to_string(result.solve.total_remote_shadow_state_dofs));
        if (session.rank() == 0) {
            compare_reference(reference_path, result.solve.state,
                              field_split || block_jacobi || hypre ? 1.0e-7 : 1.0e-10);
            std::cout << "steady_maximum_shadow_state_dofs=" << result.solve.maximum_shadow_state_dofs << '\n'
                      << "steady_total_remote_shadow_state_dofs=" << result.solve.total_remote_shadow_state_dofs
                      << '\n';
            std::cout << "[PASS] one/two-rank state equivalence"
                      << (field_split    ? " with field split\n"
                          : block_jacobi ? " with block Jacobi\n"
                          : hypre        ? " with hypre\n"
                                         : "\n");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] MPI equivalence benchmark raised: " << error.what() << '\n';
        return 1;
    }
}
