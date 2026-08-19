#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/petsc_solver.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <exodusII.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool expect_failure(const std::function<void()>& function, const std::string& expected, const std::string& message) {
    try {
        function();
    } catch (const std::exception& error) {
        return check(std::string(error.what()).find(expected) != std::string::npos, message);
    }
    return check(false, message);
}

bool nearly_equal(double left, double right) {
    const double scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= 2.0e-13 * scale;
}

bool compare_committed_states(
    const fuelsim::TransientCommittedState& left, const fuelsim::TransientCommittedState& right) {
    const std::array<double, 18> left_conservation = {left.conservation.generated_heat_rate,
        left.conservation.stored_heat_rate, left.conservation.convection_heat_rate,
        left.conservation.interface_heat_imbalance, left.conservation.dirichlet_heat_input_rate,
        left.conservation.global_thermal_balance, left.conservation.relative_thermal_balance,
        left.conservation.unconstrained_thermal_residual_l2, left.conservation.internal_mechanical_work_increment,
        left.conservation.pressure_traction_work_increment, left.conservation.dirichlet_reaction_work_increment,
        left.conservation.contact_work_increment, left.conservation.mechanical_work_balance,
        left.conservation.relative_mechanical_work_balance, left.conservation.unconstrained_mechanical_residual_l2,
        left.conservation.elastic_energy_change, left.conservation.plastic_dissipation_increment,
        left.conservation.creep_dissipation_increment};
    const std::array<double, 18> right_conservation = {right.conservation.generated_heat_rate,
        right.conservation.stored_heat_rate, right.conservation.convection_heat_rate,
        right.conservation.interface_heat_imbalance, right.conservation.dirichlet_heat_input_rate,
        right.conservation.global_thermal_balance, right.conservation.relative_thermal_balance,
        right.conservation.unconstrained_thermal_residual_l2, right.conservation.internal_mechanical_work_increment,
        right.conservation.pressure_traction_work_increment, right.conservation.dirichlet_reaction_work_increment,
        right.conservation.contact_work_increment, right.conservation.mechanical_work_balance,
        right.conservation.relative_mechanical_work_balance, right.conservation.unconstrained_mechanical_residual_l2,
        right.conservation.elastic_energy_change, right.conservation.plastic_dissipation_increment,
        right.conservation.creep_dissipation_increment};
    bool passed =
        check(nearly_equal(left.time, right.time) && nearly_equal(left.load_factor, right.load_factor),
            "restart preserves committed time and load") &&
        check(std::equal(left_conservation.begin(), left_conservation.end(), right_conservation.begin(), nearly_equal),
            "restart preserves the last conservation summary") &&
        check(left.solution.size() == right.solution.size(), "restart preserves nodal-state layout") &&
        check(left.material_histories.size() == right.material_histories.size(),
            "restart preserves material-region layout") &&
        check(left.contact_histories.size() == right.contact_histories.size(),
            "restart preserves contact-history layout");
    if (!passed) return false;
    for (std::size_t dof = 0; dof < left.solution.size(); ++dof)
        passed = check(nearly_equal(left.solution[dof], right.solution[dof]),
                     "restart reproduces the uninterrupted nodal state") &&
                 passed;
    for (std::size_t region = 0; region < left.material_histories.size(); ++region) {
        if (!check(left.material_histories[region].size() == right.material_histories[region].size(),
                "restart preserves material-element layout"))
            return false;
        for (std::size_t element = 0; element < left.material_histories[region].size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const fuelsim::MaterialPointState& a = left.material_histories[region][element][q];
                const fuelsim::MaterialPointState& b = right.material_histories[region][element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    passed = check(nearly_equal(a.elastic_strain[component], b.elastic_strain[component]) &&
                                       nearly_equal(a.plastic_strain[component], b.plastic_strain[component]) &&
                                       nearly_equal(a.creep_strain[component], b.creep_strain[component]),
                                 "restart reproduces committed tensor history") &&
                             passed;
                }
                passed = check(nearly_equal(a.equivalent_plastic_strain, b.equivalent_plastic_strain) &&
                                   nearly_equal(a.equivalent_creep_strain, b.equivalent_creep_strain),
                             "restart reproduces committed scalar history") &&
                         passed;
                const fuelsim::AxisymmetricStressValues& stress_a = a.stress;
                const fuelsim::AxisymmetricStressValues& stress_b = b.stress;
                passed = check(nearly_equal(stress_a.rr, stress_b.rr) && nearly_equal(stress_a.zz, stress_b.zz) &&
                                   nearly_equal(stress_a.hoop, stress_b.hoop) && nearly_equal(stress_a.rz, stress_b.rz),
                             "restart reproduces committed stresses") &&
                         passed;
            }
        }
    }
    for (std::size_t contact = 0; contact < left.contact_histories.size(); ++contact) {
        if (!check(left.contact_histories[contact].size() == right.contact_histories[contact].size(),
                "restart preserves contact-node history layout"))
            return false;
        for (std::size_t node = 0; node < left.contact_histories[contact].size(); ++node) {
            const fuelsim::ContactPointHistory& a = left.contact_histories[contact][node];
            const fuelsim::ContactPointHistory& b = right.contact_histories[contact][node];
            passed = check(nearly_equal(a.elastic_tangential_slip, b.elastic_tangential_slip) &&
                               a.sliding == b.sliding && nearly_equal(a.normal_multiplier, b.normal_multiplier),
                         "restart reproduces committed friction and normal "
                         "multiplier history") &&
                     passed;
        }
    }
    return passed;
}

fuelsim::TransientTimeOptions time_options(double end_time) { return {end_time, 1.0, 0.125, 1.0, 1.0, 0.5, 3, 20.0}; }

fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& input, double end_time) {
    const fuelsim::TransientTimeOptions& execution = input.transient_execution;
    return {end_time, execution.initial_time_step, execution.minimum_time_step, execution.maximum_time_step,
        execution.growth_factor, execution.cutback_factor, execution.maximum_cutbacks_per_step,
        execution.load_ramp_time, execution.target_nonlinear_iterations, execution.iteration_window,
        execution.time_error_relative_tolerance, execution.temperature_time_absolute_tolerance,
        execution.displacement_time_absolute_tolerance, execution.time_error_safety_factor,
        execution.strain_history_time_absolute_tolerance, execution.stress_history_time_absolute_tolerance};
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& input) {
    fuelsim::SolverOptions options{input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    options.backtracking_fallback = input.solver.backtracking_fallback;
    options.field_residual_scaling = input.solver.field_residual_scaling;
    options.residual_reduction_tolerance = input.solver.residual_reduction_tolerance;
    options.temperature_residual_absolute_tolerance = input.solver.temperature_residual_absolute_tolerance;
    options.mechanical_residual_absolute_tolerance = input.solver.mechanical_residual_absolute_tolerance;
    options.temperature_residual_scale = input.solver.temperature_residual_scale;
    options.mechanical_residual_scale = input.solver.mechanical_residual_scale;
    return options;
}

bool test_friction_history_checkpoint(const std::string& input_path, const std::string& checkpoint_path) {
    fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    input.spatial.contacts.at(0).friction_coefficient = 0.3;
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::TransientProblem source(input.spatial, mesh);
    fuelsim::TransientCommittedState state = fuelsim::rz::ProblemAccess::committed_state(source);
    if (state.contact_histories.empty() || state.contact_histories.front().empty())
        return check(false, "friction checkpoint fixture has contact-node history");
    state.contact_histories.front().front() = {2.5e-7, true};
    fuelsim::rz::ProblemAccess::restore_committed_state(source, state);
    const fuelsim::TransientCommittedState before_rollback = fuelsim::rz::ProblemAccess::committed_state(source);
    source.begin_time_step({1.0, 0.05});
    source.rollback_time_step();
    bool passed = compare_committed_states(before_rollback, fuelsim::rz::ProblemAccess::committed_state(source));
    fuelsim::write_transient_checkpoint(checkpoint_path, source, 0.5);
    fuelsim::TransientProblem restored(input.spatial, mesh);
    const double next_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    passed = check(next_time_step == 0.5, "friction checkpoint preserves the controller time step") &&
             compare_committed_states(fuelsim::rz::ProblemAccess::committed_state(source),
                 fuelsim::rz::ProblemAccess::committed_state(restored)) &&
             passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, source, 0.5);
    {
        std::fstream file(checkpoint_path, std::ios::binary | std::ios::in | std::ios::out);
        if (!file) return check(false, "friction checkpoint opens for version test");
        const std::array<unsigned char, 4> old_version = {5U, 0U, 0U, 0U};
        file.seekp(16, std::ios::beg);
        file.write(reinterpret_cast<const char*>(old_version.data()), static_cast<std::streamsize>(old_version.size()));
    }
    fuelsim::TransientProblem old_version_target(input.spatial, mesh);
    passed = expect_failure([&]() { (void)fuelsim::restore_transient_checkpoint(checkpoint_path, old_version_target); },
                 "version is not supported", "checkpoint version 6 rejects the previous format") &&
             passed;
    return check(std::remove(checkpoint_path.c_str()) == 0, "friction checkpoint artifact is removed") && passed;
}

std::size_t contact_secondary_global_node(const fuelsim::TransientProblem& problem, std::size_t source_node) {
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const std::vector<std::size_t>& source_nodes =
            fuelsim::rz::ProblemAccess::region_mesh(problem, region).source_node_ids();
        const auto found = std::find(source_nodes.begin(), source_nodes.end(), source_node);
        if (found != source_nodes.end())
            return fuelsim::rz::ProblemAccess::region_node_offset(problem, region) +
                   static_cast<std::size_t>(found - source_nodes.begin());
    }
    throw std::logic_error("Augmented-contact checkpoint secondary node mapping failed");
}

bool test_augmented_contact_transaction(const std::string& input_path, const std::string& checkpoint_path) {
    fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    fuelsim::ContactDefinition& contact = input.spatial.contacts.at(0);
    contact.mechanical_formulation = fuelsim::MechanicalContactFormulation::augmented_lagrangian;
    contact.automatic_penalty = false;
    contact.penalty = 1.0e14;
    contact.penetration_tolerance = 1.0e-9;
    contact.maximum_augmented_iterations = 10;
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    fuelsim::TransientProblem source(input.spatial, mesh);
    const fuelsim::TransientCommittedState initial = fuelsim::rz::ProblemAccess::committed_state(source);
    std::vector<double> penetrated = source.committed_solution();
    const std::vector<fuelsim::ContactNodeSummary> initial_nodes =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(source, 0, penetrated);
    const std::vector<std::size_t> secondary_sources =
        fuelsim::rz::ProblemAccess::contact_secondary_source_nodes(source, 0);
    if (initial_nodes.size() != secondary_sources.size() || initial_nodes.empty())
        return check(false, "augmented transaction fixture has contact-node history");
    constexpr double prescribed_penetration = 2.0e-9;
    for (std::size_t node = 0; node < initial_nodes.size(); ++node) {
        const std::size_t global = contact_secondary_global_node(source, secondary_sources[node]);
        penetrated[fuelsim::rz::ProblemAccess::dof_map(source).dof(fuelsim::Field::radial_displacement, global)] +=
            initial_nodes[node].gap + prescribed_penetration;
    }
    source.begin_time_step({1.0, 0.05});
    const fuelsim::AugmentedContactUpdate update = source.update_augmented_contact_multipliers(penetrated, 0);
    const fuelsim::TransientCommittedState trial = fuelsim::rz::ProblemAccess::committed_state(source);
    bool active_multiplier = false;
    for (const fuelsim::ContactPointHistory& history : trial.contact_histories.at(0))
        active_multiplier = active_multiplier || history.normal_multiplier > 0.0;
    bool passed = check(!update.converged && update.update_allowed &&
                            nearly_equal(update.maximum_penetration, prescribed_penetration) && active_multiplier,
                      "augmented outer update creates a positive trial multiplier") &&
                  check(source.committed_time() == 0.0, "augmented outer update does not advance committed time");
    source.rollback_time_step();
    passed = compare_committed_states(initial, fuelsim::rz::ProblemAccess::committed_state(source)) && passed;
    source.begin_time_step({1.0, 0.05});
    (void)source.update_augmented_contact_multipliers(penetrated, 0);
    source.commit_time_step(penetrated);
    const fuelsim::TransientCommittedState committed = fuelsim::rz::ProblemAccess::committed_state(source);
    bool committed_multiplier = false;
    for (const fuelsim::ContactPointHistory& history : committed.contact_histories.at(0))
        committed_multiplier = committed_multiplier || history.normal_multiplier > 0.0;
    passed = check(committed_multiplier, "accepted augmented state commits a positive normal "
                                         "multiplier") &&
             passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, source, 0.25);
    fuelsim::TransientProblem restored(input.spatial, mesh);
    const double next_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    passed = check(next_time_step == 0.25, "augmented checkpoint preserves the controller time step") &&
             compare_committed_states(committed, fuelsim::rz::ProblemAccess::committed_state(restored)) && passed;
    return check(std::remove(checkpoint_path.c_str()) == 0, "augmented checkpoint artifact is removed") && passed;
}

bool test_finite_strain_restart(const std::string& input_path, const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    fuelsim::TransientProblem uninterrupted(input.spatial, mesh);
    const fuelsim::TransientResult full =
        fuelsim::solve_transient(uninterrupted, time_options(input, input.transient_execution.end_time), solver);
    bool passed = check(full.completed, "uninterrupted finite-strain solve completes");
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first = fuelsim::solve_transient(split, time_options(input, 2.5), solver);
    passed = check(first.completed && nearly_equal(split.committed_time(), 2.5),
                 "finite-strain restart split reaches the deformed state") &&
             passed;
    bool active_rotated_history = false;
    for (std::size_t element = 0; element < fuelsim::rz::ProblemAccess::region_mesh(split, 0).elements().size();
        ++element) {
        for (const fuelsim::MaterialPointState& point :
            fuelsim::rz::ProblemAccess::material_history(split, 0, element)) {
            active_rotated_history =
                active_rotated_history ||
                (std::abs(point.plastic_strain[3]) > 1.0e-3 && std::abs(point.creep_strain[3]) > 1.0e-8 &&
                    point.equivalent_plastic_strain > 0.0 && point.equivalent_creep_strain > 0.0);
        }
    }
    passed = check(active_rotated_history, "finite-strain restart state contains rotated plastic and "
                                           "creep histories") &&
             passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step);
    const fuelsim::TransientCommittedState split_state = fuelsim::rz::ProblemAccess::committed_state(split);
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    passed = compare_committed_states(split_state, fuelsim::rz::ProblemAccess::committed_state(restarted)) && passed;
    fuelsim::TransientTimeOptions restart_options = time_options(input, input.transient_execution.end_time);
    restart_options.initial_time_step = restored_time_step;
    const fuelsim::TransientResult second = fuelsim::solve_transient(restarted, restart_options, solver);
    passed = check(second.completed, "restarted finite-strain solve reaches end time") &&
             compare_committed_states(fuelsim::rz::ProblemAccess::committed_state(uninterrupted),
                 fuelsim::rz::ProblemAccess::committed_state(restarted)) &&
             passed;
    return check(std::remove(checkpoint_path.c_str()) == 0, "finite-strain restart artifact is removed") && passed;
}

class ResultsObserver final : public fuelsim::TransientStepObserver {
  public:
    explicit ResultsObserver(fuelsim::ExodusTransientResultsWriter& writer) : _writer(writer), _steps(0) {}

    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep&) override {
        _writer.append(problem);
        ++_steps;
    }

    std::size_t steps() const noexcept { return _steps; }

  private:
    fuelsim::ExodusTransientResultsWriter& _writer;
    std::size_t _steps;
};

bool verify_exodus(const std::string& path, const fuelsim::UnstructuredQuad4Mesh& mesh,
    const fuelsim::TransientProblem& problem, std::size_t expected_steps) {
    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = 0;
    float version = 0.0F;
    const int exoid = ex_open(path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
    if (!check(exoid >= 0, "transient Exodus result can be reopened")) return false;
    ex_set_int64_status(exoid, EX_ALL_INT64_API);
    int nodal_variables = 0;
    int element_variables = 0;
    int global_variables = 0;
    bool passed = check(ex_inquire_int(exoid, EX_INQ_TIME) == static_cast<std::int64_t>(expected_steps),
                      "Exodus stores the initial and every accepted committed step") &&
                  check(ex_get_variable_param(exoid, EX_NODAL, &nodal_variables) == 0 && nodal_variables == 8,
                      "Exodus defines temperature, displacement, gap and pressure") &&
                  check(ex_get_variable_param(exoid, EX_ELEM_BLOCK, &element_variables) == 0 && element_variables == 56,
                      "Exodus defines stress and inelastic integration-point fields") &&
                  check(ex_get_variable_param(exoid, EX_GLOBAL, &global_variables) == 0 && global_variables == 4,
                      "Exodus defines load and conservative interface totals");
    const int last_step = static_cast<int>(expected_steps);
    double time = 0.0;
    std::vector<double> temperatures(mesh.nodes().size(), 0.0);
    passed = check(ex_get_time(exoid, last_step, &time) == 0 && nearly_equal(time, problem.committed_time()),
                 "Exodus last time equals the committed physical time") &&
             check(ex_get_var(exoid, last_step, EX_NODAL, 1, 1, static_cast<std::int64_t>(temperatures.size()),
                       temperatures.data()) == 0,
                 "Exodus temperature field is readable") &&
             passed;
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const fuelsim::RegionMesh& region_mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::rz::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t node = 0; node < region_mesh.nodes().size(); ++node) {
            const std::size_t source = region_mesh.source_node_ids()[node];
            const double expected = problem.committed_solution().at(
                fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, offset + node));
            passed = check(nearly_equal(temperatures.at(source), expected),
                         "Exodus nodal temperature uses source-mesh mapping") &&
                     passed;
        }
    }
    passed = check(ex_close(exoid) == 0, "Exodus result closes cleanly") && passed;
    return passed;
}

bool verify_steady_results(const fuelsim::FuelSimCaseDefinition& input, const fuelsim::UnstructuredQuad4Mesh& mesh,
    const std::string& results_path) {
    fuelsim::SteadyProblem problem(input.spatial, mesh);
    const fuelsim::SolverOptions solver{input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem,
        {input.steady_execution.load_steps, input.steady_execution.cutback_factor,
            input.steady_execution.maximum_cutbacks_per_step, input.steady_execution.minimum_load_increment},
        solver);
    if (!check(result.completed && result.solve.converged, "steady result fixture converges")) return false;
    fuelsim::write_steady_results(results_path, mesh, problem, result.solve.state);
    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = 0;
    float version = 0.0F;
    const int exoid = ex_open(results_path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
    if (!check(exoid >= 0, "steady Exodus result can be reopened")) return false;
    int nodal_variables = 0;
    int element_variables = 0;
    int global_variables = 0;
    const bool passed =
        check(ex_inquire_int(exoid, EX_INQ_TIME) == 1, "steady Exodus result contains one final state") &&
        check(ex_get_variable_param(exoid, EX_NODAL, &nodal_variables) == 0 && nodal_variables == 8,
            "steady Exodus result contains nodal contact fields") &&
        check(ex_get_variable_param(exoid, EX_ELEM_BLOCK, &element_variables) == 0 && element_variables == 16,
            "steady Exodus result contains four-point stresses") &&
        check(ex_get_variable_param(exoid, EX_GLOBAL, &global_variables) == 0 && global_variables == 4,
            "steady Exodus result contains interface totals") &&
        check(ex_close(exoid) == 0, "steady Exodus result closes cleanly");
    return passed;
}

bool run_tests(const std::string& steady_input_path, const std::string& transient_input_path,
    const std::string& finite_strain_input_path, const std::string& checkpoint_path, const std::string& results_path) {
    const fuelsim::FuelSimCaseDefinition steady_input = fuelsim::read_case_input(steady_input_path);
    const std::filesystem::path configured_results(results_path);
    const std::filesystem::path first_segment =
        configured_results.parent_path() /
        (configured_results.stem().string() + ".part1" + configured_results.extension().string());
    {
        std::ofstream occupied(first_segment);
        occupied << "previous segment";
    }
    const std::filesystem::path second_segment =
        configured_results.parent_path() /
        (configured_results.stem().string() + ".part2" + configured_results.extension().string());
    bool passed = check(fuelsim::next_results_segment_path(results_path) == second_segment.string(),
        "restart result segmentation preserves occupied earlier files");
    std::filesystem::remove(first_segment);
    const fuelsim::UnstructuredQuad4Mesh steady_mesh = fuelsim::read_exodus_quad4(steady_input.mesh_file);
    passed = verify_steady_results(steady_input, steady_mesh, results_path) && passed;
    passed = test_friction_history_checkpoint(transient_input_path, checkpoint_path + ".friction") && passed;
    passed = test_augmented_contact_transaction(transient_input_path, checkpoint_path + ".augmented") && passed;
    passed = test_finite_strain_restart(finite_strain_input_path, checkpoint_path + ".finite") && passed;
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(transient_input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    const fuelsim::SolverOptions solver{input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    fuelsim::TransientProblem uninterrupted(input.spatial, mesh);
    const fuelsim::TransientResult full = fuelsim::solve_transient(uninterrupted, time_options(20.0), solver);
    passed = check(full.completed, "uninterrupted PCMI solve completes") && passed;
    fuelsim::TransientProblem split(input.spatial, mesh);
    fuelsim::ExodusTransientResultsWriter writer(results_path, mesh, split);
    writer.append(split);
    ResultsObserver observer(writer);
    const fuelsim::TransientResult first = fuelsim::solve_transient(split, time_options(10.0), solver, &observer);
    passed = check(first.completed && observer.steps() == first.accepted_steps.size(),
                 "first restart segment observes every accepted step") &&
             passed;
    const std::string history_path = results_path + ".history.csv";
    {
        fuelsim::EngineeringHistoryWriter history(history_path, split);
        history.append(split, 1.0, 1.0, first.last_attempt.nonlinear_iterations);
    }
    {
        std::ifstream history(history_path);
        std::string header;
        std::string values;
        std::getline(history, header);
        std::getline(history, values);
        passed = check(header.find("region_cladding_maximum_temperature") != std::string::npos &&
                           header.find("contact_fuel_cladding_minimum_gap") != std::string::npos && !values.empty(),
                     "engineering history uses named region/contact columns") &&
                 passed;
    }
    std::filesystem::remove(history_path);
    fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step);
    const fuelsim::TransientCommittedState split_state = fuelsim::rz::ProblemAccess::committed_state(split);
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    passed =
        check(restored_time_step == first.next_time_step, "restart preserves the committed controller step") && passed;
    passed = compare_committed_states(split_state, fuelsim::rz::ProblemAccess::committed_state(restarted)) && passed;
    fuelsim::TransientTimeOptions restart_options = time_options(20.0);
    restart_options.initial_time_step = restored_time_step;
    const fuelsim::TransientResult second = fuelsim::solve_transient(restarted, restart_options, solver);
    passed = check(second.completed, "restarted PCMI solve reaches end time") &&
             compare_committed_states(fuelsim::rz::ProblemAccess::committed_state(uninterrupted),
                 fuelsim::rz::ProblemAccess::committed_state(restarted)) &&
             passed;
    passed = verify_exodus(results_path, mesh, split, observer.steps() + 1) && passed;
    fuelsim::TransientProblem mismatch(input.spatial, mesh);
    fuelsim::SpatialDefinition changed = input.spatial;
    auto changed_functions = std::make_shared<fuelsim::MaterialFunctionSet>(*changed.regions[1].material.functions);
    ++changed_functions->plasticity.version;
    changed.regions[1].material.functions = std::move(changed_functions);
    fuelsim::TransientProblem changed_problem(std::move(changed), mesh);
    passed = expect_failure([&]() { (void)fuelsim::restore_transient_checkpoint(checkpoint_path, changed_problem); },
                 "signature", "checkpoint rejects a changed registered material function version") &&
             passed;
    mismatch.begin_time_step({1.0, 0.05});
    passed = expect_failure([&]() { fuelsim::write_transient_checkpoint(checkpoint_path, mismatch, 1.0); },
                 "active time step", "checkpoint cannot capture uncommitted trial state") &&
             passed;
    mismatch.rollback_time_step();
    passed = expect_failure([&]() { fuelsim::write_transient_checkpoint(checkpoint_path, mismatch, 0.0); },
                 "next time step", "checkpoint rejects invalid controller state") &&
             passed;
    {
        std::fstream file(checkpoint_path, std::ios::binary | std::ios::in | std::ios::out);
        if (!file) return check(false, "checkpoint can be opened for corruption test");
        file.seekg(48, std::ios::beg);
        char byte = 0;
        file.read(&byte, 1);
        file.clear();
        file.seekp(48, std::ios::beg);
        byte = static_cast<char>(byte ^ 0x5A);
        file.write(&byte, 1);
    }
    passed = expect_failure([&]() { (void)fuelsim::restore_transient_checkpoint(checkpoint_path, mismatch); },
                 "checksum", "checkpoint detects payload corruption") &&
             passed;
    const int checkpoint_remove = std::remove(checkpoint_path.c_str());
    const int results_remove = std::remove(results_path.c_str());
    return check(checkpoint_remove == 0 && results_remove == 0, "M3.0 test artifacts are removed") && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: checkpoint_results_tests <steady.fsi> "
                     "<transient.fsi> <finite-strain.fsi> <checkpoint> "
                     "<results.e>\n";
        return 2;
    }
    try {
        fuelsim::PetscSession session(argc, argv, "fuelsim M3.0 checkpoint/results tests\n");
        if (!run_tests(argv[1], argv[2], argv[3], argv[4], argv[5])) return 1;
        std::cout << "[PASS] fuelsim M3.0 checkpoint and Exodus results\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.0 tests raised: " << error.what() << '\n';
        return 1;
    }
}
