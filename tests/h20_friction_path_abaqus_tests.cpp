#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct DisplacementReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> displacement;
};

struct ContactReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2;
};

struct ReactionReference final {
    std::array<double, 3> reaction;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) result.push_back(field);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("H20.33 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("H20.33 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("H20.33 CSV contains an invalid integer index: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<std::vector<DisplacementReference>> read_displacement_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.33 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,time,id,x,y,z,disp_x,disp_y,disp_z")
        throw std::invalid_argument("Unexpected H20.33 displacement header in " + path);
    std::vector<std::vector<DisplacementReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.33 displacement step is invalid: " + path);
        result[step - 1].push_back(
            {index_value(values, 2, path), {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 6, path), number(values, 7, path), number(values, 8, path)}});
    }
    return result;
}

std::vector<std::vector<ContactReference>> read_contact_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.33 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line !=
        "step,time,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected H20.33 contact header in " + path);
    std::vector<std::vector<ContactReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.33 contact step is invalid: " + path);
        result[step - 1].push_back(
            {index_value(values, 2, path), {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
                {number(values, 9, path), number(values, 10, path), number(values, 11, path)}, number(values, 12, path),
                number(values, 13, path)});
    }
    return result;
}

std::vector<ReactionReference> read_reaction_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.33 reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,time,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected H20.33 reaction header in " + path);
    std::vector<ReactionReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step != result.size() + 1 || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.33 reaction steps are not consecutive: " + path);
        result.push_back({{number(values, 2, path), number(values, 3, path), number(values, 4, path)}});
    }
    return result;
}

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.direct_factorization = definition.solver.direct_factorization;
    options.field_residual_scaling = definition.solver.field_residual_scaling;
    options.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    options.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    return options;
}

fuelsim::TransientTimeOptions time_options(
    const fuelsim::FuelSimCaseDefinition& definition, double end_time, double initial_time_step) {
    fuelsim::TransientTimeOptions options = definition.transient_execution;
    options.end_time = end_time;
    options.initial_time_step = initial_time_step;
    return options;
}

struct StepState final {
    double time;
    std::vector<double> solution;
    std::vector<fuelsim::CartesianContactNodeSummary> contact;
    std::vector<fuelsim::ContactPointHistory> histories;
};

class Recorder final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        const std::vector<double>& solution = problem.committed_solution();
        const std::vector<fuelsim::CartesianContactNodeSummary> contact =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solution);
        std::size_t sticking = 0, sliding = 0;
        for (const fuelsim::CartesianContactNodeSummary& node : contact) {
            if (!(node.pressure > 0.0)) continue;
            if (node.sliding)
                ++sliding;
            else
                ++sticking;
        }
        double normal_x = 0.0, tangential_y = 0.0, tangential_z = 0.0;
        for (const fuelsim::CartesianContactNodeSummary& node : contact) {
            normal_x += node.normal_contact_force[0];
            tangential_y += node.tangential_contact_force[1];
            tangential_z += node.tangential_contact_force[2];
        }
        std::cout << "h20_friction_path_time=" << step.time << " sticking=" << sticking << " sliding=" << sliding
                  << " normal_x=" << normal_x << " tangential_y=" << tangential_y << " tangential_z=" << tangential_z
                  << '\n';
        states.push_back({step.time, solution, contact,
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0)});
    }

    std::vector<StepState> states;
};

bool compare_abaqus_steps(const fuelsim::TransientProblem& problem, const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::vector<StepState>& states, const std::string& displacement_path, const std::string& contact_path,
    const std::string& reaction_path) {
    const auto displacement = read_displacement_reference(displacement_path);
    const auto contact = read_contact_reference(contact_path);
    const auto reaction = read_reaction_reference(reaction_path);
    if (states.size() != displacement.size() || states.size() != contact.size() || states.size() != reaction.size())
        throw std::invalid_argument("H20.33 Fuelsim and Abaqus step counts differ");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    const std::vector<std::size_t> contact_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    constexpr double tolerance = 1.0e-2;
    constexpr double nonmatching_displacement_pointwise_tolerance = 1.2e-1;
    constexpr double nonmatching_force_pointwise_tolerance = 1.5e-2;
    constexpr double zero_tolerance = 1.2e-8;
    bool passed = true;
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_metrics;
        fuelsim::test::FieldErrorMetrics normal_force, tangential_force, slip_1, slip_2;
        std::vector<bool> present(mesh.nodes().size(), false);
        double maximum_coordinate_difference = 0.0;
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& region_mesh = spatial.hex20_region_mesh(region);
            for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
                const std::size_t source = region_mesh.source_node_ids()[local];
                const auto found = std::find_if(displacement[step].begin(), displacement[step].end(),
                    [source](const DisplacementReference& value) { return value.id == source; });
                if (found == displacement[step].end() || present[source])
                    throw std::invalid_argument("H20.33 displacement source-node mapping is incomplete or repeated");
                present[source] = true;
                maximum_coordinate_difference =
                    std::max({maximum_coordinate_difference, std::abs(mesh.nodes()[source].x - found->point.x),
                        std::abs(mesh.nodes()[source].y - found->point.y),
                        std::abs(mesh.nodes()[source].z - found->point.z)});
                const std::size_t global = spatial.global_node(region, local);
                for (std::size_t component = 0; component < 3; ++component)
                    displacement_metrics[component].add(
                        states[step].solution[fields[component + 1].begin + global], found->displacement[component]);
            }
        }
        if (contact_nodes.size() != states[step].contact.size() || contact_nodes.size() != contact[step].size())
            throw std::invalid_argument("H20.33 contact-node counts differ");
        for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
            const auto found = std::find_if(contact[step].begin(), contact[step].end(),
                [&](const ContactReference& value) { return value.id == contact_nodes[node]; });
            if (found == contact[step].end()) throw std::invalid_argument("H20.33 contact-node mapping is incomplete");
            const fuelsim::CartesianContactNodeSummary& actual = states[step].contact[node];
            maximum_coordinate_difference =
                std::max({maximum_coordinate_difference, std::abs(actual.x - found->point.x),
                    std::abs(actual.y - found->point.y), std::abs(actual.z - found->point.z)});
            // Abaqus contact-force output uses the opposite contact-pair action
            // convention from the Fuelsim secondary-side summary.
            normal_force.add(-actual.normal_contact_force[0], found->normal_force[0]);
            tangential_force.add(-actual.tangential_contact_force[1], found->tangential_force[1]);
            slip_1.add(actual.tangential_slip[2], found->slip_1);
            slip_2.add(actual.tangential_slip[1], -found->slip_2);
        }
        const std::string prefix = "h20_33_step_" + std::to_string(step + 1) + "_";
        for (std::size_t component = 0; component < 3; ++component) {
            const std::string name = prefix + "displacement_" + std::string(1, "xyz"[component]);
            if (displacement_metrics[component].has_relative_norm()) {
                fuelsim::test::print_relative_metrics(name, displacement_metrics[component]);
                const bool metric_passed =
                    component == 0
                        ? fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                              displacement_metrics[component], tolerance, nonmatching_displacement_pointwise_tolerance)
                        : fuelsim::test::relative_metrics_below(displacement_metrics[component], tolerance);
                passed =
                    check(metric_passed,
                        "H20.33 step " + std::to_string(step + 1) +
                            (component == 0
                                    ? " normal displacement passes 1 percent aggregate metrics and the "
                                      "explicit 12 percent nonmatching-transition pointwise gate"
                                    : " tangential displacement passes all three Abaqus metrics below 1 percent")) &&
                    passed;
            } else {
                fuelsim::test::print_absolute_metrics(name, displacement_metrics[component]);
                passed = check(displacement_metrics[component].maximum_zero_reference_difference < zero_tolerance,
                             "H20.33 step " + std::to_string(step + 1) +
                                 " theoretical-zero displacement component passes its absolute check") &&
                         passed;
            }
        }
        fuelsim::test::print_relative_metrics(prefix + "signed_normal_force_x", normal_force);
        fuelsim::test::print_relative_metrics(prefix + "signed_tangential_force_y", tangential_force);
        if (slip_1.has_relative_norm())
            fuelsim::test::print_relative_metrics(prefix + "tangential_slip_1", slip_1);
        else
            fuelsim::test::print_absolute_metrics(prefix + "tangential_slip_1", slip_1);
        fuelsim::test::print_relative_metrics(prefix + "tangential_slip_2", slip_2);
        passed = check(maximum_coordinate_difference < 1.0e-7,
                     "H20.33 step " + std::to_string(step + 1) + " uses the tracked nonmatching Exodus mesh") &&
                 check(fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                           normal_force, tolerance, nonmatching_force_pointwise_tolerance),
                     "H20.33 step " + std::to_string(step + 1) +
                         " signed normal nodal force passes 1 percent aggregate metrics and the explicit 1.5 percent "
                         "nonmatching-transfer pointwise gate") &&
                 check(fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                           tangential_force, tolerance, nonmatching_force_pointwise_tolerance),
                     "H20.33 step " + std::to_string(step + 1) +
                         " signed tangential nodal force passes 1 percent aggregate metrics and the explicit 1.5 "
                         "percent nonmatching-transfer pointwise gate") &&
                 check((slip_1.has_relative_norm() ? fuelsim::test::relative_metrics_below(slip_1, tolerance)
                                                   : slip_1.maximum_zero_reference_difference < zero_tolerance) &&
                           fuelsim::test::relative_metrics_below(slip_2, tolerance),
                     "H20.33 step " + std::to_string(step + 1) +
                         " both Abaqus tangential-slip components pass all three metrics below 1 percent") &&
                 passed;
        double actual_normal = 0.0, actual_tangential = 0.0;
        for (const fuelsim::CartesianContactNodeSummary& node : states[step].contact) {
            actual_normal += node.normal_contact_force[0];
            actual_tangential += node.tangential_contact_force[1];
        }
        const double normal_resultant_error =
                         std::abs(actual_normal - reaction[step].reaction[0]) / std::abs(reaction[step].reaction[0]),
                     tangential_resultant_error = std::abs(actual_tangential - reaction[step].reaction[1]) /
                                                  std::abs(reaction[step].reaction[1]);
        std::cout << prefix << "normal_resultant_relative_error=" << normal_resultant_error << '\n'
                  << prefix << "tangential_resultant_relative_error=" << tangential_resultant_error << '\n';
        passed = check(normal_resultant_error < tolerance && tangential_resultant_error < tolerance,
                     "H20.33 step " + std::to_string(step + 1) +
                         " signed normal and tangential resultants agree with Abaqus below 1 percent") &&
                 passed;
    }
    return passed;
}

bool histories_identical(const std::vector<fuelsim::ContactPointHistory>& actual,
    const std::vector<fuelsim::ContactPointHistory>& expected) {
    if (actual.size() != expected.size()) return false;
    for (std::size_t point = 0; point < actual.size(); ++point)
        if (actual[point].elastic_tangential_slip != expected[point].elastic_tangential_slip ||
            actual[point].sliding != expected[point].sliding ||
            actual[point].normal_multiplier != expected[point].normal_multiplier ||
            actual[point].cartesian_elastic_tangential_slip != expected[point].cartesian_elastic_tangential_slip)
            return false;
    return true;
}

bool recorded_states_identical(
    const std::vector<StepState>& actual, const std::vector<StepState>& expected, std::size_t expected_begin) {
    if (actual.size() + expected_begin > expected.size()) return false;
    for (std::size_t step = 0; step < actual.size(); ++step)
        if (actual[step].time != expected[expected_begin + step].time ||
            actual[step].solution != expected[expected_begin + step].solution ||
            !histories_identical(actual[step].histories, expected[expected_begin + step].histories))
            return false;
    return true;
}

std::array<std::size_t, 2> stick_slide_counts(const StepState& state) {
    std::array<std::size_t, 2> result{};
    for (const fuelsim::CartesianContactNodeSummary& node : state.contact) {
        if (!(node.pressure > 0.0)) continue;
        ++result[node.sliding ? 1 : 0];
    }
    return result;
}

double tangential_resultant_y(const StepState& state) {
    double result = 0.0;
    for (const fuelsim::CartesianContactNodeSummary& node : state.contact) result += node.tangential_contact_force[1];
    return result;
}

bool run_path(const std::string& input_path, const std::string& displacement_path, const std::string& contact_path,
    const std::string& reaction_path, const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    Recorder full_recorder;
    fuelsim::TransientProblem full(input.spatial, mesh);
    const fuelsim::TransientResult full_result = fuelsim::solve_transient(
        full, time_options(input, input.transient_execution.end_time, 1.0), solver, &full_recorder);
    bool passed = check(full_result.completed && full_recorder.states.size() == 7,
        "H20 friction path completes all seven prescribed load states");
    if (full_recorder.states.size() == 7) {
        const auto first = stick_slide_counts(full_recorder.states[0]);
        const auto second = stick_slide_counts(full_recorder.states[1]);
        const auto mixed = stick_slide_counts(full_recorder.states[2]);
        const auto forward = stick_slide_counts(full_recorder.states[3]);
        const auto unload = stick_slide_counts(full_recorder.states[4]);
        const auto reverse = stick_slide_counts(full_recorder.states[5]);
        const auto restick = stick_slide_counts(full_recorder.states[6]);
        passed = check(first[0] == 13 && first[1] == 0 && second[0] == 13 && second[1] == 0,
                     "H20 friction path starts with two fully sticking states") &&
                 check(mixed[0] > 0 && mixed[1] > 0 && forward[1] > 0,
                     "H20 friction path contains simultaneous sticking and sliding before forward sliding") &&
                 check(unload[1] > 0 && reverse[1] > 0 && tangential_resultant_y(full_recorder.states[4]) < 0.0 &&
                           tangential_resultant_y(full_recorder.states[5]) < 0.0,
                     "H20 friction path covers unloading and reverse sliding with reversed signed force") &&
                 check(restick[0] == 13 && restick[1] == 0,
                     "H20 friction path finishes with all thirteen constraints restuck") &&
                 passed;
    }
    passed = compare_abaqus_steps(full, mesh, full_recorder.states, displacement_path, contact_path, reaction_path) &&
             passed;

    Recorder first_recorder;
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(input, 4.0, 1.0), solver, &first_recorder);
    passed = check(first.completed && first_recorder.states.size() == 4,
                 "H20 friction path reaches its nonzero-slip checkpoint state") &&
             check(recorded_states_identical(first_recorder.states, full_recorder.states, 0),
                 "H20 friction split run reproduces the first four uninterrupted nodal and contact-history states") &&
             passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step);
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    Recorder second_recorder;
    const fuelsim::TransientResult second = fuelsim::solve_transient(
        restarted, time_options(input, input.transient_execution.end_time, restored_step), solver, &second_recorder);
    passed = check(second.completed && second_recorder.states.size() == 3,
                 "H20 friction checkpoint restores and completes the reversed path") &&
             check(recorded_states_identical(second_recorder.states, full_recorder.states, 4),
                 "H20 friction restart reproduces every remaining nodal and contact-history state exactly") &&
             passed;
    double maximum_restart_difference = 0.0;
    const std::vector<double>& expected = full.committed_solution();
    const std::vector<double>& actual = restarted.committed_solution();
    if (expected.size() != actual.size())
        passed = check(false, "H20 friction restart preserves the global state layout") && passed;
    else
        for (std::size_t dof = 0; dof < expected.size(); ++dof)
            maximum_restart_difference = std::max(maximum_restart_difference, std::abs(expected[dof] - actual[dof]));
    passed =
        check(maximum_restart_difference < 1.0e-13, "H20 friction restart reproduces the uninterrupted final state") &&
        check(std::remove(checkpoint_path.c_str()) == 0, "H20 friction checkpoint artifact is removed") && passed;
    std::cout << "h20_friction_restart_maximum_absolute_difference=" << maximum_restart_difference << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: fuelsim_h20_friction_path_abaqus_tests <case.fsi> <displacement.csv> <contact.csv> "
                     "<reaction.csv> <checkpoint>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim HEX20 Abaqus friction-path comparison\n");
        if (!run_path(argv[1], argv[2], argv[3], argv[4], argv[5])) return 1;
        std::cout << "[PASS] fuelsim HEX20 friction path and restart\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20 friction path raised: " << error.what() << '\n';
        return 1;
    }
}
