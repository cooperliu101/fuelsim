#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
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
    options.minimum_time_step = std::min(options.minimum_time_step, initial_time_step);
    options.maximum_time_step = initial_time_step;
    options.growth_factor = 1.0;
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
        if (std::abs(step.time - std::round(step.time)) > 1.0e-12) return;
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

bool histories_identical(const std::vector<fuelsim::ContactPointHistory>& actual,
    const std::vector<fuelsim::ContactPointHistory>& expected) {
    if (actual.size() != expected.size()) return false;
    for (std::size_t point = 0; point < actual.size(); ++point)
        if (actual[point].elastic_tangential_slip != expected[point].elastic_tangential_slip ||
            actual[point].sliding != expected[point].sliding ||
            actual[point].normal_multiplier != expected[point].normal_multiplier ||
            actual[point].cartesian_elastic_tangential_slip != expected[point].cartesian_elastic_tangential_slip ||
            actual[point].cartesian_total_tangential_slip != expected[point].cartesian_total_tangential_slip ||
            actual[point].cartesian_tangent_basis_initialized != expected[point].cartesian_tangent_basis_initialized ||
            actual[point].cartesian_contact_normal != expected[point].cartesian_contact_normal ||
            actual[point].cartesian_contact_tangent_first != expected[point].cartesian_contact_tangent_first)
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

bool run_path(const std::string& input_path, const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    Recorder full_recorder;
    fuelsim::TransientProblem full(input.spatial, mesh);
    const fuelsim::TransientResult full_result = fuelsim::solve_transient(
        full, time_options(input, input.transient_execution.end_time, 0.1), solver, &full_recorder);
    bool passed = check(full_result.completed && full_recorder.states.size() == 7,
        "H20 friction path completes all seven prescribed load states");
    if (full_recorder.states.size() == 7) {
        const auto first = stick_slide_counts(full_recorder.states[0]);
        const auto second = stick_slide_counts(full_recorder.states[1]);
        const auto forward_entry = stick_slide_counts(full_recorder.states[2]);
        const auto forward = stick_slide_counts(full_recorder.states[3]);
        const auto unload = stick_slide_counts(full_recorder.states[4]);
        const auto reverse = stick_slide_counts(full_recorder.states[5]);
        const auto restick = stick_slide_counts(full_recorder.states[6]);
        passed = check(first[0] == 13 && first[1] == 0 && second[0] == 13 && second[1] == 0,
                     "H20 friction path starts with two fully sticking states") &&
                 check(forward_entry[0] == 0 && forward_entry[1] == 13 && forward[0] == 0 && forward[1] == 13,
                     "H20 friction path enters and remains in committed forward sliding") &&
                 check(unload[0] == 0 && unload[1] == 13 && reverse[0] == 0 && reverse[1] == 13 &&
                           tangential_resultant_y(full_recorder.states[4]) < 0.0 &&
                           tangential_resultant_y(full_recorder.states[5]) < 0.0,
                     "H20 friction path covers committed unloading and reverse sliding with reversed signed force") &&
                 check(restick[0] == 13 && restick[1] == 0,
                     "H20 friction path finishes with all thirteen constraints restuck") &&
                 passed;
    }

    Recorder first_recorder;
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(input, 4.0, 0.1), solver, &first_recorder);
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
    if (argc != 3) return 2;
    try {
        fuelsim::PetscSession session(argc, argv, "HEX20 internal friction transaction contract");
        return run_path(argv[1], argv[2]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
