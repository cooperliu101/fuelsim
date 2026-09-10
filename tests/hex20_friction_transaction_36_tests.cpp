#include "contact_types.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
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
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double dot(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

std::array<double, 3> subtract(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

std::array<double, 3> cross(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return {first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

std::array<double, 3> unit(const std::array<double, 3>& value) {
    const double measure = std::sqrt(dot(value, value));
    if (!(measure > 0.0))
        throw std::invalid_argument("H20.36 contact face is degenerate");
    return {value[0] / measure, value[1] / measure, value[2] / measure};
}

double maximum_secondary_face_nonplanarity(const fuelsim::UnstructuredHex20Mesh& mesh) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {{{{0, 1, 5, 4, 8, 13, 16, 12}},
        {{1, 2, 6, 5, 9, 14, 17, 13}},
        {{2, 3, 7, 6, 10, 15, 18, 14}},
        {{3, 0, 4, 7, 11, 12, 19, 15}},
        {{0, 3, 2, 1, 11, 10, 9, 8}},
        {{4, 5, 6, 7, 16, 17, 18, 19}}}};
    double result = 0.0;
    for (const fuelsim::ElementSide& side : mesh.side_set("secondary_contact").sides) {
        const auto& element = mesh.elements().at(side.element);
        const auto& local_nodes = face_nodes.at(side.local_side);
        const fuelsim::CartesianPoint3& origin = mesh.nodes().at(element.nodes[local_nodes[0]]);
        const std::array<double, 3> normal =
            unit(cross(subtract(mesh.nodes().at(element.nodes[local_nodes[1]]), origin),
                subtract(mesh.nodes().at(element.nodes[local_nodes[3]]), origin)));
        for (std::size_t local : local_nodes)
            result = std::max(result, std::abs(dot(subtract(mesh.nodes().at(element.nodes[local]), origin), normal)));
    }
    return result;
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

fuelsim::TransientTimeOptions
time_options(const fuelsim::FuelSimCaseDefinition& definition, double end_time, double initial_time_step) {
    fuelsim::TransientTimeOptions options = definition.transient_execution;
    options.end_time = end_time;
    options.initial_time_step = initial_time_step;
    options.minimum_time_step = std::min(options.minimum_time_step, initial_time_step);
    options.maximum_time_step = initial_time_step;
    options.growth_factor = 1.0;
    return options;
}

struct ContactNumericalEvidence final {
    std::array<double, 3> residual_balance{};
    double maximum_directional_error = 0.0;
};

ContactNumericalEvidence inspect_contact_transition(fuelsim::TransientProblem& problem,
    fuelsim::TransientCommittedState previous,
    double target_time,
    const std::vector<double>& state) {
    fuelsim::cartesian::ProblemAccess::restore_committed_state(problem, std::move(previous));
    problem.begin_time_step({target_time, 1.0, false});
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    ContactNumericalEvidence result;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size()), direction(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) {
            local[index] = state[dofs[index]];
            direction[index] = std::sin(static_cast<double>(index + 1));
        }
        std::vector<double> residual, jacobian;
        problem.compute_contribution(contribution, local, residual, &jacobian);
        constexpr double perturbation = 1.0e-10;
        std::vector<double> plus = local, minus = local;
        for (std::size_t index = 0; index < local.size(); ++index) {
            plus[index] += perturbation * direction[index];
            minus[index] -= perturbation * direction[index];
        }
        std::vector<double> plus_residual, minus_residual;
        problem.compute_contribution(contribution, plus, plus_residual, nullptr);
        problem.compute_contribution(contribution, minus, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                analytic += jacobian[row * local.size() + column] * direction[column];
            const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += (analytic - reference) * (analytic - reference);
            reference_squared += reference * reference;
        }
        if (!(reference_squared > 0.0))
            throw std::logic_error("H20.36 contact Jacobian direction is zero");
        result.maximum_directional_error =
            std::max(result.maximum_directional_error, std::sqrt(difference_squared / reference_squared));
        for (std::size_t row = 0; row < dofs.size(); ++row)
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (dofs[row] >= field.begin && dofs[row] < field.end)
                    result.residual_balance[component] += residual[row];
            }
    }
    problem.rollback_time_step();
    return result;
}

struct StepState final {
    double time;
    std::vector<double> solution;
    std::vector<fuelsim::CartesianContactNodeSummary> contact;
    std::vector<fuelsim::ContactPointHistory> histories;
    fuelsim::TransientCommittedState committed;
};

class Recorder final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        if (std::abs(step.time - std::round(step.time)) > 1.0e-12)
            return;
        const std::vector<double>& solution = problem.committed_solution();
        const std::vector<fuelsim::CartesianContactNodeSummary> contact =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solution);
        std::size_t sticking = 0, sliding = 0;
        for (const fuelsim::CartesianContactNodeSummary& node : contact) {
            if (!(node.pressure > 0.0))
                continue;
            if (node.sliding)
                ++sliding;
            else
                ++sticking;
        }
        double normal_radial = 0.0, tangential_circumferential = 0.0, tangential_z = 0.0;
        for (const fuelsim::CartesianContactNodeSummary& node : contact) {
            const double radius = std::hypot(node.x, node.y);
            normal_radial += (node.x * node.normal_contact_force[0] + node.y * node.normal_contact_force[1]) / radius;
            tangential_circumferential +=
                (-node.y * node.tangential_contact_force[0] + node.x * node.tangential_contact_force[1]) / radius;
            tangential_z += node.tangential_contact_force[2];
        }
        std::cout << "h20_36_time=" << step.time << " sticking=" << sticking << " sliding=" << sliding
                  << " normal_radial=" << normal_radial << " tangential_circumferential=" << tangential_circumferential
                  << " tangential_axial=" << tangential_z << '\n';
        states.push_back({step.time,
            solution,
            contact,
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0),
            fuelsim::cartesian::ProblemAccess::committed_state(problem)});
    }

    std::vector<StepState> states;
};

bool histories_identical(const std::vector<fuelsim::ContactPointHistory>& actual,
    const std::vector<fuelsim::ContactPointHistory>& expected) {
    if (actual.size() != expected.size())
        return false;
    for (std::size_t point = 0; point < actual.size(); ++point)
        if (actual[point].elastic_tangential_slip != expected[point].elastic_tangential_slip
            || actual[point].sliding != expected[point].sliding
            || actual[point].normal_multiplier != expected[point].normal_multiplier
            || actual[point].cartesian_elastic_tangential_slip != expected[point].cartesian_elastic_tangential_slip
            || actual[point].cartesian_total_tangential_slip != expected[point].cartesian_total_tangential_slip
            || actual[point].cartesian_tangent_basis_initialized != expected[point].cartesian_tangent_basis_initialized
            || actual[point].cartesian_contact_normal != expected[point].cartesian_contact_normal
            || actual[point].cartesian_contact_tangent_first != expected[point].cartesian_contact_tangent_first)
            return false;
    return true;
}

bool recorded_states_identical(const std::vector<StepState>& actual,
    const std::vector<StepState>& expected,
    std::size_t expected_begin) {
    if (actual.size() + expected_begin > expected.size())
        return false;
    for (std::size_t step = 0; step < actual.size(); ++step)
        if (actual[step].time != expected[expected_begin + step].time
            || actual[step].solution != expected[expected_begin + step].solution
            || !histories_identical(actual[step].histories, expected[expected_begin + step].histories))
            return false;
    return true;
}

std::array<std::size_t, 2> stick_slide_counts(const StepState& state) {
    std::array<std::size_t, 2> result{};
    for (const fuelsim::CartesianContactNodeSummary& node : state.contact) {
        if (!(node.pressure > 0.0))
            continue;
        ++result[node.sliding ? 1 : 0];
    }
    return result;
}

std::array<double, 2> tangential_resultant(const StepState& state) {
    std::array<double, 2> result{};
    for (const fuelsim::CartesianContactNodeSummary& node : state.contact) {
        const double radius = std::hypot(node.x, node.y);
        result[0] += (-node.y * node.tangential_contact_force[0] + node.x * node.tangential_contact_force[1]) / radius;
        result[1] += node.tangential_contact_force[2];
    }
    return result;
}

bool run_path(const std::string& input_path, const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    Recorder full_recorder;
    fuelsim::TransientProblem full(input.spatial, mesh);
    const fuelsim::TransientCommittedState initial_state = fuelsim::cartesian::ProblemAccess::committed_state(full);
    const fuelsim::TransientResult full_result = fuelsim::solve_transient(full,
        time_options(input, input.transient_execution.end_time, 0.1),
        solver,
        &full_recorder);
    bool passed = check(full_result.completed && full_recorder.states.size() == 7,
        "H20.36 curved friction path completes all seven prescribed load states");
    if (full_recorder.states.size() == 7) {
        const auto first = stick_slide_counts(full_recorder.states[0]);
        const auto second = stick_slide_counts(full_recorder.states[1]);
        const auto mixed = stick_slide_counts(full_recorder.states[2]);
        const auto forward = stick_slide_counts(full_recorder.states[3]);
        const auto unload = stick_slide_counts(full_recorder.states[4]);
        const auto reverse = stick_slide_counts(full_recorder.states[5]);
        const auto restick = stick_slide_counts(full_recorder.states[6]);
        const std::array<double, 2> forward_resultant = tangential_resultant(full_recorder.states[3]);
        const std::array<double, 2> reverse_resultant = tangential_resultant(full_recorder.states[5]);
        passed = check(std::all_of(full_recorder.states.begin(),
                           full_recorder.states.end(),
                           [](const StepState& state) {
                               const auto count = stick_slide_counts(state);
                               return count[0] + count[1] == 37;
                           }),
                     "H20.36 keeps all 37 curved constraints closed throughout the path")
                 && check(first[0] == 37 && first[1] == 0 && second[0] == 37 && second[1] == 0,
                     "H20.36 starts with two fully sticking states")
                 && check(mixed[0] > 0 && mixed[1] > 0 && forward[1] > 0,
                     "H20.36 contains simultaneous sticking and sliding before forward sliding")
                 && check(unload[0] > 0 && unload[1] > 0 && reverse[0] == 0 && reverse[1] == 37
                              && forward_resultant[1] * reverse_resultant[1] < 0.0,
                     "H20.36 covers mixed unloading, committed full reverse sliding, and reverses the driven axial "
                     "tangential resultant")
                 && check(restick[0] == 37 && restick[1] == 0, "H20.36 finishes with all 37 constraints restuck")
                 && check(maximum_secondary_face_nonplanarity(mesh) > 1.0e-3,
                     "H20.36 uses genuinely quadratic contact faces rather than planar facets")
                 && passed;
    }

    if (full_recorder.states.size() == 7) {
        fuelsim::TransientProblem sticking_probe(input.spatial, mesh);
        const ContactNumericalEvidence sticking =
            inspect_contact_transition(sticking_probe, initial_state, 1.0, full_recorder.states[0].solution);
        fuelsim::TransientProblem sliding_probe(input.spatial, mesh);
        const ContactNumericalEvidence sliding = inspect_contact_transition(sliding_probe,
            full_recorder.states[4].committed,
            6.0,
            full_recorder.states[5].solution);
        std::cout << "h20_36_sticking_contact_jacobian_directional_error=" << sticking.maximum_directional_error << '\n'
                  << "h20_36_sliding_contact_jacobian_directional_error=" << sliding.maximum_directional_error << '\n'
                  << "h20_36_sticking_contact_residual_balance=" << sticking.residual_balance[0] << ','
                  << sticking.residual_balance[1] << ',' << sticking.residual_balance[2] << '\n'
                  << "h20_36_sliding_contact_residual_balance=" << sliding.residual_balance[0] << ','
                  << sliding.residual_balance[1] << ',' << sliding.residual_balance[2] << '\n';
        passed = check(sticking.maximum_directional_error < 1.0e-7 && sliding.maximum_directional_error < 1.0e-7,
                     "H20.36 curved sticking and sliding contact Jacobians match centered directional differences")
                 && check(std::all_of(sticking.residual_balance.begin(),
                              sticking.residual_balance.end(),
                              [](double value) { return std::abs(value) < 1.0e-8; })
                              && std::all_of(sliding.residual_balance.begin(),
                                  sliding.residual_balance.end(),
                                  [](double value) { return std::abs(value) < 1.0e-8; }),
                     "H20.36 curved sticking and sliding residuals are action-reaction conservative")
                 && passed;
    }

    Recorder first_recorder;
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(input, 4.0, 0.1), solver, &first_recorder);
    passed =
        check(first.completed && first_recorder.states.size() == 4, "H20.36 reaches its nonzero-slip checkpoint state")
        && check(recorded_states_identical(first_recorder.states, full_recorder.states, 0),
            "H20.36 split run reproduces the first four uninterrupted nodal and contact-history states")
        && passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step);
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    Recorder second_recorder;
    const fuelsim::TransientResult second = fuelsim::solve_transient(restarted,
        time_options(input, input.transient_execution.end_time, restored_step),
        solver,
        &second_recorder);
    passed = check(second.completed && second_recorder.states.size() == 3,
                 "H20.36 checkpoint restores and completes the reversed path")
             && check(recorded_states_identical(second_recorder.states, full_recorder.states, 4),
                 "H20.36 restart reproduces every remaining nodal and contact-history state exactly")
             && passed;
    double maximum_restart_difference = 0.0;
    const std::vector<double>& expected = full.committed_solution();
    const std::vector<double>& actual = restarted.committed_solution();
    if (expected.size() != actual.size())
        passed = check(false, "H20.36 restart preserves the global state layout") && passed;
    else
        for (std::size_t dof = 0; dof < expected.size(); ++dof)
            maximum_restart_difference = std::max(maximum_restart_difference, std::abs(expected[dof] - actual[dof]));
    passed = check(maximum_restart_difference < 1.0e-13, "H20.36 restart reproduces the uninterrupted final state")
             && check(std::remove(checkpoint_path.c_str()) == 0, "H20.36 checkpoint artifact is removed") && passed;
    std::cout << "h20_36_restart_maximum_absolute_difference=" << maximum_restart_difference << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3)
        return 2;
    try {
        fuelsim::PetscSession session(argc, argv, "HEX20 internal friction transaction contract");
        return run_path(argv[1], argv[2]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
