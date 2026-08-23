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
    double slip_1, slip_2, opening, pressure;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) result.push_back(field);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("H20.36 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("H20.36 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("H20.36 CSV contains an invalid integer index: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<std::vector<DisplacementReference>> read_displacement_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.36 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,time,id,x,y,z,disp_x,disp_y,disp_z")
        throw std::invalid_argument("Unexpected H20.36 displacement header in " + path);
    std::vector<std::vector<DisplacementReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.36 displacement step is invalid: " + path);
        result[step - 1].push_back(
            {index_value(values, 2, path), {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 6, path), number(values, 7, path), number(values, 8, path)}});
    }
    return result;
}

std::vector<std::vector<ContactReference>> read_contact_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.36 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line !=
        "step,time,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected H20.36 contact header in " + path);
    std::vector<std::vector<ContactReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.36 contact step is invalid: " + path);
        result[step - 1].push_back(
            {index_value(values, 2, path), {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
                {number(values, 9, path), number(values, 10, path), number(values, 11, path)}, number(values, 12, path),
                number(values, 13, path), number(values, 14, path), number(values, 15, path)});
    }
    return result;
}

bool check(bool condition, const std::string& message) {
    if (condition) return true;
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
    return {first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

std::array<double, 3> unit(const std::array<double, 3>& value) {
    const double measure = std::sqrt(dot(value, value));
    if (!(measure > 0.0)) throw std::invalid_argument("H20.36 contact face is degenerate");
    return {value[0] / measure, value[1] / measure, value[2] / measure};
}

double maximum_secondary_face_nonplanarity(const fuelsim::UnstructuredHex20Mesh& mesh) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {
        {{{0, 1, 5, 4, 8, 13, 16, 12}}, {{1, 2, 6, 5, 9, 14, 17, 13}}, {{2, 3, 7, 6, 10, 15, 18, 14}},
            {{3, 0, 4, 7, 11, 12, 19, 15}}, {{0, 3, 2, 1, 11, 10, 9, 8}}, {{4, 5, 6, 7, 16, 17, 18, 19}}}};
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

fuelsim::TransientTimeOptions time_options(
    const fuelsim::FuelSimCaseDefinition& definition, double end_time, double initial_time_step) {
    fuelsim::TransientTimeOptions options = definition.transient_execution;
    options.end_time = end_time;
    options.initial_time_step = initial_time_step;
    return options;
}

struct ContactNumericalEvidence final {
    std::array<double, 3> residual_balance{};
    double maximum_directional_error = 0.0;
};

ContactNumericalEvidence inspect_contact_transition(fuelsim::TransientProblem& problem,
    fuelsim::TransientCommittedState previous, double target_time, const std::vector<double>& state) {
    fuelsim::cartesian::ProblemAccess::restore_committed_state(problem, std::move(previous));
    problem.begin_time_step({target_time, 1.0, false});
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    ContactNumericalEvidence result;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
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
        if (!(reference_squared > 0.0)) throw std::logic_error("H20.36 contact Jacobian direction is zero");
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
        states.push_back({step.time, solution, contact,
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0),
            fuelsim::cartesian::ProblemAccess::committed_state(problem)});
    }

    std::vector<StepState> states;
};

bool compare_abaqus_steps(const fuelsim::TransientProblem& problem, const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::vector<StepState>& states, const std::string& displacement_path, const std::string& contact_path) {
    const auto displacement = read_displacement_reference(displacement_path);
    const auto contact = read_contact_reference(contact_path);
    if (states.size() != displacement.size() || states.size() != contact.size())
        throw std::invalid_argument("H20.36 Fuelsim and Abaqus step counts differ");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    const std::vector<std::size_t> contact_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    constexpr double tolerance = 1.0e-2;
    constexpr double radial_displacement_pointwise_tolerance = 1.4e-2;
    constexpr double axial_force_pointwise_tolerance = 2.0e-2;
    constexpr double axial_slip_pointwise_tolerance = 6.0e-2;
    constexpr double zero_tolerance = 1.2e-8;
    bool passed = true;
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_metrics;
        fuelsim::test::FieldErrorMetrics normal_force, circumferential_force, axial_force, slip_1, slip_2;
        std::vector<bool> present(mesh.nodes().size(), false);
        double maximum_coordinate_difference = 0.0;
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& region_mesh = spatial.hex20_region_mesh(region);
            for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
                const std::size_t source = region_mesh.source_node_ids()[local];
                const auto found = std::find_if(displacement[step].begin(), displacement[step].end(),
                    [source](const DisplacementReference& value) { return value.id == source; });
                if (found == displacement[step].end() || present[source])
                    throw std::invalid_argument("H20.36 displacement source-node mapping is incomplete or repeated");
                present[source] = true;
                maximum_coordinate_difference =
                    std::max({maximum_coordinate_difference, std::abs(mesh.nodes()[source].x - found->point.x),
                        std::abs(mesh.nodes()[source].y - found->point.y),
                        std::abs(mesh.nodes()[source].z - found->point.z)});
                const std::size_t global = spatial.global_node(region, local);
                const double radius = std::hypot(found->point.x, found->point.y);
                const double actual_x = states[step].solution[fields[1].begin + global];
                const double actual_y = states[step].solution[fields[2].begin + global];
                const double actual_z = states[step].solution[fields[3].begin + global];
                displacement_metrics[0].add((found->point.x * actual_x + found->point.y * actual_y) / radius,
                    (found->point.x * found->displacement[0] + found->point.y * found->displacement[1]) / radius);
                displacement_metrics[1].add((-found->point.y * actual_x + found->point.x * actual_y) / radius,
                    (-found->point.y * found->displacement[0] + found->point.x * found->displacement[1]) / radius);
                displacement_metrics[2].add(actual_z, found->displacement[2]);
            }
        }
        if (contact_nodes.size() != states[step].contact.size() || contact_nodes.size() != contact[step].size())
            throw std::invalid_argument("H20.36 contact-node counts differ");
        double actual_normal_resultant = 0.0, reference_normal_resultant = 0.0;
        double actual_circumferential_resultant = 0.0, reference_circumferential_resultant = 0.0;
        double actual_axial_resultant = 0.0, reference_axial_resultant = 0.0;
        bool abaqus_contact_closed = true;
        for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
            const auto found = std::find_if(contact[step].begin(), contact[step].end(),
                [&](const ContactReference& value) { return value.id == contact_nodes[node]; });
            if (found == contact[step].end()) throw std::invalid_argument("H20.36 contact-node mapping is incomplete");
            const fuelsim::CartesianContactNodeSummary& actual = states[step].contact[node];
            maximum_coordinate_difference =
                std::max({maximum_coordinate_difference, std::abs(actual.x - found->point.x),
                    std::abs(actual.y - found->point.y), std::abs(actual.z - found->point.z)});
            const double radius = std::hypot(found->point.x, found->point.y);
            const std::array<double, 3> radial = {found->point.x / radius, found->point.y / radius, 0.0};
            const std::array<double, 3> circumferential = {-found->point.y / radius, found->point.x / radius, 0.0};
            const std::array<double, 3> actual_normal = {
                -actual.normal_contact_force[0], -actual.normal_contact_force[1], -actual.normal_contact_force[2]};
            const std::array<double, 3> actual_tangent = {-actual.tangential_contact_force[0],
                -actual.tangential_contact_force[1], -actual.tangential_contact_force[2]};
            const double actual_normal_value = dot(actual_normal, radial);
            const double reference_normal_value = dot(found->normal_force, radial);
            const double actual_circumferential = dot(actual_tangent, circumferential);
            const double reference_circumferential = dot(found->tangential_force, circumferential);
            normal_force.add(actual_normal_value, reference_normal_value);
            circumferential_force.add(actual_circumferential, reference_circumferential);
            axial_force.add(actual_tangent[2], found->tangential_force[2]);
            // For these C3D20 S6 faces, Abaqus local direction 1 opposes the
            // positive circumferential direction and direction 2 opposes the
            // positive global axial direction.
            slip_1.add(-dot(actual.tangential_slip, circumferential), found->slip_1);
            slip_2.add(-actual.tangential_slip[2], found->slip_2);
            actual_normal_resultant += actual_normal_value;
            reference_normal_resultant += reference_normal_value;
            actual_circumferential_resultant += actual_circumferential;
            reference_circumferential_resultant += reference_circumferential;
            actual_axial_resultant += actual_tangent[2];
            reference_axial_resultant += found->tangential_force[2];
            abaqus_contact_closed = abaqus_contact_closed && found->pressure > 0.0 && found->opening < 0.0;
        }
        const std::string prefix = "h20_36_step_" + std::to_string(step + 1) + "_";
        const double circumferential_displacement_difference_on_axial_scale =
            displacement_metrics[1].maximum_absolute_difference / displacement_metrics[2].maximum_reference;
        const double circumferential_force_difference_on_axial_scale =
            circumferential_force.maximum_absolute_difference / axial_force.maximum_reference;
        const double circumferential_slip_difference_on_axial_scale =
            slip_1.maximum_absolute_difference / slip_2.maximum_reference;
        fuelsim::test::print_relative_metrics(prefix + "displacement_radial", displacement_metrics[0]);
        if (displacement_metrics[1].has_relative_norm())
            fuelsim::test::print_relative_metrics(prefix + "displacement_circumferential", displacement_metrics[1]);
        fuelsim::test::print_absolute_metrics(prefix + "displacement_circumferential", displacement_metrics[1]);
        std::cout << prefix << "displacement_circumferential_maximum_difference_on_axial_scale="
                  << circumferential_displacement_difference_on_axial_scale << '\n';
        fuelsim::test::print_relative_metrics(prefix + "displacement_axial", displacement_metrics[2]);
        passed = check(fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                           displacement_metrics[0], tolerance, radial_displacement_pointwise_tolerance) &&
                           displacement_metrics[0].maximum_zero_reference_difference < zero_tolerance &&
                           fuelsim::test::relative_metrics_below(displacement_metrics[2], tolerance) &&
                           displacement_metrics[2].maximum_zero_reference_difference < zero_tolerance,
                     "H20.36 step " + std::to_string(step + 1) +
                         "driven displacement fields pass 1 percent aggregate metrics, the explicit 1.4 percent "
                         "radial pointwise gate, and their separate zero-reference checks") &&
                 check(circumferential_displacement_difference_on_axial_scale < tolerance,
                     "H20.36 step " + std::to_string(step + 1) +
                         "undriven circumferential displacement difference is below 1 percent of the driven "
                         "axial-displacement peak") &&
                 passed;
        fuelsim::test::print_relative_metrics(prefix + "signed_normal_radial_force", normal_force);
        fuelsim::test::print_relative_metrics(
            prefix + "signed_tangential_circumferential_force", circumferential_force);
        fuelsim::test::print_relative_metrics(prefix + "signed_tangential_axial_force", axial_force);
        if (slip_1.has_relative_norm()) fuelsim::test::print_relative_metrics(prefix + "tangential_slip_1", slip_1);
        fuelsim::test::print_absolute_metrics(
            prefix + "signed_tangential_circumferential_force", circumferential_force);
        fuelsim::test::print_absolute_metrics(prefix + "tangential_slip_1", slip_1);
        std::cout << prefix << "signed_tangential_circumferential_force_maximum_difference_on_axial_scale="
                  << circumferential_force_difference_on_axial_scale << '\n'
                  << prefix << "tangential_slip_1_maximum_difference_on_axial_scale="
                  << circumferential_slip_difference_on_axial_scale << '\n';
        fuelsim::test::print_relative_metrics(prefix + "tangential_slip_2", slip_2);
        passed = check(maximum_coordinate_difference < 1.0e-12,
                     "H20.36 step " + std::to_string(step + 1) + " uses the tracked quadratic Exodus mesh") &&
                 check(abaqus_contact_closed,
                     "H20.36 step " + std::to_string(step + 1) + " keeps all Abaqus contact nodes closed") &&
                 check(fuelsim::test::relative_metrics_below(normal_force, tolerance) &&
                           fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                               axial_force, tolerance, axial_force_pointwise_tolerance),
                     "H20.36 step " + std::to_string(step + 1) +
                         "signed normal and driven axial nodal forces pass 1 percent aggregate metrics and the "
                         "explicit 2 percent small-value axial-force pointwise gate") &&
                 check(circumferential_force_difference_on_axial_scale < tolerance,
                     "H20.36 step " + std::to_string(step + 1) +
                         "undriven circumferential nodal-force difference is below 1 percent of the driven "
                         "axial-force peak") &&
                 check(fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                           slip_2, tolerance, axial_slip_pointwise_tolerance),
                     "H20.36 step " + std::to_string(step + 1) +
                         "driven axial tangential slip passes 1 percent aggregate metrics and the explicit 6 percent "
                         "near-zero pointwise gate") &&
                 check(circumferential_slip_difference_on_axial_scale < tolerance,
                     "H20.36 step " + std::to_string(step + 1) +
                         "undriven circumferential slip difference is below 1 percent of the driven axial-slip "
                         "peak") &&
                 passed;
        const double normal_resultant_error =
            std::abs(actual_normal_resultant - reference_normal_resultant) / std::abs(reference_normal_resultant);
        const double circumferential_resultant_difference =
            std::abs(actual_circumferential_resultant - reference_circumferential_resultant);
        const double circumferential_resultant_error =
            circumferential_resultant_difference / std::abs(reference_circumferential_resultant);
        const double axial_resultant_error =
            std::abs(actual_axial_resultant - reference_axial_resultant) / std::abs(reference_axial_resultant);
        const double circumferential_resultant_difference_on_axial_scale =
            circumferential_resultant_difference / std::abs(reference_axial_resultant);
        std::cout << prefix << "normal_resultant_relative_error=" << normal_resultant_error << '\n'
                  << prefix << "circumferential_resultant_relative_error=" << circumferential_resultant_error << '\n'
                  << prefix << "circumferential_resultant_absolute_difference=" << circumferential_resultant_difference
                  << '\n'
                  << prefix << "circumferential_resultant_difference_on_axial_scale="
                  << circumferential_resultant_difference_on_axial_scale << '\n'
                  << prefix << "axial_resultant_relative_error=" << axial_resultant_error << '\n';
        passed = check(normal_resultant_error < tolerance && axial_resultant_error < tolerance,
                     "H20.36 step " + std::to_string(step + 1) +
                         "signed normal and driven axial resultants agree with Abaqus below 1 percent") &&
                 check(circumferential_resultant_difference_on_axial_scale < tolerance,
                     "H20.36 step " + std::to_string(step + 1) +
                         "undriven circumferential resultant difference is below 1 percent of the driven axial "
                         "resultant") &&
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
            actual[point].cartesian_elastic_tangential_slip != expected[point].cartesian_elastic_tangential_slip ||
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

std::array<double, 2> tangential_resultant(const StepState& state) {
    std::array<double, 2> result{};
    for (const fuelsim::CartesianContactNodeSummary& node : state.contact) {
        const double radius = std::hypot(node.x, node.y);
        result[0] += (-node.y * node.tangential_contact_force[0] + node.x * node.tangential_contact_force[1]) / radius;
        result[1] += node.tangential_contact_force[2];
    }
    return result;
}

bool run_path(const std::string& input_path, const std::string& displacement_path, const std::string& contact_path,
    const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    Recorder full_recorder;
    fuelsim::TransientProblem full(input.spatial, mesh);
    const fuelsim::TransientCommittedState initial_state = fuelsim::cartesian::ProblemAccess::committed_state(full);
    const fuelsim::TransientResult full_result = fuelsim::solve_transient(
        full, time_options(input, input.transient_execution.end_time, 1.0), solver, &full_recorder);
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
        passed = check(std::all_of(full_recorder.states.begin(), full_recorder.states.end(),
                           [](const StepState& state) {
                               const auto count = stick_slide_counts(state);
                               return count[0] + count[1] == 37;
                           }),
                     "H20.36 keeps all 37 curved constraints closed throughout the path") &&
                 check(first[0] == 37 && first[1] == 0 && second[0] == 37 && second[1] == 0,
                     "H20.36 starts with two fully sticking states") &&
                 check(mixed[0] > 0 && mixed[1] > 0 && forward[1] > 0,
                     "H20.36 contains simultaneous sticking and sliding before forward sliding") &&
                 check(unload[0] > 0 && unload[1] > 0 && reverse[0] > 0 && reverse[1] > 0 &&
                           forward_resultant[1] * reverse_resultant[1] < 0.0,
                     "H20.36 covers unloading and reverses the driven axial tangential resultant") &&
                 check(restick[0] == 37 && restick[1] == 0, "H20.36 finishes with all 37 constraints restuck") &&
                 check(maximum_secondary_face_nonplanarity(mesh) > 1.0e-3,
                     "H20.36 uses genuinely quadratic contact faces rather than planar facets") &&
                 passed;
    }
    passed = compare_abaqus_steps(full, mesh, full_recorder.states, displacement_path, contact_path) && passed;

    if (full_recorder.states.size() == 7) {
        fuelsim::TransientProblem sticking_probe(input.spatial, mesh);
        const ContactNumericalEvidence sticking =
            inspect_contact_transition(sticking_probe, initial_state, 1.0, full_recorder.states[0].solution);
        fuelsim::TransientProblem sliding_probe(input.spatial, mesh);
        const ContactNumericalEvidence sliding = inspect_contact_transition(
            sliding_probe, full_recorder.states[4].committed, 6.0, full_recorder.states[5].solution);
        std::cout << "h20_36_sticking_contact_jacobian_directional_error=" << sticking.maximum_directional_error << '\n'
                  << "h20_36_sliding_contact_jacobian_directional_error=" << sliding.maximum_directional_error << '\n'
                  << "h20_36_sticking_contact_residual_balance=" << sticking.residual_balance[0] << ','
                  << sticking.residual_balance[1] << ',' << sticking.residual_balance[2] << '\n'
                  << "h20_36_sliding_contact_residual_balance=" << sliding.residual_balance[0] << ','
                  << sliding.residual_balance[1] << ',' << sliding.residual_balance[2] << '\n';
        passed = check(sticking.maximum_directional_error < 1.0e-7 && sliding.maximum_directional_error < 1.0e-7,
                     "H20.36 curved sticking and sliding contact Jacobians match centered directional differences") &&
                 check(std::all_of(sticking.residual_balance.begin(), sticking.residual_balance.end(),
                           [](double value) { return std::abs(value) < 1.0e-8; }) &&
                           std::all_of(sliding.residual_balance.begin(), sliding.residual_balance.end(),
                               [](double value) { return std::abs(value) < 1.0e-8; }),
                     "H20.36 curved sticking and sliding residuals are action-reaction conservative") &&
                 passed;
    }

    Recorder first_recorder;
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(input, 4.0, 1.0), solver, &first_recorder);
    passed = check(first.completed && first_recorder.states.size() == 4,
                 "H20.36 reaches its nonzero-slip checkpoint state") &&
             check(recorded_states_identical(first_recorder.states, full_recorder.states, 0),
                 "H20.36 split run reproduces the first four uninterrupted nodal and contact-history states") &&
             passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step);
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    Recorder second_recorder;
    const fuelsim::TransientResult second = fuelsim::solve_transient(
        restarted, time_options(input, input.transient_execution.end_time, restored_step), solver, &second_recorder);
    passed = check(second.completed && second_recorder.states.size() == 3,
                 "H20.36 checkpoint restores and completes the reversed path") &&
             check(recorded_states_identical(second_recorder.states, full_recorder.states, 4),
                 "H20.36 restart reproduces every remaining nodal and contact-history state exactly") &&
             passed;
    double maximum_restart_difference = 0.0;
    const std::vector<double>& expected = full.committed_solution();
    const std::vector<double>& actual = restarted.committed_solution();
    if (expected.size() != actual.size())
        passed = check(false, "H20.36 restart preserves the global state layout") && passed;
    else
        for (std::size_t dof = 0; dof < expected.size(); ++dof)
            maximum_restart_difference = std::max(maximum_restart_difference, std::abs(expected[dof] - actual[dof]));
    passed = check(maximum_restart_difference < 1.0e-13, "H20.36 restart reproduces the uninterrupted final state") &&
             check(std::remove(checkpoint_path.c_str()) == 0, "H20.36 checkpoint artifact is removed") && passed;
    std::cout << "h20_36_restart_maximum_absolute_difference=" << maximum_restart_difference << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_h20_curved_friction_path_abaqus_tests <case.fsi> <displacement.csv> "
                     "<contact.csv> <checkpoint>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.36 curved friction-path Abaqus comparison\n");
        if (!run_path(argv[1], argv[2], argv[3], argv[4])) return 1;
        std::cout << "[PASS] fuelsim H20.36 curved friction path and restart\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.36 curved friction path raised: " << error.what() << '\n';
        return 1;
    }
}
