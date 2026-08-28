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
struct NodeReference final {
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

struct ReactionReference final {
    std::array<double, 3> reaction;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("B3.4 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("B3.4 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("B3.4 CSV contains an invalid integer: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.4 Abaqus nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "id,x,y,z,u1,u2,u3")
        throw std::invalid_argument("Unexpected B3.4 nodal header: " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)}});
    }
    if (result.size() != 16) throw std::invalid_argument("B3.4 Abaqus reference must contain sixteen nodes");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.4 Abaqus contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected B3.4 contact header: " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
                {number(values, 7, path), number(values, 8, path), number(values, 9, path)}, number(values, 10, path),
                number(values, 11, path), number(values, 12, path), number(values, 13, path)});
    }
    if (result.size() != 4) throw std::invalid_argument("B3.4 Abaqus reference must contain four contact nodes");
    return result;
}

ReactionReference read_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.4 Abaqus reaction reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "reaction_x,reaction_y,reaction_z" || !std::getline(input, line))
        throw std::invalid_argument("Unexpected B3.4 reaction reference: " + path);
    const std::vector<std::string> values = split_csv(line);
    return {{{number(values, 0, path), number(values, 1, path), number(values, 2, path)}}};
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& input) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = input.solver.absolute_tolerance;
    result.relative_tolerance = input.solver.relative_tolerance;
    result.step_tolerance = input.solver.step_tolerance;
    result.maximum_iterations = input.solver.maximum_iterations;
    result.temperature_residual_scale = input.solver.temperature_residual_scale;
    result.mechanical_residual_scale = input.solver.mechanical_residual_scale;
    return result;
}

std::array<fuelsim::test::FieldErrorMetrics, 3> compare_nodes(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::SteadyProblem& problem, const std::vector<double>& state,
    const std::vector<NodeReference>& reference, double& maximum_coordinate_difference) {
    std::array<fuelsim::test::FieldErrorMetrics, 3> result;
    std::vector<bool> present(mesh.nodes().size(), false);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto found = std::find_if(
                reference.begin(), reference.end(), [&](const NodeReference& value) { return value.id == source; });
            if (found == reference.end() || present[source])
                throw std::invalid_argument("B3.4 Abaqus source-node mapping is not unique");
            present[source] = true;
            maximum_coordinate_difference =
                std::max(maximum_coordinate_difference, std::max({std::abs(mesh.nodes()[source].x - found->point.x),
                                                            std::abs(mesh.nodes()[source].y - found->point.y),
                                                            std::abs(mesh.nodes()[source].z - found->point.z)}));
            result[0].add(state[dofs.dof(fuelsim::Field::displacement_x, offset + local)], found->displacement[0]);
            result[1].add(state[dofs.dof(fuelsim::Field::displacement_y, offset + local)], found->displacement[1]);
            result[2].add(state[dofs.dof(fuelsim::Field::displacement_z, offset + local)], found->displacement[2]);
        }
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("B3.4 comparison did not visit every node");
    return result;
}

bool histories_equal(
    const std::vector<fuelsim::ContactPointHistory>& first, const std::vector<fuelsim::ContactPointHistory>& second) {
    if (first.size() != second.size()) return false;
    for (std::size_t point = 0; point < first.size(); ++point)
        if (first[point].elastic_tangential_slip != second[point].elastic_tangential_slip ||
            first[point].sliding != second[point].sliding ||
            first[point].normal_multiplier != second[point].normal_multiplier ||
            first[point].cartesian_elastic_tangential_slip != second[point].cartesian_elastic_tangential_slip ||
            first[point].cartesian_total_tangential_slip != second[point].cartesian_total_tangential_slip ||
            first[point].cartesian_tangent_basis_initialized != second[point].cartesian_tangent_basis_initialized ||
            first[point].cartesian_contact_normal != second[point].cartesian_contact_normal ||
            first[point].cartesian_contact_tangent_first != second[point].cartesian_contact_tangent_first)
            return false;
    return true;
}

bool compare(const fuelsim::FuelSimCaseDefinition& input, const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::SteadyProblem& problem, const fuelsim::SteadyResult& solve, const std::string& nodes_path,
    const std::string& contact_path, const std::string& reaction_path) {
    double maximum_coordinate_difference = 0.0;
    const auto displacement =
        compare_nodes(mesh, problem, solve.solve.state, read_nodes(nodes_path), maximum_coordinate_difference);
    const std::vector<ContactReference> reference = read_contact(contact_path);
    const std::vector<fuelsim::CartesianContactNodeSummary> actual =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const std::vector<std::size_t> source_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    if (actual.size() != source_nodes.size()) throw std::invalid_argument("B3.4 contact-node counts differ");

    fuelsim::test::FieldErrorMetrics normal_force, tangential_force_y, slip_y, opening, pressure;
    double maximum_theoretical_zero_force = 0.0;
    std::size_t active = 0, sliding = 0;
    for (std::size_t node = 0; node < actual.size(); ++node) {
        const auto found = std::find_if(reference.begin(), reference.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[node]; });
        if (found == reference.end()) throw std::invalid_argument("B3.4 Abaqus contact-node mapping is incomplete");
        maximum_coordinate_difference = std::max(
            maximum_coordinate_difference, std::max({std::abs(mesh.nodes()[source_nodes[node]].x - found->point.x),
                                               std::abs(mesh.nodes()[source_nodes[node]].y - found->point.y),
                                               std::abs(mesh.nodes()[source_nodes[node]].z - found->point.z)}));
        normal_force.add(-actual[node].normal_contact_force[0], found->normal_force[0]);
        tangential_force_y.add(-actual[node].tangential_contact_force[1], found->tangential_force[1]);
        slip_y.add(actual[node].tangential_slip[1], -found->slip_2);
        opening.add(actual[node].gap, found->opening);
        pressure.add(actual[node].pressure, found->pressure);
        maximum_theoretical_zero_force = std::max(maximum_theoretical_zero_force,
            std::max({std::abs(actual[node].normal_contact_force[1]), std::abs(actual[node].normal_contact_force[2]),
                std::abs(actual[node].tangential_contact_force[0]), std::abs(actual[node].tangential_contact_force[2]),
                std::abs(found->normal_force[1]), std::abs(found->normal_force[2]),
                std::abs(found->tangential_force[0]), std::abs(found->tangential_force[2])}));
        if (actual[node].pressure > 0.0) ++active;
        if (actual[node].sliding) ++sliding;
    }

    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    const ReactionReference reaction = read_reaction(reaction_path);
    const double normal_resultant_error =
        std::abs(interface.total_contact_force - std::abs(reaction.reaction[0])) / std::abs(reaction.reaction[0]);
    const double tangential_resultant_error =
        std::abs(interface.total_tangential_force - std::abs(reaction.reaction[1])) / std::abs(reaction.reaction[1]);

    fuelsim::test::print_relative_metrics("b34_displacement_x", displacement[0]);
    fuelsim::test::print_relative_metrics("b34_displacement_y", displacement[1]);
    fuelsim::test::print_absolute_metrics("b34_displacement_z", displacement[2]);
    fuelsim::test::print_relative_metrics("b34_signed_normal_force_x", normal_force);
    fuelsim::test::print_relative_metrics("b34_signed_tangential_force_y", tangential_force_y);
    fuelsim::test::print_relative_metrics("b34_tangential_slip_y", slip_y);
    fuelsim::test::print_relative_metrics("b34_opening", opening);
    fuelsim::test::print_relative_metrics("b34_pressure", pressure);
    std::cout << "b34_displacement_z_zero_reference_count=" << displacement[2].zero_reference_count << '\n'
              << "b34_displacement_z_maximum_absolute_difference=" << displacement[2].maximum_absolute_difference
              << '\n'
              << "b34_theoretical_zero_contact_force_component_value_count=32\n"
              << "b34_theoretical_zero_contact_force_maximum_absolute_difference=" << maximum_theoretical_zero_force
              << '\n'
              << "b34_normal_resultant_relative_error=" << normal_resultant_error << '\n'
              << "b34_tangential_resultant_relative_error=" << tangential_resultant_error << '\n'
              << "b34_maximum_mesh_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << "b34_active_contact_nodes=" << active << '\n'
              << "b34_sliding_contact_nodes=" << sliding << '\n';

    constexpr double tolerance = 5.0e-3;
    return check(active == 4 && sliding > 0 && interface.active_contact_nodes == 4,
               "B3.4 keeps four active constraints and reaches the Coulomb sliding branch") &&
           check(fuelsim::test::relative_metrics_below(displacement[0], tolerance) &&
                     fuelsim::test::relative_metrics_below(displacement[1], tolerance),
               "B3.4 in-plane nodal displacements pass all three Abaqus metrics below 0.5 percent") &&
           check(displacement[2].maximum_absolute_difference < 5.0e-10,
               "B3.4 out-of-plane symmetry displacement agrees with Abaqus in absolute value") &&
           check(fuelsim::test::relative_metrics_below(normal_force, tolerance) &&
                     fuelsim::test::relative_metrics_below(tangential_force_y, tolerance),
               "B3.4 signed in-plane nodal contact-force fields pass all three Abaqus metrics below 0.5 percent") &&
           check(fuelsim::test::relative_metrics_below(slip_y, tolerance) &&
                     fuelsim::test::relative_metrics_below(opening, tolerance) &&
                     fuelsim::test::relative_metrics_below(pressure, tolerance),
               "B3.4 slip, opening, and pressure pass all three Abaqus metrics below 0.5 percent") &&
           check(maximum_theoretical_zero_force < 1.0e-4,
               "B3.4 theoretical-zero transverse contact-force components pass their absolute check") &&
           check(normal_resultant_error < tolerance && tangential_resultant_error < tolerance,
               "B3.4 normal and tangential resultants agree with Abaqus below 0.5 percent") &&
           check(maximum_coordinate_difference < 5.0e-10,
               "B3.4 Abaqus and Fuelsim use matching tracked HEX8 coordinates") &&
           check(std::abs(interface.total_tangential_force -
                          input.spatial.contacts[0].friction_coefficient * interface.total_contact_force) <
                     1.0e-12 * interface.total_contact_force,
               "B3.4 Fuelsim sliding resultant lies on the Coulomb cap");
}

bool restart(const fuelsim::FuelSimCaseDefinition& input, const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::string& checkpoint_path) {
    fuelsim::TransientProblem full(input.spatial, mesh);
    const fuelsim::TransientResult solve =
        fuelsim::solve_transient(full, {1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0, 0.0}, solver_options(input));
    bool passed = check(solve.completed && solve.accepted_steps.size() == 1,
        "B3.4 transient sliding state commits one physical time step");
    fuelsim::write_transient_checkpoint(checkpoint_path, full, 0.25);
    fuelsim::TransientProblem restored(input.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    const bool histories_match =
        histories_equal(fuelsim::cartesian::ProblemAccess::committed_contact_histories(restored).at(0),
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(full).at(0));
    passed =
        check(restored_step == 0.25 && restored.committed_solution() == full.committed_solution() && histories_match,
            "B3.4 checkpoint restores the exact three-dimensional friction transaction") &&
        passed;
    passed = check(std::remove(checkpoint_path.c_str()) == 0, "B3.4 removes its checkpoint artifact") && passed;
    return passed;
}

bool run(const std::string& input_path, const std::string& nodes_path, const std::string& contact_path,
    const std::string& reaction_path, const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    if (input.problem != fuelsim::CaseProblem::steady || input.geometry != fuelsim::CaseGeometry::cartesian_3d ||
        input.spatial.contacts.size() != 1 || input.spatial.contacts[0].thermal ||
        !input.spatial.contacts[0].mechanical || input.spatial.contacts[0].friction_coefficient != 0.001 ||
        input.spatial.contacts[0].mechanical_discretization !=
            fuelsim::MechanicalContactDiscretization::surface_to_surface ||
        input.spatial.contacts[0].mechanical_sliding != fuelsim::MechanicalContactSliding::small ||
        input.spatial.contacts[0].friction_slip_tolerance != 1.0e-8)
        throw std::invalid_argument(
            "B3.4 requires isolated small-sliding surface contact with explicit Abaqus slip tolerance");
    const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(input.mesh_file);
    fuelsim::SteadyProblem problem(input.spatial, mesh);
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem,
        {input.steady_execution.load_steps, input.steady_execution.cutback_factor,
            input.steady_execution.maximum_cutbacks_per_step, input.steady_execution.minimum_load_increment},
        solver_options(input));
    bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
        "B3.4 sliding-friction path completes one full load step");
    passed = compare(input, mesh, problem, solve, nodes_path, contact_path, reaction_path) && passed;
    passed = restart(input, mesh, checkpoint_path) && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: fuelsim_b34_hex8_sliding_contact_abaqus_tests "
                     "<case.fsi> <nodes.csv> <contact.csv> <reaction.csv> <checkpoint.bin>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B3.4 HEX8 sliding-contact Abaqus comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3], argv[4], argv[5]);
        if (passed && session.rank() == 0)
            std::cout << "[PASS] B3.4 HEX8 sliding-contact Abaqus comparison and restart\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B3.4 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
