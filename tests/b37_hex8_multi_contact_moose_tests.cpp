#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
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
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> fields{};
    bool present = false;
};

struct ContactReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double pressure;
    double nodal_area;
    double tangential_force;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

std::size_t column(const std::vector<std::string>& header, const std::string& name, const std::string& path) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE CSV is missing column " + name + ": " + path);
    return static_cast<std::size_t>(found - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete B3.7 MOOSE row: " + path);
    const double value = std::stod(values[index]);
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite B3.7 MOOSE value: " + path);
    return value;
}

std::size_t node_id(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    const double value = number(values, index, path);
    if (value < 0.0 || std::floor(value) != value) throw std::invalid_argument("Invalid B3.7 node identifier");
    return static_cast<std::size_t>(value);
}

std::vector<NodeReference> read_nodes(const std::string& path, std::size_t node_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.7 MOOSE nodes: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("Empty B3.7 MOOSE nodal output: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::size_t t = column(header, "T", path), dx = column(header, "disp_x", path),
                      dy = column(header, "disp_y", path), dz = column(header, "disp_z", path),
                      id = column(header, "id", path), x = column(header, "x", path), y = column(header, "y", path),
                      z = column(header, "z", path);
    std::vector<NodeReference> result(node_count);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        const std::size_t source = node_id(values, id, path);
        if (source >= node_count) throw std::invalid_argument("B3.7 MOOSE node identifier is outside the mesh");
        NodeReference candidate{{number(values, x, path), number(values, y, path), number(values, z, path)},
            {number(values, t, path), number(values, dx, path), number(values, dy, path), number(values, dz, path)},
            true};
        if (result[source].present) throw std::invalid_argument("B3.7 MOOSE emitted duplicate node values");
        result[source] = candidate;
    }
    if (std::any_of(result.begin(), result.end(), [](const NodeReference& value) { return !value.present; }))
        throw std::invalid_argument("B3.7 MOOSE nodal output does not contain every source node");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path, const std::string& tangential_name) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.7 MOOSE contact output: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("Empty B3.7 MOOSE contact output: " + path);
    const auto header = split_csv(line);
    const std::size_t pressure = column(header, "contact_pressure", path), id = column(header, "id", path),
                      area = column(header, "nodal_area", path), tangential = column(header, tangential_name, path),
                      x = column(header, "x", path), y = column(header, "y", path), z = column(header, "z", path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        result.push_back(
            {node_id(values, id, path), {number(values, x, path), number(values, y, path), number(values, z, path)},
                number(values, pressure, path), number(values, area, path), number(values, tangential, path)});
    }
    if (result.empty()) throw std::invalid_argument("B3.7 MOOSE contact output has no nodes: " + path);
    return result;
}

std::array<fuelsim::test::FieldErrorMetrics, 4> compare_nodes(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::SteadyProblem& problem, const std::vector<double>& state,
    const std::vector<NodeReference>& reference, double& coordinate_error) {
    if (reference.size() != mesh.nodes().size()) throw std::invalid_argument("B3.7 node counts differ");
    std::array<fuelsim::test::FieldErrorMetrics, 4> result;
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    std::vector<bool> seen(mesh.nodes().size(), false);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (source >= mesh.nodes().size() || seen[source]) throw std::invalid_argument("B3.7 node mapping repeats");
            seen[source] = true;
            const auto& actual = mesh.nodes()[source];
            const auto& expected = reference[source];
            coordinate_error = std::max({coordinate_error, std::abs(actual.x - expected.point.x),
                std::abs(actual.y - expected.point.y), std::abs(actual.z - expected.point.z)});
            result[0].add(state[dofs.dof(fuelsim::Field::temperature, offset + local)], expected.fields[0]);
            result[1].add(state[dofs.dof(fuelsim::Field::displacement_x, offset + local)], expected.fields[1]);
            result[2].add(state[dofs.dof(fuelsim::Field::displacement_y, offset + local)], expected.fields[2]);
            result[3].add(state[dofs.dof(fuelsim::Field::displacement_z, offset + local)], expected.fields[3]);
        }
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end())
        throw std::invalid_argument("B3.7 comparison did not visit every source node");
    return result;
}

bool compare_contact(const fuelsim::SpatialDefinition& spatial, const fuelsim::SteadyProblem& problem,
    const std::vector<double>& state, std::size_t contact, const std::vector<ContactReference>& reference,
    double& coordinate_error, const std::string& name) {
    const auto actual = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, contact, state);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, contact);
    if (actual.size() != source_nodes.size() || reference.size() != source_nodes.size())
        throw std::invalid_argument("B3.7 contact node counts differ for " + name);
    fuelsim::test::FieldErrorMetrics pressure, tangential;
    double actual_force = 0.0, reference_force = 0.0, actual_tangential_force = 0.0, reference_tangential_force = 0.0;
    std::size_t active = 0, sliding = 0;
    for (std::size_t index = 0; index < actual.size(); ++index) {
        const auto found = std::find_if(reference.begin(), reference.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[index]; });
        if (found == reference.end()) throw std::invalid_argument("B3.7 MOOSE contact node is missing for " + name);
        coordinate_error = std::max({coordinate_error, std::abs(actual[index].x - found->point.x),
            std::abs(actual[index].y - found->point.y), std::abs(actual[index].z - found->point.z)});
        pressure.add(actual[index].pressure, found->pressure);
        tangential.add(actual[index].tangential_traction, found->tangential_force / found->nodal_area);
        actual_force += actual[index].contact_force;
        actual_tangential_force += actual[index].tangential_force;
        reference_force += found->pressure * found->nodal_area;
        reference_tangential_force += found->tangential_force;
        if (actual[index].pressure > 0.0) {
            ++active;
            if (actual[index].sliding) ++sliding;
        }
    }
    const auto interface = fuelsim::cartesian::ProblemAccess::summarize_interface(problem, contact, state);
    if (!(std::abs(reference_force) > 0.0) || !(std::abs(reference_tangential_force) > 0.0))
        throw std::invalid_argument("B3.7 contact reference resultants must be nonzero for " + name);
    const double force_error = std::abs(actual_force - reference_force) / std::abs(reference_force);
    const double tangential_force_error =
        std::abs(actual_tangential_force - reference_tangential_force) / std::abs(reference_tangential_force);
    fuelsim::test::print_relative_metrics("b37_" + name + "_pressure", pressure);
    fuelsim::test::print_relative_metrics("b37_" + name + "_tangential_traction", tangential);
    std::cout << "b37_" << name << "_active_contact_nodes=" << active << '\n'
              << "b37_" << name << "_sliding_contact_nodes=" << sliding << '\n'
              << "b37_" << name << "_normal_force_relative_error=" << force_error << '\n'
              << "b37_" << name << "_tangential_force_relative_error=" << tangential_force_error << '\n';
    const double maximum_capacity_excess = [&]() {
        double value = 0.0;
        for (const auto& node : actual)
            value = std::max(
                value, node.tangential_traction - spatial.contacts[contact].friction_coefficient * node.pressure);
        return value;
    }();
    return check(fuelsim::test::relative_metrics_below(pressure, 5.0e-3),
               "B3.7 " + name + " pressure three MOOSE metrics pass") &&
           check(fuelsim::test::relative_metrics_below(tangential, 5.0e-3),
               "B3.7 " + name + " tangential traction three MOOSE metrics pass") &&
           check(active > 0 && interface.active_contact_nodes == active,
               "B3.7 " + name + " has active projected contact nodes") &&
           check(sliding > 0, "B3.7 " + name + " activates Coulomb sliding") &&
           check(force_error < 5.0e-3 && tangential_force_error < 5.0e-3,
               "B3.7 " + name + " normal and tangential resultants pass") &&
           check(maximum_capacity_excess <= 1.0e-12 * interface.maximum_contact_pressure,
               "B3.7 " + name + " respects the Coulomb cap");
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_b37_hex8_multi_contact_moose_tests <case.fsi> <all-nodes.csv> "
                     "<pair-a.csv> <pair-b.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B3.7 three-dimensional multi-contact MOOSE comparison\n");
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(argv[1]);
        if (definition.problem != fuelsim::CaseProblem::steady ||
            definition.geometry != fuelsim::CaseGeometry::cartesian_3d || definition.spatial.regions.size() != 4 ||
            definition.spatial.contacts.size() != 2)
            throw std::invalid_argument("B3.7 requires four three-dimensional regions and two contact pairs");
        for (const auto& contact : definition.spatial.contacts)
            if (!contact.mechanical || contact.thermal || contact.friction_coefficient != 0.001)
                throw std::invalid_argument("B3.7 requires mechanical Coulomb friction on both contact pairs");
        const fuelsim::UnstructuredHex8Mesh source = fuelsim::read_exodus_hex8(definition.mesh_file);
        fuelsim::SteadyProblem problem(definition.spatial, source);
        if (fuelsim::cartesian::ProblemAccess::region_mesh(problem, 0).elements().size() ==
                fuelsim::cartesian::ProblemAccess::region_mesh(problem, 1).elements().size() ||
            fuelsim::cartesian::ProblemAccess::region_mesh(problem, 2).elements().size() ==
                fuelsim::cartesian::ProblemAccess::region_mesh(problem, 3).elements().size())
            throw std::invalid_argument("B3.7 requires nonmatching primary and secondary face partitions");
        const fuelsim::SolverOptions options = {definition.solver.absolute_tolerance,
            definition.solver.relative_tolerance, definition.solver.step_tolerance,
            definition.solver.maximum_iterations};
        const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem,
            {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
                definition.steady_execution.maximum_cutbacks_per_step,
                definition.steady_execution.minimum_load_increment},
            options);
        bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
            "B3.7 completes the two-pair three-dimensional load step");
        double coordinate_error = 0.0;
        const auto fields = compare_nodes(
            source, problem, solve.solve.state, read_nodes(argv[2], source.nodes().size()), coordinate_error);
        const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            if (field == 3 || fields[field].maximum_reference < 1.0e-8) {
                fuelsim::test::print_absolute_metrics("b37_" + names[field], fields[field]);
                passed = check(fuelsim::test::absolute_metrics_below(fields[field], 1.0e-8),
                             "B3.7 " + names[field] + " physical-zero absolute gate passes") &&
                         passed;
            } else if (field == 2) {
                fuelsim::test::print_relative_metrics("b37_" + names[field], fields[field]);
                passed =
                    check(fuelsim::test::relative_metrics_below_with_pointwise_tolerance(fields[field], 5.0e-3, 5.0e-2),
                        "B3.7 displacement_y aggregate metrics pass with the near-zero pointwise gate") &&
                    passed;
            } else if (fields[field].has_relative_norm()) {
                fuelsim::test::print_relative_metrics("b37_" + names[field], fields[field]);
                passed = check(fuelsim::test::relative_metrics_below(fields[field], 5.0e-3),
                             "B3.7 " + names[field] + " three MOOSE metrics pass") &&
                         passed;
            } else {
                fuelsim::test::print_absolute_metrics("b37_" + names[field], fields[field]);
                passed = check(fuelsim::test::absolute_metrics_below(fields[field], 1.0e-12),
                             "B3.7 " + names[field] + " zero-reference absolute gate passes") &&
                         passed;
            }
        }
        passed = check(coordinate_error < 1.0e-12, "B3.7 compares MOOSE values at matching coordinates") && passed;
        passed = compare_contact(definition.spatial, problem, solve.solve.state, 0,
                     read_contact(argv[3], "pair_a_tangential_force_y"), coordinate_error, "pair_a") &&
                 passed;
        passed = compare_contact(definition.spatial, problem, solve.solve.state, 1,
                     read_contact(argv[4], "pair_b_tangential_force_y"), coordinate_error, "pair_b") &&
                 passed;
        if (!passed) return 1;
        std::cout << "[PASS] B3.7 three-dimensional multi-contact friction MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B3.7 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
