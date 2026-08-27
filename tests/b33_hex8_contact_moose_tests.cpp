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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> fields{};
};

struct ContactReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double pressure;
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

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete three-dimensional contact row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_nodes(const std::string& thermal_path, const std::string& mechanical_path) {
    std::ifstream thermal(thermal_path);
    if (!thermal) throw std::runtime_error("Could not read three-dimensional thermal-contact nodes: " + thermal_path);
    std::string line;
    std::getline(thermal, line);
    if (line != "T,id,x,y,z") throw std::invalid_argument("Unexpected thermal-contact header in " + thermal_path);
    std::vector<NodeReference> result;
    while (std::getline(thermal, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, thermal_path)),
            {number(values, 2, thermal_path), number(values, 3, thermal_path), number(values, 4, thermal_path)},
            {number(values, 0, thermal_path), 0.0, 0.0, 0.0}});
    }
    std::ifstream mechanical(mechanical_path);
    if (!mechanical)
        throw std::runtime_error("Could not read three-dimensional friction-contact nodes: " + mechanical_path);
    std::getline(mechanical, line);
    if (line != "disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected friction-contact header in " + mechanical_path);
    std::size_t row = 0;
    while (std::getline(mechanical, line)) {
        const auto values = split_csv(line);
        if (row >= result.size() || result[row].id != static_cast<std::size_t>(number(values, 3, mechanical_path)))
            throw std::invalid_argument("Thermal and mechanical MOOSE node order differs");
        const fuelsim::CartesianPoint3 point = {
            number(values, 4, mechanical_path), number(values, 5, mechanical_path), number(values, 6, mechanical_path)};
        if (std::max({std::abs(point.x - result[row].point.x), std::abs(point.y - result[row].point.y),
                std::abs(point.z - result[row].point.z)}) > 1.0e-12)
            throw std::invalid_argument("Thermal and mechanical MOOSE node coordinates differ");
        result[row].fields[1] = number(values, 0, mechanical_path);
        result[row].fields[2] = number(values, 1, mechanical_path);
        result[row].fields[3] = number(values, 2, mechanical_path);
        ++row;
    }
    if (row != result.size()) throw std::invalid_argument("Thermal and mechanical MOOSE node counts differ");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional contact pressure: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "contact_pressure,id,nodal_area,penetration,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional contact-pressure header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}, number(values, 0, path)});
    }
    return result;
}

std::array<fuelsim::test::FieldErrorMetrics, 4> compare_nodes(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::SteadyProblem& problem, const std::vector<double>& state,
    const std::vector<NodeReference>& reference, double& maximum_coordinate_difference) {
    if (reference.size() != mesh.nodes().size())
        throw std::invalid_argument("Three-dimensional contact node counts differ");
    std::array<fuelsim::test::FieldErrorMetrics, 4> result;
    std::vector<bool> present(mesh.nodes().size(), false);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (source >= reference.size() || reference[source].id != source || present[source])
                throw std::invalid_argument("Three-dimensional contact source-node mapping is not unique");
            present[source] = true;
            const auto& actual_point = mesh.nodes()[source];
            const auto& expected = reference[source];
            maximum_coordinate_difference = std::max(maximum_coordinate_difference,
                std::max({std::abs(actual_point.x - expected.point.x), std::abs(actual_point.y - expected.point.y),
                    std::abs(actual_point.z - expected.point.z)}));
            result[0].add(state[dofs.dof(fuelsim::Field::temperature, offset + local)], expected.fields[0]);
            result[1].add(state[dofs.dof(fuelsim::Field::displacement_x, offset + local)], expected.fields[1]);
            result[2].add(state[dofs.dof(fuelsim::Field::displacement_y, offset + local)], expected.fields[2]);
            result[3].add(state[dofs.dof(fuelsim::Field::displacement_z, offset + local)], expected.fields[3]);
        }
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Three-dimensional contact comparison did not visit every node");
    return result;
}

fuelsim::test::FieldErrorMetrics compare_pressure(const fuelsim::SteadyProblem& problem,
    const std::vector<fuelsim::CartesianContactNodeSummary>& actual, const std::vector<ContactReference>& reference,
    double& maximum_coordinate_difference) {
    const std::vector<std::size_t> source_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    if (actual.size() != source_nodes.size() || reference.size() != source_nodes.size())
        throw std::invalid_argument("Three-dimensional contact-pressure node counts differ");
    fuelsim::test::FieldErrorMetrics result;
    for (std::size_t index = 0; index < actual.size(); ++index) {
        const auto match = std::find_if(reference.begin(), reference.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[index]; });
        if (match == reference.end()) throw std::invalid_argument("MOOSE contact node is missing");
        maximum_coordinate_difference = std::max(maximum_coordinate_difference,
            std::max({std::abs(actual[index].x - match->point.x), std::abs(actual[index].y - match->point.y),
                std::abs(actual[index].z - match->point.z)}));
        result.add(actual[index].pressure, match->pressure);
    }
    return result;
}

bool run(const std::string& input_path, const std::string& thermal_path, const std::string& mechanical_path,
    const std::string& contact_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d || definition.spatial.contacts.size() != 1 ||
        !definition.spatial.contacts[0].thermal || !definition.spatial.contacts[0].mechanical ||
        definition.spatial.contacts[0].friction_coefficient != 0.2)
        throw std::invalid_argument("B3.3 comparison requires coupled three-dimensional contact and friction");
    const fuelsim::UnstructuredHex8Mesh source = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, source);
    const fuelsim::SolverOptions options = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks_per_step, definition.steady_execution.minimum_load_increment},
        options);
    bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 10,
                      "B3.3 coupled contact path completes ten load steps") &&
                  check(solve.aggregate_timing.workspace_setups == 1, "B3.3 contact path reuses one PETSc workspace");
    double coordinate_error = 0.0;
    const auto fields =
        compare_nodes(source, problem, solve.solve.state, read_nodes(thermal_path, mechanical_path), coordinate_error);
    const std::vector<fuelsim::CartesianContactNodeSummary> contact =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const fuelsim::test::FieldErrorMetrics pressure =
        compare_pressure(problem, contact, read_contact(contact_path), coordinate_error);
    constexpr double tolerance = 5.0e-3;
    const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t field = 0; field < fields.size(); ++field) {
        fuelsim::test::print_relative_metrics("b33_" + names[field], fields[field]);
        passed = check(fuelsim::test::relative_metrics_below(fields[field], tolerance),
                     "B3.3 " + names[field] + " three MOOSE errors are below 0.5 percent") &&
                 passed;
    }
    fuelsim::test::print_relative_metrics("b33_contact_pressure", pressure);
    passed = check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                 "B3.3 contact-pressure three MOOSE errors are below 0.5 percent") &&
             check(coordinate_error < 1.0e-12, "B3.3 compares MOOSE values at matching coordinates") && passed;
    std::size_t active = 0;
    std::size_t sliding = 0;
    double maximum_capacity_excess = 0.0;
    for (const auto& node : contact) {
        if (!(node.pressure > 0.0)) continue;
        ++active;
        if (node.sliding) ++sliding;
        maximum_capacity_excess = std::max(maximum_capacity_excess,
            node.tangential_traction - definition.spatial.contacts[0].friction_coefficient * node.pressure);
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    passed = check(active == 4 && sliding == 0 && interface.active_contact_nodes == 4,
                 "B3.3 has four active projected nodes in the Coulomb sticking branch") &&
             check(maximum_capacity_excess <= 1.0e-12 * interface.maximum_contact_pressure,
                 "B3.3 tangential traction respects the Coulomb cap") &&
             check(interface.total_heat_rate > 0.0 && interface.total_contact_force > 0.0 &&
                       interface.total_tangential_force > 0.0,
                 "B3.3 produces nonzero heat, normal-force, and friction-force transfer") &&
             passed;
    fuelsim::SpatialDefinition frictionless_definition = definition.spatial;
    frictionless_definition.contacts[0].friction_coefficient = 0.0;
    fuelsim::SteadyProblem frictionless(std::move(frictionless_definition), source);
    const fuelsim::SteadyResult frictionless_solve = fuelsim::solve_steady(frictionless,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks_per_step, definition.steady_execution.minimum_load_increment},
        options);
    double difference_squared = 0.0;
    double scale_squared = 0.0;
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t node = 0; node < dofs.node_count(); ++node) {
        const std::size_t dof = dofs.dof(fuelsim::Field::displacement_y, node);
        const double difference = solve.solve.state[dof] - frictionless_solve.solve.state[dof];
        difference_squared += difference * difference;
        scale_squared += solve.solve.state[dof] * solve.solve.state[dof];
    }
    const double friction_effect = std::sqrt(difference_squared / scale_squared);
    passed = check(frictionless_solve.completed && friction_effect > 1.0e-3,
                 "B3.3 friction measurably changes the three-dimensional tangential field") &&
             passed;
    std::cout << "b33_active_contact_nodes=" << active << '\n'
              << "b33_sliding_contact_nodes=" << sliding << '\n'
              << "b33_total_heat_rate=" << interface.total_heat_rate << '\n'
              << "b33_total_contact_force=" << interface.total_contact_force << '\n'
              << "b33_total_tangential_force=" << interface.total_tangential_force << '\n'
              << "b33_frictional_tangential_field_change=" << friction_effect << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_b33_hex8_contact_moose_tests "
                     "<case.fsi> <thermal-nodes.csv> <mechanical-nodes.csv> <contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B3.3 three-dimensional contact MOOSE comparison\n");
        if (!run(argv[1], argv[2], argv[3], argv[4])) return 1;
        std::cout << "[PASS] B3.3 three-dimensional contact MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B3.3 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
