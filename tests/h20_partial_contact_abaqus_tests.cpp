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
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Reference final {
    std::size_t node;
    fuelsim::CartesianPoint3 point;
    double gap, pressure;
    std::array<double, 3> normal_force;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete H20.41 reference row: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid H20.41 reference value: " + path);
    return result;
}

std::vector<Reference> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.41 contact reference: " + path);
    std::string line;
    const std::string expected = "id,x,y,z,gap,pressure,normal_x,normal_y,normal_z";
    if (!std::getline(input, line) || line != expected)
        throw std::invalid_argument("Unexpected H20.41 reference header: " + path);
    std::vector<Reference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({static_cast<std::size_t>(number(values, 0, path)),
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}, number(values, 4, path),
            number(values, 5, path), {number(values, 6, path), number(values, 7, path), number(values, 8, path)}});
    }
    if (result.size() != 29) throw std::invalid_argument("H20.41 requires 29 secondary contact nodes");
    return result;
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = definition.solver.absolute_tolerance;
    result.relative_tolerance = definition.solver.relative_tolerance;
    result.step_tolerance = definition.solver.step_tolerance;
    result.maximum_iterations = definition.solver.maximum_iterations;
    result.linear_solver = definition.solver.linear_solver;
    result.direct_factorization = definition.solver.direct_factorization;
    result.field_residual_scaling = definition.solver.field_residual_scaling;
    result.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    return result;
}

double coordinate_difference(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return std::max({std::abs(first.x - second.x), std::abs(first.y - second.y), std::abs(first.z - second.z)});
}

bool run(const std::string& case_path, const std::string& reference_path) {
    const std::vector<Reference> references = read_reference(reference_path);
    std::map<std::size_t, Reference> reference_by_node;
    for (const Reference& reference : references)
        if (!reference_by_node.emplace(reference.node, reference).second)
            throw std::invalid_argument("H20.41 contains a repeated reference node");

    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(
        solve.completed && solve.solve.converged && solve.completed_steps == definition.steady_execution.load_steps,
        "H20.41 completes every Fuelsim load step");
    if (!solve.completed) return false;

    const auto summaries = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    if (summaries.size() != source_nodes.size() || summaries.size() != references.size())
        throw std::logic_error("H20.41 contact summary and reference sizes differ");

    const double penalty = definition.spatial.contacts.at(0).penalty;
    fuelsim::test::FieldErrorMetrics constraint_pressure, recovered_pressure;
    fuelsim::test::GroupedFieldErrorMetrics normal_force;
    std::size_t fuelsim_constraint_active = 0, fuelsim_recovered_positive = 0, abaqus_constraint_active = 0,
                abaqus_recovered_positive = 0, abaqus_recovery_extension = 0;
    double coordinate_error = 0.0, abaqus_resultant = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const Reference& reference = reference_by_node.at(source_nodes[node]);
        coordinate_error =
            std::max(coordinate_error, coordinate_difference(mesh.nodes().at(reference.node), reference.point));
        const double actual_constraint = std::max(-penalty * summaries[node].gap, 0.0);
        const double reference_constraint = std::max(-penalty * reference.gap, 0.0);
        constraint_pressure.add(actual_constraint, reference_constraint);
        recovered_pressure.add(summaries[node].pressure, reference.pressure);
        std::array<double, 3> actual_normal_force{};
        for (std::size_t component = 0; component < actual_normal_force.size(); ++component)
            actual_normal_force[component] = -summaries[node].normal_contact_force[component];
        normal_force.add(actual_normal_force.data(), reference.normal_force.data(), actual_normal_force.size());
        fuelsim_constraint_active += actual_constraint > 0.0 ? 1U : 0U;
        fuelsim_recovered_positive += summaries[node].pressure > 0.0 ? 1U : 0U;
        abaqus_constraint_active += reference_constraint > 0.0 ? 1U : 0U;
        abaqus_recovered_positive += reference.pressure > 0.0 ? 1U : 0U;
        abaqus_recovery_extension += reference.gap > 0.0 && reference.pressure > 0.0 ? 1U : 0U;
        abaqus_resultant += reference.normal_force[0];
    }

    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const fuelsim::Hex20RegionMesh& secondary = spatial.hex20_region_mesh(1);
    const fuelsim::Hex20RegionBoundary contact = secondary.map_side_set(mesh, "secondary_contact");
    std::size_t partial_faces = 0;
    for (const fuelsim::Quad8FaceElement& face : contact.faces) {
        std::size_t active = 0;
        for (std::size_t local_node : face.nodes) {
            const std::size_t source = secondary.source_node_ids().at(local_node);
            active += reference_by_node.at(source).gap < 0.0 ? 1U : 0U;
        }
        partial_faces += active > 0 && active < face.nodes.size() ? 1U : 0U;
    }

    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    const double resultant_error =
        std::abs(interface.total_contact_force - std::abs(abaqus_resultant)) / std::abs(abaqus_resultant);
    fuelsim::test::print_relative_metrics("h20_41_constraint_pressure", constraint_pressure);
    fuelsim::test::print_relative_metrics("h20_41_recovered_pressure", recovered_pressure);
    fuelsim::test::print_grouped_relative_metrics("h20_41_normal_force", normal_force);
    std::cout << "h20_41_fuelsim_constraint_active_nodes=" << fuelsim_constraint_active << '\n'
              << "h20_41_fuelsim_recovered_positive_nodes=" << fuelsim_recovered_positive << '\n'
              << "h20_41_abaqus_constraint_active_nodes=" << abaqus_constraint_active << '\n'
              << "h20_41_abaqus_recovered_positive_nodes=" << abaqus_recovered_positive << '\n'
              << "h20_41_abaqus_recovery_extension_nodes=" << abaqus_recovery_extension << '\n'
              << "h20_41_partial_faces=" << partial_faces << '\n'
              << "h20_41_fuelsim_resultant=" << interface.total_contact_force << '\n'
              << "h20_41_abaqus_resultant=" << std::abs(abaqus_resultant) << '\n'
              << "h20_41_resultant_relative_error=" << resultant_error << '\n';

    constexpr double tolerance = 1.0e-2;
    passed = check(coordinate_error < 1.0e-14, "H20.41 uses identical Abaqus and Fuelsim coordinates") && passed;
    passed =
        check(fuelsim_constraint_active == 11 && fuelsim_recovered_positive == 11 && abaqus_constraint_active == 11 &&
                  abaqus_recovered_positive == 16 && abaqus_recovery_extension == 5 && partial_faces == 4,
            "H20.41 places the contact front inside quadratic secondary faces") &&
        passed;
    passed = check(interface.projected_contact_nodes == references.size() && interface.unprojected_contact_nodes == 0 &&
                       interface.active_contact_nodes == fuelsim_constraint_active,
                 "H20.41 projects every constraint and reports physical active-node topology") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(constraint_pressure, tolerance),
                 "H20.41 constraint-pressure metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::grouped_relative_metrics_below(normal_force, tolerance),
                 "H20.41 nodal normal-force metrics are below one percent") &&
             passed;
    passed = check(recovered_pressure.maximum_zero_reference_difference == 0.0,
                 "H20.41 reports zero pressure at every Abaqus zero-pressure node") &&
             passed;
    passed =
        check(resultant_error < tolerance, "H20.41 total normal force differs from Abaqus by less than one percent") &&
        passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_h20_41_hex20_partial_contact_abaqus_tests <case.fsi> <contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.41 HEX20 partial-contact Abaqus comparison\n");
        const bool passed = run(argv[1], argv[2]);
        if (passed && session.rank() == 0) std::cout << "[PASS] H20.41 HEX20 partial-contact Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.41 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
