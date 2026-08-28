#include "fuelsim/core/cartesian3d_hex20.hpp"
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

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) result.push_back(field);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("H20.35 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("H20.35 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("H20.35 CSV contains an invalid node index: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<DisplacementReference> read_displacements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.35 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "id,x,y,z,disp_x,disp_y,disp_z")
        throw std::invalid_argument("Unexpected H20.35 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.35 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2")
        throw std::invalid_argument("Unexpected H20.35 contact header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
                {number(values, 7, path), number(values, 8, path), number(values, 9, path)}, number(values, 10, path),
                number(values, 11, path)});
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

double coordinate_difference(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return std::max({std::abs(first.x - second.x), std::abs(first.y - second.y), std::abs(first.z - second.z)});
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
    if (!(measure > 0.0)) throw std::invalid_argument("H20.35 contact face is degenerate");
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

struct ContactNumericalEvidence final {
    std::array<double, 3> residual_balance{};
    double maximum_directional_error = 0.0;
};

ContactNumericalEvidence inspect_contact_numerics(const fuelsim::SteadyProblem& problem,
    const fuelsim::cartesian::SpatialAssembly& spatial, const std::vector<double>& state) {
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
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        constexpr double perturbation = 1.0e-10;
        std::vector<double> plus = local, minus = local;
        for (std::size_t index = 0; index < local.size(); ++index) {
            plus[index] += perturbation * direction[index];
            minus[index] -= perturbation * direction[index];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                analytic += jacobian[row * local.size() + column] * direction[column];
            const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += (analytic - reference) * (analytic - reference);
            reference_squared += reference * reference;
        }
        if (!(reference_squared > 0.0)) throw std::logic_error("H20.35 contact Jacobian direction is zero");
        result.maximum_directional_error =
            std::max(result.maximum_directional_error, std::sqrt(difference_squared / reference_squared));
        for (std::size_t row = 0; row < dofs.size(); ++row)
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (dofs[row] >= field.begin && dofs[row] < field.end)
                    result.residual_balance[component] += residual[row];
            }
    }
    return result;
}

bool metric_passes(const std::string& name, const fuelsim::test::FieldErrorMetrics& metric, double relative_tolerance,
    double zero_tolerance) {
    if (metric.has_relative_norm()) {
        fuelsim::test::print_relative_metrics(name, metric);
        return check(fuelsim::test::relative_metrics_below(metric, relative_tolerance) &&
                         metric.maximum_zero_reference_difference < zero_tolerance,
            name + " passes all three relative metrics and its separate zero-reference check");
    }
    fuelsim::test::print_absolute_metrics(name, metric);
    return check(metric.maximum_zero_reference_difference < zero_tolerance,
        name + " passes its theoretical-zero absolute check");
}

bool run_case(const std::string& case_path, const std::string& displacement_path, const std::string& contact_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult result =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(
        result.completed && result.solve.converged && result.completed_steps == definition.steady_execution.load_steps,
        "H20.35 quadratic-cylinder friction solve completes every load step");
    const std::vector<DisplacementReference> displacement = read_displacements(displacement_path);
    const std::vector<ContactReference> contact = read_contact(contact_path);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_error;
    fuelsim::test::GroupedFieldErrorMetrics displacement_vector_error;
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_error = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto reference = std::find_if(displacement.begin(), displacement.end(),
                [source](const DisplacementReference& value) { return value.id == source; });
            if (reference == displacement.end() || present[source])
                throw std::invalid_argument("H20.35 displacement source-node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_error =
                std::max(maximum_coordinate_error, coordinate_difference(mesh.nodes()[source], reference->point));
            const std::size_t global = spatial.global_node(region, local);
            std::array<double, 3> actual_displacement{}, reference_displacement{};
            for (std::size_t component = 0; component < 3; ++component) {
                actual_displacement[component] = result.solve.state[fields[component + 1].begin + global];
                reference_displacement[component] = reference->displacement[component];
                displacement_error[component].add(actual_displacement[component], reference_displacement[component]);
            }
            displacement_vector_error.add(actual_displacement.data(), reference_displacement.data(), 3);
        }
    }

    const auto summaries = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, result.solve.state);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    if (summaries.size() != contact.size() || source_nodes.size() != contact.size())
        throw std::invalid_argument("H20.35 contact-node counts differ");
    fuelsim::test::FieldErrorMetrics normal_force, circumferential_force, axial_force, slip_1, slip_2;
    double actual_normal_resultant = 0.0, reference_normal_resultant = 0.0;
    double actual_circumferential_resultant = 0.0, reference_circumferential_resultant = 0.0;
    double actual_axial_resultant = 0.0, reference_axial_resultant = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const auto reference = std::find_if(contact.begin(), contact.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[node]; });
        if (reference == contact.end()) throw std::invalid_argument("H20.35 contact-node mapping is incomplete");
        maximum_coordinate_error = std::max(
            maximum_coordinate_error, coordinate_difference(mesh.nodes()[source_nodes[node]], reference->point));
        const double radius = std::hypot(reference->point.x, reference->point.y);
        const std::array<double, 3> radial = {reference->point.x / radius, reference->point.y / radius, 0.0};
        const std::array<double, 3> circumferential = {-reference->point.y / radius, reference->point.x / radius, 0.0};
        const std::array<double, 3> actual_normal = {-summaries[node].normal_contact_force[0],
            -summaries[node].normal_contact_force[1], -summaries[node].normal_contact_force[2]};
        const std::array<double, 3> actual_tangent = {-summaries[node].tangential_contact_force[0],
            -summaries[node].tangential_contact_force[1], -summaries[node].tangential_contact_force[2]};
        const double actual_normal_value = dot(actual_normal, radial);
        const double reference_normal_value = dot(reference->normal_force, radial);
        const double actual_circumferential = dot(actual_tangent, circumferential);
        const double reference_circumferential = dot(reference->tangential_force, circumferential);
        normal_force.add(actual_normal_value, reference_normal_value);
        circumferential_force.add(actual_circumferential, reference_circumferential);
        axial_force.add(actual_tangent[2], reference->tangential_force[2]);
        // On the C3D20 S6 faces in this tracked mesh, Abaqus local direction 1
        // is opposite the positive circumferential direction and local
        // direction 2 is opposite the global axial direction.
        slip_1.add(-dot(summaries[node].tangential_slip, circumferential), reference->slip_1);
        slip_2.add(-summaries[node].tangential_slip[2], reference->slip_2);
        actual_normal_resultant += actual_normal_value;
        reference_normal_resultant += reference_normal_value;
        actual_circumferential_resultant += actual_circumferential;
        reference_circumferential_resultant += reference_circumferential;
        actual_axial_resultant += actual_tangent[2];
        reference_axial_resultant += reference->tangential_force[2];
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, result.solve.state);
    const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    const ContactNumericalEvidence numerical = inspect_contact_numerics(problem, spatial, result.solve.state);
    const double face_nonplanarity = maximum_secondary_face_nonplanarity(mesh);
    constexpr double relative_tolerance = 1.0e-2;
    constexpr double zero_tolerance = 1.0e-10;
    for (std::size_t component = 0; component < 3; ++component)
        fuelsim::test::print_relative_metrics(
            "h20_35_displacement_" + std::string(1, "xyz"[component]), displacement_error[component]);
    fuelsim::test::print_grouped_relative_metrics("h20_35_displacement_vector", displacement_vector_error);
    passed = check(fuelsim::test::grouped_relative_metrics_below(displacement_vector_error, relative_tolerance) &&
                       displacement_vector_error.maximum_zero_reference_difference < zero_tolerance,
                 "H20.35 complete displacement-vector metrics and its separate zero-reference check are below 1 "
                 "percent") &&
             passed;
    passed =
        metric_passes("h20_35_signed_normal_radial_force", normal_force, relative_tolerance, zero_tolerance) && passed;
    passed = metric_passes("h20_35_signed_tangential_circumferential_force", circumferential_force, relative_tolerance,
                 zero_tolerance) &&
             passed;
    passed = metric_passes("h20_35_signed_tangential_axial_force", axial_force, relative_tolerance, zero_tolerance) &&
             passed;
    passed = metric_passes("h20_35_tangential_slip_1", slip_1, relative_tolerance, zero_tolerance) && passed;
    passed = metric_passes("h20_35_tangential_slip_2", slip_2, relative_tolerance, zero_tolerance) && passed;
    const double normal_resultant_error =
        std::abs(actual_normal_resultant - reference_normal_resultant) / std::abs(reference_normal_resultant);
    const double circumferential_resultant_error =
        std::abs(actual_circumferential_resultant - reference_circumferential_resultant) /
        std::abs(reference_circumferential_resultant);
    const double axial_resultant_error =
        std::abs(actual_axial_resultant - reference_axial_resultant) / std::abs(reference_axial_resultant);
    std::cout << "h20_35_normal_resultant_relative_error=" << normal_resultant_error << '\n'
              << "h20_35_circumferential_resultant_relative_error=" << circumferential_resultant_error << '\n'
              << "h20_35_axial_resultant_relative_error=" << axial_resultant_error << '\n'
              << "h20_35_contact_residual_balance=" << numerical.residual_balance[0] << ','
              << numerical.residual_balance[1] << ',' << numerical.residual_balance[2] << '\n'
              << "h20_35_contact_jacobian_directional_error=" << numerical.maximum_directional_error << '\n'
              << "h20_35_maximum_secondary_face_nonplanarity=" << face_nonplanarity << '\n';
    return check(displacement.size() == mesh.nodes().size() &&
                     std::all_of(present.begin(), present.end(), [](bool value) { return value; }),
               "H20.35 compares every tracked Exodus mesh node") &&
           check(maximum_coordinate_error < 1.0e-12, "H20.35 references preserve the tracked Exodus coordinates") &&
           check(interface.active_contact_nodes == summaries.size() && interface.unprojected_contact_nodes == 0,
               "H20.35 keeps all 37 quadratic-surface contact nodes active and projected") &&
           check(
               face_nonplanarity > 1.0e-3, "H20.35 uses genuinely quadratic contact faces rather than planar facets") &&
           check(histories.size() == summaries.size() &&
                     std::all_of(histories.begin(), histories.end(),
                         [](const fuelsim::ContactPointHistory& history) { return !history.sliding; }),
               "H20.35 keeps every constraint in the two-direction sticking branch") &&
           check(std::abs(reference_circumferential_resultant) > 1.0 && std::abs(reference_axial_resultant) > 1.0,
               "H20.35 activates nonzero circumferential and axial tangential resultants") &&
           check(normal_resultant_error < relative_tolerance && circumferential_resultant_error < relative_tolerance &&
                     axial_resultant_error < relative_tolerance,
               "H20.35 signed normal and both tangential resultants agree with Abaqus below 1 percent") &&
           check(std::abs(numerical.residual_balance[0]) < 1.0e-8 && std::abs(numerical.residual_balance[1]) < 1.0e-8 &&
                     std::abs(numerical.residual_balance[2]) < 1.0e-8,
               "H20.35 frictional contact residual is action-reaction conservative") &&
           check(numerical.maximum_directional_error < 1.0e-7,
               "H20.35 frictional contact Jacobian matches a centered directional difference") &&
           passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_h20_curved_friction_abaqus_tests <case.fsi> <displacement.csv> <contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.35 curved friction Abaqus comparison\n");
        if (!run_case(argv[1], argv[2], argv[3])) return 1;
        std::cout << "[PASS] fuelsim H20.35 quadratic curved friction comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.35 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
