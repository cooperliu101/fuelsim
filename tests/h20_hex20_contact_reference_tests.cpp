#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/field_error_metrics.hpp"
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
struct TemperatureReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double temperature;
};

struct DisplacementReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> displacement{};
};

struct ContactReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double pressure;
    double nodal_area;
};

struct ReactionReference final {
    double normal_force;
    double tangential_force;
};

struct AreaRuleResult final {
    fuelsim::test::FieldErrorMetrics displacement_x;
    double normal_force = 0.0;
    double normal_force_relative_error = 0.0;
    std::size_t active_nodes = 0;
    std::size_t positive_area_nodes = 0;
    std::size_t negative_area_nodes = 0;
    bool passed = false;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        values.push_back(value);
    return values;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete HEX20 contact row in " + path);
    return std::stod(values[index]);
}

std::vector<TemperatureReference> read_temperature(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read HEX20 thermal-contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "T,id,x,y,z")
        throw std::invalid_argument("Unexpected thermal-contact header in " + path);
    std::vector<TemperatureReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            number(values, 0, path)});
    }
    return result;
}

std::vector<DisplacementReference> read_displacement(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read HEX20 mechanical-contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected mechanical-contact header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 3, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read HEX20 contact-pressure reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "contact_pressure,id,nodal_area,penetration,x,y,z")
        throw std::invalid_argument("Unexpected contact-pressure header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            number(values, 0, path),
            number(values, 2, path)});
    }
    return result;
}

ReactionReference read_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read HEX20 sliding-contact reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected sliding-contact reaction header in " + path);
    std::vector<std::string> final_values;
    while (std::getline(input, line))
        if (!line.empty())
            final_values = split_csv(line);
    if (final_values.empty() || std::abs(number(final_values, 0, path) - 1.0) > 1.0e-12)
        throw std::invalid_argument("HEX20 sliding-contact reference does not end at unit load");
    return {std::abs(number(final_values, 1, path)),
        std::hypot(number(final_values, 2, path), number(final_values, 3, path))};
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.preconditioner = definition.solver.preconditioner;
    options.field_residual_scaling = definition.solver.field_residual_scaling;
    return options;
}

bool run_positive_lumped_path(const fuelsim::FuelSimCaseDefinition& comparison_definition,
    const fuelsim::UnstructuredHex20Mesh& mesh,
    std::size_t expected_steps,
    const std::string& label) {
    fuelsim::FuelSimCaseDefinition definition = comparison_definition;
    definition.spatial.contacts[0].mechanical_discretization =
        fuelsim::MechanicalContactDiscretization::node_to_surface;
    definition.spatial.contacts[0].quad8_nodal_area_rule = fuelsim::Quad8NodalAreaRule::positive_lumped;
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    const auto contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    std::size_t active = 0, sliding = 0, positive_areas = 0;
    double area_sum = 0.0;
    for (const auto& node : contact) {
        if (node.tributary_area > 0.0)
            ++positive_areas;
        area_sum += node.tributary_area;
        if (!(node.pressure > 0.0))
            continue;
        ++active;
        if (node.sliding)
            ++sliding;
    }
    bool passed =
        check(solve.completed && solve.solve.converged && solve.completed_steps == expected_steps,
            label + " positive-lumped path completes its load steps")
        && check(contact.size() == 8 && positive_areas == 8 && area_sum > 0.0,
            label + " positive-lumped path keeps all eight nodal areas positive")
        && check(active > 0 && interface.active_contact_nodes == active && interface.total_contact_force > 0.0,
            label + " positive-lumped path activates projected contact nodes");
    if (definition.spatial.contacts[0].friction_coefficient > 0.0)
        passed = check(sliding > 0 && interface.total_tangential_force > 0.0,
                     label + " positive-lumped path enters Coulomb sliding")
                 && passed;
    std::cout << label << "_positive_lumped_area_sum=" << area_sum << '\n'
              << label << "_positive_lumped_active_contact_nodes=" << active << '\n'
              << label << "_positive_lumped_sliding_contact_nodes=" << sliding << '\n';
    return passed;
}

double coordinate_difference(const fuelsim::CartesianPoint3& actual, const fuelsim::CartesianPoint3& expected) {
    return std::max(
        {std::abs(actual.x - expected.x), std::abs(actual.y - expected.y), std::abs(actual.z - expected.z)});
}

bool run_thermal(const std::string& case_path, const std::string& reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || !definition.spatial.contacts[0].thermal
        || definition.spatial.contacts[0].mechanical)
        throw std::invalid_argument("H20.16 requires isolated HEX20 thermal contact");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
                      "H20.16 thermal-contact solve completes one load step")
                  && check(solve.aggregate_timing.workspace_setups == 1,
                      "H20.16 thermal-contact solve constructs one PETSc workspace");
    const auto reference = read_temperature(reference_path);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    fuelsim::test::FieldErrorMetrics metrics;
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_difference = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            if (!region_mesh.temperature_nodes()[local])
                continue;
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto found = std::find_if(reference.begin(), reference.end(), [&](const TemperatureReference& value) {
                return value.id == source;
            });
            if (found == reference.end() || present[source])
                throw std::invalid_argument("HEX20 thermal source-node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_difference =
                std::max(maximum_coordinate_difference, coordinate_difference(mesh.nodes().at(source), found->point));
            metrics.add(solve.solve.state.at(fields[0].begin + spatial.global_temperature_node(region, local)),
                found->temperature);
        }
    }
    constexpr double tolerance = 5.0e-3;
    fuelsim::test::print_relative_metrics("h20_16_temperature", metrics);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    passed =
        check(metrics.value_count == reference.size(), "H20.16 compares all sixteen first-order temperature nodes")
        && check(fuelsim::test::relative_metrics_below(metrics, tolerance),
            "H20.16 temperature three MOOSE errors are below 0.5 percent")
        && check(maximum_coordinate_difference < 1.0e-12, "H20.16 compares temperature on the tracked MOOSE HEX20 mesh")
        && check(interface.total_heat_rate > 0.0, "H20.16 transfers nonzero heat across the interface") && passed;
    std::cout << "h20_16_total_heat_rate=" << interface.total_heat_rate << '\n';
    return passed;
}

bool run_mechanical(const std::string& case_path,
    const std::string& displacement_path,
    const std::string& contact_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || definition.spatial.contacts[0].thermal
        || !definition.spatial.contacts[0].mechanical || definition.spatial.contacts[0].friction_coefficient != 0.0
        || definition.spatial.contacts[0].quad8_nodal_area_rule != fuelsim::Quad8NodalAreaRule::consistent_shape)
        throw std::invalid_argument("H20.17 requires isolated frictionless HEX20 mechanical contact");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 10,
                      "H20.17 mechanical-contact solve completes ten load steps")
                  && check(solve.aggregate_timing.workspace_setups == 1,
                      "H20.17 mechanical-contact solve reuses one PETSc workspace");
    const auto displacement = read_displacement(displacement_path);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_metrics;
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_difference = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto found = std::find_if(displacement.begin(),
                displacement.end(),
                [&](const DisplacementReference& value) { return value.id == source; });
            if (found == displacement.end() || present[source])
                throw std::invalid_argument("HEX20 displacement source-node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_difference =
                std::max(maximum_coordinate_difference, coordinate_difference(mesh.nodes().at(source), found->point));
            const std::size_t global = spatial.global_node(region, local);
            for (std::size_t component = 0; component < displacement_metrics.size(); ++component)
                displacement_metrics[component].add(solve.solve.state.at(fields[component + 1].begin + global),
                    found->displacement[component]);
        }
    }
    const auto contact_reference = read_contact(contact_path);
    const auto contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    if (contact.size() != source_nodes.size() || contact_reference.size() != source_nodes.size())
        throw std::invalid_argument("HEX20 contact-pressure node counts differ");
    fuelsim::test::FieldErrorMetrics pressure_metrics;
    std::size_t active = 0;
    std::size_t sliding = 0;
    for (std::size_t index = 0; index < contact.size(); ++index) {
        const auto found = std::find_if(contact_reference.begin(),
            contact_reference.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[index]; });
        if (found == contact_reference.end())
            throw std::invalid_argument("MOOSE HEX20 contact node is missing");
        maximum_coordinate_difference = std::max(maximum_coordinate_difference,
            coordinate_difference({contact[index].x, contact[index].y, contact[index].z}, found->point));
        pressure_metrics.add(contact[index].pressure, found->pressure);
        std::cout << "h20_17_contact_node=" << source_nodes[index] << ",gap=" << contact[index].gap
                  << ",pressure=" << contact[index].pressure << ",reference_pressure=" << found->pressure
                  << ",area=" << contact[index].tributary_area << '\n';
        if (!(contact[index].pressure > 0.0))
            continue;
        ++active;
        if (contact[index].sliding)
            ++sliding;
    }
    constexpr double tolerance = 5.0e-3;
    const std::array<std::string, 3> names = {"displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t component = 0; component < displacement_metrics.size(); ++component) {
        fuelsim::test::print_relative_metrics("h20_17_" + names[component], displacement_metrics[component]);
        passed = check(fuelsim::test::relative_metrics_below(displacement_metrics[component], tolerance),
                     "H20.17 " + names[component] + " three MOOSE errors are below 0.5 percent")
                 && passed;
    }
    fuelsim::test::print_relative_metrics("h20_17_contact_pressure", pressure_metrics);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    passed =
        check(displacement.size() == mesh.nodes().size(), "H20.17 compares all forty second-order displacement nodes")
        && check(fuelsim::test::relative_metrics_below(pressure_metrics, tolerance),
            "H20.17 contact-pressure three MOOSE errors are below 0.5 percent")
        && check(maximum_coordinate_difference < 1.0e-12,
            "H20.17 compares displacement and pressure on the tracked MOOSE HEX20 mesh")
        && check(active > 0 && interface.active_contact_nodes == active,
            "H20.17 activates projected quadratic-face contact nodes")
        && check(interface.total_contact_force > 0.0 && interface.total_tangential_force == 0.0,
            "H20.17 transfers normal force without friction force")
        && passed;
    std::cout << "h20_17_active_contact_nodes=" << active << '\n'
              << "h20_17_sliding_contact_nodes=" << sliding << '\n'
              << "h20_17_total_contact_force=" << interface.total_contact_force << '\n'
              << "h20_17_total_tangential_force=" << interface.total_tangential_force << '\n';
    return run_positive_lumped_path(definition, mesh, 10, "h20_17") && passed;
}

bool run_sliding(const std::string& case_path, const std::string& reaction_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || definition.spatial.contacts[0].thermal
        || !definition.spatial.contacts[0].mechanical || definition.spatial.contacts[0].friction_coefficient != 0.001
        || definition.spatial.contacts[0].quad8_nodal_area_rule != fuelsim::Quad8NodalAreaRule::consistent_shape)
        throw std::invalid_argument("H20.18 requires isolated HEX20 Coulomb sliding contact");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
        "H20.18 sliding-contact solve completes one load step");
    const auto contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    std::size_t active = 0;
    std::size_t sliding = 0;
    for (const auto& node : contact) {
        if (!(node.pressure > 0.0))
            continue;
        ++active;
        if (node.sliding)
            ++sliding;
    }
    const ReactionReference reference = read_reaction(reaction_path);
    const double normal_error =
        std::abs(interface.total_contact_force - reference.normal_force) / reference.normal_force;
    const double tangential_error =
        std::abs(interface.total_tangential_force - reference.tangential_force) / reference.tangential_force;
    constexpr double tolerance = 5.0e-3;
    passed = check(active == 8 && sliding > 0 && interface.active_contact_nodes == 8,
                 "H20.18 keeps all eight quadratic-face nodes active and enters Coulomb sliding")
             && check(normal_error < tolerance && tangential_error < tolerance,
                 "H20.18 normal and tangential resultants agree with MOOSE below 0.5 percent")
             && check(interface.total_tangential_force
                          > definition.spatial.contacts[0].friction_coefficient * interface.total_contact_force,
                 "H20.18 signed quadratic nodal areas retain the MOOSE edge-midpoint friction resultant")
             && passed;
    std::cout << "h20_18_active_contact_nodes=" << active << '\n'
              << "h20_18_sliding_contact_nodes=" << sliding << '\n'
              << "h20_18_normal_resultant=" << interface.total_contact_force << '\n'
              << "h20_18_reference_normal_resultant=" << reference.normal_force << '\n'
              << "h20_18_tangential_resultant=" << interface.total_tangential_force << '\n'
              << "h20_18_reference_tangential_resultant=" << reference.tangential_force << '\n'
              << "h20_18_normal_resultant_relative_error=" << normal_error << '\n'
              << "h20_18_tangential_resultant_relative_error=" << tangential_error << '\n';
    return run_positive_lumped_path(definition, mesh, 1, "h20_18") && passed;
}

AreaRuleResult run_area_rule(const fuelsim::FuelSimCaseDefinition& comparison_definition,
    const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::vector<DisplacementReference>& displacement,
    const ReactionReference& reaction,
    fuelsim::Quad8NodalAreaRule area_rule,
    double secondary_reference_x_shift,
    const std::string& label,
    fuelsim::MechanicalContactDiscretization discretization =
        fuelsim::MechanicalContactDiscretization::node_to_surface) {
    fuelsim::FuelSimCaseDefinition definition = comparison_definition;
    definition.spatial.contacts[0].mechanical_discretization = discretization;
    definition.spatial.contacts[0].quad8_nodal_area_rule = area_rule;
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    AreaRuleResult result;
    result.passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
                        label + " completes the pure-normal load step")
                    && check(solve.aggregate_timing.workspace_setups == 1, label + " constructs one PETSc workspace");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_difference = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto found = std::find_if(displacement.begin(),
                displacement.end(),
                [&](const DisplacementReference& value) { return value.id == source; });
            if (found == displacement.end() || present[source])
                throw std::invalid_argument("HEX20 contact displacement source-node mapping is incomplete or repeated");
            present[source] = true;
            fuelsim::CartesianPoint3 expected_point = mesh.nodes().at(source);
            if (region == 1)
                expected_point.x += secondary_reference_x_shift;
            maximum_coordinate_difference =
                std::max(maximum_coordinate_difference, coordinate_difference(expected_point, found->point));
            const std::size_t global = spatial.global_node(region, local);
            result.displacement_x.add(solve.solve.state.at(fields[1].begin + global), found->displacement[0]);
        }
    }
    const auto contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    for (const auto& node : contact) {
        if (node.pressure > 0.0)
            ++result.active_nodes;
        if (node.tributary_area > 0.0)
            ++result.positive_area_nodes;
        else if (node.tributary_area < 0.0)
            ++result.negative_area_nodes;
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    result.normal_force = interface.total_contact_force;
    result.normal_force_relative_error = std::abs(result.normal_force - reaction.normal_force) / reaction.normal_force;
    result.passed =
        check(displacement.size() == mesh.nodes().size(), label + " compares all forty displacement nodes")
        && check(maximum_coordinate_difference < 1.0e-9, label + " uses the expected HEX20 reference coordinates")
        && check(result.active_nodes > 0 && interface.active_contact_nodes == result.active_nodes,
            label + " activates its projected secondary face nodes")
        && check(interface.total_tangential_force < 1.0e-10, label + " remains frictionless under pure-normal loading")
        && result.passed;
    fuelsim::test::print_relative_metrics(label + "_displacement_x", result.displacement_x);
    std::cout << label << "_active_contact_nodes=" << result.active_nodes << '\n'
              << label << "_positive_area_nodes=" << result.positive_area_nodes << '\n'
              << label << "_negative_area_nodes=" << result.negative_area_nodes << '\n'
              << label << "_normal_resultant=" << result.normal_force << '\n'
              << label << "_reference_normal_resultant=" << reaction.normal_force << '\n'
              << label << "_normal_resultant_relative_error=" << result.normal_force_relative_error << '\n';
    return result;
}

bool run_mortar_area_comparison(const std::string& case_path,
    const std::string& displacement_path,
    const std::string& reaction_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || definition.spatial.contacts[0].thermal
        || !definition.spatial.contacts[0].mechanical || definition.spatial.contacts[0].friction_coefficient != 0.0
        || definition.spatial.contacts[0].quad8_nodal_area_rule != fuelsim::Quad8NodalAreaRule::positive_lumped)
        throw std::invalid_argument("H20.19 requires isolated pure-normal HEX20 mechanical contact");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    const auto displacement = read_displacement(displacement_path);
    const ReactionReference reaction = read_reaction(reaction_path);
    const AreaRuleResult positive = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::positive_lumped,
        0.0,
        "h20_19_positive_lumped");
    const AreaRuleResult consistent = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::consistent_shape,
        0.0,
        "h20_19_consistent_shape");
    const AreaRuleResult surface = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::positive_lumped,
        0.0,
        "h20_19_surface_to_surface",
        fuelsim::MechanicalContactDiscretization::surface_to_surface);
    constexpr double consistent_tolerance = 1.0e-2;
    bool passed =
        check(positive.positive_area_nodes == 8 && positive.negative_area_nodes == 0,
            "H20.19 positive lumping gives eight positive secondary nodal areas")
        && check(consistent.positive_area_nodes == 4 && consistent.negative_area_nodes == 4,
            "H20.19 consistent Quad8 integration retains four negative corner areas")
        && check(fuelsim::test::relative_metrics_below(consistent.displacement_x, consistent_tolerance)
                     && consistent.normal_force_relative_error < consistent_tolerance,
            "H20.19 consistent-shape displacement and resultant errors are below 1 percent")
        && check(consistent.displacement_x.relative_l2() < positive.displacement_x.relative_l2()
                     && consistent.displacement_x.relative_absolute_peak()
                            < positive.displacement_x.relative_absolute_peak() + 1.0e-15
                     && consistent.displacement_x.maximum_pointwise_relative
                            < positive.displacement_x.maximum_pointwise_relative
                     && consistent.normal_force_relative_error < positive.normal_force_relative_error,
            "H20.19 MOOSE mortar ranks consistent-shape integration no worse in relative peak and closer in "
            "the other field metrics and normal resultant")
        && check(positive.displacement_x.relative_l2() > 5.0e-2 && positive.normal_force_relative_error > 5.0e-2,
            "H20.19 distinguishes the positive-lumped response from the mortar reference")
        && check(fuelsim::test::relative_metrics_below(surface.displacement_x, consistent_tolerance)
                     && surface.normal_force_relative_error < consistent_tolerance,
            "H20.19 surface-to-surface displacement and resultant errors are below 1 percent")
        && positive.passed && consistent.passed && surface.passed;
    return passed;
}

bool run_abaqus_sliding_comparison(const std::string& case_path,
    const std::string& displacement_path,
    const std::string& reaction_path) {
    fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || definition.spatial.contacts[0].thermal
        || !definition.spatial.contacts[0].mechanical)
        throw std::invalid_argument("H20.23 requires isolated HEX20 mechanical contact");
    definition.spatial.contacts[0].mechanical_discretization =
        fuelsim::MechanicalContactDiscretization::surface_to_surface;
    definition.spatial.contacts[0].quad8_nodal_area_rule = fuelsim::Quad8NodalAreaRule::positive_lumped;
    definition.spatial.contacts[0].friction_coefficient = 0.001;
    definition.spatial.contacts[0].friction_slip_tolerance = 1.0e-6;
    const auto secondary_y = std::find_if(definition.spatial.boundary_conditions.begin(),
        definition.spatial.boundary_conditions.end(),
        [](const fuelsim::BoundaryConditionDefinition& boundary) { return boundary.name == "secondary_y"; });
    if (secondary_y == definition.spatial.boundary_conditions.end())
        throw std::invalid_argument("H20.23 secondary_y boundary condition is missing");
    secondary_y->value = 2.8e-6;
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
        "H20.23 surface-to-surface sliding path completes one load step");
    const auto reference = read_displacement(displacement_path);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    std::array<fuelsim::test::FieldErrorMetrics, 3> metrics;
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_difference = 0.0, symmetry_plane_z_maximum_absolute_difference = 0.0;
    std::size_t symmetry_plane_z_count = 0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto found = std::find_if(reference.begin(),
                reference.end(),
                [&](const DisplacementReference& value) { return value.id == source; });
            if (found == reference.end() || present[source])
                throw std::invalid_argument("H20.23 displacement source-node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_difference =
                std::max(maximum_coordinate_difference, coordinate_difference(mesh.nodes()[source], found->point));
            const std::size_t global = spatial.global_node(region, local);
            for (std::size_t component = 0; component < 2; ++component)
                metrics[component].add(solve.solve.state[fields[component + 1].begin + global],
                    found->displacement[component]);
            const double actual_z = solve.solve.state[fields[3].begin + global];
            if (std::abs(found->point.z - 5.0e-3) > 1.0e-8) {
                metrics[2].add(actual_z, found->displacement[2]);
            } else {
                ++symmetry_plane_z_count;
                symmetry_plane_z_maximum_absolute_difference =
                    std::max(symmetry_plane_z_maximum_absolute_difference, std::abs(actual_z - found->displacement[2]));
            }
        }
    }
    constexpr double tolerance = 1.0e-2;
    const std::array<std::string, 3> names = {"displacement_x", "displacement_y", "displacement_z_off_symmetry_plane"};
    for (std::size_t component = 0; component < metrics.size(); ++component) {
        fuelsim::test::print_relative_metrics("h20_23_" + names[component], metrics[component]);
        passed = check(fuelsim::test::relative_metrics_below(metrics[component], tolerance),
                     "H20.23 " + names[component] + " three Abaqus errors are below 1 percent")
                 && passed;
    }
    const ReactionReference reaction = read_reaction(reaction_path);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    const double normal_error = std::abs(interface.total_contact_force - reaction.normal_force) / reaction.normal_force,
                 tangential_error =
                     std::abs(interface.total_tangential_force - reaction.tangential_force) / reaction.tangential_force;
    const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    fuelsim::SteadyProblem numerical_problem(definition.spatial, mesh);
    const auto& numerical_spatial = fuelsim::cartesian::ProblemAccess::view(numerical_problem);
    std::array<double, 3> contact_balance{};
    double maximum_jacobian_directional_error = 0.0;
    for (std::size_t contribution = 0; contribution < numerical_spatial.contribution_count(); ++contribution) {
        if (numerical_spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        numerical_problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local)
            local_state[local] = solve.solve.state[dofs[local]];
        std::vector<double> residual, jacobian;
        numerical_spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local_state.size()), plus = local_state, minus = local_state;
        constexpr double perturbation = 1.0e-10;
        for (std::size_t local = 0; local < local_state.size(); ++local) {
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] += perturbation * direction[local];
            minus[local] -= perturbation * direction[local];
        }
        std::vector<double> plus_residual, minus_residual;
        numerical_spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        numerical_spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local_state.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[row * local_state.size() + column] * direction[column];
            const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - finite_difference, 2);
            reference_squared += finite_difference * finite_difference;
            for (std::size_t component = 0; component < contact_balance.size(); ++component) {
                const auto& field = fields[component + 1];
                if (dofs[row] >= field.begin && dofs[row] < field.end)
                    contact_balance[component] += residual[row];
            }
        }
        maximum_jacobian_directional_error =
            std::max(maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
    }
    passed =
        check(reference.size() == mesh.nodes().size() && maximum_coordinate_difference < 1.0e-9,
            "H20.23 compares all forty normal-displacement nodes on the tracked Abaqus mesh")
        && check(symmetry_plane_z_count == 8 && symmetry_plane_z_maximum_absolute_difference < 1.0e-16,
            "H20.23 checks the theoretical-zero symmetry plane by absolute difference")
        && check(normal_error < tolerance && tangential_error < tolerance,
            "H20.23 normal and tangential resultants agree with Abaqus below 1 percent")
        && check(
            interface.total_tangential_force > 0.0 && histories.size() == 8
                && std::all_of(histories.begin(), histories.end(), [](const auto& history) { return history.sliding; }),
            "H20.23 keeps all eight Abaqus-style averaged constraints in Coulomb sliding")
        && check(std::all_of(contact_balance.begin(),
                     contact_balance.end(),
                     [](double value) { return std::abs(value) < 1.0e-10; }),
            "H20.23 averaged friction residual is action-reaction conservative")
        && check(maximum_jacobian_directional_error < 1.0e-7,
            "H20.23 averaged sliding-friction Jacobian matches a centered directional difference")
        && passed;
    std::cout << "h20_23_normal_resultant=" << interface.total_contact_force << '\n'
              << "h20_23_reference_normal_resultant=" << reaction.normal_force << '\n'
              << "h20_23_normal_resultant_relative_error=" << normal_error << '\n'
              << "h20_23_tangential_resultant=" << interface.total_tangential_force << '\n'
              << "h20_23_reference_tangential_resultant=" << reaction.tangential_force << '\n'
              << "h20_23_tangential_resultant_relative_error=" << tangential_error << '\n'
              << "h20_23_symmetry_plane_z_count=" << symmetry_plane_z_count << '\n'
              << "h20_23_symmetry_plane_z_maximum_absolute_difference=" << symmetry_plane_z_maximum_absolute_difference
              << '\n'
              << "h20_23_contact_jacobian_directional_error=" << maximum_jacobian_directional_error << '\n';
    return passed;
}

bool run_abaqus_area_comparison(const std::string& case_path,
    const std::string& displacement_path,
    const std::string& reaction_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || definition.spatial.contacts[0].thermal
        || !definition.spatial.contacts[0].mechanical || definition.spatial.contacts[0].friction_coefficient != 0.0
        || definition.spatial.contacts[0].quad8_nodal_area_rule != fuelsim::Quad8NodalAreaRule::positive_lumped)
        throw std::invalid_argument("H20.21 requires isolated pure-normal HEX20 mechanical contact");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    const auto displacement = read_displacement(displacement_path);
    const ReactionReference reaction = read_reaction(reaction_path);
    const AreaRuleResult positive = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::positive_lumped,
        0.0,
        "h20_21_positive_lumped");
    const AreaRuleResult consistent = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::consistent_shape,
        0.0,
        "h20_21_consistent_shape");
    const AreaRuleResult surface = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::positive_lumped,
        0.0,
        "h20_21_surface_to_surface",
        fuelsim::MechanicalContactDiscretization::surface_to_surface);
    constexpr double consistent_tolerance = 5.0e-3;
    return check(fuelsim::test::relative_metrics_below(consistent.displacement_x, consistent_tolerance)
                     && consistent.normal_force_relative_error < consistent_tolerance,
               "H20.21 consistent-shape displacement and resultant errors are below 0.5 percent")
           && check(consistent.displacement_x.relative_l2() < positive.displacement_x.relative_l2()
                        && consistent.displacement_x.relative_absolute_peak()
                               <= positive.displacement_x.relative_absolute_peak()
                        && consistent.displacement_x.maximum_pointwise_relative
                               < positive.displacement_x.maximum_pointwise_relative
                        && consistent.normal_force_relative_error < positive.normal_force_relative_error,
               "H20.21 Abaqus surface contact ranks consistent-shape integration no worse in relative peak and closer "
               "in the other field metrics and normal resultant")
           && check(positive.displacement_x.relative_l2() > 5.0e-2 && positive.normal_force_relative_error > 5.0e-2,
               "H20.21 distinguishes the positive-lumped response from the Abaqus reference")
           && check(fuelsim::test::relative_metrics_below(surface.displacement_x, consistent_tolerance)
                        && surface.normal_force_relative_error < consistent_tolerance,
               "H20.21 surface-to-surface displacement and resultant errors are below 0.5 percent")
           && positive.passed && consistent.passed && surface.passed;
}

bool run_unnormalized_area_comparison(const std::string& case_path,
    const std::string& displacement_path,
    const std::string& contact_path,
    const std::string& reaction_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    if (definition.problem != fuelsim::CaseProblem::steady || definition.geometry != fuelsim::CaseGeometry::cartesian_3d
        || definition.spatial.contacts.size() != 1 || definition.spatial.contacts[0].thermal
        || !definition.spatial.contacts[0].mechanical || definition.spatial.contacts[0].friction_coefficient != 0.0
        || definition.spatial.contacts[0].quad8_nodal_area_rule != fuelsim::Quad8NodalAreaRule::positive_lumped)
        throw std::invalid_argument("H20.22 requires isolated pure-normal HEX20 mechanical contact");
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    const auto displacement = read_displacement(displacement_path);
    const auto contact = read_contact(contact_path);
    const ReactionReference reaction = read_reaction(reaction_path);
    std::size_t active_nodes = 0;
    std::size_t positive_area_nodes = 0;
    std::size_t negative_area_nodes = 0;
    for (const ContactReference& node : contact) {
        if (node.pressure > 0.0)
            ++active_nodes;
        if (node.nodal_area > 0.0)
            ++positive_area_nodes;
        else if (node.nodal_area < 0.0)
            ++negative_area_nodes;
    }
    const AreaRuleResult positive = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::positive_lumped,
        0.0,
        "h20_22_positive_lumped");
    const AreaRuleResult consistent = run_area_rule(definition,
        mesh,
        displacement,
        reaction,
        fuelsim::Quad8NodalAreaRule::consistent_shape,
        0.0,
        "h20_22_consistent_shape");
    std::cout << "h20_22_moose_active_contact_nodes=" << active_nodes << '\n'
              << "h20_22_moose_positive_area_nodes=" << positive_area_nodes << '\n'
              << "h20_22_moose_negative_area_nodes=" << negative_area_nodes << '\n';
    return check(contact.size() == 8 && active_nodes == 4,
               "H20.22 unnormalized MOOSE contact leaves only the four edge-midpoint nodes active")
           && check(positive_area_nodes == 4 && negative_area_nodes == 4,
               "H20.22 MOOSE diagnostics retain the signed consistent Quad8 nodal areas")
           && check(positive.displacement_x.relative_l2() < consistent.displacement_x.relative_l2()
                        && positive.displacement_x.maximum_pointwise_relative
                               < consistent.displacement_x.maximum_pointwise_relative
                        && positive.normal_force_relative_error < consistent.normal_force_relative_error,
               "H20.22 unnormalized MOOSE response ranks positive lumping closer than consistent-shape integration")
           && positive.passed && consistent.passed;
}
} // namespace

int main(int argc, char** argv) {
    const bool thermal_only = argc == 4 && std::string(argv[1]) == "--thermal-only";
    if (!thermal_only && argc != 18) {
        std::cerr << "Usage: fuelsim_h20_hex20_contact_reference_tests --thermal-only "
                     "<thermal.fsi> <temperature.csv>\n"
                     "   or: fuelsim_h20_hex20_contact_reference_tests <thermal.fsi> <temperature.csv> "
                     "<mechanical.fsi> <displacement.csv> <contact.csv> <sliding.fsi> <reaction.csv> "
                     "<mortar.fsi> <mortar_displacement.csv> <mortar_reaction.csv> "
                     "<abaqus_friction_displacement.csv> <abaqus_friction_reaction.csv> "
                     "<abaqus_displacement.csv> <abaqus_reaction.csv> "
                     "<unnormalized_displacement.csv> <unnormalized_contact.csv> <unnormalized_reaction.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim HEX20 contact external-reference comparisons\n");
        if (thermal_only) {
            const bool passed = run_thermal(argv[2], argv[3]);
            if (passed && session.rank() == 0)
                std::cout << "[PASS] HEX20 thermal-contact external comparison\n";
            return passed ? 0 : 1;
        }
        // Retain the legacy node-to-surface comparison calls for source-level
        // reference while the production C3D20T branch is disabled.
        const bool passed = run_thermal(argv[1], argv[2]) && run_mechanical(argv[3], argv[4], argv[5])
                            && run_sliding(argv[6], argv[7]) && run_mortar_area_comparison(argv[8], argv[9], argv[10])
                            && run_abaqus_sliding_comparison(argv[8], argv[11], argv[12])
                            && run_abaqus_area_comparison(argv[8], argv[13], argv[14])
                            && run_unnormalized_area_comparison(argv[8], argv[15], argv[16], argv[17]);
        if (passed && session.rank() == 0)
            std::cout << "[PASS] HEX20 contact external-reference comparisons\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] HEX20 contact external-reference comparison raised: " << error.what() << '\n';
        return 1;
    }
}
