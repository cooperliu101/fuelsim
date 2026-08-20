#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::UnstructuredHex8Mesh two_element_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0}) {
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0}) nodes.push_back({x, y, z});
    }
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) { return z * 6 + y * 3 + x; };
    std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0), node(1, 0, 0), node(1, 1, 0), node(0, 1, 0), node(0, 0, 1), node(1, 0, 1), node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0), node(2, 0, 0), node(2, 1, 0), node(1, 1, 0), node(1, 0, 1), node(2, 0, 1), node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    const std::vector<fuelsim::ElementSide> all_faces = {
        {0, 3}, {1, 1}, {0, 0}, {1, 0}, {0, 2}, {1, 2}, {0, 4}, {1, 4}, {0, 5}, {1, 5}};
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), {1, 1}, {{1, "solid"}}, {},
        {{10, "all", all_faces}, {11, "x0", {{0, 3}}}, {12, "y0", {{0, 0}, {1, 0}}}, {13, "z0", {{0, 4}, {1, 4}}},
            {14, "x2", {{1, 1}}}});
}

fuelsim::UnstructuredHex8Mesh two_region_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double origin : {0.0, 2.0}) {
        for (double z : {0.0, 1.0}) {
            for (double y : {0.0, 1.0})
                for (double x : {0.0, 1.0}) nodes.push_back({origin + x, y, z});
        }
    }
    const auto cube = [](std::size_t offset) {
        return fuelsim::Hex8Element{
            {{offset, offset + 1, offset + 3, offset + 2, offset + 4, offset + 5, offset + 7, offset + 6}}};
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), {cube(0), cube(8)}, {1, 2}, {{1, "first"}, {2, "second"}},
        {},
        {{10, "first_all", {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}}, {11, "first_x0", {{0, 3}}},
            {12, "first_y0", {{0, 0}}}, {13, "first_z0", {{0, 4}}},
            {20, "second_all", {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}}, {21, "second_x0", {{1, 3}}},
            {22, "second_y0", {{1, 0}}}, {23, "second_z0", {{1, 4}}}});
}

fuelsim::UnstructuredHex8Mesh shared_two_region_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0}) {
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0}) nodes.push_back({x, y, z});
    }
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) { return z * 6 + y * 3 + x; };
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0), node(1, 0, 0), node(1, 1, 0), node(0, 1, 0), node(0, 0, 1), node(1, 0, 1), node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0), node(2, 0, 0), node(2, 1, 0), node(1, 1, 0), node(1, 0, 1), node(2, 0, 1), node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), elements, {1, 2}, {{1, "meat"}, {2, "clad"}}, {},
        {{10, "meat_x0", {{0, 3}}}, {11, "clad_x2", {{1, 1}}}, {12, "meat_y0", {{0, 0}}}, {13, "meat_z0", {{0, 4}}},
            {14, "clad_y0", {{1, 0}}}, {15, "clad_z0", {{1, 4}}}});
}

fuelsim::UnstructuredHex8Mesh contact_projection_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0})
        for (double y : {-1.0, 0.0, 1.0})
            for (double x : {0.0, 1.0}) nodes.push_back({x, y, z});
    const auto primary_node = [](std::size_t x, std::size_t y, std::size_t z) { return z * 6 + y * 2 + x; };
    std::vector<fuelsim::Hex8Element> elements = {
        {{{primary_node(0, 0, 0), primary_node(1, 0, 0), primary_node(1, 1, 0), primary_node(0, 1, 0),
            primary_node(0, 0, 1), primary_node(1, 0, 1), primary_node(1, 1, 1), primary_node(0, 1, 1)}}},
        {{{primary_node(0, 1, 0), primary_node(1, 1, 0), primary_node(1, 2, 0), primary_node(0, 2, 0),
            primary_node(0, 1, 1), primary_node(1, 1, 1), primary_node(1, 2, 1), primary_node(0, 2, 1)}}},
    };
    const std::size_t secondary_offset = nodes.size();
    for (double z : {0.25, 0.75})
        for (double y : {-0.25, 0.25})
            for (double x : {1.0, 2.0}) nodes.push_back({x, y, z});
    elements.push_back({{{secondary_offset, secondary_offset + 1, secondary_offset + 3, secondary_offset + 2,
        secondary_offset + 4, secondary_offset + 5, secondary_offset + 7, secondary_offset + 6}}});
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), {1, 1, 2},
        {{1, "primary"}, {2, "secondary"}}, {},
        {{10, "primary_right", {{0, 1}, {1, 1}}}, {20, "secondary_left", {{2, 3}}}});
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.25, 1.0e-5, 300.0, 0.0, 0.0, 0.0, 6000.0, 1000.0);
}

fuelsim::SpatialDefinition steady_definition() {
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", material(), 0.0, 300.0});
    definition.boundary_conditions = {
        {"temperature", fuelsim::BoundaryConditionType::dirichlet, "all", fuelsim::Field::temperature, 400.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
    };
    return definition;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = 1.0e-9;
    options.relative_tolerance = 1.0e-11;
    options.maximum_iterations = 12;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::gmres;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::field_split;
    options.field_residual_scaling = true;
    options.linear_relative_tolerance = 1.0e-11;
    options.maximum_linear_iterations = 200;
    return options;
}

bool check_uniform_solution(const fuelsim::spatial_detail::SpatialLayout& dofs,
    const fuelsim::UnstructuredHex8Mesh& mesh, const std::vector<double>& state, double temperature) {
    const double thermal_strain = 1.0e-5 * (temperature - 300.0);
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node)
        if (!check(std::abs(state[dofs.dof(fuelsim::Field::temperature, node)] - temperature) < 2.0e-9,
                "three-dimensional temperature matches the uniform analytic solution") ||
            !check(std::abs(state[dofs.dof(fuelsim::Field::displacement_x, node)] -
                            thermal_strain * mesh.nodes()[node].x) < 2.0e-11 &&
                       std::abs(state[dofs.dof(fuelsim::Field::displacement_y, node)] -
                                thermal_strain * mesh.nodes()[node].y) < 2.0e-11 &&
                       std::abs(state[dofs.dof(fuelsim::Field::displacement_z, node)] -
                                thermal_strain * mesh.nodes()[node].z) < 2.0e-11,
                "three-dimensional free thermal expansion matches the analytic solution"))
            return false;
    return true;
}

bool test_steady(
    const fuelsim::PetscSession& session, const fuelsim::UnstructuredHex8Mesh& mesh, const std::string& results_path) {
    fuelsim::SteadyProblem problem(steady_definition(), mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    bool passed =
        check(result.completed && result.solve.converged, "three-dimensional steady solve converges") &&
        check(result.solve.field_names ==
                  std::vector<std::string>{"temperature", "displacement_x", "displacement_y", "displacement_z"},
            "three-dimensional solve reports four field-major fields") &&
        check(problem.contribution_count() == 2, "two HEX8 elements expose two independent contributions") &&
        check_uniform_solution(fuelsim::cartesian::ProblemAccess::dof_map(problem), mesh, result.solve.state, 400.0);
    const auto expected_partition = problem.contribution_partition(
        static_cast<std::size_t>(session.rank()), static_cast<std::size_t>(session.size()));
    const std::size_t expected_begin = expected_partition.first;
    const std::size_t expected_end = expected_partition.second;
    passed = check(result.solve.local_contribution_begin == expected_begin &&
                       result.solve.local_contribution_end == expected_end,
                 "each message-passing rank owns only its exact HEX8 contribution interval") &&
             passed;
    session.collective_root_action(
        [&]() { fuelsim::write_steady_results(results_path, mesh, problem, result.solve.state); });
    return passed;
}

bool test_convection_boundary(const fuelsim::UnstructuredHex8Mesh& mesh) {
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", material(), 0.0, 300.0});
    fuelsim::BoundaryConditionDefinition convection{
        "right_coolant", fuelsim::BoundaryConditionType::convection, "x2", fuelsim::Field::temperature, 0.0};
    convection.heat_transfer_coefficient = 10.0;
    convection.ambient_temperature = 400.0;
    definition.boundary_conditions = {
        {"left_temperature", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::temperature, 300.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
    };
    definition.boundary_conditions.push_back(convection);
    fuelsim::SteadyProblem problem(definition, mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    bool passed = check(result.completed && result.solve.converged,
        "three-dimensional input-style convection boundary solve converges");
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t node = 0; passed && node < mesh.nodes().size(); ++node) {
        const double x = mesh.nodes()[node].x;
        const double expected = 300.0 + (100.0 / 3.0) * x;
        passed = check(std::abs(result.solve.state[dofs.dof(fuelsim::Field::temperature, node)] - expected) < 2.0e-8,
                     "three-dimensional convection boundary matches the one-dimensional conduction solution") &&
                 passed;
    }
    return passed;
}

bool test_transient(const fuelsim::PetscSession& session, const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::string& checkpoint_path, const std::string& results_path) {
    fuelsim::SpatialDefinition spatial = steady_definition();
    spatial.boundary_conditions.erase(spatial.boundary_conditions.begin());
    spatial.regions[0].volumetric_heat_source = 6.0e6;
    const fuelsim::SpatialDefinition definition = spatial;
    fuelsim::TransientProblem problem(definition, mesh);
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, {1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 2, 0.0}, solver_options());
    bool passed = check(result.completed && result.accepted_steps.size() == 1,
                      "three-dimensional Backward Euler transient solve accepts one physical time step") &&
                  check_uniform_solution(
                      fuelsim::cartesian::ProblemAccess::dof_map(problem), mesh, problem.committed_solution(), 301.0);
    session.collective_root_action([&]() {
        fuelsim::write_transient_checkpoint(checkpoint_path, problem, 0.25);
        fuelsim::ExodusTransientResultsWriter writer(results_path, mesh, problem);
        writer.append(problem);
    });
    fuelsim::TransientProblem restored(definition, mesh);
    const double next_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    passed = check(next_time_step == 0.25 && restored.committed_time() == 1.0 &&
                       restored.committed_solution() == problem.committed_solution(),
                 "three-dimensional checkpoint restores the exact committed nodal state and controller step") &&
             passed;
    bool stresses_match = true;
    for (std::size_t region = 0; stresses_match && region < definition.regions.size(); ++region) {
        for (std::size_t element = 0;
            stresses_match &&
            element < fuelsim::cartesian::ProblemAccess::region_mesh(problem, region).elements().size();
            ++element) {
            const auto before = fuelsim::cartesian::ProblemAccess::stress(problem, region, element);
            const auto after = fuelsim::cartesian::ProblemAccess::stress(restored, region, element);
            for (std::size_t q = 0; stresses_match && q < 8; ++q) {
                const fuelsim::SymmetricTensor3Values& first = before[q];
                const fuelsim::SymmetricTensor3Values& second = after[q];
                stresses_match = first.xx == second.xx && first.yy == second.yy && first.zz == second.zz &&
                                 first.xy == second.xy && first.yz == second.yz && first.xz == second.xz;
            }
        }
    }
    passed = check(stresses_match,
                 "three-dimensional checkpoint restores all eight integration points and six stress components") &&
             passed;
    return passed;
}

bool test_multiple_regions() {
    const fuelsim::UnstructuredHex8Mesh mesh = two_region_mesh();
    fuelsim::SpatialDefinition definition;
    definition.regions = {{"first", "first", material(), 0.0, 300.0}, {"second", "second", material(), 0.0, 300.0}};
    for (const std::string& prefix : {std::string("first"), std::string("second")}) {
        definition.boundary_conditions.push_back({prefix + "_temperature", fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_all", fuelsim::Field::temperature, 325.0});
        definition.boundary_conditions.push_back({prefix + "_x", fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_x0", fuelsim::Field::displacement_x, 0.0});
        definition.boundary_conditions.push_back({prefix + "_y", fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_y0", fuelsim::Field::displacement_y, 0.0});
        definition.boundary_conditions.push_back({prefix + "_z", fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_z0", fuelsim::Field::displacement_z, 0.0});
    }
    fuelsim::SteadyProblem problem(definition, mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    const auto& first = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 0);
    const auto& second = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 1);
    bool passed =
        check(result.completed && result.solve.converged, "two-region three-dimensional solve converges") &&
        check(first.nodes().size() == 8 && second.nodes().size() == 8 &&
                  fuelsim::cartesian::ProblemAccess::region_node_offset(problem, 1) == 8 && problem.dof_count() == 64,
            "three-dimensional regions keep independent nodes and exact four-field offsets");
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    const double strain = 1.0e-5 * 25.0;
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const double origin = node < 8 ? 0.0 : 2.0;
        const auto& point = mesh.nodes()[node];
        passed =
            check(std::abs(result.solve.state[dofs.dof(fuelsim::Field::temperature, node)] - 325.0) < 2.0e-9 &&
                      std::abs(result.solve.state[dofs.dof(fuelsim::Field::displacement_x, node)] -
                               strain * (point.x - origin)) < 2.0e-11 &&
                      std::abs(result.solve.state[dofs.dof(fuelsim::Field::displacement_y, node)] - strain * point.y) <
                          2.0e-11 &&
                      std::abs(result.solve.state[dofs.dof(fuelsim::Field::displacement_z, node)] - strain * point.z) <
                          2.0e-11,
                "each independent three-dimensional region matches free thermal expansion") &&
            passed;
    }
    return passed;
}

bool test_shared_nodes(const fuelsim::PetscSession& session) {
    const fuelsim::UnstructuredHex8Mesh mesh = shared_two_region_mesh();
    const std::filesystem::path input = std::filesystem::temp_directory_path() / "fuelsim_shared_hex8_input.e";
    try {
        session.collective_root_action([&]() { fuelsim::write_exodus_hex8(input.string(), mesh); });
    } catch (const std::exception& error) {
        return check(false, std::string("shared-node Exodus input write succeeds: ") + error.what());
    }
    const fuelsim::UnstructuredHex8Mesh exodus_mesh = fuelsim::read_exodus_hex8(input.string());
    fuelsim::SpatialDefinition definition;
    definition.regions = {{"meat", "meat", material(), 0.0, 300.0}, {"clad", "clad", material(), 0.0, 300.0}};
    definition.boundary_conditions = {
        {"temperature_left", fuelsim::BoundaryConditionType::dirichlet, "meat_x0", fuelsim::Field::temperature, 325.0},
        {"temperature_right", fuelsim::BoundaryConditionType::dirichlet, "clad_x2", fuelsim::Field::temperature, 325.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "meat_x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "meat_y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "meat_z0", fuelsim::Field::displacement_z, 0.0},
    };
    fuelsim::SteadyProblem problem(definition, exodus_mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    const auto& first = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 0);
    const auto& second = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 1);
    const auto first_shared = std::find(first.source_node_ids().begin(), first.source_node_ids().end(), 1);
    const auto second_shared = std::find(second.source_node_ids().begin(), second.source_node_ids().end(), 1);
    bool passed =
        check(result.completed && result.solve.converged, "shared-node three-dimensional solve converges") &&
        check(problem.dof_count() == 48 && first.nodes().size() == 8 && second.nodes().size() == 8,
            "shared source nodes reduce the global four-field system to twelve unique nodes") &&
        check(first_shared != first.source_node_ids().end() && second_shared != second.source_node_ids().end(),
            "shared source node is present in both material regions");
    if (passed) {
        const std::size_t first_local = static_cast<std::size_t>(first_shared - first.source_node_ids().begin());
        const std::size_t second_local = static_cast<std::size_t>(second_shared - second.source_node_ids().begin());
        passed = check(fuelsim::cartesian::ProblemAccess::dof_map(problem).global_node(0, first_local) ==
                           fuelsim::cartesian::ProblemAccess::dof_map(problem).global_node(1, second_local),
            "shared source node maps to one global node across material regions");
    }
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t node = 0; passed && node < exodus_mesh.nodes().size(); ++node)
        passed = check(std::abs(result.solve.state[dofs.dof(fuelsim::Field::temperature, node)] - 325.0) < 2.0e-9,
            "shared-node solve preserves the continuous uniform temperature");
    const std::filesystem::path output = std::filesystem::temp_directory_path() / "fuelsim_shared_hex8.e";
    try {
        session.collective_root_action(
            [&]() { fuelsim::write_steady_results(output.string(), exodus_mesh, problem, result.solve.state); });
    } catch (const std::exception& error) {
        passed = check(false, std::string("shared-node Exodus output succeeds: ") + error.what());
    }
    session.collective_root_action([&]() {
        (void)std::remove(input.c_str());
        (void)std::remove(output.c_str());
    });
    return passed;
}

bool test_contact_projection_transfer() {
    const fuelsim::UnstructuredHex8Mesh mesh = contact_projection_mesh();
    fuelsim::SpatialDefinition definition;
    definition.regions = {
        {"primary", "primary", material(), 0.0, 300.0}, {"secondary", "secondary", material(), 0.0, 400.0}};
    definition.contacts.push_back({"interface", "primary_right", "secondary_left", true, true, 1.0, 0.01, 1.0e6, 0.1});
    fuelsim::SteadyProblem problem(definition, mesh);
    const std::vector<std::size_t> secondary_sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    const std::size_t source = 12;
    const auto source_position = std::find(secondary_sources.begin(), secondary_sources.end(), source);
    if (source_position == secondary_sources.end())
        throw std::logic_error("Three-dimensional contact test could not locate its secondary node");
    const std::size_t summary_index = static_cast<std::size_t>(source_position - secondary_sources.begin());
    const auto& secondary_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 1);
    const auto local_position =
        std::find(secondary_mesh.source_node_ids().begin(), secondary_mesh.source_node_ids().end(), source);
    if (local_position == secondary_mesh.source_node_ids().end())
        throw std::logic_error("Three-dimensional contact test could not map its secondary node");
    const std::size_t global_node = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, 1) +
                                    static_cast<std::size_t>(local_position - secondary_mesh.source_node_ids().begin());
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    const auto transferred = [&](double current_y) {
        std::vector<double> state = problem.initial_state();
        state[dofs.dof(fuelsim::Field::displacement_x, global_node)] = -1.0e-4;
        state[dofs.dof(fuelsim::Field::displacement_y, global_node)] = current_y - mesh.nodes()[source].y;
        problem.validate_state(state);
        return fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state).at(summary_index);
    };
    constexpr double offset = 1.0e-9;
    const fuelsim::CartesianContactNodeSummary before = transferred(-offset);
    const fuelsim::CartesianContactNodeSummary at = transferred(0.0);
    const fuelsim::CartesianContactNodeSummary after = transferred(offset);
    const double force_scale = std::max({1.0, std::abs(before.contact_force), std::abs(after.contact_force)});
    bool passed =
        check(before.projected && at.projected && after.projected && before.primary_face != after.primary_face &&
                  (at.primary_face == before.primary_face || at.primary_face == after.primary_face),
            "three-dimensional internal primary edge has one owner and transfers ownership once") &&
        check(std::abs(before.contact_force - after.contact_force) < 1.0e-8 * force_scale,
            "three-dimensional contact force is continuous across primary-face ownership transfer");
    std::vector<double> lost = problem.initial_state();
    for (std::size_t source_node : secondary_sources) {
        const auto local =
            std::find(secondary_mesh.source_node_ids().begin(), secondary_mesh.source_node_ids().end(), source_node);
        const std::size_t node = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, 1) +
                                 static_cast<std::size_t>(local - secondary_mesh.source_node_ids().begin());
        lost[dofs.dof(fuelsim::Field::displacement_y, node)] = 2.0;
    }
    bool rejected = false;
    try {
        problem.validate_state(lost);
    } catch (const std::domain_error&) { rejected = true; }
    return check(rejected, "three-dimensional contact rejects a thermal point or mechanical node outside the complete "
                           "primary surface") &&
           passed;
}

bool test_pressure_configuration_selection(const fuelsim::UnstructuredHex8Mesh& mesh) {
    const auto warning_count = [&](bool finite_strain, bool current_configuration) {
        fuelsim::SpatialDefinition definition;
        definition.regions.push_back({"solid", "solid", material(), 0.0, 300.0});
        definition.regions.front().strain_formulation =
            finite_strain ? fuelsim::StrainFormulation::finite : fuelsim::StrainFormulation::small;
        fuelsim::BoundaryConditionDefinition pressure{
            "pressure", fuelsim::BoundaryConditionType::pressure, "x2", fuelsim::Field::displacement_x, 1.0e6};
        pressure.use_displaced_geometry = current_configuration;
        definition.boundary_conditions.push_back(pressure);
        fuelsim::SteadyProblem problem(std::move(definition), mesh);
        return fuelsim::cartesian::ProblemAccess::dof_map(problem).configuration_warnings().size();
    };
    const std::size_t reference_small = warning_count(false, false);
    const std::size_t current_small = warning_count(false, true);
    const std::size_t reference_finite = warning_count(true, false);
    const std::size_t current_finite = warning_count(true, true);
    return check(reference_small == 0 && current_small == 1 && reference_finite == 1 && current_finite == 0,
        "three-dimensional pressure accepts both configurations and warns for non-recommended choices");
}

fuelsim::SpatialDefinition inelastic_definition(bool creep, bool plasticity) {
    fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 10.0, 2.0e11, 0.3, 0.0, 600.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    if (creep) properties = fuelsim::test::with_norton(std::move(properties), 1.0e-4, 1.0e8, 3.0);
    if (plasticity) properties = fuelsim::test::with_plasticity(std::move(properties), 2.0e8, 2.0e9);
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", std::move(properties), 0.0, 600.0});
    definition.boundary_conditions = {
        {"temperature", fuelsim::BoundaryConditionType::dirichlet, "all", fuelsim::Field::temperature, 600.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
        {"traction", fuelsim::BoundaryConditionType::traction, "x2", fuelsim::Field::displacement_x, 2.01e8, true},
    };
    return definition;
}

bool test_inelastic_branches(const fuelsim::UnstructuredHex8Mesh& mesh) {
    bool passed = true;
    for (const std::array<bool, 2> branch : {std::array<bool, 2>{false, true}, {true, false}, {true, true}}) {
        fuelsim::TransientProblem problem(inelastic_definition(branch[0], branch[1]), mesh);
        fuelsim::SolverOptions options = solver_options();
        options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
        options.maximum_iterations = 30;
        const fuelsim::TransientResult result =
            fuelsim::solve_transient(problem, {1.0, 0.1, 0.1, 0.1, 1.0, 0.5, 0, 1.0}, options);
        double maximum_plastic = 0.0, maximum_creep = 0.0;
        bool finite = true;
        for (std::size_t element = 0; element < 2; ++element)
            for (const fuelsim::CartesianMaterialPointState& point :
                fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element)) {
                maximum_plastic = std::max(maximum_plastic, point.equivalent_plastic_strain);
                maximum_creep = std::max(maximum_creep, point.equivalent_creep_strain);
                finite = finite && std::isfinite(point.stress.xx) && std::isfinite(point.stress.yy) &&
                         std::isfinite(point.stress.zz) && std::isfinite(point.stress.xy) &&
                         std::isfinite(point.stress.yz) && std::isfinite(point.stress.xz);
            }
        passed = check(result.completed && result.accepted_steps.size() == 10 && finite &&
                           (branch[1] ? maximum_plastic > 0.0 : maximum_plastic == 0.0) &&
                           (branch[0] ? maximum_creep > 0.0 : maximum_creep == 0.0),
                     "three-dimensional plastic, creep, and coupled transient branches solve and commit once") &&
                 passed;
    }
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_hex8_solver_tests <steady.e> <transient.e> <checkpoint.bin>\n";
        return 2;
    }
    fuelsim::PetscSession session(argc, argv, "fuelsim HEX8 solver tests\n");
    const fuelsim::UnstructuredHex8Mesh mesh = two_element_mesh();
    bool passed = test_steady(session, mesh, argv[1]);
    passed = test_convection_boundary(mesh) && passed;
    passed = test_transient(session, mesh, argv[3], argv[2]) && passed;
    passed = test_multiple_regions() && passed;
    passed = test_shared_nodes(session) && passed;
    passed = test_contact_projection_transfer() && passed;
    passed = test_pressure_configuration_selection(mesh) && passed;
    passed = test_inelastic_branches(mesh) && passed;
    session.collective_root_action([&]() {
        (void)std::remove(argv[1]);
        (void)std::remove(argv[2]);
        (void)std::remove(argv[3]);
    });
    if (passed && session.rank() == 0)
        std::cout << "HEX8 steady, transient, message-passing, results, and checkpoint tests passed\n";
    return passed ? 0 : 1;
}
