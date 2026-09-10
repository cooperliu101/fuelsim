#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <map>
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

bool same_material_point(const fuelsim::CartesianMaterialPointState& first,
    const fuelsim::CartesianMaterialPointState& second) {
    return first.elastic_strain == second.elastic_strain && first.plastic_strain == second.plastic_strain
           && first.creep_strain == second.creep_strain
           && first.equivalent_plastic_strain == second.equivalent_plastic_strain
           && first.equivalent_creep_strain == second.equivalent_creep_strain && first.stress.xx == second.stress.xx
           && first.stress.yy == second.stress.yy && first.stress.zz == second.stress.zz
           && first.stress.xy == second.stress.xy && first.stress.yz == second.stress.yz
           && first.stress.xz == second.stress.xz;
}

fuelsim::UnstructuredHex8Mesh two_element_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0}) {
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0})
                nodes.push_back({x, y, z});
    }
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) {
        return z * 6 + y * 3 + x;
    };
    std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0),
            node(1, 0, 0),
            node(1, 1, 0),
            node(0, 1, 0),
            node(0, 0, 1),
            node(1, 0, 1),
            node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0),
            node(2, 0, 0),
            node(2, 1, 0),
            node(1, 1, 0),
            node(1, 0, 1),
            node(2, 0, 1),
            node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    const std::vector<fuelsim::ElementSide> all_faces =
        {{0, 3}, {1, 1}, {0, 0}, {1, 0}, {0, 2}, {1, 2}, {0, 4}, {1, 4}, {0, 5}, {1, 5}};
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        std::move(elements),
        {1, 1},
        {{1, "solid"}},
        {},
        {{10, "all", all_faces},
            {11, "x0", {{0, 3}}},
            {12, "y0", {{0, 0}, {1, 0}}},
            {13, "z0", {{0, 4}, {1, 4}}},
            {14, "x2", {{1, 1}}}});
}

fuelsim::UnstructuredHex8Mesh two_region_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double origin : {0.0, 2.0}) {
        for (double z : {0.0, 1.0}) {
            for (double y : {0.0, 1.0})
                for (double x : {0.0, 1.0})
                    nodes.push_back({origin + x, y, z});
        }
    }
    const auto cube = [](std::size_t offset) {
        return fuelsim::Hex8Element{
            {{offset, offset + 1, offset + 3, offset + 2, offset + 4, offset + 5, offset + 7, offset + 6}}};
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        {cube(0), cube(8)},
        {1, 2},
        {{1, "first"}, {2, "second"}},
        {},
        {{10, "first_all", {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}},
            {11, "first_x0", {{0, 3}}},
            {12, "first_y0", {{0, 0}}},
            {13, "first_z0", {{0, 4}}},
            {20, "second_all", {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}},
            {21, "second_x0", {{1, 3}}},
            {22, "second_y0", {{1, 0}}},
            {23, "second_z0", {{1, 4}}}});
}

fuelsim::UnstructuredHex8Mesh shared_two_region_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0}) {
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0})
                nodes.push_back({x, y, z});
    }
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) {
        return z * 6 + y * 3 + x;
    };
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0),
            node(1, 0, 0),
            node(1, 1, 0),
            node(0, 1, 0),
            node(0, 0, 1),
            node(1, 0, 1),
            node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0),
            node(2, 0, 0),
            node(2, 1, 0),
            node(1, 1, 0),
            node(1, 0, 1),
            node(2, 0, 1),
            node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        elements,
        {1, 2},
        {{1, "meat"}, {2, "clad"}},
        {},
        {{10, "meat_x0", {{0, 3}}},
            {11, "clad_x2", {{1, 1}}},
            {12, "meat_y0", {{0, 0}}},
            {13, "meat_z0", {{0, 4}}},
            {14, "clad_y0", {{1, 0}}},
            {15, "clad_z0", {{1, 4}}}});
}

fuelsim::UnstructuredHex8Mesh contact_projection_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0})
        for (double y : {-1.0, 0.0, 1.0})
            for (double x : {0.0, 1.0})
                nodes.push_back({x, y, z});
    const auto primary_node = [](std::size_t x, std::size_t y, std::size_t z) {
        return z * 6 + y * 2 + x;
    };
    std::vector<fuelsim::Hex8Element> elements = {
        {{{primary_node(0, 0, 0),
            primary_node(1, 0, 0),
            primary_node(1, 1, 0),
            primary_node(0, 1, 0),
            primary_node(0, 0, 1),
            primary_node(1, 0, 1),
            primary_node(1, 1, 1),
            primary_node(0, 1, 1)}}},
        {{{primary_node(0, 1, 0),
            primary_node(1, 1, 0),
            primary_node(1, 2, 0),
            primary_node(0, 2, 0),
            primary_node(0, 1, 1),
            primary_node(1, 1, 1),
            primary_node(1, 2, 1),
            primary_node(0, 2, 1)}}},
    };
    const std::size_t secondary_offset = nodes.size();
    for (double z : {0.25, 0.75})
        for (double y : {-0.25, 0.25})
            for (double x : {1.0, 2.0})
                nodes.push_back({x, y, z});
    elements.push_back({{{secondary_offset,
        secondary_offset + 1,
        secondary_offset + 3,
        secondary_offset + 2,
        secondary_offset + 4,
        secondary_offset + 5,
        secondary_offset + 7,
        secondary_offset + 6}}});
    std::vector<std::size_t> primary_all(secondary_offset), secondary_all(8);
    for (std::size_t node = 0; node < primary_all.size(); ++node)
        primary_all[node] = node;
    for (std::size_t node = 0; node < secondary_all.size(); ++node)
        secondary_all[node] = secondary_offset + node;
    std::vector<fuelsim::ElementSide> primary_all_faces, secondary_all_faces;
    for (std::size_t element = 0; element < 2; ++element)
        for (std::size_t side = 0; side < 6; ++side)
            primary_all_faces.push_back({element, side});
    for (std::size_t side = 0; side < 6; ++side)
        secondary_all_faces.push_back({2, side});
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        std::move(elements),
        {1, 1, 2},
        {{1, "primary"}, {2, "secondary"}},
        {{30, "primary_all", std::move(primary_all)}, {40, "secondary_all", std::move(secondary_all)}},
        {{10, "primary_right", {{0, 1}, {1, 1}}},
            {20, "secondary_left", {{2, 3}}},
            {50, "primary_all", std::move(primary_all_faces)},
            {60, "secondary_all", std::move(secondary_all_faces)}});
}

fuelsim::Hex8Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map,
    double x0,
    double x1,
    double y0,
    double y1,
    double z0,
    double z1) {
    const std::array<fuelsim::CartesianPoint3, 8> points = {{{x0, y0, z0},
        {x1, y0, z0},
        {x1, y1, z0},
        {x0, y1, z0},
        {x0, y0, z1},
        {x1, y0, z1},
        {x1, y1, z1},
        {x0, y1, z1}}};
    fuelsim::Hex8Element element{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const std::array<double, 3> key = {points[local].x, points[local].y, points[local].z};
        const auto inserted = node_map.emplace(key, nodes.size());
        if (inserted.second)
            nodes.push_back(points[local]);
        element.nodes[local] = inserted.first->second;
    }
    return element;
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

double mechanical_contact_directional_error(fuelsim::SteadyProblem& problem,
    const std::vector<double>& global_state,
    double perturbation) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum_error = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local)
            local_state[local] = global_state[dofs[local]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local_state.size()), plus = local_state, minus = local_state;
        for (std::size_t local = 0; local < local_state.size(); ++local) {
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] += perturbation * direction[local];
            minus[local] -= perturbation * direction[local];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local_state.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[row * local_state.size() + column] * direction[column];
            const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - finite_difference, 2);
            reference_squared += finite_difference * finite_difference;
        }
        maximum_error = std::max(maximum_error, std::sqrt(difference_squared / reference_squared));
    }
    return maximum_error;
}

bool check_uniform_solution(const fuelsim::spatial_detail::SpatialLayout& dofs,
    const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::vector<double>& state,
    double temperature) {
    const double thermal_strain = 1.0e-5 * (temperature - 300.0);
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node)
        if (!check(std::abs(state[dofs.dof(fuelsim::Field::temperature, node)] - temperature) < 2.0e-9,
                "three-dimensional temperature matches the uniform analytic solution")
            || !check(
                std::abs(state[dofs.dof(fuelsim::Field::displacement_x, node)] - thermal_strain * mesh.nodes()[node].x)
                        < 2.0e-11
                    && std::abs(state[dofs.dof(fuelsim::Field::displacement_y, node)]
                                - thermal_strain * mesh.nodes()[node].y)
                           < 2.0e-11
                    && std::abs(state[dofs.dof(fuelsim::Field::displacement_z, node)]
                                - thermal_strain * mesh.nodes()[node].z)
                           < 2.0e-11,
                "three-dimensional free thermal expansion matches the analytic solution"))
            return false;
    return true;
}

bool test_steady(const fuelsim::PetscSession& session,
    const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::string& results_path) {
    fuelsim::SteadyProblem problem(steady_definition(), mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    bool passed =
        check(result.completed && result.solve.converged, "three-dimensional steady solve converges")
        && check(result.solve.field_names
                     == std::vector<std::string>{"temperature", "displacement_x", "displacement_y", "displacement_z"},
            "three-dimensional solve reports four field-major fields")
        && check(problem.contribution_count() == 2, "two HEX8 elements expose two independent contributions")
        && check_uniform_solution(fuelsim::cartesian::ProblemAccess::dof_map(problem), mesh, result.solve.state, 400.0);
    const auto expected_partition = problem.contribution_partition(static_cast<std::size_t>(session.rank()),
        static_cast<std::size_t>(session.size()));
    const std::size_t expected_begin = expected_partition.first;
    const std::size_t expected_end = expected_partition.second;
    passed = check(result.solve.local_contribution_begin == expected_begin
                       && result.solve.local_contribution_end == expected_end,
                 "each message-passing rank owns only its exact HEX8 contribution interval")
             && passed;
    session.collective_root_action(
        [&]() { fuelsim::write_steady_results(results_path, mesh, problem, result.solve.state); });
    return passed;
}

bool test_small_strain_steady_predictor(const fuelsim::UnstructuredHex8Mesh& mesh) {
    fuelsim::SpatialDefinition definition = steady_definition();
    definition.regions[0].strain_formulation = fuelsim::StrainFormulation::finite;
    definition.regions[0].hex8_element_formulation = fuelsim::Hex8ElementFormulation::c3d8rt;
    fuelsim::SteadyProblem problem(std::move(definition), mesh);
    fuelsim::SteadyLoadOptions load_options{1, 0.5, 4, 1.0e-6, true};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, load_options, solver_options());
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    return check(result.completed && result.solve.converged && result.used_small_strain_predictor
                     && result.predictor_nonlinear_iterations > 0 && result.solve.nonlinear_iterations > 0
                     && result.total_nonlinear_iterations
                            == result.predictor_nonlinear_iterations + result.solve.nonlinear_iterations
                     && result.aggregate_timing.workspace_setups == 1 && result.predictor_timing.workspace_setups == 1
                     && result.solve.timing.workspace_setups == 0
                     && spatial.region(0).strain_formulation == fuelsim::StrainFormulation::finite,
        "small-strain steady predictor reuses one PETSc workspace and restores the finite-strain formulation");
}

bool test_convection_boundary(const fuelsim::UnstructuredHex8Mesh& mesh) {
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", material(), 0.0, 300.0});
    fuelsim::BoundaryConditionDefinition convection{"right_coolant",
        fuelsim::BoundaryConditionType::convection,
        "x2",
        fuelsim::Field::temperature,
        0.0};
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
                     "three-dimensional convection boundary matches the one-dimensional conduction solution")
                 && passed;
    }
    return passed;
}

bool test_transient(const fuelsim::PetscSession& session,
    const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::string& checkpoint_path,
    const std::string& results_path) {
    fuelsim::SpatialDefinition spatial = steady_definition();
    spatial.boundary_conditions.erase(spatial.boundary_conditions.begin());
    spatial.regions[0].volumetric_heat_source = 6.0e6;
    const fuelsim::SpatialDefinition definition = spatial;
    fuelsim::TransientProblem problem(definition, mesh);
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, {1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 2, 0.0}, solver_options());
    bool passed = check(result.completed && result.accepted_steps.size() == 1,
                      "three-dimensional Backward Euler transient solve accepts one physical time step")
                  && check_uniform_solution(fuelsim::cartesian::ProblemAccess::dof_map(problem),
                      mesh,
                      problem.committed_solution(),
                      301.0);
    session.collective_root_action([&]() {
        fuelsim::write_transient_checkpoint(checkpoint_path, problem, 0.25);
        fuelsim::ExodusTransientResultsWriter writer(results_path, mesh, problem);
        writer.append(problem);
    });
    fuelsim::TransientProblem restored(definition, mesh);
    const double next_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    passed = check(next_time_step == 0.25 && restored.committed_time() == 1.0
                       && restored.committed_solution() == problem.committed_solution(),
                 "three-dimensional checkpoint restores the exact committed nodal state and controller step")
             && passed;
    bool stresses_match = true;
    for (std::size_t region = 0; stresses_match && region < definition.regions.size(); ++region) {
        for (std::size_t element = 0;
            stresses_match
            && element < fuelsim::cartesian::ProblemAccess::region_mesh(problem, region).elements().size();
            ++element) {
            const auto before = fuelsim::cartesian::ProblemAccess::stress(problem, region, element);
            const auto after = fuelsim::cartesian::ProblemAccess::stress(restored, region, element);
            for (std::size_t q = 0; stresses_match && q < 8; ++q) {
                const fuelsim::SymmetricTensor3Values& first = before[q];
                const fuelsim::SymmetricTensor3Values& second = after[q];
                stresses_match = first.xx == second.xx && first.yy == second.yy && first.zz == second.zz
                                 && first.xy == second.xy && first.yz == second.yz && first.xz == second.xz;
            }
        }
    }
    passed = check(stresses_match,
                 "three-dimensional checkpoint restores all eight integration points and six stress components")
             && passed;
    return passed;
}

bool test_multiple_regions() {
    const fuelsim::UnstructuredHex8Mesh mesh = two_region_mesh();
    fuelsim::SpatialDefinition definition;
    definition.regions = {{"first", "first", material(), 0.0, 300.0}, {"second", "second", material(), 0.0, 300.0}};
    for (const std::string& prefix : {std::string("first"), std::string("second")}) {
        definition.boundary_conditions.push_back({prefix + "_temperature",
            fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_all",
            fuelsim::Field::temperature,
            325.0});
        definition.boundary_conditions.push_back({prefix + "_x",
            fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_x0",
            fuelsim::Field::displacement_x,
            0.0});
        definition.boundary_conditions.push_back({prefix + "_y",
            fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_y0",
            fuelsim::Field::displacement_y,
            0.0});
        definition.boundary_conditions.push_back({prefix + "_z",
            fuelsim::BoundaryConditionType::dirichlet,
            prefix + "_z0",
            fuelsim::Field::displacement_z,
            0.0});
    }
    fuelsim::SteadyProblem problem(definition, mesh);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    const auto& first = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 0);
    const auto& second = fuelsim::cartesian::ProblemAccess::region_mesh(problem, 1);
    bool passed = check(result.completed && result.solve.converged, "two-region three-dimensional solve converges")
                  && check(first.nodes().size() == 8 && second.nodes().size() == 8
                               && fuelsim::cartesian::ProblemAccess::region_node_offset(problem, 1) == 8
                               && problem.dof_count() == 64,
                      "three-dimensional regions keep independent nodes and exact four-field offsets");
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    const double strain = 1.0e-5 * 25.0;
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const double origin = node < 8 ? 0.0 : 2.0;
        const auto& point = mesh.nodes()[node];
        passed =
            check(std::abs(result.solve.state[dofs.dof(fuelsim::Field::temperature, node)] - 325.0) < 2.0e-9
                      && std::abs(result.solve.state[dofs.dof(fuelsim::Field::displacement_x, node)]
                                  - strain * (point.x - origin))
                             < 2.0e-11
                      && std::abs(result.solve.state[dofs.dof(fuelsim::Field::displacement_y, node)] - strain * point.y)
                             < 2.0e-11
                      && std::abs(result.solve.state[dofs.dof(fuelsim::Field::displacement_z, node)] - strain * point.z)
                             < 2.0e-11,
                "each independent three-dimensional region matches free thermal expansion")
            && passed;
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
        check(result.completed && result.solve.converged, "shared-node three-dimensional solve converges")
        && check(problem.dof_count() == 48 && first.nodes().size() == 8 && second.nodes().size() == 8,
            "shared source nodes reduce the global four-field system to twelve unique nodes")
        && check(first_shared != first.source_node_ids().end() && second_shared != second.source_node_ids().end(),
            "shared source node is present in both material regions");
    if (passed) {
        const std::size_t first_local = static_cast<std::size_t>(first_shared - first.source_node_ids().begin());
        const std::size_t second_local = static_cast<std::size_t>(second_shared - second.source_node_ids().begin());
        passed = check(fuelsim::cartesian::ProblemAccess::dof_map(problem).global_node(0, first_local)
                           == fuelsim::cartesian::ProblemAccess::dof_map(problem).global_node(1, second_local),
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
    definition.regions = {{"primary", "primary", material(), 0.0, 300.0},
        {"secondary", "secondary", material(), 0.0, 400.0}};
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
    const std::size_t global_node =
        fuelsim::cartesian::ProblemAccess::region_node_offset(problem, 1)
        + static_cast<std::size_t>(local_position - secondary_mesh.source_node_ids().begin());
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
    bool passed = check(before.projected && at.projected && after.projected && before.primary_face != after.primary_face
                            && (at.primary_face == before.primary_face || at.primary_face == after.primary_face),
                      "three-dimensional internal primary edge has one owner and transfers ownership once")
                  && check(std::abs(before.contact_force - after.contact_force) < 1.0e-8 * force_scale,
                      "three-dimensional contact force is continuous across primary-face ownership transfer");
    std::vector<double> lost = problem.initial_state();
    for (std::size_t source_node : secondary_sources) {
        const auto local =
            std::find(secondary_mesh.source_node_ids().begin(), secondary_mesh.source_node_ids().end(), source_node);
        const std::size_t node = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, 1)
                                 + static_cast<std::size_t>(local - secondary_mesh.source_node_ids().begin());
        lost[dofs.dof(fuelsim::Field::displacement_y, node)] = 2.0;
    }
    bool rejected = false;
    try {
        problem.validate_state(lost);
    } catch (const std::domain_error&) {
        rejected = true;
    }
    return check(rejected,
               "three-dimensional contact rejects a thermal point or mechanical node outside the complete "
               "primary surface")
           && passed;
}

bool test_surface_contact_finite_sliding() {
    const fuelsim::UnstructuredHex8Mesh mesh = contact_projection_mesh();
    bool passed = true;
    for (const fuelsim::StrainFormulation strain :
        {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        fuelsim::SpatialDefinition definition;
        definition.regions = {{"primary", "primary", material(), 0.0, 300.0, -1, "", strain},
            {"secondary", "secondary", material(), 0.0, 300.0, -1, "", strain}};
        fuelsim::ContactDefinition contact;
        contact.name = "finite_sliding_interface";
        contact.primary = "primary_right";
        contact.secondary = "secondary_left";
        contact.mechanical = true;
        contact.penalty = 1.0e8;
        contact.friction_coefficient = 0.2;
        contact.friction_slip_tolerance = 1.0e-5;
        contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
        contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
        definition.contacts.push_back(contact);
        fuelsim::SteadyProblem problem(definition, mesh);
        const auto& view = fuelsim::cartesian::ProblemAccess::view(problem);
        const auto distinctive_dof = [&](double y) {
            const auto& region = view.region_mesh(0);
            const auto found = std::find_if(region.nodes().begin(), region.nodes().end(), [&](const auto& point) {
                return point.x == 1.0 && point.y == y && point.z == 0.0;
            });
            if (found == region.nodes().end())
                throw std::logic_error("HEX8 finite-sliding primary node was not found");
            return view.dof(fuelsim::Field::displacement_x,
                view.global_node(0, static_cast<std::size_t>(found - region.nodes().begin())));
        };
        const std::size_t lower_dof = distinctive_dof(-1.0), upper_dof = distinctive_dof(1.0);
        const auto owner_counts = [&](const std::vector<double>& state) {
            problem.validate_state(state);
            std::array<std::size_t, 3> result{};
            for (std::size_t contribution = view.volume_contribution_count(); contribution < view.contribution_count();
                ++contribution) {
                if (view.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
                    continue;
                std::vector<std::size_t> dofs;
                problem.contribution_dofs(contribution, dofs);
                const bool lower = std::find(dofs.begin(), dofs.end(), lower_dof) != dofs.end();
                const bool upper = std::find(dofs.begin(), dofs.end(), upper_dof) != dofs.end();
                result[0] += lower ? 1U : 0U;
                result[1] += upper ? 1U : 0U;
                result[2] += lower && upper ? 1U : 0U;
            }
            return result;
        };
        const fuelsim::ProblemStateSnapshot initial_snapshot = problem.capture_internal_state();
        std::vector<double> lower_state = problem.initial_state(), upper_state = problem.initial_state();
        for (std::size_t local = 0; local < view.region_mesh(1).nodes().size(); ++local) {
            const std::size_t global = view.global_node(1, local);
            lower_state[view.dof(fuelsim::Field::displacement_x, global)] = -1.0e-4;
            lower_state[view.dof(fuelsim::Field::displacement_y, global)] = -0.5;
            upper_state[view.dof(fuelsim::Field::displacement_x, global)] = -1.0e-4;
            upper_state[view.dof(fuelsim::Field::displacement_y, global)] = 0.5;
            upper_state[view.dof(fuelsim::Field::displacement_z, global)] = 2.0e-4;
        }
        const std::array<std::size_t, 3> lower_owners = owner_counts(lower_state),
                                         upper_owners = owner_counts(upper_state);
        passed = check(lower_owners == std::array<std::size_t, 3>{4, 0, 0}
                           && upper_owners == std::array<std::size_t, 3>{0, 4, 0},
                     "HEX8 finite sliding uniquely transfers all four node-centered points across a primary-face edge")
                 && passed;
        std::vector<double> edge_state = lower_state;
        for (std::size_t local = 0; local < view.region_mesh(1).nodes().size(); ++local) {
            const std::size_t global = view.global_node(1, local);
            edge_state[view.dof(fuelsim::Field::displacement_y, global)] = 0.125;
        }
        const std::array<std::size_t, 3> edge_owners = owner_counts(edge_state);
        passed = check(edge_owners[0] + edge_owners[1] == 4 && edge_owners[2] == 0,
                     "HEX8 finite-sliding points exactly on the internal primary edge have one owner")
                 && passed;
        const fuelsim::InterfaceSummary upper_summary =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, upper_state);
        const auto upper_nodes = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, upper_state);
        passed = check(upper_summary.projected_contact_nodes == 4 && upper_summary.unprojected_contact_nodes == 0
                           && upper_summary.total_contact_force > 0.0
                           && std::all_of(upper_nodes.begin(),
                               upper_nodes.end(),
                               [](const auto& node) { return node.projected && node.primary_face == 1; }),
                     "HEX8 finite-sliding output recovery follows the current upper primary face")
                 && passed;
        const double jacobian_error = mechanical_contact_directional_error(problem, upper_state, 1.0e-8);
        passed = check(jacobian_error < 2.0e-5,
                     "HEX8 finite-sliding contact Jacobian matches a centered directional difference")
                 && passed;
        problem.commit_internal_state(upper_state);
        const auto& upper_histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        passed = check(upper_histories.size() == 4
                           && std::all_of(upper_histories.begin(),
                               upper_histories.end(),
                               [](const auto& history) {
                                   return history.sliding && history.cartesian_tangent_basis_initialized
                                          && std::abs(history.cartesian_elastic_tangential_slip[1]) > 0.0
                                          && std::abs(history.cartesian_elastic_tangential_slip[2]) > 0.0;
                               }),
                     "HEX8 finite sliding commits four biaxial friction histories on the new primary face")
                 && passed;
        const std::array<std::size_t, 3> reverse_owners = owner_counts(lower_state);
        const double transported_jacobian_error = mechanical_contact_directional_error(problem, lower_state, 1.0e-8);
        passed = check(reverse_owners == std::array<std::size_t, 3>{4, 0, 0},
                     "HEX8 finite sliding uniquely transfers ownership back across the primary-face edge")
                 && check(transported_jacobian_error < 2.0e-5,
                     "HEX8 finite-sliding Jacobian includes the transported biaxial friction history")
                 && passed;
        std::vector<double> outside_state = lower_state;
        for (std::size_t local = 0; local < view.region_mesh(1).nodes().size(); ++local) {
            const std::size_t global = view.global_node(1, local);
            outside_state[view.dof(fuelsim::Field::displacement_y, global)] = 2.0;
        }
        bool outside_rejected = false;
        try {
            problem.validate_state(outside_state);
        } catch (const std::domain_error&) {
            outside_rejected = true;
        }
        passed = check(outside_rejected,
                     "HEX8 finite sliding rejects a state after its points leave the complete primary surface")
                 && passed;
        problem.commit_internal_state(lower_state);
        problem.restore_internal_state(initial_snapshot, problem.initial_state());
        const auto& restored = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        passed = check(std::all_of(restored.begin(),
                           restored.end(),
                           [](const auto& history) {
                               return !history.sliding && !history.cartesian_tangent_basis_initialized
                                      && history.cartesian_elastic_tangential_slip == std::array<double, 3>{};
                           }),
                     "HEX8 finite-sliding rollback restores friction and tangent-basis histories")
                 && passed;
    }
    return passed;
}

bool test_finite_sliding_search_tree() {
    constexpr std::size_t primary_face_count = 65;
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_faces;
    for (std::size_t face = 0; face < primary_face_count; ++face) {
        const double lower = static_cast<double>(face);
        elements.push_back(append_cuboid(nodes, primary_nodes, 0.0, 1.0, lower, lower + 1.0, 0.0, 1.0));
        blocks.push_back(1);
        primary_faces.push_back({face, 1});
    }
    const std::size_t secondary_element = elements.size();
    elements.push_back(append_cuboid(nodes, secondary_nodes, 1.0, 2.0, 0.1, 0.9, 0.0, 1.0));
    blocks.push_back(2);
    fuelsim::UnstructuredHex8Mesh mesh(std::move(nodes),
        std::move(elements),
        std::move(blocks),
        {{1, "primary"}, {2, "secondary"}},
        {},
        {{10, "primary_right", std::move(primary_faces)}, {20, "secondary_left", {{secondary_element, 3}}}});
    fuelsim::SpatialDefinition spatial;
    spatial.regions = {{"primary", "primary", material(), 0.0, 300.0},
        {"secondary", "secondary", material(), 0.0, 300.0}};
    fuelsim::ContactDefinition contact;
    contact.name = "finite_sliding_search_tree_interface";
    contact.primary = "primary_right";
    contact.secondary = "secondary_left";
    contact.mechanical = true;
    contact.penalty = 1.0e8;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    spatial.contacts.push_back(contact);
    fuelsim::SteadyProblem problem(spatial, mesh);
    const auto& view = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> state = problem.initial_state();
    for (std::size_t local = 0; local < view.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = view.global_node(1, local);
        state[view.dof(fuelsim::Field::displacement_x, global)] = -1.0e-4;
        state[view.dof(fuelsim::Field::displacement_y, global)] = 64.0;
    }
    problem.validate_state(state);
    const std::vector<fuelsim::CartesianContactNodeSummary> summaries =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, state);
    std::size_t active_contributions = 0;
    for (std::size_t contribution = view.volume_contribution_count(); contribution < view.contribution_count();
        ++contribution) {
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        if (!dofs.empty())
            ++active_contributions;
    }
    return check(summaries.size() == 4 && interface.total_contact_force > 0.0 && active_contributions == 4
                     && std::all_of(summaries.begin(),
                         summaries.end(),
                         [](const auto& summary) { return summary.projected && summary.primary_face == 64; }),
        "HEX8 finite sliding uses the search tree to find the last of 65 primary faces");
}

bool test_finite_sliding_end_to_end() {
    const fuelsim::UnstructuredHex8Mesh mesh = contact_projection_mesh();
    bool passed = true;
    for (const fuelsim::StrainFormulation strain :
        {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        fuelsim::SpatialDefinition spatial;
        spatial.regions = {{"primary", "primary", material(), 0.0, 300.0, -1, "", strain},
            {"secondary", "secondary", material(), 0.0, 300.0, -1, "", strain}};
        fuelsim::ContactDefinition contact;
        contact.name = "finite_sliding_end_to_end_interface";
        contact.primary = "primary_right";
        contact.secondary = "secondary_left";
        contact.mechanical = true;
        contact.penalty = 1.0e8;
        contact.friction_coefficient = 0.2;
        contact.friction_slip_tolerance = 1.0e-5;
        contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
        contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
        spatial.contacts.push_back(contact);
        for (const fuelsim::Field field :
            {fuelsim::Field::displacement_x, fuelsim::Field::displacement_y, fuelsim::Field::displacement_z})
            spatial.boundary_conditions.push_back(
                {"fix_primary", fuelsim::BoundaryConditionType::dirichlet, "primary_all", field, 0.0});
        spatial.boundary_conditions.push_back({"move_secondary_x",
            fuelsim::BoundaryConditionType::dirichlet,
            "secondary_all",
            fuelsim::Field::displacement_x,
            -1.0e-4,
            true});
        spatial.boundary_conditions.push_back({"move_secondary_y",
            fuelsim::BoundaryConditionType::dirichlet,
            "secondary_all",
            fuelsim::Field::displacement_y,
            0.5,
            true});
        spatial.boundary_conditions.push_back({"move_secondary_z",
            fuelsim::BoundaryConditionType::dirichlet,
            "secondary_all",
            fuelsim::Field::displacement_z,
            2.0e-4,
            true});
        spatial.boundary_conditions.push_back({"primary_temperature",
            fuelsim::BoundaryConditionType::dirichlet,
            "primary_all",
            fuelsim::Field::temperature,
            300.0});
        spatial.boundary_conditions.push_back({"secondary_temperature",
            fuelsim::BoundaryConditionType::dirichlet,
            "secondary_all",
            fuelsim::Field::temperature,
            300.0});
        fuelsim::SteadyProblem problem(spatial, mesh);
        fuelsim::SolverOptions options = solver_options();
        options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
        const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {4, 0.5, 4, 1.0e-6}, options);
        const std::vector<fuelsim::CartesianContactNodeSummary> summaries =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, result.solve.state);
        const fuelsim::InterfaceSummary interface =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, result.solve.state);
        const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        passed = check(result.completed && result.solve.converged && result.rejected_steps.empty(),
                     "HEX8 finite-sliding end-to-end load path converges without a rejected load step")
                 && check(summaries.size() == 4 && interface.total_contact_force > 0.0
                              && std::all_of(summaries.begin(),
                                  summaries.end(),
                                  [](const auto& summary) { return summary.projected && summary.primary_face == 1; }),
                     "HEX8 finite-sliding end-to-end solve transfers the complete secondary face")
                 && check(histories.size() == 4
                              && std::all_of(histories.begin(),
                                  histories.end(),
                                  [](const auto& history) {
                                      return history.sliding && history.cartesian_tangent_basis_initialized
                                             && std::abs(history.cartesian_elastic_tangential_slip[1]) > 0.0
                                             && std::abs(history.cartesian_elastic_tangential_slip[2]) > 0.0;
                                  }),
                     "HEX8 finite-sliding end-to-end solve commits two-component friction histories")
                 && passed;
    }
    return passed;
}

bool test_boundary_configuration_selection(const fuelsim::UnstructuredHex8Mesh& mesh) {
    const auto warning_count = [&](fuelsim::BoundaryConditionType type,
                                   fuelsim::Field field,
                                   bool finite_strain,
                                   bool current_configuration,
                                   fuelsim::Hex8ElementFormulation element_formulation,
                                   bool configuration_explicit = true) {
        fuelsim::SpatialDefinition definition;
        definition.regions.push_back({"solid", "solid", material(), 0.0, 300.0});
        definition.regions.front().strain_formulation =
            finite_strain ? fuelsim::StrainFormulation::finite : fuelsim::StrainFormulation::small;
        definition.regions.front().hex8_element_formulation = element_formulation;
        fuelsim::BoundaryConditionDefinition boundary{"boundary", type, "x2", field, 1.0e6};
        if (type == fuelsim::BoundaryConditionType::convection) {
            boundary.heat_transfer_coefficient = 1000.0;
            boundary.ambient_temperature = 300.0;
        }
        boundary.use_displaced_geometry = current_configuration;
        boundary.configuration_explicit = configuration_explicit;
        definition.boundary_conditions.push_back(boundary);
        fuelsim::SteadyProblem problem(std::move(definition), mesh);
        return fuelsim::cartesian::ProblemAccess::dof_map(problem).configuration_warnings().size();
    };
    const auto check_type = [&](fuelsim::BoundaryConditionType type,
                                fuelsim::Field field,
                                fuelsim::Hex8ElementFormulation element_formulation,
                                const char* element_name,
                                const char* boundary_name) {
        const std::size_t reference_small = warning_count(type, field, false, false, element_formulation);
        const std::size_t current_small = warning_count(type, field, false, true, element_formulation);
        const std::size_t reference_finite = warning_count(type, field, true, false, element_formulation);
        const std::size_t current_finite = warning_count(type, field, true, true, element_formulation);
        bool result = check(reference_small == 0 && current_small == 1 && reference_finite == 1 && current_finite == 0,
            std::string(element_name) + " " + boundary_name
                + " accepts both configurations and warns for non-recommended choices");
        const std::size_t default_small = warning_count(type, field, false, false, element_formulation, false);
        const std::size_t default_finite = warning_count(type, field, true, false, element_formulation, false);
        result = check(default_small == 0 && default_finite == 0,
                     std::string(element_name) + " " + boundary_name
                         + " omits configuration without warning for the strain-dependent recommendation")
                 && result;
        return result;
    };
    const auto convection_has_geometry_columns = [&](bool finite_strain,
                                                     bool current_configuration,
                                                     fuelsim::Hex8ElementFormulation element_formulation,
                                                     bool configuration_explicit) {
        fuelsim::SpatialDefinition definition;
        definition.regions.push_back({"solid", "solid", material(), 0.0, 300.0});
        definition.regions.front().strain_formulation =
            finite_strain ? fuelsim::StrainFormulation::finite : fuelsim::StrainFormulation::small;
        definition.regions.front().hex8_element_formulation = element_formulation;
        fuelsim::BoundaryConditionDefinition boundary{"convection",
            fuelsim::BoundaryConditionType::convection,
            "x2",
            fuelsim::Field::temperature,
            0.0};
        boundary.heat_transfer_coefficient = 1000.0;
        boundary.ambient_temperature = 300.0;
        boundary.use_displaced_geometry = current_configuration;
        boundary.configuration_explicit = configuration_explicit;
        definition.boundary_conditions.push_back(boundary);
        fuelsim::SteadyProblem problem(std::move(definition), mesh);
        const fuelsim::cartesian::SpatialAssembly& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
            if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::convection)
                continue;
            std::vector<unsigned char> pattern;
            problem.contribution_jacobian_pattern(contribution, pattern);
            for (std::size_t row = 0; row < 4; ++row)
                for (std::size_t column = 4; column < 16; ++column)
                    if (pattern[row * 16 + column] != 0U)
                        return true;
            return false;
        }
        throw std::logic_error("convection contribution is missing from the configuration test");
    };
    bool passed = true;
    for (const auto element : {fuelsim::Hex8ElementFormulation::c3d8t, fuelsim::Hex8ElementFormulation::c3d8rt}) {
        const char* element_name = element == fuelsim::Hex8ElementFormulation::c3d8t ? "C3D8T" : "C3D8RT";
        passed = check_type(fuelsim::BoundaryConditionType::pressure,
                     fuelsim::Field::displacement_x,
                     element,
                     element_name,
                     "pressure")
                 && check_type(fuelsim::BoundaryConditionType::traction,
                     fuelsim::Field::displacement_x,
                     element,
                     element_name,
                     "traction")
                 && check_type(fuelsim::BoundaryConditionType::heat_flux,
                     fuelsim::Field::temperature,
                     element,
                     element_name,
                     "surface heat flux")
                 && check_type(fuelsim::BoundaryConditionType::convection,
                     fuelsim::Field::temperature,
                     element,
                     element_name,
                     "convection")
                 && passed;
        const bool explicit_reference = convection_has_geometry_columns(true, false, element, true);
        const bool explicit_current = convection_has_geometry_columns(false, true, element, true);
        const bool default_small = convection_has_geometry_columns(false, false, element, false);
        const bool default_finite = convection_has_geometry_columns(true, false, element, false);
        passed =
            check(!explicit_reference && explicit_current && !default_small && default_finite,
                std::string(element_name) + " convection honors explicit configuration and strain-dependent defaults")
            && passed;
    }
    return passed;
}

fuelsim::SpatialDefinition inelastic_definition(bool creep,
    bool plasticity,
    fuelsim::Hex8ElementFormulation element_formulation,
    fuelsim::StrainFormulation strain_formulation = fuelsim::StrainFormulation::small) {
    fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 10.0, 2.0e11, 0.3, 0.0, 600.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    if (creep)
        properties = fuelsim::test::with_norton(std::move(properties), 1.0e-4, 1.0e8, 3.0);
    if (plasticity)
        properties = fuelsim::test::with_plasticity(std::move(properties), 2.0e8, 2.0e9);
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"solid", "solid", std::move(properties), 0.0, 600.0});
    definition.regions.back().hex8_element_formulation = element_formulation;
    definition.regions.back().strain_formulation = strain_formulation;
    definition.boundary_conditions = {
        {"temperature", fuelsim::BoundaryConditionType::dirichlet, "all", fuelsim::Field::temperature, 600.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
        {"traction", fuelsim::BoundaryConditionType::traction, "x2", fuelsim::Field::displacement_x, 2.01e8, true},
    };
    return definition;
}

bool test_inelastic_branches(const fuelsim::PetscSession& session,
    const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::string& checkpoint_path) {
    struct RepresentativeCase final {
        fuelsim::Hex8ElementFormulation formulation;
        fuelsim::StrainFormulation strain;
        std::array<bool, 2> branch;
    };

    const std::array<RepresentativeCase, 4> cases = {
        RepresentativeCase{fuelsim::Hex8ElementFormulation::c3d8t, fuelsim::StrainFormulation::small, {false, true}},
        RepresentativeCase{fuelsim::Hex8ElementFormulation::c3d8t, fuelsim::StrainFormulation::small, {true, false}},
        RepresentativeCase{fuelsim::Hex8ElementFormulation::c3d8rt, fuelsim::StrainFormulation::small, {true, true}},
        RepresentativeCase{fuelsim::Hex8ElementFormulation::c3d8rt, fuelsim::StrainFormulation::finite, {true, true}},
    };
    bool passed = true;
    for (const RepresentativeCase& current : cases) {
        const fuelsim::SpatialDefinition definition =
            inelastic_definition(current.branch[0], current.branch[1], current.formulation, current.strain);
        fuelsim::TransientProblem problem(definition, mesh);
        fuelsim::SolverOptions options = solver_options();
        options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
        options.maximum_iterations = 30;
        const fuelsim::TransientResult result =
            fuelsim::solve_transient(problem, {1.0, 0.1, 0.1, 0.1, 1.0, 0.5, 0, 1.0}, options);
        double maximum_plastic = 0.0, maximum_creep = 0.0;
        bool finite = true, material_point_count_matches = true;
        const std::size_t expected_points = current.formulation == fuelsim::Hex8ElementFormulation::c3d8rt ? 1 : 8;
        for (std::size_t element = 0; element < 2; ++element) {
            const fuelsim::CartesianMaterialHistory& history =
                fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element);
            material_point_count_matches = material_point_count_matches && history.size() == expected_points;
            for (const fuelsim::CartesianMaterialPointState& point : history) {
                maximum_plastic = std::max(maximum_plastic, point.equivalent_plastic_strain);
                maximum_creep = std::max(maximum_creep, point.equivalent_creep_strain);
                finite = finite && std::isfinite(point.stress.xx) && std::isfinite(point.stress.yy)
                         && std::isfinite(point.stress.zz) && std::isfinite(point.stress.xy)
                         && std::isfinite(point.stress.yz) && std::isfinite(point.stress.xz);
            }
        }
        passed = check(result.completed && result.accepted_steps.size() == 10 && finite && material_point_count_matches
                           && (current.branch[1] ? maximum_plastic > 0.0 : maximum_plastic == 0.0)
                           && (current.branch[0] ? maximum_creep > 0.0 : maximum_creep == 0.0),
                     "representative full- and reduced-integration three-dimensional inelastic paths solve "
                     "and commit the formulation-specific material-point count")
                 && passed;
        if (current.formulation == fuelsim::Hex8ElementFormulation::c3d8rt
            && current.strain == fuelsim::StrainFormulation::finite && current.branch[0] && current.branch[1]) {
            session.collective_root_action(
                [&]() { fuelsim::write_transient_checkpoint(checkpoint_path, problem, 0.025); });
            fuelsim::TransientProblem restored(definition, mesh);
            const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
            bool history_matches = restored.committed_solution() == problem.committed_solution();
            for (std::size_t element = 0; element < 2; ++element) {
                const fuelsim::CartesianMaterialHistory& before =
                    fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element);
                const fuelsim::CartesianMaterialHistory& after =
                    fuelsim::cartesian::ProblemAccess::material_history(restored, 0, element);
                history_matches = history_matches && before.size() == 1 && after.size() == 1
                                  && same_material_point(before.front(), after.front());
            }
            passed = check(restored_step == 0.025 && restored.committed_time() == problem.committed_time()
                               && history_matches,
                         "C3D8RT checkpoint restores the exact nodal state and single coupled material point")
                     && passed;
        }
    }
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_hex8_solver_tests <steady.e> <transient.e> <checkpoint.bin>\n";
        return 2;
    }
    const bool mpi_only = argc > 4 && std::string(argv[4]) == "--mpi-only";
    fuelsim::PetscSession session(argc, argv, "fuelsim HEX8 solver tests\n");
    const fuelsim::UnstructuredHex8Mesh mesh = two_element_mesh();
    bool passed = test_steady(session, mesh, argv[1]);
    passed = test_transient(session, mesh, argv[3], argv[2]) && passed;
    passed = test_shared_nodes(session) && passed;
    passed = test_finite_sliding_end_to_end() && passed;
    if (!mpi_only) {
        passed = test_small_strain_steady_predictor(mesh) && passed;
        passed = test_convection_boundary(mesh) && passed;
        passed = test_multiple_regions() && passed;
        passed = test_contact_projection_transfer() && passed;
        passed = test_surface_contact_finite_sliding() && passed;
        passed = test_finite_sliding_search_tree() && passed;
        passed = test_boundary_configuration_selection(mesh) && passed;
        passed = test_inelastic_branches(session, mesh, argv[3]) && passed;
    }
    session.collective_root_action([&]() {
        (void)std::remove(argv[1]);
        (void)std::remove(argv[2]);
        (void)std::remove(argv[3]);
    });
    if (passed && session.rank() == 0)
        std::cout << "HEX8 steady, transient, message-passing, results, and checkpoint tests passed\n";
    return passed ? 0 : 1;
}
