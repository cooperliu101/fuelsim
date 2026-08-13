#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
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
        {{10, "all", all_faces}, {11, "x0", {{0, 3}}}, {12, "y0", {{0, 0}, {1, 0}}}, {13, "z0", {{0, 4}, {1, 4}}}});
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
    const std::size_t expected_begin = problem.contribution_count() * static_cast<std::size_t>(session.rank()) /
                                       static_cast<std::size_t>(session.size());
    const std::size_t expected_end = problem.contribution_count() * static_cast<std::size_t>(session.rank() + 1) /
                                     static_cast<std::size_t>(session.size());
    passed = check(result.solve.local_contribution_begin == expected_begin &&
                       result.solve.local_contribution_end == expected_end,
                 "each message-passing rank owns only its exact HEX8 contribution interval") &&
             passed;
    session.collective_root_action(
        [&]() { fuelsim::write_steady_results(results_path, mesh, problem, result.solve.state); });
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
} // namespace
int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_hex8_solver_tests <steady.e> <transient.e> <checkpoint.bin>\n";
        return 2;
    }
    fuelsim::PetscSession session(argc, argv, "fuelsim HEX8 solver tests\n");
    const fuelsim::UnstructuredHex8Mesh mesh = two_element_mesh();
    bool passed = test_steady(session, mesh, argv[1]);
    passed = test_transient(session, mesh, argv[3], argv[2]) && passed;
    passed = test_multiple_regions() && passed;
    session.collective_root_action([&]() {
        (void)std::remove(argv[1]);
        (void)std::remove(argv[2]);
        (void)std::remove(argv[3]);
    });
    if (passed && session.rank() == 0)
        std::cout << "HEX8 steady, transient, message-passing, results, and checkpoint tests passed\n";
    return passed ? 0 : 1;
}
