#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::UnstructuredHex20Mesh unit_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {0.5, 0.0, 0.0}, {1.0, 0.5, 0.0},
        {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5},
        {0.5, 0.0, 1.0}, {1.0, 0.5, 1.0}, {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0}};
    fuelsim::Hex20Element element{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node) element.nodes[node] = node;
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes), {element}, {1}, {{1, "solid"}}, {},
        {{10, "all", {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}}, {11, "x0", {{0, 3}}}, {12, "y0", {{0, 0}}},
            {13, "z0", {{0, 4}}}, {14, "x1", {{0, 1}}}});
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.25, 1.0e-5, 300.0, 0.0, 0.0, 0.0, 6000.0, 1000.0);
}

fuelsim::SpatialDefinition definition(bool fixed_temperature) {
    fuelsim::SpatialDefinition value;
    value.regions.push_back({"solid", "solid", material(), fixed_temperature ? 0.0 : 6.0e6, 300.0});
    value.boundary_conditions = {
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0}};
    if (fixed_temperature)
        value.boundary_conditions.push_back(
            {"temperature", fuelsim::BoundaryConditionType::dirichlet, "all", fuelsim::Field::temperature, 400.0});
    return value;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions value;
    value.absolute_tolerance = 1.0e-8;
    value.relative_tolerance = 1.0e-11;
    value.maximum_iterations = 12;
    value.linear_solver = fuelsim::SolverOptions::LinearSolver::gmres;
    value.preconditioner = fuelsim::SolverOptions::Preconditioner::field_split;
    value.field_residual_scaling = true;
    value.linear_relative_tolerance = 1.0e-11;
    value.maximum_linear_iterations = 300;
    return value;
}

bool check_uniform(const fuelsim::cartesian::SpatialAssembly& spatial, const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::vector<double>& state, double temperature) {
    const auto& fields = spatial.field_layout();
    const double strain = 1.0e-5 * (temperature - 300.0);
    bool passed = true;
    for (std::size_t node = 0; node < 8; ++node)
        passed =
            check(std::abs(state[fields[0].begin + spatial.global_temperature_node(0, node)] - temperature) < 3.0e-8,
                "HEX20 corner temperature matches the uniform analytic solution") &&
            passed;
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const std::size_t global = spatial.global_node(0, node);
        passed = check(std::abs(state[fields[1].begin + global] - strain * mesh.nodes()[node].x) < 3.0e-10 &&
                           std::abs(state[fields[2].begin + global] - strain * mesh.nodes()[node].y) < 3.0e-10 &&
                           std::abs(state[fields[3].begin + global] - strain * mesh.nodes()[node].z) < 3.0e-10,
                     "HEX20 quadratic displacement field matches uniform free thermal expansion") &&
                 passed;
    }
    return passed;
}

bool test_steady_and_io(const fuelsim::PetscSession& session, const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::string& mesh_path, const std::string& results_path) {
    session.collective_root_action([&]() { fuelsim::write_exodus_hex20(mesh_path, mesh); });
    const fuelsim::UnstructuredHex20Mesh restored = fuelsim::read_exodus_hex20(mesh_path);
    bool passed =
        check(fuelsim::exodus_uses_hex20(mesh_path) && restored.elements()[0].nodes == mesh.elements()[0].nodes,
            "Exodus HEX20 topology detection and connectivity round trip are exact");
    fuelsim::SteadyProblem problem(definition(true), restored);
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto partition = problem.contribution_partition(
        static_cast<std::size_t>(session.rank()), static_cast<std::size_t>(session.size()));
    passed = check(result.completed && result.solve.converged && problem.dof_count() == 68,
                 "non-contact HEX20 steady problem converges with 8 temperature and 60 displacement DOFs") &&
             check(result.solve.local_contribution_begin == partition.first &&
                       result.solve.local_contribution_end == partition.second,
                 "each message-passing rank owns only its exact HEX20 contribution interval") &&
             check_uniform(spatial, restored, result.solve.state, 400.0) && passed;
    session.collective_root_action(
        [&]() { fuelsim::write_steady_results(results_path, restored, problem, result.solve.state); });
    return passed;
}

bool test_transient_restart(const fuelsim::PetscSession& session, const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::string& checkpoint_path, const std::string& results_path) {
    const fuelsim::SpatialDefinition spatial = definition(false);
    fuelsim::TransientProblem problem(spatial, mesh);
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, {1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 2, 0.0}, solver_options());
    bool passed =
        check(result.completed && result.accepted_steps.size() == 1,
            "HEX20 Backward Euler transient accepts one physical time step") &&
        check_uniform(fuelsim::cartesian::ProblemAccess::view(problem), mesh, problem.committed_solution(), 301.0);
    const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, 0, 0);
    passed = check(history.size() == 27, "HEX20 transient stores 27 committed material-point states") && passed;
    session.collective_root_action([&]() {
        fuelsim::write_transient_checkpoint(checkpoint_path, problem, 0.25);
        fuelsim::ExodusTransientResultsWriter writer(results_path, mesh, problem);
        writer.append(problem);
    });
    fuelsim::TransientProblem restored(spatial, mesh);
    const double next_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    passed = check(next_time_step == 0.25 && restored.committed_solution() == problem.committed_solution() &&
                       fuelsim::cartesian::ProblemAccess::material_history(restored, 0, 0).size() == 27,
                 "HEX20 checkpoint restores the exact nodal state and all 27 material points") &&
             passed;
    return passed;
}

bool test_contact_projection(const fuelsim::UnstructuredHex20Mesh& mesh) {
    (void)mesh;
    std::vector<fuelsim::CartesianPoint3> nodes;
    const auto append_cube = [&nodes](double origin) {
        const std::array<fuelsim::CartesianPoint3, 8> corners = {
            {{origin, 0.0, 0.0}, {origin + 1.0, 0.0, 0.0}, {origin + 1.0, 1.0, 0.0}, {origin, 1.0, 0.0},
                {origin, 0.0, 1.0}, {origin + 1.0, 0.0, 1.0}, {origin + 1.0, 1.0, 1.0}, {origin, 1.0, 1.0}}};
        for (const auto& point : corners) nodes.push_back(point);
        const std::array<std::pair<std::size_t, std::size_t>, 12> edges = {
            {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}}};
        for (const auto& edge : edges) {
            const auto& first = corners[edge.first];
            const auto& second = corners[edge.second];
            nodes.push_back({0.5 * (first.x + second.x), 0.5 * (first.y + second.y), 0.5 * (first.z + second.z)});
        }
    };
    append_cube(0.0);
    append_cube(1.0);
    fuelsim::Hex20Element first{}, second{};
    for (std::size_t node = 0; node < 20; ++node) {
        first.nodes[node] = node;
        second.nodes[node] = 20 + node;
    }
    fuelsim::UnstructuredHex20Mesh contact_mesh(std::move(nodes), {first, second}, {1, 2},
        {{1, "primary"}, {2, "secondary"}}, {}, {{10, "primary_right", {{0, 1}}}, {20, "secondary_left", {{1, 3}}}});
    fuelsim::SpatialDefinition spatial;
    spatial.regions = {
        {"primary", "primary", material(), 0.0, 300.0}, {"secondary", "secondary", material(), 0.0, 400.0}};
    fuelsim::ContactDefinition contact;
    contact.name = "interface";
    contact.primary = "primary_right";
    contact.secondary = "secondary_left";
    contact.thermal = true;
    contact.mechanical = true;
    contact.gap_conductivity = 1.0;
    contact.minimum_gap = 1.0e-6;
    contact.penalty = 1.0e8;
    contact.friction_coefficient = 0.1;
    spatial.contacts.push_back(contact);
    fuelsim::SteadyProblem problem(spatial, contact_mesh);
    const auto& view = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto summary =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, problem.initial_state());
    bool passed = check(summary.size() == 8, "HEX20 contact exposes all eight quadratic secondary face nodes");
    for (const auto& node : summary)
        passed = check(node.projected && std::abs(node.tributary_area) > 0.0 && std::abs(node.gap) < 1.0e-10,
                     "HEX20 Q8 mechanical contact projects every secondary node with nonzero consistent area") &&
                 passed;
    const auto interface = fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, problem.initial_state());
    passed = check(interface.projected_contact_nodes == 8 && interface.unprojected_contact_nodes == 0,
                 "HEX20 contact validation preserves unique primary projection") &&
             passed;
    passed = check(view.contribution_count() > view.volume_contribution_count(),
                 "HEX20 contact contributes thermal and mechanical surface kernels") &&
             passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_hex20_solver_tests <mesh.e> <steady.e> <checkpoint.bin>\n";
        return 2;
    }
    fuelsim::PetscSession session(argc, argv, "fuelsim HEX20-U2/T1 solver tests\n");
    const fuelsim::UnstructuredHex20Mesh mesh = unit_mesh();
    const std::string transient_results = std::string(argv[2]) + ".transient.e";
    const bool passed = test_steady_and_io(session, mesh, argv[1], argv[2]) &&
                        test_transient_restart(session, mesh, argv[3], transient_results) &&
                        test_contact_projection(mesh);
    session.collective_root_action([&]() {
        (void)std::remove(argv[1]);
        (void)std::remove(argv[2]);
        (void)std::remove(argv[3]);
        (void)std::remove(transient_results.c_str());
    });
    if (passed && session.rank() == 0) std::cout << "All HEX20-U2/T1 solver tests passed\n";
    return passed ? 0 : 1;
}
