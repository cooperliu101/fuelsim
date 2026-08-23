#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
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

fuelsim::Hex20Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map, double x0, double x1, double y0, double y1, double z0,
    double z1) {
    const std::array<fuelsim::CartesianPoint3, 8> corners = {{{x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
        {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}}};
    const std::array<std::pair<std::size_t, std::size_t>, 12> edges = {
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}}};
    std::array<fuelsim::CartesianPoint3, 20> points{};
    std::copy(corners.begin(), corners.end(), points.begin());
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const auto& first = corners[edges[edge].first];
        const auto& second = corners[edges[edge].second];
        points[8 + edge] = {0.5 * (first.x + second.x), 0.5 * (first.y + second.y), 0.5 * (first.z + second.z)};
    }
    fuelsim::Hex20Element element{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const std::array<double, 3> key = {points[local].x, points[local].y, points[local].z};
        const auto inserted = node_map.emplace(key, nodes.size());
        if (inserted.second) nodes.push_back(points[local]);
        element.nodes[local] = inserted.first->second;
    }
    return element;
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

double mechanical_contact_directional_error(
    fuelsim::SteadyProblem& problem, const std::vector<double>& global_state, double perturbation) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum_error = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local) local_state[local] = global_state[dofs[local]];
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
    contact.friction_elastic_slip = 1.0e-5;
    spatial.contacts.push_back(contact);
    fuelsim::SteadyProblem problem(spatial, contact_mesh);
    const auto& view = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto summary =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, problem.initial_state());
    bool passed = check(summary.size() == 8, "HEX20 contact exposes all eight quadratic secondary face nodes");
    double area_sum = 0.0;
    std::size_t positive_areas = 0, negative_areas = 0;
    for (const auto& node : summary) {
        if (node.tributary_area > 0.0)
            ++positive_areas;
        else if (node.tributary_area < 0.0)
            ++negative_areas;
        passed = check(node.projected && node.tributary_area != 0.0 && std::abs(node.gap) < 1.0e-10,
                     "HEX20 surface-to-surface contact projects every secondary integration point to the face") &&
                 passed;
        area_sum += node.tributary_area;
    }
    passed = check(positive_areas == 8 && negative_areas == 0 && std::abs(area_sum - 1.0) < 1.0e-12,
                 "HEX20 Abaqus-style averaged contact stores eight positive constraint areas") &&
             passed;
    const auto interface = fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, problem.initial_state());
    passed = check(interface.projected_contact_nodes == 8 && interface.unprojected_contact_nodes == 0,
                 "HEX20 contact validation preserves unique primary projection") &&
             passed;
    passed = check(view.contribution_count() > view.volume_contribution_count(),
                 "HEX20 contact contributes thermal and mechanical surface kernels") &&
             passed;
    const fuelsim::ProblemStateSnapshot snapshot = problem.capture_internal_state();
    std::vector<double> sticking_state = problem.initial_state();
    for (std::size_t local_node = 0; local_node < view.hex20_region_mesh(1).nodes().size(); ++local_node) {
        const std::size_t global = view.global_node(1, local_node);
        sticking_state[view.dof(fuelsim::Field::displacement_x, global)] = -1.0e-4;
        sticking_state[view.dof(fuelsim::Field::displacement_y, global)] = 1.0e-6;
    }
    problem.validate_state(sticking_state);
    const double sticking_jacobian_error = mechanical_contact_directional_error(problem, sticking_state, 1.0e-8);
    problem.commit_internal_state(sticking_state);
    const auto& sticking_history = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    passed = check(sticking_history.size() == 8 &&
                       std::all_of(sticking_history.begin(), sticking_history.end(),
                           [](const auto& history) {
                               return !history.sliding &&
                                      std::abs(history.cartesian_elastic_tangential_slip[1] - 1.0e-6) < 1.0e-14;
                           }),
                 "HEX20 averaged friction commits all eight elastic sticking histories") &&
             passed;
    passed = check(sticking_jacobian_error < 1.0e-7,
                 "HEX20 averaged sticking-friction Jacobian matches a centered directional difference") &&
             passed;
    problem.restore_internal_state(snapshot, problem.initial_state());
    std::vector<double> sliding_state = problem.initial_state();
    for (std::size_t local_node = 0; local_node < view.hex20_region_mesh(1).nodes().size(); ++local_node) {
        const std::size_t global = view.global_node(1, local_node);
        sliding_state[view.dof(fuelsim::Field::displacement_x, global)] = -1.0e-4;
        sliding_state[view.dof(fuelsim::Field::displacement_y, global)] = 1.0e-4;
    }
    problem.validate_state(sliding_state);
    problem.commit_internal_state(sliding_state);
    const auto& sliding_history = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    passed = check(sliding_history.size() == 8 &&
                       std::all_of(sliding_history.begin(), sliding_history.end(),
                           [](const auto& history) {
                               return history.sliding && std::abs(history.cartesian_elastic_tangential_slip[1]) > 0.0;
                           }),
                 "HEX20 surface-to-surface friction commits one sliding history at each averaged constraint") &&
             passed;
    problem.restore_internal_state(snapshot, problem.initial_state());
    const auto& restored_history = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    passed =
        check(std::all_of(restored_history.begin(), restored_history.end(),
                  [](const auto& history) {
                      return !history.sliding && history.cartesian_elastic_tangential_slip == std::array<double, 3>{};
                  }),
            "HEX20 surface-to-surface friction rollback restores every averaged-constraint history") &&
        passed;
    return passed;
}

bool test_surface_contact_fixed_reference_graph() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    const fuelsim::Hex20Element primary_lower = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
                                primary_upper = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 1.0, 2.0, 0.0, 1.0),
                                secondary = append_cuboid(nodes, secondary_nodes, 1.0, 2.0, 0.0, 2.0, 0.0, 1.0);
    fuelsim::UnstructuredHex20Mesh mesh(std::move(nodes), {primary_lower, primary_upper, secondary}, {1, 1, 2},
        {{1, "primary"}, {2, "secondary"}}, {},
        {{10, "primary_right", {{0, 1}, {1, 1}}}, {20, "secondary_left", {{2, 3}}}});
    fuelsim::SpatialDefinition spatial;
    spatial.regions = {
        {"primary", "primary", material(), 0.0, 300.0}, {"secondary", "secondary", material(), 0.0, 300.0}};
    fuelsim::ContactDefinition contact;
    contact.name = "split_interface";
    contact.primary = "primary_right";
    contact.secondary = "secondary_left";
    contact.mechanical = true;
    contact.penalty = 1.0e8;
    spatial.contacts.push_back(contact);
    fuelsim::SteadyProblem problem(spatial, mesh);
    const auto& view = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto distinctive_dof = [&](double y) {
        const auto& region = view.hex20_region_mesh(0);
        const auto found = std::find_if(region.nodes().begin(), region.nodes().end(),
            [&](const auto& point) { return point.x == 1.0 && point.y == y && point.z == 0.0; });
        if (found == region.nodes().end()) throw std::logic_error("HEX20 split-primary node was not found");
        const std::size_t local = static_cast<std::size_t>(found - region.nodes().begin());
        return view.dof(fuelsim::Field::displacement_x, view.global_node(0, local));
    };
    const std::size_t lower_dof = distinctive_dof(0.0), upper_dof = distinctive_dof(2.0);
    const auto owner_counts = [&](const std::vector<double>& state) {
        problem.validate_state(state);
        std::array<std::size_t, 3> result{};
        for (std::size_t contribution = view.volume_contribution_count(); contribution < view.contribution_count();
            ++contribution) {
            std::vector<std::size_t> dofs;
            problem.contribution_dofs(contribution, dofs);
            const bool lower = std::find(dofs.begin(), dofs.end(), lower_dof) != dofs.end();
            const bool upper = std::find(dofs.begin(), dofs.end(), upper_dof) != dofs.end();
            if (!lower && !upper)
                throw std::logic_error("HEX20 surface integration point has no selected primary face");
            result[0] += lower ? 1U : 0U;
            result[1] += upper ? 1U : 0U;
            result[2] += lower && upper ? 1U : 0U;
        }
        return result;
    };
    std::vector<double> state = problem.initial_state();
    const std::array<std::size_t, 3> reference_owners = owner_counts(state);
    bool passed = check(reference_owners[0] > 0 && reference_owners[1] > 0 && reference_owners[2] > 0,
        "HEX20 averaged small-sliding constraints include both primary faces and cross-face support");
    for (std::size_t local = 0; local < view.hex20_region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = view.global_node(1, local);
        state[view.dof(fuelsim::Field::displacement_y, global)] = 0.05;
    }
    passed = check(owner_counts(state) == reference_owners,
                 "HEX20 averaged small-sliding reference graph remains fixed after a tangential displacement") &&
             passed;
    std::vector<double> tilted = problem.initial_state();
    for (std::size_t region = 0; region < 2; ++region)
        for (std::size_t local = 0; local < view.hex20_region_mesh(region).nodes().size(); ++local) {
            const auto& point = view.hex20_region_mesh(region).nodes()[local];
            const std::size_t global = view.global_node(region, local);
            tilted[view.dof(fuelsim::Field::displacement_x, global)] = 0.05 * point.y - (region == 1 ? 1.0e-4 : 0.0);
        }
    problem.validate_state(tilted);
    const fuelsim::InterfaceSummary tilted_summary =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, tilted);
    passed = check(tilted_summary.total_contact_force > 0.0 && tilted_summary.unprojected_contact_nodes == 0,
                 "HEX20 averaged small-sliding contact remains active under a compatible interface tilt") &&
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
                        test_contact_projection(mesh) && test_surface_contact_fixed_reference_graph();
    session.collective_root_action([&]() {
        (void)std::remove(argv[1]);
        (void)std::remove(argv[2]);
        (void)std::remove(argv[3]);
        (void)std::remove(transient_results.c_str());
    });
    if (passed && session.rank() == 0) std::cout << "All HEX20-U2/T1 solver tests passed\n";
    return passed ? 0 : 1;
}
