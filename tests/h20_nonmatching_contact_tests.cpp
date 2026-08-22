#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
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

fuelsim::UnstructuredHex20Mesh nonmatching_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    std::vector<fuelsim::Hex20Element> elements;
    for (std::size_t segment = 0; segment < 4; ++segment) {
        const double lower = 0.5 * static_cast<double>(segment);
        elements.push_back(append_cuboid(nodes, primary_nodes, 0.0, 1.0, lower, lower + 0.5, 0.0, 1.0));
    }
    for (const std::array<double, 2> interval : {std::array<double, 2>{0.0, 1.2}, std::array<double, 2>{1.2, 2.0}})
        elements.push_back(append_cuboid(nodes, secondary_nodes, 1.0, 2.0, interval[0], interval[1], 0.0, 1.0));
    std::vector<std::int64_t> element_blocks(6, 1);
    element_blocks[4] = 2;
    element_blocks[5] = 2;
    const auto sides = [](std::initializer_list<fuelsim::ElementSide> values) {
        return std::vector<fuelsim::ElementSide>(values);
    };
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes), std::move(elements), std::move(element_blocks),
        {{1, "primary"}, {2, "secondary"}}, {},
        {{10, "primary_contact", sides({{0, 1}, {1, 1}, {2, 1}, {3, 1}})},
            {11, "primary_x0", sides({{0, 3}, {1, 3}, {2, 3}, {3, 3}})}, {12, "primary_y0", sides({{0, 0}})},
            {13, "primary_z0", sides({{0, 4}, {1, 4}, {2, 4}, {3, 4}})},
            {20, "secondary_contact", sides({{4, 3}, {5, 3}})}, {21, "secondary_x2", sides({{4, 1}, {5, 1}})},
            {22, "secondary_y0", sides({{4, 0}})}, {23, "secondary_z0", sides({{4, 4}, {5, 4}})}});
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition value;
    value.regions = {fuelsim::RegionDefinition{"primary", "primary", material(), 0.0, 300.0},
        fuelsim::RegionDefinition{"secondary", "secondary", material(), 0.0, 300.0}};
    value.boundary_conditions = {
        {"primary_x", fuelsim::BoundaryConditionType::dirichlet, "primary_x0", fuelsim::Field::displacement_x, 0.0},
        {"primary_y", fuelsim::BoundaryConditionType::dirichlet, "primary_y0", fuelsim::Field::displacement_y, 0.0},
        {"primary_z", fuelsim::BoundaryConditionType::dirichlet, "primary_z0", fuelsim::Field::displacement_z, 0.0},
        {"secondary_x", fuelsim::BoundaryConditionType::dirichlet, "secondary_x2", fuelsim::Field::displacement_x,
            -1.0e-5},
        {"secondary_y", fuelsim::BoundaryConditionType::dirichlet, "secondary_y0", fuelsim::Field::displacement_y, 0.0},
        {"secondary_z", fuelsim::BoundaryConditionType::dirichlet, "secondary_z0", fuelsim::Field::displacement_z, 0.0},
        {"primary_temperature", fuelsim::BoundaryConditionType::dirichlet, "primary_z0", fuelsim::Field::temperature,
            300.0},
        {"secondary_temperature", fuelsim::BoundaryConditionType::dirichlet, "secondary_z0",
            fuelsim::Field::temperature, 300.0}};
    fuelsim::ContactDefinition contact;
    contact.name = "nonmatching_surface_contact";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.thermal = false;
    contact.mechanical = true;
    contact.penalty = 1.0e11;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    value.contacts.push_back(contact);
    return value;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions value;
    value.absolute_tolerance = 1.0e-8;
    value.relative_tolerance = 1.0e-11;
    value.maximum_iterations = 30;
    value.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    value.direct_factorization = fuelsim::SolverOptions::DirectFactorization::mumps;
    value.field_residual_scaling = true;
    value.linear_relative_tolerance = 1.0e-11;
    value.maximum_linear_iterations = 400;
    return value;
}
} // namespace

int main(int argc, char** argv) {
    fuelsim::PetscSession session(argc, argv, "fuelsim HEX20 nonmatching contact tests\n");
    const fuelsim::UnstructuredHex20Mesh mesh = nonmatching_mesh();
    fuelsim::SteadyProblem problem(definition(), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::array<double, 4> primary_midpoints = {0.25, 0.75, 1.25, 1.75};
    std::array<std::size_t, 4> primary_midpoint_dofs{};
    for (std::size_t segment = 0; segment < primary_midpoints.size(); ++segment) {
        const auto& region = spatial.hex20_region_mesh(0);
        const auto found = std::find_if(region.nodes().begin(), region.nodes().end(), [&](const auto& point) {
            return std::abs(point.x - 1.0) < 1.0e-14 && std::abs(point.y - primary_midpoints[segment]) < 1.0e-14 &&
                   std::abs(point.z) < 1.0e-14;
        });
        if (found == region.nodes().end()) {
            std::cerr << "[FAIL] missing primary face midpoint node\n";
            return 1;
        }
        primary_midpoint_dofs[segment] = spatial.dof(fuelsim::Field::displacement_x,
            spatial.global_node(0, static_cast<std::size_t>(found - region.nodes().begin())));
    }
    const auto owner_counts = [&](const std::vector<double>& state) {
        problem.validate_state(state);
        std::array<std::size_t, 4> counts{};
        for (std::size_t contribution = spatial.volume_contribution_count();
            contribution < spatial.contribution_count(); ++contribution) {
            std::vector<std::size_t> dofs;
            problem.contribution_dofs(contribution, dofs);
            const auto owner = std::find_if(primary_midpoint_dofs.begin(), primary_midpoint_dofs.end(),
                [&](std::size_t dof) { return std::find(dofs.begin(), dofs.end(), dof) != dofs.end(); });
            if (owner == primary_midpoint_dofs.end()) {
                std::cerr << "[FAIL] nonmatching contact contribution has no unique primary owner\n";
                return std::array<std::size_t, 4>{};
            }
            ++counts[static_cast<std::size_t>(owner - primary_midpoint_dofs.begin())];
        }
        return counts;
    };
    const std::array<std::size_t, 4> owners = owner_counts(problem.initial_state());
    bool passed = check(owners == std::array<std::size_t, 4>{3, 3, 6, 6},
        "nonmatching HEX20 surface integration points select all four primary faces exactly once");
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {1, 0.5, 4, 1.0e-6}, solver_options());
    passed = check(result.completed && result.solve.converged && result.completed_steps == 1,
                 "nonmatching HEX20 surface contact converges under compression") &&
             passed;
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, result.solve.state);
    std::array<double, 3> contact_balance{};
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        std::vector<double> residual;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local) local_state[local] = result.solve.state[dofs[local]];
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, nullptr);
        for (std::size_t local = 0; local < dofs.size(); ++local) {
            const std::size_t global = dofs[local];
            for (std::size_t component = 0; component < contact_balance.size(); ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (global >= field.begin && global < field.end) contact_balance[component] += residual[local];
            }
        }
    }
    passed = check(interface.projected_contact_nodes > 0 && interface.unprojected_contact_nodes == 0,
                 "nonmatching HEX20 contact projects every secondary boundary node") &&
             check(interface.active_contact_nodes > 0 && interface.total_contact_force > 0.0,
                 "nonmatching HEX20 contact develops positive compressive interface force") &&
             check(std::abs(contact_balance[0]) < 1.0e-8 && std::abs(contact_balance[1]) < 1.0e-8 &&
                       std::abs(contact_balance[2]) < 1.0e-8,
                 "nonmatching HEX20 contact residual is action-reaction conservative") &&
             passed;
    const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    passed = check(histories.size() == 18, "nonmatching HEX20 contact stores one committed history for each of "
                                           "eighteen secondary quadrature points") &&
             passed;
    if (session.rank() == 0) {
        std::cout << "h20_24_primary_face_owner_counts=" << owners[0] << ',' << owners[1] << ',' << owners[2] << ','
                  << owners[3] << '\n'
                  << "h20_24_projected_secondary_nodes=" << interface.projected_contact_nodes << '\n'
                  << "h20_24_active_contact_nodes=" << interface.active_contact_nodes << '\n'
                  << "h20_24_total_contact_force=" << interface.total_contact_force << '\n';
    }
    return passed ? 0 : 1;
}
