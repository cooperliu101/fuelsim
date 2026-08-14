#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"
#include "support/jacobian_check.hpp"
#include "support/material_factory.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
fuelsim::UnstructuredQuad4Mesh three_region_mesh() {
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {0.0, 0.0},
            {1.0, 0.0},
            {1.0, 1.0},
            {0.0, 1.0},
            {1.1, 0.0},
            {1.2, 0.0},
            {1.2, 1.0},
            {1.1, 1.0},
            {1.3, 0.0},
            {1.4, 0.0},
            {1.4, 1.0},
            {1.3, 1.0},
        },
        {
            {{{0, 1, 2, 3}}},
            {{{4, 5, 6, 7}}},
            {{{8, 9, 10, 11}}},
        },
        {1, 2, 3},
        {
            {1, "pellet"},
            {2, "clad_1"},
            {3, "clad_2"},
        },
        {},
        {
            {11, "pellet_axis", {{{0, 3}}}},
            {12, "pellet_outer", {{{0, 1}}}},
            {13, "pellet_bottom", {{{0, 0}}}},
            {14, "pellet_top", {{{0, 2}}}},
            {21, "clad_1_inner", {{{1, 3}}}},
            {22, "clad_1_outer", {{{1, 1}}}},
            {23, "clad_1_bottom", {{{1, 0}}}},
            {24, "clad_1_top", {{{1, 2}}}},
            {31, "clad_2_inner", {{{2, 3}}}},
            {32, "clad_2_outer", {{{2, 1}}}},
            {33, "clad_2_bottom", {{{2, 0}}}},
            {34, "clad_2_top", {{{2, 2}}}},
        });
}
fuelsim::UnstructuredQuad4Mesh annular_boundary_mesh() {
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {1.0, 0.0},
            {2.0, 0.0},
            {2.0, 1.0},
            {1.0, 1.0},
        },
        {{{{0, 1, 2, 3}}}}, {1}, {{1, "solid"}}, {},
        {
            {1, "bottom", {{{0, 0}}}},
            {2, "right", {{{0, 1}}}},
            {3, "top", {{{0, 2}}}},
            {4, "left", {{{0, 3}}}},
        });
}
fuelsim::UnstructuredQuad4Mesh two_pellet_nonmatching_mesh() {
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {0.0, 0.0},
            {0.002, 0.0},
            {0.004, 0.0},
            {0.0, 0.001},
            {0.002, 0.001},
            {0.004, 0.001},
            {0.0, 0.001002},
            {0.001, 0.001002},
            {0.0025, 0.001002},
            {0.0041, 0.001002},
            {0.0, 0.002},
            {0.001, 0.002},
            {0.0025, 0.002},
            {0.0041, 0.002},
        },
        {
            {{{0, 1, 4, 3}}},
            {{{1, 2, 5, 4}}},
            {{{6, 7, 11, 10}}},
            {{{7, 8, 12, 11}}},
            {{{8, 9, 13, 12}}},
        },
        {10, 10, 20, 20, 20}, {{10, "lower_pellet"}, {20, "upper_pellet"}}, {},
        {
            {101, "lower_axis", {{{0, 3}}}},
            {102, "lower_bottom", {{{0, 0}, {1, 0}}}},
            {103, "lower_top", {{{0, 2}, {1, 2}}}},
            {104, "lower_outer", {{{1, 1}}}},
            {201, "upper_axis", {{{2, 3}}}},
            {202, "upper_bottom", {{{2, 0}, {3, 0}, {4, 0}}}},
            {203, "upper_top", {{{2, 2}, {3, 2}, {4, 2}}}},
            {204, "upper_outer", {{{4, 1}}}},
        });
}
fuelsim::UnstructuredQuad4Mesh l_shaped_primary_mesh() {
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {2.0, 0.0},
            {3.0, 0.0},
            {3.0, 1.0},
            {2.0, 1.0},
            {1.0, 0.5},
            {1.5, 0.5},
            {1.5, 1.0},
            {1.0, 1.0},
        },
        {
            {{{0, 1, 2, 3}}},
            {{{4, 5, 6, 7}}},
        },
        {10, 20}, {{10, "tool"}, {20, "slug"}}, {},
        {
            {11, "tool_corner", {{{0, 1}, {0, 2}}}},
            {22, "slug_face", {{{1, 1}}}},
        });
}
fuelsim::UnstructuredQuad4Mesh coincident_fuel_clad_mesh(double clad_inner_radius) {
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {0.0, 0.0},
            {1.0, 0.0},
            {1.0, 1.0},
            {0.0, 1.0},
            {clad_inner_radius, 0.0},
            {1.2, 0.0},
            {1.2, 1.0},
            {clad_inner_radius, 1.0},
        },
        {
            {{{0, 1, 2, 3}}},
            {{{4, 5, 6, 7}}},
        },
        {1, 2}, {{1, "fuel"}, {2, "clad"}}, {},
        {
            {11, "fuel_axis", {{{0, 3}}}},
            {12, "fuel_outer", {{{0, 1}}}},
            {13, "fuel_bottom", {{{0, 0}}}},
            {21, "clad_inner", {{{1, 3}}}},
            {22, "clad_outer", {{{1, 1}}}},
            {23, "clad_bottom", {{{1, 0}}}},
        });
}
fuelsim::UnstructuredQuad4Mesh overlapping_material_mesh() {
    // Both blocks occupy the lower half of the tall block, so the two
    // materials sit on the same side of their shared faces.
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {0.0, 0.0},
            {1.0, 0.0},
            {1.0, 1.0},
            {0.0, 1.0},
            {0.0, 0.0},
            {1.0, 0.0},
            {1.0, 0.5},
            {0.0, 0.5},
        },
        {
            {{{0, 1, 2, 3}}},
            {{{4, 5, 6, 7}}},
        },
        {10, 20}, {{10, "tall"}, {20, "short"}}, {},
        {
            {11, "tall_bottom", {{{0, 0}}}},
            {12, "tall_top", {{{0, 2}}}},
            {21, "short_bottom", {{{1, 0}}}},
        });
}
fuelsim::ThermoelasticProperties thermoelastic(double conductivity) {
    return fuelsim::test::thermoelastic(0.0, conductivity, 2.0e11, 0.3, 1.0e-5, 300.0, 0.0, 0.0, 0.0, 10.0, 20.0);
}
fuelsim::RegionDefinition region(
    const std::string& name, const std::string& block, double initial_temperature, double heat_source) {
    return {name, block, thermoelastic(10.0), heat_source, initial_temperature};
}
fuelsim::ContactDefinition contact(const std::string& name, const std::string& primary, const std::string& secondary) {
    return {name, primary, secondary, true, true, 0.2, 1.0e-5, 1.0e14};
}
fuelsim::BoundaryConditionDefinition dirichlet(
    const std::string& name, const std::string& boundary, fuelsim::Field field, double value) {
    return {name, fuelsim::BoundaryConditionType::dirichlet, boundary, field, value};
}
fuelsim::SpatialDefinition single_region_definition() {
    return {{region("pellet", "pellet", 500.0, 2.0e5)}, {},
        {
            dirichlet("axis", "pellet_axis", fuelsim::Field::radial_displacement, 0.0),
            dirichlet("bottom", "pellet_bottom", fuelsim::Field::axial_displacement, 0.0),
            dirichlet("outer_temperature", "pellet_outer", fuelsim::Field::temperature, 300.0),
        }};
}
fuelsim::SpatialDefinition three_region_definition() {
    return {{region("pellet", "pellet", 500.0, 2.0e5), region("inner_clad", "clad_1", 400.0, 0.0),
                region("outer_clad", "clad_2", 300.0, 0.0)},
        {contact("pellet_to_inner", "clad_1_inner", "pellet_outer"),
            contact("inner_to_outer", "clad_2_inner", "clad_1_outer")},
        {
            dirichlet("axis", "pellet_axis", fuelsim::Field::radial_displacement, 0.0),
            dirichlet("pellet_bottom", "pellet_bottom", fuelsim::Field::axial_displacement, 0.0),
            dirichlet("inner_bottom", "clad_1_bottom", fuelsim::Field::axial_displacement, 0.0),
            dirichlet("outer_bottom", "clad_2_bottom", fuelsim::Field::axial_displacement, 0.0),
            dirichlet("outer_temperature", "clad_2_outer", fuelsim::Field::temperature, 300.0),
        }};
}
bool test_single_region(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    const fuelsim::SteadyProblem problem(single_region_definition(), mesh);
    const std::vector<double> state = problem.initial_state();
    bool passed =
        check(fuelsim::rz::ProblemAccess::region_count(problem) == 1 &&
                  fuelsim::rz::ProblemAccess::contact_count(problem) == 0,
            "one Exodus block forms a standalone steady problem") &&
        check(fuelsim::rz::ProblemAccess::region_index(problem, "pellet") == 0 &&
                  fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) == 0 &&
                  fuelsim::rz::ProblemAccess::region_element_offset(problem, 0) == 0,
            "single-region indices start at zero") &&
        check(problem.dof_count() == 12 && fuelsim::rz::ProblemAccess::volume_contribution_count(problem) == 1 &&
                  problem.contribution_count() == 1,
            "single Quad4 region has one 12-DOF contribution") &&
        check(state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, 0)] == 500.0 &&
                  state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, 1)] == 300.0,
            "region initial temperature and boundary value are applied");
    const std::vector<fuelsim::FieldDescriptor>& fields = problem.field_layout();
    const fuelsim::LocalDofs rz_dofs = fuelsim::rz::ProblemAccess::contribution_dofs(problem, 0);
    bool rz_layout = fields.size() == 3 && fields[0].name == "temperature" && fields[0].begin == 0 &&
                     fields[0].end == fuelsim::rz::ProblemAccess::dof_map(problem).node_count() &&
                     fields[0].category == fuelsim::FieldCategory::thermal && fields[1].name == "radial" &&
                     fields[1].begin == fuelsim::rz::ProblemAccess::dof_map(problem).node_count() &&
                     fields[1].end == 2 * fuelsim::rz::ProblemAccess::dof_map(problem).node_count() &&
                     fields[1].category == fuelsim::FieldCategory::mechanical && fields[2].name == "axial" &&
                     fields[2].begin == 2 * fuelsim::rz::ProblemAccess::dof_map(problem).node_count() &&
                     fields[2].end == problem.dof_count() && fields[2].category == fuelsim::FieldCategory::mechanical &&
                     rz_dofs.size() == fuelsim::local_dof_count;
    const fuelsim::Quad4Element& element = fuelsim::rz::ProblemAccess::region_mesh(problem, 0).elements().front();
    for (std::size_t local_node = 0; local_node < 4; ++local_node) {
        const std::size_t global_node =
            fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) + element.nodes[local_node];
        rz_layout = rz_layout &&
                    rz_dofs[local_node] ==
                        fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, global_node) &&
                    rz_dofs[4 + local_node] == fuelsim::rz::ProblemAccess::dof_map(problem).dof(
                                                   fuelsim::Field::radial_displacement, global_node) &&
                    rz_dofs[8 + local_node] == fuelsim::rz::ProblemAccess::dof_map(problem).dof(
                                                   fuelsim::Field::axial_displacement, global_node);
    }
    passed = check(rz_layout, "RZ adapter preserves [T0..T3, ur0..ur3, uz0..uz3] and field metadata") && passed;
    const fuelsim::LocalValues local = fuelsim::rz::ProblemAccess::contribution_state(problem, 0, state);
    const fuelsim::rz::LocalLinearization system =
        fuelsim::rz::ProblemAccess::linearize_contribution(problem, 0, local);
    passed = check(std::all_of(system.residual.begin(), system.residual.end(),
                       [](double value) { return std::isfinite(value); }),
                 "single-region AD residual is finite") &&
             passed;
    return passed;
}
bool test_time_controlled_pressure(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SpatialDefinition definition = single_region_definition();
    definition.time_tables.emplace_back(
        "pressure_history", std::vector<double>{0.0, 1.0}, std::vector<double>{1.0, 2.0});
    fuelsim::BoundaryConditionDefinition pressure{"outer_pressure", fuelsim::BoundaryConditionType::pressure,
        "pellet_outer", fuelsim::Field::radial_displacement, 10.0, false};
    pressure.function = "pressure_history";
    definition.boundary_conditions.push_back(std::move(pressure));
    fuelsim::SteadyProblem problem(std::move(definition), mesh);
    problem.set_time(0.0);
    const std::size_t pressure_contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
    const fuelsim::LocalValues local =
        fuelsim::rz::ProblemAccess::contribution_state(problem, pressure_contribution, problem.initial_state());
    const fuelsim::LocalResidual first =
        fuelsim::rz::ProblemAccess::contribution_residual(problem, pressure_contribution, local);
    problem.set_time(1.0);
    const fuelsim::LocalResidual second =
        fuelsim::rz::ProblemAccess::contribution_residual(problem, pressure_contribution, local);
    bool passed = true;
    for (std::size_t dof = 0; dof < first.size(); ++dof)
        passed = check(std::abs(second[dof] - 2.0 * first[dof]) < 1.0e-12,
                     "pressure time table scales the assembled load") &&
                 passed;
    return passed;
}
bool test_pressure_parent_edge_orientation() {
    const fuelsim::UnstructuredQuad4Mesh mesh = annular_boundary_mesh();
    constexpr double pressure_value = 3.0;
    constexpr double pi = 3.141592653589793238462643383279502884;
    const auto resultant = [&](const std::string& boundary) {
        fuelsim::SpatialDefinition definition = {{region("solid", "solid", 600.0, 0.0)}, {}, {}};
        definition.boundary_conditions.push_back({"pressure", fuelsim::BoundaryConditionType::pressure, boundary,
            fuelsim::Field::radial_displacement, pressure_value});
        const fuelsim::SteadyProblem problem(std::move(definition), mesh);
        const std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
        const fuelsim::LocalValues local =
            fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, problem.initial_state());
        const fuelsim::LocalResidual residual =
            fuelsim::rz::ProblemAccess::contribution_residual(problem, contribution, local);
        std::array<double, 2> force = {0.0, 0.0};
        for (std::size_t node = 0; node < 4; ++node) {
            force[0] += residual[4 + node];
            force[1] += residual[8 + node];
        }
        return force;
    };
    const std::array<double, 2> left = resultant("left");
    const std::array<double, 2> right = resultant("right");
    const std::array<double, 2> bottom = resultant("bottom");
    const std::array<double, 2> top = resultant("top");
    const double left_expected = -2.0 * pi * pressure_value;
    const double right_expected = 4.0 * pi * pressure_value;
    const double bottom_expected = -3.0 * pi * pressure_value;
    const double top_expected = 3.0 * pi * pressure_value;
    const double tolerance = 1.0e-12;
    std::cout << "pressure_resultants_left=" << left[0] << ',' << left[1] << '\n';
    std::cout << "pressure_resultants_right=" << right[0] << ',' << right[1] << '\n';
    std::cout << "pressure_resultants_bottom=" << bottom[0] << ',' << bottom[1] << '\n';
    std::cout << "pressure_resultants_top=" << top[0] << ',' << top[1] << '\n';
    return check(std::abs(left[0] - left_expected) < tolerance && std::abs(left[1]) < tolerance &&
                     std::abs(right[0] - right_expected) < tolerance && std::abs(right[1]) < tolerance &&
                     std::abs(bottom[0]) < tolerance && std::abs(bottom[1] - bottom_expected) < tolerance &&
                     std::abs(top[0]) < tolerance && std::abs(top[1] - top_expected) < tolerance,
        "pressure uses the parent Quad4 outward normal on left, "
        "right, bottom, and top boundaries");
}
bool test_global_field_diagnostics(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SteadyProblem problem(single_region_definition(), mesh);
    const std::vector<double> state = problem.initial_state();
    std::vector<double> direction(problem.dof_count(), 0.0);
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::dof_map(problem).node_count(); ++node) {
        direction[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, node)] = 0.2;
        direction[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, node)] = 1.0e-6;
        direction[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, node)] = -0.7e-6;
    }
    const fuelsim::test::DirectionalJacobianCheck diagnostic =
        fuelsim::test::check_directional_jacobian(problem, state, direction, 1.0e-4);
    bool passed = true;
    for (std::size_t field = 0; field < 3; ++field) {
        const double reference = diagnostic.finite_difference_directional_derivative.l2[field];
        passed = check(diagnostic.difference.l2[field] <= 1.0e-7 * (1.0 + reference),
                     "global field Jacobian matches centered differences") &&
                 passed;
    }
    return passed;
}
bool test_three_regions(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SteadyProblem problem(three_region_definition(), mesh);
    const std::vector<double> state = problem.initial_state();
    bool passed =
        check(mesh.side_set_block_id("pellet_outer") == 1 && mesh.side_set_block_id("clad_1_inner") == 2 &&
                  mesh.side_set_block_id("clad_2_inner") == 3,
            "side-set ownership is inferred from adjacent Exodus elements") &&
        check(fuelsim::rz::ProblemAccess::region_count(problem) == 3 &&
                  fuelsim::rz::ProblemAccess::contact_count(problem) == 2,
            "three regions and two named contact pairs are composed") &&
        check(fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) == 0 &&
                  fuelsim::rz::ProblemAccess::region_node_offset(problem, 1) == 4 &&
                  fuelsim::rz::ProblemAccess::region_node_offset(problem, 2) == 8,
            "arbitrary regions receive independent node ranges") &&
        check(fuelsim::rz::ProblemAccess::region_element_offset(problem, 0) == 0 &&
                  fuelsim::rz::ProblemAccess::region_element_offset(problem, 1) == 1 &&
                  fuelsim::rz::ProblemAccess::region_element_offset(problem, 2) == 2,
            "arbitrary regions receive independent element ranges") &&
        check(problem.dof_count() == 36 && fuelsim::rz::ProblemAccess::volume_contribution_count(problem) == 3 &&
                  problem.contribution_count() == 11,
            "three volumes, four STS integration points, and four NTS nodes are assembled") &&
        check(fuelsim::rz::ProblemAccess::contact(problem, 0).primary == "clad_1_inner" &&
                  fuelsim::rz::ProblemAccess::contact(problem, 0).secondary == "pellet_outer" &&
                  fuelsim::rz::ProblemAccess::contact(problem, 1).primary == "clad_2_inner" &&
                  fuelsim::rz::ProblemAccess::contact(problem, 1).secondary == "clad_1_outer",
            "contacts contain only primary and secondary side-set names");
    for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
        contribution < problem.contribution_count(); ++contribution) {
        const fuelsim::LocalValues local = fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, state);
        const fuelsim::rz::LocalLinearization system =
            fuelsim::rz::ProblemAccess::linearize_contribution(problem, contribution, local);
        passed = check(std::all_of(system.residual.begin(), system.residual.end(),
                           [](double value) { return std::isfinite(value); }),
                     "multi-contact AD residual is finite") &&
                 passed;
    }
    fuelsim::SpatialDefinition solve_definition = three_region_definition();
    for (fuelsim::RegionDefinition& region_value : solve_definition.regions) {
        region_value.initial_temperature = 300.0;
        region_value.volumetric_heat_source = 0.0;
    }
    for (fuelsim::ContactDefinition& contact_value : solve_definition.contacts) {
        contact_value.penalty = 1.0e6;
        // This subcase isolates simultaneous mechanical activation. Its large
        // prescribed radial closure gives the three one-element bodies very
        // different Poisson axial contractions, so their thermal faces no
        // longer overlap in the converged configuration.
        contact_value.thermal = false;
    }
    solve_definition.boundary_conditions.push_back(
        dirichlet("pellet_temperature", "pellet_axis", fuelsim::Field::temperature, 300.0));
    solve_definition.boundary_conditions.push_back(
        dirichlet("inner_clad_temperature", "clad_1_inner", fuelsim::Field::temperature, 300.0));
    fuelsim::BoundaryConditionDefinition pellet_closure =
        dirichlet("pellet_contact_closure", "pellet_outer", fuelsim::Field::radial_displacement, 0.101);
    pellet_closure.scale_with_load = true;
    solve_definition.boundary_conditions.push_back(pellet_closure);
    fuelsim::BoundaryConditionDefinition inner_clad_closure =
        dirichlet("inner_clad_contact_closure", "clad_1_outer", fuelsim::Field::radial_displacement, 0.101);
    inner_clad_closure.scale_with_load = true;
    solve_definition.boundary_conditions.push_back(inner_clad_closure);
    solve_definition.boundary_conditions.push_back(
        dirichlet("inner_clad_contact_anchor", "clad_1_inner", fuelsim::Field::radial_displacement, 0.0));
    solve_definition.boundary_conditions.push_back(
        dirichlet("outer_clad_contact_anchor", "clad_2_inner", fuelsim::Field::radial_displacement, 0.0));
    fuelsim::SteadyProblem solve_problem(std::move(solve_definition), mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(solve_problem, {4, 0.5, 12, 1.0e-6}, fuelsim::SolverOptions{});
    if (!solve.completed)
        std::cerr << "multi-contact solve failure: "
                  << fuelsim::solve_failure_category_name(solve.solve.failure_category) << ": "
                  << solve.solve.failure_message << '\n';
    passed =
        check(solve.completed && solve.solve.converged, "three-region two-contact PETSc solve converges") && passed;
    if (solve.completed && solve.solve.converged) {
        for (std::size_t contact_value = 0; contact_value < fuelsim::rz::ProblemAccess::contact_count(solve_problem);
            ++contact_value) {
            const fuelsim::InterfaceSummary summary =
                fuelsim::rz::ProblemAccess::summarize_interface(solve_problem, contact_value, solve.solve.state);
            passed = check(summary.projected_contact_nodes > 0 && summary.active_contact_nodes > 0 &&
                               summary.unprojected_contact_nodes == 0,
                         "every converged contact pair remains projected "
                         "and mechanically active") &&
                     passed;
        }
    }
    return passed;
}
bool test_nonmatching_pellet_faces() {
    const fuelsim::UnstructuredQuad4Mesh mesh = two_pellet_nonmatching_mesh();
    fuelsim::SpatialDefinition definition = {
        {region("lower", "lower_pellet", 700.0, 0.0), region("upper", "upper_pellet", 500.0, 0.0)},
        {contact("pellet_stack", "upper_bottom", "lower_top")}, {}};
    const fuelsim::SteadyProblem problem(std::move(definition), mesh);
    const std::vector<double> state = problem.initial_state();
    bool passed = check(fuelsim::rz::ProblemAccess::region_count(problem) == 2 &&
                            fuelsim::rz::ProblemAccess::contact_count(problem) == 1,
        "two pellet blocks form one axial contact pair");
    for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
        contribution < problem.contribution_count(); ++contribution) {
        const fuelsim::LocalResidual residual = fuelsim::rz::ProblemAccess::contribution_residual(
            problem, contribution, fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, state));
        passed = check(std::all_of(residual.begin(), residual.end(), [](double value) { return std::isfinite(value); }),
                     "nonmatching pellet-face contribution is finite") &&
                 passed;
    }
    const fuelsim::InterfaceSummary summary = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
    const double expected_area = 3.141592653589793238462643383279502884 * 0.004 * 0.004;
    const double expected_heat_rate = expected_area * (0.2 / 1.0e-5) * (700.0 - 500.0);
    std::cout << "nonmatching_heat_rate=" << summary.total_heat_rate << " expected=" << expected_heat_rate << '\n';
    passed = check(std::abs(summary.total_heat_rate - expected_heat_rate) < 1.0e-12 * expected_heat_rate,
                 "split nonmatching STS integration covers the full pellet face") &&
             check(summary.projected_contact_nodes == 3 && summary.active_contact_nodes == 0,
                 "all nonmatching pellet-face NTS nodes project uniquely") &&
             passed;
    return passed;
}
bool test_l_shaped_primary_collinear_candidate() {
    // The ordered primary chain bends at (3,1): segment (2,1)->(3,1) then
    // (3,1)->(3,0). The secondary node (1.5,1.0) lies exactly on the
    // extension of the first segment but its reference projection fraction
    // is -0.5, far outside the segment. Construction must assign that
    // candidate a deterministic normal orientation instead of rejecting
    // the whole problem with a zero reference-gap error.
    const fuelsim::UnstructuredQuad4Mesh mesh = l_shaped_primary_mesh();
    const fuelsim::SpatialDefinition definition = {
        {region("tool", "tool", 400.0, 0.0), region("slug", "slug", 500.0, 0.0)},
        {{"corner_contact", "tool_corner", "slug_face", false, true, 0.2, 1.0e-5, 1.0e14}}, {}};
    try {
        const fuelsim::SteadyProblem problem(definition, mesh);
        const std::vector<double> state = problem.initial_state();
        bool passed =
            check(fuelsim::rz::ProblemAccess::region_count(problem) == 2 &&
                      fuelsim::rz::ProblemAccess::contact_count(problem) == 1,
                "L-shaped primary chain with a far collinear candidate "
                "constructs successfully") &&
            check(problem.contribution_count() == fuelsim::rz::ProblemAccess::volume_contribution_count(problem) + 4,
                "two secondary nodes times two primary segments form four "
                "mechanical candidates");
        for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
            contribution < problem.contribution_count(); ++contribution) {
            const fuelsim::LocalResidual residual = fuelsim::rz::ProblemAccess::contribution_residual(
                problem, contribution, fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, state));
            passed =
                check(std::all_of(residual.begin(), residual.end(), [](double value) { return std::isfinite(value); }),
                    "L-shaped primary candidate residual is finite") &&
                passed;
        }
        const fuelsim::InterfaceSummary summary = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
        passed =
            check(summary.active_contact_nodes == 0, "far collinear candidate keeps the open gap inactive") && passed;
        return passed;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] L-shaped primary chain construction raised: " << error.what() << '\n';
        return false;
    }
}
bool test_zero_initial_gap_construction() {
    // Coincident fuel and cladding surfaces (zero reference normal gap) are
    // a legal initial condition: the normal orientation comes from the
    // material side of each boundary edge's parent element.
    const fuelsim::UnstructuredQuad4Mesh mesh = coincident_fuel_clad_mesh(1.0);
    const fuelsim::SpatialDefinition definition = {
        {region("fuel", "fuel", 500.0, 0.0), region("clad", "clad", 300.0, 0.0)},
        {contact("fuel_clad", "clad_inner", "fuel_outer")}, {}};
    try {
        const fuelsim::SteadyProblem problem(definition, mesh);
        const std::vector<double> state = problem.initial_state();
        bool passed =
            check(fuelsim::rz::ProblemAccess::region_count(problem) == 2 &&
                      fuelsim::rz::ProblemAccess::contact_count(problem) == 1,
                "coincident fuel-cladding surfaces construct one contact") &&
            check(problem.contribution_count() == fuelsim::rz::ProblemAccess::volume_contribution_count(problem) + 4,
                "zero-gap contact assembles two STS integration points and two NTS candidates");
        for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
            contribution < problem.contribution_count(); ++contribution) {
            const fuelsim::rz::LocalLinearization system = fuelsim::rz::ProblemAccess::linearize_contribution(
                problem, contribution, fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, state));
            passed = check(std::all_of(system.residual.begin(), system.residual.end(),
                               [](double value) { return std::isfinite(value); }),
                         "zero-gap contact contribution is finite at the "
                         "initial state") &&
                     passed;
        }
        const fuelsim::InterfaceSummary summary = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
        passed = check(summary.projected_contact_nodes == 2 && summary.unprojected_contact_nodes == 0 &&
                           summary.active_contact_nodes == 0,
                     "zero-gap secondary nodes project at exactly zero gap "
                     "with no pressure") &&
                 passed;
        return passed;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] zero initial gap construction raised: " << error.what() << '\n';
        return false;
    }
}
bool test_overlapping_material_rejected() {
    // The two blocks occupy the same space next to the interface, so both
    // parent-element centroids lie on the same side of the primary line and
    // no meaningful zero-gap orientation exists.
    const fuelsim::UnstructuredQuad4Mesh mesh = overlapping_material_mesh();
    const fuelsim::SpatialDefinition degenerate = {
        {region("tall", "tall", 400.0, 0.0), region("short", "short", 400.0, 0.0)},
        {{"overlap", "tall_bottom", "short_bottom", false, true, 0.2, 1.0e-5, 1.0e14}}, {}};
    bool rejected = false;
    try {
        const fuelsim::SteadyProblem problem(degenerate, mesh);
        (void)problem;
    } catch (const std::invalid_argument&) { rejected = true; }
    bool passed = check(rejected, "coincident surfaces with both materials on the same "
                                  "side remain an explicit error");
    // Control: the same overlapping blocks with a positive 1 m gap between
    // the faces constructs normally, because the material-side check only
    // applies when a secondary node rides exactly on the segment.
    const fuelsim::SpatialDefinition open = {{region("tall", "tall", 400.0, 0.0), region("short", "short", 400.0, 0.0)},
        {{"overlap_open", "tall_top", "short_bottom", false, true, 0.2, 1.0e-5, 1.0e14}}, {}};
    try {
        const fuelsim::SteadyProblem problem(open, mesh);
        const std::vector<double> state = problem.initial_state();
        const fuelsim::InterfaceSummary summary = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
        passed = check(summary.projected_contact_nodes == 2 && summary.active_contact_nodes == 0,
                     "positive-gap overlapping blocks construct with the "
                     "material-side check dormant") &&
                 passed;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] positive-gap overlapping construction raised: " << error.what() << '\n';
        return false;
    }
    return passed;
}
fuelsim::SpatialDefinition zero_gap_solve_definition(double closure_displacement) {
    fuelsim::SpatialDefinition definition = {{region("fuel", "fuel", 500.0, 2.0e2), region("clad", "clad", 300.0, 0.0)},
        {contact("fuel_clad", "clad_inner", "fuel_outer")},
        {
            dirichlet("axis", "fuel_axis", fuelsim::Field::radial_displacement, 0.0),
            dirichlet("fuel_bottom", "fuel_bottom", fuelsim::Field::axial_displacement, 0.0),
            dirichlet("clad_bottom", "clad_bottom", fuelsim::Field::axial_displacement, 0.0),
            dirichlet("outer_temperature", "clad_outer", fuelsim::Field::temperature, 300.0),
            dirichlet("clad_inner_anchor", "clad_inner", fuelsim::Field::radial_displacement, 0.0),
        }};
    definition.contacts[0].penalty = 1.0e6;
    fuelsim::BoundaryConditionDefinition closure =
        dirichlet("closure", "fuel_outer", fuelsim::Field::radial_displacement, closure_displacement);
    closure.scale_with_load = true;
    definition.boundary_conditions.push_back(std::move(closure));
    return definition;
}
bool test_zero_initial_gap_solve() {
    // End-to-end oracle: a coincident-interface problem pressed 0.01 m into
    // contact must behave like the same problem with a 1e-9 m initial gap
    // pressed 0.01+1e-9 m, because both reach the same current interface
    // state from consistent normal orientations.
    constexpr double penetration = 0.01;
    constexpr double epsilon_gap = 1.0e-9;
    constexpr double expected_pressure = 1.0e6 * penetration;
    double pressures[2] = {0.0, 0.0};
    bool passed = true;
    for (std::size_t variant = 0; variant < 2; ++variant) {
        const double clad_inner = 1.0 + (variant == 0 ? 0.0 : epsilon_gap);
        const double closure_value = penetration + (variant == 0 ? 0.0 : epsilon_gap);
        const fuelsim::UnstructuredQuad4Mesh mesh = coincident_fuel_clad_mesh(clad_inner);
        fuelsim::SteadyProblem problem(zero_gap_solve_definition(closure_value), mesh);
        const fuelsim::SteadyResult solve =
            fuelsim::solve_steady(problem, {4, 0.5, 12, 1.0e-6}, fuelsim::SolverOptions{});
        if (!solve.completed)
            std::cerr << "zero initial gap solve failure (variant " << variant
                      << "): " << fuelsim::solve_failure_category_name(solve.solve.failure_category) << ": "
                      << solve.solve.failure_message << '\n';
        passed =
            check(solve.completed && solve.solve.converged,
                variant == 0 ? "coincident-interface PETSc solve converges" : "1e-9 m opened twin solve converges") &&
            passed;
        if (!solve.completed || !solve.solve.converged) return false;
        const fuelsim::InterfaceSummary summary =
            fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
        pressures[variant] = summary.maximum_contact_pressure;
        std::cout << (variant == 0 ? "zero_initial_gap" : "zero_initial_gap_opened_twin")
                  << "_maximum_contact_pressure=" << summary.maximum_contact_pressure << '\n'
                  << (variant == 0 ? "zero_initial_gap" : "zero_initial_gap_opened_twin")
                  << "_total_contact_force=" << summary.total_contact_force << '\n'
                  << (variant == 0 ? "zero_initial_gap" : "zero_initial_gap_opened_twin")
                  << "_total_heat_rate=" << summary.total_heat_rate << '\n';
        passed = check(summary.projected_contact_nodes == 2 && summary.unprojected_contact_nodes == 0 &&
                           summary.active_contact_nodes == 2,
                     "every zero-gap secondary node stays projected and "
                     "mechanically active under the 0.01 m closure") &&
                 check(std::abs(summary.maximum_contact_pressure - expected_pressure) < 1.0e-4 * expected_pressure,
                     "contact pressure matches the penalty times the 0.01 m "
                     "closure within the solver tolerance") &&
                 check(summary.total_contact_force > 0.0, "active zero-gap contact carries a positive total force") &&
                 passed;
        if (variant == 0) {
            // All fuel heat leaves through the interface: the source is
            // 2e2 W/m3 over the unit-length half-cross-section of radius 1 m.
            constexpr double expected_heat_rate = 2.0e2 * 3.141592653589793238462643383279502884;
            passed = check(std::abs(summary.total_heat_rate - expected_heat_rate) < 1.0e-3 * expected_heat_rate,
                         "zero-gap thermal contact conducts the full fuel "
                         "heat generation to the cladding") &&
                     passed;
        }
    }
    const double oracle_error = std::abs(pressures[0] - pressures[1]) / expected_pressure;
    std::cout << "zero_initial_gap_twin_pressure_relative_difference=" << oracle_error << '\n';
    passed = check(oracle_error < 1.0e-6, "coincident and 1e-9 m opened solves agree on the contact "
                                          "pressure") &&
             passed;
    return passed;
}
bool test_transient_regions(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SpatialDefinition definition;
    definition = three_region_definition();
    fuelsim::TransientProblem problem(std::move(definition), mesh);
    const std::vector<double> initial = problem.committed_solution();
    problem.begin_time_step({1.0, 0.5});
    const fuelsim::LocalValues local = fuelsim::rz::ProblemAccess::contribution_state(problem, 2, initial);
    (void)fuelsim::rz::ProblemAccess::linearize_contribution(problem, 2, local);
    problem.rollback_time_step();
    bool passed = check(fuelsim::rz::ProblemAccess::region_count(problem) == 3 && problem.committed_time() == 0.0 &&
                            problem.committed_solution() == initial && !problem.time_step_active(),
        "three-region transient rollback preserves committed state");
    problem.begin_time_step({1.0, 1.0});
    problem.commit_time_step(initial);
    passed =
        check(problem.committed_time() == 1.0 && problem.committed_load_factor() == 1.0 && !problem.time_step_active(),
            "three-region transient state commits once") &&
        passed;
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const fuelsim::RegionStateSummary summary =
            fuelsim::rz::ProblemAccess::summarize_region_history(problem, region);
        passed =
            check(summary.maximum_equivalent_plastic_strain == 0.0 && summary.maximum_equivalent_creep_strain == 0.0,
                "elastic region keeps zero inelastic history") &&
            passed;
    }
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "fuelsim multi-region contact solve tests\n");
        const fuelsim::UnstructuredQuad4Mesh mesh = three_region_mesh();
        const bool passed = test_single_region(mesh) && test_time_controlled_pressure(mesh) &&
                            test_pressure_parent_edge_orientation() && test_global_field_diagnostics(mesh) &&
                            test_three_regions(mesh) && test_nonmatching_pellet_faces() &&
                            test_l_shaped_primary_collinear_candidate() && test_zero_initial_gap_construction() &&
                            test_overlapping_material_rejected() && test_zero_initial_gap_solve() &&
                            test_transient_regions(mesh);
        if (!passed) return 1;
        std::cout << "[PASS] single- and multi-region problem tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] multi-region tests raised: " << error.what() << '\n';
        return 1;
    }
}
