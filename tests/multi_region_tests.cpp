#include "fuelsim/diagnostics.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
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
            {0.0, 0.0},         {0.002, 0.0},        {0.004, 0.0},
            {0.0, 0.001},       {0.002, 0.001},      {0.004, 0.001},
            {0.0, 0.001002},    {0.001, 0.001002},   {0.0025, 0.001002},
            {0.0041, 0.001002}, {0.0, 0.002},        {0.001, 0.002},
            {0.0025, 0.002},    {0.0041, 0.002},
        },
        {
            {{{0, 1, 4, 3}}},
            {{{1, 2, 5, 4}}},
            {{{6, 7, 11, 10}}},
            {{{7, 8, 12, 11}}},
            {{{8, 9, 13, 12}}},
        },
        {10, 10, 20, 20, 20},
        {{10, "lower_pellet"}, {20, "upper_pellet"}}, {},
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

fuelsim::ThermoelasticProperties thermoelastic(double conductivity) {
    return {0.0, conductivity, 2.0e11, 0.3, 1.0e-5, 300.0};
}

fuelsim::RegionDefinition region(const std::string& name,
                                 const std::string& block,
                                 double initial_temperature,
                                 double heat_source) {
    return {name, block, thermoelastic(10.0), heat_source, initial_temperature};
}

fuelsim::ContactDefinition contact(const std::string& name,
                                   const std::string& primary,
                                   const std::string& secondary) {
    return {name, primary, secondary, true, true, 0.2, 1.0e-5, 1.0e14};
}

fuelsim::BoundaryConditionDefinition dirichlet(const std::string& name,
                                               const std::string& boundary,
                                               fuelsim::Field field,
                                               double value) {
    return {name, fuelsim::BoundaryConditionType::dirichlet, boundary, field,
            value};
}

fuelsim::SteadyProblemDefinition single_region_definition() {
    return {{region("pellet", "pellet", 500.0, 2.0e5)},
            {},
            {
                dirichlet("axis", "pellet_axis",
                          fuelsim::Field::radial_displacement, 0.0),
                dirichlet("bottom", "pellet_bottom",
                          fuelsim::Field::axial_displacement, 0.0),
                dirichlet("outer_temperature", "pellet_outer",
                          fuelsim::Field::temperature, 300.0),
            }};
}

fuelsim::SteadyProblemDefinition three_region_definition() {
    return {{region("pellet", "pellet", 500.0, 2.0e5),
             region("inner_clad", "clad_1", 400.0, 0.0),
             region("outer_clad", "clad_2", 300.0, 0.0)},
            {contact("pellet_to_inner", "clad_1_inner", "pellet_outer"),
             contact("inner_to_outer", "clad_2_inner", "clad_1_outer")},
            {
                dirichlet("axis", "pellet_axis",
                          fuelsim::Field::radial_displacement, 0.0),
                dirichlet("pellet_bottom", "pellet_bottom",
                          fuelsim::Field::axial_displacement, 0.0),
                dirichlet("inner_bottom", "clad_1_bottom",
                          fuelsim::Field::axial_displacement, 0.0),
                dirichlet("outer_bottom", "clad_2_bottom",
                          fuelsim::Field::axial_displacement, 0.0),
                dirichlet("outer_temperature", "clad_2_outer",
                          fuelsim::Field::temperature, 300.0),
            }};
}

fuelsim::TransientInelasticProperties elastic_transient_material() {
    return {10.0,
            20.0,
            fuelsim::InelasticBehavior::elastic,
            {0.0, 1.0, 1.0},
            {1.0, 0.0}};
}

bool test_single_region(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    const fuelsim::SteadyProblem problem(single_region_definition(), mesh);
    const std::vector<double> state = problem.initial_state();
    bool passed =
        check(problem.region_count() == 1 && problem.contact_count() == 0,
              "one Exodus block forms a standalone steady problem") &&
        check(problem.region_index("pellet") == 0 &&
                  problem.region_node_offset(0) == 0 &&
                  problem.region_element_offset(0) == 0,
              "single-region indices start at zero") &&
        check(problem.dof_count() == 12 &&
                  problem.volume_contribution_count() == 1 &&
                  problem.contribution_count() == 1,
              "single Quad4 region has one 12-DOF contribution") &&
        check(state[problem.dof_map().temperature(0)] == 500.0 &&
                  state[problem.dof_map().temperature(1)] == 300.0,
              "region initial temperature and boundary value are applied");

    const fuelsim::LocalValues local = problem.contribution_state(0, state);
    const fuelsim::LocalSystem system =
        problem.linearize_contribution(0, local);
    passed =
        check(std::all_of(system.residual.begin(), system.residual.end(),
                          [](double value) { return std::isfinite(value); }),
              "single-region AD residual is finite") &&
        passed;
    return passed;
}

bool test_time_controlled_pressure(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SteadyProblemDefinition definition = single_region_definition();
    definition.time_tables.emplace_back("pressure_history",
                                        std::vector<double>{0.0, 1.0},
                                        std::vector<double>{1.0, 2.0});
    fuelsim::BoundaryConditionDefinition pressure{
        "outer_pressure",
        fuelsim::BoundaryConditionType::pressure,
        "pellet_outer",
        fuelsim::Field::radial_displacement,
        10.0,
        false};
    pressure.function = "pressure_history";
    definition.boundary_conditions.push_back(std::move(pressure));
    fuelsim::SteadyProblem problem(std::move(definition), mesh);
    problem.set_time(0.0);
    const std::size_t pressure_contribution =
        problem.volume_contribution_count();
    const fuelsim::LocalValues local = problem.contribution_state(
        pressure_contribution, problem.initial_state());
    const fuelsim::LocalResidual first = problem.contribution_residual(
        pressure_contribution, local);
    problem.set_time(1.0);
    const fuelsim::LocalResidual second = problem.contribution_residual(
        pressure_contribution, local);
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
        fuelsim::SteadyProblemDefinition definition = {
            {region("solid", "solid", 600.0, 0.0)}, {}, {}};
        definition.boundary_conditions.push_back(
            {"pressure", fuelsim::BoundaryConditionType::pressure, boundary,
             fuelsim::Field::radial_displacement, pressure_value});
        const fuelsim::SteadyProblem problem(std::move(definition), mesh);
        const std::size_t contribution = problem.volume_contribution_count();
        const fuelsim::LocalValues local = problem.contribution_state(
            contribution, problem.initial_state());
        const fuelsim::LocalResidual residual =
            problem.contribution_residual(contribution, local);
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
    std::cout << "pressure_resultants_left=" << left[0] << ',' << left[1]
              << '\n';
    std::cout << "pressure_resultants_right=" << right[0] << ',' << right[1]
              << '\n';
    std::cout << "pressure_resultants_bottom=" << bottom[0] << ','
              << bottom[1] << '\n';
    std::cout << "pressure_resultants_top=" << top[0] << ',' << top[1]
              << '\n';
    return check(std::abs(left[0] - left_expected) < tolerance &&
                     std::abs(left[1]) < tolerance &&
                     std::abs(right[0] - right_expected) < tolerance &&
                     std::abs(right[1]) < tolerance &&
                     std::abs(bottom[0]) < tolerance &&
                     std::abs(bottom[1] - bottom_expected) < tolerance &&
                     std::abs(top[0]) < tolerance &&
                     std::abs(top[1] - top_expected) < tolerance,
                 "pressure uses the parent Quad4 outward normal on left, "
                 "right, bottom, and top boundaries");
}

bool test_global_field_diagnostics(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SteadyProblem problem(single_region_definition(), mesh);
    const std::vector<double> state = problem.initial_state();
    std::vector<double> direction(problem.dof_count(), 0.0);
    for (std::size_t node = 0; node < problem.dof_map().node_count(); ++node) {
        direction[problem.dof_map().temperature(node)] = 0.2;
        direction[problem.dof_map().radial_displacement(node)] = 1.0e-6;
        direction[problem.dof_map().axial_displacement(node)] = -0.7e-6;
    }
    const fuelsim::DirectionalJacobianCheck diagnostic =
        fuelsim::check_directional_jacobian(problem, problem.dof_map(), state,
                                            direction, 1.0e-4);
    bool passed = true;
    for (std::size_t field = 0; field < 3; ++field) {
        const double reference =
            diagnostic.finite_difference_directional_derivative.l2[field];
        passed =
            check(diagnostic.difference.l2[field] <= 1.0e-7 * (1.0 + reference),
                  "global field Jacobian matches centered differences") &&
            passed;
    }
    return passed;
}

bool test_three_regions(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::SteadyProblem problem(three_region_definition(), mesh);
    const std::vector<double> state = problem.initial_state();
    bool passed =
        check(mesh.side_set_block_id("pellet_outer") == 1 &&
                  mesh.side_set_block_id("clad_1_inner") == 2 &&
                  mesh.side_set_block_id("clad_2_inner") == 3,
              "side-set ownership is inferred from adjacent Exodus elements") &&
        check(problem.region_count() == 3 && problem.contact_count() == 2,
              "three regions and two named contact pairs are composed") &&
        check(problem.region_node_offset(0) == 0 &&
                  problem.region_node_offset(1) == 4 &&
                  problem.region_node_offset(2) == 8,
              "arbitrary regions receive independent node ranges") &&
        check(problem.region_element_offset(0) == 0 &&
                  problem.region_element_offset(1) == 1 &&
                  problem.region_element_offset(2) == 2,
              "arbitrary regions receive independent element ranges") &&
        check(
            problem.dof_count() == 36 &&
                problem.volume_contribution_count() == 3 &&
                problem.contribution_count() == 9,
            "three volumes, two STS edges, and four NTS nodes are assembled") &&
        check(problem.contact(0).primary == "clad_1_inner" &&
                  problem.contact(0).secondary == "pellet_outer" &&
                  problem.contact(1).primary == "clad_2_inner" &&
                  problem.contact(1).secondary == "clad_1_outer",
              "contacts contain only primary and secondary side-set names");

    for (std::size_t contribution = problem.volume_contribution_count();
         contribution < problem.contribution_count(); ++contribution) {
        const fuelsim::LocalValues local =
            problem.contribution_state(contribution, state);
        const fuelsim::LocalSystem system =
            problem.linearize_contribution(contribution, local);
        passed = check(std::all_of(
                           system.residual.begin(), system.residual.end(),
                           [](double value) { return std::isfinite(value); }),
                       "multi-contact AD residual is finite") &&
                 passed;
    }
    fuelsim::SteadyProblemDefinition solve_definition =
        three_region_definition();
    for (fuelsim::RegionDefinition& region_value : solve_definition.regions) {
        region_value.initial_temperature = 300.0;
        region_value.volumetric_heat_source = 0.0;
    }
    for (fuelsim::ContactDefinition& contact_value :
         solve_definition.contacts)
        contact_value.penalty = 1.0e6;
    fuelsim::BoundaryConditionDefinition pellet_closure = dirichlet(
        "pellet_contact_closure", "pellet_outer",
        fuelsim::Field::radial_displacement, 0.101);
    pellet_closure.scale_with_load = true;
    solve_definition.boundary_conditions.push_back(pellet_closure);
    fuelsim::BoundaryConditionDefinition inner_clad_closure = dirichlet(
        "inner_clad_contact_closure", "clad_1_outer",
        fuelsim::Field::radial_displacement, 0.101);
    inner_clad_closure.scale_with_load = true;
    solve_definition.boundary_conditions.push_back(inner_clad_closure);
    solve_definition.boundary_conditions.push_back(dirichlet(
        "inner_clad_contact_anchor", "clad_1_inner",
        fuelsim::Field::radial_displacement, 0.0));
    solve_definition.boundary_conditions.push_back(dirichlet(
        "outer_clad_contact_anchor", "clad_2_inner",
        fuelsim::Field::radial_displacement, 0.0));
    fuelsim::SteadyProblem solve_problem(std::move(solve_definition), mesh);
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(
        solve_problem, {4, 0.5, 12, 1.0e-6}, fuelsim::SolverOptions{});
    if (!solve.completed)
        std::cerr << "multi-contact solve failure: "
                  << fuelsim::solve_failure_category_name(
                         solve.solve.failure_category)
                  << ": " << solve.solve.failure_message << '\n';
    passed =
        check(solve.completed && solve.solve.converged,
              "three-region two-contact PETSc solve converges") &&
        passed;
    if (solve.completed && solve.solve.converged) {
        for (std::size_t contact_value = 0;
             contact_value < solve_problem.contact_count(); ++contact_value) {
            const fuelsim::InterfaceSummary summary =
                solve_problem.summarize_interface(contact_value,
                                                  solve.solve.state);
            passed = check(summary.projected_contact_nodes > 0 &&
                               summary.active_contact_nodes > 0 &&
                               summary.unprojected_contact_nodes == 0,
                           "every converged contact pair remains projected "
                           "and mechanically active") &&
                     passed;
        }
    }
    return passed;
}

bool test_nonmatching_pellet_faces() {
    const fuelsim::UnstructuredQuad4Mesh mesh =
        two_pellet_nonmatching_mesh();
    fuelsim::SteadyProblemDefinition definition = {
        {region("lower", "lower_pellet", 700.0, 0.0),
         region("upper", "upper_pellet", 500.0, 0.0)},
        {contact("pellet_stack", "upper_bottom", "lower_top")},
        {}};
    const fuelsim::SteadyProblem problem(std::move(definition), mesh);
    const std::vector<double> state = problem.initial_state();
    bool passed =
        check(problem.region_count() == 2 && problem.contact_count() == 1,
              "two pellet blocks form one axial contact pair");
    for (std::size_t contribution = problem.volume_contribution_count();
         contribution < problem.contribution_count(); ++contribution) {
        const fuelsim::LocalResidual residual = problem.contribution_residual(
            contribution,
            problem.contribution_state(contribution, state));
        passed =
            check(std::all_of(residual.begin(), residual.end(),
                              [](double value) { return std::isfinite(value); }),
                  "nonmatching pellet-face contribution is finite") &&
            passed;
    }

    const fuelsim::InterfaceSummary summary =
        problem.summarize_interface(0, state);
    const double expected_area =
        3.141592653589793238462643383279502884 * 0.004 * 0.004;
    const double expected_heat_rate =
        expected_area * (0.2 / 1.0e-5) * (700.0 - 500.0);
    std::cout << "nonmatching_heat_rate=" << summary.total_heat_rate
              << " expected=" << expected_heat_rate << '\n';
    passed =
        check(std::abs(summary.total_heat_rate - expected_heat_rate) <
                  1.0e-12 * expected_heat_rate,
              "split nonmatching STS integration covers the full pellet face") &&
        check(summary.projected_contact_nodes == 3 &&
                  summary.active_contact_nodes == 0,
              "all nonmatching pellet-face NTS nodes project uniquely") &&
        passed;
    return passed;
}

bool test_transient_regions(const fuelsim::UnstructuredQuad4Mesh& mesh) {
    fuelsim::TransientProblemDefinition definition;
    definition.spatial = three_region_definition();
    for (const fuelsim::RegionDefinition& region : definition.spatial.regions)
        definition.regions.push_back(
            {region.name, elastic_transient_material()});

    fuelsim::TransientProblem problem(std::move(definition), mesh);
    const std::vector<double> initial = problem.committed_solution();
    problem.begin_time_step({1.0, 0.5});
    const fuelsim::LocalValues local = problem.contribution_state(2, initial);
    (void)problem.linearize_contribution(2, local);
    problem.rollback_time_step();

    bool passed =
        check(problem.region_count() == 3 && problem.committed_time() == 0.0 &&
                  problem.committed_solution() == initial &&
                  !problem.time_step_active(),
              "three-region transient rollback preserves committed state");

    problem.begin_time_step({1.0, 1.0});
    problem.commit_time_step(initial);
    passed = check(problem.committed_time() == 1.0 &&
                       problem.committed_load_factor() == 1.0 &&
                       !problem.time_step_active(),
                   "three-region transient state commits once") &&
             passed;
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const fuelsim::RegionInelasticSummary summary =
            problem.summarize_region_history(region);
        passed = check(summary.maximum_equivalent_plastic_strain == 0.0 &&
                           summary.maximum_equivalent_creep_strain == 0.0,
                       "elastic region keeps zero inelastic history") &&
                 passed;
    }
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(
            argc, argv, "fuelsim multi-region contact solve tests\n");
        const fuelsim::UnstructuredQuad4Mesh mesh = three_region_mesh();
        const bool passed =
            test_single_region(mesh) && test_time_controlled_pressure(mesh) &&
            test_pressure_parent_edge_orientation() &&
            test_global_field_diagnostics(mesh) && test_three_regions(mesh) &&
            test_nonmatching_pellet_faces() && test_transient_regions(mesh);
        if (!passed)
            return 1;
        std::cout << "[PASS] single- and multi-region problem tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] multi-region tests raised: " << error.what()
                  << '\n';
        return 1;
    }
}
