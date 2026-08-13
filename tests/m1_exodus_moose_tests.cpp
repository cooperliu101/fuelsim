#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
bool run_comparison(const std::string& input_path, const std::string& nodal_reference_path,
    const std::string& pressure_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument("M1 comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    bool passed = check(source.nodes().size() == 528 && source.elements().size() == 460,
        "M1 input card reads all MOOSE Exodus entities");
    passed =
        check(source.element_block("fuel").id == 0 && source.element_block("clad").id == 1 &&
                  source.side_set("fuel_right").sides.size() == 10 && source.side_set("clad_left").sides.size() == 10,
            "M1 MOOSE block and side-set metadata are preserved") &&
        passed;
    fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
    const fuelsim::SolverOptions options = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        options);
    passed = check(result.completed && result.solve.converged, "M1 input-card load path converged") && passed;
    passed =
        check(result.aggregate_timing.workspace_setups == 1, "M1 input-card path reuses one PETSc workspace") && passed;
    const std::vector<double>& state = result.solve.state;
    const std::vector<fuelsim::test::NodalFieldReference> nodal_reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, state, nodal_reference);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, state);
    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(pressure_reference_path, pressure_coordinates);
    const fuelsim::test::FieldErrorMetrics pressure =
        fuelsim::test::compare_moose_contact_pressure(contact_nodes, pressure_reference, pressure_coordinates, 1.0e-12);
    passed = check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                 "M1 compares every MOOSE node at matching coordinates") &&
             passed;
    constexpr double tolerance = 1.0e-2;
    passed = check(fuelsim::test::relative_metrics_below(fields.temperature, tolerance),
                 "M1 full-field temperature three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.radial_displacement, tolerance),
                 "M1 full-field radial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.axial_displacement, tolerance),
                 "M1 full-field axial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                 "M1 full-field contact pressure three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m1_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m1_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m1_axial_displacement", fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m1_contact_pressure", pressure);
    const fuelsim::InterfaceSummary interface = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
    fuelsim::test::FieldErrorMetrics total_force;
    total_force.add(interface.total_contact_force, 663.8896691615588);
    passed = check(fuelsim::test::relative_metrics_below(total_force, tolerance),
                 "M1 total contact force three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m1_total_contact_force", total_force);
    passed = check(interface.projected_contact_nodes == 11 && interface.active_contact_nodes == 11,
                 "M1 projects and activates all fuel-surface nodes") &&
             passed;
    const auto solve_penalty = [&](double penalty) {
        fuelsim::SpatialDefinition modified = definition.spatial_definition();
        modified.contacts.at(0).penalty = penalty;
        fuelsim::SteadyProblem penalty_problem(std::move(modified), source);
        const fuelsim::SteadyResult penalty_result = fuelsim::solve_steady(penalty_problem,
            {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
                definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
            options);
        if (!penalty_result.completed || !penalty_result.solve.converged)
            throw std::runtime_error("M4 penalty-convergence solve did not converge");
        return fuelsim::rz::ProblemAccess::summarize_interface(penalty_problem, 0, penalty_result.solve.state);
    };
    const fuelsim::InterfaceSummary low_penalty = solve_penalty(2.5e13);
    const fuelsim::InterfaceSummary medium_penalty = solve_penalty(5.0e13);
    const double low_penetration = -low_penalty.minimum_contact_gap;
    const double medium_penetration = -medium_penalty.minimum_contact_gap;
    const double high_penetration = -interface.minimum_contact_gap;
    const double low_to_medium_force_change =
        std::abs(medium_penalty.total_contact_force - low_penalty.total_contact_force);
    const double medium_to_high_force_change =
        std::abs(interface.total_contact_force - medium_penalty.total_contact_force);
    std::cout << "penalty_convergence_penetrations=" << low_penetration << ',' << medium_penetration << ','
              << high_penetration << '\n';
    std::cout << "penalty_convergence_force_changes=" << low_to_medium_force_change << ','
              << medium_to_high_force_change << '\n';
    passed = check(low_penetration > medium_penetration && medium_penetration > high_penetration &&
                       high_penetration > 0.0 && medium_to_high_force_change < low_to_medium_force_change,
                 "penalty refinement reduces penetration and contact-force "
                 "increments") &&
             passed;
    fuelsim::SpatialDefinition automatic_definition = definition.spatial_definition();
    automatic_definition.contacts[0].automatic_penalty = true;
    automatic_definition.contacts[0].penalty = 0.0;
    automatic_definition.contacts[0].penalty_factor = 1.0;
    fuelsim::SteadyProblem automatic_problem(std::move(automatic_definition), source);
    const double fuel_normal_length = 0.00412 / 40.0;
    const double clad_normal_length = (0.004692 - 0.004122) / 6.0;
    const double interface_stiffness = 1.0 / (fuel_normal_length / 2.0e11 + clad_normal_length / 7.5e10);
    const double expected_automatic_penalty = interface_stiffness;
    const fuelsim::SteadyResult automatic_result = fuelsim::solve_steady(automatic_problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        options);
    const double automatic_penalty = fuelsim::rz::ProblemAccess::contact(automatic_problem, 0).penalty;
    std::cout << "automatic_penalty=" << automatic_penalty << '\n';
    std::cout << "automatic_penalty_interface_stiffness=" << interface_stiffness << '\n';
    passed = check(automatic_result.completed && automatic_result.solve.converged &&
                       std::abs(automatic_penalty - expected_automatic_penalty) < 1.0e-12 * expected_automatic_penalty,
                 "automatic penalty uses the two-sided normal compliance and "
                 "converges end to end") &&
             passed;
    fuelsim::SpatialDefinition augmented_definition = definition.spatial_definition();
    augmented_definition.contacts[0].mechanical_formulation =
        fuelsim::MechanicalContactFormulation::augmented_lagrangian;
    augmented_definition.contacts[0].automatic_penalty = false;
    augmented_definition.contacts[0].penalty = 0.25 * interface_stiffness;
    augmented_definition.contacts[0].penetration_tolerance = 1.0e-9;
    augmented_definition.contacts[0].maximum_augmented_iterations = 50;
    fuelsim::SteadyProblem augmented_problem(std::move(augmented_definition), source);
    const fuelsim::SteadyResult augmented_result = fuelsim::solve_steady(augmented_problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        options);
    const fuelsim::InterfaceSummary augmented_interface =
        fuelsim::rz::ProblemAccess::summarize_interface(augmented_problem, 0, augmented_result.solve.state);
    const double augmented_penetration = std::max(-augmented_interface.minimum_contact_gap, 0.0);
    const double high_penalty_condition_proxy = 1.0 + 2.0 * (10.0 * interface_stiffness) / interface_stiffness;
    const double augmented_condition_proxy =
        1.0 + 2.0 * fuelsim::rz::ProblemAccess::contact(augmented_problem, 0).penalty / interface_stiffness;
    std::cout << "augmented_penetration=" << augmented_penetration << '\n';
    std::cout << "augmented_multiplier_updates=" << augmented_result.solve.augmented_lagrangian_iterations << '\n';
    std::cout << "contact_condition_proxies=" << high_penalty_condition_proxy << ',' << augmented_condition_proxy
              << '\n';
    passed =
        check(augmented_result.completed && augmented_result.solve.converged &&
                  augmented_result.solve.augmented_lagrangian_iterations > 0 && augmented_penetration <= 1.0e-9 &&
                  augmented_penetration < high_penetration && augmented_result.aggregate_timing.workspace_setups == 1 &&
                  augmented_condition_proxy < 0.1 * high_penalty_condition_proxy,
            "augmented contact reaches the penetration tolerance with "
            "one PETSc workspace and a lower two-body tangent condition "
            "proxy than the high automatic penalty") &&
        passed;
    fuelsim::SpatialDefinition failing_definition = definition.spatial_definition();
    failing_definition.contacts[0].mechanical_formulation = fuelsim::MechanicalContactFormulation::augmented_lagrangian;
    failing_definition.contacts[0].automatic_penalty = false;
    failing_definition.contacts[0].penalty = 0.25 * interface_stiffness;
    failing_definition.contacts[0].penetration_tolerance = 1.0e-20;
    failing_definition.contacts[0].maximum_augmented_iterations = 1;
    fuelsim::SteadyProblem failing_problem(std::move(failing_definition), source);
    const fuelsim::SteadyResult failing_result = fuelsim::solve_steady(failing_problem, {1, 0.5, 0, 1.0e-6}, options);
    bool multiplier_rolled_back = true;
    for (const fuelsim::ContactPointHistory& history :
        fuelsim::rz::ProblemAccess::committed_contact_histories(failing_problem).at(0))
        multiplier_rolled_back = multiplier_rolled_back && history.normal_multiplier == 0.0;
    passed = check(!failing_result.completed && !failing_result.solve.converged &&
                       failing_result.solve.failure_category == fuelsim::SolveFailureCategory::contact_constraint &&
                       multiplier_rolled_back && failing_problem.load_factor() == 0.0,
                 "failed augmented load step restores the accepted load and all "
                 "normal multipliers") &&
             passed;
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m1_exodus_moose_tests <m1.fsi> "
                     "<all-nodes.csv> <contact-pressure.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim input-card M1 MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3])) return 1;
        std::cout << "[PASS] input-card M1 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M1 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
