#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
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

bool run_comparison(const std::string& input_path,
                    const std::string& nodal_reference_path,
                    const std::string& pressure_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument(
            "M1 comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    bool passed =
        check(source.nodes().size() == 528 && source.elements().size() == 460,
              "M1 input card reads all MOOSE Exodus entities");
    passed = check(source.element_block("fuel").id == 0 &&
                       source.element_block("clad").id == 1 &&
                       source.side_set("fuel_right").sides.size() == 10 &&
                       source.side_set("clad_left").sides.size() == 10,
                   "M1 MOOSE block and side-set metadata are preserved") &&
             passed;

    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    const fuelsim::SolverOptions options = {
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance, definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem,
        {definition.steady_execution.load_steps,
         definition.steady_execution.cutback_factor,
         definition.steady_execution.maximum_cutbacks,
         definition.steady_execution.minimum_load_increment},
        options);
    passed = check(result.completed && result.solve.converged,
                   "M1 input-card load path converged") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "M1 input-card path reuses one PETSc workspace") &&
             passed;

    const std::vector<double>& state = result.solve.state;
    const std::vector<fuelsim::test::NodalFieldReference> nodal_reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, state,
                                                  nodal_reference);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(0, state);
    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(
            pressure_reference_path, pressure_coordinates);
    const fuelsim::test::FieldErrorMetrics pressure =
        fuelsim::test::compare_moose_contact_pressure(
            contact_nodes, pressure_reference, pressure_coordinates, 1.0e-12);

    passed = check(fields.node_count == source.nodes().size() &&
                       fields.maximum_coordinate_difference < 1.0e-12,
                   "M1 compares every MOOSE node at matching coordinates") &&
             passed;

    constexpr double tolerance = 1.0e-2;
    passed = check(fuelsim::test::relative_metrics_below(fields.temperature,
                                                         tolerance),
                   "M1 full-field temperature three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       fields.radial_displacement, tolerance),
                   "M1 full-field radial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       fields.axial_displacement, tolerance),
                   "M1 full-field axial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                   "M1 full-field contact pressure three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m1_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m1_radial_displacement",
                                          fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m1_axial_displacement",
                                          fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m1_contact_pressure", pressure);

    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(0, state);
    fuelsim::test::FieldErrorMetrics total_force;
    total_force.add(interface.total_contact_force, 663.8896691615588);
    passed =
        check(fuelsim::test::relative_metrics_below(total_force, tolerance),
              "M1 total contact force three errors pass") &&
        passed;
    fuelsim::test::print_relative_metrics("m1_total_contact_force",
                                          total_force);
    passed = check(interface.projected_contact_nodes == 11 &&
                       interface.active_contact_nodes == 11,
                   "M1 projects and activates all fuel-surface nodes") &&
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
        fuelsim::PetscSession session(
            argc, argv, "fuelsim input-card M1 MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] input-card M1 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M1 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
