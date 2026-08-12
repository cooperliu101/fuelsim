#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/rz_problem_access.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
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

bool run_comparison(const std::string& input_path, const std::string& nodal_reference_path,
                    const std::string& pressure_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument("M5.4 comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
    const fuelsim::SolverOptions solver = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
                                           definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
         definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        solver);
    bool passed =
        check(result.completed && result.solve.converged, "M5.4 augmented-contact load path converges") &&
        check(result.solve.augmented_lagrangian_iterations > 0, "M5.4 final load step performs a multiplier update") &&
        check(result.aggregate_timing.workspace_setups == 1,
              "M5.4 reuses one PETSc workspace across all outer iterations");

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
    constexpr double field_tolerance = 1.0e-2;
    passed = check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                   "M5.4 compares every MOOSE node at matching coordinates") &&
             check(fuelsim::test::relative_metrics_below(fields.temperature, field_tolerance),
                   "M5.4 temperature three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(fields.radial_displacement, field_tolerance),
                   "M5.4 radial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(fields.axial_displacement, field_tolerance),
                   "M5.4 axial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(pressure, field_tolerance),
                   "M5.4 contact-pressure three full-field errors pass") &&
             passed;

    const fuelsim::InterfaceSummary interface = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
    const double maximum_penetration = std::max(-interface.minimum_contact_gap, 0.0);
    passed = check(interface.projected_contact_nodes == contact_nodes.size() &&
                       interface.active_contact_nodes == contact_nodes.size() &&
                       maximum_penetration <= definition.contacts.at(0).penetration_tolerance,
                   "M5.4 activates every contact node and satisfies the 1 nm "
                   "penetration tolerance") &&
             passed;

    fuelsim::test::print_relative_metrics("m54_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m54_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m54_axial_displacement", fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m54_contact_pressure", pressure);
    std::cout << "m54_maximum_penetration=" << maximum_penetration << '\n'
              << "m54_penetration_tolerance=" << definition.contacts.at(0).penetration_tolerance << '\n'
              << "m54_multiplier_updates=" << result.solve.augmented_lagrangian_iterations << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m54_augmented_contact_moose_tests "
                     "<case.fsi> <all-nodes.csv> <contact-pressure.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M5.4 augmented-contact MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] M5.4 augmented-contact MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.4 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
