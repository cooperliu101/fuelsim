#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <cmath>
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
    if (definition.problem != fuelsim::CaseProblem::steady ||
        definition.contacts.size() != 1 ||
        definition.contacts[0].friction_coefficient != 0.3)
        throw std::invalid_argument(
            "M5.1 comparison requires the Coulomb friction input card");
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
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
    bool passed =
        check(result.completed && result.solve.converged,
              "M5.1 twenty-step frictional load path converges") &&
        check(result.aggregate_timing.workspace_setups == 1,
              "M5.1 frictional path reuses one PETSc workspace");

    const std::vector<double>& state = result.solve.state;
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, state, reference);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(0, state);
    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(
            pressure_reference_path, pressure_coordinates);
    const fuelsim::test::FieldErrorMetrics pressure =
        fuelsim::test::compare_moose_contact_pressure(
            contact_nodes, pressure_reference, pressure_coordinates, 1.0e-12);

    constexpr double tolerance = 1.0e-2;
    passed = check(fields.node_count == source.nodes().size() &&
                       fields.maximum_coordinate_difference < 1.0e-12,
                   "M5.1 compares every MOOSE node at matching coordinates") &&
             check(fuelsim::test::relative_metrics_below(fields.temperature,
                                                          tolerance),
                   "M5.1 temperature three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(
                       fields.radial_displacement, tolerance),
                   "M5.1 radial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(
                       fields.axial_displacement, tolerance),
                   "M5.1 axial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                   "M5.1 contact-pressure three full-field errors pass") &&
             passed;

    std::size_t active = 0;
    std::size_t sliding = 0;
    double maximum_capacity_excess = 0.0;
    for (const fuelsim::ContactNodeSummary& node : contact_nodes) {
        if (!(node.pressure > 0.0))
            continue;
        ++active;
        if (node.sliding)
            ++sliding;
        maximum_capacity_excess = std::max(
            maximum_capacity_excess,
            std::abs(node.tangential_traction) -
                definition.contacts[0].friction_coefficient * node.pressure);
    }
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(0, state);
    passed = check(active == contact_nodes.size() && sliding > 0 &&
                       maximum_capacity_excess <=
                           1.0e-12 * interface.maximum_contact_pressure &&
                       std::abs(interface.total_tangential_force) > 0.0,
                   "M5.1 activates friction, reaches sliding, respects every "
                   "Coulomb cap, and produces a nonzero shear resultant") &&
             passed;

    fuelsim::SpatialDefinition frictionless_definition =
        definition.spatial_definition();
    frictionless_definition.contacts[0].friction_coefficient = 0.0;
    fuelsim::SteadyProblem frictionless(std::move(frictionless_definition),
                                        source);
    const fuelsim::SteadyResult frictionless_result = fuelsim::solve_steady(
        frictionless,
        {definition.steady_execution.load_steps,
         definition.steady_execution.cutback_factor,
         definition.steady_execution.maximum_cutbacks,
         definition.steady_execution.minimum_load_increment},
        options);
    double axial_difference_squared = 0.0;
    double axial_scale_squared = 0.0;
    for (std::size_t node = 0; node < problem.dof_map().node_count(); ++node) {
        const std::size_t dof = problem.dof_map().axial_displacement(node);
        const double difference =
            state[dof] - frictionless_result.solve.state[dof];
        axial_difference_squared += difference * difference;
        axial_scale_squared += state[dof] * state[dof];
    }
    const double friction_effect =
        std::sqrt(axial_difference_squared / axial_scale_squared);
    passed = check(frictionless_result.completed && friction_effect > 1.0e-3,
                   "M5.1 MOOSE comparison exercises a measurable frictional "
                   "change rather than the mu-zero baseline") &&
             passed;

    fuelsim::test::print_relative_metrics("m51_temperature",
                                          fields.temperature);
    fuelsim::test::print_relative_metrics("m51_radial_displacement",
                                          fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m51_axial_displacement",
                                          fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m51_contact_pressure", pressure);
    std::cout << "m51_active_contact_nodes=" << active << '\n'
              << "m51_sliding_contact_nodes=" << sliding << '\n'
              << "m51_total_tangential_force="
              << interface.total_tangential_force << '\n'
              << "m51_frictional_axial_field_change=" << friction_effect
              << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m51_coulomb_friction_moose_tests "
                     "<case.fsi> <all-nodes.csv> <contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M5.1 Coulomb friction MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] M5.1 Coulomb friction MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.1 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
