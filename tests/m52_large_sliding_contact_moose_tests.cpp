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
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    const std::vector<fuelsim::ContactNodeSummary> initial_contact_nodes =
        problem.summarize_contact_nodes(0, problem.initial_state());
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
              "M5.2 twenty-step large-sliding path converges") &&
        check(result.aggregate_timing.workspace_setups == 1,
              "M5.2 keeps one PETSc workspace and locked sparsity");

    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, result.solve.state,
                                                  reference);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(0, result.solve.state);
    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(
            pressure_reference_path, pressure_coordinates);
    const fuelsim::test::FieldErrorMetrics pressure =
        fuelsim::test::compare_moose_contact_pressure(
            contact_nodes, pressure_reference, pressure_coordinates, 1.0e-12);

    constexpr double tolerance = 1.0e-2;
    passed = check(fields.node_count == 402 &&
                       fields.maximum_coordinate_difference < 1.0e-12,
                   "M5.2 compares all 402 MOOSE nodes on the tracked mesh") &&
             check(fuelsim::test::relative_metrics_below(fields.temperature,
                                                          tolerance),
                   "M5.2 temperature three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(
                       fields.radial_displacement, tolerance),
                   "M5.2 radial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(
                       fields.axial_displacement, tolerance),
                   "M5.2 axial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                   "M5.2 contact-pressure three full-field errors pass") &&
             passed;

    std::size_t active = 0;
    std::size_t maximum_segment_change = 0;
    double maximum_current_radius = 0.0;
    const std::vector<std::size_t> secondary_sources =
        problem.contact_secondary_source_nodes(0);
    for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
        const fuelsim::ContactNodeSummary& contact = contact_nodes[node];
        if (contact.pressure > 0.0)
            ++active;
        const std::size_t initial_segment =
            initial_contact_nodes.at(node).primary_segment;
        if (contact.primary_segment >= initial_segment)
            maximum_segment_change =
                std::max(maximum_segment_change,
                         contact.primary_segment - initial_segment);
        const std::size_t source_node = secondary_sources[node];
        std::size_t global_node = problem.dof_map().node_count();
        for (std::size_t region = 0; region < problem.region_count(); ++region) {
            const std::vector<std::size_t>& source_nodes =
                problem.region_mesh(region).source_node_ids();
            const auto found =
                std::find(source_nodes.begin(), source_nodes.end(), source_node);
            if (found == source_nodes.end())
                continue;
            global_node = problem.region_node_offset(region) +
                          static_cast<std::size_t>(found - source_nodes.begin());
            break;
        }
        if (global_node == problem.dof_map().node_count())
            throw std::logic_error("M5.2 secondary node mapping failed");
        maximum_current_radius = std::max(
            maximum_current_radius,
            contact.r + result.solve.state[problem.dof_map().radial_displacement(
                            global_node)]);
    }
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(0, result.solve.state);
    passed = check(active == contact_nodes.size() &&
                       maximum_segment_change >= 2 &&
                       maximum_current_radius > 4.2e-3 &&
                       maximum_current_radius < 8.0e-3 &&
                       interface.total_contact_force > 0.0,
                   "M5.2 keeps all five nodes active while crossing at least "
                   "two primary segments inside the full chain") &&
             passed;

    fuelsim::test::print_relative_metrics("m52_temperature",
                                          fields.temperature);
    fuelsim::test::print_relative_metrics("m52_radial_displacement",
                                          fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m52_axial_displacement",
                                          fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m52_contact_pressure", pressure);
    std::cout << "m52_active_contact_nodes=" << active << '\n'
              << "m52_maximum_primary_segment_change="
              << maximum_segment_change << '\n'
              << "m52_maximum_secondary_current_radius="
              << maximum_current_radius << '\n'
              << "m52_total_contact_force="
              << interface.total_contact_force << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m52_large_sliding_contact_moose_tests "
                     "<case.fsi> <all-nodes.csv> <contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M5.2 large-sliding MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] M5.2 large-sliding MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.2 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
