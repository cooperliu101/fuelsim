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
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    bool passed = true;

    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem,
        {definition.steady_execution.load_steps,
         definition.steady_execution.cutback_factor,
         definition.steady_execution.maximum_cutbacks,
         definition.steady_execution.minimum_load_increment},
        options);
    passed = check(result.completed && result.solve.converged,
                   "non-tensor MOOSE mesh solve converged") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "non-tensor load path reuses one PETSc workspace") &&
             passed;

    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    passed = check(reference.size() == source.nodes().size(),
                   "MOOSE reference covers every Exodus node") &&
             passed;
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, result.solve.state,
                                                  reference);
    passed = check(fields.node_count == source.nodes().size() &&
                       fields.maximum_coordinate_difference < 1.0e-12,
                   "every MOOSE node is compared at matching coordinates") &&
             passed;

    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(
            pressure_reference_path, pressure_coordinates);
    const std::vector<fuelsim::ContactNodeSummary> pressure_values =
        problem.summarize_contact_nodes(0, result.solve.state);
    const fuelsim::test::FieldErrorMetrics pressure =
        fuelsim::test::compare_moose_contact_pressure(
            pressure_values, pressure_reference, pressure_coordinates, 1.0e-12);

    constexpr double tolerance = 1.0e-2;
    passed =
        check(fuelsim::test::relative_metrics_below(fields.temperature,
                                                    tolerance) &&
                  fuelsim::test::relative_metrics_below(
                      fields.radial_displacement, tolerance) &&
                  fuelsim::test::relative_metrics_below(
                      fields.axial_displacement, tolerance) &&
                  fuelsim::test::relative_metrics_below(pressure, tolerance),
              "all MOOSE full-field three-metric errors pass 1%") &&
        passed;

    fuelsim::test::print_relative_metrics("unstructured_temperature",
                                          fields.temperature);
    fuelsim::test::print_relative_metrics("unstructured_radial_displacement",
                                          fields.radial_displacement);
    fuelsim::test::print_relative_metrics("unstructured_axial_displacement",
                                          fields.axial_displacement);
    fuelsim::test::print_relative_metrics("unstructured_contact_pressure",
                                          pressure);
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m1_unstructured_moose_tests "
                     "<case.fsi> <all_nodes.csv> <fuel_surface.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim non-tensor Quad4 MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] non-tensor Quad4 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
