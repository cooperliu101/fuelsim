#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool test_moose_mesh_backward_euler_heat_source(
    const std::string& input_path, const std::string& nodal_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient)
        throw std::invalid_argument(
            "M2.1 comparison requires a transient input card");
    const fuelsim::UnstructuredQuad4Mesh imported =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.transient_definition(),
                                      imported);
    const fuelsim::TransientTimeOptions time_options = {
        definition.transient_execution.end_time,
        definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step,
        definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor,
        definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks,
        definition.transient_execution.load_ramp_time};
    const fuelsim::SolverOptions solver_options = {
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance, definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, time_options, solver_options);

    bool passed = check(problem.dof_map().node_count() == 15 &&
                            problem.contribution_count() == 8,
                        "M2.1 uses all 8 elements from the MOOSE Exodus mesh");
    passed = check(result.completed && result.accepted_steps.size() == 10,
                   "M2.1 input-card transient completes ten steps") &&
             passed;

    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(
            problem, result.committed_state, reference);
    passed = check(fields.node_count == imported.nodes().size() &&
                       fields.maximum_coordinate_difference < 1.0e-12,
                   "M2.1 compares every MOOSE node at matching coordinates") &&
             passed;
    passed =
        check(fuelsim::test::relative_metrics_below(fields.temperature, 1.0e-3),
              "M2.1 full-field temperature three errors pass") &&
        passed;
    passed = check(fuelsim::test::absolute_metrics_below(
                       fields.radial_displacement, 1.0e-12),
                   "M2.1 zero radial field three absolute errors pass") &&
             passed;
    passed = check(fuelsim::test::absolute_metrics_below(
                       fields.axial_displacement, 1.0e-12),
                   "M2.1 zero axial field three absolute errors pass") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "ten backward-Euler steps create one PETSc workspace") &&
             passed;
    passed = check(result.aggregate_timing.solve_calls == 10,
                   "one PETSc solve is issued per backward-Euler step") &&
             passed;

    fuelsim::test::print_relative_metrics("m21_temperature",
                                          fields.temperature);
    fuelsim::test::print_absolute_metrics("m21_radial_displacement",
                                          fields.radial_displacement);
    fuelsim::test::print_absolute_metrics("m21_axial_displacement",
                                          fields.axial_displacement);
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_m2_solver_tests <m21.fsi> "
                     "<all-nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim input-card M2.1 MOOSE comparison test\n");
        if (!test_moose_mesh_backward_euler_heat_source(argv[1], argv[2]))
            return 1;
        std::cout << "[PASS] input-card M2.1 MOOSE comparison test\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M2.1 test raised: " << error.what() << '\n';
        return 1;
    }
}
