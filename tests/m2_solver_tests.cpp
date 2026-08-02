#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool test_moose_mesh_backward_euler_heat_source(const std::string& input_path) {
    constexpr double expected_temperature = 610.0;
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

    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_actual = 0.0;
    double maximum_pointwise_relative = 0.0;
    for (std::size_t node = 0; node < problem.dof_map().node_count(); ++node) {
        const double actual =
            result.committed_state.at(problem.dof_map().temperature(node));
        const double difference = actual - expected_temperature;
        difference_squared += difference * difference;
        reference_squared += expected_temperature * expected_temperature;
        maximum_actual = std::max(maximum_actual, std::abs(actual));
        maximum_pointwise_relative =
            std::max(maximum_pointwise_relative,
                     std::abs(difference) / expected_temperature);
    }
    const double relative_l2 =
        std::sqrt(difference_squared / reference_squared);
    const double relative_absolute_peak =
        std::abs(maximum_actual - expected_temperature) / expected_temperature;
    passed = check(relative_l2 < 1.0e-3 && relative_absolute_peak < 1.0e-3 &&
                       maximum_pointwise_relative < 1.0e-3,
                   "M2.1 three MOOSE temperature metrics are below 0.1%") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "ten backward-Euler steps create one PETSc workspace") &&
             passed;
    passed = check(result.aggregate_timing.solve_calls == 10,
                   "one PETSc solve is issued per backward-Euler step") &&
             passed;

    std::cout << "m21_temperature_relative_l2=" << relative_l2 << '\n';
    std::cout << "m21_temperature_relative_absolute_peak="
              << relative_absolute_peak << '\n';
    std::cout << "m21_temperature_maximum_pointwise_relative="
              << maximum_pointwise_relative << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m2_solver_tests <m21.fsi>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim input-card M2.1 MOOSE comparison test\n");
        if (!test_moose_mesh_backward_euler_heat_source(argv[1]))
            return 1;
        std::cout << "[PASS] input-card M2.1 MOOSE comparison test\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M2.1 test raised: " << error.what() << '\n';
        return 1;
    }
}
