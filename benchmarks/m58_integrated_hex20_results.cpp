#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
class ProgressObserver final : public fuelsim::TransientStepObserver {
  public:
    explicit ProgressObserver(bool active) : _active(active) {}

    void accepted_step(const fuelsim::TransientProblem&, const fuelsim::TransientAcceptedStep& step) override {
        if (!_active) return;
        ++_accepted_steps;
        std::cout << "m58_hex20_progress_accepted_steps=" << _accepted_steps << '\n'
                  << "m58_hex20_progress_time=" << step.time << '\n'
                  << "m58_hex20_progress_time_step=" << step.time_step << '\n'
                  << "m58_hex20_progress_nonlinear_iterations=" << step.nonlinear_iterations << '\n'
                  << "m58_hex20_progress_linear_iterations=" << step.linear_iterations << '\n'
                  << "m58_hex20_progress_cutbacks=" << step.cutbacks << '\n'
                  << std::flush;
    }

  private:
    bool _active;
    std::size_t _accepted_steps = 0;
};

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = definition.solver.absolute_tolerance;
    result.relative_tolerance = definition.solver.relative_tolerance;
    result.step_tolerance = definition.solver.step_tolerance;
    result.maximum_iterations = definition.solver.maximum_iterations;
    result.linear_solver = definition.solver.linear_solver;
    result.preconditioner = definition.solver.preconditioner;
    result.direct_factorization = definition.solver.direct_factorization;
    result.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    result.jacobian_lag = definition.solver.jacobian_lag;
    result.line_search = definition.solver.line_search;
    result.backtracking_fallback = definition.solver.backtracking_fallback;
    result.field_residual_scaling = definition.solver.field_residual_scaling;
    result.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
    result.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
    result.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
    result.temperature_residual_scale = definition.solver.temperature_residual_scale;
    result.mechanical_residual_scale = definition.solver.mechanical_residual_scale;
    return result;
}

fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    const auto& input = definition.transient_execution;
    return {input.end_time, input.initial_time_step, input.minimum_time_step, input.maximum_time_step,
        input.growth_factor, input.cutback_factor, input.maximum_cutbacks_per_step, input.load_ramp_time,
        input.target_nonlinear_iterations, input.iteration_window, input.time_error_relative_tolerance,
        input.temperature_time_absolute_tolerance, input.displacement_time_absolute_tolerance,
        input.time_error_safety_factor, input.strain_history_time_absolute_tolerance,
        input.stress_history_time_absolute_tolerance, input.include_thermal_time_term, input.use_linear_time_predictor};
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_m58_integrated_hex20_results <case.fsi> <final-results.e>\n";
        return 2;
    }
    try {
        fuelsim::PetscSession session(argc, argv, "fuelsim M5.8 integrated HEX20 final-result writer\n");
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(argv[1]);
        if (definition.problem != fuelsim::CaseProblem::transient ||
            definition.geometry != fuelsim::CaseGeometry::cartesian_3d)
            throw std::invalid_argument("M5.8 HEX20 result writing requires a three-dimensional transient case");
        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
        fuelsim::TransientProblem problem(definition.spatial, mesh);
        std::vector<std::size_t> local_dofs;
        std::size_t maximum_contribution_dofs = 0, maximum_sparsity_dofs = 0;
        for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
            problem.contribution_dofs(contribution, local_dofs);
            maximum_contribution_dofs = std::max(maximum_contribution_dofs, local_dofs.size());
        }
        for (std::size_t contribution = 0; contribution < problem.sparsity_contribution_count(); ++contribution) {
            problem.sparsity_contribution_dofs(contribution, local_dofs);
            maximum_sparsity_dofs = std::max(maximum_sparsity_dofs, local_dofs.size());
        }
        if (session.rank() == 0)
            std::cout << "m58_hex20_dofs=" << problem.dof_count() << '\n'
                      << "m58_hex20_contributions=" << problem.contribution_count() << '\n'
                      << "m58_hex20_sparsity_contributions=" << problem.sparsity_contribution_count() << '\n'
                      << "m58_hex20_maximum_contribution_dofs=" << maximum_contribution_dofs << '\n'
                      << "m58_hex20_maximum_sparsity_dofs=" << maximum_sparsity_dofs << '\n'
                      << std::flush;
        ProgressObserver observer(session.rank() == 0);
        const fuelsim::TransientResult result =
            fuelsim::solve_transient(problem, time_options(definition), solver_options(definition), &observer);
        if (session.rank() == 0)
            std::cout << "m58_hex20_accepted_steps=" << result.accepted_steps.size() << '\n'
                      << "m58_hex20_rejected_steps=" << result.rejected_steps.size() << '\n'
                      << "m58_hex20_nonlinear_iterations=" << result.total_nonlinear_iterations << '\n'
                      << "m58_hex20_linear_iterations=" << result.total_linear_iterations << '\n'
                      << "m58_hex20_residual_evaluations=" << result.aggregate_timing.residual_evaluations << '\n'
                      << "m58_hex20_jacobian_evaluations=" << result.aggregate_timing.jacobian_evaluations << '\n'
                      << "m58_hex20_workspace_setups=" << result.aggregate_timing.workspace_setups << '\n'
                      << "m58_hex20_setup_seconds=" << result.aggregate_timing.setup_seconds << '\n'
                      << "m58_hex20_residual_callback_seconds=" << result.aggregate_timing.residual_callback_seconds
                      << '\n'
                      << "m58_hex20_jacobian_callback_seconds=" << result.aggregate_timing.jacobian_callback_seconds
                      << '\n'
                      << "m58_hex20_nonlinear_solve_seconds=" << result.aggregate_timing.nonlinear_solve_seconds << '\n'
                      << "m58_hex20_total_seconds=" << result.aggregate_timing.total_seconds << '\n'
                      << "m58_hex20_dofs=" << problem.dof_count() << '\n';
        if (!result.completed) {
            std::cerr << "M5.8 HEX20 solve failed: " << result.last_attempt.failure_message << '\n';
            return 1;
        }
        fuelsim::ExodusTransientResultsWriter writer(argv[2], mesh, problem);
        writer.append(problem);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "M5.8 HEX20 result writing failed: " << error.what() << '\n';
        return 1;
    }
}
