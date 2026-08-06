#include "support/steady_fuel_cladding_solver.hpp"

#include <chrono>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {

using SteadyClock = std::chrono::steady_clock;

double seconds_since(const SteadyClock::time_point& start) {
    return std::chrono::duration<double>(SteadyClock::now() - start).count();
}

void accumulate_timing(SolveTiming& total, const SolveTiming& step) {
    total.setup_seconds += step.setup_seconds;
    total.nonlinear_solve_seconds += step.nonlinear_solve_seconds;
    total.residual_callback_seconds += step.residual_callback_seconds;
    total.jacobian_callback_seconds += step.jacobian_callback_seconds;
    total.total_seconds += step.total_seconds;
    total.residual_evaluations += step.residual_evaluations;
    total.jacobian_evaluations += step.jacobian_evaluations;
    total.workspace_setups += step.workspace_setups;
    total.solve_calls += step.solve_calls;
}

} // namespace

SteadyFuelCladdingLoadResult SteadyFuelCladdingLoadStepper::solve(
    const SteadyFuelCladdingParameters& target_parameters,
    std::size_t load_steps, const SolverOptions& options) const {
    return solve(
        target_parameters,
        StructuredRzMesh::make_annulus(0.0, target_parameters.fuel_radius,
                                       target_parameters.fuel_length,
                                       target_parameters.fuel_radial_elements,
                                       target_parameters.axial_elements),
        StructuredRzMesh::make_annulus(
            target_parameters.cladding_inner_radius,
            target_parameters.cladding_outer_radius,
            target_parameters.cladding_length,
            target_parameters.cladding_radial_elements,
            target_parameters.axial_elements),
        load_steps, options);
}

SteadyFuelCladdingLoadResult SteadyFuelCladdingLoadStepper::solve(
    const SteadyFuelCladdingParameters& target_parameters,
    StructuredRzMesh fuel_mesh, StructuredRzMesh cladding_mesh,
    std::size_t load_steps, const SolverOptions& options) const {
    if (load_steps == 0)
        throw std::invalid_argument(
            "SteadyFuelCladdingLoadStepper load_steps must be positive");

    const SteadyClock::time_point total_start = SteadyClock::now();
    const SteadyClock::time_point problem_setup_start = SteadyClock::now();
    SteadyFuelCladdingProblem problem(target_parameters, std::move(fuel_mesh),
                                      std::move(cladding_mesh));

    SteadyFuelCladdingLoadResult result;
    result.problem_setup_seconds = seconds_since(problem_setup_start);

    PetscSolver nonlinear_solver;
    std::vector<double> state = problem.initial_state();

    for (std::size_t step = 1; step <= load_steps; ++step) {
        const double step_heat_source =
            target_parameters.volumetric_heat_source *
            static_cast<double>(step) / static_cast<double>(load_steps);
        problem.set_volumetric_heat_source(step_heat_source);

        SolveResult step_result =
            nonlinear_solver.solve(problem, state, options);
        accumulate_timing(result.aggregate_timing, step_result.timing);
        result.total_nonlinear_iterations += step_result.nonlinear_iterations;
        result.total_linear_iterations += step_result.linear_iterations;
        if (!step_result.converged) {
            result.solve = std::move(step_result);
            result.completed_steps = step - 1;
            result.total_seconds = seconds_since(total_start);
            return result;
        }

        result.completed_steps = step;
        if (step == load_steps)
            result.solve = std::move(step_result);
        else
            state = std::move(step_result.state);
    }

    result.completed = true;
    result.total_seconds = seconds_since(total_start);
    return result;
}

} // namespace fuelsim
