#include "fuelsim/m1_solver.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace fuelsim {

M1LoadStepResult M1LoadStepper::solve(const M1Parameters& target_parameters,
                                      std::size_t load_steps,
                                      const SolverOptions& options) const {
    if (load_steps == 0)
        throw std::invalid_argument(
            "M1LoadStepper load_steps must be positive");

    const PetscSequentialSolver nonlinear_solver;
    std::vector<double> state;
    SolveResult step_result{};

    for (std::size_t step = 1; step <= load_steps; ++step) {
        M1Parameters step_parameters = target_parameters;
        step_parameters.volumetric_heat_source =
            target_parameters.volumetric_heat_source *
            static_cast<double>(step) / static_cast<double>(load_steps);
        const M1Problem step_problem(step_parameters);
        if (state.empty())
            state = step_problem.initial_state();

        step_result = nonlinear_solver.solve(step_problem, state, options);
        if (!step_result.converged)
            return {std::move(step_result), step - 1, false};
        state = step_result.state;
    }

    return {std::move(step_result), load_steps, true};
}

} // namespace fuelsim
