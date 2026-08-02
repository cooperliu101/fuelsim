#include "fuelsim/petsc_solver.hpp"
#include "support/steady_fuel_cladding_solver.hpp"

#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>

namespace {

constexpr std::size_t fuel_radial_elements = 100;
constexpr std::size_t cladding_radial_elements = 16;
constexpr std::size_t axial_elements = 64;
constexpr std::size_t load_steps = 20;

fuelsim::SteadyFuelCladdingParameters make_medium_case() {
    return {
        0.004120,
        0.004122,
        0.004692,
        0.010,
        0.010020,
        fuel_radial_elements,
        cladding_radial_elements,
        axial_elements,
        {
            3824.0,
            0.61,
            2.0e11,
            0.316,
            10.0e-6,
            600.0,
        },
        {
            0.0,
            16.0,
            75.0e9,
            0.3,
            5.0e-6,
            600.0,
        },
        2.0e8,
        600.0,
        600.0,
        0.4,
        1.0e-6,
        1.0e14,
    };
}

} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim M1 benchmark: 23,010 DOFs and 20 load "
            "steps\n");

        const fuelsim::SteadyFuelCladdingParameters parameters =
            make_medium_case();
        fuelsim::SolverOptions options;
        options.maximum_iterations = 50;
        const fuelsim::SteadyFuelCladdingLoadStepper load_stepper;
        const fuelsim::SteadyFuelCladdingLoadResult result =
            load_stepper.solve(parameters, load_steps, options);

        const std::size_t nodes =
            (fuel_radial_elements + 1 + cladding_radial_elements + 1) *
            (axial_elements + 1);
        const std::size_t elements =
            (fuel_radial_elements + cladding_radial_elements) * axial_elements;

        if (session.rank() != 0)
            return result.completed && result.solve.converged ? 0 : 1;

        std::cout << std::boolalpha << std::scientific << std::setprecision(12);
        std::cout << "case=m1-medium\n";
        std::cout << "mpi_ranks=" << session.size() << '\n';
        std::cout << "nodes=" << nodes << '\n';
        std::cout << "elements=" << elements << '\n';
        std::cout << "dofs=" << 3 * nodes << '\n';
        std::cout << "load_steps=" << load_steps << '\n';
        std::cout << "converged="
                  << (result.completed && result.solve.converged) << '\n';
        std::cout << "load_steps_completed=" << result.completed_steps << '\n';
        std::cout << "nonlinear_iterations_total="
                  << result.total_nonlinear_iterations << '\n';
        std::cout << "residual_norm_last_step=" << result.solve.residual_norm
                  << '\n';
        std::cout << "timing_problem_setup=" << result.problem_setup_seconds
                  << '\n';
        std::cout << "timing_solver_setup="
                  << result.aggregate_timing.setup_seconds << '\n';
        std::cout << "timing_nonlinear_solve="
                  << result.aggregate_timing.nonlinear_solve_seconds << '\n';
        std::cout << "timing_residual_callbacks="
                  << result.aggregate_timing.residual_callback_seconds << '\n';
        std::cout << "timing_jacobian_callbacks="
                  << result.aggregate_timing.jacobian_callback_seconds << '\n';
        std::cout << "timing_load_path_total=" << result.total_seconds << '\n';
        std::cout << "residual_evaluations="
                  << result.aggregate_timing.residual_evaluations << '\n';
        std::cout << "jacobian_evaluations="
                  << result.aggregate_timing.jacobian_evaluations << '\n';
        std::cout << "petsc_workspace_setups="
                  << result.aggregate_timing.workspace_setups << '\n';

        return result.completed && result.solve.converged ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
