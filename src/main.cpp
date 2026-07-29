#include "fuelsim/m1_problem.hpp"
#include "fuelsim/m1_solver.hpp"
#include "fuelsim/petsc_solver.hpp"

#include <exception>
#include <iomanip>
#include <iostream>

namespace {

fuelsim::M1Parameters make_reference_case() {
    return {
        0.004120,
        0.004122,
        0.004692,
        0.010,
        40,
        6,
        10,
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
            "fuelsim M1: axisymmetric fuel, cladding, gap heat transfer, "
            "and frictionless normal contact\n");

        const fuelsim::M1Parameters parameters = make_reference_case();
        const fuelsim::M1Problem problem(parameters);
        const fuelsim::M1LoadStepper load_stepper;
        fuelsim::SolverOptions options;
        options.maximum_iterations = 50;
        constexpr std::size_t load_steps = 20;
        const fuelsim::M1LoadStepResult continuation =
            load_stepper.solve(parameters, load_steps, options);
        const fuelsim::SolveResult& result = continuation.solve;

        const fuelsim::StructuredRzMesh& fuel = problem.fuel_mesh();
        const fuelsim::StructuredRzMesh& cladding = problem.cladding_mesh();
        const fuelsim::DofMap& dofs = problem.dof_map();
        const std::size_t mid_z = fuel.axial_elements() / 2;
        const std::size_t fuel_axis_mid = fuel.node_id(0, mid_z);
        const std::size_t fuel_outer_mid =
            fuel.node_id(fuel.radial_elements(), mid_z);
        const std::size_t fuel_axis_top =
            fuel.node_id(0, fuel.axial_elements());
        const std::size_t cladding_inner_mid = cladding.node_id(0, mid_z);
        const fuelsim::InterfaceSummary interface =
            problem.summarize_interface(result.state);

        std::cout << std::boolalpha << std::scientific << std::setprecision(12);
        std::cout << "converged="
                  << (continuation.completed && result.converged) << '\n';
        std::cout << "convergence_reason="
                  << fuelsim::petsc_convergence_reason_name(
                         result.convergence_reason)
                  << '\n';
        std::cout << "load_steps_completed=" << continuation.completed_steps
                  << '\n';
        std::cout << "nonlinear_iterations_last_step="
                  << result.nonlinear_iterations << '\n';
        std::cout << "residual_norm_last_step=" << result.residual_norm << '\n';
        std::cout << "temperature_center="
                  << result.state[dofs.temperature(
                         problem.fuel_global_node(fuel_axis_mid))]
                  << '\n';
        std::cout << "temperature_fuel_surface="
                  << result.state[dofs.temperature(
                         problem.fuel_global_node(fuel_outer_mid))]
                  << '\n';
        std::cout << "temperature_cladding_inner="
                  << result.state[dofs.temperature(
                         problem.cladding_global_node(cladding_inner_mid))]
                  << '\n';
        std::cout << "radial_displacement_fuel_surface="
                  << result.state[dofs.radial_displacement(
                         problem.fuel_global_node(fuel_outer_mid))]
                  << '\n';
        std::cout << "radial_displacement_cladding_inner="
                  << result.state[dofs.radial_displacement(
                         problem.cladding_global_node(cladding_inner_mid))]
                  << '\n';
        std::cout << "axial_displacement_fuel_top="
                  << result.state[dofs.axial_displacement(
                         problem.fuel_global_node(fuel_axis_top))]
                  << '\n';
        std::cout << "minimum_gap=" << interface.minimum_gap << '\n';
        std::cout << "maximum_gap=" << interface.maximum_gap << '\n';
        std::cout << "maximum_contact_pressure="
                  << interface.maximum_contact_pressure << '\n';
        std::cout << "total_heat_rate=" << interface.total_heat_rate << '\n';
        std::cout << "total_contact_force=" << interface.total_contact_force
                  << '\n';

        return continuation.completed && result.converged ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim failed: " << error.what() << '\n';
        return 1;
    }
}
