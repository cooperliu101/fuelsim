#include "fuelsim/m2_problem.hpp"
#include "fuelsim/m2_solver.hpp"
#include "fuelsim/petsc_solver.hpp"

#include <exception>
#include <iomanip>
#include <iostream>

namespace {

fuelsim::M2Parameters make_generic_demonstration_case() {
    const fuelsim::M1Parameters base = {
        0.004120,
        0.004122,
        0.004692,
        0.010,
        0.010020,
        12,
        4,
        4,
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

    const fuelsim::TransientInelasticProperties fuel = {
        10970.0,
        300.0,
        fuelsim::InelasticBehavior::norton_creep,
        {
            1.0e-5,
            1.0e8,
            3.0,
        },
        {
            1.0,
            0.0,
        },
    };
    const fuelsim::TransientInelasticProperties cladding = {
        6500.0,
        330.0,
        fuelsim::InelasticBehavior::norton_creep_j2_plasticity,
        {
            1.0e-5,
            5.0e6,
            3.0,
        },
        {
            5.0e6,
            2.0e9,
        },
    };
    return {base, fuel, cladding};
}

} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim M2 generic demonstration: transient heat capacity, "
            "Norton creep, J2 plasticity, and their coupled update\n");

        fuelsim::M2Problem problem(make_generic_demonstration_case());
        const fuelsim::M2TimeOptions time_options = {
            20.0, 1.0, 0.125, 1.0, 1.0, 0.5, 3, 20.0,
        };
        fuelsim::SolverOptions solver_options;
        solver_options.maximum_iterations = 50;

        const fuelsim::M2TimeStepper time_stepper;
        const fuelsim::M2TransientResult result =
            time_stepper.solve(problem, time_options, solver_options);
        const fuelsim::RegionInelasticSummary fuel_history =
            problem.summarize_fuel_history();
        const fuelsim::RegionInelasticSummary cladding_history =
            problem.summarize_cladding_history();
        const fuelsim::InterfaceSummary interface =
            problem.summarize_interface(result.committed_state);

        std::cout << std::boolalpha << std::scientific << std::setprecision(12);
        std::cout << "completed=" << result.completed << '\n';
        std::cout << "committed_time=" << result.committed_time << '\n';
        std::cout << "accepted_steps=" << result.accepted_steps.size() << '\n';
        std::cout << "total_cutbacks=" << result.total_cutbacks << '\n';
        std::cout << "nonlinear_iterations_total="
                  << result.total_nonlinear_iterations << '\n';
        std::cout << "petsc_workspace_setups="
                  << result.aggregate_timing.workspace_setups << '\n';
        std::cout << "fuel_maximum_equivalent_creep_strain="
                  << fuel_history.maximum_equivalent_creep_strain << '\n';
        std::cout << "cladding_maximum_equivalent_plastic_strain="
                  << cladding_history.maximum_equivalent_plastic_strain << '\n';
        std::cout << "cladding_maximum_equivalent_creep_strain="
                  << cladding_history.maximum_equivalent_creep_strain << '\n';
        std::cout << "minimum_gap=" << interface.minimum_gap << '\n';
        std::cout << "maximum_contact_pressure="
                  << interface.maximum_contact_pressure << '\n';
        std::cout << "projected_contact_nodes="
                  << interface.projected_contact_nodes << '\n';
        std::cout << "active_contact_nodes=" << interface.active_contact_nodes
                  << '\n';

        const bool accepted =
            result.completed && result.accepted_steps.size() == 20 &&
            result.aggregate_timing.workspace_setups == 1 &&
            fuel_history.maximum_equivalent_creep_strain > 0.0 &&
            cladding_history.maximum_equivalent_plastic_strain > 0.0 &&
            cladding_history.maximum_equivalent_creep_strain > 0.0 &&
            interface.projected_contact_nodes ==
                problem.parameters().base.axial_elements + 1;
        return accepted ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim M2 failed: " << error.what() << '\n';
        return 1;
    }
}
