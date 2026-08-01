#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "support/steady_fuel_cladding_problem.hpp"
#include "support/steady_fuel_cladding_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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

double relative_error(double actual, double expected) {
    return std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-30);
}

fuelsim::SteadyFuelCladdingParameters reference_parameters() {
    return {
        0.004120,
        0.004122,
        0.004692,
        0.010000,
        0.010020,
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

bool run_comparison(const char* mesh_path) {
    const fuelsim::UnstructuredQuad4Mesh imported =
        fuelsim::ExodusMeshIo::read_quad4(mesh_path);
    bool passed = check(imported.nodes().size() == 528,
                        "MOOSE Exodus mesh has 528 nodes");
    passed = check(imported.elements().size() == 460,
                   "MOOSE Exodus mesh has 460 Quad4 elements") &&
             passed;
    passed = check(imported.element_block("fuel").id == 0 &&
                       imported.element_block("clad").id == 1,
                   "MOOSE Exodus fuel/clad block IDs are preserved") &&
             passed;
    passed = check(imported.node_set("fuel_right").nodes.size() == 11 &&
                       imported.node_set("clad_left").nodes.size() == 11,
                   "MOOSE Exodus contact node sets are preserved") &&
             passed;
    passed = check(imported.side_set("fuel_right").sides.size() == 10 &&
                       imported.side_set("clad_left").sides.size() == 10,
                   "MOOSE Exodus contact side sets are preserved") &&
             passed;

    const fuelsim::StructuredRzMesh fuel =
        fuelsim::StructuredRzMesh::from_unstructured_block(
            imported, "fuel",
            {"fuel_left", "fuel_right", "fuel_bottom", "fuel_top"});
    const fuelsim::StructuredRzMesh cladding =
        fuelsim::StructuredRzMesh::from_unstructured_block(
            imported, "clad",
            {"clad_left", "clad_right", "clad_bottom", "clad_top"});
    passed =
        check(fuel.nodes().size() == 451 && fuel.elements().size() == 400 &&
                  cladding.nodes().size() == 77 &&
                  cladding.elements().size() == 60,
              "MOOSE blocks convert to the expected fuelsim meshes") &&
        passed;

    const fuelsim::SteadyFuelCladdingParameters parameters =
        reference_parameters();
    fuelsim::SolverOptions options;
    options.maximum_iterations = 50;
    constexpr std::size_t load_steps = 20;
    const fuelsim::SteadyFuelCladdingLoadStepper load_stepper;
    const fuelsim::SteadyFuelCladdingLoadResult continuation =
        load_stepper.solve(parameters, fuel, cladding, load_steps, options);
    passed = check(continuation.completed && continuation.solve.converged,
                   "M1 solve on the MOOSE Exodus mesh converged") &&
             passed;

    const fuelsim::SteadyFuelCladdingProblem problem(parameters, fuel,
                                                     cladding);
    const fuelsim::SolveResult& result = continuation.solve;
    const std::size_t axial_mid = fuel.axial_elements() / 2;
    const std::size_t fuel_center = fuel.node_id(0, axial_mid);
    const std::size_t fuel_surface =
        fuel.node_id(fuel.radial_elements(), axial_mid);
    const std::size_t cladding_inner = cladding.node_id(0, axial_mid);
    const std::size_t fuel_axis_top = fuel.node_id(0, fuel.axial_elements());
    const fuelsim::DofMap& dofs = problem.dof_map();

    const double temperature_center =
        result.state[dofs.temperature(problem.fuel_global_node(fuel_center))];
    const double temperature_fuel_surface =
        result.state[dofs.temperature(problem.fuel_global_node(fuel_surface))];
    const double temperature_cladding_inner = result.state[dofs.temperature(
        problem.cladding_global_node(cladding_inner))];
    const double radial_fuel_surface = result.state[dofs.radial_displacement(
        problem.fuel_global_node(fuel_surface))];
    const double radial_cladding_inner = result.state[dofs.radial_displacement(
        problem.cladding_global_node(cladding_inner))];
    const double axial_fuel_top = result.state[dofs.axial_displacement(
        problem.fuel_global_node(fuel_axis_top))];
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(result.state);
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.state);

    constexpr double expected_temperature_center = 751.10972999158;
    constexpr double expected_temperature_fuel_surface = 614.75750803533;
    constexpr double expected_temperature_cladding_inner = 613.72910879010;
    constexpr double expected_radial_fuel_surface = 3.2141740355603e-6;
    constexpr double expected_radial_cladding_inner = 1.1901972438575e-6;
    constexpr double expected_axial_fuel_top = 1.0531979661574e-5;
    constexpr double expected_total_contact_force = 663.8896691615588;
    constexpr std::array<double, 11> expected_contact_pressure = {
        2547100.2219749, 2545172.8286605, 2532649.5016160, 2507500.3516350,
        2465844.2361049, 2413364.3996226, 2413608.4076008, 2557256.2070568,
        2243329.5849293, 1643992.4525044, 6029728.9198301,
    };

    const std::array<double, 6> point_errors = {
        relative_error(temperature_center, expected_temperature_center),
        relative_error(temperature_fuel_surface,
                       expected_temperature_fuel_surface),
        relative_error(temperature_cladding_inner,
                       expected_temperature_cladding_inner),
        relative_error(radial_fuel_surface, expected_radial_fuel_surface),
        relative_error(radial_cladding_inner, expected_radial_cladding_inner),
        relative_error(axial_fuel_top, expected_axial_fuel_top),
    };
    const double maximum_point_field_error =
        *std::max_element(point_errors.begin(), point_errors.end());

    double pressure_difference_squared = 0.0;
    double pressure_reference_squared = 0.0;
    double maximum_pressure_point_error = 0.0;
    for (std::size_t node = 0; node < expected_contact_pressure.size();
         ++node) {
        const double difference =
            contact_nodes.at(node).pressure - expected_contact_pressure[node];
        pressure_difference_squared += difference * difference;
        pressure_reference_squared +=
            expected_contact_pressure[node] * expected_contact_pressure[node];
        maximum_pressure_point_error =
            std::max(maximum_pressure_point_error,
                     relative_error(contact_nodes.at(node).pressure,
                                    expected_contact_pressure[node]));
    }
    const double pressure_relative_l2 =
        std::sqrt(pressure_difference_squared / pressure_reference_squared);
    const double contact_force_error = relative_error(
        interface.total_contact_force, expected_total_contact_force);

    passed = check(maximum_point_field_error < 1.0e-2,
                   "MOOSE Exodus temperature/displacement points are within "
                   "1%") &&
             passed;
    passed = check(pressure_relative_l2 < 1.0e-2,
                   "MOOSE Exodus contact-pressure relative L2 is within 1%") &&
             passed;
    passed = check(maximum_pressure_point_error < 1.0e-2,
                   "MOOSE Exodus pointwise contact pressure is within 1%") &&
             passed;
    passed = check(contact_force_error < 1.0e-2,
                   "MOOSE Exodus total contact force is within 1%") &&
             passed;
    passed = check(interface.projected_contact_nodes == 11 &&
                       interface.active_contact_nodes == 11,
                   "All MOOSE Exodus contact nodes are projected and active") &&
             passed;

    std::cout << "m1_exodus_maximum_temperature_displacement_point_error="
              << maximum_point_field_error << '\n';
    std::cout << "m1_exodus_contact_pressure_relative_l2="
              << pressure_relative_l2 << '\n';
    std::cout << "m1_exodus_contact_pressure_maximum_point_error="
              << maximum_pressure_point_error << '\n';
    std::cout << "m1_exodus_total_contact_force_error=" << contact_force_error
              << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m1_exodus_moose_tests <mesh.e>\n";
        return 2;
    }

    try {
        const std::string mesh_path = argv[1];
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M1 MOOSE Exodus mesh comparison\n");
        if (!run_comparison(mesh_path.c_str()))
            return 1;
        std::cout << "[PASS] M1 MOOSE Exodus mesh comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
