#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/m1_problem.hpp"
#include "fuelsim/m1_solver.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem.hpp"

#include <algorithm>
#include <array>
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

double metric_relative_error(double actual, double expected) {
    return std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-30);
}

fuelsim::ThermoelasticProperties constant_material(double conductivity,
                                                   double thermal_expansion) {
    return {
        0.0, conductivity, 75.0e9, 0.3, thermal_expansion, 600.0,
    };
}

bool test_thermal_cylinder() {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double conductivity = 4.0;
    constexpr double heat_source = 2.0e8;
    constexpr double outer_temperature = 600.0;

    const fuelsim::M0Parameters parameters = {
        0.0,
        radius,
        length,
        32,
        2,
        constant_material(conductivity, 0.0),
        heat_source,
        outer_temperature,
        outer_temperature,
        0.0,
        0.0,
    };

    fuelsim::M0Problem problem(parameters);
    fuelsim::PetscSequentialSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    bool passed = check(result.converged, "thermal cylinder SNES converged");
    double maximum_scaled_error = 0.0;
    const double center_rise =
        heat_source * radius * radius / (4.0 * conductivity);
    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const double r = problem.mesh().nodes()[node].r;
        const double expected =
            outer_temperature +
            heat_source * (radius * radius - r * r) / (4.0 * conductivity);
        const double actual = result.state[problem.dof_map().temperature(node)];
        maximum_scaled_error = std::max(
            maximum_scaled_error, std::abs(actual - expected) / center_rise);
    }
    passed = check(maximum_scaled_error < 1.0e-3,
                   "thermal cylinder temperature error is below 0.1%; actual=" +
                       std::to_string(maximum_scaled_error)) &&
             passed;
    std::cout << "thermal_cylinder_maximum_scaled_error="
              << maximum_scaled_error << '\n';
    return passed;
}

bool test_free_thermal_expansion() {
    constexpr double radius = 0.004;
    constexpr double length = 0.01;
    constexpr double alpha = 1.0e-5;
    constexpr double temperature = 700.0;
    constexpr double temperature_change = 100.0;

    const fuelsim::M0Parameters parameters = {
        0.0, radius,      length, 8,   4,   constant_material(4.0, alpha),
        0.0, temperature, 600.0,  0.0, 0.0,
    };

    const fuelsim::M0Problem problem(parameters);
    fuelsim::PetscSequentialSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    bool passed = check(result.converged, "free expansion SNES converged");
    double maximum_temperature_error = 0.0;
    double maximum_displacement_error = 0.0;
    const double displacement_scale = alpha * temperature_change * length;
    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const fuelsim::RzPoint& point = problem.mesh().nodes()[node];
        const double actual_temperature =
            result.state[problem.dof_map().temperature(node)];
        const double actual_radial =
            result.state[problem.dof_map().radial_displacement(node)];
        const double actual_axial =
            result.state[problem.dof_map().axial_displacement(node)];
        maximum_temperature_error =
            std::max(maximum_temperature_error,
                     std::abs(actual_temperature - temperature));
        maximum_displacement_error = std::max(
            maximum_displacement_error,
            std::abs(actual_radial - alpha * temperature_change * point.r));
        maximum_displacement_error = std::max(
            maximum_displacement_error,
            std::abs(actual_axial - alpha * temperature_change * point.z));
    }

    double maximum_stress = 0.0;
    for (std::size_t element = 0; element < problem.element_count();
         ++element) {
        const fuelsim::LocalValues state =
            problem.element_state(element, result.state);
        const auto stresses = problem.kernel().stress_values(
            problem.element_geometry(element), state);
        for (const fuelsim::AxisymmetricStressValues& stress : stresses) {
            maximum_stress = std::max(maximum_stress, std::abs(stress.rr));
            maximum_stress = std::max(maximum_stress, std::abs(stress.zz));
            maximum_stress = std::max(maximum_stress, std::abs(stress.hoop));
            maximum_stress = std::max(maximum_stress, std::abs(stress.rz));
        }
    }

    passed = check(maximum_temperature_error < 1.0e-9,
                   "free expansion temperature is uniform") &&
             passed;
    passed = check(maximum_displacement_error / displacement_scale < 1.0e-8,
                   "free expansion displacement matches analytic field") &&
             passed;
    passed = check(maximum_stress < 100.0,
                   "free expansion stress is below 100 Pa") &&
             passed;
    std::cout << "free_expansion_maximum_temperature_error="
              << maximum_temperature_error << '\n';
    std::cout << "free_expansion_maximum_displacement_relative_error="
              << maximum_displacement_error / displacement_scale << '\n';
    std::cout << "free_expansion_maximum_stress=" << maximum_stress << '\n';
    return passed;
}

bool test_lame_open_ended_cylinder() {
    constexpr double inner_radius = 0.004;
    constexpr double outer_radius = 0.005;
    constexpr double length = 0.01;
    constexpr double pressure = 1.0e6;
    constexpr double young_modulus = 75.0e9;
    constexpr double poisson_ratio = 0.3;

    fuelsim::ThermoelasticProperties material = constant_material(4.0, 0.0);
    material.young_modulus = young_modulus;
    material.poisson_ratio = poisson_ratio;

    const fuelsim::M0Parameters parameters = {
        inner_radius, outer_radius, length, 48,       2,   material,
        0.0,          600.0,        600.0,  pressure, 0.0,
    };

    const fuelsim::M0Problem problem(parameters);
    fuelsim::PetscSequentialSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    bool passed = check(result.converged, "Lame cylinder SNES converged");
    const double A =
        pressure * inner_radius * inner_radius /
        (outer_radius * outer_radius - inner_radius * inner_radius);
    const double B =
        pressure * inner_radius * inner_radius * outer_radius * outer_radius /
        (outer_radius * outer_radius - inner_radius * inner_radius);
    const double axial_strain = -2.0 * poisson_ratio * A / young_modulus;

    double maximum_radial_relative_error = 0.0;
    double maximum_axial_relative_error = 0.0;
    const double radial_scale = ((1.0 - poisson_ratio) * A * inner_radius +
                                 (1.0 + poisson_ratio) * B / inner_radius) /
                                young_modulus;
    const double axial_scale = std::abs(axial_strain * length);

    for (std::size_t node = 0; node < problem.mesh().nodes().size(); ++node) {
        const fuelsim::RzPoint& point = problem.mesh().nodes()[node];
        const double expected_radial = ((1.0 - poisson_ratio) * A * point.r +
                                        (1.0 + poisson_ratio) * B / point.r) /
                                       young_modulus;
        const double expected_axial = axial_strain * point.z;
        const double actual_radial =
            result.state[problem.dof_map().radial_displacement(node)];
        const double actual_axial =
            result.state[problem.dof_map().axial_displacement(node)];
        maximum_radial_relative_error =
            std::max(maximum_radial_relative_error,
                     std::abs(actual_radial - expected_radial) / radial_scale);
        if (axial_scale > 0.0) {
            maximum_axial_relative_error =
                std::max(maximum_axial_relative_error,
                         std::abs(actual_axial - expected_axial) / axial_scale);
        }
    }

    passed = check(maximum_radial_relative_error < 3.0e-3,
                   "Lame radial displacement error is below 0.3%") &&
             passed;
    passed = check(maximum_axial_relative_error < 3.0e-3,
                   "Lame axial displacement error is below 0.3%") &&
             passed;
    std::cout << "lame_maximum_radial_relative_error="
              << maximum_radial_relative_error << '\n';
    std::cout << "lame_maximum_axial_relative_error="
              << maximum_axial_relative_error << '\n';
    return passed;
}

bool test_moose_reference(const std::string& mesh_path) {
    const fuelsim::M0Parameters parameters = {
        0.0,
        0.00412,
        0.010,
        40,
        10,
        {
            3824.0,
            0.61,
            2.0e11,
            0.316,
            10.0e-6,
            600.0,
        },
        2.0e8,
        600.0,
        600.0,
        0.0,
        0.0,
    };

    const fuelsim::UnstructuredQuad4Mesh imported =
        fuelsim::ExodusMeshIo::read_quad4(mesh_path);
    const fuelsim::StructuredRzMesh imported_mesh =
        fuelsim::StructuredRzMesh::from_unstructured_block(
            imported, 0,
            {"fuel_left", "fuel_right", "fuel_bottom", "fuel_top"});
    const fuelsim::M0Problem problem(parameters, imported_mesh);
    fuelsim::PetscSequentialSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    const fuelsim::StructuredRzMesh& mesh = problem.mesh();
    const fuelsim::DofMap& dofs = problem.dof_map();
    const std::size_t mid_z = mesh.axial_elements() / 2;
    const std::size_t axis_mid = mesh.node_id(0, mid_z);
    const std::size_t outer_mid = mesh.node_id(mesh.radial_elements(), mid_z);
    const std::size_t outer_top =
        mesh.node_id(mesh.radial_elements(), mesh.axial_elements());
    const std::size_t axis_top = mesh.node_id(0, mesh.axial_elements());

    // Generated by verification/moose/m0_simple_fuel_rz.i with July/MOOSE
    // commit 93b11698be and PETSc 3.25.2.
    constexpr double expected_temperature_center = 733.4201407826;
    constexpr double expected_radial_outer_mid = 2.5827190582385e-6;
    constexpr double expected_radial_outer_top = 3.5772306219751e-6;
    constexpr double expected_axial_axis_top = 8.7881251481059e-6;

    const double temperature_center = result.state[dofs.temperature(axis_mid)];
    const double radial_outer_mid =
        result.state[dofs.radial_displacement(outer_mid)];
    const double radial_outer_top =
        result.state[dofs.radial_displacement(outer_top)];
    const double axial_axis_top =
        result.state[dofs.axial_displacement(axis_top)];
    const double temperature_error =
        metric_relative_error(temperature_center, expected_temperature_center);
    const double radial_mid_error =
        metric_relative_error(radial_outer_mid, expected_radial_outer_mid);
    const double radial_top_error =
        metric_relative_error(radial_outer_top, expected_radial_outer_top);
    const double axial_top_error =
        metric_relative_error(axial_axis_top, expected_axial_axis_top);

    bool passed =
        check(result.converged, "MOOSE reference case SNES converged");
    passed = check(temperature_error < 1.0e-10,
                   "MOOSE temperature_center relative error is below 1e-10") &&
             passed;
    passed = check(radial_mid_error < 1.0e-10,
                   "MOOSE radial_displacement_outer_mid relative error is "
                   "below 1e-10") &&
             passed;
    passed = check(radial_top_error < 1.0e-10,
                   "MOOSE radial_displacement_outer_top relative error is "
                   "below 1e-10") &&
             passed;
    passed = check(axial_top_error < 1.0e-10,
                   "MOOSE axial_displacement_axis_top relative error is below "
                   "1e-10") &&
             passed;
    std::cout << "moose_temperature_center_relative_error=" << temperature_error
              << '\n';
    std::cout << "moose_radial_outer_mid_relative_error=" << radial_mid_error
              << '\n';
    std::cout << "moose_radial_outer_top_relative_error=" << radial_top_error
              << '\n';
    std::cout << "moose_axial_axis_top_relative_error=" << axial_top_error
              << '\n';
    return passed;
}

bool test_m1_open_gap_analytic_thermal() {
    constexpr double fuel_radius = 0.004;
    constexpr double cladding_inner_radius = 0.0041;
    constexpr double cladding_outer_radius = 0.0046;
    constexpr double length = 0.010;
    constexpr double fuel_conductivity = 4.0;
    constexpr double cladding_conductivity = 16.0;
    constexpr double gap_conductivity = 0.4;
    constexpr double heat_source = 1.0e8;
    constexpr double outer_temperature = 600.0;

    const fuelsim::M1Parameters parameters = {
        fuel_radius,
        cladding_inner_radius,
        cladding_outer_radius,
        length,
        length,
        32,
        8,
        2,
        constant_material(fuel_conductivity, 0.0),
        constant_material(cladding_conductivity, 0.0),
        heat_source,
        outer_temperature,
        outer_temperature,
        gap_conductivity,
        1.0e-6,
        1.0e14,
    };
    const fuelsim::M1Problem problem(parameters);
    fuelsim::PetscSequentialSolver solver;
    const fuelsim::SolveResult result =
        solver.solve(problem, problem.initial_state());

    const double cladding_rise =
        heat_source * fuel_radius * fuel_radius /
        (2.0 * cladding_conductivity) *
        std::log(cladding_outer_radius / cladding_inner_radius);
    const double gap_rise = heat_source * fuel_radius *
                            (cladding_inner_radius - fuel_radius) /
                            (2.0 * gap_conductivity);
    const double fuel_rise =
        heat_source * fuel_radius * fuel_radius / (4.0 * fuel_conductivity);
    const double expected_cladding_inner = outer_temperature + cladding_rise;
    const double expected_fuel_surface = expected_cladding_inner + gap_rise;
    const double expected_center = expected_fuel_surface + fuel_rise;

    const std::size_t axial_mid = problem.fuel_mesh().axial_elements() / 2;
    const std::size_t fuel_center_local =
        problem.fuel_mesh().node_id(0, axial_mid);
    const std::size_t fuel_surface_local = problem.fuel_mesh().node_id(
        problem.fuel_mesh().radial_elements(), axial_mid);
    const std::size_t cladding_inner_local =
        problem.cladding_mesh().node_id(0, axial_mid);
    const fuelsim::DofMap& dofs = problem.dof_map();

    const double actual_center = result.state[dofs.temperature(
        problem.fuel_global_node(fuel_center_local))];
    const double actual_fuel_surface = result.state[dofs.temperature(
        problem.fuel_global_node(fuel_surface_local))];
    const double actual_cladding_inner = result.state[dofs.temperature(
        problem.cladding_global_node(cladding_inner_local))];
    const double temperature_scale = expected_center - outer_temperature;
    const double center_error =
        std::abs(actual_center - expected_center) / temperature_scale;
    const double fuel_surface_error =
        std::abs(actual_fuel_surface - expected_fuel_surface) /
        temperature_scale;
    const double cladding_inner_error =
        std::abs(actual_cladding_inner - expected_cladding_inner) /
        temperature_scale;

    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.state);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double expected_heat_rate =
        heat_source * pi * fuel_radius * fuel_radius * length;
    const double heat_balance_error =
        metric_relative_error(interface.total_heat_rate, expected_heat_rate);

    bool passed = check(result.converged, "M1 open-gap thermal SNES converged");
    passed = check(center_error < 1.0e-3,
                   "M1 analytic center temperature error is below 0.1%") &&
             passed;
    passed = check(fuel_surface_error < 1.0e-3,
                   "M1 analytic fuel-surface temperature error is below "
                   "0.1%") &&
             passed;
    passed = check(cladding_inner_error < 1.0e-3,
                   "M1 analytic cladding-inner temperature error is below "
                   "0.1%") &&
             passed;
    passed = check(interface.minimum_gap > 0.0 &&
                       interface.maximum_contact_pressure == 0.0,
                   "M1 analytic thermal case remains out of contact") &&
             passed;
    passed = check(heat_balance_error < 1.0e-9,
                   "M1 interface heat rate balances generated power") &&
             passed;

    std::cout << "m1_analytic_center_temperature_scaled_error=" << center_error
              << '\n';
    std::cout << "m1_analytic_fuel_surface_scaled_error=" << fuel_surface_error
              << '\n';
    std::cout << "m1_analytic_cladding_inner_scaled_error="
              << cladding_inner_error << '\n';
    std::cout << "m1_analytic_heat_balance_relative_error="
              << heat_balance_error << '\n';
    return passed;
}

fuelsim::M1Parameters m1_reference_parameters() {
    return {
        0.004120,
        0.004122,
        0.004692,
        0.010,
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

bool test_m1_closed_gap_end_to_end(const std::string& mesh_path) {
    const fuelsim::M1Parameters target_parameters = m1_reference_parameters();
    const fuelsim::UnstructuredQuad4Mesh imported =
        fuelsim::ExodusMeshIo::read_quad4(mesh_path);
    const fuelsim::StructuredRzMesh fuel =
        fuelsim::StructuredRzMesh::from_unstructured_block(
            imported, "fuel",
            {"fuel_left", "fuel_right", "fuel_bottom", "fuel_top"});
    const fuelsim::StructuredRzMesh cladding =
        fuelsim::StructuredRzMesh::from_unstructured_block(
            imported, "clad",
            {"clad_left", "clad_right", "clad_bottom", "clad_top"});
    fuelsim::SolverOptions options;
    options.maximum_iterations = 50;

    constexpr std::size_t load_steps = 20;
    const fuelsim::M1LoadStepper load_stepper;
    const fuelsim::M1LoadStepResult continuation = load_stepper.solve(
        target_parameters, fuel, cladding, load_steps, options);
    const fuelsim::SolveResult& result = continuation.solve;

    const fuelsim::M1Problem problem(target_parameters, fuel, cladding);

    const std::size_t axial_mid = problem.fuel_mesh().axial_elements() / 2;
    const std::size_t fuel_center_local =
        problem.fuel_mesh().node_id(0, axial_mid);
    const std::size_t fuel_surface_local = problem.fuel_mesh().node_id(
        problem.fuel_mesh().radial_elements(), axial_mid);
    const std::size_t cladding_inner_local =
        problem.cladding_mesh().node_id(0, axial_mid);
    const std::size_t fuel_axis_top_local =
        problem.fuel_mesh().node_id(0, problem.fuel_mesh().axial_elements());
    const fuelsim::DofMap& dofs = problem.dof_map();

    const double temperature_center = result.state[dofs.temperature(
        problem.fuel_global_node(fuel_center_local))];
    const double temperature_fuel_surface = result.state[dofs.temperature(
        problem.fuel_global_node(fuel_surface_local))];
    const double temperature_cladding_inner = result.state[dofs.temperature(
        problem.cladding_global_node(cladding_inner_local))];
    const double radial_fuel_surface = result.state[dofs.radial_displacement(
        problem.fuel_global_node(fuel_surface_local))];
    const double radial_cladding_inner = result.state[dofs.radial_displacement(
        problem.cladding_global_node(cladding_inner_local))];
    const double axial_fuel_top = result.state[dofs.axial_displacement(
        problem.fuel_global_node(fuel_axis_top_local))];
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(result.state);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(result.state);

    constexpr double pi = 3.141592653589793238462643383279502884;
    const fuelsim::M1Parameters& parameters = problem.parameters();
    const double expected_heat_rate =
        parameters.volumetric_heat_source * pi * parameters.fuel_radius *
        parameters.fuel_radius * parameters.fuel_length;
    const double heat_balance_error =
        metric_relative_error(interface.total_heat_rate, expected_heat_rate);

    // Generated by verification/moose/m1_fuel_cladding_gap_rz.i with
    // area-normalized node-face penalty contact and a cladding height of
    // 10.020 mm.
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

    const double moose_temperature_center_error =
        metric_relative_error(temperature_center, expected_temperature_center);
    const double moose_temperature_fuel_surface_error = metric_relative_error(
        temperature_fuel_surface, expected_temperature_fuel_surface);
    const double moose_temperature_cladding_inner_error = metric_relative_error(
        temperature_cladding_inner, expected_temperature_cladding_inner);
    const double moose_radial_fuel_surface_error = metric_relative_error(
        radial_fuel_surface, expected_radial_fuel_surface);
    const double moose_radial_cladding_inner_error = metric_relative_error(
        radial_cladding_inner, expected_radial_cladding_inner);
    const double moose_axial_fuel_top_error =
        metric_relative_error(axial_fuel_top, expected_axial_fuel_top);
    double pressure_difference_squared = 0.0;
    double pressure_reference_squared = 0.0;
    double pressure_maximum_point_error = 0.0;
    for (std::size_t node = 0; node < expected_contact_pressure.size();
         ++node) {
        const double difference =
            contact_nodes.at(node).pressure - expected_contact_pressure[node];
        pressure_difference_squared += difference * difference;
        pressure_reference_squared +=
            expected_contact_pressure[node] * expected_contact_pressure[node];
        pressure_maximum_point_error =
            std::max(pressure_maximum_point_error,
                     metric_relative_error(contact_nodes.at(node).pressure,
                                           expected_contact_pressure[node]));
    }
    const double moose_contact_pressure_relative_l2 =
        std::sqrt(pressure_difference_squared / pressure_reference_squared);
    const double moose_total_contact_force_error = metric_relative_error(
        interface.total_contact_force, expected_total_contact_force);

    bool passed = check(continuation.completed && result.converged,
                        "M1 closed-gap load continuation converged");
    passed = check(interface.minimum_gap < 0.0,
                   "M1 reference case closes the fuel-cladding gap") &&
             passed;
    passed = check(interface.maximum_contact_pressure > 0.0,
                   "M1 reference case develops contact pressure") &&
             passed;
    passed = check(interface.projected_contact_nodes ==
                       target_parameters.axial_elements + 1,
                   "M1 taller cladding contains every NTS projection") &&
             passed;
    passed = check(continuation.aggregate_timing.workspace_setups == 1,
                   "M1 load path creates one reusable PETSc workspace") &&
             passed;
    passed = check(continuation.aggregate_timing.solve_calls == load_steps,
                   "M1 load path reuses the workspace for all load steps") &&
             passed;
    passed =
        check(continuation.aggregate_timing.residual_evaluations > 0 &&
                  continuation.aggregate_timing.jacobian_evaluations > 0 &&
                  continuation.aggregate_timing.nonlinear_solve_seconds > 0.0 &&
                  continuation.total_seconds > 0.0,
              "M1 load path reports internal solver timing and evaluations") &&
        passed;
    passed = check(heat_balance_error < 1.0e-8,
                   "M1 closed-gap interface heat balances generated power") &&
             passed;
    passed =
        check(temperature_center > temperature_fuel_surface &&
                  temperature_fuel_surface > temperature_cladding_inner &&
                  temperature_cladding_inner > parameters.outer_temperature,
              "M1 temperatures decrease from fuel center to cladding "
              "outer surface") &&
        passed;
    passed = check(moose_temperature_center_error < 1.0e-2,
                   "M1 MOOSE center temperature error is below 1%") &&
             passed;
    passed = check(moose_temperature_fuel_surface_error < 1.0e-2,
                   "M1 MOOSE fuel-surface temperature error is below 1%") &&
             passed;
    passed = check(moose_temperature_cladding_inner_error < 1.0e-2,
                   "M1 MOOSE cladding-inner temperature error is below 1%") &&
             passed;
    passed = check(moose_radial_fuel_surface_error < 1.0e-2,
                   "M1 MOOSE fuel-surface displacement error is below 1%") &&
             passed;
    passed = check(moose_radial_cladding_inner_error < 1.0e-2,
                   "M1 MOOSE cladding-inner displacement error is below 1%") &&
             passed;
    passed = check(moose_axial_fuel_top_error < 1.0e-2,
                   "M1 MOOSE fuel-top axial displacement error is below 1%") &&
             passed;
    passed = check(moose_contact_pressure_relative_l2 < 1.0e-2,
                   "M1 MOOSE nodal contact-pressure L2 error is below 1%") &&
             passed;
    passed = check(pressure_maximum_point_error < 1.0e-2,
                   "M1 MOOSE maximum nodal pressure error is below 1%") &&
             passed;
    passed = check(moose_total_contact_force_error < 1.0e-2,
                   "M1 MOOSE total contact-force error is below 1%") &&
             passed;

    std::cout << "m1_temperature_center=" << temperature_center << '\n';
    std::cout << "m1_temperature_fuel_surface=" << temperature_fuel_surface
              << '\n';
    std::cout << "m1_temperature_cladding_inner=" << temperature_cladding_inner
              << '\n';
    std::cout << "m1_radial_fuel_surface=" << radial_fuel_surface << '\n';
    std::cout << "m1_radial_cladding_inner=" << radial_cladding_inner << '\n';
    std::cout << "m1_axial_fuel_top=" << axial_fuel_top << '\n';
    std::cout << "m1_minimum_gap=" << interface.minimum_gap << '\n';
    std::cout << "m1_maximum_gap=" << interface.maximum_gap << '\n';
    std::cout << "m1_minimum_contact_gap=" << interface.minimum_contact_gap
              << '\n';
    std::cout << "m1_maximum_contact_pressure="
              << interface.maximum_contact_pressure << '\n';
    std::cout << "m1_total_heat_rate=" << interface.total_heat_rate << '\n';
    std::cout << "m1_total_contact_force=" << interface.total_contact_force
              << '\n';
    std::cout << "m1_projected_contact_nodes="
              << interface.projected_contact_nodes << '\n';
    std::cout << "m1_active_contact_nodes=" << interface.active_contact_nodes
              << '\n';
    std::cout << "m1_total_nonlinear_iterations="
              << continuation.total_nonlinear_iterations << '\n';
    std::cout << "m1_timing_problem_setup="
              << continuation.problem_setup_seconds << '\n';
    std::cout << "m1_timing_solver_setup="
              << continuation.aggregate_timing.setup_seconds << '\n';
    std::cout << "m1_timing_nonlinear_solve="
              << continuation.aggregate_timing.nonlinear_solve_seconds << '\n';
    std::cout << "m1_timing_residual_callbacks="
              << continuation.aggregate_timing.residual_callback_seconds
              << '\n';
    std::cout << "m1_timing_jacobian_callbacks="
              << continuation.aggregate_timing.jacobian_callback_seconds
              << '\n';
    std::cout << "m1_timing_total=" << continuation.total_seconds << '\n';
    std::cout << "m1_heat_balance_relative_error=" << heat_balance_error
              << '\n';
    std::cout << "m1_moose_temperature_center_relative_error="
              << moose_temperature_center_error << '\n';
    std::cout << "m1_moose_temperature_fuel_surface_relative_error="
              << moose_temperature_fuel_surface_error << '\n';
    std::cout << "m1_moose_temperature_cladding_inner_relative_error="
              << moose_temperature_cladding_inner_error << '\n';
    std::cout << "m1_moose_radial_fuel_surface_relative_error="
              << moose_radial_fuel_surface_error << '\n';
    std::cout << "m1_moose_radial_cladding_inner_relative_error="
              << moose_radial_cladding_inner_error << '\n';
    std::cout << "m1_moose_axial_fuel_top_relative_error="
              << moose_axial_fuel_top_error << '\n';
    std::cout << "m1_moose_contact_pressure_relative_l2="
              << moose_contact_pressure_relative_l2 << '\n';
    std::cout << "m1_moose_contact_pressure_maximum_point_relative_error="
              << pressure_maximum_point_error << '\n';
    std::cout << "m1_moose_total_contact_force_relative_error="
              << moose_total_contact_force_error << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_solver_tests <m0_mesh.e> <m1_mesh.e>\n";
        return 2;
    }

    try {
        const std::string m0_mesh_path = argv[1];
        const std::string m1_mesh_path = argv[2];
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M0 and M1 numerical acceptance tests\n");

        bool passed = true;
        passed = test_thermal_cylinder() && passed;
        passed = test_free_thermal_expansion() && passed;
        passed = test_lame_open_ended_cylinder() && passed;
        passed = test_moose_reference(m0_mesh_path) && passed;
        passed = test_m1_open_gap_analytic_thermal() && passed;
        passed = test_m1_closed_gap_end_to_end(m1_mesh_path) && passed;
        if (!passed)
            return 1;

        std::cout << "[PASS] fuelsim M0 and M1 solver acceptance tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] solver tests raised: " << error.what() << '\n';
        return 1;
    }
}
