#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ErrorMetrics final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_actual = 0.0;
    double maximum_reference = 0.0;
    double maximum_pointwise_relative = 0.0;

    void add(double actual, double reference) {
        const double difference = actual - reference;
        difference_squared += difference * difference;
        reference_squared += reference * reference;
        maximum_actual = std::max(maximum_actual, std::abs(actual));
        maximum_reference = std::max(maximum_reference, std::abs(reference));
        maximum_pointwise_relative =
            std::max(maximum_pointwise_relative,
                     std::abs(difference) / std::abs(reference));
    }

    double relative_l2() const {
        return std::sqrt(difference_squared / reference_squared);
    }

    double relative_absolute_peak() const {
        return std::abs(maximum_actual - maximum_reference) / maximum_reference;
    }
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool check_metrics(const std::string& name, const ErrorMetrics& metrics,
                   double tolerance) {
    std::cout << name << "_relative_l2=" << metrics.relative_l2() << '\n';
    std::cout << name
              << "_relative_absolute_peak=" << metrics.relative_absolute_peak()
              << '\n';
    std::cout << name << "_maximum_pointwise_relative="
              << metrics.maximum_pointwise_relative << '\n';
    return check(metrics.relative_l2() < tolerance &&
                     metrics.relative_absolute_peak() < tolerance &&
                     metrics.maximum_pointwise_relative < tolerance,
                 name + " three MOOSE error metrics pass");
}

bool same_coordinate(double lhs, double rhs) {
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= 1.0e-12 * scale;
}

std::size_t find_node(const fuelsim::RegionMesh& mesh, double radius,
                      double axial_coordinate) {
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        if (same_coordinate(mesh.nodes()[node].r, radius) &&
            same_coordinate(mesh.nodes()[node].z, axial_coordinate))
            return node;
    }
    throw std::invalid_argument("M1 comparison point is not in the mesh");
}

bool run_comparison(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument(
            "M1 comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    bool passed =
        check(source.nodes().size() == 528 && source.elements().size() == 460,
              "M1 input card reads all MOOSE Exodus entities");
    passed = check(source.element_block("fuel").id == 0 &&
                       source.element_block("clad").id == 1 &&
                       source.side_set("fuel_right").sides.size() == 10 &&
                       source.side_set("clad_left").sides.size() == 10,
                   "M1 MOOSE block and side-set metadata are preserved") &&
             passed;

    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    const fuelsim::SolverOptions options = {
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance, definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem, definition.steady_execution.load_steps, options);
    passed = check(result.completed && result.solve.converged,
                   "M1 input-card load path converged") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "M1 input-card path reuses one PETSc workspace") &&
             passed;

    const std::size_t fuel_region = problem.region_index("fuel");
    const std::size_t clad_region = problem.region_index("cladding");
    const fuelsim::RegionMesh& fuel = problem.region_mesh(fuel_region);
    const fuelsim::RegionMesh& clad = problem.region_mesh(clad_region);
    const std::size_t fuel_center = find_node(fuel, 0.0, 0.005);
    const std::size_t fuel_surface = find_node(fuel, 0.00412, 0.005);
    const std::size_t fuel_top = find_node(fuel, 0.0, 0.010);
    const std::size_t clad_inner = find_node(clad, 0.004122, 0.00501);
    const std::size_t fuel_offset = problem.region_node_offset(fuel_region);
    const std::size_t clad_offset = problem.region_node_offset(clad_region);
    const fuelsim::DofMap& dofs = problem.dof_map();
    const std::vector<double>& state = result.solve.state;

    ErrorMetrics temperature;
    ErrorMetrics radial_displacement;
    ErrorMetrics axial_displacement;
    temperature.add(state[dofs.temperature(fuel_offset + fuel_center)],
                    751.10972999158);
    temperature.add(state[dofs.temperature(fuel_offset + fuel_surface)],
                    614.75750803533);
    temperature.add(state[dofs.temperature(clad_offset + clad_inner)],
                    613.72910879010);
    radial_displacement.add(
        state[dofs.radial_displacement(fuel_offset + fuel_surface)],
        3.2141740355603e-6);
    radial_displacement.add(
        state[dofs.radial_displacement(clad_offset + clad_inner)],
        1.1901972438575e-6);
    axial_displacement.add(
        state[dofs.axial_displacement(fuel_offset + fuel_top)],
        1.0531979661574e-5);

    constexpr std::array<double, 11> expected_pressure = {
        2547100.2219749, 2545172.8286605, 2532649.5016160, 2507500.3516350,
        2465844.2361049, 2413364.3996226, 2413608.4076008, 2557256.2070568,
        2243329.5849293, 1643992.4525044, 6029728.9198301,
    };
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(0, state);
    ErrorMetrics pressure;
    if (contact_nodes.size() == expected_pressure.size()) {
        for (std::size_t node = 0; node < expected_pressure.size(); ++node)
            pressure.add(contact_nodes[node].pressure, expected_pressure[node]);
    }
    passed = check(contact_nodes.size() == expected_pressure.size(),
                   "M1 pressure vectors have matching node counts") &&
             passed;

    constexpr double tolerance = 1.0e-2;
    passed = check_metrics("m1_temperature", temperature, tolerance) && passed;
    passed = check_metrics("m1_radial_displacement", radial_displacement,
                           tolerance) &&
             passed;
    passed =
        check_metrics("m1_axial_displacement", axial_displacement, tolerance) &&
        passed;
    passed =
        check_metrics("m1_contact_pressure", pressure, tolerance) && passed;

    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(0, state);
    ErrorMetrics total_force;
    total_force.add(interface.total_contact_force, 663.8896691615588);
    passed = check_metrics("m1_total_contact_force", total_force, tolerance) &&
             passed;
    passed = check(interface.projected_contact_nodes == 11 &&
                       interface.active_contact_nodes == 11,
                   "M1 projects and activates all fuel-surface nodes") &&
             passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m1_exodus_moose_tests <m1.fsi>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim input-card M1 MOOSE comparison\n");
        if (!run_comparison(argv[1]))
            return 1;
        std::cout << "[PASS] input-card M1 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M1 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
