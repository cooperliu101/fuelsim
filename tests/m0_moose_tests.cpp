#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
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
    throw std::invalid_argument("M0 comparison point is not in the mesh");
}

bool run_comparison(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument(
            "M0 comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    const fuelsim::SolverOptions options = {
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance, definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem, definition.steady_execution.load_steps, options);
    bool passed = check(result.completed && result.solve.converged,
                        "M0 input-card solve converged");
    passed = check(problem.region_count() == 1 &&
                       problem.region_mesh(0).elements().size() == 400,
                   "M0 input card selects the complete MOOSE block") &&
             passed;

    const fuelsim::RegionMesh& mesh = problem.region_mesh(0);
    const fuelsim::DofMap& dofs = problem.dof_map();
    const std::vector<double>& state = result.solve.state;
    const std::size_t axis_mid = find_node(mesh, 0.0, 0.005);
    const std::size_t outer_mid = find_node(mesh, 0.00412, 0.005);
    const std::size_t outer_top = find_node(mesh, 0.00412, 0.010);
    const std::size_t axis_top = find_node(mesh, 0.0, 0.010);

    ErrorMetrics temperature;
    ErrorMetrics radial_displacement;
    ErrorMetrics axial_displacement;
    temperature.add(state[dofs.temperature(axis_mid)], 733.4201407826);
    radial_displacement.add(state[dofs.radial_displacement(outer_mid)],
                            2.5827190582385e-6);
    radial_displacement.add(state[dofs.radial_displacement(outer_top)],
                            3.5772306219751e-6);
    axial_displacement.add(state[dofs.axial_displacement(axis_top)],
                           8.7881251481059e-6);
    constexpr double tolerance = 1.0e-10;
    passed = check_metrics("m0_temperature", temperature, tolerance) && passed;
    passed = check_metrics("m0_radial_displacement", radial_displacement,
                           tolerance) &&
             passed;
    passed =
        check_metrics("m0_axial_displacement", axial_displacement, tolerance) &&
        passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_m0_moose_tests <m0.fsi>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim input-card M0 MOOSE comparison\n");
        if (!run_comparison(argv[1]))
            return 1;
        std::cout << "[PASS] input-card M0 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M0 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
