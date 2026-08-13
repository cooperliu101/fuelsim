#include "fuelsim/problem_solver.hpp"
#include "support/material_factory.hpp"
#include "support/mesh_fixture.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double density = 2.0;
constexpr double specific_heat = 3.0;
constexpr double conductivity = 5.0;
constexpr double end_time = 1.0;
struct ManufacturedParameters final {
    double base_temperature;
    double radial_quadratic;
    double linear_time;
    double quadratic_time;
};
struct ErrorMetrics final {
    double absolute_l2 = 0.0;
    double relative_l2 = 0.0;
    double maximum_absolute = 0.0;
};
struct ManufacturedResult final {
    ErrorMetrics temperature;
    double maximum_displacement = 0.0;
    std::size_t accepted_steps = 0;
    std::size_t workspace_setups = 0;
};
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
double exact_temperature(const ManufacturedParameters& parameters, double r, double time) {
    return parameters.base_temperature + parameters.radial_quadratic * r * r + parameters.linear_time * time +
           parameters.quadratic_time * time * time;
}
double exact_heat_source(const ManufacturedParameters& parameters, double time) {
    return density * specific_heat * (parameters.linear_time + 2.0 * parameters.quadratic_time * time) -
           4.0 * conductivity * parameters.radial_quadratic;
}
fuelsim::UnstructuredQuad4Mesh make_mesh(std::size_t radial_elements, std::size_t axial_elements) {
    return fuelsim::test::make_disconnected_annular_mesh(
        {{1, "solid", 0.0, 1.0, 1.0, radial_elements, axial_elements}});
}
fuelsim::BoundaryConditionDefinition dirichlet(const std::string& name, const std::string& boundary,
    fuelsim::Field field, double value, const std::string& function = {}) {
    fuelsim::BoundaryConditionDefinition result{};
    result.name = name;
    result.type = fuelsim::BoundaryConditionType::dirichlet;
    result.boundary = boundary;
    result.field = field;
    result.value = value;
    result.function = function;
    return result;
}
fuelsim::SpatialDefinition make_definition(const ManufacturedParameters& parameters, std::size_t time_steps) {
    std::vector<double> times;
    std::vector<double> outer_temperatures;
    std::vector<double> heat_sources;
    times.reserve(time_steps + 1);
    outer_temperatures.reserve(time_steps + 1);
    heat_sources.reserve(time_steps + 1);
    for (std::size_t step = 0; step <= time_steps; ++step) {
        const double time = end_time * static_cast<double>(step) / static_cast<double>(time_steps);
        times.push_back(time);
        outer_temperatures.push_back(exact_temperature(parameters, 1.0, time));
        heat_sources.push_back(exact_heat_source(parameters, time));
    }
    fuelsim::SpatialDefinition spatial;
    spatial.regions.push_back({"solid", "solid",
        fuelsim::test::thermoelastic(0.0, conductivity, 1.0e6, 0.3, 0.0, 300.0, 0.0, 0.0, 0.0, density, specific_heat),
        1.0, parameters.base_temperature, -1, "source"});
    spatial.boundary_conditions.push_back(
        dirichlet("axis_radial", "solid_inner", fuelsim::Field::radial_displacement, 0.0));
    spatial.boundary_conditions.push_back(
        dirichlet("bottom_axial", "solid_bottom", fuelsim::Field::axial_displacement, 0.0));
    spatial.boundary_conditions.push_back(
        dirichlet("outer_temperature", "solid_outer", fuelsim::Field::temperature, 1.0, "outer_temperature"));
    spatial.time_tables.emplace_back("outer_temperature", times, outer_temperatures);
    spatial.time_tables.emplace_back("source", times, heat_sources);
    return spatial;
}
void set_exact_initial_state(fuelsim::TransientProblem& problem, const ManufacturedParameters& parameters) {
    fuelsim::TransientCommittedState state = fuelsim::rz::ProblemAccess::committed_state(problem);
    const fuelsim::RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, 0);
    const std::size_t offset = fuelsim::rz::ProblemAccess::region_node_offset(problem, 0);
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const std::size_t dof =
            fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, offset + node);
        state.solution[dof] = exact_temperature(parameters, mesh.nodes()[node].r, 0.0);
    }
    fuelsim::rz::ProblemAccess::restore_committed_state(problem, std::move(state));
}
ErrorMetrics temperature_error(const fuelsim::TransientProblem& problem, const ManufacturedParameters& parameters) {
    const std::array<double, 3> points = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
    const std::array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    const fuelsim::RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, 0);
    const std::size_t offset = fuelsim::rz::ProblemAccess::region_node_offset(problem, 0);
    const std::vector<double>& state = problem.committed_solution();
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    ErrorMetrics result;
    for (const fuelsim::Quad4Element& element : mesh.elements()) {
        for (std::size_t i = 0; i < points.size(); ++i) {
            const double xi = points[i];
            for (std::size_t j = 0; j < points.size(); ++j) {
                const double eta = points[j];
                const std::array<double, 4> shape = {0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta),
                    0.25 * (1.0 + xi) * (1.0 + eta), 0.25 * (1.0 - xi) * (1.0 + eta)};
                const std::array<double, 4> dshape_dxi = {
                    -0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)};
                const std::array<double, 4> dshape_deta = {
                    -0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)};
                double r = 0.0;
                double actual = 0.0;
                double dr_dxi = 0.0;
                double dr_deta = 0.0;
                double dz_dxi = 0.0;
                double dz_deta = 0.0;
                for (std::size_t node = 0; node < shape.size(); ++node) {
                    const fuelsim::RzPoint& point = mesh.nodes().at(element.nodes[node]);
                    r += shape[node] * point.r;
                    actual += shape[node] * state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(
                                                fuelsim::Field::temperature, offset + element.nodes[node])];
                    dr_dxi += dshape_dxi[node] * point.r;
                    dr_deta += dshape_deta[node] * point.r;
                    dz_dxi += dshape_dxi[node] * point.z;
                    dz_deta += dshape_deta[node] * point.z;
                }
                const double determinant = dr_dxi * dz_deta - dr_deta * dz_dxi;
                if (!(determinant > 0.0) || !(r >= 0.0))
                    throw std::runtime_error("manufactured heat error quadrature has invalid "
                                             "geometry");
                const double exact = exact_temperature(parameters, r, end_time);
                const double difference = actual - exact;
                const double measure = 2.0 * pi * r * determinant * weights[i] * weights[j];
                difference_squared += difference * difference * measure;
                reference_squared += exact * exact * measure;
                result.maximum_absolute = std::max(result.maximum_absolute, std::abs(difference));
            }
        }
    }
    result.absolute_l2 = std::sqrt(difference_squared);
    result.relative_l2 = std::sqrt(difference_squared / reference_squared);
    return result;
}
ManufacturedResult solve_manufactured(std::size_t radial_elements, std::size_t axial_elements, std::size_t time_steps,
    const ManufacturedParameters& parameters) {
    const fuelsim::UnstructuredQuad4Mesh mesh = make_mesh(radial_elements, axial_elements);
    fuelsim::TransientProblem problem(make_definition(parameters, time_steps), mesh);
    set_exact_initial_state(problem, parameters);
    const double time_step = end_time / static_cast<double>(time_steps);
    fuelsim::SolverOptions solver_options;
    solver_options.absolute_tolerance = 1.0e-11;
    solver_options.relative_tolerance = 1.0e-12;
    solver_options.step_tolerance = 1.0e-14;
    solver_options.linear_relative_tolerance = 1.0e-12;
    solver_options.maximum_iterations = 20;
    solver_options.backtracking_fallback = false;
    const fuelsim::TransientTimeOptions time_options = {
        end_time, time_step, time_step, time_step, 1.0, 0.5, 0, end_time};
    const fuelsim::TransientResult solve = fuelsim::solve_transient(problem, time_options, solver_options);
    if (!solve.completed) throw std::runtime_error("manufactured transient heat solve did not complete");
    ManufacturedResult result;
    result.temperature = temperature_error(problem, parameters);
    result.accepted_steps = solve.accepted_steps.size();
    result.workspace_setups = solve.aggregate_timing.workspace_setups;
    const std::size_t node_count = fuelsim::rz::ProblemAccess::dof_map(problem).node_count();
    for (std::size_t node = 0; node < node_count; ++node) {
        result.maximum_displacement = std::max(result.maximum_displacement,
            std::abs(problem.committed_solution()[fuelsim::rz::ProblemAccess::dof_map(problem).dof(
                fuelsim::Field::radial_displacement, node)]));
        result.maximum_displacement = std::max(result.maximum_displacement,
            std::abs(problem.committed_solution()[fuelsim::rz::ProblemAccess::dof_map(problem).dof(
                fuelsim::Field::axial_displacement, node)]));
    }
    return result;
}
double observed_order(double coarse, double fine) { return std::log(coarse / fine) / std::log(2.0); }
bool test_spatial_order() {
    const ManufacturedParameters parameters = {300.0, -50.0, 10.0, 0.0};
    const std::array<std::size_t, 4> elements = {4, 8, 16, 32};
    std::array<ManufacturedResult, 4> results{};
    for (std::size_t level = 0; level < elements.size(); ++level)
        results[level] = solve_manufactured(elements[level], 2, 1, parameters);
    std::array<double, 3> orders{};
    bool passed = true;
    for (std::size_t level = 0; level < orders.size(); ++level) {
        orders[level] =
            observed_order(results[level].temperature.absolute_l2, results[level + 1].temperature.absolute_l2);
        passed = check(orders[level] > 1.9 && orders[level] < 2.1, "manufactured transient heat spatial L2 order is "
                                                                   "second order") &&
                 passed;
    }
    for (const ManufacturedResult& result : results) {
        passed =
            check(result.accepted_steps == 1 && result.workspace_setups == 1 && result.maximum_displacement < 1.0e-13,
                "spatial manufactured solve uses one step, one "
                "workspace, and zero mechanics") &&
            passed;
    }
    std::cout << "m21_manufactured_spatial_reference_variation_K=50\n"
              << "m21_manufactured_spatial_relative_l2=";
    for (const ManufacturedResult& result : results) std::cout << result.temperature.relative_l2 << ',';
    std::cout << '\n' << "m21_manufactured_spatial_maximum_absolute_K=";
    for (const ManufacturedResult& result : results) std::cout << result.temperature.maximum_absolute << ',';
    std::cout << '\n'
              << "m21_manufactured_spatial_orders=" << orders[0] << ',' << orders[1] << ',' << orders[2] << '\n';
    return passed;
}
bool test_temporal_order() {
    const ManufacturedParameters parameters = {300.0, -1.0, 0.0, 10.0};
    const std::array<std::size_t, 4> steps = {5, 10, 20, 40};
    std::array<ManufacturedResult, 4> results{};
    for (std::size_t level = 0; level < steps.size(); ++level)
        results[level] = solve_manufactured(64, 2, steps[level], parameters);
    std::array<double, 3> orders{};
    bool passed = true;
    for (std::size_t level = 0; level < orders.size(); ++level) {
        orders[level] =
            observed_order(results[level].temperature.absolute_l2, results[level + 1].temperature.absolute_l2);
        passed = check(orders[level] > 0.85 && orders[level] < 1.15, "manufactured transient heat temporal L2 order is "
                                                                     "first order") &&
                 passed;
    }
    for (std::size_t level = 0; level < results.size(); ++level) {
        passed = check(results[level].accepted_steps == steps[level] && results[level].workspace_setups == 1 &&
                           results[level].maximum_displacement < 1.0e-13,
                     "temporal manufactured solve uses every fixed step, "
                     "one workspace, and zero mechanics") &&
                 passed;
    }
    std::cout << "m21_manufactured_temporal_reference_change_K=10\n"
              << "m21_manufactured_temporal_relative_l2=";
    for (const ManufacturedResult& result : results) std::cout << result.temperature.relative_l2 << ',';
    std::cout << '\n' << "m21_manufactured_temporal_maximum_absolute_K=";
    for (const ManufacturedResult& result : results) std::cout << result.temperature.maximum_absolute << ',';
    std::cout << '\n'
              << "m21_manufactured_temporal_orders=" << orders[0] << ',' << orders[1] << ',' << orders[2] << '\n';
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M2.1 nonuniform manufactured heat convergence test\n");
        bool passed = test_spatial_order();
        passed = test_temporal_order() && passed;
        if (!passed) return 1;
        std::cout << "[PASS] nonuniform manufactured transient heat spatial "
                     "and temporal orders\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] manufactured heat test raised: " << error.what() << '\n';
        return 1;
    }
}
