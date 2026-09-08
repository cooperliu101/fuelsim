#include "fuelsim/core/spatial_definition.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/material_factory.hpp"
#include "support/mesh_fixture.hpp"
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
constexpr std::size_t fuel_radial_elements = 100;
constexpr std::size_t cladding_radial_elements = 16;
constexpr std::size_t load_steps = 20;

struct BenchmarkCase final {
    fuelsim::UnstructuredQuad4Mesh mesh;
    fuelsim::SpatialDefinition definition;
};

fuelsim::BoundaryConditionDefinition
dirichlet(const std::string& name, const std::string& boundary, fuelsim::Field field, double value) {
    fuelsim::BoundaryConditionDefinition condition{};
    condition.name = name;
    condition.type = fuelsim::BoundaryConditionType::dirichlet;
    condition.boundary = boundary;
    condition.field = field;
    condition.value = value;
    return condition;
}

BenchmarkCase make_case(std::size_t requested_fuel_radial_elements,
    std::size_t requested_cladding_radial_elements,
    std::size_t axial_elements) {
    BenchmarkCase result{
        fuelsim::test::make_disconnected_annular_mesh(
            {{1, "fuel", 0.0, 0.004120, 0.010, requested_fuel_radial_elements, axial_elements},
                {2, "clad", 0.004122, 0.004692, 0.010020, requested_cladding_radial_elements, axial_elements}}),
        {}};
    result.definition.regions.push_back(
        {"fuel", "fuel", fuelsim::test::thermoelastic(3824.0, 0.61, 2.0e11, 0.316, 10.0e-6, 600.0), 2.0e8, 600.0});
    result.definition.regions.push_back(
        {"clad", "clad", fuelsim::test::thermoelastic(0.0, 16.0, 75.0e9, 0.3, 5.0e-6, 600.0), 0.0, 600.0});
    result.definition.contacts.push_back({"fuel_clad", "clad_inner", "fuel_outer", true, true, 0.4, 1.0e-6, 1.0e14});
    result.definition.boundary_conditions.push_back(
        dirichlet("fuel_axis", "fuel_inner", fuelsim::Field::radial_displacement, 0.0));
    result.definition.boundary_conditions.push_back(
        dirichlet("fuel_bottom", "fuel_bottom", fuelsim::Field::axial_displacement, 0.0));
    result.definition.boundary_conditions.push_back(
        dirichlet("clad_bottom", "clad_bottom", fuelsim::Field::axial_displacement, 0.0));
    result.definition.boundary_conditions.push_back(
        dirichlet("clad_temperature", "clad_outer", fuelsim::Field::temperature, 600.0));
    return result;
}

std::size_t parse_positive_size(const char* text, const char* name) {
    const std::string value(text);
    std::size_t consumed = 0;
    const unsigned long long parsed = std::stoull(value, &consumed);
    if (value.empty() || value.front() == '-' || consumed != value.size() || parsed == 0
        || parsed > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    return static_cast<std::size_t>(parsed);
}

void configure_solver(const std::string& name, fuelsim::SolverOptions& options) {
    if (name == "direct") {
        options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
        return;
    }
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::gmres;
    if (name == "block_jacobi")
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::block_jacobi;
    else if (name == "field_split")
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::field_split;
    else if (name == "hypre")
        options.preconditioner = fuelsim::SolverOptions::Preconditioner::hypre;
    else
        throw std::invalid_argument("solver must be direct, block_jacobi, field_split, or hypre");
}
} // namespace

int main(int argc, char** argv) {
    try {
        const std::string case_name = argc > 1 ? argv[1] : "medium";
        const std::string solver_name = argc > 2 ? argv[2] : "direct";
        const std::size_t requested_load_steps = argc > 3 ? parse_positive_size(argv[3], "load steps") : load_steps;
        const std::string scaling_name = argc > 4 ? argv[4] : "unscaled";
        if (scaling_name != "unscaled" && scaling_name != "scaled")
            throw std::invalid_argument("field scaling must be unscaled or scaled");
        if (argc > 5)
            throw std::invalid_argument("usage: fuelsim_m1_single_core_benchmark "
                                        "[medium|large] [direct|block_jacobi|field_split|hypre] "
                                        "[load_steps] [unscaled|scaled]");
        const std::size_t radial_multiplier = case_name == "medium" ? 1U : case_name == "large" ? 2U : 0U;
        if (radial_multiplier == 0)
            throw std::invalid_argument("case must be medium or large");
        const std::size_t requested_fuel_radial_elements = radial_multiplier * fuel_radial_elements;
        const std::size_t requested_cladding_radial_elements = radial_multiplier * cladding_radial_elements;
        const std::size_t axial_elements = 64U;
        fuelsim::PetscSession session(argc, argv, "fuelsim M1 engineering-scale linear-solver benchmark\n");
        const auto problem_setup_start = std::chrono::steady_clock::now();
        BenchmarkCase benchmark_case =
            make_case(requested_fuel_radial_elements, requested_cladding_radial_elements, axial_elements);
        fuelsim::SteadyProblem problem(std::move(benchmark_case.definition), benchmark_case.mesh);
        const double problem_setup_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - problem_setup_start).count();
        fuelsim::SolverOptions options;
        options.maximum_iterations = 50;
        options.field_residual_scaling = scaling_name == "scaled";
        configure_solver(solver_name, options);
        const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, {requested_load_steps}, options);
        const std::size_t nodes =
            (requested_fuel_radial_elements + 1 + requested_cladding_radial_elements + 1) * (axial_elements + 1);
        const std::size_t elements =
            (requested_fuel_radial_elements + requested_cladding_radial_elements) * axial_elements;
        if (session.rank() != 0)
            return result.completed && result.solve.converged ? 0 : 1;
        std::cout << std::boolalpha << std::scientific << std::setprecision(12);
        std::cout << "case=m1-" << case_name << '\n';
        std::cout << "linear_solver=" << solver_name << '\n';
        std::cout << "field_residual_scaling=" << scaling_name << '\n';
        std::cout << "mpi_ranks=" << session.size() << '\n';
        std::cout << "nodes=" << nodes << '\n';
        std::cout << "elements=" << elements << '\n';
        std::cout << "dofs=" << 3 * nodes << '\n';
        std::cout << "load_steps=" << requested_load_steps << '\n';
        std::cout << "converged=" << (result.completed && result.solve.converged) << '\n';
        std::cout << "load_steps_completed=" << result.completed_steps << '\n';
        std::cout << "nonlinear_iterations_total=" << result.total_nonlinear_iterations << '\n';
        std::cout << "linear_iterations_total=" << result.total_linear_iterations << '\n';
        std::cout << "residual_norm_last_step=" << result.solve.residual_norm << '\n';
        std::cout << "timing_problem_setup=" << problem_setup_seconds << '\n';
        std::cout << "timing_solver_setup=" << result.aggregate_timing.setup_seconds << '\n';
        std::cout << "timing_nonlinear_solve=" << result.aggregate_timing.nonlinear_solve_seconds << '\n';
        std::cout << "timing_residual_callbacks=" << result.aggregate_timing.residual_callback_seconds << '\n';
        std::cout << "timing_jacobian_callbacks=" << result.aggregate_timing.jacobian_callback_seconds << '\n';
        std::cout << "timing_load_path_total=" << result.total_seconds << '\n';
        std::cout << "residual_evaluations=" << result.aggregate_timing.residual_evaluations << '\n';
        std::cout << "jacobian_evaluations=" << result.aggregate_timing.jacobian_evaluations << '\n';
        std::cout << "petsc_workspace_setups=" << result.aggregate_timing.workspace_setups << '\n';
        std::cout << "global_state_dofs=" << result.solve.global_state_dofs << '\n';
        std::cout << "maximum_shadow_state_dofs=" << result.solve.maximum_shadow_state_dofs << '\n';
        std::cout << "total_shadow_state_dofs=" << result.solve.total_shadow_state_dofs << '\n';
        std::cout << "total_remote_shadow_state_dofs=" << result.solve.total_remote_shadow_state_dofs << '\n';
        std::cout << "maximum_shadow_state_bytes=" << result.solve.maximum_shadow_state_dofs * sizeof(double) << '\n';
        std::cout << "maximum_shadow_workspace_bytes="
                  << result.solve.maximum_shadow_state_dofs * (2 * sizeof(double) + sizeof(std::uint32_t))
                         + result.solve.global_state_dofs * sizeof(std::uint32_t)
                  << '\n';
        std::cout << "global_to_shadow_lookup_bytes=" << result.solve.global_state_dofs * sizeof(std::uint32_t) << '\n';
        std::cout << "replicated_callback_state_workspace_bytes=" << result.solve.global_state_dofs * 2 * sizeof(double)
                  << '\n';
        std::cout << "remote_shadow_bytes_per_callback=" << result.solve.total_remote_shadow_state_dofs * sizeof(double)
                  << '\n';
        return result.completed && result.solve.converged ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "fuelsim benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
