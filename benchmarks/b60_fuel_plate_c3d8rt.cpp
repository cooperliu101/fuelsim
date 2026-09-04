#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

void write_nodes(const std::string& path, const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::cartesian::SpatialAssembly& spatial, const std::vector<double>& state) {
    const auto& dofs = spatial;
    std::vector<std::size_t> source_global(mesh.nodes().size(), 0U);
    std::vector<bool> present(mesh.nodes().size(), false);
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids().at(local);
            const std::size_t global = dofs.global_node(region, local);
            if (present[source] && source_global[source] != global)
                throw std::runtime_error("B6.0 interface node maps to multiple global nodes");
            source_global[source] = global;
            present[source] = true;
        }
    }
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B6.0 Fuelsim nodal output: " + path);
    output << "id,x,y,z,temperature,displacement_x,displacement_y,displacement_z\n";
    output << std::scientific << std::setprecision(17);
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        if (!present[source]) throw std::runtime_error("B6.0 mesh node was not assigned a global degree of freedom");
        const std::size_t global = source_global[source];
        const auto& point = mesh.nodes()[source];
        output << (source + 1U) << ',' << point.x << ',' << point.y << ',' << point.z << ','
               << state[dofs.dof(fuelsim::Field::temperature, global)] << ','
               << state[dofs.dof(fuelsim::Field::displacement_x, global)] << ','
               << state[dofs.dof(fuelsim::Field::displacement_y, global)] << ','
               << state[dofs.dof(fuelsim::Field::displacement_z, global)] << '\n';
    }
}

void write_nodes(const std::string& path, const fuelsim::UnstructuredHex20Mesh& mesh,
    const fuelsim::cartesian::SpatialAssembly& spatial, const std::vector<double>& state) {
    std::vector<std::size_t> source_global(mesh.nodes().size(), 0U);
    std::vector<bool> present(mesh.nodes().size(), false);
    std::vector<bool> temperature_active(mesh.nodes().size(), false);
    std::vector<std::size_t> source_temperature_global(mesh.nodes().size(), 0U);
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids().at(local);
            const std::size_t global = spatial.global_node(region, local);
            if (present[source] && source_global[source] != global)
                throw std::runtime_error("B6.0 HEX20 interface node maps to multiple global nodes");
            source_global[source] = global;
            present[source] = true;
            if (region_mesh.temperature_nodes().at(local)) {
                const std::size_t temperature_global = spatial.global_temperature_node(region, local);
                if (temperature_active[source] && source_temperature_global[source] != temperature_global)
                    throw std::runtime_error("B6.0 HEX20 temperature node maps to multiple global nodes");
                source_temperature_global[source] = temperature_global;
                temperature_active[source] = true;
            }
        }
    }
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B6.0 Fuelsim nodal output: " + path);
    output << "id,x,y,z,temperature,displacement_x,displacement_y,displacement_z\n";
    output << std::scientific << std::setprecision(17);
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        if (!present[source])
            throw std::runtime_error("B6.0 HEX20 mesh node was not assigned a global degree of freedom");
        const std::size_t global = source_global[source];
        const auto& point = mesh.nodes()[source];
        const double temperature =
            temperature_active[source]
                ? state[spatial.dof(fuelsim::Field::temperature, source_temperature_global[source])]
                : 0.0;
        output << (source + 1U) << ',' << point.x << ',' << point.y << ',' << point.z << ',' << temperature << ','
               << state[spatial.dof(fuelsim::Field::displacement_x, global)] << ','
               << state[spatial.dof(fuelsim::Field::displacement_y, global)] << ','
               << state[spatial.dof(fuelsim::Field::displacement_z, global)] << '\n';
    }
}

void write_transient_timing(const std::string& path, const fuelsim::TransientResult& result, double setup_seconds,
    double solver_seconds, double total_seconds, const fuelsim::TransientProblem& problem) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B6.0 timing output: " + path);
    output << "metric\tvalue\n"
           << std::setprecision(17) << "setup_seconds\t" << setup_seconds << '\n'
           << "solver_seconds\t" << solver_seconds << '\n'
           << "total_seconds\t" << total_seconds << '\n'
           << "global_state_dofs\t" << problem.dof_count() << '\n'
           << "accepted_steps\t" << result.accepted_steps.size() << '\n'
           << "rejected_steps\t" << result.rejected_steps.size() << '\n'
           << "nonlinear_iterations\t" << result.total_nonlinear_iterations << '\n'
           << "linear_iterations\t" << result.total_linear_iterations << '\n'
           << "residual_evaluations\t" << result.aggregate_timing.residual_evaluations << '\n'
           << "jacobian_evaluations\t" << result.aggregate_timing.jacobian_evaluations << '\n'
           << "residual_callback_seconds\t" << result.aggregate_timing.residual_callback_seconds << '\n'
           << "jacobian_callback_seconds\t" << result.aggregate_timing.jacobian_callback_seconds << '\n'
           << "local_residual_assembly_seconds\t" << result.aggregate_timing.local_residual_assembly_seconds << '\n'
           << "local_jacobian_assembly_seconds\t" << result.aggregate_timing.local_jacobian_assembly_seconds << '\n'
           << "workspace_setups\t" << result.aggregate_timing.workspace_setups << '\n';
    for (std::size_t step = 0; step < result.accepted_steps.size(); ++step) {
        output << "step_" << (step + 1U) << "_nonlinear_iterations\t"
               << result.accepted_steps[step].nonlinear_iterations << '\n'
               << "step_" << (step + 1U) << "_linear_iterations\t" << result.accepted_steps[step].linear_iterations
               << '\n';
    }
}

void write_steady_timing(const std::string& path, const fuelsim::SteadyResult& result, double setup_seconds,
    double solver_seconds, double total_seconds, const fuelsim::SteadyProblem& problem) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B6.0 timing output: " + path);
    output << "metric\tvalue\n"
           << std::setprecision(17) << "setup_seconds\t" << setup_seconds << '\n'
           << "solver_seconds\t" << solver_seconds << '\n'
           << "total_seconds\t" << total_seconds << '\n'
           << "global_state_dofs\t" << problem.dof_count() << '\n'
           << "completed_load_steps\t" << result.completed_steps << '\n'
           << "rejected_load_steps\t" << result.rejected_steps.size() << '\n'
           << "load_cutbacks\t" << result.total_cutbacks << '\n'
           << "used_small_strain_predictor\t" << result.used_small_strain_predictor << '\n'
           << "predictor_nonlinear_iterations\t" << result.predictor_nonlinear_iterations << '\n'
           << "predictor_linear_iterations\t" << result.predictor_linear_iterations << '\n'
           << "predictor_jacobian_evaluations\t" << result.predictor_timing.jacobian_evaluations << '\n'
           << "finite_corrector_nonlinear_iterations\t" << result.solve.nonlinear_iterations << '\n'
           << "finite_corrector_linear_iterations\t" << result.solve.linear_iterations << '\n'
           << "finite_corrector_jacobian_evaluations\t" << result.solve.timing.jacobian_evaluations << '\n'
           << "nonlinear_iterations\t" << result.total_nonlinear_iterations << '\n'
           << "linear_iterations\t" << result.total_linear_iterations << '\n'
           << "residual_evaluations\t" << result.aggregate_timing.residual_evaluations << '\n'
           << "jacobian_evaluations\t" << result.aggregate_timing.jacobian_evaluations << '\n'
           << "residual_callback_seconds\t" << result.aggregate_timing.residual_callback_seconds << '\n'
           << "jacobian_callback_seconds\t" << result.aggregate_timing.jacobian_callback_seconds << '\n'
           << "local_residual_assembly_seconds\t" << result.aggregate_timing.local_residual_assembly_seconds << '\n'
           << "local_jacobian_assembly_seconds\t" << result.aggregate_timing.local_jacobian_assembly_seconds << '\n'
           << "workspace_setups\t" << result.aggregate_timing.workspace_setups << '\n';
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.preconditioner = definition.solver.preconditioner;
    options.direct_factorization = definition.solver.direct_factorization;
    options.mumps_ordering = definition.solver.mumps_ordering;
    options.jacobian_lag = definition.solver.jacobian_lag;
    options.predictor_jacobian_lag = definition.solver.predictor_jacobian_lag;
    options.line_search = definition.solver.line_search;
    options.field_residual_scaling = definition.solver.field_residual_scaling;
    options.field_residual_convergence = definition.solver.field_residual_convergence;
    options.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
    options.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
    options.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
    return options;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_b60_fuel_plate_c3d8rt_benchmark <case.fsi> <nodal.csv> <timing.tsv>\n";
        return 2;
    }
    try {
        fuelsim::PetscSession session(argc, argv, "fuelsim B6.0 C3D8RT fuel-plate bending benchmark\n");
        const auto total_start = Clock::now();
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(argv[1]);
        const auto setup_end = Clock::now();
        const fuelsim::SolverOptions options = solver_options(definition);
        if (fuelsim::exodus_uses_hex20(definition.mesh_file)) {
            const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
            if (definition.problem == fuelsim::CaseProblem::steady) {
                fuelsim::SteadyProblem problem(definition.spatial, mesh);
                const auto problem_end = Clock::now();
                const auto solve_start = Clock::now();
                const fuelsim::SteadyResult result =
                    fuelsim::solve_steady(problem, definition.steady_execution, options);
                const auto solve_end = Clock::now();
                if (!result.completed || !result.solve.converged)
                    throw std::runtime_error(
                        "B6.0 HEX20 steady Fuelsim solve did not complete: " + result.solve.failure_message);
                if (session.rank() == 0) {
                    write_nodes(argv[2], mesh, fuelsim::cartesian::ProblemAccess::view(problem), result.solve.state);
                    const double setup_seconds = std::chrono::duration<double>(problem_end - setup_end).count();
                    const double solver_seconds = std::chrono::duration<double>(solve_end - solve_start).count();
                    const double total_seconds = std::chrono::duration<double>(solve_end - total_start).count();
                    write_steady_timing(argv[3], result, setup_seconds, solver_seconds, total_seconds, problem);
                    std::cout << std::scientific << std::setprecision(12) << "b60_dofs=" << problem.dof_count() << '\n'
                              << "b60_completed_load_steps=" << result.completed_steps << '\n'
                              << "b60_setup_seconds=" << setup_seconds << '\n'
                              << "b60_solver_seconds=" << solver_seconds << '\n'
                              << "b60_total_seconds=" << total_seconds << '\n'
                              << "b60_workspace_setups=" << result.aggregate_timing.workspace_setups << '\n';
                }
                return 0;
            }
            fuelsim::TransientProblem problem(definition.spatial, mesh);
            const auto problem_end = Clock::now();
            const auto& execution = definition.transient_execution;
            const fuelsim::TransientTimeOptions time_options = {execution.end_time, execution.initial_time_step,
                execution.minimum_time_step, execution.maximum_time_step, execution.growth_factor,
                execution.cutback_factor, execution.maximum_cutbacks_per_step, execution.load_ramp_time,
                execution.target_nonlinear_iterations, execution.iteration_window,
                execution.time_error_relative_tolerance, execution.temperature_time_absolute_tolerance,
                execution.displacement_time_absolute_tolerance, execution.time_error_safety_factor,
                execution.strain_history_time_absolute_tolerance, execution.stress_history_time_absolute_tolerance,
                execution.include_thermal_time_term, execution.use_linear_time_predictor};
            const auto solve_start = Clock::now();
            const fuelsim::TransientResult result = fuelsim::solve_transient(problem, time_options, options);
            const auto solve_end = Clock::now();
            if (!result.completed)
                throw std::runtime_error(
                    "B6.0 HEX20 Fuelsim solve did not complete: " + result.last_attempt.failure_message);
            if (session.rank() == 0) {
                write_nodes(
                    argv[2], mesh, fuelsim::cartesian::ProblemAccess::view(problem), problem.committed_solution());
                const double setup_seconds = std::chrono::duration<double>(problem_end - setup_end).count();
                const double solver_seconds = std::chrono::duration<double>(solve_end - solve_start).count();
                const double total_seconds = std::chrono::duration<double>(solve_end - total_start).count();
                write_transient_timing(argv[3], result, setup_seconds, solver_seconds, total_seconds, problem);
                std::cout << std::scientific << std::setprecision(12) << "b60_dofs=" << problem.dof_count() << '\n'
                          << "b60_accepted_steps=" << result.accepted_steps.size() << '\n'
                          << "b60_setup_seconds=" << setup_seconds << '\n'
                          << "b60_solver_seconds=" << solver_seconds << '\n'
                          << "b60_total_seconds=" << total_seconds << '\n'
                          << "b60_workspace_setups=" << result.aggregate_timing.workspace_setups << '\n';
            }
            return 0;
        }
        const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
        if (definition.problem == fuelsim::CaseProblem::steady) {
            fuelsim::SteadyProblem problem(definition.spatial, mesh);
            const auto problem_end = Clock::now();
            const auto solve_start = Clock::now();
            const fuelsim::SteadyResult result = fuelsim::solve_steady(problem, definition.steady_execution, options);
            const auto solve_end = Clock::now();
            if (!result.completed || !result.solve.converged)
                throw std::runtime_error("B6.0 steady Fuelsim solve did not complete: " + result.solve.failure_message);
            if (session.rank() == 0) {
                write_nodes(argv[2], mesh, fuelsim::cartesian::ProblemAccess::view(problem), result.solve.state);
                const double setup_seconds = std::chrono::duration<double>(problem_end - setup_end).count();
                const double solver_seconds = std::chrono::duration<double>(solve_end - solve_start).count();
                const double total_seconds = std::chrono::duration<double>(solve_end - total_start).count();
                write_steady_timing(argv[3], result, setup_seconds, solver_seconds, total_seconds, problem);
                std::cout << std::scientific << std::setprecision(12) << "b60_dofs=" << problem.dof_count() << '\n'
                          << "b60_completed_load_steps=" << result.completed_steps << '\n'
                          << "b60_setup_seconds=" << setup_seconds << '\n'
                          << "b60_solver_seconds=" << solver_seconds << '\n'
                          << "b60_total_seconds=" << total_seconds << '\n'
                          << "b60_workspace_setups=" << result.aggregate_timing.workspace_setups << '\n';
            }
            return 0;
        }
        fuelsim::TransientProblem problem(definition.spatial, mesh);
        const auto problem_end = Clock::now();
        const auto& execution = definition.transient_execution;
        const fuelsim::TransientTimeOptions time_options = {execution.end_time, execution.initial_time_step,
            execution.minimum_time_step, execution.maximum_time_step, execution.growth_factor, execution.cutback_factor,
            execution.maximum_cutbacks_per_step, execution.load_ramp_time, execution.target_nonlinear_iterations,
            execution.iteration_window, execution.time_error_relative_tolerance,
            execution.temperature_time_absolute_tolerance, execution.displacement_time_absolute_tolerance,
            execution.time_error_safety_factor, execution.strain_history_time_absolute_tolerance,
            execution.stress_history_time_absolute_tolerance, execution.include_thermal_time_term,
            execution.use_linear_time_predictor};
        const auto solve_start = Clock::now();
        const fuelsim::TransientResult result = fuelsim::solve_transient(problem, time_options, options);
        const auto solve_end = Clock::now();
        if (!result.completed)
            throw std::runtime_error("B6.0 Fuelsim solve did not complete: " + result.last_attempt.failure_message);
        if (session.rank() == 0) {
            write_nodes(argv[2], mesh, fuelsim::cartesian::ProblemAccess::view(problem), problem.committed_solution());
            const double setup_seconds = std::chrono::duration<double>(problem_end - setup_end).count();
            const double solver_seconds = std::chrono::duration<double>(solve_end - solve_start).count();
            const double total_seconds = std::chrono::duration<double>(solve_end - total_start).count();
            write_transient_timing(argv[3], result, setup_seconds, solver_seconds, total_seconds, problem);
            std::cout << std::scientific << std::setprecision(12) << "b60_dofs=" << problem.dof_count() << '\n'
                      << "b60_accepted_steps=" << result.accepted_steps.size() << '\n'
                      << "b60_setup_seconds=" << setup_seconds << '\n'
                      << "b60_solver_seconds=" << solver_seconds << '\n'
                      << "b60_total_seconds=" << total_seconds << '\n'
                      << "b60_workspace_setups=" << result.aggregate_timing.workspace_setups << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B6.0 benchmark: " << error.what() << '\n';
        return 1;
    }
}
