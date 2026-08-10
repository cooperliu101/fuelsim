#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void write_reference(const std::string& path, const std::vector<double>& state) {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output)
        throw std::runtime_error("Could not write MPI reference: " + path);
    output << state.size() << '\n' << std::setprecision(17);
    for (double value : state)
        output << value << '\n';
    if (!output)
        throw std::runtime_error("Could not complete MPI reference: " + path);
}

std::vector<double> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MPI reference: " + path);
    std::size_t count = 0;
    input >> count;
    std::vector<double> result(count, 0.0);
    for (double& value : result)
        input >> value;
    if (!input)
        throw std::runtime_error("MPI reference is incomplete: " + path);
    return result;
}

std::vector<double> flatten_transient_state(const fuelsim::TransientProblem& problem) {
    std::vector<double> result = {problem.committed_time(), problem.committed_load_factor()};
    result.insert(result.end(), problem.committed_solution().begin(), problem.committed_solution().end());
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        for (std::size_t element = 0; element < problem.region_mesh(region).elements().size(); ++element) {
            const fuelsim::Quad4MaterialHistory& history = problem.material_history(region, element);
            const auto& stresses = problem.material_stress(region, element);
            for (std::size_t q = 0; q < history.size(); ++q) {
                result.insert(result.end(), history[q].plastic_strain.begin(), history[q].plastic_strain.end());
                result.insert(result.end(), history[q].creep_strain.begin(), history[q].creep_strain.end());
                result.push_back(history[q].equivalent_plastic_strain);
                result.push_back(history[q].equivalent_creep_strain);
                result.push_back(stresses[q].rr);
                result.push_back(stresses[q].zz);
                result.push_back(stresses[q].hoop);
                result.push_back(stresses[q].rz);
            }
        }
    }
    for (const auto& contact : problem.committed_state().contact_histories) {
        for (const fuelsim::ContactPointHistory& history : contact) {
            result.push_back(history.elastic_tangential_slip);
            result.push_back(history.normal_multiplier);
            result.push_back(history.sliding ? 1.0 : 0.0);
        }
    }
    return result;
}

void compare_reference(const std::string& path, const std::vector<double>& state, double tolerance) {
    const std::vector<double> reference = read_reference(path);
    if (reference.size() != state.size())
        throw std::runtime_error("MPI reference state size differs");
    double maximum_absolute = 0.0;
    double maximum_scaled = 0.0;
    std::size_t maximum_scaled_index = 0;
    for (std::size_t value = 0; value < reference.size(); ++value) {
        const double difference = std::abs(state[value] - reference[value]);
        maximum_absolute = std::max(maximum_absolute, difference);
        const double scaled = difference / (1.0 + std::abs(reference[value]));
        if (scaled > maximum_scaled) {
            maximum_scaled = scaled;
            maximum_scaled_index = value;
        }
    }
    if (!(maximum_scaled < tolerance)) {
        std::ostringstream message;
        message << std::scientific << std::setprecision(12) << "one/multi-rank state difference exceeds tolerance: "
                << "maximum absolute=" << maximum_absolute << ", maximum scaled=" << maximum_scaled
                << ", index=" << maximum_scaled_index << ", one-rank=" << reference[maximum_scaled_index]
                << ", two-rank=" << state[maximum_scaled_index];
        throw std::runtime_error(message.str());
    }
    std::cout << std::scientific << std::setprecision(12) << "mpi_equivalence_maximum_absolute=" << maximum_absolute
              << '\n'
              << "mpi_equivalence_maximum_scaled=" << maximum_scaled << '\n';
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_mpi_equivalence_tests "
                     "<write|compare|compare_field_split|compare_block_jacobi|"
                     "compare_hypre|write_transient|compare_transient|"
                     "test_io_failure> "
                     "<reference> "
                     "<case.fsi> [PETSc options]\n";
        return 2;
    }
    try {
        const std::string mode = argv[1];
        const std::string reference_path = argv[2];
        const std::string input_path = argv[3];
        for (int index = 4; index < argc; ++index)
            argv[index - 3] = argv[index];
        argc -= 3;
        argv[argc] = nullptr;

        fuelsim::PetscSession session(argc, argv, "fuelsim one/two-rank equivalence test\n");
        if (mode == "test_io_failure") {
            if (session.size() != 2)
                throw std::invalid_argument("Collective I/O failure test requires two ranks");
            bool caught = false;
            try {
                session.collective_root_action([]() { throw std::runtime_error("intentional root I/O failure"); });
            } catch (const std::runtime_error& error) {
                caught = std::string(error.what()).find("collective root-rank I/O failed") != std::string::npos;
            }
            if (!caught)
                throw std::runtime_error("Collective root I/O failure did not reach every rank");
            session.collective_root_action([]() {});
            if (session.rank() == 0)
                std::cout << "[PASS] root I/O failure reached every rank\n";
            return 0;
        }
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::CaseInputReader::read(input_path);
        const fuelsim::UnstructuredQuad4Mesh source = fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
        fuelsim::SolverOptions options;
        options.absolute_tolerance = definition.solver.absolute_tolerance;
        options.relative_tolerance = definition.solver.relative_tolerance;
        options.step_tolerance = definition.solver.step_tolerance;
        options.maximum_iterations = definition.solver.maximum_iterations;
        options.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
        options.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
        options.backtracking_fallback = definition.solver.backtracking_fallback;
        options.field_residual_scaling = definition.solver.field_residual_scaling;
        options.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
        options.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
        options.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
        options.temperature_residual_scale = definition.solver.temperature_residual_scale;
        options.mechanical_residual_scale = definition.solver.mechanical_residual_scale;
        if (definition.problem == fuelsim::CaseProblem::transient) {
            if (mode != "write_transient" && mode != "compare_transient")
                throw std::invalid_argument("Transient MPI case requires write_transient or "
                                            "compare_transient mode");
            fuelsim::TransientProblem problem(definition.transient_definition(), source);
            const fuelsim::TransientExecutionInput& execution = definition.transient_execution;
            const fuelsim::TransientResult result = fuelsim::solve_transient(
                problem,
                {execution.end_time, execution.initial_time_step, execution.minimum_time_step,
                 execution.maximum_time_step, execution.growth_factor, execution.cutback_factor,
                 execution.maximum_cutbacks, execution.load_ramp_time, execution.target_nonlinear_iterations,
                 execution.iteration_window, execution.time_error_relative_tolerance,
                 execution.temperature_time_absolute_tolerance, execution.displacement_time_absolute_tolerance,
                 execution.time_error_safety_factor, execution.strain_history_time_absolute_tolerance,
                 execution.stress_history_time_absolute_tolerance},
                options);
            if (!result.completed || result.aggregate_timing.workspace_setups != 1) {
                std::ostringstream message;
                message << "Transient MPI equivalence solve failed or rebuilt "
                           "its workspace: completed="
                        << result.completed << ", setups=" << result.aggregate_timing.workspace_setups
                        << ", category=" << fuelsim::solve_failure_category_name(result.last_attempt.failure_category)
                        << ", message=" << result.last_attempt.failure_message;
                throw std::runtime_error(message.str());
            }
            const fuelsim::SolveResult& shadow = result.last_attempt;
            if (shadow.global_state_dofs != problem.dof_count())
                throw std::runtime_error("Transient shadow-state global size is incorrect");
            if (session.size() == 1 && shadow.maximum_shadow_state_dofs != problem.dof_count())
                throw std::runtime_error("One-rank transient solve does not cover the full state");
            if (session.size() > 1 &&
                (shadow.maximum_shadow_state_dofs > problem.dof_count() ||
                 shadow.total_shadow_state_dofs > static_cast<std::size_t>(session.size()) * problem.dof_count() ||
                 shadow.total_remote_shadow_state_dofs == 0))
                throw std::runtime_error(
                    "Multi-rank transient shadow-state bounds failed: global=" + std::to_string(problem.dof_count()) +
                    ", maximum=" + std::to_string(shadow.maximum_shadow_state_dofs) +
                    ", total=" + std::to_string(shadow.total_shadow_state_dofs) +
                    ", remote=" + std::to_string(shadow.total_remote_shadow_state_dofs));
            const std::vector<double> state = flatten_transient_state(problem);
            if (mode == "write_transient") {
                if (session.size() != 1)
                    throw std::invalid_argument("Transient MPI reference requires one rank");
                write_reference(reference_path, state);
                std::cout << "[PASS] wrote transient one-rank MPI reference\n";
                return 0;
            }
            if (session.size() < 2)
                throw std::invalid_argument("Transient MPI comparison requires at least two ranks");
            if (session.rank() == 0) {
                // Distributed assembly changes the summation order. The mixed
                // absolute/relative gate still checks every stored value while
                // allowing sub-micro-Pascal roundoff in nominally zero stress
                // components.
                compare_reference(reference_path, state, 1.0e-7);
                std::cout << "transient_maximum_shadow_state_dofs=" << shadow.maximum_shadow_state_dofs << '\n'
                          << "transient_total_remote_shadow_state_dofs=" << shadow.total_remote_shadow_state_dofs
                          << '\n';
                std::cout << "[PASS] transient nodal, integration-point, and contact-history MPI equivalence\n";
            }
            return 0;
        }

        fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
        const bool field_split = mode == "compare_field_split";
        const bool block_jacobi = mode == "compare_block_jacobi";
        const bool hypre = mode == "compare_hypre";
        options.linear_solver = field_split || block_jacobi || hypre ? fuelsim::SolverOptions::LinearSolver::gmres
                                                                     : fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner = field_split    ? fuelsim::SolverOptions::Preconditioner::field_split
                                 : block_jacobi ? fuelsim::SolverOptions::Preconditioner::block_jacobi
                                 : hypre        ? fuelsim::SolverOptions::Preconditioner::hypre
                                                : fuelsim::SolverOptions::Preconditioner::lu;
        const fuelsim::SteadyResult result = fuelsim::solve_steady(
            problem,
            {field_split ? 2U : definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
             definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
            options);
        if (!result.completed || !result.solve.converged)
            throw std::runtime_error("MPI equivalence solve did not converge");
        if (result.aggregate_timing.workspace_setups != 1)
            throw std::runtime_error("MPI equivalence solve did not reuse one workspace");
        if (result.solve.global_state_dofs != problem.dof_count())
            throw std::runtime_error("Steady shadow-state global size is incorrect");
        if (session.size() == 1 && result.solve.maximum_shadow_state_dofs != problem.dof_count())
            throw std::runtime_error("One-rank steady solve does not cover the full state");

        if (mode == "write") {
            if (session.size() != 1)
                throw std::invalid_argument("MPI reference must be written with one rank");
            write_reference(reference_path, result.solve.state);
            std::cout << "[PASS] wrote one-rank MPI reference\n";
            return 0;
        }
        if (mode != "compare" && !field_split && !block_jacobi && !hypre)
            throw std::invalid_argument("Unknown MPI equivalence mode: " + mode);
        if (session.size() != 2)
            throw std::invalid_argument("MPI comparison must run with exactly two ranks");
        const std::size_t expected_begin =
            problem.contribution_count() * static_cast<std::size_t>(result.solve.mpi_rank) / 2U;
        const std::size_t expected_end =
            problem.contribution_count() * static_cast<std::size_t>(result.solve.mpi_rank + 1) / 2U;
        if (result.solve.local_contribution_begin != expected_begin ||
            result.solve.local_contribution_end != expected_end)
            throw std::runtime_error("MPI contribution partition differs from ownership contract");
        if (!(result.solve.total_shadow_state_dofs < 2 * problem.dof_count()) ||
            result.solve.total_remote_shadow_state_dofs == 0)
            throw std::runtime_error(
                "Two-rank steady shadow-state bounds failed: global=" + std::to_string(problem.dof_count()) +
                ", maximum=" + std::to_string(result.solve.maximum_shadow_state_dofs) +
                ", total=" + std::to_string(result.solve.total_shadow_state_dofs) +
                ", remote=" + std::to_string(result.solve.total_remote_shadow_state_dofs));
        if (session.rank() == 0) {
            compare_reference(reference_path, result.solve.state,
                              field_split || block_jacobi || hypre ? 1.0e-7 : 1.0e-10);
            std::cout << "steady_maximum_shadow_state_dofs=" << result.solve.maximum_shadow_state_dofs << '\n'
                      << "steady_total_remote_shadow_state_dofs=" << result.solve.total_remote_shadow_state_dofs
                      << '\n';
            std::cout << "[PASS] one/two-rank state equivalence"
                      << (field_split    ? " with field split\n"
                          : block_jacobi ? " with block Jacobi\n"
                          : hypre        ? " with hypre\n"
                                         : "\n");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] MPI equivalence test raised: " << error.what() << '\n';
        return 1;
    }
}
