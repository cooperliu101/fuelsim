#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void write_reference(const std::string& path,
                     const std::vector<double>& state) {
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

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: fuelsim_mpi_equivalence_tests "
                     "<write|compare|compare_field_split> <reference> "
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

        fuelsim::PetscSession session(
            argc, argv, "fuelsim one/two-rank equivalence test\n");
        const fuelsim::FuelSimCaseDefinition definition =
            fuelsim::CaseInputReader::read(input_path);
        const fuelsim::UnstructuredQuad4Mesh source =
            fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
        fuelsim::SteadyProblem problem(definition.steady_definition(), source);
        fuelsim::SolverOptions options;
        options.absolute_tolerance = definition.solver.absolute_tolerance;
        options.relative_tolerance = definition.solver.relative_tolerance;
        options.step_tolerance = definition.solver.step_tolerance;
        options.maximum_iterations = definition.solver.maximum_iterations;
        const bool field_split = mode == "compare_field_split";
        options.linear_solver =
            field_split ? fuelsim::SolverOptions::LinearSolver::gmres
                        : fuelsim::SolverOptions::LinearSolver::direct;
        options.preconditioner =
            field_split
                ? fuelsim::SolverOptions::Preconditioner::field_split
                : fuelsim::SolverOptions::Preconditioner::lu;
        const fuelsim::SteadyResult result = fuelsim::solve_steady(
            problem, definition.steady_execution.load_steps, options);
        if (!result.completed || !result.solve.converged)
            throw std::runtime_error("MPI equivalence solve did not converge");
        if (result.aggregate_timing.workspace_setups != 1)
            throw std::runtime_error(
                "MPI equivalence solve did not reuse one workspace");

        if (mode == "write") {
            if (session.size() != 1)
                throw std::invalid_argument(
                    "MPI reference must be written with one rank");
            write_reference(reference_path, result.solve.state);
            std::cout << "[PASS] wrote one-rank MPI reference\n";
            return 0;
        }
        if (mode != "compare" && !field_split)
            throw std::invalid_argument("Unknown MPI equivalence mode: " + mode);
        if (session.size() != 2)
            throw std::invalid_argument(
                "MPI comparison must run with exactly two ranks");
        const std::vector<double> reference = read_reference(reference_path);
        if (reference.size() != result.solve.state.size())
            throw std::runtime_error("MPI reference state size differs");
        double maximum_absolute = 0.0;
        double maximum_scaled = 0.0;
        for (std::size_t dof = 0; dof < reference.size(); ++dof) {
            const double difference =
                std::abs(result.solve.state[dof] - reference[dof]);
            maximum_absolute = std::max(maximum_absolute, difference);
            maximum_scaled = std::max(
                maximum_scaled, difference / (1.0 + std::abs(reference[dof])));
        }
        const std::size_t expected_begin =
            problem.contribution_count() *
            static_cast<std::size_t>(result.solve.mpi_rank) / 2U;
        const std::size_t expected_end =
            problem.contribution_count() *
            static_cast<std::size_t>(result.solve.mpi_rank + 1) / 2U;
        if (result.solve.local_contribution_begin != expected_begin ||
            result.solve.local_contribution_end != expected_end)
            throw std::runtime_error(
                "MPI contribution partition differs from ownership contract");
        const double tolerance = field_split ? 1.0e-7 : 1.0e-10;
        if (!(maximum_scaled < tolerance))
            throw std::runtime_error(
                "one/two-rank state difference exceeds tolerance");
        if (session.rank() == 0) {
            std::cout << std::scientific << std::setprecision(12)
                      << "mpi_equivalence_maximum_absolute="
                      << maximum_absolute << '\n'
                      << "mpi_equivalence_maximum_scaled=" << maximum_scaled
                      << '\n'
                      << "[PASS] one/two-rank state equivalence"
                      << (field_split ? " with field split\n" : "\n");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] MPI equivalence test raised: " << error.what()
                  << '\n';
        return 1;
    }
}
