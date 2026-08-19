#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::TransientTimeOptions time_options(double end_time) { return {end_time, 0.05, 0.05, 0.05, 1.0, 0.5, 2, 1.0}; }

bool run_test(const std::string& input_path, const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh = fuelsim::read_exodus_quad4(input.mesh_file);
    const fuelsim::SolverOptions solver{input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    fuelsim::TransientProblem uninterrupted(input.spatial, mesh);
    const fuelsim::TransientResult full = fuelsim::solve_transient(uninterrupted, time_options(1.0), solver);
    bool passed = check(full.completed, "M5.2 uninterrupted dynamic-search path completes");
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first = fuelsim::solve_transient(split, time_options(0.5), solver);
    passed = check(first.completed, "M5.2 pre-checkpoint dynamic-search path completes") && passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, split, 0.05);
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    const fuelsim::TransientResult second = fuelsim::solve_transient(restarted, time_options(1.0), solver);
    passed = check(restored_step == 0.05 && second.completed, "M5.2 checkpoint restores its controller step and "
                                                              "continues") &&
             passed;
    const std::vector<double>& expected = uninterrupted.committed_solution();
    const std::vector<double>& actual = restarted.committed_solution();
    double maximum_absolute = 0.0;
    double maximum_scaled = 0.0;
    for (std::size_t dof = 0; dof < expected.size(); ++dof) {
        const double difference = std::abs(actual[dof] - expected[dof]);
        maximum_absolute = std::max(maximum_absolute, difference);
        maximum_scaled = std::max(maximum_scaled, difference / (1.0 + std::abs(expected[dof])));
    }
    const std::vector<fuelsim::ContactNodeSummary> expected_contact =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(uninterrupted, 0, expected);
    const std::vector<fuelsim::ContactNodeSummary> actual_contact =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(restarted, 0, actual);
    bool contact_equal = expected_contact.size() == actual_contact.size();
    for (std::size_t node = 0; contact_equal && node < expected_contact.size(); ++node) {
        contact_equal = expected_contact[node].projected == actual_contact[node].projected &&
                        expected_contact[node].primary_segment == actual_contact[node].primary_segment &&
                        expected_contact[node].pressure == actual_contact[node].pressure;
    }
    passed =
        check(maximum_scaled < 1.0e-13 && contact_equal, "M5.2 restarted nodal state, active primary segments, and "
                                                         "pressures reproduce the uninterrupted path") &&
        passed;
    std::cout << "m52_restart_maximum_absolute=" << maximum_absolute << '\n'
              << "m52_restart_maximum_scaled=" << maximum_scaled << '\n';
    return check(std::remove(checkpoint_path.c_str()) == 0, "M5.2 checkpoint artifact is removed") && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_m52_large_sliding_restart_tests "
                     "<case.fsi> <checkpoint>\n";
        return 2;
    }
    try {
        fuelsim::PetscSession session(argc, argv, "fuelsim M5.2 large-sliding restart test\n");
        if (!run_test(argv[1], argv[2])) return 1;
        std::cout << "[PASS] M5.2 large-sliding checkpoint restart\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.2 restart raised: " << error.what() << '\n';
        return 1;
    }
}
