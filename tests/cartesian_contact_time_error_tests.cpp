#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using namespace fuelsim;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void advance(TransientProblem& problem, PetscSolver& solver, const SolverOptions& options, double end) {
    problem.begin_time_step({end, 1.0});
    const auto solved = solver.solve(problem, problem.committed_solution(), options);
    require(solved.converged, "Contact history refinement solve did not converge");
    problem.commit_time_step(solved.state);
}

void check(const FuelSimCaseDefinition& input, const std::string& output) {
    TransientProblem problem(input.spatial, read_exodus_hex8(input.mesh_file));
    PetscSolver solver;
    for (int step = 1; step <= 10; ++step)
        advance(problem, solver, input.solver, step * 0.02);
    const auto baseline = problem.capture_state();
    const auto initial = BackendAccess::committed_state(problem);
    require(!initial.contact_histories.empty() && !initial.contact_histories[0].empty(), "No real contact history");
    auto options = input.transient_execution;
    options.time_error_relative_tolerance = 1e-3;
    options.displacement_time_absolute_tolerance = 1e-8;
    options.strain_history_time_absolute_tolerance = 1e-8;
    options.stress_history_time_absolute_tolerance = 1.0;
    require(problem.step_doubling_error(baseline, baseline, options).maximum == 0.0, "Identical snapshots differ");
    auto altered = initial;
    altered.contact_histories[0][0].cartesian_elastic_tangential_slip[1] += 1e-5;
    BackendAccess::restore_committed_state(problem, altered);
    require(problem.step_doubling_error(baseline, problem.capture_state(), options).contact_friction > 0.0,
        "Elastic slip missing from time error");
    altered = initial;
    altered.contact_histories[0][0].cartesian_total_tangential_slip[1] += 1e-5;
    BackendAccess::restore_committed_state(problem, altered);
    require(problem.step_doubling_error(baseline, problem.capture_state(), options).contact_friction > 0.0,
        "Accumulated slip missing from time error");
    altered = initial;
    altered.contact_histories[0][0].normal_multiplier += 100.0;
    BackendAccess::restore_committed_state(problem, altered);
    require(problem.step_doubling_error(baseline, problem.capture_state(), options).contact_normal_multiplier > 0.0,
        "Normal multiplier missing from time error");
    altered = initial;
    altered.contact_histories[0][0].sliding = !altered.contact_histories[0][0].sliding;
    BackendAccess::restore_committed_state(problem, altered);
    require(std::isinf(problem.step_doubling_error(baseline, problem.capture_state(), options).maximum),
        "Sliding state mismatch must reject the time step");
    altered = initial;
    auto& contact = altered.contact_histories[0][0];
    require(contact.cartesian_tangent_basis_initialized, "Expected initialized finite-sliding basis");
    for (double& value : contact.cartesian_contact_tangent_first)
        value = -value;
    BackendAccess::restore_committed_state(problem, altered);
    require(problem.step_doubling_error(baseline, problem.capture_state(), options).maximum == 0.0,
        "Equivalent tangent sign caused a false history error");
    const auto n = contact.cartesian_contact_normal, t = contact.cartesian_contact_tangent_first;
    contact.cartesian_contact_tangent_first = {n[1] * t[2] - n[2] * t[1],
        n[2] * t[0] - n[0] * t[2],
        n[0] * t[1] - n[1] * t[0]};
    BackendAccess::restore_committed_state(problem, altered);
    require(problem.step_doubling_error(baseline, problem.capture_state(), options).contact_friction > 0.0,
        "Tangent orientation change was ignored");
    altered = initial;
    for (double& value : altered.contact_histories[0][0].cartesian_contact_normal)
        value = -value;
    BackendAccess::restore_committed_state(problem, altered);
    require(problem.step_doubling_error(baseline, problem.capture_state(), options).contact_friction > 0.0,
        "Oriented contact normal change was ignored");
    altered = initial;
    altered.contact_histories[0][0].cartesian_tangent_basis_initialized = false;
    BackendAccess::restore_committed_state(problem, altered);
    require(std::isinf(problem.step_doubling_error(baseline, problem.capture_state(), options).maximum),
        "Initialization state mismatch must reject the time step");
    std::ofstream table(output);
    table << "dt,contact_friction,contact_normal_multiplier,maximum\n" << std::setprecision(17);
    double previous = 0.0;
    for (double dt : {0.04, 0.02, 0.01}) {
        problem.restore_state(baseline);
        advance(problem, solver, input.solver, 0.2 + dt);
        const auto full = problem.capture_state();
        problem.restore_state(baseline);
        advance(problem, solver, input.solver, 0.2 + dt / 2.0);
        advance(problem, solver, input.solver, 0.2 + dt);
        const auto error = problem.step_doubling_error(full, problem.capture_state(), options);
        require(std::isfinite(error.contact_friction) && error.contact_friction > 0.0,
            "Real contact refinement must detect a finite positive history error");
        if (dt == 0.04)
            require(error.contact_friction > 1.0, "Coarse step must exceed the contact acceptance threshold");
        if (dt == 0.01)
            require(error.contact_friction < 1.0, "Refinement must resolve the contact history error");
        if (previous > 0.0)
            require(error.contact_friction < previous, "Contact history error did not decrease with time refinement");
        previous = error.contact_friction;
        table << dt << ',' << error.contact_friction << ',' << error.contact_normal_multiplier << ',' << error.maximum
              << '\n';
        std::cout << std::setprecision(17) << "dt=" << dt << " contact_friction=" << error.contact_friction
                  << " contact_normal_multiplier=" << error.contact_normal_multiplier << " maximum=" << error.maximum
                  << '\n';
    }
    require(static_cast<bool>(table), "Cannot write contact refinement evidence");
    std::cout << "history_contracts=passed\n";
}
} // namespace

int main(int argc, char** argv) {
    fuelsim::PetscSession session(argc, argv, "Cartesian contact time error contracts");
    try {
        require(argc == 3 && session.size() == 1, "Expected input and refinement CSV, one rank");
        check(fuelsim::read_case_input(argv[1]), argv[2]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
