#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
using namespace fuelsim;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

// Internal transaction coverage: reuse a complete input solely to obtain a valid
// mesh and model. Inspect solver rollback/history directly; this is not an
// external-solver acceptance case or a replacement for production verification.
std::unique_ptr<TransientProblem> make_problem(const FuelSimCaseDefinition& input) {
    if (input.geometry == CaseGeometry::axisymmetric_1d)
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_bar2(input.mesh_file));
    if (input.geometry == CaseGeometry::generalized_plane_strain)
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_plane_quad8(input.mesh_file));
    if (input.geometry == CaseGeometry::axisymmetric_rz) {
        const auto kind = input.spatial.regions.front().rz_element_formulation;
        if (kind == RzElementFormulation::cax8t || kind == RzElementFormulation::cax8rt)
            return std::make_unique<TransientProblem>(input.spatial, read_exodus_quad8(input.mesh_file));
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_quad4(input.mesh_file));
    }
    if (exodus_uses_hex20(input.mesh_file))
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_hex20(input.mesh_file));
    return std::make_unique<TransientProblem>(input.spatial, read_exodus_hex8(input.mesh_file));
}

void exercise(const char* path) {
    const auto input = read_case_input(path);
    auto storage = make_problem(input);
    auto& problem = *storage;
    auto options = input.transient_execution;
    options.time_error_relative_tolerance = 0;
    options.adaptive_algorithm = fuelsim::AdaptiveTimeAlgorithm::creep_rate;
    options.creep_strain_time_tolerance = 1e-5;
    options.minimum_time_step = 1e-12;
    options.maximum_cutbacks_per_step = 30;
    options.growth_factor = 2;
    const auto initial = problem.capture_state();
    const auto state = BackendAccess::committed_state(problem);
    const auto rates = problem.committed_creep_rates();
    require(!rates.empty(), "Every supported topology exposes active material points");
    problem.begin_time_step({options.initial_time_step, 1, false});
    bool active_rejected = false;
    try {
        (void)problem.committed_creep_rates();
    } catch (const std::logic_error&) {
        active_rejected = true;
    }
    require(active_rejected, "Rate sampling must reject active trial state");
    problem.rollback_time_step();
    problem.restore_state(initial);
    auto impossible = options;
    impossible.creep_strain_time_tolerance = 1e-30;
    impossible.maximum_cutbacks_per_step = 0;
    const auto failed = solve_transient(problem, impossible, input.solver);
    require(!failed.completed && failed.time_error_rejections == 1,
        "Converged oversized step must fail the material-point error criterion");
    const auto restored = BackendAccess::committed_state(problem);
    require(restored.solution == state.solution && restored.previous_solution == state.previous_solution
                && restored.time == state.time && restored.load_factor == state.load_factor
                && restored.previous_time == state.previous_time && restored.raw_residual == state.raw_residual
                && restored.external_load_residual == state.external_load_residual
                && problem.committed_creep_rates() == rates && !problem.time_step_active(),
        "Failed control step must restore nodal fields, rates, clock, loads, residuals and predictor state");
    auto error_options = options;
    error_options.time_error_relative_tolerance = 1e-6;
    require(problem.step_doubling_error(initial, problem.capture_state(), error_options).maximum == 0,
        "Failed control step must restore every material and contact history");
    for (const auto& field : transient_conservation_fields)
        require(restored.conservation.*field.member == state.conservation.*field.member,
            "Failed control step must restore energy and thermal diagnostics");
    const auto result = solve_transient(problem, options, input.solver);
    require(result.completed, "Adaptive controller must finish the supported element case");
    for (const auto& step : result.accepted_steps)
        require(step.time_error_estimate <= 1 && step.time_step <= options.maximum_time_step * (1 + 1e-12),
            "Accepted steps obey rate tolerance and time-step cap");
    const auto final_rates = problem.committed_creep_rates();
    require(std::any_of(final_rates.begin(), final_rates.end(), [](double rate) { return rate > 0; }),
        "Creep must remain active after retries");
    std::cout << path << ": " << result.accepted_steps.size() << " accepted steps, " << result.time_error_rejections
              << " error retries\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Creep-rate adaptive state transactions");
        require(argc > 1, "Expected input cards");
        for (int i = 1; i < argc; ++i)
            exercise(argv[i]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
