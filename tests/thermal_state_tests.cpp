#include "core/problem_backend_access.hpp"
#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::runtime_error("Expected thermal card and checkpoint path");
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_quad8(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        if (problem.field_layout().size() != 1 || problem.dof_count() != mesh.nodes().size())
            throw std::runtime_error("Thermal model must contain exactly all node temperatures");
        const auto original = problem.committed_solution();
        const auto snapshot = problem.capture_state();
        problem.begin_time_step({2, 1, true});
        auto perturbed = original;
        for (double& value : perturbed)
            value += 1;
        problem.commit_time_step(perturbed);
        const auto accepted = problem.capture_state();
        problem.restore_state(snapshot);
        if (problem.committed_solution() != original || problem.committed_time() != 0)
            throw std::runtime_error("Thermal snapshot did not restore temperature and time");
        problem.begin_time_step({1, 0.5, true});
        problem.rollback_time_step();
        if (problem.time_step_active() || problem.committed_solution() != original)
            throw std::runtime_error("Thermal rollback changed committed temperature");
        const auto state = fuelsim::BackendAccess::committed_state(problem);
        if (!state.material_histories.empty() || !state.quad8_material_histories.empty()
            || !state.cartesian_material_histories.empty() || !state.contact_histories.empty())
            throw std::runtime_error("Thermal state allocated mechanical history");
        fuelsim::write_transient_checkpoint(argv[2], problem, 2.0);
        problem.restore_state(accepted);
        const double next = fuelsim::restore_transient_checkpoint(argv[2], problem);
        if (next != 2.0 || problem.committed_solution() != original || problem.committed_time() != 0.0)
            throw std::runtime_error("Thermal checkpoint did not restore accepted temperature and time");
        const auto error = problem.step_doubling_error(snapshot, accepted, input.transient_execution);
        if (!std::isfinite(error.maximum) || error.maximum <= 0)
            throw std::runtime_error("Thermal time-error estimator ignored temperature");
        std::cout << "thermal temperature layout, snapshot, rollback, and time error passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
