#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/exodus_result_reader.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
std::string checkpoint_bytes(const std::string& path, const fuelsim::TransientProblem& problem) {
    fuelsim::write_transient_checkpoint(path, problem, 0.1);
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot read friction checkpoint bytes");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    try {
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_quad4(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        const double step = fuelsim::restore_transient_checkpoint(argv[2], problem);
        std::cout << std::hexfloat << "restored_step=" << step << " restored_time=" << problem.committed_time() << '\n';
        // Landing on the final amplitude knot can shorten 0.1 by a few representable intervals.
        if (std::abs(step - 0.1) > 4 * std::numeric_limits<double>::epsilon() * 0.1 ||
            std::abs(problem.committed_time() - 1.0) > 2e-15)
            throw std::runtime_error("B8.1 production checkpoint did not preserve final time and controller step");
        const auto output = fuelsim::test::read_final_exodus_results(argv[3]);
        const auto committed = fuelsim::rz::ProblemAccess::committed_state(problem);
        const auto& history = committed.contact_histories.at(0);
        if (history.size() != 3) throw std::runtime_error("B8.1 must commit three unique contact histories");
        for (std::size_t i = 0; i < 3; ++i) {
            const auto n = 2 * i + 1;
            const double elastic = output.nodal("contact_elastic_tangential_slip_interface")[n];
            const double bound = 8 * std::numeric_limits<double>::epsilon() *
                                 std::max(std::abs(elastic), std::abs(history[i].elastic_tangential_slip));
            if (output.nodal("contact_total_tangential_slip_interface")[n] != history[i].total_tangential_slip ||
                std::abs(elastic - history[i].elastic_tangential_slip) > bound ||
                output.nodal("contact_sliding_interface")[n] != (history[i].sliding ? 1.0 : 0.0))
                throw std::runtime_error("Production contact output differs from committed friction history");
        }
        const std::string temporary = std::string(argv[2]) + ".contract";
        const auto original = checkpoint_bytes(temporary, problem);
        fuelsim::TransientProblem restored(input.spatial, mesh);
        (void)fuelsim::restore_transient_checkpoint(temporary, restored);
        if (checkpoint_bytes(temporary, restored) != original)
            throw std::runtime_error("Friction checkpoint roundtrip changed serialized state");
        const auto full = problem.capture_state();
        auto changed = committed;
        changed.contact_histories[0][0].total_tangential_slip += 1e-4;
        fuelsim::rz::ProblemAccess::restore_committed_state(problem, changed);
        const auto half = problem.capture_state();
        fuelsim::TransientTimeOptions options{};
        options.displacement_time_absolute_tolerance = 1e-8;
        const auto estimate = problem.step_doubling_error(full, half, options);
        if (estimate.maximum != 0.0)
            throw std::runtime_error("Output-only accumulated slip must not change physical-state time error");
        changed.contact_histories[0][0].elastic_tangential_slip += 1e-4;
        fuelsim::rz::ProblemAccess::restore_committed_state(problem, changed);
        if (!(problem.step_doubling_error(full, problem.capture_state(), options).contact_friction > 1.0))
            throw std::runtime_error("Step-doubling must detect elastic slip that changes contact forces");
        problem.restore_state(full);
        if (checkpoint_bytes(temporary, problem) != original)
            throw std::runtime_error("Snapshot restore changed friction or material state");
        problem.begin_time_step({1.1, 1.0});
        auto trial = problem.committed_solution();
        trial[2 * mesh.nodes().size() + 1] += 1e-9;
        problem.validate_state(trial);
        fuelsim::ContributionWorkspace workspace;
        for (std::size_t i = 0; i < problem.contribution_count(); ++i)
            problem.evaluate_contribution(i, trial, workspace, true);
        problem.rollback_time_step();
        if (checkpoint_bytes(temporary, problem) != original)
            throw std::runtime_error("Trial evaluation or rollback changed complete committed state");
        auto wrong_input = input.spatial;
        for (auto& region : wrong_input.regions) region.rz_element_formulation = fuelsim::RzElementFormulation::quad4;
        fuelsim::TransientProblem wrong_element(wrong_input, mesh);
        bool rejected = false;
        try {
            (void)fuelsim::restore_transient_checkpoint(temporary, wrong_element);
        } catch (const std::exception&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Checkpoint signature did not distinguish quad4 and cax4t");
        if (std::remove(temporary.c_str()) != 0)
            throw std::runtime_error("Cannot remove temporary friction checkpoint");
        std::cout << "rz_friction_checkpoint_roundtrip=exact\nrz_friction_trial_rollback=exact\n"
                     "rz_friction_step_doubling=passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
