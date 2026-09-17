#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/problem_signature.hpp"
#include "io/results_io.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<double> values(const fuelsim::TransientProblem& problem) {
    const auto state = fuelsim::BackendAccess::committed_state(problem);
    std::vector<double> result{state.time, state.load_factor, state.previous_time};
    for (const auto* field :
        {&state.solution, &state.previous_solution, &state.raw_residual, &state.external_load_residual})
        result.insert(result.end(), field->begin(), field->end());
    for (const auto& field : fuelsim::transient_conservation_fields)
        result.push_back(state.conservation.*field.member);
    for (const auto& region : state.cartesian_material_histories)
        for (const auto& element : region)
            for (const auto& point : element) {
                for (const auto* tensor : {&point.elastic_strain, &point.plastic_strain, &point.creep_strain})
                    result.insert(result.end(), tensor->begin(), tensor->end());
                result.insert(result.end(),
                    {point.stress.xx,
                        point.stress.yy,
                        point.stress.zz,
                        point.stress.xy,
                        point.stress.yz,
                        point.stress.xz,
                        point.equivalent_plastic_strain,
                        point.equivalent_creep_strain});
            }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    using namespace fuelsim;
    try {
        require(argc == 3, "Expected source input and isolated test directory");
        const auto input = read_case_input(argv[1]);
        const auto mesh = read_exodus_plane_quad8(input.mesh_file);
        std::filesystem::create_directories(argv[2]);
        const auto mesh_path = (std::filesystem::path(argv[2]) / "roundtrip.e").string();
        write_exodus_plane_quad8(mesh_path, mesh);
        const auto restored_mesh = read_exodus_plane_quad8(mesh_path);
        require(restored_mesh.nodes() == mesh.nodes(), "Plane mesh roundtrip preserves coordinates");
        for (std::size_t e = 0; e < mesh.elements().size(); ++e)
            require(restored_mesh.elements()[e].nodes == mesh.elements()[e].nodes,
                "Plane mesh roundtrip preserves connectivity");
        TransientProblem problem(input.spatial, mesh);
        const auto initial = problem.capture_state();
        const auto before = values(problem);
        problem.begin_time_step({1, 1, false});
        auto trial = problem.committed_solution();
        for (const auto& condition : problem.dirichlet_conditions())
            trial[condition.dof] = condition.value;
        problem.validate_state(trial);
        std::vector<double> complete(problem.dof_count(), 0), partitioned(problem.dof_count(), 0);
        std::vector<unsigned> owners(problem.contribution_count(), 0);
        ContributionWorkspace workspace;
        for (std::size_t index = 0; index < problem.contribution_count(); ++index) {
            problem.evaluate_contribution(index, trial, workspace, true);
            for (std::size_t i = 0; i < workspace.dofs.size(); ++i)
                complete[workspace.dofs[i]] += workspace.residual[i];
        }
        for (std::size_t rank = 0; rank < 2; ++rank) {
            const auto range = problem.contribution_partition(rank, 2);
            auto shadow = std::vector<double>(trial.size(), std::numeric_limits<double>::quiet_NaN());
            for (auto dof : problem.required_state_dofs(range.first, range.second))
                shadow[dof] = trial[dof];
            problem.validate_local_state(range.first, range.second, shadow);
            for (auto index = range.first; index < range.second; ++index) {
                ++owners[index];
                problem.evaluate_contribution(index, shadow, workspace, true);
                for (std::size_t i = 0; i < workspace.dofs.size(); ++i)
                    partitioned[workspace.dofs[i]] += workspace.residual[i];
            }
        }
        for (auto count : owners)
            require(count == 1, "Each volume and contact contribution belongs to exactly one rank");
        require(complete == partitioned, "Owned contributions with only required shadow fields recover full residual");
        problem.commit_time_step(trial);
        const auto committed = values(problem);
        require(committed != before, "Commit must update nodal state, time and material history");
        const auto checkpoint = (std::filesystem::path(argv[2]) / "plane.chk").string();
        write_transient_checkpoint(checkpoint, problem, .25);
        TransientProblem restored(input.spatial, restored_mesh);
        require(restore_transient_checkpoint(checkpoint, restored) == .25, "Checkpoint restores next time step");
        require(values(restored) == committed,
            "Checkpoint restores all control fields and nine-point histories exactly");
        problem.begin_time_step({2, 1, false});
        trial = problem.committed_solution();
        const auto& spatial = BackendAccess::plane_spatial(problem);
        trial[spatial.section_dof(0, 0)] = -1;
        bool rejected = false;
        try {
            problem.validate_state(trial);
        } catch (const std::domain_error&) {
            rejected = true;
        }
        require(rejected, "Negative plane thickness rejects the trial state");
        problem.rollback_time_step();
        require(values(problem) == committed, "Rejected time step preserves the complete accepted transaction");
        problem.restore_state(initial);
        require(values(problem) == before, "Snapshot restore recovers fields, controls, history and time");
        auto changed = input.spatial;
        changed.contacts[0].gap_conductance += 1;
        TransientProblem different(changed, mesh);
        rejected = false;
        try {
            (void)restore_transient_checkpoint(checkpoint, different);
        } catch (const std::exception&) {
            rejected = true;
        }
        require(rejected, "Checkpoint rejects a changed thermal contact law");
        changed = input.spatial;
        changed.generalized_plane_strain[0].prescribed[0] = .001;
        TransientProblem changed_constraint(changed, mesh);
        require(transient_problem_signature(changed_constraint) != transient_problem_signature(problem),
            "Checkpoint signature includes section constraint values");
        changed = input.spatial;
        changed.generalized_plane_strain[0].blocks.push_back("upper");
        changed.generalized_plane_strain.pop_back();
        TransientProblem shared(changed, mesh);
        const auto& shared_layout = BackendAccess::plane_spatial(shared);
        std::vector<std::size_t> lower, upper;
        shared_layout.contribution_dofs(0, lower);
        shared_layout.contribution_dofs(1, upper);
        require(shared_layout.section_count() == 1 && shared.dof_count() + 3 == problem.dof_count(),
            "Multiple blocks in one section share exactly three controls");
        for (std::size_t i = 20; i < 23; ++i)
            require(lower[i] == upper[i], "Every block uses its section controls");
        require(std::abs(shared_layout.section_origin(0)[1] - .00005) < 1e-14,
            "Section origin is the area centroid of all assigned blocks");
        for (unsigned failure = 0; failure < 4; ++failure) {
            changed = input.spatial;
            if (failure == 0)
                changed.generalized_plane_strain[0].blocks.push_back("upper");
            else if (failure == 1)
                changed.generalized_plane_strain.pop_back();
            else if (failure == 2)
                changed.generalized_plane_strain[0].blocks = {"missing"};
            else
                changed.generalized_plane_strain[0].initial_thickness = 0;
            rejected = false;
            try {
                TransientProblem invalid(changed, mesh);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            require(rejected, "Invalid section ownership or thickness must be rejected");
        }
        auto connected_elements = mesh.elements();
        connected_elements[1].nodes[0] = connected_elements[0].nodes[3];
        UnstructuredPlaneQuad8Mesh connected_mesh(mesh.nodes(),
            connected_elements,
            mesh.element_block_ids(),
            mesh.element_blocks(),
            mesh.node_sets(),
            mesh.side_sets());
        rejected = false;
        try {
            TransientProblem invalid(input.spatial, connected_mesh);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "Shared field nodes cannot belong to different sections");
        std::cout << "Plane mesh, snapshot, checkpoint and rejected-step transaction checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
