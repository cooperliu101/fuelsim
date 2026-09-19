#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace {
using namespace fuelsim;

std::unique_ptr<TransientProblem> make_problem(const FuelSimCaseDefinition& input) {
    if (input.geometry == CaseGeometry::generalized_plane_strain)
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_plane_quad8(input.mesh_file));
    if (input.geometry == CaseGeometry::axisymmetric_1d)
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_bar2(input.mesh_file));
    if (input.geometry == CaseGeometry::cartesian_3d) {
        if (exodus_uses_hex20(input.mesh_file))
            return std::make_unique<TransientProblem>(input.spatial, read_exodus_hex20(input.mesh_file));
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_hex8(input.mesh_file));
    }
    if (exodus_uses_quad8(input.mesh_file))
        return std::make_unique<TransientProblem>(input.spatial, read_exodus_quad8(input.mesh_file));
    return std::make_unique<TransientProblem>(input.spatial, read_exodus_quad4(input.mesh_file));
}

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<char> checkpoint(const std::string& path, const TransientProblem& problem) {
    write_transient_checkpoint(path, problem, 0.25);
    std::ifstream file(path, std::ios::binary);
    require(static_cast<bool>(file), "Cannot read contract checkpoint");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void check(const FuelSimCaseDefinition& input, const std::string& path) {
    auto instance = make_problem(input);
    auto& problem = *instance;
    const auto identity = problem.discretization_identity();
    const auto snapshot = problem.capture_state();
    const auto initial = checkpoint(path, problem);
    const bool fixed = problem.contribution_metadata_is_fixed();
    const bool dynamic = problem.jacobian_sparsity_is_state_dependent();
    const auto candidate = problem.committed_solution();
    problem.begin_time_step({0.25, 0.0, true});
    problem.validate_state(candidate);
    ContributionWorkspace workspace;
    std::vector<unsigned char> pattern;
    std::vector<unsigned> owners(problem.contribution_count(), 0);
    for (std::size_t rank = 0; rank < 2; ++rank) {
        const auto range = problem.contribution_partition(rank, 2);
        const auto required = problem.required_state_dofs(range.first, range.second);
        std::vector<double> shadow(problem.dof_count(), std::numeric_limits<double>::quiet_NaN());
        for (const auto dof : required) {
            require(dof < problem.dof_count(), "Shadow degree of freedom exceeds the global layout");
            shadow[dof] = candidate[dof];
        }
        problem.validate_local_state(range.first, range.second, shadow);
        for (std::size_t index = range.first; index < range.second; ++index) {
            ++owners[index];
            problem.evaluate_contribution(index, shadow, workspace, false);
            const auto first = workspace.residual;
            const auto dofs = workspace.dofs;
            problem.contribution_jacobian_pattern(index, pattern);
            require(pattern.size() == dofs.size() * dofs.size(), "Local Jacobian pattern has the wrong size");
            problem.evaluate_contribution(index, candidate, workspace, true);
            problem.evaluate_contribution(index, candidate, workspace, false);
            require(first == workspace.residual, "Intervening Jacobian evaluation changed the repeated residual");
            require(dofs == workspace.dofs, "Repeated evaluation changed the contribution map");
        }
    }
    for (const auto count : owners)
        require(count == 1, "A contribution must have exactly one partition owner");
    for (std::size_t index = 0; index < problem.sparsity_contribution_count(); ++index) {
        std::vector<std::size_t> dofs;
        problem.sparsity_contribution_dofs(index, dofs);
        problem.sparsity_contribution_jacobian_pattern(index, pattern);
        require(pattern.size() == dofs.size() * dofs.size(), "Sparsity contribution pattern has the wrong size");
        for (const auto dof : dofs)
            require(dof < problem.dof_count(), "Sparsity degree of freedom exceeds the global layout");
    }
    problem.rollback_time_step();
    require(initial == checkpoint(path, problem), "Evaluation or rollback contaminated committed state");
    require(identity == problem.discretization_identity(), "Time controls changed discretization identity");
    require(fixed == problem.contribution_metadata_is_fixed(), "Time controls changed metadata lifetime");
    require(dynamic == problem.jacobian_sparsity_is_state_dependent(), "Time controls changed graph policy");

    auto other = make_problem(input);
    bool rejected = false;
    try {
        other->restore_state(snapshot);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "A snapshot from another problem was accepted");
    require(initial == checkpoint(path, *other), "Rejected snapshot changed its destination");

    if (!problem.uses_radial_gps()) {
        // Remote rejection of staged materials, contacts and diagnostics must
        // leave exactly the same checkpoint payload after rollback.
        {
            problem.begin_time_step({0.25, 0.0, true});
            unsigned phase = 0;
            rejected = false;
            try {
                problem.commit_time_step(candidate, 0, problem.contribution_count(), [&](std::vector<double>& values) {
                    ++phase;
                    values.at(0) = 1.0;
                });
            } catch (const std::domain_error&) {
                rejected = true;
            }
            require(rejected && phase == 1, "Remote commit failure did not follow the shared protocol");
            problem.rollback_time_step();
            require(initial == checkpoint(path, problem), "Failed publication left partially committed state");
        }
    }
    auto invalid = BackendAccess::committed_state(problem);
    if (problem.uses_radial_gps())
        invalid.material_histories.emplace_back();
    else
        invalid.radial_material_histories.emplace_back();
    rejected = false;
    try {
        BackendAccess::restore_committed_state(problem, std::move(invalid));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "An incompatible backend history was accepted");
    require(initial == checkpoint(path, problem), "Rejected history restore changed committed state");

    auto updated = candidate;
    for (const auto& field : problem.field_layout())
        if (field.category == FieldCategory::thermal)
            for (std::size_t dof = field.begin; dof < field.end; ++dof)
                updated[dof] += 1.0;
    problem.begin_time_step({0.25, 0.0, true});
    problem.commit_time_step(updated);
    const auto accepted = problem.capture_state();
    const auto accepted_file = checkpoint(path, problem);
    require(accepted_file != initial, "A successful commit did not publish new state");
    problem.restore_state(snapshot);
    require(initial == checkpoint(path, problem), "An old snapshot was overwritten by trial or accepted buffers");
    problem.restore_state(accepted);
    require(accepted_file == checkpoint(path, problem), "Accepted snapshot did not restore complete state");
    if (BackendAccess::uses_thermal(problem)) {
        const auto state = BackendAccess::committed_state(problem);
        require(problem.field_layout().size() == 1 && problem.field_layout().front().category == FieldCategory::thermal,
            "Pure thermal backend allocated mechanical degrees of freedom");
        require(state.material_histories.empty() && state.quad8_material_histories.empty()
                    && state.radial_material_histories.empty() && state.cartesian_material_histories.empty()
                    && state.contact_histories.empty(),
            "Pure thermal backend allocated mechanical histories");
    }
    std::cout << "Contribution metadata, callback purity, identity and publication transactions passed\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Expected input card and scratch checkpoint");
        check(fuelsim::read_case_input(argv[1]), argv[2]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
