#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char** argv) {
    using namespace fuelsim;
    try {
        if (argc != 2)
            throw std::invalid_argument("Expected the static sliding contact input");
        const auto input = read_case_input(argv[1]);
        const auto mesh = read_exodus_plane_quad8(input.mesh_file);
        TransientProblem problem(input.spatial, mesh);
        problem.begin_time_step({1.0, 1.0, false});
        auto state = problem.committed_solution();
        for (const auto& condition : problem.dirichlet_conditions())
            state[condition.dof] = condition.value;
        const auto& spatial = BackendAccess::plane_spatial(problem);
        std::vector<double> direction(state.size(), 0.0);
        for (std::size_t r = 0; r < spatial.region_count(); ++r)
            for (std::size_t e = 0; e < spatial.region_element_count(r); ++e) {
                std::vector<std::size_t> dofs;
                spatial.contribution_dofs(spatial.region_element_offset(r) + e, dofs);
                const auto& geometry = spatial.geometry(r, e);
                for (std::size_t n = 0; n < 8; ++n) {
                    const auto& x = geometry.coordinates[n];
                    // Affine motion keeps the shared primary edge tangents continuous.
                    direction[dofs[4 + n]] = .0001 * (r == 0 ? 1.0 : -2.0) + .003 * x[0];
                    direction[dofs[12 + n]] = .0001 * (r == 0 ? -1.0 : 1.0) + .002 * x[0];
                    if (n < 4)
                        direction[dofs[n]] = 5.0 * (r == 0 ? 1.0 : -1.0) + 20.0 * x[0];
                }
                for (std::size_t n = 0; n < 3; ++n)
                    direction[dofs[20 + n]] = .001 * static_cast<double>(n + 1);
            }
        const auto assemble = [&](const std::vector<double>& values, std::size_t ranks, bool tangent) {
            std::vector<double> result(2 * state.size(), 0.0);
            std::vector<unsigned> owners(problem.contribution_count(), 0);
            for (std::size_t rank = 0; rank < ranks; ++rank) {
                const auto range = problem.contribution_partition(rank, ranks);
                std::vector<double> shadow(state.size(), std::numeric_limits<double>::quiet_NaN());
                for (auto dof : problem.required_state_dofs(range.first, range.second))
                    shadow[dof] = values[dof];
                problem.validate_local_state(range.first, range.second, shadow);
                ContributionWorkspace work;
                for (auto index = std::max(range.first, spatial.contact_offset()); index < range.second; ++index) {
                    ++owners[index];
                    problem.evaluate_contribution(index, shadow, work, tangent);
                    for (std::size_t i = 0; i < work.dofs.size(); ++i) {
                        result[work.dofs[i]] += work.residual[i];
                        if (tangent)
                            for (std::size_t j = 0; j < work.dofs.size(); ++j)
                                result[state.size() + work.dofs[i]] +=
                                    work.jacobian[i * work.dofs.size() + j] * direction[work.dofs[j]];
                    }
                }
            }
            for (auto index = spatial.contact_offset(); index < owners.size(); ++index)
                if (owners[index] != 1)
                    throw std::runtime_error("Segment contribution ownership must be unique");
            return result;
        };
        const auto exact = assemble(state, 1, true);
        for (std::size_t ranks : {std::size_t{2}, std::size_t{4}})
            if (assemble(state, ranks, true) != exact)
                throw std::runtime_error("Segment residual and tangent require only declared shadow fields");
        auto plus = state, minus = state;
        constexpr double epsilon = 1e-5;
        for (std::size_t i = 0; i < state.size(); ++i) {
            plus[i] += epsilon * direction[i];
            minus[i] -= epsilon * direction[i];
        }
        const auto upper = assemble(plus, 1, false), lower = assemble(minus, 1, false);
        for (std::size_t i = 0; i < state.size(); ++i) {
            const double derivative = exact[state.size() + i];
            const double difference = (upper[i] - lower[i]) / (2 * epsilon);
            if (std::abs(difference - derivative) > 2e-6 * std::max(1.0, std::abs(derivative)))
                throw std::runtime_error("Moving segment boundary tangent disagrees with central differences");
        }
        std::cout << "Moving contact segments: tangent and two/four partition ownership passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
