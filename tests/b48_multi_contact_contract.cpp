#include "core/problem_backend_access.hpp"
#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
bool histories_equal(const std::vector<std::vector<fuelsim::ContactPointHistory>>& first,
    const std::vector<std::vector<fuelsim::ContactPointHistory>>& second) {
    if (first.size() != second.size())
        return false;
    for (std::size_t pair = 0; pair < first.size(); ++pair) {
        if (first[pair].size() != second[pair].size())
            return false;
        for (std::size_t point = 0; point < first[pair].size(); ++point) {
            const auto& left = first[pair][point];
            const auto& right = second[pair][point];
            if (left.elastic_tangential_slip != right.elastic_tangential_slip || left.sliding != right.sliding
                || left.normal_multiplier != right.normal_multiplier
                || left.cartesian_elastic_tangential_slip != right.cartesian_elastic_tangential_slip
                || left.cartesian_total_tangential_slip != right.cartesian_total_tangential_slip
                || left.cartesian_tangent_basis_initialized != right.cartesian_tangent_basis_initialized
                || left.cartesian_contact_normal != right.cartesian_contact_normal
                || left.cartesian_contact_tangent_first != right.cartesian_contact_tangent_first)
                return false;
        }
    }
    return true;
}

double jacobian_error(const fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double perturbation = 1.0e-7;
        for (std::size_t index = 0; index < local.size(); ++index) {
            direction[index] = std::sin(static_cast<double>(index + 1));
            plus[index] += perturbation * direction[index];
            minus[index] -= perturbation * direction[index];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                analytic += jacobian[row * local.size() + column] * direction[column];
            const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - reference, 2);
            reference_squared += reference * reference;
        }
        if (reference_squared > 0.0)
            maximum = std::max(maximum, std::sqrt(difference_squared / reference_squared));
    }
    return maximum;
}

double contact_action_reaction_error(const fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
        std::vector<double> residual;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, nullptr);
        std::array<double, 3> resultant{};
        const auto& fields = spatial.field_layout();
        for (std::size_t row = 0; row < dofs.size(); ++row)
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = fields.at(component + 1);
                if (dofs[row] >= field.begin && dofs[row] < field.end)
                    resultant[component] += residual[row];
            }
        for (double value : resultant)
            maximum = std::max(maximum, std::abs(value));
    }
    return maximum;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 5)
        return 2;
    try {
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex8(input.mesh_file);
        fuelsim::TransientProblem full(input.spatial, mesh), restarted(input.spatial, mesh);
        const double full_step = fuelsim::restore_transient_checkpoint(argv[2], full);
        const double restart_step = fuelsim::restore_transient_checkpoint(argv[3], restarted);
        if (full.committed_time() != 4 || restarted.committed_time() != 4 || full_step != restart_step
            || full.committed_solution() != restarted.committed_solution()
            || !histories_equal(fuelsim::cartesian::ProblemAccess::committed_contact_histories(full),
                fuelsim::cartesian::ProblemAccess::committed_contact_histories(restarted)))
            throw std::runtime_error("B4.8 restart must exactly preserve the complete nodal and friction histories");
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(full);
        auto trial = full.committed_solution();
        for (std::size_t pair = 0; pair < 2; ++pair) {
            const std::size_t region = pair == 0 ? 1 : 3;
            for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local) {
                const auto global = spatial.global_node(region, local);
                trial[spatial.dof(fuelsim::Field::displacement_y, global)] += pair == 0 ? 1.2e-3 : -1.2e-3;
                trial[spatial.dof(fuelsim::Field::displacement_z, global)] += 1.4e-3;
            }
        }
        spatial.validate_state(trial);
        const double derivative = jacobian_error(full, trial),
                     conservation = contact_action_reaction_error(full, trial);
        std::cout << "b48_contact_jacobian_directional_error=" << derivative << '\n'
                  << "b48_contact_action_reaction_maximum_absolute=" << conservation << '\n';
        if (!(derivative < 1e-5 && conservation < 1e-10))
            throw std::runtime_error("B4.8 local contact derivative or conservation check failed");
        // Reproduce MPI's partially refreshed projection cache without MPI or
        // another solve. Both commits start from the same production checkpoint.
        for (std::size_t rank = 0; rank < 2; ++rank) {
            fuelsim::TransientProblem complete(input.spatial, mesh), partial(input.spatial, mesh);
            fuelsim::restore_transient_checkpoint(argv[4], complete);
            fuelsim::restore_transient_checkpoint(argv[4], partial);
            if (complete.committed_time() != 2 || partial.committed_time() != 2)
                throw std::runtime_error("B4.8 partition regression requires the second-step checkpoint");
            complete.begin_time_step({4.0, 1.0, true});
            partial.begin_time_step({4.0, 1.0, true});
            const auto& state = full.committed_solution();
            complete.validate_state(state);
            const auto interval = partial.contribution_partition(rank, 2);
            partial.validate_local_state(interval.first, interval.second, state);
            complete.commit_time_step(state);
            partial.commit_time_step(state);
            if (fuelsim::BackendAccess::committed_raw_residual(complete)
                    != fuelsim::BackendAccess::committed_raw_residual(partial)
                || !histories_equal(fuelsim::cartesian::ProblemAccess::committed_contact_histories(complete),
                    fuelsim::cartesian::ProblemAccess::committed_contact_histories(partial)))
                throw std::runtime_error("B4.8 partial projection validation changed committed reactions or history");
        }
        std::cout << "[PASS] B4.8 commit refreshes complete contact projections after either local partition\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
