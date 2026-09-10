#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

// Internal directional-derivative and conservation contracts at the production solution.
// This executable deliberately has no solver dependency and never solves the case.
void verify(const std::string& card, const std::string& output) {
    const auto definition = fuelsim::read_case_input(card);
    const auto mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const auto result = fuelsim::test::read_final_exodus_results(output);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    auto state = problem.initial_state();
    const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& local_mesh = spatial.region_mesh(region);
        for (std::size_t local = 0; local < local_mesh.nodes().size(); ++local)
            for (std::size_t field = 0; field < names.size(); ++field)
                state.at(spatial.field_layout()[field].begin + spatial.global_node(region, local)) =
                    result.nodal(names[field]).at(local_mesh.source_node_ids()[local]);
    }
    problem.validate_state(state);
    std::size_t mechanical_contributions = 0;
    std::array<double, 3> residual_balance{};
    double maximum_jacobian_directional_error = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        ++mechanical_contributions;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double perturbation = 1.0e-10;
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
            const std::size_t global = dofs[row];
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (global >= field.begin && global < field.end)
                    residual_balance[component] += residual[row];
            }
        }
        maximum_jacobian_directional_error =
            std::max(maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
    }

    const auto nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    std::cout << "mechanical_constraint_count=" << mechanical_contributions
              << " contact_jacobian_directional_error=" << maximum_jacobian_directional_error
              << " residual_balance=" << residual_balance[0] << ',' << residual_balance[1] << ',' << residual_balance[2]
              << '\n';
    if (mechanical_contributions != nodes.size() || nodes.empty() || !(maximum_jacobian_directional_error < 1.0e-7)
        || !(std::abs(residual_balance[0]) < 1.0e-8) || !(std::abs(residual_balance[1]) < 1.0e-8)
        || !(std::abs(residual_balance[2]) < 1.0e-8))
        throw std::runtime_error("B3.9 unique constraint, directional derivative or action-reaction contract failed");
}

int main(int argc, char** argv) {
    if (argc != 3)
        return 2;
    try {
        verify(argv[1], argv[2]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
