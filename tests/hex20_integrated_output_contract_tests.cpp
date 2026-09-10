#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {
double mechanical_contact_directional_error(fuelsim::TransientProblem& problem,
    const std::vector<double>& global_state,
    double perturbation) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> trial_state = global_state;
    const std::vector<std::size_t> source_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    const std::vector<fuelsim::CartesianContactNodeSummary> summaries =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, global_state);
    if (source_nodes.size() != summaries.size())
        throw std::logic_error("C3D20T contact Jacobian check has inconsistent secondary-node maps");
    std::map<std::size_t, std::size_t> source_to_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const fuelsim::Hex20RegionMesh& mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local)
            source_to_global.emplace(mesh.source_node_ids()[local], spatial.global_node(region, local));
    }
    constexpr double sliding_trial_increment = 1.0e-5;
    for (std::size_t node = 0; node < source_nodes.size(); ++node) {
        const auto& slip = summaries[node].tangential_slip;
        const double magnitude = std::hypot(slip[0], slip[1], slip[2]);
        if (!(magnitude > 0.0))
            continue;
        const std::size_t global = source_to_global.at(source_nodes[node]);
        trial_state[spatial.dof(fuelsim::Field::displacement_x, global)] +=
            sliding_trial_increment * slip[0] / magnitude;
        trial_state[spatial.dof(fuelsim::Field::displacement_y, global)] +=
            sliding_trial_increment * slip[1] / magnitude;
        trial_state[spatial.dof(fuelsim::Field::displacement_z, global)] +=
            sliding_trial_increment * slip[2] / magnitude;
    }
    spatial.validate_state(trial_state);
    double maximum_error = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size()), direction(dofs.size()), plus(dofs.size()), minus(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local) {
            local_state[local] = trial_state[dofs[local]];
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] = local_state[local] + perturbation * direction[local];
            minus[local] = local_state[local] - perturbation * direction[local];
        }
        std::vector<double> residual, jacobian, plus_residual, minus_residual;
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        const auto evaluate_perturbed = [&](const std::vector<double>& local, std::vector<double>& result) {
            auto perturbed = trial_state;
            for (std::size_t i = 0; i < dofs.size(); ++i)
                perturbed[dofs[i]] = local[i];
            spatial.validate_state(perturbed);
            std::vector<std::size_t> perturbed_dofs;
            problem.contribution_dofs(contribution, perturbed_dofs);
            if (perturbed_dofs != dofs)
                throw std::runtime_error("Integrated contact derivative crossed a candidate support transition");
            spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, result, nullptr);
        };
        evaluate_perturbed(plus, plus_residual);
        evaluate_perturbed(minus, minus_residual);
        spatial.validate_state(trial_state);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local_state.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[row * local_state.size() + column] * direction[column];
            const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - finite_difference, 2);
            reference_squared += finite_difference * finite_difference;
        }
        if (reference_squared > 0.0) {
            const double error = std::sqrt(difference_squared / reference_squared);
            maximum_error = std::max(maximum_error, error);
        }
    }
    spatial.validate_state(global_state);
    return maximum_error;
}

void run(const std::string& input_path, const std::string& checkpoint_path, const std::string& output_path) {
    const auto definition = fuelsim::read_case_input(input_path);
    const auto mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    fuelsim::restore_transient_checkpoint(checkpoint_path, problem);
    const auto output = fuelsim::test::read_final_exodus_results(output_path);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& state = problem.committed_solution();
    if (output.time != problem.committed_time() || output.step_count != 2)
        throw std::runtime_error("C3D20T output and checkpoint physical times differ");
    std::map<std::size_t, std::size_t> displacement_nodes, temperature_nodes;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const auto source = region_mesh.source_node_ids()[local];
            const auto global = spatial.global_node(region, local);
            displacement_nodes.emplace(source, global);
            for (std::size_t component = 0; component < 3; ++component)
                if (output.nodal("displacement_" + std::string(1, "xyz"[component])).at(source)
                    != state.at(spatial.field_layout()[component + 1].begin + global))
                    throw std::runtime_error("C3D20T displacement output differs from committed state");
            if (region_mesh.temperature_nodes()[local]) {
                const auto temperature_node = spatial.global_temperature_node(region, local);
                temperature_nodes.emplace(source, temperature_node);
                if (output.nodal("temperature").at(source) != state.at(temperature_node))
                    throw std::runtime_error("C3D20T corner-temperature output differs from committed state");
            }
        }
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const auto source = region_mesh.source_element_ids()[element];
            const auto& geometry = spatial.hex20_region_element_geometry(region, element);
            for (std::size_t q = 0; q < 27; ++q) {
                const auto& point = geometry.mechanical_points[q];
                std::array<double, 3> current = {point.position.x, point.position.y, point.position.z};
                for (std::size_t local = 0; local < 20; ++local) {
                    const auto global = spatial.global_node(region, region_mesh.elements()[element].nodes[local]);
                    for (std::size_t component = 0; component < 3; ++component)
                        current[component] += point.displacement_shape[local]
                                              * state.at(spatial.field_layout()[component + 1].begin + global);
                }
                for (std::size_t component = 0; component < 3; ++component)
                    if (output.element("current_" + std::string(1, "xyz"[component]) + "_q" + std::to_string(q))
                            .at(source)
                        != current[component])
                        throw std::runtime_error(
                            "C3D20T material-point output coordinate differs from committed geometry");
            }
        }
    }
    if (displacement_nodes.size() != 112 || temperature_nodes.size() != 40)
        throw std::runtime_error("C3D20T mixed-order source-node maps changed");
    if (!definition.spatial.contacts.empty()) {
        const bool frictional = definition.spatial.contacts.front().friction_coefficient > 0.0;
        const double error = mechanical_contact_directional_error(problem, state, frictional ? 1e-9 : 1e-8);
        std::cout << "c3d20t_contact_jacobian_directional_error=" << error << '\n';
        if (!(error < 2e-5))
            throw std::runtime_error("C3D20T finite-sliding contact Jacobian failed centered differences");
    }
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4)
        return 2;
    try {
        run(argv[1], argv[2], argv[3]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
