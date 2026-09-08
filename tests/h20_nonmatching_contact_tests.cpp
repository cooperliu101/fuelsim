#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (!condition)
        std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

bool test_friction_length(fuelsim::SpatialDefinition definition, const fuelsim::UnstructuredHex20Mesh& mesh) {
    definition.boundary_conditions.clear();
    definition.contacts[0].thermal = false;
    definition.contacts[0].mechanical = true;
    definition.contacts[0].friction_coefficient = 0.3;
    definition.contacts[0].friction_slip_tolerance = 1e-5;
    fuelsim::SteadyProblem problem(definition, mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    auto state = problem.initial_state();
    for (std::size_t node = 0; node < spatial.hex20_region_mesh(1).nodes().size(); ++node) {
        const auto global = spatial.global_node(1, node);
        state[spatial.dof(fuelsim::Field::displacement_x, global)] = -1e-5;
        state[spatial.dof(fuelsim::Field::displacement_y, global)] = 1e-6;
    }
    problem.validate_state(state);
    const auto summary = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
    double tangential_force = 0.0, normal_force = 0.0, maximum_gap_error = 0.0;
    for (const auto& node : summary) {
        tangential_force += node.tangential_contact_force[1];
        normal_force += node.normal_contact_force[0];
        maximum_gap_error = std::max(maximum_gap_error, std::abs(node.gap + 1e-5));
    }
    // Independently recovered at all 13 Abaqus nodes; see diagnose_c3d20rt_friction_length.py.
    constexpr double native_elastic_slip_length = 9.94936153005124e-6;
    const double expected = std::abs(normal_force) * 0.3 * 1e-6 / native_elastic_slip_length;
    const double error = std::abs(std::abs(tangential_force) - expected) / expected;
    std::cout << "h20_nonmatching_friction_length_relative_error=" << error << '\n';
    return check(error < 1e-10, "nonuniform face areas reproduce the independently measured elastic slip length")
           && check(maximum_gap_error < 1e-16, "coincident reference faces do not acquire a cancellation gap");
}

bool test_thermal_partition(fuelsim::SpatialDefinition definition, const fuelsim::UnstructuredHex20Mesh& mesh) {
    definition.boundary_conditions.clear();
    definition.contacts[0].thermal = true;
    definition.contacts[0].mechanical = false;
    definition.contacts[0].gap_conductivity = 0.001;
    definition.contacts[0].minimum_gap = 0.001;
    definition.regions[0].initial_temperature = 300.0;
    definition.regions[1].initial_temperature = 301.0;
    fuelsim::SteadyProblem problem(definition, mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto state = problem.initial_state();
    problem.validate_state(state);
    std::vector<double> assembled(problem.dof_count());
    double maximum_directional_error = 0.0;
    std::size_t thermal_constraint_count = 0;
    bool unique_partition_ownership = true;
    const auto first_partition = spatial.contribution_partition(0, 2);
    const auto second_partition = spatial.contribution_partition(1, 2);
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::thermal_contact)
            continue;
        ++thermal_constraint_count;
        const unsigned owners =
            static_cast<unsigned>(contribution >= first_partition.first && contribution < first_partition.second)
            + static_cast<unsigned>(contribution >= second_partition.first && contribution < second_partition.second);
        unique_partition_ownership = unique_partition_ownership && owners == 1;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size()), direction(dofs.size());
        for (std::size_t i = 0; i < dofs.size(); ++i) {
            local[i] = state[dofs[i]];
            direction[i] = (dofs[i] < spatial.dof(fuelsim::Field::displacement_x, 0) ? 1.0 : 0.001)
                           * std::sin(static_cast<double>(i + 1));
        }
        std::vector<double> residual, jacobian, plus = local, minus = local, rp, rm;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        constexpr double epsilon = 1e-4;
        for (std::size_t i = 0; i < local.size(); ++i) {
            assembled[dofs[i]] += residual[i];
            plus[i] += epsilon * direction[i];
            minus[i] -= epsilon * direction[i];
        }
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, rp, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, rm, nullptr);
        double difference = 0.0, reference = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double exact = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                exact += jacobian[row * local.size() + column] * direction[column];
            const double numerical = (rp[row] - rm[row]) / (2.0 * epsilon);
            difference += (exact - numerical) * (exact - numerical);
            reference += numerical * numerical;
        }
        maximum_directional_error = std::max(maximum_directional_error, std::sqrt(difference / reference));
    }
    double maximum_distribution_error = 0.0, balance = 0.0;
    for (std::size_t region = 0; region < 2; ++region) {
        const auto& local_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t node = 0; node < local_mesh.nodes().size(); ++node) {
            if (!local_mesh.temperature_nodes()[node])
                continue;
            const auto& coordinate = local_mesh.nodes()[node];
            double exact = 0.0;
            if (coordinate.x == 1.0) {
                if (region == 0) {
                    // Native constant-flux probe; primary transfer uses the identified contact sampling rule.
                    const std::array<double, 5> weights = {0.1254903810568766,
                        0.24925220819632693,
                        0.24976556426742502,
                        0.2507507226492291,
                        0.1247411238302501};
                    exact = -weights.at(static_cast<std::size_t>(std::lround(2.0 * coordinate.y)));
                } else
                    exact = coordinate.y == 0.0 ? 0.3 : (coordinate.y == 2.0 ? 0.2 : 0.5);
            }
            const double actual =
                assembled[spatial.dof(fuelsim::Field::temperature, spatial.global_temperature_node(region, node))];
            maximum_distribution_error = std::max(maximum_distribution_error, std::abs(actual - exact));
            balance += actual;
        }
    }
    std::cout << "h20_nonmatching_unit_heat_distribution_error=" << maximum_distribution_error
              << "\nh20_nonmatching_thermal_directional_error=" << maximum_directional_error << '\n';
    return check(thermal_constraint_count == 6 && unique_partition_ownership,
               "six shared-corner thermal constraints each belong to exactly one MPI contribution interval")
           && check(maximum_distribution_error < 1e-6,
               "nonmatching contact matches independently measured native constant heat at each node")
           && check(std::abs(balance) < 1e-12, "nonmatching thermal contact conserves heat")
           && check(maximum_directional_error < 1e-6, "partitioned thermal contact retains its geometric Jacobian");
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3)
        return 2;
    const auto definition = fuelsim::read_case_input(argv[1]);
    const auto mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::array<double, 4> primary_midpoints = {0.25, 0.75, 1.25, 1.75};
    std::array<std::size_t, 4> primary_midpoint_dofs{};
    for (std::size_t segment = 0; segment < primary_midpoints.size(); ++segment) {
        const auto& region = spatial.hex20_region_mesh(0);
        const auto found = std::find_if(region.nodes().begin(), region.nodes().end(), [&](const auto& point) {
            return std::abs(point.x - 1.0) < 1.0e-14 && std::abs(point.y - primary_midpoints[segment]) < 1.0e-14
                   && std::abs(point.z) < 1.0e-14;
        });
        if (found == region.nodes().end()) {
            std::cerr << "[FAIL] missing primary face midpoint node\n";
            return 1;
        }
        primary_midpoint_dofs[segment] = spatial.dof(fuelsim::Field::displacement_x,
            spatial.global_node(0, static_cast<std::size_t>(found - region.nodes().begin())));
    }
    const auto primary_participation = [&](const std::vector<double>& state) {
        problem.validate_state(state);
        std::array<std::size_t, 4> counts{};
        std::size_t constraint_count = 0, maximum_faces = 0;
        for (std::size_t contribution = spatial.volume_contribution_count();
            contribution < spatial.contribution_count();
            ++contribution) {
            if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
                continue;
            std::vector<std::size_t> dofs;
            problem.contribution_dofs(contribution, dofs);
            std::size_t faces = 0;
            for (std::size_t primary = 0; primary < primary_midpoint_dofs.size(); ++primary)
                if (std::find(dofs.begin(), dofs.end(), primary_midpoint_dofs[primary]) != dofs.end()) {
                    ++counts[primary];
                    ++faces;
                }
            if (faces == 0) {
                std::cerr << "[FAIL] nonmatching averaged constraint has no primary participation\n";
                return std::make_tuple(std::array<std::size_t, 4>{}, std::size_t{0}, std::size_t{0});
            }
            ++constraint_count;
            maximum_faces = std::max(maximum_faces, faces);
        }
        return std::make_tuple(counts, constraint_count, maximum_faces);
    };
    const auto participation = primary_participation(problem.initial_state());
    const std::array<std::size_t, 4>& primary_counts = std::get<0>(participation);
    const std::size_t constraint_count = std::get<1>(participation), maximum_faces = std::get<2>(participation);
    bool passed = check(
        constraint_count == 13
            && std::all_of(primary_counts.begin(), primary_counts.end(), [](std::size_t count) { return count > 0; })
            && maximum_faces > 1,
        "nonmatching HEX20 contact constructs thirteen averaged constraints spanning all four primary faces");
    std::vector<double> state(problem.dof_count());
    const auto output = fuelsim::test::read_final_exodus_results(argv[2]);
    if (output.nodes.size() != mesh.nodes().size())
        throw std::runtime_error("HEX20 result node count differs");
    const std::array<std::string, 3> displacement_names = {"displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& local_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < local_mesh.nodes().size(); ++local) {
            const auto source = local_mesh.source_node_ids()[local];
            const auto global = spatial.global_node(region, local);
            for (std::size_t component = 0; component < 3; ++component)
                state.at(spatial.field_layout()[component + 1].begin + global) =
                    output.nodal(displacement_names[component]).at(source);
        }
        for (std::size_t local = 0; local < local_mesh.nodes().size(); ++local) {
            if (!local_mesh.temperature_nodes()[local])
                continue;
            const auto source = local_mesh.source_node_ids()[local];
            state.at(spatial.dof(fuelsim::Field::temperature, spatial.global_temperature_node(region, local))) =
                output.nodal("temperature").at(source);
        }
    }
    problem.validate_state(state);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, state);
    std::array<double, 3> contact_balance{};
    double maximum_jacobian_directional_error = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        std::vector<double> residual;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local)
            local_state[local] = state[dofs[local]];
        std::vector<double> jacobian;
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local_state.size()), plus = local_state, minus = local_state;
        constexpr double perturbation = 1.0e-10;
        for (std::size_t local = 0; local < direction.size(); ++local) {
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] += perturbation * direction[local];
            minus[local] -= perturbation * direction[local];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local_state.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[row * local_state.size() + column] * direction[column];
            const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - finite_difference, 2);
            reference_squared += finite_difference * finite_difference;
        }
        maximum_jacobian_directional_error =
            std::max(maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
        for (std::size_t local = 0; local < dofs.size(); ++local) {
            const std::size_t global = dofs[local];
            for (std::size_t component = 0; component < contact_balance.size(); ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (global >= field.begin && global < field.end)
                    contact_balance[component] += residual[local];
            }
        }
    }
    passed = check(interface.projected_contact_nodes > 0 && interface.unprojected_contact_nodes == 0,
                 "nonmatching HEX20 contact projects every secondary boundary node")
             && check(interface.active_contact_nodes > 0 && interface.total_contact_force > 0.0,
                 "nonmatching HEX20 contact develops positive compressive interface force")
             && check(std::abs(contact_balance[0]) < 1.0e-8 && std::abs(contact_balance[1]) < 1.0e-8
                          && std::abs(contact_balance[2]) < 1.0e-8,
                 "nonmatching HEX20 contact residual is action-reaction conservative")
             && check(maximum_jacobian_directional_error < 1.0e-7,
                 "nonmatching HEX20 averaged-contact Jacobian matches a centered directional difference")
             && passed;
    const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
    passed = check(histories.size() == constraint_count,
                 "nonmatching HEX20 contact stores one transaction slot for every averaged constraint")
             && passed;
    {
        std::cout << "h20_24_primary_face_participation_counts=" << primary_counts[0] << ',' << primary_counts[1] << ','
                  << primary_counts[2] << ',' << primary_counts[3] << '\n'
                  << "h20_24_averaged_constraint_count=" << constraint_count << '\n'
                  << "h20_24_maximum_primary_faces_per_constraint=" << maximum_faces << '\n'
                  << "h20_24_contact_jacobian_directional_error=" << maximum_jacobian_directional_error << '\n'
                  << "h20_24_projected_secondary_nodes=" << interface.projected_contact_nodes << '\n'
                  << "h20_24_active_contact_nodes=" << interface.active_contact_nodes << '\n'
                  << "h20_24_total_contact_force=" << interface.total_contact_force << '\n';
    }
    passed = test_thermal_partition(definition.spatial, mesh) && passed;
    passed = test_friction_length(definition.spatial, mesh) && passed;
    return passed ? 0 : 1;
}
