#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool check_projection(const std::string& input_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::read_exodus_quad4(definition.mesh_file);
    bool passed = check(source.side_set("fuel_right").sides.size() != source.side_set("clad_left").sides.size(),
        "shared model retains nonmatching contact segmentation");
    fuelsim::SteadyProblem problem(definition.spatial, source);
    const std::vector<std::size_t> secondary_sources =
        fuelsim::rz::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    std::vector<std::size_t> ordered_secondary_sources = secondary_sources;
    std::sort(ordered_secondary_sources.begin(),
        ordered_secondary_sources.end(),
        [&source](std::size_t lhs, std::size_t rhs) { return source.nodes().at(lhs).z < source.nodes().at(rhs).z; });
    const std::size_t sliding_source = ordered_secondary_sources.at(1);
    std::size_t sliding_global_node = problem.dof_count();
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const std::vector<std::size_t>& source_nodes =
            fuelsim::rz::ProblemAccess::region_mesh(problem, region).source_node_ids();
        const auto found = std::find(source_nodes.begin(), source_nodes.end(), sliding_source);
        if (found == source_nodes.end())
            continue;
        sliding_global_node = fuelsim::rz::ProblemAccess::region_node_offset(problem, region)
                              + static_cast<std::size_t>(found - source_nodes.begin());
        break;
    }
    if (sliding_global_node == problem.dof_count())
        throw std::logic_error("M3.3 could not locate a secondary contact node");
    std::vector<double> lost_projection_state = problem.initial_state();
    lost_projection_state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
        sliding_global_node)] += 1.0e-2;
    bool lost_projection_rejected = false;
    try {
        problem.validate_state(lost_projection_state);
    } catch (const std::domain_error&) {
        lost_projection_rejected = true;
    }
    passed = check(lost_projection_rejected,
                 "M3.3 rejects a secondary node outside the complete "
                 "primary chain")
             && passed;
    const std::size_t secondary_index = static_cast<std::size_t>(
        std::find(secondary_sources.begin(), secondary_sources.end(), sliding_source) - secondary_sources.begin());
    const std::size_t primary_region = fuelsim::rz::ProblemAccess::region_index(problem, "cladding");
    const fuelsim::RegionBoundary primary_boundary = fuelsim::rz::ProblemAccess::region_mesh(problem, primary_region)
                                                         .map_side_set(source, definition.spatial.contacts[0].primary);
    std::vector<double> primary_positions;
    primary_positions.reserve(primary_boundary.nodes.size());
    for (const std::size_t node : primary_boundary.nodes)
        primary_positions.push_back(
            fuelsim::rz::ProblemAccess::region_mesh(problem, primary_region).nodes().at(node).z);
    std::sort(primary_positions.begin(), primary_positions.end());
    const auto initial_contact =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, problem.initial_state());
    passed = check(primary_positions.size() >= 5
                       && primary_positions.size() - 2 >= initial_contact.at(secondary_index).primary_segment + 2,
                 "projection test crosses beyond adjacent primary segments")
             && passed;
    const double target_position = 0.5 * (primary_positions[primary_positions.size() - 2] + primary_positions.back());
    std::vector<double> large_sliding_state = problem.initial_state();
    large_sliding_state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
        sliding_global_node)] = target_position - source.nodes().at(sliding_source).z;
    large_sliding_state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement,
        sliding_global_node)] = 3.0e-5;
    problem.validate_state(large_sliding_state);
    const std::vector<fuelsim::ContactNodeSummary> large_sliding_summary =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, large_sliding_state);
    passed = check(large_sliding_summary.at(secondary_index).projected
                       && large_sliding_summary.at(secondary_index).primary_segment == primary_positions.size() - 2
                       && large_sliding_summary.at(secondary_index).pressure > 0.0,
                 "M3.3 dynamically transfers a secondary node across more "
                 "than the former three-segment window")
             && passed;
    const double transfer_position = primary_positions[primary_positions.size() - 2];
    constexpr double transfer_offset = 1.0e-10;
    const auto transferred_node = [&](double position) {
        std::vector<double> state = problem.initial_state();
        state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
            sliding_global_node)] = position - source.nodes().at(sliding_source).z;
        state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement,
            sliding_global_node)] = 3.0e-5;
        problem.validate_state(state);
        return fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, state).at(secondary_index);
    };
    const fuelsim::ContactNodeSummary before_transfer = transferred_node(transfer_position - transfer_offset);
    const fuelsim::ContactNodeSummary at_transfer = transferred_node(transfer_position);
    const fuelsim::ContactNodeSummary after_transfer = transferred_node(transfer_position + transfer_offset);
    const double transfer_force_scale =
        std::max({1.0, std::abs(before_transfer.contact_force), std::abs(after_transfer.contact_force)});
    passed =
        check(before_transfer.primary_segment + 1 == after_transfer.primary_segment
                  && at_transfer.primary_segment == after_transfer.primary_segment,
            "M3.3 internal primary vertex has one owner and transfers "
            "ownership exactly once")
        && check(std::abs(before_transfer.contact_force - after_transfer.contact_force) < 1.0e-5 * transfer_force_scale,
            "M3.3 contact force is continuous across dynamic segment "
            "ownership transfer")
        && passed;
    problem.validate_state(large_sliding_state);
    double radial_contact_sum = 0.0;
    double axial_contact_sum = 0.0;
    double contact_force_scale = 0.0;
    for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
        contribution < problem.contribution_count();
        ++contribution) {
        if (fuelsim::rz::ProblemAccess::contribution_type(problem, contribution)
            != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        const fuelsim::LocalResidual local = fuelsim::rz::ProblemAccess::contribution_residual(problem,
            contribution,
            fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, large_sliding_state));
        for (std::size_t row = 4; row < 8; ++row) {
            radial_contact_sum += local[row];
            contact_force_scale += std::abs(local[row]);
        }
        for (std::size_t row = 8; row < 12; ++row) {
            axial_contact_sum += local[row];
            contact_force_scale += std::abs(local[row]);
        }
    }
    passed = check(std::abs(radial_contact_sum) < 1.0e-13 * (1.0 + contact_force_scale)
                       && std::abs(axial_contact_sum) < 1.0e-13 * (1.0 + contact_force_scale),
                 "M3.3 dynamically selected contact reactions remain "
                 "discretely conservative")
             && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2)
        return 2;
    try {
        return check_projection(argv[1]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
