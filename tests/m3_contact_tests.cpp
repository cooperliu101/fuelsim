#include "fuelsim/case_input.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"
#include "support/moose_field_comparison.hpp"
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
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) result.push_back(field);
    return result;
}
std::size_t column(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE contact CSV is missing column: " + name);
    return static_cast<std::size_t>(found - header.begin());
}
std::vector<std::pair<double, double>> read_radial_pressure(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE contact CSV: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE contact CSV is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::size_t radius = column(header, "x");
    const std::size_t pressure = column(header, "contact_pressure");
    std::vector<std::pair<double, double>> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (radius >= values.size() || pressure >= values.size())
            throw std::invalid_argument("MOOSE contact CSV row is incomplete");
        result.emplace_back(std::stod(values[radius]), std::stod(values[pressure]));
    }
    std::sort(result.begin(), result.end());
    return result;
}
bool run_comparison(const std::string& input_path, const std::string& nodal_reference_path,
    const std::string& pressure_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::read_exodus_quad4(definition.mesh_file);
    bool passed =
        check(source.nodes().size() == 36 && source.elements().size() == 20,
            "M3.3 reads the complete tracked two-pellet MOOSE mesh") &&
        check(source.side_set("lower_top").sides.size() == 4 && source.side_set("upper_bottom").sides.size() == 6,
            "M3.3 keeps nonmatching four-to-six contact segmentation");
    fuelsim::SteadyProblem problem(definition.spatial_definition(), source);
    const std::vector<std::size_t> secondary_sources =
        fuelsim::rz::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    std::vector<std::size_t> ordered_secondary_sources = secondary_sources;
    std::sort(ordered_secondary_sources.begin(), ordered_secondary_sources.end(),
        [&source](std::size_t lhs, std::size_t rhs) { return source.nodes().at(lhs).r < source.nodes().at(rhs).r; });
    const std::size_t sliding_source = ordered_secondary_sources[ordered_secondary_sources.size() / 2];
    std::size_t sliding_global_node = problem.dof_count();
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const std::vector<std::size_t>& source_nodes =
            fuelsim::rz::ProblemAccess::region_mesh(problem, region).source_node_ids();
        const auto found = std::find(source_nodes.begin(), source_nodes.end(), sliding_source);
        if (found == source_nodes.end()) continue;
        sliding_global_node = fuelsim::rz::ProblemAccess::region_node_offset(problem, region) +
                              static_cast<std::size_t>(found - source_nodes.begin());
        break;
    }
    if (sliding_global_node == problem.dof_count())
        throw std::logic_error("M3.3 could not locate a secondary contact node");
    std::vector<double> lost_projection_state = problem.initial_state();
    lost_projection_state[fuelsim::rz::ProblemAccess::dof_map(problem).radial_displacement(sliding_global_node)] +=
        1.0e-2;
    bool lost_projection_rejected = false;
    try {
        problem.validate_state(lost_projection_state);
    } catch (const std::domain_error&) { lost_projection_rejected = true; }
    passed = check(lost_projection_rejected, "M3.3 rejects a secondary node outside the complete "
                                             "primary chain") &&
             passed;
    const std::size_t secondary_index = static_cast<std::size_t>(
        std::find(secondary_sources.begin(), secondary_sources.end(), sliding_source) - secondary_sources.begin());
    const std::size_t primary_region = fuelsim::rz::ProblemAccess::region_index(problem, "upper");
    const fuelsim::RegionBoundary primary_boundary = fuelsim::rz::ProblemAccess::region_mesh(problem, primary_region)
                                                         .map_side_set(source, definition.contacts[0].primary);
    std::vector<double> primary_radii;
    primary_radii.reserve(primary_boundary.nodes.size());
    for (const std::size_t node : primary_boundary.nodes)
        primary_radii.push_back(fuelsim::rz::ProblemAccess::region_mesh(problem, primary_region).nodes().at(node).r);
    std::sort(primary_radii.begin(), primary_radii.end());
    const double target_radius = 0.5 * (primary_radii[primary_radii.size() - 2] + primary_radii.back());
    std::vector<double> large_sliding_state = problem.initial_state();
    large_sliding_state[fuelsim::rz::ProblemAccess::dof_map(problem).radial_displacement(sliding_global_node)] =
        target_radius - source.nodes().at(sliding_source).r;
    large_sliding_state[fuelsim::rz::ProblemAccess::dof_map(problem).axial_displacement(sliding_global_node)] = 3.0e-6;
    problem.validate_state(large_sliding_state);
    const std::vector<fuelsim::ContactNodeSummary> large_sliding_summary =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, large_sliding_state);
    passed = check(large_sliding_summary.at(secondary_index).projected &&
                       large_sliding_summary.at(secondary_index).primary_segment == primary_radii.size() - 2 &&
                       large_sliding_summary.at(secondary_index).pressure > 0.0,
                 "M3.3 dynamically transfers a secondary node across more "
                 "than the former three-segment window") &&
             passed;
    const double transfer_radius = primary_radii[primary_radii.size() - 2];
    constexpr double transfer_offset = 1.0e-10;
    const auto transferred_node = [&](double radius) {
        std::vector<double> state = problem.initial_state();
        state[fuelsim::rz::ProblemAccess::dof_map(problem).radial_displacement(sliding_global_node)] =
            radius - source.nodes().at(sliding_source).r;
        state[fuelsim::rz::ProblemAccess::dof_map(problem).axial_displacement(sliding_global_node)] = 3.0e-6;
        problem.validate_state(state);
        return fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, state).at(secondary_index);
    };
    const fuelsim::ContactNodeSummary before_transfer = transferred_node(transfer_radius - transfer_offset);
    const fuelsim::ContactNodeSummary at_transfer = transferred_node(transfer_radius);
    const fuelsim::ContactNodeSummary after_transfer = transferred_node(transfer_radius + transfer_offset);
    const double transfer_force_scale =
        std::max({1.0, std::abs(before_transfer.contact_force), std::abs(after_transfer.contact_force)});
    passed =
        check(before_transfer.primary_segment + 1 == after_transfer.primary_segment &&
                  at_transfer.primary_segment == after_transfer.primary_segment,
            "M3.3 internal primary vertex has one owner and transfers "
            "ownership exactly once") &&
        check(std::abs(before_transfer.contact_force - after_transfer.contact_force) < 1.0e-5 * transfer_force_scale,
            "M3.3 contact force is continuous across dynamic segment "
            "ownership transfer") &&
        passed;
    problem.validate_state(large_sliding_state);
    double radial_contact_sum = 0.0;
    double axial_contact_sum = 0.0;
    double contact_force_scale = 0.0;
    for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
        contribution < problem.contribution_count(); ++contribution) {
        if (fuelsim::rz::ProblemAccess::contribution_type(problem, contribution) !=
            fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        const fuelsim::LocalResidual local = fuelsim::rz::ProblemAccess::contribution_residual(problem, contribution,
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
    passed = check(std::abs(radial_contact_sum) < 1.0e-13 * (1.0 + contact_force_scale) &&
                       std::abs(axial_contact_sum) < 1.0e-13 * (1.0 + contact_force_scale),
                 "M3.3 dynamically selected contact reactions remain "
                 "discretely conservative") &&
             passed;
    const fuelsim::SolverOptions options = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks, definition.steady_execution.minimum_load_increment},
        options);
    passed = check(result.completed && result.solve.converged, "M3.3 two-pellet contact solve converges") && passed;
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, result.solve.state, reference);
    constexpr double tolerance = 1.0e-2;
    passed = check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                 "M3.3 compares every MOOSE node at matching coordinates") &&
             check(fuelsim::test::relative_metrics_below(fields.temperature, tolerance),
                 "M3.3 temperature three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(fields.radial_displacement, tolerance),
                 "M3.3 radial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(fields.axial_displacement, tolerance),
                 "M3.3 axial-displacement three full-field errors pass") &&
             passed;
    const std::vector<fuelsim::ContactNodeSummary> actual =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, result.solve.state);
    const std::vector<std::pair<double, double>> pressure_reference = read_radial_pressure(pressure_reference_path);
    if (actual.size() != pressure_reference.size())
        throw std::invalid_argument("M3.3 fuelsim and MOOSE contact node counts differ");
    fuelsim::test::FieldErrorMetrics pressure;
    for (std::size_t node = 0; node < actual.size(); ++node) {
        if (std::abs(actual[node].r - pressure_reference[node].first) > 1.0e-12)
            throw std::invalid_argument("M3.3 fuelsim and MOOSE contact radii differ");
        pressure.add(actual[node].pressure, pressure_reference[node].second);
    }
    passed = check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                 "M3.3 contact-pressure three full-field errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m33_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m33_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m33_axial_displacement", fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m33_contact_pressure", pressure);
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m3_contact_tests <case.fsi> "
                     "<all-nodes.csv> <contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M3.3 contact comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3])) return 1;
        std::cout << "[PASS] M3.3 nonmatching two-pellet MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.3 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
