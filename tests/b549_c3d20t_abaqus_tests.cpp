#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct TemperatureReference final {
    std::size_t node;
    fuelsim::CartesianPoint3 point;
    double value;
};

struct DisplacementReference final {
    std::size_t node;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> value;
};

struct MaterialReference final {
    std::size_t element, point;
    fuelsim::CartesianPoint3 position;
    double volume, equivalent_stress, equivalent_plastic_strain, equivalent_creep_strain;
};

struct ContactReference final {
    std::size_t node;
    fuelsim::CartesianPoint3 point;
    double gap, pressure;
    std::array<double, 3> normal_force, tangential_force, tangential_slip;
    bool has_slip;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete B5.49 CSV row: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B5.49 CSV number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double value = number(values, column, path);
    if (value < 1.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        value != std::floor(value))
        throw std::invalid_argument("Invalid B5.49 CSV index: " + path);
    return static_cast<std::size_t>(value);
}

std::vector<TemperatureReference> read_temperature(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B5.49 temperature reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "id,x,y,z,temperature")
        throw std::invalid_argument("Unexpected B5.49 temperature header: " + path);
    std::vector<TemperatureReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, path) - 1,
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}, number(values, 4, path)});
    }
    if (result.size() != 40) throw std::invalid_argument("B5.49 requires 40 corner-temperature rows");
    return result;
}

std::vector<DisplacementReference> read_displacement(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B5.49 displacement reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "id,x,y,z,displacement_x,displacement_y,displacement_z")
        throw std::invalid_argument("Unexpected B5.49 displacement header: " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, path) - 1,
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}});
    }
    if (result.size() != 112) throw std::invalid_argument("B5.49 requires 112 displacement rows");
    return result;
}

std::vector<MaterialReference> read_material(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B5.49 material reference: " + path);
    std::string line;
    const std::string expected = "element,integration_point,current_x,current_y,current_z,ivol,vonmises_stress,"
                                 "effective_plastic_strain,effective_creep_strain";
    if (!std::getline(input, line) || line != expected)
        throw std::invalid_argument("Unexpected B5.49 material header: " + path);
    std::vector<MaterialReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, path), index_value(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)}, number(values, 5, path),
            number(values, 6, path), number(values, 7, path), number(values, 8, path)});
    }
    if (result.size() != 216) throw std::invalid_argument("B5.49 requires 216 material-point rows");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read C3D20T contact reference: " + path);
    std::string line;
    const std::string expected = "id,x,y,z,gap,pressure,normal_x,normal_y,normal_z,shear_x,shear_y,shear_z";
    const std::string friction_expected = expected + ",slip_1,slip_2,tangent_1_x,tangent_1_y,tangent_1_z,tangent_2_x,"
                                                     "tangent_2_y,tangent_2_z";
    if (!std::getline(input, line) || (line != expected && line != friction_expected))
        throw std::invalid_argument("Unexpected C3D20T contact header: " + path);
    const bool has_slip = line == friction_expected;
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        ContactReference reference{index_value(values, 0, path) - 1,
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}, number(values, 4, path),
            number(values, 5, path), {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
            {number(values, 9, path), number(values, 10, path), number(values, 11, path)}, {}, has_slip};
        if (has_slip) {
            const double slip_first = number(values, 12, path), slip_second = number(values, 13, path);
            for (std::size_t component = 0; component < 3; ++component)
                reference.tangential_slip[component] = slip_first * number(values, 14 + component, path) +
                                                       slip_second * number(values, 17 + component, path);
        } else
            for (const double force : reference.tangential_force)
                if (force != 0.0)
                    throw std::invalid_argument("C3D20T frictionless contact reference contains shear force: " + path);
        result.push_back(reference);
    }
    if (result.size() != 8) throw std::invalid_argument("C3D20T contact reference requires eight secondary nodes");
    return result;
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = definition.solver.absolute_tolerance;
    result.relative_tolerance = definition.solver.relative_tolerance;
    result.step_tolerance = definition.solver.step_tolerance;
    result.maximum_iterations = definition.solver.maximum_iterations;
    result.linear_solver = definition.solver.linear_solver;
    result.preconditioner = definition.solver.preconditioner;
    result.direct_factorization = definition.solver.direct_factorization;
    result.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    result.jacobian_lag = definition.solver.jacobian_lag;
    result.line_search = definition.solver.line_search;
    result.backtracking_fallback = definition.solver.backtracking_fallback;
    result.field_residual_scaling = definition.solver.field_residual_scaling;
    result.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
    result.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
    result.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
    return result;
}

fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    const auto& input = definition.transient_execution;
    return {input.end_time, input.initial_time_step, input.minimum_time_step, input.maximum_time_step,
        input.growth_factor, input.cutback_factor, input.maximum_cutbacks_per_step, input.load_ramp_time,
        input.target_nonlinear_iterations, input.iteration_window, input.time_error_relative_tolerance,
        input.temperature_time_absolute_tolerance, input.displacement_time_absolute_tolerance,
        input.time_error_safety_factor, input.strain_history_time_absolute_tolerance,
        input.stress_history_time_absolute_tolerance, input.include_thermal_time_term, input.use_linear_time_predictor};
}

double mechanical_contact_directional_error(
    fuelsim::TransientProblem& problem, const std::vector<double>& global_state, double perturbation) {
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
        if (!(magnitude > 0.0)) continue;
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
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
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
        if (reference_squared > 0.0) {
            const double error = std::sqrt(difference_squared / reference_squared);
            maximum_error = std::max(maximum_error, error);
        }
    }
    spatial.validate_state(global_state);
    return maximum_error;
}

double distance(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return std::sqrt(
        std::pow(first.x - second.x, 2) + std::pow(first.y - second.y, 2) + std::pow(first.z - second.z, 2));
}

double equivalent_stress(const fuelsim::SymmetricTensor3Values& stress) {
    const double mean = (stress.xx + stress.yy + stress.zz) / 3.0;
    return std::sqrt(
        1.5 * (std::pow(stress.xx - mean, 2) + std::pow(stress.yy - mean, 2) + std::pow(stress.zz - mean, 2) +
                  2.0 * (stress.xy * stress.xy + stress.yz * stress.yz + stress.xz * stress.xz)));
}

bool run(const std::string& case_path, const std::string& temperature_path, const std::string& displacement_path,
    const std::string& material_path, const std::string& contact_path) {
    const auto temperature_reference = read_temperature(temperature_path);
    const auto displacement_reference = read_displacement(displacement_path);
    const auto material_reference = read_material(material_path);
    const std::vector<ContactReference> contact_reference =
        contact_path.empty() ? std::vector<ContactReference>{} : read_contact(contact_path);
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const bool frictional =
        !definition.spatial.contacts.empty() && definition.spatial.contacts.front().friction_coefficient > 0.0;
    const std::string prefix = frictional ? "b555_" : contact_reference.empty() ? "b549_" : "b550_";
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    const fuelsim::TransientResult solve =
        fuelsim::solve_transient(problem, time_options(definition), solver_options(definition));
    bool passed = check(solve.completed && solve.accepted_steps.size() == 20 && solve.rejected_steps.empty(),
        "B5.49 completes twenty fixed increments without reducing the time step");
    if (!solve.completed) return false;
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::vector<double>& state = problem.committed_solution();

    std::map<std::size_t, std::size_t> displacement_nodes, temperature_nodes;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const fuelsim::Hex20RegionMesh& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            displacement_nodes.emplace(source, spatial.global_node(region, local));
            if (region_mesh.temperature_nodes()[local])
                temperature_nodes.emplace(source, spatial.global_temperature_node(region, local));
        }
    }
    if (displacement_nodes.size() != 112 || temperature_nodes.size() != 40)
        throw std::logic_error("B5.49 mixed-order source-node maps have changed");

    fuelsim::test::FieldErrorMetrics temperature_metrics;
    fuelsim::test::GroupedFieldErrorMetrics displacement_metrics;
    double maximum_reference_coordinate_difference = 0.0;
    for (const auto& reference : temperature_reference) {
        maximum_reference_coordinate_difference = std::max(
            maximum_reference_coordinate_difference, distance(mesh.nodes().at(reference.node), reference.point));
        temperature_metrics.add(
            state[spatial.dof(fuelsim::Field::temperature, temperature_nodes.at(reference.node))], reference.value);
    }
    for (const auto& reference : displacement_reference) {
        maximum_reference_coordinate_difference = std::max(
            maximum_reference_coordinate_difference, distance(mesh.nodes().at(reference.node), reference.point));
        std::array<double, 3> actual{};
        const std::size_t global = displacement_nodes.at(reference.node);
        for (std::size_t component = 0; component < 3; ++component) {
            const fuelsim::Field field = component == 0   ? fuelsim::Field::displacement_x
                                         : component == 1 ? fuelsim::Field::displacement_y
                                                          : fuelsim::Field::displacement_z;
            actual[component] = state[spatial.dof(field, global)];
        }
        displacement_metrics.add(actual.data(), reference.value.data(), actual.size());
    }

    std::map<std::size_t, std::vector<MaterialReference>> references_by_element;
    for (const auto& reference : material_reference) references_by_element[reference.element].push_back(reference);
    fuelsim::test::FieldErrorMetrics stress_metrics, plastic_metrics, creep_metrics;
    double maximum_material_coordinate_difference = 0.0,
           minimum_material_second_to_first_distance_ratio = std::numeric_limits<double>::infinity();
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const fuelsim::Hex20RegionMesh& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const std::size_t source_element = region_mesh.source_element_ids()[element] + 1;
            const auto& references = references_by_element.at(source_element);
            const auto& mesh_element = region_mesh.elements()[element];
            const fuelsim::Hex20Geometry& geometry = spatial.hex20_region_element_geometry(region, element);
            const fuelsim::CartesianMaterialHistory& history =
                fuelsim::cartesian::ProblemAccess::material_history(problem, region, element);
            std::array<fuelsim::CartesianPoint3, 27> positions{};
            for (std::size_t q = 0; q < positions.size(); ++q) {
                positions[q] = geometry.mechanical_points[q].position;
                for (std::size_t local = 0; local < 20; ++local) {
                    const std::size_t global = spatial.global_node(region, mesh_element.nodes[local]);
                    const double shape = geometry.mechanical_points[q].displacement_shape[local];
                    positions[q].x += shape * state[spatial.dof(fuelsim::Field::displacement_x, global)];
                    positions[q].y += shape * state[spatial.dof(fuelsim::Field::displacement_y, global)];
                    positions[q].z += shape * state[spatial.dof(fuelsim::Field::displacement_z, global)];
                }
            }
            std::vector<std::tuple<double, std::size_t, std::size_t>> candidates;
            for (std::size_t actual = 0; actual < 27; ++actual)
                for (std::size_t reference = 0; reference < 27; ++reference)
                    candidates.emplace_back(
                        distance(positions[actual], references[reference].position), actual, reference);
            for (std::size_t actual = 0; actual < 27; ++actual) {
                std::array<double, 27> distances{};
                for (std::size_t reference = 0; reference < 27; ++reference)
                    distances[reference] = distance(positions[actual], references[reference].position);
                std::sort(distances.begin(), distances.end());
                if (distances[0] > 0.0)
                    minimum_material_second_to_first_distance_ratio =
                        std::min(minimum_material_second_to_first_distance_ratio, distances[1] / distances[0]);
            }
            std::sort(candidates.begin(), candidates.end());
            std::array<bool, 27> actual_used{}, reference_used{};
            std::size_t pair_count = 0;
            for (const auto& candidate : candidates) {
                const std::size_t actual = std::get<1>(candidate), reference = std::get<2>(candidate);
                if (actual_used[actual] || reference_used[reference]) continue;
                actual_used[actual] = reference_used[reference] = true;
                ++pair_count;
                maximum_material_coordinate_difference =
                    std::max(maximum_material_coordinate_difference, std::get<0>(candidate));
                stress_metrics.add(equivalent_stress(history[actual].stress), references[reference].equivalent_stress);
                plastic_metrics.add(
                    history[actual].equivalent_plastic_strain, references[reference].equivalent_plastic_strain);
                creep_metrics.add(
                    history[actual].equivalent_creep_strain, references[reference].equivalent_creep_strain);
            }
            if (pair_count != 27) throw std::logic_error("B5.49 material-point association is incomplete");
        }
    }

    fuelsim::test::print_relative_metrics(prefix + "temperature", temperature_metrics);
    fuelsim::test::print_grouped_relative_metrics(prefix + "displacement_vector", displacement_metrics);
    fuelsim::test::print_relative_metrics(prefix + "equivalent_stress", stress_metrics);
    fuelsim::test::print_relative_metrics(prefix + "equivalent_plastic_strain", plastic_metrics);
    fuelsim::test::print_relative_metrics(prefix + "equivalent_creep_strain", creep_metrics);
    double abaqus_contact_force = 0.0, fuelsim_contact_force = 0.0;
    double abaqus_nodal_tangential_force_magnitude_sum = 0.0, fuelsim_integrated_tangential_traction_magnitude = 0.0;
    double contact_jacobian_directional_error = 0.0;
    std::size_t active_contact_nodes = 0, sticking_contact_nodes = 0, sliding_contact_nodes = 0;
    fuelsim::test::FieldErrorMetrics contact_gap_metrics, constraint_contact_pressure_metrics,
        recovered_contact_pressure_metrics;
    fuelsim::test::GroupedFieldErrorMetrics contact_normal_force_metrics, contact_tangential_force_metrics,
        contact_tangential_slip_metrics, contact_tangential_resultant_metrics;
    if (!contact_reference.empty()) {
        std::map<std::size_t, ContactReference> contact_by_node;
        std::array<double, 3> resultant{}, tangential_resultant{};
        for (const ContactReference& reference : contact_reference) {
            if (reference.has_slip != frictional)
                throw std::invalid_argument("C3D20T contact reference friction fields do not match the input");
            contact_by_node.emplace(reference.node, reference);
            maximum_reference_coordinate_difference = std::max(
                maximum_reference_coordinate_difference, distance(mesh.nodes().at(reference.node), reference.point));
            for (std::size_t component = 0; component < 3; ++component) {
                resultant[component] += reference.normal_force[component];
                tangential_resultant[component] += reference.tangential_force[component];
            }
            abaqus_nodal_tangential_force_magnitude_sum +=
                std::hypot(reference.tangential_force[0], reference.tangential_force[1], reference.tangential_force[2]);
        }
        abaqus_contact_force = std::hypot(resultant[0], resultant[1], resultant[2]);
        const auto summaries = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        if (summaries.size() != source_nodes.size())
            throw std::logic_error("C3D20T contact summary and source-node maps have different sizes");
        for (std::size_t node = 0; node < summaries.size(); ++node) {
            const auto& summary = summaries[node];
            const ContactReference& reference = contact_by_node.at(source_nodes[node]);
            if (summary.pressure > 0.0) {
                ++active_contact_nodes;
                ++(summary.sliding ? sliding_contact_nodes : sticking_contact_nodes);
            }
            const double penalty = definition.spatial.contacts.at(0).penalty;
            contact_gap_metrics.add(summary.gap, reference.gap);
            constraint_contact_pressure_metrics.add(
                std::max(-penalty * summary.gap, 0.0), std::max(-penalty * reference.gap, 0.0));
            recovered_contact_pressure_metrics.add(summary.pressure, reference.pressure);
            std::array<double, 3> fuelsim_secondary_force{}, fuelsim_secondary_tangential_force{};
            for (std::size_t component = 0; component < fuelsim_secondary_force.size(); ++component)
                fuelsim_secondary_force[component] = -summary.normal_contact_force[component];
            for (std::size_t component = 0; component < fuelsim_secondary_tangential_force.size(); ++component)
                fuelsim_secondary_tangential_force[component] = -summary.tangential_contact_force[component];
            contact_normal_force_metrics.add(
                fuelsim_secondary_force.data(), reference.normal_force.data(), reference.normal_force.size());
            if (frictional) {
                contact_tangential_force_metrics.add(fuelsim_secondary_tangential_force.data(),
                    reference.tangential_force.data(), reference.tangential_force.size());
                contact_tangential_slip_metrics.add(
                    summary.tangential_slip.data(), reference.tangential_slip.data(), reference.tangential_slip.size());
            }
        }
        const fuelsim::InterfaceSummary interface =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, state);
        fuelsim_contact_force = interface.total_contact_force;
        fuelsim_integrated_tangential_traction_magnitude = interface.total_tangential_force;
        std::array<double, 3> fuelsim_tangential_resultant{};
        for (const auto& summary : summaries)
            for (std::size_t component = 0; component < fuelsim_tangential_resultant.size(); ++component)
                fuelsim_tangential_resultant[component] -= summary.tangential_contact_force[component];
        contact_tangential_resultant_metrics.add(
            fuelsim_tangential_resultant.data(), tangential_resultant.data(), tangential_resultant.size());
        contact_jacobian_directional_error =
            mechanical_contact_directional_error(problem, state, frictional ? 1.0e-9 : 1.0e-8);
        const double contact_tolerance = frictional ? 5.0e-3 : 1.0e-2;
        passed = check(active_contact_nodes > 0 && interface.active_contact_nodes == active_contact_nodes,
                     "C3D20T contact comparison has active Fuelsim contact constraints") &&
                 passed;
        passed =
            check(abaqus_contact_force > 0.0 &&
                      std::abs(fuelsim_contact_force - abaqus_contact_force) / abaqus_contact_force < contact_tolerance,
                "C3D20T total contact force passes its Abaqus tolerance") &&
            passed;
        passed = check(fuelsim::test::grouped_relative_metrics_below(contact_normal_force_metrics, contact_tolerance),
                     "C3D20T secondary nodal normal-force metrics pass their Abaqus tolerance") &&
                 passed;
        constexpr double pressure_tolerance = 5.0e-3;
        passed = check(fuelsim::test::relative_metrics_below(constraint_contact_pressure_metrics, pressure_tolerance),
                     "C3D20T node-centered constraint-pressure metrics are below 0.5 percent") &&
                 passed;
        passed = check(fuelsim::test::relative_metrics_below(recovered_contact_pressure_metrics, pressure_tolerance),
                     "C3D20T recovered nodal contact-pressure metrics are below 0.5 percent") &&
                 passed;
        if (frictional) {
            passed = check(abaqus_nodal_tangential_force_magnitude_sum > 0.0 &&
                               fuelsim_integrated_tangential_traction_magnitude > 0.0 && sliding_contact_nodes > 0,
                         "C3D20T friction comparison activates nonzero tangential force and sliding") &&
                     passed;
            passed =
                check(fuelsim::test::relative_metrics_below(contact_gap_metrics, contact_tolerance) &&
                          fuelsim::test::grouped_relative_metrics_below(
                              contact_tangential_force_metrics, contact_tolerance) &&
                          fuelsim::test::grouped_relative_metrics_below(
                              contact_tangential_slip_metrics, contact_tolerance) &&
                          fuelsim::test::grouped_relative_metrics_below(
                              contact_tangential_resultant_metrics, contact_tolerance),
                    "C3D20T gap, tangential-force, tangential-slip, and resultant metrics are below 0.5 percent") &&
                passed;
        }
        passed = check(contact_jacobian_directional_error < 2.0e-5,
                     "C3D20T finite-sliding contact Jacobian matches a centered directional difference") &&
                 passed;
        fuelsim::test::print_relative_metrics(prefix + "contact_gap", contact_gap_metrics);
        fuelsim::test::print_relative_metrics(
            prefix + "constraint_contact_pressure", constraint_contact_pressure_metrics);
        fuelsim::test::print_relative_metrics(
            prefix + "recovered_contact_pressure", recovered_contact_pressure_metrics);
        fuelsim::test::print_grouped_relative_metrics(prefix + "contact_normal_force", contact_normal_force_metrics);
        if (frictional) {
            fuelsim::test::print_grouped_relative_metrics(
                prefix + "contact_tangential_force", contact_tangential_force_metrics);
            fuelsim::test::print_grouped_relative_metrics(
                prefix + "contact_tangential_slip", contact_tangential_slip_metrics);
            fuelsim::test::print_grouped_relative_metrics(
                prefix + "contact_tangential_resultant", contact_tangential_resultant_metrics);
        }
    }
    std::cout << prefix << "reference_coordinate_maximum_difference=" << maximum_reference_coordinate_difference << '\n'
              << prefix << "material_coordinate_maximum_difference=" << maximum_material_coordinate_difference << '\n'
              << prefix
              << "material_minimum_second_to_first_distance_ratio=" << minimum_material_second_to_first_distance_ratio
              << '\n'
              << prefix << "active_contact_nodes=" << active_contact_nodes << '\n'
              << prefix << "sticking_contact_nodes=" << sticking_contact_nodes << '\n'
              << prefix << "sliding_contact_nodes=" << sliding_contact_nodes << '\n'
              << prefix << "fuelsim_total_contact_force=" << fuelsim_contact_force << '\n'
              << prefix << "abaqus_total_contact_force=" << abaqus_contact_force << '\n'
              << prefix
              << "fuelsim_integrated_tangential_traction_magnitude=" << fuelsim_integrated_tangential_traction_magnitude
              << '\n'
              << prefix << "abaqus_nodal_tangential_force_magnitude_sum=" << abaqus_nodal_tangential_force_magnitude_sum
              << '\n'
              << prefix << "contact_jacobian_directional_error=" << contact_jacobian_directional_error << '\n'
              << prefix << "accepted_steps=" << solve.accepted_steps.size() << '\n'
              << prefix << "rejected_steps=" << solve.rejected_steps.size() << '\n'
              << prefix << "nonlinear_iterations=" << solve.total_nonlinear_iterations << '\n'
              << prefix << "dofs=" << problem.dof_count() << '\n';
    const double tolerance = frictional ? 5.0e-3 : 1.0e-2;
    passed = check(maximum_reference_coordinate_difference < 1.0e-14,
                 "B5.49 uses identical tracked Fuelsim and Abaqus reference coordinates") &&
             passed;
    passed = check(minimum_material_second_to_first_distance_ratio > 10.0,
                 "B5.49 material-point association remains unambiguous in the current configuration") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(temperature_metrics, tolerance),
                 "C3D20T temperature metrics pass their Abaqus tolerance") &&
             passed;
    passed = check(fuelsim::test::grouped_relative_metrics_below(displacement_metrics, tolerance) &&
                       displacement_metrics.maximum_zero_reference_difference < 1.0e-12,
                 "C3D20T displacement-vector metrics pass their Abaqus tolerance") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(stress_metrics, tolerance),
                 "C3D20T equivalent-stress metrics pass their Abaqus tolerance") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(plastic_metrics, tolerance) &&
                       plastic_metrics.nonzero_reference_count > 0 &&
                       plastic_metrics.maximum_zero_reference_difference < 1.0e-14,
                 "C3D20T active equivalent-plastic-strain metrics pass their Abaqus tolerance") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(creep_metrics, tolerance) &&
                       creep_metrics.nonzero_reference_count > 0 &&
                       creep_metrics.maximum_zero_reference_difference < 1.0e-14,
                 "C3D20T active equivalent-creep-strain metrics pass their Abaqus tolerance") &&
             passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5 && argc != 6) {
        std::cerr << "Usage: fuelsim_b549_c3d20t_abaqus_tests <case.fsi> <temperature.csv> <displacement.csv> "
                     "<material.csv> [contact.csv]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.49 C3D20T Abaqus comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3], argv[4], argc == 6 ? argv[5] : "");
        if (passed && session.rank() == 0) std::cout << "[PASS] B5.49 C3D20T Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.49 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
