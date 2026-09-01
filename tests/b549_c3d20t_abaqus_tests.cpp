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
    const std::string& material_path) {
    const auto temperature_reference = read_temperature(temperature_path);
    const auto displacement_reference = read_displacement(displacement_path);
    const auto material_reference = read_material(material_path);
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
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
    double maximum_material_coordinate_difference = 0.0;
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

    fuelsim::test::print_relative_metrics("b549_temperature", temperature_metrics);
    fuelsim::test::print_grouped_relative_metrics("b549_displacement_vector", displacement_metrics);
    fuelsim::test::print_relative_metrics("b549_equivalent_stress", stress_metrics);
    fuelsim::test::print_relative_metrics("b549_equivalent_plastic_strain", plastic_metrics);
    fuelsim::test::print_relative_metrics("b549_equivalent_creep_strain", creep_metrics);
    std::cout << "b549_reference_coordinate_maximum_difference=" << maximum_reference_coordinate_difference << '\n'
              << "b549_material_coordinate_maximum_difference=" << maximum_material_coordinate_difference << '\n'
              << "b549_accepted_steps=" << solve.accepted_steps.size() << '\n'
              << "b549_rejected_steps=" << solve.rejected_steps.size() << '\n'
              << "b549_nonlinear_iterations=" << solve.total_nonlinear_iterations << '\n'
              << "b549_dofs=" << problem.dof_count() << '\n';
    constexpr double tolerance = 1.0e-2;
    passed = check(maximum_reference_coordinate_difference < 1.0e-14,
                 "B5.49 uses identical tracked Fuelsim and Abaqus reference coordinates") &&
             passed;
    passed = check(maximum_material_coordinate_difference < 1.0e-6,
                 "B5.49 material-point association remains unambiguous in the current configuration") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(temperature_metrics, tolerance),
                 "B5.49 temperature metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::grouped_relative_metrics_below(displacement_metrics, tolerance) &&
                       displacement_metrics.maximum_zero_reference_difference < 1.0e-12,
                 "B5.49 displacement-vector metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(stress_metrics, tolerance),
                 "B5.49 equivalent-stress metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(plastic_metrics, tolerance) &&
                       plastic_metrics.nonzero_reference_count > 0 &&
                       plastic_metrics.maximum_zero_reference_difference < 1.0e-14,
                 "B5.49 active equivalent-plastic-strain metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(creep_metrics, tolerance) &&
                       creep_metrics.nonzero_reference_count > 0 &&
                       creep_metrics.maximum_zero_reference_difference < 1.0e-14,
                 "B5.49 active equivalent-creep-strain metrics are below one percent") &&
             passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_b549_c3d20t_abaqus_tests <case.fsi> <temperature.csv> <displacement.csv> "
                     "<material.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.49 C3D20T Abaqus comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3], argv[4]);
        if (passed && session.rank() == 0) std::cout << "[PASS] B5.49 C3D20T Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.49 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
