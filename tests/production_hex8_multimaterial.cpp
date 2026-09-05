#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
struct Point {
    double x, y, z;
};

struct Tensor {
    double xx, yy, zz, xy, yz, xz;
};

constexpr double length = 2.0, half_thickness = 0.1, width = 0.25, step_time = 1.0e7;

struct CaseSpec final {
    const char* name;
    std::size_t nx, subdivisions_per_layer, nz;
    bool distorted;
};

constexpr std::array<CaseSpec, 3> cases = {
    CaseSpec{"coarse", 2, 1, 1, false}, CaseSpec{"refined", 4, 2, 2, false}, CaseSpec{"distorted", 4, 2, 2, true}};

struct NodeReference final {
    std::string case_name;
    std::size_t stage, node;
    double time;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::string case_name;
    std::size_t stage, element;
    double time;
    Point position;
    double temperature;
    std::array<double, 3> heat_flux;
    Tensor stress, strain;
    double integration_volume;
};

struct BendingSummary {
    double peak, cooled;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.9 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.9 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

double canonical_zero(double value, double tolerance) { return std::abs(value) < tolerance ? 0.0 : value; }

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.9 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "case,stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.9 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 12) throw std::invalid_argument("Unexpected Abaqus B5.9 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 4, path);
            const double tolerance = field >= 5 ? 1.0e-5 : (field == 4 ? 1.0e-8 : 1.0e-16);
            if (field != 0) fields[field] = canonical_zero(fields[field], tolerance);
        }
        result.push_back({values[0], positive_integer(number(values, 1, path), path),
            positive_integer(number(values, 3, path), path), number(values, 2, path), fields});
    }
    if (result.size() != 672) throw std::invalid_argument("Abaqus B5.9 nodal reference must contain 672 rows");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.9 integration points: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "case,stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,"
                "hfl_z_w_m2,s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,"
                "e13_engineering,e23_engineering,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.9 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 25)
            throw std::invalid_argument("Unexpected Abaqus B5.9 integration-point column count in " + path);
        std::array<double, 3> heat_flux = {number(values, 9, path), number(values, 10, path), number(values, 11, path)};
        for (double& value : heat_flux) value = canonical_zero(value, 1.0e-8);
        Tensor stress = {number(values, 12, path), number(values, 13, path), number(values, 14, path),
            number(values, 15, path), number(values, 17, path), number(values, 16, path)};
        Tensor strain = {number(values, 18, path), number(values, 19, path), number(values, 20, path),
            0.5 * number(values, 21, path), 0.5 * number(values, 23, path), 0.5 * number(values, 22, path)};
        for (double* value : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *value = canonical_zero(*value, 1.0e-4);
        for (double* value : {&strain.xx, &strain.yy, &strain.zz, &strain.xy, &strain.yz, &strain.xz})
            *value = canonical_zero(*value, 1.0e-15);
        result.push_back({values[0], positive_integer(number(values, 1, path), path),
            positive_integer(number(values, 3, path), path), number(values, 2, path),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)}, number(values, 8, path),
            heat_flux, stress, strain, number(values, 24, path)});
    }
    if (result.size() != 2176)
        throw std::invalid_argument("Abaqus B5.9 integration-point reference must contain 2176 rows");
    return result;
}

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance) {
    if (metrics.has_relative_norm() && !fuelsim::test::relative_metrics_below(metrics, relative_tolerance))
        return false;
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

void print_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    if (metrics.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metrics);
    else
        fuelsim::test::print_absolute_metrics(name, metrics);
}

BendingSummary run_case(const CaseSpec& spec, const std::string& output_path,
    const std::vector<NodeReference>& node_reference, const std::vector<IntegrationReference>& integration_reference,
    bool& passed) {
    const auto history = fuelsim::test::read_exodus_history(output_path);
    if (history.size() != 5) throw std::invalid_argument("B5.9 requires four output increments");
    const auto& final = history.back();
    const std::size_t node_count = (spec.nx + 1) * (2 * spec.subdivisions_per_layer + 1) * (spec.nz + 1);
    const std::size_t element_count = spec.nx * 2 * spec.subdivisions_per_layer * spec.nz;
    if (final.nodes.size() != node_count ||
        final.block_element_counts != std::vector<std::size_t>{element_count / 2, element_count / 2})
        throw std::invalid_argument("B5.9 output mesh must preserve all shared nodes and both material blocks");
    std::set<std::pair<std::size_t, std::size_t>> seen_nodes;
    std::set<std::tuple<std::size_t, std::size_t, std::size_t>> seen_points;
    std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
    const std::array<std::string, 8> fields = {"temperature", "displacement_x", "displacement_y", "displacement_z",
        "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
    std::size_t compared_nodes = 0;
    for (std::size_t stage = 1; stage <= 4; ++stage) {
        const auto& snapshot = history.at(stage);
        if (std::abs(snapshot.time - static_cast<double>(stage) * step_time) > 1e-6)
            throw std::invalid_argument("B5.9 output time differs");
        for (const NodeReference& reference : node_reference) {
            if (reference.case_name != spec.name || reference.stage != stage) continue;
            if (std::abs(reference.time - snapshot.time) > 1.0e-6 || reference.node < 1 ||
                reference.node > node_count || !seen_nodes.emplace(stage, reference.node).second)
                throw std::invalid_argument("Abaqus B5.9 nodal label or time lies outside the Fuelsim case");
            for (std::size_t field = 0; field < fields.size(); ++field)
                nodal_metrics[field].add(snapshot.nodal(fields[field]).at(reference.node - 1), reference.fields[field]);
            ++compared_nodes;
        }
        std::cout << "b59_" << spec.name << "_stage_" << stage
                  << "_relative_thermal_balance=" << snapshot.global("conservation_relative_thermal_balance")
                  << " stored_heat_rate=" << snapshot.global("conservation_stored_heat_rate")
                  << " dirichlet_heat_input_rate=" << snapshot.global("conservation_dirichlet_heat_input_rate")
                  << " global_thermal_balance=" << snapshot.global("conservation_global_thermal_balance")
                  << " unconstrained_thermal_residual_l2="
                  << snapshot.global("conservation_unconstrained_thermal_residual_l2") << '\n';
        const double thermal_scale = std::abs(snapshot.global("conservation_stored_heat_rate")) +
                                     std::abs(snapshot.global("conservation_dirichlet_heat_input_rate"));
        const bool relative_balance_required = thermal_scale >= 1.0e-3;
        const bool thermal_balance_passed =
            std::abs(snapshot.global("conservation_global_thermal_balance")) < 1.0e-10 &&
            snapshot.global("conservation_unconstrained_thermal_residual_l2") < 1.0e-10 &&
            (!relative_balance_required || snapshot.global("conservation_relative_thermal_balance") < 1.0e-10);
        passed =
            check(thermal_balance_passed,
                std::string("B5.9 ") + spec.name + " thermal conservation closes at stage " + std::to_string(stage) +
                    (relative_balance_required ? " by relative and absolute metrics"
                                               : " by the separately gated near-zero absolute metric")) &&
            passed;
    }
    passed = check(compared_nodes == 4 * node_count,
                 std::string("B5.9 ") + spec.name + " compares every node at every stage") &&
             passed;

    const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y", "displacement_z",
        "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
    const std::array<double, 8> zero_tolerances = {1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 2.0e-4, 5.0e-2, 5.0e-2, 5.0e-2};
    for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
        print_metrics(std::string("b59_") + spec.name + '_' + nodal_names[field], nodal_metrics[field]);
        passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, zero_tolerances[field]),
                     std::string("B5.9 ") + spec.name + ' ' + nodal_names[field] +
                         " metrics are below the acceptance limits") &&
                 passed;
    }

    std::array<fuelsim::test::FieldErrorMetrics, 17> integration_metrics;
    std::size_t compared_points = 0;
    double maximum_coordinate_error = 0.0;
    for (const IntegrationReference& reference : integration_reference) {
        if (reference.case_name != spec.name) continue;
        if (reference.element < 1 || reference.element > element_count)
            throw std::invalid_argument("Abaqus B5.9 element label lies outside the Fuelsim mesh");
        const auto& snapshot = history.at(reference.stage);
        if (std::abs(reference.time - snapshot.time) > 1.0e-6)
            throw std::invalid_argument("Abaqus B5.9 integration-point time does not match the Fuelsim snapshot");
        const std::size_t element = reference.element - 1;
        std::size_t closest = 0;
        double closest_squared = std::numeric_limits<double>::max();
        for (std::size_t q = 0; q < 8; ++q) {
            const auto suffix = "_q" + std::to_string(q);
            const Point point = {snapshot.element("reference_x" + suffix).at(element),
                snapshot.element("reference_y" + suffix).at(element),
                snapshot.element("reference_z" + suffix).at(element)};
            const double distance_squared = std::pow(point.x - reference.position.x, 2) +
                                            std::pow(point.y - reference.position.y, 2) +
                                            std::pow(point.z - reference.position.z, 2);
            if (distance_squared < closest_squared) {
                closest = q;
                closest_squared = distance_squared;
            }
        }
        maximum_coordinate_error = std::max(maximum_coordinate_error, std::sqrt(closest_squared));
        if (!seen_points.emplace(reference.stage, element, closest).second)
            throw std::invalid_argument("B5.9 integration mapping is not unique");
        const auto suffix = "_q" + std::to_string(closest);
        const double material_temperature = snapshot.element("material_temperature" + suffix).at(element);
        std::array<double, 3> actual_heat_flux = {snapshot.element("heat_flux_x" + suffix).at(element),
            snapshot.element("heat_flux_y" + suffix).at(element), snapshot.element("heat_flux_z" + suffix).at(element)};
        const auto tensor = [&](const std::string& prefix) {
            return Tensor{snapshot.element(prefix + "xx" + suffix).at(element),
                snapshot.element(prefix + "yy" + suffix).at(element),
                snapshot.element(prefix + "zz" + suffix).at(element),
                snapshot.element(prefix + "xy" + suffix).at(element),
                snapshot.element(prefix + "yz" + suffix).at(element),
                snapshot.element(prefix + "xz" + suffix).at(element)};
        };
        const Tensor actual_stress = tensor("stress_");
        const Tensor actual_strain = tensor("infinitesimal_strain_");
        for (std::size_t component = 0; component < 3; ++component)
            integration_metrics[component].add(actual_heat_flux[component], reference.heat_flux[component]);
        const std::array<double, 6> actual_stress_components = {
            actual_stress.xx, actual_stress.yy, actual_stress.zz, actual_stress.xy, actual_stress.yz, actual_stress.xz};
        const std::array<double, 6> reference_stress_components = {reference.stress.xx, reference.stress.yy,
            reference.stress.zz, reference.stress.xy, reference.stress.yz, reference.stress.xz};
        const std::array<double, 6> actual_strain_components = {
            actual_strain.xx, actual_strain.yy, actual_strain.zz, actual_strain.xy, actual_strain.yz, actual_strain.xz};
        const std::array<double, 6> reference_strain_components = {reference.strain.xx, reference.strain.yy,
            reference.strain.zz, reference.strain.xy, reference.strain.yz, reference.strain.xz};
        for (std::size_t component = 0; component < 6; ++component) {
            integration_metrics[3 + component].add(
                actual_stress_components[component], reference_stress_components[component]);
            integration_metrics[9 + component].add(
                actual_strain_components[component], reference_strain_components[component]);
        }
        integration_metrics[15].add(material_temperature, reference.temperature);
        integration_metrics[16].add(
            snapshot.element("integration_measure" + suffix).at(element), reference.integration_volume);
        ++compared_points;
    }
    passed = check(compared_points == 32 * element_count,
                 std::string("B5.9 ") + spec.name + " compares all eight points in every element and stage") &&
             passed;
    const std::array<std::string, 17> integration_names = {"heat_flux_x", "heat_flux_y", "heat_flux_z", "stress_xx",
        "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz",
        "strain_xy", "strain_yz", "strain_xz", "material_temperature", "integration_volume"};
    const std::array<double, 17> integration_zero_tolerances = {1.0e-8, 1.0e-8, 1.0e-8, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3,
        1.0e-3, 1.0e-3, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-14};
    for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
        print_metrics(std::string("b59_") + spec.name + '_' + integration_names[field], integration_metrics[field]);
        passed =
            check(metrics_pass(integration_metrics[field], 1.0e-3, integration_zero_tolerances[field]),
                std::string("B5.9 ") + spec.name + ' ' + integration_names[field] + " metrics are below 0.1 percent") &&
            passed;
    }
    passed = check(maximum_coordinate_error < 1.0e-12,
                 std::string("B5.9 ") + spec.name + " maps all Abaqus integration points uniquely") &&
             passed;
    const auto bending = [&](std::size_t stage) {
        const auto& values = history.at(stage).nodal("displacement_x");
        double top = 0, bottom = 0;
        const std::size_t ny = 2 * spec.subdivisions_per_layer;
        for (std::size_t z = 0; z <= spec.nz; ++z) {
            top += values.at(z * (ny + 1) * (spec.nx + 1) + ny * (spec.nx + 1) + spec.nx);
            bottom += values.at(z * (ny + 1) * (spec.nx + 1) + spec.nx);
        }
        return (top - bottom) / static_cast<double>(spec.nz + 1);
    };
    const BendingSummary result{bending(2), bending(4)};
    std::cout << "b59_" << spec.name << "_peak_axial_bending_indicator=" << result.peak << " b59_" << spec.name
              << "_cooled_axial_bending_indicator=" << result.cooled << '\n';
    passed =
        check(std::abs(result.peak) > 1.0e-4,
            std::string("B5.9 ") + spec.name + " develops resolved two-layer thermal bending") &&
        check(std::abs(result.cooled) > 1.0e-5 && std::abs(result.cooled) < 0.2 * std::abs(result.peak),
            std::string("B5.9 ") + spec.name + " retains the smaller deformation caused by the cooled 330 K state") &&
        passed;
    return result;
}
} // namespace

namespace fuelsim::test {
bool check_hex8_multimaterial(const std::string& coarse, const std::string& refined, const std::string& distorted,
    const std::string& nodal_path, const std::string& integration_path) {
    const auto nodes = read_nodes(nodal_path);
    const auto points = read_integration(integration_path);
    const std::array<std::string, 3> paths = {coarse, refined, distorted};
    std::array<BendingSummary, 3> bending{};
    bool passed = true;
    for (std::size_t i = 0; i < 3; ++i) bending[i] = run_case(cases[i], paths[i], nodes, points, passed);
    const double sensitivity = std::abs(bending[2].peak - bending[1].peak) / std::abs(bending[1].peak);
    std::cout << "b59_refined_distorted_bending_relative_difference=" << sensitivity << '\n';
    return check(sensitivity < 0.1, "B5.9 distorted/refined bending difference stays below 10 percent") && passed;
}
} // namespace fuelsim::test
