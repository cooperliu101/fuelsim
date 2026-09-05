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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Point {
    double x, y, z;
};

struct Tensor {
    double xx, yy, zz, xy, yz, xz;
};

constexpr double preload_time = 1.0e-9;
constexpr double time_step = 0.1;

const char* case_id(bool finite_strain, bool reduced_integration) {
    return reduced_integration ? (finite_strain ? "B5.36" : "B5.42") : (finite_strain ? "B5.16" : "B5.11");
}

const char* metric_prefix(bool finite_strain, bool reduced_integration) {
    return reduced_integration ? (finite_strain ? "b536_" : "b542_") : (finite_strain ? "b516_" : "b511_");
}

struct NodeReference final {
    std::size_t stage, node;
    double time;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t stage;
    double time;
    Point position;
    double temperature;
    Tensor stress, strain, elastic_strain, creep_strain;
    double equivalent_creep_strain, integration_volume;
};

struct EnergyReference final {
    std::size_t stage;
    double time, internal_energy, creep_dissipation, external_work;
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
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.11 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.11 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

double zero_noise(double value, double threshold) { return std::abs(value) < threshold ? 0.0 : value; }

Tensor tensor(
    const std::vector<std::string>& values, std::size_t first, const std::string& path, double zero_threshold) {
    Tensor result = {number(values, first, path), number(values, first + 1, path), number(values, first + 2, path),
        0.5 * number(values, first + 3, path), 0.5 * number(values, first + 5, path),
        0.5 * number(values, first + 4, path)};
    for (double* value : {&result.xx, &result.yy, &result.zz, &result.xy, &result.yz, &result.xz})
        *value = zero_noise(*value, zero_threshold);
    return result;
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.11 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.11 nodal header in " + path);
    std::vector<NodeReference> result;
    double maximum_active_reaction = 0.0, maximum_raw_transverse_reaction = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11) throw std::invalid_argument("Unexpected Abaqus B5.11 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 3, path);
            if (field == 5) maximum_active_reaction = std::max(maximum_active_reaction, std::abs(fields[field]));
            if (field == 6 || field == 7) {
                maximum_raw_transverse_reaction = std::max(maximum_raw_transverse_reaction, std::abs(fields[field]));
                fields[field] = 0.0;
            } else if (field != 0)
                fields[field] = zero_noise(fields[field], field >= 5 ? 1.0e-4 : (field == 4 ? 1.0e-10 : 1.0e-16));
        }
        result.push_back({positive_integer(number(values, 0, path), path),
            positive_integer(number(values, 2, path), path), number(values, 1, path), fields});
    }
    if (result.size() != 80) throw std::invalid_argument("Abaqus B5.11 nodal reference must contain 80 rows");
    const double relative_residual = maximum_raw_transverse_reaction / maximum_active_reaction;
    std::cout << "abaqus_raw_transverse_reaction_maximum_absolute=" << maximum_raw_transverse_reaction << '\n'
              << "abaqus_raw_transverse_reaction_relative_active_scale=" << relative_residual << '\n';
    if (!(maximum_active_reaction > 0.0 && relative_residual < 1.0e-7))
        throw std::invalid_argument("Abaqus transverse reaction residual is too large to canonicalize as zero");
    return result;
}

std::vector<IntegrationReference> read_integration(
    const std::string& path, bool finite_strain, bool reduced_integration) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.11 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string strain_header = finite_strain
                                          ? "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,"
                                          : "e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,";
    if (line != "stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
                "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa," +
                    strain_header +
                    "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
                    "ce11,ce22,ce33,ce12_engineering,ce13_engineering,ce23_engineering,ceeq,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.11 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    double maximum_active_stress = 0.0, maximum_raw_transverse_stress = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 34)
            throw std::invalid_argument("Unexpected Abaqus B5.11 integration-point column count in " + path);
        if (positive_integer(number(values, 2, path), path) != 1)
            throw std::invalid_argument("Abaqus B5.11 reference must contain one element");
        Tensor stress = {number(values, 8, path), number(values, 9, path), number(values, 10, path),
            number(values, 11, path), number(values, 13, path), number(values, 12, path)};
        maximum_active_stress = std::max(maximum_active_stress, std::abs(stress.xx));
        maximum_raw_transverse_stress =
            std::max({maximum_raw_transverse_stress, std::abs(stress.yy), std::abs(stress.zz)});
        for (double* value : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *value = zero_noise(*value, 1.0e-3);
        stress.yy = 0.0;
        stress.zz = 0.0;
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}, number(values, 7, path),
            stress, tensor(values, 14, path, 1.0e-15), tensor(values, 20, path, 1.0e-15),
            tensor(values, 26, path, 1.0e-15), number(values, 32, path), number(values, 33, path)});
    }
    if (result.size() != (reduced_integration ? 10 : 80))
        throw std::invalid_argument("Abaqus Norton integration-point reference has an unexpected row count");
    const double relative_residual = maximum_raw_transverse_stress / maximum_active_stress;
    std::cout << "abaqus_raw_transverse_stress_maximum_absolute=" << maximum_raw_transverse_stress << '\n'
              << "abaqus_raw_transverse_stress_relative_active_scale=" << relative_residual << '\n';
    if (!(maximum_active_stress > 0.0 && relative_residual < 1.0e-7))
        throw std::invalid_argument("Abaqus transverse stress residual is too large to canonicalize as zero");
    return result;
}

std::vector<EnergyReference> read_energy(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.11 energy: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,internal_energy_j,creep_dissipation_j,external_work_j")
        throw std::invalid_argument("Unexpected Abaqus B5.11 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5) throw std::invalid_argument("Unexpected Abaqus B5.11 energy column count in " + path);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            number(values, 2, path), number(values, 3, path), number(values, 4, path)});
    }
    if (result.size() != 10) throw std::invalid_argument("Abaqus B5.11 energy reference must contain ten rows");
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

std::array<double, 6> components(const Tensor& value) {
    return {value.xx, value.yy, value.zz, value.xy, value.yz, value.xz};
}

std::array<double, 6> output_tensor(
    const fuelsim::test::ExodusResults& output, const std::string& prefix, std::size_t q) {
    std::array<double, 6> values{};
    const std::array<std::string, 6> names = {"xx", "yy", "zz", "xy", "yz", "xz"};
    for (std::size_t component = 0; component < 6; ++component)
        values[component] = output.element(prefix + names[component] + "_q" + std::to_string(q)).at(0);
    return values;
}
} // namespace

namespace fuelsim::test {
bool check_hex8_norton_abaqus(const std::string& output_path, const std::string& branch, const std::string& nodal_path,
    const std::string& integration_path, const std::string& energy_path) {
    if (branch != "small" && branch != "finite" && branch != "reduced" && branch != "finite_reduced")
        throw std::invalid_argument("Unknown Norton reference branch");
    const bool finite_strain = branch == "finite" || branch == "finite_reduced";
    const bool reduced_integration = branch == "reduced" || branch == "finite_reduced";
    const auto node_reference = read_nodes(nodal_path);
    const auto integration_reference = read_integration(integration_path, finite_strain, reduced_integration);
    const auto energy_reference = read_energy(energy_path);
    auto history = read_exodus_history(output_path);
    if (history.size() != 12 || history.at(1).time != preload_time)
        throw std::invalid_argument("Norton output requires one preload and ten hold steps");
    std::map<std::size_t, ExodusResults> storage;
    std::map<std::size_t, const ExodusResults*> snapshots;
    bool passed = true;
    for (std::size_t stage = 1; stage <= 10; ++stage) {
        auto frame = std::move(history.at(stage + 1));
        if (frame.nodes.size() != 8 ||
            std::abs(frame.time - preload_time - time_step * static_cast<double>(stage)) > 1e-12)
            throw std::invalid_argument("Norton output does not contain every prescribed hold");
        storage.emplace(stage, std::move(frame));
        snapshots.emplace(stage, &storage.at(stage));
        passed = check(storage.at(stage).global("conservation_creep_dissipation_increment") >= -1e-8,
                     "Norton creep dissipation is nonnegative at every hold") &&
                 passed;
    }
    std::array<FieldErrorMetrics, 8> nodal_metrics;
    const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y", "displacement_z",
        "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
    std::set<std::pair<std::size_t, std::size_t>> nodes_seen;
    for (const auto& reference : node_reference) {
        const auto& frame = storage.at(reference.stage);
        if (reference.node < 1 || reference.node > 8 ||
            std::abs(reference.time - (frame.time - preload_time)) > 1e-12 ||
            !nodes_seen.emplace(reference.stage, reference.node).second)
            throw std::invalid_argument("Inelastic reference node or time mapping is not unique");
        for (std::size_t field = 0; field < 8; ++field)
            nodal_metrics[field].add(frame.nodal(nodal_names[field]).at(reference.node - 1), reference.fields[field]);
    }
    const std::array<double, 8> nodal_zero_tolerances = {
        1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-10, 1.0e-2, 1.0e-2, 1.0e-2};
    for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
        print_metrics(metric_prefix(finite_strain, reduced_integration) + nodal_names[field], nodal_metrics[field]);
        passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, nodal_zero_tolerances[field]),
                     std::string(case_id(finite_strain, reduced_integration)) + " " + nodal_names[field] +
                         " metrics are below the acceptance limits") &&
                 passed;
    }

    std::array<FieldErrorMetrics, 27> integration_metrics;
    std::set<std::pair<std::size_t, std::size_t>> points_seen;
    double maximum_coordinate_error = 0;
    for (const auto& reference : integration_reference) {
        const auto& frame = storage.at(reference.stage);
        if (std::abs(frame.time - preload_time - reference.time) > 1e-12)
            throw std::invalid_argument("Inelastic integration reference time differs");
        std::size_t closest = 0;
        double distance = std::numeric_limits<double>::infinity();
        for (std::size_t q = 0; q < (reduced_integration ? 1u : 8u); ++q) {
            const std::string prefix = finite_strain ? "current_" : "reference_";
            const auto suffix = "_q" + std::to_string(q);
            const double candidate = std::hypot(frame.element(prefix + "x" + suffix).at(0) - reference.position.x,
                frame.element(prefix + "y" + suffix).at(0) - reference.position.y,
                frame.element(prefix + "z" + suffix).at(0) - reference.position.z);
            if (candidate < distance) {
                closest = q;
                distance = candidate;
            }
        }
        if (!std::isfinite(distance) || !points_seen.emplace(reference.stage, closest).second)
            throw std::invalid_argument("Inelastic integration point association is not unique");
        maximum_coordinate_error = std::max(maximum_coordinate_error, distance);
        const auto actual_stress = output_tensor(frame, "stress_", closest);
        const auto actual_total =
            output_tensor(frame, finite_strain ? "logarithmic_strain_" : "infinitesimal_strain_", closest);
        const auto actual_elastic = output_tensor(frame, "elastic_", closest);
        const auto actual_creep = output_tensor(frame, "creep_", closest);
        const auto expected_stress = components(reference.stress);
        const auto expected_total = components(reference.strain);
        const auto expected_elastic = components(reference.elastic_strain);
        const auto expected_creep = components(reference.creep_strain);
        for (std::size_t component = 0; component < 6; ++component) {
            integration_metrics[component].add(actual_stress[component], expected_stress[component]);
            integration_metrics[6 + component].add(actual_total[component], expected_total[component]);
            integration_metrics[12 + component].add(actual_elastic[component], expected_elastic[component]);
            integration_metrics[18 + component].add(actual_creep[component], expected_creep[component]);
        }
        const auto suffix = "_q" + std::to_string(closest);
        integration_metrics[24].add(frame.element("equiv_creep" + suffix).at(0), reference.equivalent_creep_strain);
        integration_metrics[25].add(frame.element("material_temperature" + suffix).at(0), reference.temperature);
        integration_metrics[26].add(frame.element("current_measure" + suffix).at(0), reference.integration_volume);
    }
    const std::array<std::string, 27> integration_names = {"stress_xx", "stress_yy", "stress_zz", "stress_xy",
        "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz",
        "elastic_strain_xx", "elastic_strain_yy", "elastic_strain_zz", "elastic_strain_xy", "elastic_strain_yz",
        "elastic_strain_xz", "creep_strain_xx", "creep_strain_yy", "creep_strain_zz", "creep_strain_xy",
        "creep_strain_yz", "creep_strain_xz", "equivalent_creep_strain", "material_temperature", "integration_volume"};
    for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
        print_metrics(
            metric_prefix(finite_strain, reduced_integration) + integration_names[field], integration_metrics[field]);
        const double zero_tolerance = field < 6 ? 1.0e-2 : 1.0e-14;
        passed = check(metrics_pass(integration_metrics[field], 1.0e-3, zero_tolerance),
                     std::string(case_id(finite_strain, reduced_integration)) + " " + integration_names[field] +
                         " metrics are below 0.1 percent") &&
                 passed;
    }
    passed = check(maximum_coordinate_error < (finite_strain ? 1.0e-10 : 1.0e-12),
                 std::string(case_id(finite_strain, reduced_integration)) +
                     " maps every integration point at all ten accepted times") &&
             passed;
    std::cout << metric_prefix(finite_strain, reduced_integration)
              << "maximum_integration_coordinate_absolute_difference=" << maximum_coordinate_error << '\n';

    std::array<fuelsim::test::FieldErrorMetrics, 3> energy_metrics;
    double cumulative_elastic = history.at(1).global("conservation_elastic_energy_change");
    double cumulative_creep = history.at(1).global("conservation_creep_dissipation_increment");
    std::size_t expected_stage = 0;
    for (const EnergyReference& reference : energy_reference) {
        ++expected_stage;
        if (reference.stage != expected_stage ||
            std::abs(reference.time - time_step * static_cast<double>(expected_stage)) > 1e-12)
            throw std::invalid_argument("Norton energy history must cover the ten holds in time order");
        const ExodusResults& snapshot = *snapshots.at(reference.stage);
        cumulative_elastic += snapshot.global("conservation_elastic_energy_change");
        cumulative_creep += snapshot.global("conservation_creep_dissipation_increment");
        energy_metrics[0].add(cumulative_elastic + cumulative_creep, reference.internal_energy);
        energy_metrics[1].add(cumulative_creep, reference.creep_dissipation);
        energy_metrics[2].add(cumulative_elastic + cumulative_creep, reference.external_work);
    }
    const std::array<std::string, 3> energy_names = {"internal_energy", "creep_dissipation", "external_work"};
    for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
        print_metrics(metric_prefix(finite_strain, reduced_integration) + energy_names[field], energy_metrics[field]);
        passed = check(metrics_pass(energy_metrics[field], 1.0e-3, 1.0e-8),
                     std::string(case_id(finite_strain, reduced_integration)) + " " + energy_names[field] +
                         " metrics are below 0.1 percent") &&
                 passed;
    }
    passed = check(snapshots.at(1)->element("equiv_creep_q0").at(0) > 0.0 &&
                       snapshots.at(10)->element("equiv_creep_q0").at(0) >
                           9.9 * snapshots.at(1)->element("equiv_creep_q0").at(0),
                 std::string(case_id(finite_strain, reduced_integration)) +
                     " accumulates positive Norton creep throughout the ten constant-force holds") &&
             passed;
    return passed;
}
} // namespace fuelsim::test
