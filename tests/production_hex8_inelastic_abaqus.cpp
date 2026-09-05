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
#include <utility>
#include <vector>

namespace {
struct Point {
    double x, y, z;
};

struct Tensor {
    double xx, yy, zz, xy, yz, xz;
};

constexpr double standard_time_step = 0.1;

enum class Branch {
    plastic,
    coupled,
    noncoaxial,
    finite_plastic,
    finite_coupled,
    finite_noncoaxial,
    reduced_plastic,
    reduced_coupled,
    finite_reduced_plastic,
    finite_reduced_coupled,
    finite_reduced_noncoaxial
};

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
    Tensor stress, strain, elastic_strain, plastic_strain, creep_strain;
    double equivalent_plastic_strain, equivalent_creep_strain, integration_volume;
};

struct EnergyReference final {
    std::size_t stage;
    double time, internal_energy, plastic_dissipation, creep_dissipation, external_work;
};

const char* case_id(Branch branch) {
    switch (branch) {
    case Branch::plastic: return "B5.10";
    case Branch::coupled: return "B5.12";
    case Branch::noncoaxial: return "B5.13";
    case Branch::finite_plastic: return "B5.15";
    case Branch::finite_coupled: return "B5.17";
    case Branch::finite_noncoaxial: return "B5.18";
    case Branch::reduced_plastic: return "B5.41";
    case Branch::reduced_coupled: return "B5.43";
    case Branch::finite_reduced_plastic: return "B5.35";
    case Branch::finite_reduced_coupled: return "B5.37";
    case Branch::finite_reduced_noncoaxial: return "B5.45";
    }
    throw std::logic_error("Unknown Abaqus inelastic branch");
}

const char* metric_prefix(Branch branch) {
    switch (branch) {
    case Branch::plastic: return "b510_";
    case Branch::coupled: return "b512_";
    case Branch::noncoaxial: return "b513_";
    case Branch::finite_plastic: return "b515_";
    case Branch::finite_coupled: return "b517_";
    case Branch::finite_noncoaxial: return "b518_";
    case Branch::reduced_plastic: return "b541_";
    case Branch::reduced_coupled: return "b543_";
    case Branch::finite_reduced_plastic: return "b535_";
    case Branch::finite_reduced_coupled: return "b537_";
    case Branch::finite_reduced_noncoaxial: return "b545_";
    }
    throw std::logic_error("Unknown Abaqus inelastic branch");
}

bool finite_strain(Branch branch) {
    return branch == Branch::finite_plastic || branch == Branch::finite_coupled ||
           branch == Branch::finite_noncoaxial || branch == Branch::finite_reduced_plastic ||
           branch == Branch::finite_reduced_coupled || branch == Branch::finite_reduced_noncoaxial;
}

bool reduced_integration(Branch branch) {
    return branch == Branch::reduced_plastic || branch == Branch::reduced_coupled ||
           branch == Branch::finite_reduced_plastic || branch == Branch::finite_reduced_coupled ||
           branch == Branch::finite_reduced_noncoaxial;
}

bool plastic_only(Branch branch) {
    return branch == Branch::plastic || branch == Branch::reduced_plastic || branch == Branch::finite_plastic ||
           branch == Branch::finite_reduced_plastic;
}

bool noncoaxial(Branch branch) {
    return branch == Branch::noncoaxial || branch == Branch::finite_noncoaxial ||
           branch == Branch::finite_reduced_noncoaxial;
}

bool long_noncoaxial(Branch branch) { return branch == Branch::finite_reduced_noncoaxial; }

bool monotonic_coupled(Branch branch) {
    return branch == Branch::coupled || branch == Branch::reduced_coupled || branch == Branch::finite_coupled ||
           branch == Branch::finite_reduced_coupled;
}

std::size_t stage_count(Branch branch) { return long_noncoaxial(branch) ? 100 : noncoaxial(branch) ? 20 : 10; }

double time_step(Branch branch) { return noncoaxial(branch) ? 0.001 : standard_time_step; }

Branch parse_branch(const std::string& value) {
    if (value == "plastic") return Branch::plastic;
    if (value == "coupled") return Branch::coupled;
    if (value == "noncoaxial") return Branch::noncoaxial;
    if (value == "finite_plastic") return Branch::finite_plastic;
    if (value == "finite_coupled") return Branch::finite_coupled;
    if (value == "finite_noncoaxial") return Branch::finite_noncoaxial;
    if (value == "reduced_plastic") return Branch::reduced_plastic;
    if (value == "reduced_coupled") return Branch::reduced_coupled;
    if (value == "finite_reduced_plastic") return Branch::finite_reduced_plastic;
    if (value == "finite_reduced_coupled") return Branch::finite_reduced_coupled;
    if (value == "finite_reduced_noncoaxial") return Branch::finite_reduced_noncoaxial;
    throw std::invalid_argument("Abaqus inelastic branch is not recognized");
}

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
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.10 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.10 index is not a positive integer in " + path);
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

std::vector<NodeReference> read_nodes(const std::string& path, Branch branch) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.10 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.10 nodal header in " + path);
    std::vector<NodeReference> result;
    double maximum_active_reaction = 0.0;
    double maximum_raw_theoretical_zero_reaction = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11) throw std::invalid_argument("Unexpected Abaqus B5.10 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 3, path);
            if (field == 5 || field == 6)
                maximum_active_reaction = std::max(maximum_active_reaction, std::abs(fields[field]));
            const bool theoretical_zero_reaction =
                (noncoaxial(branch) && field == 7) || (!noncoaxial(branch) && (field == 6 || field == 7));
            if (theoretical_zero_reaction) {
                maximum_raw_theoretical_zero_reaction =
                    std::max(maximum_raw_theoretical_zero_reaction, std::abs(fields[field]));
                fields[field] = 0.0;
            } else if (field != 0)
                fields[field] = zero_noise(fields[field], field >= 5 ? 1.0e-3 : (field == 4 ? 1.0e-10 : 1.0e-16));
        }
        result.push_back({positive_integer(number(values, 0, path), path),
            positive_integer(number(values, 2, path), path), number(values, 1, path), fields});
    }
    if (result.size() != 8 * stage_count(branch))
        throw std::invalid_argument("Abaqus small-strain nodal reference has an unexpected row count");
    {
        const double relative_residual = maximum_raw_theoretical_zero_reaction / maximum_active_reaction;
        std::cout << metric_prefix(branch)
                  << "abaqus_raw_theoretical_zero_reaction_maximum_absolute=" << maximum_raw_theoretical_zero_reaction
                  << '\n';
        std::cout << metric_prefix(branch)
                  << "abaqus_raw_theoretical_zero_reaction_relative_active_scale=" << relative_residual << '\n';
        if (!(maximum_active_reaction > 0.0 && relative_residual < 1.0e-7))
            throw std::invalid_argument(
                "Abaqus theoretical-zero reaction residual is too large to canonicalize as zero");
    }
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path, Branch branch) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.10 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string strain_header = finite_strain(branch)
                                          ? "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,"
                                          : "e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,";
    const std::string common = "stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
                               "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa," +
                               strain_header +
                               "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
                               "pe11,pe22,pe33,pe12_engineering,pe13_engineering,pe23_engineering,peeq,";
    const std::string expected_header =
        common + (plastic_only(branch)
                         ? "ivol_m3"
                         : "ce11,ce22,ce33,ce12_engineering,ce13_engineering,ce23_engineering,ceeq,ivol_m3");
    if (line != expected_header)
        throw std::invalid_argument("Unexpected Abaqus B5.10 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    double maximum_active_stress = 0.0;
    double maximum_raw_free_z_stress = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != (plastic_only(branch) ? 34 : 41))
            throw std::invalid_argument("Unexpected Abaqus B5.10 integration-point column count in " + path);
        if (positive_integer(number(values, 2, path), path) != 1)
            throw std::invalid_argument("Abaqus B5.10 reference must contain one element");
        Tensor stress = {number(values, 8, path), number(values, 9, path), number(values, 10, path),
            number(values, 11, path), number(values, 13, path), number(values, 12, path)};
        maximum_active_stress =
            std::max({maximum_active_stress, std::abs(stress.xx), std::abs(stress.yy), std::abs(stress.xy)});
        if (noncoaxial(branch)) maximum_raw_free_z_stress = std::max(maximum_raw_free_z_stress, std::abs(stress.zz));
        for (double* value : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *value = zero_noise(*value, 1.0);
        if (noncoaxial(branch)) stress.zz = 0.0;
        const Tensor creep = plastic_only(branch) ? Tensor{} : tensor(values, 33, path, 1.0e-15);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}, number(values, 7, path),
            stress, tensor(values, 14, path, 1.0e-15), tensor(values, 20, path, 1.0e-15),
            tensor(values, 26, path, 1.0e-15), creep, number(values, 32, path),
            plastic_only(branch) ? 0.0 : number(values, 39, path),
            number(values, plastic_only(branch) ? 33 : 40, path)});
    }
    if (result.size() != (reduced_integration(branch) ? 1 : 8) * stage_count(branch))
        throw std::invalid_argument("Abaqus small-strain integration-point reference has an unexpected row count");
    if (noncoaxial(branch)) {
        const double relative_residual = maximum_raw_free_z_stress / maximum_active_stress;
        std::cout << metric_prefix(branch) << "abaqus_raw_free_z_stress_maximum_absolute=" << maximum_raw_free_z_stress
                  << '\n';
        std::cout << metric_prefix(branch) << "abaqus_raw_free_z_stress_relative_active_scale=" << relative_residual
                  << '\n';
        if (!(maximum_active_stress > 0.0 && relative_residual < 1.0e-7))
            throw std::invalid_argument("Abaqus B5.13 free-z stress residual is too large to canonicalize as zero");
    }
    return result;
}

std::vector<EnergyReference> read_energy(const std::string& path, Branch branch) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.10 energy: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected_header =
        plastic_only(branch)
            ? "stage,time_s,internal_energy_j,plastic_dissipation_j,external_work_j"
            : "stage,time_s,internal_energy_j,plastic_dissipation_j,creep_dissipation_j,external_work_j";
    if (line != expected_header) throw std::invalid_argument("Unexpected Abaqus B5.10 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != (plastic_only(branch) ? 5 : 6))
            throw std::invalid_argument("Unexpected Abaqus B5.10 energy column count in " + path);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            number(values, 2, path), number(values, 3, path), plastic_only(branch) ? 0.0 : number(values, 4, path),
            number(values, plastic_only(branch) ? 4 : 5, path)});
    }
    if (result.size() != stage_count(branch))
        throw std::invalid_argument("Abaqus small-strain energy reference has an unexpected row count");
    return result;
}

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance,
    double qualified_pointwise_absolute_tolerance = 0.0) {
    if (metrics.has_relative_norm()) {
        const bool aggregate_passed =
            metrics.relative_l2() < relative_tolerance && metrics.relative_absolute_peak() < relative_tolerance;
        const double pointwise_absolute_difference =
            std::abs(metrics.maximum_pointwise_relative_actual - metrics.maximum_pointwise_relative_reference);
        const bool pointwise_passed = metrics.maximum_pointwise_relative_error() < relative_tolerance ||
                                      (qualified_pointwise_absolute_tolerance > 0.0 &&
                                          pointwise_absolute_difference < qualified_pointwise_absolute_tolerance);
        if (!aggregate_passed || !pointwise_passed) return false;
    }
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
bool check_hex8_inelastic_abaqus(const std::string& output_path, const std::string& branch_name,
    const std::string& nodal_path, const std::string& integration_path, const std::string& energy_path) {
    const Branch branch = parse_branch(branch_name);
    const auto node_reference = read_nodes(nodal_path, branch);
    const auto integration_reference = read_integration(integration_path, branch);
    const auto energy_reference = read_energy(energy_path, branch);
    const std::size_t stages = stage_count(branch);
    auto history = read_exodus_history(output_path);
    std::map<std::size_t, ExodusResults> storage;
    std::map<std::size_t, const ExodusResults*> snapshots;
    bool passed = true;
    for (std::size_t stage = 1; stage <= stages; ++stage) {
        auto frame = std::move(history.at(stage));
        if (frame.step_count != stages + 1 || frame.nodes.size() != 8 ||
            std::abs(frame.time - time_step(branch) * static_cast<double>(stage)) > 1e-12)
            throw std::invalid_argument("Inelastic output does not contain the prescribed full path");
        storage.emplace(stage, std::move(frame));
        snapshots.emplace(stage, &storage.at(stage));
        passed = check(storage.at(stage).global("conservation_plastic_dissipation_increment") >= -1e-8 &&
                           storage.at(stage).global("conservation_creep_dissipation_increment") >= -1e-8,
                     "Inelastic dissipation is nonnegative at every stage") &&
                 passed;
    }
    std::array<FieldErrorMetrics, 8> nodal_metrics;
    const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y", "displacement_z",
        "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
    std::set<std::pair<std::size_t, std::size_t>> nodes_seen;
    for (const auto& reference : node_reference) {
        const auto& frame = storage.at(reference.stage);
        if (reference.node < 1 || reference.node > 8 || std::abs(reference.time - frame.time) > 1e-12 ||
            !nodes_seen.emplace(reference.stage, reference.node).second)
            throw std::invalid_argument("Inelastic reference node or time mapping is not unique");
        for (std::size_t field = 0; field < 8; ++field)
            nodal_metrics[field].add(frame.nodal(nodal_names[field]).at(reference.node - 1), reference.fields[field]);
    }
    const std::array<double, 8> nodal_zero_tolerances = {
        1.0e-12, 1.0e-15, 1.0e-15, 1.0e-15, 1.0e-10, 1.0e-2, 1.0e-2, 1.0e-2};
    for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
        print_metrics(metric_prefix(branch) + nodal_names[field], nodal_metrics[field]);
        passed =
            check(metrics_pass(nodal_metrics[field], 1.0e-3, nodal_zero_tolerances[field]),
                std::string(case_id(branch)) + " " + nodal_names[field] + " metrics are below the acceptance limits") &&
            passed;
    }

    std::array<FieldErrorMetrics, 34> integration_metrics;
    std::set<std::pair<std::size_t, std::size_t>> points_seen;
    double maximum_coordinate_error = 0;
    for (const auto& reference : integration_reference) {
        const auto& frame = storage.at(reference.stage);
        if (std::abs(frame.time - reference.time) > 1e-12)
            throw std::invalid_argument("Inelastic integration reference time differs");
        std::size_t closest = 0;
        double distance = std::numeric_limits<double>::infinity();
        for (std::size_t q = 0; q < (reduced_integration(branch) ? 1u : 8u); ++q) {
            const std::string prefix = finite_strain(branch) ? "current_" : "reference_";
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
            output_tensor(frame, finite_strain(branch) ? "logarithmic_strain_" : "infinitesimal_strain_", closest);
        const auto actual_elastic = output_tensor(frame, "elastic_", closest);
        const auto actual_plastic = output_tensor(frame, "plastic_", closest);
        const auto actual_creep = output_tensor(frame, "creep_", closest);
        const auto expected_stress = components(reference.stress);
        const auto expected_total = components(reference.strain);
        const auto expected_elastic = components(reference.elastic_strain);
        const auto expected_plastic = components(reference.plastic_strain);
        const auto expected_creep = components(reference.creep_strain);
        for (std::size_t component = 0; component < 6; ++component) {
            integration_metrics[component].add(actual_stress[component], expected_stress[component]);
            integration_metrics[6 + component].add(actual_total[component], expected_total[component]);
            integration_metrics[12 + component].add(actual_elastic[component], expected_elastic[component]);
            integration_metrics[18 + component].add(actual_plastic[component], expected_plastic[component]);
            integration_metrics[24 + component].add(actual_creep[component], expected_creep[component]);
        }
        const auto suffix = "_q" + std::to_string(closest);
        integration_metrics[30].add(frame.element("equiv_plastic" + suffix).at(0), reference.equivalent_plastic_strain);
        integration_metrics[31].add(frame.element("equiv_creep" + suffix).at(0), reference.equivalent_creep_strain);
        integration_metrics[32].add(frame.element("material_temperature" + suffix).at(0), reference.temperature);
        integration_metrics[33].add(frame.element("current_measure" + suffix).at(0), reference.integration_volume);
    }
    const std::array<std::string, 34> integration_names = {"stress_xx", "stress_yy", "stress_zz", "stress_xy",
        "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz",
        "elastic_strain_xx", "elastic_strain_yy", "elastic_strain_zz", "elastic_strain_xy", "elastic_strain_yz",
        "elastic_strain_xz", "plastic_strain_xx", "plastic_strain_yy", "plastic_strain_zz", "plastic_strain_xy",
        "plastic_strain_yz", "plastic_strain_xz", "creep_strain_xx", "creep_strain_yy", "creep_strain_zz",
        "creep_strain_xy", "creep_strain_yz", "creep_strain_xz", "equivalent_plastic_strain", "equivalent_creep_strain",
        "material_temperature", "integration_volume"};
    for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
        print_metrics(metric_prefix(branch) + integration_names[field], integration_metrics[field]);
        const double zero_tolerance = long_noncoaxial(branch) && field == 2 ? 1.0e-1 : field < 6 ? 1.0e-2 : 1.0e-14;
        const double qualified_pointwise_absolute_tolerance = long_noncoaxial(branch) && field == 14 ? 1.0e-10 : 0.0;
        passed =
            check(metrics_pass(
                      integration_metrics[field], 1.0e-3, zero_tolerance, qualified_pointwise_absolute_tolerance),
                std::string(case_id(branch)) + " " + integration_names[field] + " metrics are below 0.1 percent") &&
            passed;
    }
    passed = check(maximum_coordinate_error < 1.0e-11,
                 std::string(case_id(branch)) + " maps every integration point at every accepted time") &&
             passed;
    std::cout << metric_prefix(branch)
              << "maximum_integration_coordinate_absolute_difference=" << maximum_coordinate_error << '\n';

    std::array<fuelsim::test::FieldErrorMetrics, 4> energy_metrics;
    double cumulative_elastic = 0.0, cumulative_plastic = 0.0, cumulative_creep = 0.0;
    std::size_t energy_stage = 0;
    for (const EnergyReference& reference : energy_reference) {
        const ExodusResults& snapshot = *snapshots.at(reference.stage);
        if (reference.stage != ++energy_stage || std::abs(reference.time - snapshot.time) > 1e-12)
            throw std::invalid_argument("Inelastic energy reference must contain every stage in time order");
        cumulative_elastic += snapshot.global("conservation_elastic_energy_change");
        cumulative_plastic += snapshot.global("conservation_plastic_dissipation_increment");
        cumulative_creep += snapshot.global("conservation_creep_dissipation_increment");
        energy_metrics[0].add(cumulative_elastic + cumulative_plastic + cumulative_creep, reference.internal_energy);
        energy_metrics[1].add(cumulative_plastic, reference.plastic_dissipation);
        energy_metrics[2].add(cumulative_creep, reference.creep_dissipation);
        energy_metrics[3].add(cumulative_elastic + cumulative_plastic + cumulative_creep, reference.external_work);
    }
    const std::array<std::string, 4> energy_names = {
        "internal_energy", "plastic_dissipation", "creep_dissipation", "external_work"};
    for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
        print_metrics(metric_prefix(branch) + energy_names[field], energy_metrics[field]);
        passed = check(metrics_pass(energy_metrics[field], 1.0e-3, 1.0e-8),
                     std::string(case_id(branch)) + " " + energy_names[field] + " metrics are below 0.1 percent") &&
                 passed;
    }
    if (plastic_only(branch))
        passed = check(snapshots.at(3)->element("equiv_plastic_q0").at(0) > 0.0 &&
                           snapshots.at(4)->element("equiv_plastic_q0").at(0) ==
                               snapshots.at(3)->element("equiv_plastic_q0").at(0) &&
                           snapshots.at(7)->element("equiv_plastic_q0").at(0) >
                               snapshots.at(4)->element("equiv_plastic_q0").at(0) &&
                           snapshots.at(9)->element("equiv_plastic_q0").at(0) >
                               snapshots.at(7)->element("equiv_plastic_q0").at(0),
                     std::string(case_id(branch)) +
                         " activates plasticity, preserves it during elastic unload, and accumulates it on reversal") &&
                 passed;
    else if (monotonic_coupled(branch))
        passed = check(snapshots.at(1)->element("equiv_plastic_q0").at(0) > 0.0 &&
                           snapshots.at(1)->element("equiv_creep_q0").at(0) > 0.0 &&
                           snapshots.at(stages - 1)->element("equiv_plastic_q0").at(0) >
                               snapshots.at(1)->element("equiv_plastic_q0").at(0) &&
                           snapshots.at(stages - 1)->element("equiv_creep_q0").at(0) >
                               snapshots.at(1)->element("equiv_creep_q0").at(0),
                     std::string(case_id(branch)) +
                         " keeps plasticity and creep active and accumulating throughout the path") &&
                 passed;
    else if (long_noncoaxial(branch))
        passed =
            check(snapshots.at(1)->element("equiv_plastic_q0").at(0) > 0.0 &&
                      snapshots.at(1)->element("equiv_creep_q0").at(0) > 0.0 &&
                      snapshots.at(stages - 1)->element("equiv_plastic_q0").at(0) >
                          snapshots.at(1)->element("equiv_plastic_q0").at(0) &&
                      snapshots.at(stages - 1)->element("equiv_creep_q0").at(0) >
                          snapshots.at(1)->element("equiv_creep_q0").at(0) &&
                      snapshots.at(stages - 1)->element("plastic_xy_q0").at(0) > 0.0,
                std::string(case_id(branch)) + " keeps both mechanisms active while the loading direction changes") &&
            passed;
    else
        passed = check(snapshots.at(1)->element("equiv_plastic_q0").at(0) > 0.0 &&
                           snapshots.at(1)->element("equiv_creep_q0").at(0) > 0.0 &&
                           snapshots.at(stages - 1)->element("equiv_plastic_q0").at(0) >
                               snapshots.at(1)->element("equiv_plastic_q0").at(0) &&
                           snapshots.at(stages - 1)->element("equiv_creep_q0").at(0) >
                               snapshots.at(1)->element("equiv_creep_q0").at(0) &&
                           snapshots.at(9)->element("plastic_xy_q0").at(0) *
                                   snapshots.at(stages - 1)->element("plastic_xy_q0").at(0) <
                               0.0,
                     std::string(case_id(branch)) +
                         " accumulates both mechanisms and reverses the plastic shear direction") &&
                 passed;
    if (long_noncoaxial(branch)) {
        double maximum_rotation_degrees = 0.0;
        for (std::size_t stage = 0; stage < stages; ++stage) {
            const auto& frame = storage.at(stage + 1);
            const double axial = frame.nodal("displacement_x").at(1) - frame.nodal("displacement_x").at(0);
            const double shear = frame.nodal("displacement_y").at(1) - frame.nodal("displacement_y").at(0);
            maximum_rotation_degrees =
                std::max(maximum_rotation_degrees, std::abs(std::atan2(shear, 2.0 + axial)) * 180.0 / std::acos(-1.0));
        }
        std::cout << metric_prefix(branch) << "maximum_polar_rotation_degrees=" << maximum_rotation_degrees << '\n';
        passed = check(maximum_rotation_degrees > 25.0,
                     std::string(case_id(branch)) + " prescribed path exceeds 25 degrees of polar rotation") &&
                 passed;
    }
    return passed;
}
} // namespace fuelsim::test
