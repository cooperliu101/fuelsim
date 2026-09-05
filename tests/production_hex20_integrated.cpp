#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct Point final {
    double x, y, z;
};

struct Stress final {
    double xx, yy, zz, xy, yz, xz;
};

struct TemperatureReference final {
    std::size_t node;
    Point point;
    double value;
};

struct DisplacementReference final {
    std::size_t node;
    Point point;
    std::array<double, 3> value;
};

struct MaterialReference final {
    std::size_t element, point;
    Point position;
    double volume, equivalent_stress, equivalent_plastic_strain, equivalent_creep_strain;
};

struct ContactReference final {
    std::size_t node;
    Point point;
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

double distance(const Point& first, const Point& second) {
    return std::sqrt(
        std::pow(first.x - second.x, 2) + std::pow(first.y - second.y, 2) + std::pow(first.z - second.z, 2));
}

double equivalent_stress(const Stress& stress) {
    const double mean = (stress.xx + stress.yy + stress.zz) / 3.0;
    return std::sqrt(
        1.5 * (std::pow(stress.xx - mean, 2) + std::pow(stress.yy - mean, 2) + std::pow(stress.zz - mean, 2) +
                  2.0 * (stress.xy * stress.xy + stress.yz * stress.yz + stress.xz * stress.xz)));
}

bool run(const std::string& output_path, const std::string& temperature_path, const std::string& displacement_path,
    const std::string& material_path, const std::string& contact_path, bool frictional, std::size_t summary_active) {
    const auto temperature_reference = read_temperature(temperature_path);
    const auto displacement_reference = read_displacement(displacement_path);
    const auto material_reference = read_material(material_path);
    const std::vector<ContactReference> contact_reference =
        contact_path.empty() ? std::vector<ContactReference>{} : read_contact(contact_path);
    const auto output = fuelsim::test::read_final_exodus_results(output_path);
    const std::string prefix = frictional ? "b555_" : contact_reference.empty() ? "b549_" : "b550_";
    if (output.nodes.size() != 112 || output.block_element_counts.size() != 2 ||
        output.block_element_counts[0] + output.block_element_counts[1] != 8)
        throw std::invalid_argument("B5.49 output must preserve 112 source nodes and eight two-region elements");
    bool passed = true;
    const auto point_at = [&](std::size_t node) {
        const auto& p = output.nodes.at(node);
        return Point{p[0], p[1], p[2]};
    };
    std::set<std::size_t> temperature_nodes, displacement_nodes;
    for (const auto& reference : temperature_reference)
        if (!temperature_nodes.insert(reference.node).second)
            throw std::invalid_argument("Repeated temperature reference node");
    for (const auto& reference : displacement_reference)
        if (!displacement_nodes.insert(reference.node).second)
            throw std::invalid_argument("Repeated displacement reference node");
    fuelsim::test::FieldErrorMetrics temperature_metrics;
    fuelsim::test::GroupedFieldErrorMetrics displacement_metrics;
    double maximum_reference_coordinate_difference = 0.0;
    for (const auto& reference : temperature_reference) {
        maximum_reference_coordinate_difference =
            std::max(maximum_reference_coordinate_difference, distance(point_at(reference.node), reference.point));
        temperature_metrics.add(output.nodal("temperature").at(reference.node), reference.value);
    }
    for (const auto& reference : displacement_reference) {
        maximum_reference_coordinate_difference =
            std::max(maximum_reference_coordinate_difference, distance(point_at(reference.node), reference.point));
        std::array<double, 3> actual{};
        for (std::size_t component = 0; component < 3; ++component)
            actual[component] = output.nodal("displacement_" + std::string(1, "xyz"[component])).at(reference.node);
        displacement_metrics.add(actual.data(), reference.value.data(), actual.size());
    }

    std::map<std::size_t, std::vector<MaterialReference>> references_by_element;
    for (const auto& reference : material_reference) references_by_element[reference.element].push_back(reference);
    fuelsim::test::FieldErrorMetrics stress_metrics, plastic_metrics, creep_metrics;
    double maximum_material_coordinate_difference = 0.0,
           minimum_material_second_to_first_distance_ratio = std::numeric_limits<double>::infinity();
    for (std::size_t element = 0; element < 8; ++element) {
        const auto& references = references_by_element.at(element + 1);
        if (references.size() != 27)
            throw std::invalid_argument("Material reference requires 27 unique points per element");
        std::array<Point, 27> positions{};
        for (std::size_t q = 0; q < 27; ++q) {
            const std::string suffix = "_q" + std::to_string(q);
            positions[q] = {output.element("current_x" + suffix).at(element),
                output.element("current_y" + suffix).at(element), output.element("current_z" + suffix).at(element)};
            if (!std::isfinite(positions[q].x) || !std::isfinite(positions[q].y) || !std::isfinite(positions[q].z))
                throw std::invalid_argument("C3D20T material-point current coordinate is not finite");
        }
        std::vector<std::tuple<double, std::size_t, std::size_t>> candidates;
        for (std::size_t actual = 0; actual < 27; ++actual)
            for (std::size_t reference = 0; reference < 27; ++reference)
                candidates.emplace_back(distance(positions[actual], references[reference].position), actual, reference);
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
            const std::string suffix = "_q" + std::to_string(actual);
            const Stress stress = {output.element("stress_xx" + suffix).at(element),
                output.element("stress_yy" + suffix).at(element), output.element("stress_zz" + suffix).at(element),
                output.element("stress_xy" + suffix).at(element), output.element("stress_yz" + suffix).at(element),
                output.element("stress_xz" + suffix).at(element)};
            stress_metrics.add(equivalent_stress(stress), references[reference].equivalent_stress);
            plastic_metrics.add(
                output.element("equiv_plastic" + suffix).at(element), references[reference].equivalent_plastic_strain);
            creep_metrics.add(
                output.element("equiv_creep" + suffix).at(element), references[reference].equivalent_creep_strain);
        }
        if (pair_count != 27) throw std::logic_error("B5.49 material-point association is incomplete");
    }

    fuelsim::test::print_relative_metrics(prefix + "temperature", temperature_metrics);
    fuelsim::test::print_grouped_relative_metrics(prefix + "displacement_vector", displacement_metrics);
    fuelsim::test::print_relative_metrics(prefix + "equivalent_stress", stress_metrics);
    fuelsim::test::print_relative_metrics(prefix + "equivalent_plastic_strain", plastic_metrics);
    fuelsim::test::print_relative_metrics(prefix + "equivalent_creep_strain", creep_metrics);
    double abaqus_contact_force = 0.0, fuelsim_contact_force = 0.0;
    double abaqus_nodal_tangential_force_magnitude_sum = 0.0, fuelsim_integrated_tangential_traction_magnitude = 0.0;
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
            maximum_reference_coordinate_difference =
                std::max(maximum_reference_coordinate_difference, distance(point_at(reference.node), reference.point));
            for (std::size_t component = 0; component < 3; ++component) {
                resultant[component] += reference.normal_force[component];
                tangential_resultant[component] += reference.tangential_force[component];
            }
            abaqus_nodal_tangential_force_magnitude_sum +=
                std::hypot(reference.tangential_force[0], reference.tangential_force[1], reference.tangential_force[2]);
        }
        abaqus_contact_force = std::hypot(resultant[0], resultant[1], resultant[2]);
        const auto& projected = output.nodal("contact_projected_coupled_contact");
        std::vector<std::size_t> source_nodes;
        for (std::size_t node = 0; node < projected.size(); ++node)
            if (std::isfinite(projected[node])) source_nodes.push_back(node);
        if (source_nodes.size() != contact_reference.size())
            throw std::invalid_argument("C3D20T contact output size changed");
        for (std::size_t node = 0; node < source_nodes.size(); ++node) {
            const auto scalar = [&](const std::string& field) {
                const double value = output.nodal("contact_" + field + "_coupled_contact").at(source_nodes[node]);
                if (!std::isfinite(value)) throw std::invalid_argument("C3D20T contact output is not finite");
                return value;
            };
            passed = check(scalar("projected") == 1.0, "C3D20T contact constraint remains projected") && passed;
            const std::array<double, 3> slip = {scalar("total_slip_x"), scalar("total_slip_y"), scalar("total_slip_z")};
            const ContactReference& reference = contact_by_node.at(source_nodes[node]);
            if (scalar("pressure") > 0.0) {
                ++active_contact_nodes;
                ++((scalar("sliding") == 1.0) ? sliding_contact_nodes : sticking_contact_nodes);
            }
            constexpr double penalty = 1e9; // Matched tracked Abaqus penalty.
            contact_gap_metrics.add(scalar("gap"), reference.gap);
            constraint_contact_pressure_metrics.add(
                std::max(-penalty * scalar("gap"), 0.0), std::max(-penalty * reference.gap, 0.0));
            recovered_contact_pressure_metrics.add(scalar("pressure"), reference.pressure);
            std::array<double, 3> fuelsim_secondary_force{}, fuelsim_secondary_tangential_force{};
            for (std::size_t component = 0; component < fuelsim_secondary_force.size(); ++component)
                fuelsim_secondary_force[component] = -scalar("normal_force_" + std::string(1, "xyz"[component]));
            for (std::size_t component = 0; component < fuelsim_secondary_tangential_force.size(); ++component)
                fuelsim_secondary_tangential_force[component] =
                    -scalar("tangential_force_" + std::string(1, "xyz"[component]));
            contact_normal_force_metrics.add(
                fuelsim_secondary_force.data(), reference.normal_force.data(), reference.normal_force.size());
            if (frictional) {
                contact_tangential_force_metrics.add(fuelsim_secondary_tangential_force.data(),
                    reference.tangential_force.data(), reference.tangential_force.size());
                contact_tangential_slip_metrics.add(
                    slip.data(), reference.tangential_slip.data(), reference.tangential_slip.size());
            }
        }
        fuelsim_contact_force = output.global("contact_force_coupled_contact");
        fuelsim_integrated_tangential_traction_magnitude = output.global("contact_tangential_force_coupled_contact");
        std::array<double, 3> fuelsim_tangential_resultant{};
        for (const auto source : source_nodes)
            for (std::size_t component = 0; component < 3; ++component)
                fuelsim_tangential_resultant[component] -=
                    output.nodal("contact_tangential_force_" + std::string(1, "xyz"[component]) + "_coupled_contact")
                        .at(source);
        contact_tangential_resultant_metrics.add(fuelsim_tangential_resultant.data(), tangential_resultant.data(), 3);
        const double contact_tolerance = frictional ? 5.0e-3 : 1.0e-2;
        passed = check(active_contact_nodes > 0 && summary_active == active_contact_nodes,
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
              << '\n';
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

namespace fuelsim::test {
bool check_hex20_integrated(const std::string& output_path, const std::string& temperature_path,
    const std::string& displacement_path, const std::string& material_path, const std::string& contact_path,
    bool frictional, std::size_t summary_active) {
    return run(
        output_path, temperature_path, displacement_path, material_path, contact_path, frictional, summary_active);
}
} // namespace fuelsim::test
