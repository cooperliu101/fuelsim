#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct DisplacementReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> displacement;
};

struct ContactReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) result.push_back(field);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("H20.35 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("H20.35 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("H20.35 CSV contains an invalid node index: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<DisplacementReference> read_displacements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.35 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "id,x,y,z,disp_x,disp_y,disp_z")
        throw std::invalid_argument("Unexpected H20.35 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.35 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2")
        throw std::invalid_argument("Unexpected H20.35 contact header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
                {number(values, 7, path), number(values, 8, path), number(values, 9, path)}, number(values, 10, path),
                number(values, 11, path)});
    }
    return result;
}

double coordinate_difference(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::max({std::abs(a[0] - b[0]), std::abs(a[1] - b[1]), std::abs(a[2] - b[2])});
}

double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

struct ContactValues final {
    std::array<double, 3> normal_contact_force, tangential_contact_force, tangential_slip;
};

struct ContactCounts final {
    std::size_t active_contact_nodes = 0, unprojected_contact_nodes = 0;
};

bool metric_passes(const std::string& name, const fuelsim::test::FieldErrorMetrics& metric, double relative_tolerance,
    double zero_tolerance) {
    if (metric.has_relative_norm()) {
        fuelsim::test::print_relative_metrics(name, metric);
        return check(fuelsim::test::relative_metrics_below(metric, relative_tolerance) &&
                         metric.maximum_zero_reference_difference < zero_tolerance,
            name + " passes all three relative metrics and its separate zero-reference check");
    }
    fuelsim::test::print_absolute_metrics(name, metric);
    return check(metric.maximum_zero_reference_difference < zero_tolerance,
        name + " passes its theoretical-zero absolute check");
}

} // namespace

namespace fuelsim::test {
bool check_hex20_curved_friction(
    const std::string& output_path, const std::string& displacement_path, const std::string& contact_path) {
    const auto output = read_final_exodus_results(output_path);
    const auto displacement = read_displacements(displacement_path);
    const auto contact = read_contact(contact_path);
    std::array<FieldErrorMetrics, 3> displacement_error;
    GroupedFieldErrorMetrics displacement_vector_error;
    std::vector<bool> present(output.nodes.size(), false);
    double maximum_coordinate_error = 0.0;
    const std::array<std::string, 3> names = {"displacement_x", "displacement_y", "displacement_z"};
    for (const auto& row : displacement) {
        if (row.id >= present.size() || present[row.id]) throw std::runtime_error("H20.35 invalid nodal mapping");
        present[row.id] = true;
        maximum_coordinate_error =
            std::max(maximum_coordinate_error, coordinate_difference(output.nodes[row.id], row.point));
        std::array<double, 3> actual{};
        for (std::size_t component = 0; component < 3; ++component) {
            actual[component] = output.nodal(names[component]).at(row.id);
            displacement_error[component].add(actual[component], row.displacement[component]);
        }
        displacement_vector_error.add(actual.data(), row.displacement.data(), 3);
    }
    std::vector<ContactValues> summaries;
    std::vector<std::size_t> source_nodes;
    ContactCounts interface;
    bool all_sticking = true, passed = true;
    const auto& projected = output.nodal("contact_projected_interface");
    for (std::size_t node = 0; node < projected.size(); ++node) {
        if (std::isnan(projected[node])) continue;
        source_nodes.push_back(node);
        ContactValues values;
        for (std::size_t component = 0; component < 3; ++component) {
            const std::string suffix = std::string(1, "xyz"[component]) + "_interface";
            values.normal_contact_force[component] = output.nodal("contact_normal_force_" + suffix)[node];
            values.tangential_contact_force[component] = output.nodal("contact_tangential_force_" + suffix)[node];
            values.tangential_slip[component] = output.nodal("contact_total_slip_" + suffix)[node];
        }
        summaries.push_back(values);
        if (projected[node] != 1.0) ++interface.unprojected_contact_nodes;
        if (output.nodal("contact_pressure_interface")[node] > 0.0) ++interface.active_contact_nodes;
        all_sticking = all_sticking && output.nodal("contact_sliding_interface")[node] == 0.0;
    }
    if (summaries.size() != contact.size()) throw std::runtime_error("H20.35 contact node counts differ");
    std::vector<bool> contact_seen(output.nodes.size(), false);
    for (const auto& row : contact) {
        if (row.id >= contact_seen.size() || contact_seen[row.id])
            throw std::runtime_error("H20.35 invalid contact ID");
        contact_seen[row.id] = true;
    }
    fuelsim::test::FieldErrorMetrics normal_force, circumferential_force, axial_force, slip_1, slip_2;
    double actual_normal_resultant = 0.0, reference_normal_resultant = 0.0;
    double actual_circumferential_resultant = 0.0, reference_circumferential_resultant = 0.0;
    double actual_axial_resultant = 0.0, reference_axial_resultant = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const auto reference = std::find_if(contact.begin(), contact.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[node]; });
        if (reference == contact.end()) throw std::invalid_argument("H20.35 contact-node mapping is incomplete");
        maximum_coordinate_error = std::max(
            maximum_coordinate_error, coordinate_difference(output.nodes[source_nodes[node]], reference->point));
        const double radius = std::hypot(reference->point[0], reference->point[1]);
        const std::array<double, 3> radial = {reference->point[0] / radius, reference->point[1] / radius, 0.0};
        const std::array<double, 3> circumferential = {
            -reference->point[1] / radius, reference->point[0] / radius, 0.0};
        const std::array<double, 3> actual_normal = {-summaries[node].normal_contact_force[0],
            -summaries[node].normal_contact_force[1], -summaries[node].normal_contact_force[2]};
        const std::array<double, 3> actual_tangent = {-summaries[node].tangential_contact_force[0],
            -summaries[node].tangential_contact_force[1], -summaries[node].tangential_contact_force[2]};
        const double actual_normal_value = dot(actual_normal, radial);
        const double reference_normal_value = dot(reference->normal_force, radial);
        const double actual_circumferential = dot(actual_tangent, circumferential);
        const double reference_circumferential = dot(reference->tangential_force, circumferential);
        normal_force.add(actual_normal_value, reference_normal_value);
        circumferential_force.add(actual_circumferential, reference_circumferential);
        axial_force.add(actual_tangent[2], reference->tangential_force[2]);
        // On the C3D20 S6 faces in this tracked mesh, Abaqus local direction 1
        // is opposite the positive circumferential direction and local
        // direction 2 is opposite the global axial direction.
        slip_1.add(-dot(summaries[node].tangential_slip, circumferential), reference->slip_1);
        slip_2.add(-summaries[node].tangential_slip[2], reference->slip_2);
        actual_normal_resultant += actual_normal_value;
        reference_normal_resultant += reference_normal_value;
        actual_circumferential_resultant += actual_circumferential;
        reference_circumferential_resultant += reference_circumferential;
        actual_axial_resultant += actual_tangent[2];
        reference_axial_resultant += reference->tangential_force[2];
    }
    constexpr double relative_tolerance = 1.0e-2;
    constexpr double zero_tolerance = 1.0e-10;
    for (std::size_t component = 0; component < 3; ++component)
        fuelsim::test::print_relative_metrics(
            "h20_35_displacement_" + std::string(1, "xyz"[component]), displacement_error[component]);
    fuelsim::test::print_grouped_relative_metrics("h20_35_displacement_vector", displacement_vector_error);
    passed = check(fuelsim::test::grouped_relative_metrics_below(displacement_vector_error, relative_tolerance) &&
                       displacement_vector_error.maximum_zero_reference_difference < zero_tolerance,
                 "H20.35 complete displacement-vector metrics and its separate zero-reference check are below 1 "
                 "percent") &&
             passed;
    passed =
        metric_passes("h20_35_signed_normal_radial_force", normal_force, relative_tolerance, zero_tolerance) && passed;
    passed = metric_passes("h20_35_signed_tangential_circumferential_force", circumferential_force, relative_tolerance,
                 zero_tolerance) &&
             passed;
    passed = metric_passes("h20_35_signed_tangential_axial_force", axial_force, relative_tolerance, zero_tolerance) &&
             passed;
    passed = metric_passes("h20_35_tangential_slip_1", slip_1, relative_tolerance, zero_tolerance) && passed;
    passed = metric_passes("h20_35_tangential_slip_2", slip_2, relative_tolerance, zero_tolerance) && passed;
    const double normal_resultant_error =
        std::abs(actual_normal_resultant - reference_normal_resultant) / std::abs(reference_normal_resultant);
    const double circumferential_resultant_error =
        std::abs(actual_circumferential_resultant - reference_circumferential_resultant) /
        std::abs(reference_circumferential_resultant);
    const double axial_resultant_error =
        std::abs(actual_axial_resultant - reference_axial_resultant) / std::abs(reference_axial_resultant);
    std::cout << "h20_35_normal_resultant_relative_error=" << normal_resultant_error << '\n'
              << "h20_35_circumferential_resultant_relative_error=" << circumferential_resultant_error << '\n'
              << "h20_35_axial_resultant_relative_error=" << axial_resultant_error << '\n'
              << '\n';
    return check(displacement.size() == output.nodes.size() &&
                     std::all_of(present.begin(), present.end(), [](bool value) { return value; }),
               "H20.35 compares every tracked Exodus mesh node") &&
           check(maximum_coordinate_error < 1.0e-12, "H20.35 references preserve the tracked Exodus coordinates") &&
           check(interface.active_contact_nodes == summaries.size() && interface.unprojected_contact_nodes == 0,
               "H20.35 keeps all 37 quadratic-surface contact nodes active and projected") &&
           check(summaries.size() == 37 && all_sticking,
               "H20.35 keeps every constraint in the two-direction sticking branch") &&
           check(std::abs(reference_circumferential_resultant) > 1.0 && std::abs(reference_axial_resultant) > 1.0,
               "H20.35 activates nonzero circumferential and axial tangential resultants") &&
           check(normal_resultant_error < relative_tolerance && circumferential_resultant_error < relative_tolerance &&
                     axial_resultant_error < relative_tolerance,
               "H20.35 signed normal and both tangential resultants agree with Abaqus below 1 percent") &&
           passed;
}
} // namespace fuelsim::test
