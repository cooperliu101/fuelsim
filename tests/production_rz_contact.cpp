#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace fuelsim::test {
namespace {
bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        fields.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos)
            return fields;
        begin = separator + 1;
    }
}

std::size_t column(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end())
        throw std::invalid_argument("MOOSE contact CSV is missing column: " + name);
    return static_cast<std::size_t>(found - header.begin());
}

double value(const std::vector<std::string>& fields, std::size_t index, const std::string& path) {
    if (index >= fields.size())
        throw std::invalid_argument("MOOSE contact CSV row is incomplete: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(fields[index], &parsed);
    if (parsed != fields[index].size() || !std::isfinite(result))
        throw std::invalid_argument("MOOSE contact CSV contains an invalid number: " + path);
    return result;
}

struct ContactReference final {
    std::vector<double> coordinates;
    std::vector<double> pressures;
    std::vector<double> tangential_tractions;
    std::vector<double> radii;
    std::vector<double> axial_coordinates;
    double total_force = 0.0, total_tangential_force = 0.0;
    bool has_tangential_force = false;
};

ContactReference read_contact_reference(const std::string& path) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE contact CSV: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE contact CSV is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::size_t pressure = column(header, "contact_pressure");
    const std::size_t coordinate = column(header, "y");
    const std::size_t radius = column(header, "x");
    const std::size_t radial_displacement = column(header, "disp_x");
    const std::size_t axial_displacement = column(header, "disp_y");
    const std::size_t nodal_area = column(header, "nodal_area");
    const auto tangential_found = std::find_if(header.begin(), header.end(), [](const std::string& name) {
        constexpr const char suffix[] = "tangential_force_y";
        return name == suffix
               || (name.size() > sizeof(suffix) - 1
                   && name.compare(name.size() - (sizeof(suffix) - 1), sizeof(suffix) - 1, suffix) == 0);
    });
    const bool has_tangential_force = tangential_found != header.end();
    const std::size_t tangential_force =
        has_tangential_force ? static_cast<std::size_t>(tangential_found - header.begin()) : 0;
    std::vector<std::tuple<double, double, double, double, double, double>> values;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv(line);
        const double y = value(fields, coordinate, path);
        const double p = value(fields, pressure, path);
        values.push_back({y,
            value(fields, radius, path) + value(fields, radial_displacement, path),
            y + value(fields, axial_displacement, path),
            p,
            has_tangential_force ? value(fields, tangential_force, path) : 0.0,
            value(fields, nodal_area, path)});
    }
    if (values.empty())
        throw std::invalid_argument("MOOSE contact CSV has no values: " + path);
    std::sort(values.begin(), values.end());
    ContactReference result;
    result.coordinates.reserve(values.size());
    result.pressures.reserve(values.size());
    result.tangential_tractions.reserve(values.size());
    result.radii.reserve(values.size());
    result.axial_coordinates.reserve(values.size());
    result.has_tangential_force = has_tangential_force;
    for (const auto& entry : values) {
        result.coordinates.push_back(std::get<0>(entry));
        result.radii.push_back(std::get<1>(entry));
        result.axial_coordinates.push_back(std::get<2>(entry));
        result.pressures.push_back(std::get<3>(entry));
        const double area = std::get<5>(entry);
        if (!(area > 0.0))
            throw std::invalid_argument("MOOSE contact CSV contains a nonpositive nodal area: " + path);
        result.tangential_tractions.push_back(std::get<4>(entry) / area);
    }
    for (std::size_t node = 0; node + 1 < result.coordinates.size(); ++node) {
        const double dr = result.radii[node + 1] - result.radii[node];
        const double dz = result.axial_coordinates[node + 1] - result.axial_coordinates[node];
        const double edge_length = std::hypot(dr, dz);
        const double area_node =
            2.0 * pi * 0.5 * edge_length * (2.0 * result.radii[node] + result.radii[node + 1]) / 3.0;
        const double area_next =
            2.0 * pi * 0.5 * edge_length * (2.0 * result.radii[node + 1] + result.radii[node]) / 3.0;
        result.total_force += result.pressures[node] * area_node + result.pressures[node + 1] * area_next;
        result.total_tangential_force +=
            result.tangential_tractions[node] * area_node + result.tangential_tractions[node + 1] * area_next;
    }
    return result;
}
} // namespace

bool check_rz_multi_contact(const std::string& path,
    const std::string& first_reference,
    const std::string& second_reference,
    double friction) {
    const auto result = read_final_exodus_results(path);
    bool passed = true;
    const std::array<std::string, 2> names = {"pellet_to_inner_clad", "inner_to_outer_clad"};
    const std::array<std::string, 2> references = {first_reference, second_reference};
    for (std::size_t pair = 0; pair < names.size(); ++pair) {
        const auto reference = read_contact_reference(references[pair]);
        const auto& projected = result.nodal("contact_projected_" + names[pair]);
        const auto& pressure = result.nodal("contact_pressure_" + names[pair]);
        const auto& traction = result.nodal("contact_tangential_traction_" + names[pair]);
        const auto& sliding = result.nodal("contact_sliding_" + names[pair]);
        std::vector<std::pair<double, std::size_t>> nodes;
        for (std::size_t node = 0; node < projected.size(); ++node)
            if (!std::isnan(projected[node]))
                nodes.emplace_back(result.nodes[node][1], node);
        std::sort(nodes.begin(), nodes.end());
        if (nodes.empty() || nodes.size() != reference.coordinates.size())
            throw std::invalid_argument("RZ multi-contact reference node count differs");
        FieldErrorMetrics pressure_error, traction_error;
        double maximum_pressure = 0.0, maximum_excess = 0.0;
        std::size_t sliding_count = 0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            const auto node = nodes[i].second;
            if (std::abs(nodes[i].first - reference.coordinates[i]) >= 1.0e-12)
                throw std::invalid_argument("RZ contact reference coordinates differ");
            passed = check(projected[node] == 1.0 && pressure[node] > 0.0,
                         "Every node of contact pair " + names[pair] + " is projected and active")
                     && passed;
            pressure_error.add(pressure[node], reference.pressures[i]);
            maximum_pressure = std::max(maximum_pressure, pressure[node]);
            if (friction > 0.0) {
                traction_error.add(traction[node], reference.tangential_tractions[i]);
                maximum_excess = std::max(maximum_excess, std::abs(traction[node]) - friction * pressure[node]);
                sliding_count += sliding[node] == 1.0 ? 1U : 0U;
            }
        }
        const double normal = result.global("contact_force_" + names[pair]);
        const double tangential = result.global("contact_tangential_force_" + names[pair]);
        passed = check(relative_metrics_below(pressure_error, 1.0e-3),
                     "RZ multi-contact pressure retains three 0.1 percent gates")
                 && passed;
        passed = check(std::abs(normal - reference.total_force)
                           < 1.0e-3 * std::max({1.0, std::abs(normal), reference.total_force}),
                     "RZ multi-contact total normal force agrees")
                 && passed;
        if (friction > 0.0) {
            passed = check(reference.has_tangential_force && relative_metrics_below(traction_error, 1.0e-3)
                               && sliding_count > 0 && maximum_excess <= 1.0e-12 * std::max(1.0, maximum_pressure),
                         "RZ multi-contact retains tangential field, sliding and Coulomb-cap gates")
                     && passed;
            passed =
                check(std::abs(tangential - reference.total_tangential_force)
                          < 1.0e-3 * std::max({1.0, std::abs(tangential), std::abs(reference.total_tangential_force)}),
                    "RZ multi-contact total tangential force agrees")
                && passed;
            print_relative_metrics(names[pair] + "_tangential_traction", traction_error);
        }
        print_relative_metrics(names[pair] + "_pressure", pressure_error);
        std::cout << names[pair] << "_normal_force=" << normal << " reference=" << reference.total_force
                  << " tangential_force=" << tangential << " reference=" << reference.total_tangential_force << '\n';
    }
    return passed;
}
} // namespace fuelsim::test
