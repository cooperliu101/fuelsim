#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
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
    if (index >= values.size()) throw std::invalid_argument("Incomplete reference row in " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[index], &parsed);
    if (parsed != values[index].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid reference number in " + path);
    return result;
}

std::size_t identifier(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    const double result = number(values, index, path);
    if (result < 0.0 || std::floor(result) != result)
        throw std::invalid_argument("Invalid reference identifier in " + path);
    return static_cast<std::size_t>(result);
}

using fuelsim::test::FieldErrorMetrics;
using fuelsim::test::print_relative_metrics;
using fuelsim::test::relative_metrics_below;

struct RzNodalReference final {
    double radius;
    double axial_coordinate;
    double temperature;
    double radial_displacement;
    double axial_displacement;
};

struct RzNodalComparison final {
    FieldErrorMetrics temperature;
    FieldErrorMetrics radial_displacement;
    FieldErrorMetrics axial_displacement;
    double maximum_coordinate_difference = 0.0;
    std::size_t node_count = 0;
};

std::size_t column_index(const std::vector<std::string>& header, const std::string& name, const std::string& path) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("Reference CSV is missing column '" + name + "': " + path);
    return static_cast<std::size_t>(found - header.begin());
}

std::vector<RzNodalReference> read_rz_nodal_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE nodal reference is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::array<std::size_t, 6> columns = {column_index(header, "T", path), column_index(header, "disp_x", path),
        column_index(header, "disp_y", path), column_index(header, "id", path), column_index(header, "x", path),
        column_index(header, "y", path)};
    std::vector<RzNodalReference> result;
    std::vector<bool> present;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv(line);
        const std::size_t id = identifier(fields, columns[3], path);
        if (id >= result.size()) {
            result.resize(id + 1);
            present.resize(id + 1, false);
        }
        const RzNodalReference row = {number(fields, columns[4], path), number(fields, columns[5], path),
            number(fields, columns[0], path), number(fields, columns[1], path), number(fields, columns[2], path)};
        if (present[id]) {
            const RzNodalReference& previous = result[id];
            if (previous.radius != row.radius || previous.axial_coordinate != row.axial_coordinate ||
                previous.temperature != row.temperature || previous.radial_displacement != row.radial_displacement ||
                previous.axial_displacement != row.axial_displacement)
                throw std::invalid_argument("MOOSE emitted inconsistent duplicate values for a shared RZ node");
        } else {
            result[id] = row;
            present[id] = true;
        }
    }
    if (result.empty() || std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("MOOSE nodal reference IDs must be contiguous: " + path);
    return result;
}

RzNodalComparison compare_rz_nodal_results(
    const fuelsim::test::ExodusResults& results, const std::vector<RzNodalReference>& reference) {
    if (results.nodes.size() != reference.size())
        throw std::invalid_argument("MOOSE and production-result node counts differ");
    const std::vector<double>& temperature = results.nodal("temperature");
    const std::vector<double>& radial_displacement = results.nodal("displacement_r");
    const std::vector<double>& axial_displacement = results.nodal("displacement_z");
    RzNodalComparison comparison;
    comparison.node_count = results.nodes.size();
    for (std::size_t node = 0; node < results.nodes.size(); ++node) {
        comparison.maximum_coordinate_difference = std::max(
            {comparison.maximum_coordinate_difference, std::abs(results.nodes[node][0] - reference[node].radius),
                std::abs(results.nodes[node][1] - reference[node].axial_coordinate)});
        comparison.temperature.add(temperature[node], reference[node].temperature);
        comparison.radial_displacement.add(radial_displacement[node], reference[node].radial_displacement);
        comparison.axial_displacement.add(axial_displacement[node], reference[node].axial_displacement);
    }
    return comparison;
}

std::vector<double> read_contact_pressure_reference(
    const std::string& path, std::vector<double>& coordinates, std::size_t coordinate_component) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE pressure reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE pressure reference is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::size_t pressure_column = column_index(header, "contact_pressure", path);
    const std::size_t coordinate_column = column_index(header, coordinate_component == 0 ? "x" : "y", path);
    std::vector<std::pair<double, double>> values;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv(line);
        values.emplace_back(number(fields, coordinate_column, path), number(fields, pressure_column, path));
    }
    std::sort(values.begin(), values.end());
    coordinates.clear();
    std::vector<double> result;
    for (const auto& value : values) {
        coordinates.push_back(value.first);
        result.push_back(value.second);
    }
    if (result.empty()) throw std::invalid_argument("MOOSE pressure reference has no values: " + path);
    return result;
}

bool run_m0(const std::string& results_path, const std::string& nodal_reference_path) {
    const std::vector<RzNodalReference> reference = read_rz_nodal_reference(nodal_reference_path);
    const RzNodalComparison fields =
        compare_rz_nodal_results(fuelsim::test::read_final_exodus_results(results_path), reference);
    constexpr double tolerance = 1.0e-10;
    bool passed = check(fields.node_count == reference.size() && fields.maximum_coordinate_difference < 1.0e-12,
        "M0 compares every production Exodus node at matching coordinates");
    passed =
        check(relative_metrics_below(fields.temperature, tolerance), "M0 full-field temperature three errors pass") &&
        passed;
    passed = check(relative_metrics_below(fields.radial_displacement, tolerance),
                 "M0 full-field radial displacement three errors pass") &&
             passed;
    passed = check(relative_metrics_below(fields.axial_displacement, tolerance),
                 "M0 full-field axial displacement three errors pass") &&
             passed;
    print_relative_metrics("m0_temperature", fields.temperature);
    print_relative_metrics("m0_radial_displacement", fields.radial_displacement);
    print_relative_metrics("m0_axial_displacement", fields.axial_displacement);
    return passed;
}

std::map<std::string, std::string> read_summary(const std::string& path);
double summary_number(const std::map<std::string, std::string>& summary, const std::string& name);

bool run_m21(
    const std::string& results_path, const std::string& summary_path, const std::string& nodal_reference_path) {
    const std::vector<RzNodalReference> reference = read_rz_nodal_reference(nodal_reference_path);
    const RzNodalComparison fields =
        compare_rz_nodal_results(fuelsim::test::read_final_exodus_results(results_path), reference);
    constexpr double temperature_tolerance = 1.0e-3;
    constexpr double zero_displacement_tolerance = 1.0e-12;
    bool passed = check(fields.node_count == reference.size() && fields.maximum_coordinate_difference < 1.0e-12,
        "M2.1 compares every production Exodus node at matching coordinates");
    passed = check(relative_metrics_below(fields.temperature, temperature_tolerance),
                 "M2.1 full-field temperature three errors pass 0.1 percent") &&
             passed;
    passed = check(fields.radial_displacement.maximum_absolute_difference < zero_displacement_tolerance &&
                       fields.radial_displacement.maximum_actual < zero_displacement_tolerance,
                 "M2.1 radial displacement remains zero") &&
             passed;
    passed = check(fields.axial_displacement.maximum_absolute_difference < zero_displacement_tolerance &&
                       fields.axial_displacement.maximum_actual < zero_displacement_tolerance,
                 "M2.1 axial displacement remains zero") &&
             passed;
    print_relative_metrics("m21_temperature", fields.temperature);
    std::cout << "m21_radial_displacement_maximum_absolute=" << fields.radial_displacement.maximum_absolute_difference
              << '\n';
    std::cout << "m21_axial_displacement_maximum_absolute=" << fields.axial_displacement.maximum_absolute_difference
              << '\n';
    const std::map<std::string, std::string> summary = read_summary(summary_path);
    passed = check(summary_number(summary, "accepted_steps") == 10.0,
                 "M2.1 production transient completes ten accepted steps") &&
             passed;
    passed = check(summary_number(summary, "petsc_workspace_setups") == 1.0,
                 "M2.1 production transient reuses one PETSc workspace") &&
             passed;
    return passed;
}

std::map<std::string, std::string> read_summary(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read fuelsim summary: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "metric,value") throw std::invalid_argument("Unexpected fuelsim summary header in " + path);
    std::map<std::string, std::string> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t separator = line.find(',');
        if (separator == std::string::npos || separator == 0 || separator + 1 == line.size())
            throw std::invalid_argument("Invalid fuelsim summary row in " + path);
        if (!result.emplace(line.substr(0, separator), line.substr(separator + 1)).second)
            throw std::invalid_argument("Duplicate fuelsim summary metric in " + path);
    }
    return result;
}

double summary_number(const std::map<std::string, std::string>& summary, const std::string& name) {
    const auto found = summary.find(name);
    if (found == summary.end()) throw std::invalid_argument("fuelsim summary is missing metric '" + name + "'");
    return number({found->second}, 0, "fuelsim summary metric " + name);
}

FieldErrorMetrics compare_contact_pressure(const fuelsim::test::ExodusResults& results,
    const std::string& variable_name, const std::string& reference_path, double coordinate_tolerance,
    std::size_t coordinate_component = 1) {
    std::vector<double> reference_coordinates;
    const std::vector<double> reference =
        read_contact_pressure_reference(reference_path, reference_coordinates, coordinate_component);
    const std::vector<double>& pressure = results.nodal(variable_name);
    std::vector<std::pair<double, double>> actual;
    for (std::size_t node = 0; node < pressure.size(); ++node)
        if (std::isfinite(pressure[node]))
            actual.emplace_back(results.nodes[node][coordinate_component], pressure[node]);
    std::sort(actual.begin(), actual.end());
    if (actual.size() != reference.size())
        throw std::invalid_argument("MOOSE and production-result contact-pressure counts differ");
    FieldErrorMetrics metrics;
    for (std::size_t node = 0; node < actual.size(); ++node) {
        if (std::abs(actual[node].first - reference_coordinates[node]) > coordinate_tolerance)
            throw std::invalid_argument("MOOSE and production-result contact coordinates differ");
        metrics.add(actual[node].second, reference[node]);
    }
    return metrics;
}

bool run_m1(const std::string& name, const std::string& results_path, const std::string& summary_path,
    const std::string& nodal_reference_path, const std::string& pressure_reference_path) {
    const fuelsim::test::ExodusResults results = fuelsim::test::read_final_exodus_results(results_path);
    const std::vector<RzNodalReference> reference = read_rz_nodal_reference(nodal_reference_path);
    const RzNodalComparison fields = compare_rz_nodal_results(results, reference);
    const FieldErrorMetrics pressure =
        compare_contact_pressure(results, "contact_pressure_fuel_cladding", pressure_reference_path, 1.0e-12);
    constexpr double tolerance = 1.0e-2;
    bool passed = check(fields.node_count == reference.size() && fields.maximum_coordinate_difference < 1.0e-12,
        name + " compares every production Exodus node at matching coordinates");
    passed = check(relative_metrics_below(fields.temperature, tolerance) &&
                       relative_metrics_below(fields.radial_displacement, tolerance) &&
                       relative_metrics_below(fields.axial_displacement, tolerance) &&
                       relative_metrics_below(pressure, tolerance),
                 name + " full-field three-metric errors pass 1 percent") &&
             passed;
    print_relative_metrics(name + "_temperature", fields.temperature);
    print_relative_metrics(name + "_radial_displacement", fields.radial_displacement);
    print_relative_metrics(name + "_axial_displacement", fields.axial_displacement);
    print_relative_metrics(name + "_contact_pressure", pressure);

    const std::map<std::string, std::string> summary = read_summary(summary_path);
    passed = check(summary_number(summary, "petsc_workspace_setups") == 1.0,
                 name + " production path reuses one PETSc workspace") &&
             passed;
    if (name == "m1") {
        FieldErrorMetrics total_force;
        total_force.add(summary_number(summary, "contact.fuel_cladding.total_contact_force"), 663.8896691615588);
        print_relative_metrics("m1_total_contact_force", total_force);
        passed =
            check(relative_metrics_below(total_force, tolerance), "M1 total contact force three errors pass") && passed;
        passed = check(summary_number(summary, "contact.fuel_cladding.projected_contact_nodes") == 11.0 &&
                           summary_number(summary, "contact.fuel_cladding.active_contact_nodes") == 11.0,
                     "M1 projects and activates all fuel-surface nodes") &&
                 passed;
    }
    return passed;
}

struct CartesianNodeReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 4> fields;
};

struct CartesianStressReference final {
    std::size_t id;
    std::array<double, 6> stress;
};

std::vector<CartesianNodeReference> read_cartesian_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "T,disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional MOOSE nodal header in " + path);
    std::vector<CartesianNodeReference> result;
    std::vector<bool> present;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const CartesianNodeReference row = {identifier(values, 4, path),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path), number(values, 3, path)}};
        if (row.id >= result.size()) {
            result.resize(row.id + 1);
            present.resize(row.id + 1, false);
        }
        if (present[row.id]) {
            const CartesianNodeReference& previous = result[row.id];
            if (previous.point != row.point || previous.fields != row.fields)
                throw std::invalid_argument("MOOSE emitted inconsistent duplicate values for a shared node");
        } else {
            result[row.id] = row;
            present[row.id] = true;
        }
    }
    if (result.empty() || std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Three-dimensional MOOSE node IDs must be contiguous");
    return result;
}

std::vector<CartesianStressReference> read_cartesian_stresses(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE stresses: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "id,stress_xx,stress_xy,stress_xz,stress_yy,stress_yz,stress_zz,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional MOOSE stress header in " + path);
    std::vector<CartesianStressReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back({identifier(values, 0, path),
            {number(values, 1, path), number(values, 4, path), number(values, 6, path), number(values, 2, path),
                number(values, 5, path), number(values, 3, path)}});
    }
    return result;
}

bool run_b3(const std::string& results_path, const std::string& nodal_path, const std::string& stress_path) {
    const fuelsim::test::ExodusResults results = fuelsim::test::read_final_exodus_results(results_path);
    const std::vector<CartesianNodeReference> nodes = read_cartesian_nodes(nodal_path);
    const std::vector<CartesianStressReference> stresses = read_cartesian_stresses(stress_path);
    if (nodes.size() != results.nodes.size())
        throw std::invalid_argument("Three-dimensional MOOSE and fuelsim node counts differ");
    const std::array<std::string, 4> field_names = {
        "temperature", "displacement_x", "displacement_y", "displacement_z"};
    std::array<FieldErrorMetrics, 4> field_errors;
    double coordinate_error = 0.0;
    std::vector<bool> visited(results.nodes.size(), false);
    for (const CartesianNodeReference& reference : nodes) {
        if (reference.id >= results.nodes.size() || visited[reference.id])
            throw std::invalid_argument("Three-dimensional reference node IDs must be unique and in range");
        visited[reference.id] = true;
        for (std::size_t component = 0; component < 3; ++component)
            coordinate_error = std::max(
                coordinate_error, std::abs(results.nodes[reference.id][component] - reference.point[component]));
        for (std::size_t field = 0; field < field_names.size(); ++field)
            field_errors[field].add(results.nodal(field_names[field]).at(reference.id), reference.fields[field]);
    }
    if (std::find(visited.begin(), visited.end(), false) != visited.end())
        throw std::invalid_argument("Three-dimensional reference does not cover every production result node");

    const std::array<std::string, 6> stress_names = {"xx", "yy", "zz", "xy", "yz", "xz"};
    std::array<FieldErrorMetrics, 6> stress_errors;
    for (const CartesianStressReference& reference : stresses) {
        for (std::size_t point = 0; point < 8; ++point)
            for (std::size_t component = 0; component < stress_names.size(); ++component)
                stress_errors[component].add(
                    results.element("stress_" + stress_names[component] + "_q" + std::to_string(point))
                        .at(reference.id),
                    reference.stress[component]);
    }

    bool passed = true;
    for (std::size_t field = 0; field < field_errors.size(); ++field) {
        print_relative_metrics("b3_" + field_names[field], field_errors[field]);
        passed = check(relative_metrics_below(field_errors[field], 1.0e-3) &&
                           field_errors[field].maximum_zero_reference_difference < 1.0e-10,
                     "stage B " + field_names[field] + " three metrics are below 0.1 percent") &&
                 passed;
    }
    print_relative_metrics("b3_stress_xx", stress_errors[0]);
    passed = check(relative_metrics_below(stress_errors[0], 1.0e-3),
                 "stage B nonzero stress three metrics are below 0.1 percent") &&
             passed;
    for (std::size_t component = 1; component < stress_errors.size(); ++component)
        passed = check(stress_errors[component].maximum_absolute_difference < 1.0e-6,
                     "stage B near-zero stress component satisfies its absolute tolerance") &&
                 passed;
    return check(coordinate_error < 1.0e-12, "stage B compares all nodes at matching coordinates") && passed;
}

bool run_cartesian_fields(
    const std::string& name, const std::string& results_path, const std::string& nodal_path, double tolerance) {
    const fuelsim::test::ExodusResults results = fuelsim::test::read_final_exodus_results(results_path);
    const std::vector<CartesianNodeReference> nodes = read_cartesian_nodes(nodal_path);
    if (nodes.size() != results.nodes.size())
        throw std::invalid_argument("Three-dimensional MOOSE and production-result node counts differ");
    const std::array<std::string, 4> field_names = {
        "temperature", "displacement_x", "displacement_y", "displacement_z"};
    std::array<FieldErrorMetrics, 4> errors;
    double coordinate_error = 0.0;
    for (const CartesianNodeReference& reference : nodes) {
        for (std::size_t component = 0; component < 3; ++component)
            coordinate_error = std::max(
                coordinate_error, std::abs(results.nodes[reference.id][component] - reference.point[component]));
        for (std::size_t field = 0; field < field_names.size(); ++field)
            errors[field].add(results.nodal(field_names[field]).at(reference.id), reference.fields[field]);
    }
    bool passed = true;
    for (std::size_t field = 0; field < errors.size(); ++field) {
        print_relative_metrics(name + "_" + field_names[field], errors[field]);
        passed = check(relative_metrics_below(errors[field], tolerance) &&
                           errors[field].maximum_zero_reference_difference < 1.0e-10,
                     name + " full-field three-metric errors and zero-reference absolute errors pass") &&
                 passed;
    }
    return check(coordinate_error < 1.0e-12, name + " compares every production Exodus node at matching coordinates") &&
           passed;
}

bool run_rz_fields(
    const std::string& name, const std::string& results_path, const std::string& nodal_path, double tolerance) {
    const std::vector<RzNodalReference> reference = read_rz_nodal_reference(nodal_path);
    const RzNodalComparison fields =
        compare_rz_nodal_results(fuelsim::test::read_final_exodus_results(results_path), reference);
    bool passed = check(fields.node_count == reference.size() && fields.maximum_coordinate_difference < 1.0e-12,
        name + " compares every production Exodus node at matching coordinates");
    for (const auto& field : {std::pair<std::string, const FieldErrorMetrics&>{"temperature", fields.temperature},
             {"radial_displacement", fields.radial_displacement}, {"axial_displacement", fields.axial_displacement}}) {
        print_relative_metrics(name + "_" + field.first, field.second);
        passed = check(relative_metrics_below(field.second, tolerance) &&
                           field.second.maximum_zero_reference_difference < 1.0e-10,
                     name + " full-field three-metric errors and zero-reference absolute errors pass") &&
                 passed;
    }
    return passed;
}

struct MixedOrderReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> values;
};

std::vector<MixedOrderReference> read_mixed_order_reference(const std::string& path, bool temperature) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read HEX20 MOOSE reference: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected = temperature ? "T,id,x,y,z" : "disp_x,disp_y,disp_z,id,x,y,z";
    if (line != expected) throw std::invalid_argument("Unexpected HEX20 MOOSE reference header: " + path);
    std::vector<MixedOrderReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t offset = temperature ? 1U : 3U;
        MixedOrderReference row{};
        row.id = identifier(values, offset, path);
        row.point = {
            number(values, offset + 1, path), number(values, offset + 2, path), number(values, offset + 3, path)};
        if (temperature)
            row.values[0] = number(values, 0, path);
        else
            for (std::size_t component = 0; component < 3; ++component)
                row.values[component] = number(values, component, path);
        result.push_back(row);
    }
    return result;
}

bool run_hex20_fields(const std::string& name, const std::string& results_path, const std::string& temperature_path,
    const std::string& displacement_path, double tolerance, bool compare_displacement) {
    const fuelsim::test::ExodusResults results = fuelsim::test::read_final_exodus_results(results_path);
    const std::vector<MixedOrderReference> temperature = read_mixed_order_reference(temperature_path, true);
    const std::vector<MixedOrderReference> displacement = read_mixed_order_reference(displacement_path, false);
    if (temperature.size() != 8 || displacement.size() != 20 || results.nodes.size() != 20)
        throw std::invalid_argument("HEX20 MOOSE comparison has the wrong mixed-order node counts");
    FieldErrorMetrics temperature_error;
    std::array<FieldErrorMetrics, 3> displacement_error;
    double coordinate_error = 0.0;
    for (const MixedOrderReference& row : temperature) {
        if (row.id >= results.nodes.size()) throw std::invalid_argument("HEX20 temperature node ID is out of range");
        for (std::size_t component = 0; component < 3; ++component)
            coordinate_error =
                std::max(coordinate_error, std::abs(results.nodes[row.id][component] - row.point[component]));
        temperature_error.add(results.nodal("temperature").at(row.id), row.values[0]);
    }
    const std::array<std::string, 3> displacement_names = {"displacement_x", "displacement_y", "displacement_z"};
    for (const MixedOrderReference& row : displacement) {
        if (row.id >= results.nodes.size()) throw std::invalid_argument("HEX20 displacement node ID is out of range");
        for (std::size_t component = 0; component < 3; ++component) {
            coordinate_error =
                std::max(coordinate_error, std::abs(results.nodes[row.id][component] - row.point[component]));
            displacement_error[component].add(
                results.nodal(displacement_names[component]).at(row.id), row.values[component]);
        }
    }
    print_relative_metrics(name + "_temperature", temperature_error);
    bool passed = check(relative_metrics_below(temperature_error, name == "b6_hex20" ? 1.0e-10 : tolerance),
        "HEX20 temperature relative L2, relative absolute peak, and maximum pointwise errors pass");
    if (compare_displacement) {
        for (std::size_t component = 0; component < displacement_error.size(); ++component) {
            print_relative_metrics(name + "_" + displacement_names[component], displacement_error[component]);
            passed = check(relative_metrics_below(displacement_error[component], tolerance) &&
                               displacement_error[component].maximum_zero_reference_difference < 1.0e-12,
                         "HEX20 " + displacement_names[component] +
                             " relative L2, relative absolute peak, and maximum pointwise errors pass") &&
                     passed;
        }
    }
    return check(coordinate_error < 1.0e-14, "HEX20 comparison uses identical tracked MOOSE mesh coordinates") &&
           passed;
}

bool run_b6(
    const std::string& results_path, const std::string& temperature_path, const std::string& displacement_path) {
    return run_hex20_fields("b6_hex20", results_path, temperature_path, displacement_path, 1.0e-8, true);
}

bool completed_summary(const std::string& path, const std::string& problem);

bool run_hex20_thermal_contact(
    const std::string& result_path, const std::string& summary_path, const std::string& reference_path) {
    const auto output = fuelsim::test::read_final_exodus_results(result_path);
    const auto reference = read_mixed_order_reference(reference_path, true);
    const auto& temperature = output.nodal("temperature");
    std::vector<bool> present(output.nodes.size(), false);
    FieldErrorMetrics error;
    double coordinate_error = 0.0;
    for (const auto& row : reference) {
        if (row.id >= present.size() || present[row.id])
            throw std::invalid_argument("H20.16 temperature node mapping is incomplete or repeated");
        present[row.id] = true;
        for (std::size_t component = 0; component < 3; ++component)
            coordinate_error =
                std::max(coordinate_error, std::abs(output.nodes[row.id][component] - row.point[component]));
        error.add(temperature.at(row.id), row.values[0]);
    }
    // Steady HEX20 Exodus output also interpolates temperature onto displacement
    // midpoint nodes. The tracked reference contains the sixteen solved corner
    // temperatures, not these visualization-only interpolated values.
    print_relative_metrics("h20_16_temperature", error);
    const double heat_rate = output.global("contact_heat_rate_interface");
    std::cout << "h20_16_total_heat_rate=" << heat_rate << '\n';
    return completed_summary(summary_path, "steady") &&
           check(summary_number(read_summary(summary_path), "load_steps_completed") == 1.0,
               "H20.16 completes one load step") &&
           check(output.nodes.size() == 40 && reference.size() == 16,
               "H20.16 compares all sixteen first-order temperature nodes") &&
           check(relative_metrics_below(error, 5.0e-3), "H20.16 temperature errors pass 0.5 percent") &&
           check(coordinate_error < 1.0e-12, "H20.16 tracked MOOSE coordinates match") &&
           check(heat_rate > 0.0, "H20.16 transfers nonzero heat");
}

bool run_hex20_transient(const std::string& results_path, const std::string& summary_path,
    const std::string& temperature_path, const std::string& displacement_path, double tolerance) {
    bool passed =
        run_hex20_fields("hex20_transient", results_path, temperature_path, displacement_path, tolerance, true);
    const std::map<std::string, std::string> summary = read_summary(summary_path);
    passed = check(summary_number(summary, "accepted_steps") == 10.0,
                 "HEX20 production transient accepts ten fixed steps") &&
             passed;
    passed = check(summary_number(summary, "petsc_workspace_setups") == 1.0,
                 "HEX20 production transient reuses one PETSc workspace") &&
             passed;
    return passed;
}

void require_argument_count(const std::string& mode, int argc, int expected) {
    if (argc != expected)
        throw std::invalid_argument("Production result check '" + mode + "' received the wrong number of arguments");
}

bool completed_summary(const std::string& path, const std::string& problem) {
    const auto summary = read_summary(path);
    return check(summary.at("completed") == "true" && summary.at("problem") == problem,
               "Production summary confirms completion of the expected problem") &&
           check(
               summary_number(summary, "petsc_workspace_setups") == 1.0, "Production solve reuses one PETSc workspace");
}

bool run_hex20_inelastic(const std::string& results_path, const std::string& summary_path,
    const std::string& temperature_path, const std::string& displacement_path, const std::string& material_path) {
    bool passed = run_hex20_transient(results_path, summary_path, temperature_path, displacement_path, 5.0e-3);
    passed = completed_summary(summary_path, "transient") && passed;
    const auto results = fuelsim::test::read_final_exodus_results(results_path);
    passed = check(std::abs(results.time - 1.0) < 1.0e-12, "HEX20 result reaches the specified final time") && passed;
    std::ifstream input(material_path);
    if (!input) throw std::runtime_error("Could not read material reference: " + material_path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("Empty material reference");
    const auto header = split_csv(line);
    if (!std::getline(input, line)) throw std::invalid_argument("Missing material reference row");
    const auto row = split_csv(line);
    if (identifier(row, column_index(header, "id", material_path), material_path) != 0)
        throw std::invalid_argument("Expected the single HEX20 element with source ID zero");
    while (std::getline(input, line))
        if (!line.empty()) throw std::invalid_argument("Expected exactly one HEX20 material reference row");
    const std::array<std::string, 8> reference_names = {"stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz",
        "stress_xz", "effective_plastic_strain", "effective_creep_strain"};
    const std::array<std::string, 8> output_names = {
        "stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "equiv_plastic", "equiv_creep"};
    for (std::size_t component = 0; component < reference_names.size(); ++component) {
        const double reference =
            number(row, column_index(header, reference_names[component], material_path), material_path);
        FieldErrorMetrics metrics;
        for (std::size_t q = 0; q < 27; ++q) {
            const auto& values = results.element(output_names[component] + "_q" + std::to_string(q));
            if (values.size() != 1) throw std::invalid_argument("Expected exactly one production HEX20 element");
            metrics.add(values[0], reference);
        }
        // These are the original H20.07-09 gates, including the transverse-stress absolute gates.
        const bool inactive = metrics.maximum_reference < 1.0e-14;
        const bool transverse = component >= 1 && component <= 5;
        if (inactive || transverse) {
            std::cout << output_names[component]
                      << "_maximum_absolute_difference=" << metrics.maximum_absolute_difference << '\n'
                      << output_names[component] << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
                      << output_names[component]
                      << "_maximum_zero_reference_difference=" << metrics.maximum_zero_reference_difference << '\n';
        } else {
            print_relative_metrics(output_names[component], metrics);
        }
        const bool matches = inactive ? metrics.maximum_absolute_difference < 1.0e-10
                                      : (transverse ? metrics.maximum_absolute_difference < 1.0
                                                    : relative_metrics_below(metrics, 5.0e-3));
        passed =
            check(matches, "All 27 HEX20 material points pass the original " + reference_names[component] + " gate") &&
            passed;
    }
    return passed;
}

bool run_hex8_inelastic(const std::string& result_path, const std::string& summary_path, const std::string& nodal_path,
    const std::string& material_path) {
    bool passed = completed_summary(summary_path, "transient") &&
                  run_cartesian_fields("hex8_inelastic", result_path, nodal_path, 5.0e-3);
    const auto result = fuelsim::test::read_final_exodus_results(result_path);
    passed = check(summary_number(read_summary(summary_path), "accepted_steps") == 10 &&
                       std::abs(result.time - 1.0) < 1.0e-12,
                 "HEX8 accepts ten steps and reaches the prescribed end time") &&
             passed;
    std::ifstream input(material_path);
    if (!input) throw std::runtime_error("Could not read HEX8 material reference");
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("Empty HEX8 material reference");
    const auto header = split_csv(line);
    const std::array<std::string, 8> names = {
        "stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "equiv_plastic", "equiv_creep"};
    std::array<FieldErrorMetrics, 8> metrics;
    std::vector<bool> present(result.element("stress_xx_q0").size(), false);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto row = split_csv(line);
        const auto element = identifier(row, column_index(header, "id", material_path), material_path);
        if (element >= present.size() || present[element])
            throw std::invalid_argument("Invalid or repeated material element");
        present[element] = true;
        for (std::size_t c = 0; c < 8; ++c) {
            const std::string reference_name =
                c < 6 ? names[c] : (c == 6 ? "effective_plastic_strain" : "effective_creep_strain");
            const double expected = number(row, column_index(header, reference_name, material_path), material_path);
            for (std::size_t q = 0; q < 8; ++q)
                metrics[c].add(result.element(names[c] + "_q" + std::to_string(q))[element], expected);
        }
    }
    passed = check(!present.empty() && std::all_of(present.begin(), present.end(), [](bool p) { return p; }),
                 "HEX8 compares all elements and all eight material points") &&
             passed;
    for (std::size_t c = 0; c < 8; ++c) {
        if (metrics[c].maximum_reference > 0.0)
            print_relative_metrics("hex8_" + names[c], metrics[c]);
        else
            std::cout << "hex8_" << names[c]
                      << "_maximum_absolute_difference=" << metrics[c].maximum_absolute_difference
                      << " zero_reference_count=" << metrics[c].zero_reference_count << '\n';
        const bool matches = c >= 1 && c <= 5 ? metrics[c].maximum_absolute_difference < 1.0
                                              : (c >= 6 && metrics[c].maximum_reference == 0.0
                                                        ? metrics[c].maximum_absolute_difference < 1.0e-14
                                                        : relative_metrics_below(metrics[c], 5.0e-3));
        passed = check(matches, "HEX8 material field retains its original MOOSE gate: " + names[c]) && passed;
    }
    return passed;
}

bool run_m54(const std::string& results_path, const std::string& summary_path, const std::string& nodal_path,
    const std::string& pressure_path) {
    bool passed = run_m1("m54", results_path, summary_path, nodal_path, pressure_path);
    passed = completed_summary(summary_path, "steady") && passed;
    const auto summary = read_summary(summary_path);
    const auto results = fuelsim::test::read_final_exodus_results(results_path);
    const auto& pressure = results.nodal("contact_pressure_fuel_cladding");
    const auto& gap = results.nodal("contact_gap_fuel_cladding");
    std::size_t count = 0, active = 0;
    double maximum_penetration = 0.0;
    for (std::size_t node = 0; node < pressure.size(); ++node) {
        if (std::isnan(pressure[node])) continue;
        if (!std::isfinite(pressure[node]) || !std::isfinite(gap[node]))
            throw std::invalid_argument("Invalid production contact field");
        ++count;
        if (pressure[node] > 0.0) ++active;
        maximum_penetration = std::max(maximum_penetration, -gap[node]);
    }
    passed =
        check(count > 0 && active == count &&
                  summary_number(summary, "contact.fuel_cladding.projected_contact_nodes") ==
                      static_cast<double>(count) &&
                  summary_number(summary, "contact.fuel_cladding.active_contact_nodes") == static_cast<double>(count),
            "M5.4 projects and activates every contact node") &&
        passed;
    passed = check(maximum_penetration <= 1.0e-9 && summary_number(summary, "augmented_lagrangian_iterations") > 0,
                 "M5.4 updates multipliers and satisfies the original 1 nm penetration tolerance") &&
             passed;
    std::cout << "m54_maximum_penetration=" << maximum_penetration << '\n';
    return passed;
}

bool run_planar_sts(const std::string& result_path, const std::string& summary_path,
    const std::string& displacement_path, const std::string& force_path, const std::string& resultant_path,
    bool hex8 = false) {
    const auto result = fuelsim::test::read_final_exodus_results(result_path);
    bool passed = completed_summary(summary_path, "steady");
    const auto summary = read_summary(summary_path);
    passed =
        check(summary_number(summary, "load_steps_completed") == 4, "Planar STS completes four load steps") && passed;
    const std::array<std::string, 3> components = {"displacement_x", "displacement_y", "displacement_z"};
    std::ifstream displacement(displacement_path), force(force_path), resultant(resultant_path);
    if (!displacement || !force || !resultant) throw std::runtime_error("Could not open planar STS references");
    std::string line;
    if (!std::getline(displacement, line) || line != "normal_displacement,normal_x,normal_y,normal_z,id,x,y,z")
        throw std::invalid_argument("Unexpected planar STS displacement header");
    FieldErrorMetrics displacement_error, force_error;
    std::vector<bool> visited(result.nodes.size(), false), force_visited(result.nodes.size(), false);
    std::array<double, 3> normal{};
    while (std::getline(displacement, line)) {
        if (line.empty()) continue;
        const auto row = split_csv(line);
        const auto node = identifier(row, 4, displacement_path);
        if (node >= visited.size() || visited[node])
            throw std::invalid_argument("Repeated or invalid displacement node");
        visited[node] = true;
        if (displacement_error.value_count == 0)
            for (std::size_t c = 0; c < 3; ++c) normal[c] = number(row, c + 1, displacement_path);
        double actual = 0.0;
        for (std::size_t c = 0; c < 3; ++c) {
            if (std::abs(result.nodes[node][c] - number(row, c + 5, displacement_path)) >= 1.0e-12)
                throw std::invalid_argument("Planar STS displacement coordinates differ");
            actual += normal[c] * result.nodal(components[c])[node];
        }
        displacement_error.add(actual, number(row, 0, displacement_path));
    }
    if (!std::getline(force, line) || line != "normal_force,id,x,y,z")
        throw std::invalid_argument("Unexpected planar STS force header");
    const auto& actual_force = result.nodal("contact_normal_force_interface");
    const auto& projected = result.nodal("contact_projected_interface");
    const auto& pressure = result.nodal("contact_pressure_interface");
    while (std::getline(force, line)) {
        if (line.empty()) continue;
        const auto row = split_csv(line);
        const auto node = identifier(row, 1, force_path);
        if (node >= force_visited.size() || force_visited[node])
            throw std::invalid_argument("Repeated or invalid force node");
        force_visited[node] = true;
        for (std::size_t c = 0; c < 3; ++c)
            if (std::abs(result.nodes[node][c] - number(row, c + 2, force_path)) >= 1.0e-12)
                throw std::invalid_argument("Planar STS force coordinates differ");
        passed = check(projected[node] == 1.0 && pressure[node] > 0.0,
                     "Every planar contact constraint is projected and active") &&
                 passed;
        force_error.add(actual_force[node], number(row, 0, force_path));
    }
    std::size_t contact_nodes = 0;
    std::array<double, 3> normal_resultant{};
    double maximum_tangential_force = 0.0;
    for (std::size_t node = 0; node < projected.size(); ++node)
        if (!std::isnan(projected[node])) {
            ++contact_nodes;
            if (!force_visited[node]) throw std::invalid_argument("Missing planar STS force reference");
            if (hex8) {
                const std::array<std::string, 3> axes = {"x", "y", "z"};
                for (std::size_t c = 0; c < 3; ++c)
                    normal_resultant[c] -= result.nodal("contact_normal_force_" + axes[c] + "_interface")[node];
                maximum_tangential_force =
                    std::max(maximum_tangential_force, result.nodal("contact_tangential_force_interface")[node]);
            }
        }
    if (hex8) {
        const double magnitude = std::hypot(normal_resultant[0], normal_resultant[1], normal_resultant[2]);
        double direction_error = 0.0;
        for (std::size_t c = 0; c < 3; ++c)
            direction_error = std::max(direction_error, std::abs(normal_resultant[c] / magnitude - normal[c]));
        passed = check(magnitude > 0.0 && direction_error < 1.0e-12 && maximum_tangential_force < 1.0e-10,
                     "B3.9 frictionless resultant follows the reference normal") &&
                 passed;
    }
    if (!std::getline(resultant, line) || line != "normal_contact_resultant,normal_outer_reaction" ||
        !std::getline(resultant, line))
        throw std::invalid_argument("Invalid planar STS resultant reference");
    const auto row = split_csv(line);
    const double expected = std::abs(number(row, 0, resultant_path));
    const double actual = result.global("contact_force_interface");
    print_relative_metrics("planar_sts_normal_displacement", displacement_error);
    print_relative_metrics("planar_sts_normal_force", force_error);
    std::cout << "planar_sts_normal_resultant=" << actual << " reference=" << expected
              << " reference_outer_reaction=" << std::abs(number(row, 1, resultant_path)) << '\n';
    return check(displacement_error.value_count == result.nodes.size() && contact_nodes > 0 &&
                     force_error.value_count == contact_nodes,
               "Planar STS covers every node and contact constraint") &&
           check(relative_metrics_below(displacement_error, 1.0e-2) && relative_metrics_below(force_error, 1.0e-2) &&
                     expected > 0 && std::abs(actual - expected) / expected < 1.0e-2,
               "Planar STS retains all original 1 percent gates") &&
           passed;
}

bool compare_result_files(
    const std::string& actual_path, const std::string& reference_path, double tolerance, bool uniaxial_mpi = false) {
    if (!std::isfinite(tolerance) || tolerance < 0.0)
        throw std::invalid_argument("Equivalence tolerance must be finite and nonnegative");
    const auto actual = fuelsim::test::read_final_exodus_results(actual_path);
    const auto reference = fuelsim::test::read_final_exodus_results(reference_path);
    if (actual.nodes != reference.nodes || actual.nodal_variable_names != reference.nodal_variable_names ||
        actual.element_variable_names != reference.element_variable_names ||
        std::abs(actual.time - reference.time) > 1.0e-12)
        throw std::invalid_argument("Production result mesh, variable schema or final time differs");
    bool passed = true;
    for (int category = 0; category < 2; ++category) {
        const auto& values = category == 0 ? actual.nodal_variables : actual.element_variables;
        const auto& expected = category == 0 ? reference.nodal_variables : reference.element_variables;
        const auto& names = category == 0 ? actual.nodal_variable_names : actual.element_variable_names;
        for (std::size_t variable = 0; variable < values.size(); ++variable) {
            if (values[variable].size() != expected[variable].size())
                throw std::invalid_argument("Result field sizes differ");
            double maximum_difference = 0.0, scale = 0.0;
            std::size_t compared = 0, missing = 0;
            for (std::size_t item = 0; item < values[variable].size(); ++item) {
                const double a = values[variable][item], b = expected[variable][item];
                if (std::isnan(a) && std::isnan(b)) {
                    ++missing;
                    continue;
                }
                if (!std::isfinite(a) || !std::isfinite(b))
                    throw std::invalid_argument("Result missing-value masks differ or contain infinity");
                maximum_difference = std::max(maximum_difference, std::abs(a - b));
                scale = std::max(scale, std::abs(b));
                ++compared;
            }
            // This MPI card is uniaxial. Transverse/shear stresses are physically zero;
            // distributed factorization leaves sub-micropascal roundoff, not a relative field error.
            const bool zero_stress =
                uniaxial_mpi && names[variable].rfind("stress_", 0) == 0 && names[variable].rfind("stress_xx_", 0) != 0;
            const bool zero_reaction =
                uniaxial_mpi && (names[variable] == "reaction_force_y" || names[variable] == "reaction_force_z");
            const bool zero_shear_strain =
                uniaxial_mpi &&
                (names[variable].rfind("elastic_", 0) == 0 || names[variable].rfind("plastic_", 0) == 0 ||
                    names[variable].rfind("creep_", 0) == 0) &&
                (names[variable].find("_xy_q") != std::string::npos ||
                    names[variable].find("_yz_q") != std::string::npos ||
                    names[variable].find("_xz_q") != std::string::npos);
            const double bound =
                zero_stress || zero_reaction
                    ? 1.0e-5
                    : (zero_shear_strain ? 1.0e-14
                                         : (tolerance == 0.0 ? 0.0 : (scale == 0.0 ? 1.0e-12 : tolerance * scale)));
            passed = check(maximum_difference <= bound, names[variable] + " production equivalence") && passed;
            std::cout << names[variable] << "_maximum_absolute_difference=" << maximum_difference
                      << " reference_absolute_peak=" << scale << " absolute_bound=" << bound << " compared=" << compared
                      << " missing=" << missing << '\n';
        }
    }
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: fuelsim_production_result_tests <check-mode> "
                     "<result.e> <summary-and-reference files...>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const std::string mode = argv[1];
        bool passed = false;
        if (mode == "m0") {
            require_argument_count(mode, argc, 4);
            passed = run_m0(argv[2], argv[3]);
        } else if (mode == "m21") {
            require_argument_count(mode, argc, 5);
            passed = run_m21(argv[2], argv[3], argv[4]);
        } else if (mode == "m1" || mode == "m1-unstructured") {
            require_argument_count(mode, argc, 6);
            passed = run_m1(mode == "m1" ? "m1" : "unstructured", argv[2], argv[3], argv[4], argv[5]);
        } else if (mode == "b3") {
            require_argument_count(mode, argc, 5);
            passed = run_b3(argv[2], argv[3], argv[4]);
        } else if (mode == "rz-fields") {
            require_argument_count(mode, argc, 5);
            passed = run_rz_fields("rz_case", argv[2], argv[3], std::stod(argv[4]));
        } else if (mode == "cartesian-fields") {
            require_argument_count(mode, argc, 5);
            passed = run_cartesian_fields("cartesian_case", argv[2], argv[3], std::stod(argv[4]));
        } else if (mode == "b549" || mode == "b550" || mode == "b555") {
            require_argument_count(mode, argc, mode == "b549" ? 7 : 8);
            passed = completed_summary(argv[3], "transient");
            const auto summary = read_summary(argv[3]);
            passed = check(summary_number(summary, "accepted_steps") == 20.0 &&
                               summary_number(summary, "rejected_steps") == 0.0,
                         "C3D20T integrated path completes twenty prescribed increments without rejected steps") &&
                     passed;
            const auto active =
                mode == "b549"
                    ? 0
                    : static_cast<std::size_t>(summary_number(summary, "contact.coupled_contact.active_contact_nodes"));
            passed = fuelsim::test::check_hex20_integrated(
                         argv[2], argv[4], argv[5], argv[6], mode == "b549" ? "" : argv[7], mode == "b555", active) &&
                     passed;
        } else if (mode == "hex20-friction-path-33" || mode == "hex20-friction-path-36") {
            require_argument_count(mode, argc, mode == "hex20-friction-path-33" ? 7 : 6);
            passed = completed_summary(argv[3], "transient");
            passed = check(summary_number(read_summary(argv[3]), "accepted_steps") == 70.0,
                         "HEX20 friction path completes seventy prescribed increments") &&
                     passed;
            if (mode == "hex20-friction-path-33")
                passed = fuelsim::test::check_hex20_friction_path_33(argv[2], argv[4], argv[5], argv[6]) && passed;
            else
                passed = fuelsim::test::check_hex20_friction_path_36(argv[2], argv[4], argv[5]) && passed;
        } else if (mode == "hex20-partial-contact") {
            require_argument_count(mode, argc, 5);
            passed = completed_summary(argv[3], "steady");
            passed = check(summary_number(read_summary(argv[3]), "load_steps_completed") == 4.0,
                         "H20.41 completes four load steps") &&
                     passed;
            const auto summary = read_summary(argv[3]);
            passed = check(summary_number(summary, "contact.interface.projected_contact_nodes") == 29.0 &&
                               summary_number(summary, "contact.interface.unprojected_contact_nodes") == 0.0 &&
                               summary_number(summary, "contact.interface.active_contact_nodes") == 11.0,
                         "H20.41 summary reports the physical active-constraint topology") &&
                     passed;
            passed = fuelsim::test::check_hex20_partial_contact(argv[2], argv[4]) && passed;
        } else if (mode == "b40") {
            require_argument_count(mode, argc, 6);
            passed = completed_summary(argv[3], "transient");
            passed = check(summary_number(read_summary(argv[3]), "accepted_steps") == 4.0,
                         "B4.0 completes four prescribed biaxial friction steps") &&
                     passed;
            passed = fuelsim::test::check_hex8_biaxial_friction(argv[2], argv[4], argv[5]) && passed;
        } else if (mode == "b34" || mode == "b34-transient") {
            require_argument_count(mode, argc, 7);
            const bool transient = mode == "b34-transient";
            passed = completed_summary(argv[3], transient ? "transient" : "steady");
            passed = check(summary_number(
                               read_summary(argv[3]), transient ? "accepted_steps" : "load_steps_completed") == 1.0,
                         "B3.4 completes one full loading step") &&
                     passed;
            passed = fuelsim::test::check_hex8_sliding(argv[2], argv[4], argv[5], argv[6]) && passed;
        } else if (mode == "hex20-thermal-contact") {
            require_argument_count(mode, argc, 5);
            passed = run_hex20_thermal_contact(argv[2], argv[3], argv[4]);
        } else if (mode == "hex20-curved-friction") {
            require_argument_count(mode, argc, 6);
            passed = completed_summary(argv[3], "steady");
            passed = check(summary_number(read_summary(argv[3]), "load_steps_completed") == 4.0,
                         "H20.35 completes four load steps") &&
                     passed;
            passed = fuelsim::test::check_hex20_curved_friction(argv[2], argv[4], argv[5]) && passed;
        } else if (mode == "hex20-curved") {
            require_argument_count(mode, argc, 8);
            passed = completed_summary(argv[3], "steady");
            passed = check(summary_number(read_summary(argv[3]), "load_steps_completed") == 4.0,
                         "H20.30 completes four load steps") &&
                     passed;
            passed = fuelsim::test::check_hex20_curved(argv[2], argv[4], argv[5], argv[6], argv[7]) && passed;
        } else if (mode == "b33") {
            require_argument_count(mode, argc, 8);
            passed = completed_summary(argv[3], "steady");
            passed = check(summary_number(read_summary(argv[3]), "load_steps_completed") == 10.0,
                         "B3.3 completes ten load steps") &&
                     passed;
            passed = fuelsim::test::check_hex8_sticking(argv[2], argv[4], argv[5], argv[6], argv[7]) && passed;
        } else if (mode == "b35") {
            require_argument_count(mode, argc, 4);
            passed = fuelsim::test::check_hex8_shared_plate(argv[2], argv[3], "");
        } else if (mode == "b36") {
            require_argument_count(mode, argc, 6);
            passed = completed_summary(argv[3], "transient");
            passed = check(summary_number(read_summary(argv[3]), "accepted_steps") == 40.0,
                         "B3.6 accepts forty time steps") &&
                     passed;
            passed = fuelsim::test::check_hex8_shared_plate(argv[2], argv[4], argv[5]) && passed;
        } else if (mode == "b6") {
            require_argument_count(mode, argc, 5);
            passed = run_b6(argv[2], argv[3], argv[4]);
        } else if (mode == "hex20-fields") {
            require_argument_count(mode, argc, 7);
            const std::string field_selection = argv[6];
            if (field_selection != "all" && field_selection != "temperature-only")
                throw std::invalid_argument("HEX20 field selection must be 'all' or 'temperature-only'");
            passed =
                run_hex20_fields("hex20_case", argv[2], argv[3], argv[4], std::stod(argv[5]), field_selection == "all");
        } else if (mode == "hex20-transient") {
            require_argument_count(mode, argc, 7);
            passed = run_hex20_transient(argv[2], argv[3], argv[4], argv[5], std::stod(argv[6]));
        } else if (mode == "hex20-inelastic") {
            require_argument_count(mode, argc, 7);
            passed = run_hex20_inelastic(argv[2], argv[3], argv[4], argv[5], argv[6]);
        } else if (mode == "hex8-inelastic") {
            require_argument_count(mode, argc, 6);
            passed = run_hex8_inelastic(argv[2], argv[3], argv[4], argv[5]);
        } else if (mode == "rz-steady") {
            require_argument_count(mode, argc, 6);
            passed = run_rz_fields("rz_steady", argv[2], argv[4], std::stod(argv[5]));
            passed = completed_summary(argv[3], "steady") && passed;
        } else if (mode == "frictionless") {
            require_argument_count(mode, argc, 5);
            passed = completed_summary(argv[3], "steady");
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const std::string contact = argv[4];
            const auto& projected = result.nodal("contact_projected_" + contact);
            std::size_t count = 0;
            for (std::size_t node = 0; node < projected.size(); ++node) {
                if (std::isnan(projected[node])) continue;
                ++count;
                passed = check(projected[node] == 1.0 &&
                                   result.nodal("contact_tangential_traction_" + contact)[node] == 0.0 &&
                                   result.nodal("contact_tangential_force_" + contact)[node] == 0.0,
                             "Frictionless baseline has projected contact and exactly zero tangential force") &&
                         passed;
            }
            passed = check(count > 0, "Frictionless baseline contains contact nodes") && passed;
        } else if (mode == "m51") {
            require_argument_count(mode, argc, 7);
            passed = completed_summary(argv[3], "steady");
            passed = run_rz_fields("m51", argv[2], argv[4], 1.0e-2) && passed;
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const auto baseline = fuelsim::test::read_final_exodus_results(argv[6]);
            if (result.nodes != baseline.nodes) throw std::invalid_argument("M5.1 baseline node mapping differs");
            const auto pressure_error =
                compare_contact_pressure(result, "contact_pressure_fuel_cladding", argv[5], 1.0e-12);
            print_relative_metrics("m51_contact_pressure", pressure_error);
            const auto& projected = result.nodal("contact_projected_fuel_cladding");
            double maximum_pressure = 0.0, maximum_excess = 0.0;
            std::size_t count = 0, sliding = 0;
            double difference_squared = 0.0, scale_squared = 0.0;
            for (std::size_t node = 0; node < result.nodes.size(); ++node) {
                const double axial = result.nodal("displacement_z")[node];
                const double difference = axial - baseline.nodal("displacement_z")[node];
                difference_squared += difference * difference;
                scale_squared += axial * axial;
                if (std::isnan(projected[node])) continue;
                ++count;
                const double pressure = result.nodal("contact_pressure_fuel_cladding")[node];
                maximum_pressure = std::max(maximum_pressure, pressure);
                maximum_excess = std::max(maximum_excess,
                    std::abs(result.nodal("contact_tangential_traction_fuel_cladding")[node]) - 0.3 * pressure);
                sliding += result.nodal("contact_sliding_fuel_cladding")[node] == 1.0 ? 1U : 0U;
                passed = check(projected[node] == 1.0 && pressure > 0.0,
                             "M5.1 all contact nodes are active and projected") &&
                         passed;
            }
            const double effect = std::sqrt(difference_squared / scale_squared);
            passed = check(count > 0 && sliding > 0 && maximum_excess <= 1.0e-12 * maximum_pressure &&
                               std::abs(result.global("contact_tangential_force_fuel_cladding")) > 0.0 &&
                               relative_metrics_below(pressure_error, 1.0e-2) && effect > 1.0e-3 &&
                               summary_number(read_summary(argv[3]), "load_steps_completed") == 20.0,
                         "M5.1 retains pressure, friction activation, Coulomb cap and baseline-effect gates") &&
                     passed;
            std::cout << "m51_frictional_axial_field_change=" << effect << " sliding_nodes=" << sliding << '\n';
        } else if (mode == "m52") {
            require_argument_count(mode, argc, 6);
            passed = completed_summary(argv[3], "steady");
            passed = run_rz_fields("m52", argv[2], argv[4], 1.0e-2) && passed;
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const auto pressure = compare_contact_pressure(result, "contact_pressure_pellet_stack", argv[5], 1.0e-12);
            print_relative_metrics("m52_contact_pressure", pressure);
            std::size_t active = 0;
            double maximum_radius = 0.0;
            const auto& projected = result.nodal("contact_projected_pellet_stack");
            for (std::size_t node = 0; node < projected.size(); ++node) {
                if (std::isnan(projected[node])) continue;
                passed = check(projected[node] == 1.0 && result.nodal("contact_pressure_pellet_stack")[node] > 0.0,
                             "M5.2 every contact node remains projected and active") &&
                         passed;
                ++active;
                maximum_radius = std::max(maximum_radius, result.nodal("contact_current_r_pellet_stack")[node]);
            }
            passed = check(result.nodes.size() == 402 && active == 5 && maximum_radius > 4.2e-3 &&
                               maximum_radius < 8.0e-3 && result.global("contact_force_pellet_stack") > 0.0 &&
                               relative_metrics_below(pressure, 1.0e-2) &&
                               summary_number(read_summary(argv[3]), "load_steps_completed") == 20.0 &&
                               summary_number(read_summary(argv[3]), "petsc_workspace_setups") == 1.0,
                         "M5.2 retains the field, radius, contact-force and workspace gates") &&
                     passed;
            std::cout << "m52_maximum_secondary_current_radius=" << maximum_radius << '\n';
        } else if (mode == "m31") {
            require_argument_count(mode, argc, 5);
            passed = completed_summary(argv[3], "transient");
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const auto fields = compare_rz_nodal_results(result, read_rz_nodal_reference(argv[4]));
            print_relative_metrics("m31_temperature", fields.temperature);
            passed =
                check(fields.maximum_coordinate_difference < 1.0e-12 &&
                          relative_metrics_below(fields.temperature, 1.0e-3) &&
                          fuelsim::test::absolute_metrics_below(fields.radial_displacement, 1.0e-12) &&
                          fuelsim::test::absolute_metrics_below(fields.axial_displacement, 1.0e-12) &&
                          result.step_count == 5 && summary_number(read_summary(argv[3]), "accepted_steps") == 4.0 &&
                          summary_number(read_summary(argv[3]), "contributions") == 10.0,
                    "M3.1 retains full fields, zero displacement, four steps and ten contributions") &&
                passed;
            const std::array<double, 5> times = {0.0, 2.5, 5.0, 8.0, 10.0};
            for (std::size_t step = 0; step < times.size(); ++step)
                passed = check(fuelsim::test::read_exodus_results(argv[2], step + 1).time == times[step],
                             "M3.1 lands exactly on each power-table event") &&
                         passed;
        } else if (mode == "pressure-cylinder") {
            require_argument_count(mode, argc, 4);
            passed = completed_summary(argv[3], "steady");
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            for (const auto* component : {"rr", "hoop"}) {
                double stress = 0.0;
                std::size_t count = 0;
                for (std::size_t q = 0; q < 4; ++q)
                    for (const double value :
                        result.element(std::string("stress_") + component + "_q" + std::to_string(q))) {
                        stress += value;
                        ++count;
                    }
                stress /= static_cast<double>(count);
                std::cout << "pressure_average_" << component << "=" << stress << '\n';
                passed = check(count > 0 && std::abs(stress + 1.0e6) / 1.0e6 < 1.0e-10,
                             "Radial pressure produces the analytic solid-cylinder stress") &&
                         passed;
            }
        } else if (mode == "m43") {
            require_argument_count(mode, argc, 7);
            passed = completed_summary(argv[3], "transient");
            const double tolerance = std::string(argv[5]) == "material_oracle" ? 5.0e-6 : 5.0e-3;
            passed = run_rz_fields("m43", argv[2], argv[4], tolerance) && passed;
            const auto fields = compare_rz_nodal_results(
                fuelsim::test::read_final_exodus_results(argv[2]), read_rz_nodal_reference(argv[4]));
            passed = check(fields.temperature.maximum_zero_reference_difference < 1.0e-12 &&
                               fields.radial_displacement.maximum_zero_reference_difference < 1.0e-14 &&
                               fields.axial_displacement.maximum_zero_reference_difference < 1.0e-14 &&
                               summary_number(read_summary(argv[3]), "accepted_steps") == 100.0 &&
                               summary_number(read_summary(argv[3]), "rejected_steps") == 0.0,
                         "M4.3 retains nodal zero-reference gates and one hundred accepted steps without rejection") &&
                     passed;
            passed = fuelsim::test::check_rz_noncoaxial(argv[2], argv[5], argv[6]) && passed;
        } else if (mode == "m22") {
            if (argc != 6 && argc != 7)
                throw std::invalid_argument(
                    "M2.2 checker requires result, summary, nodes, branch and optional history");
            passed = completed_summary(argv[3], "transient");
            passed = run_rz_fields("m22", argv[2], argv[4], 1.0e-3) && passed;
            const std::string branch = argv[5];
            passed = check(summary_number(read_summary(argv[3]), "accepted_steps") ==
                               (branch == "j2_unload_reload" ? 30.0 : 10.0),
                         "M2.2 commits every configured time step") &&
                     passed;
            passed = fuelsim::test::check_rz_inelastic(argv[2], branch, argc == 7 ? argv[6] : "") && passed;
        } else if (mode == "m23" || mode == "m41") {
            require_argument_count(mode, argc, 9);
            const double tolerance = mode == "m23" ? 1.0e-3 : 5.0e-3;
            passed = completed_summary(argv[3], "transient");
            passed = run_rz_fields(mode, argv[2], argv[4], tolerance) && passed;
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const auto pressure = compare_contact_pressure(result, "contact_pressure_fuel_cladding", argv[5], 1.0e-12);
            print_relative_metrics(mode + "_contact_pressure", pressure);
            passed = check(relative_metrics_below(pressure, tolerance) &&
                               summary_number(read_summary(argv[3]), "accepted_steps") == 20.0 &&
                               summary_number(read_summary(argv[3]), "total_cutbacks") == 0.0,
                         "PCMI completes twenty fixed steps without cutbacks and retains the pressure gate") &&
                     passed;
            passed = fuelsim::test::check_rz_pcmi(argv[2], argv[6], argv[7], argv[8], argv[5], tolerance) && passed;
        } else if (mode == "m57") {
            require_argument_count(mode, argc, 9);
            passed = completed_summary(argv[3], "transient");
            passed = run_rz_fields("m57", argv[2], argv[4], 5.0e-3) && passed;
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const auto pressure = compare_contact_pressure(result, "contact_pressure_fuel_cladding", argv[5], 1.0e-12);
            print_relative_metrics("m57_contact_pressure", pressure);
            const auto summary = read_summary(argv[3]);
            passed =
                check(relative_metrics_below(pressure, 5.0e-3) && summary_number(summary, "accepted_steps") == 18.0 &&
                          summary_number(summary, "time_error_rejections") == 0.0 &&
                          summary_number(summary, "maximum_accepted_time_error_estimate") > 0.0 &&
                          summary_number(summary, "maximum_accepted_time_error_estimate") < 1.0 &&
                          summary_number(summary, "minimum_accepted_time_step") <
                              summary_number(summary, "maximum_accepted_time_step") &&
                          std::abs(summary_number(summary, "conservation.interface_heat_imbalance")) < 1.0e-10,
                    "M5.7 retains adaptive-step, pressure and thermal conservation gates") &&
                passed;
            passed = fuelsim::test::check_rz_integrated(argv[2], argv[6], argv[7], argv[8]) && passed;
        } else if (mode == "m54") {
            require_argument_count(mode, argc, 6);
            passed = run_m54(argv[2], argv[3], argv[4], argv[5]);
        } else if (mode == "rz-two-pellet") {
            require_argument_count(mode, argc, 6);
            passed = completed_summary(argv[3], "steady");
            passed = run_rz_fields("m33_two_pellet", argv[2], argv[4], 1.0e-2) && passed;
            const auto result = fuelsim::test::read_final_exodus_results(argv[2]);
            const auto pressure =
                compare_contact_pressure(result, "contact_pressure_pellet_stack", argv[5], 1.0e-12, 0);
            print_relative_metrics("m33_two_pellet_pressure", pressure);
            passed =
                check(relative_metrics_below(pressure, 1.0e-2), "Two-pellet pressure retains three 1 percent gates") &&
                passed;
        } else if (mode == "rz-multi-contact") {
            require_argument_count(mode, argc, 8);
            passed = completed_summary(argv[3], "steady");
            passed =
                check(summary_number(read_summary(argv[3]), "regions") == 3.0, "RZ multi-contact has three regions") &&
                passed;
            passed = run_rz_fields("rz_multi_contact", argv[2], argv[4], 1.0e-3) && passed;
            passed = fuelsim::test::check_rz_multi_contact(argv[2], argv[5], argv[6], std::stod(argv[7])) && passed;
        } else if (mode == "hex8-multi-contact") {
            require_argument_count(mode, argc, 7);
            passed = completed_summary(argv[3], "steady");
            passed = check(summary_number(read_summary(argv[3]), "load_steps_completed") == 1.0 &&
                               summary_number(read_summary(argv[3]), "regions") == 4.0,
                         "B3.7 completes one load step for four regions") &&
                     passed;
            passed = fuelsim::test::check_hex8_multi_contact(argv[2], argv[4], argv[5], argv[6]) && passed;
        } else if (mode == "hex20-nonmatching") {
            require_argument_count(mode, argc, 9);
            passed = completed_summary(argv[3], "steady");
            passed = check(summary_number(read_summary(argv[3]), "load_steps_completed") == 1.0 &&
                               summary_number(read_summary(argv[3]), "petsc_workspace_setups") == 1.0,
                         "H20.24 completes one compression step") &&
                     passed;
            passed =
                fuelsim::test::check_hex20_nonmatching(argv[2], argv[4], argv[5], argv[6], argv[7], argv[8]) && passed;
        } else if (mode == "planar-sts" || mode == "planar-sts-hex8") {
            require_argument_count(mode, argc, 7);
            passed = run_planar_sts(argv[2], argv[3], argv[4], argv[5], argv[6], mode == "planar-sts-hex8");
        } else if (mode == "equivalence") {
            require_argument_count(mode, argc, 5);
            passed = compare_result_files(argv[2], argv[3], std::stod(argv[4]));
        } else if (mode == "b40-restart") {
            require_argument_count(mode, argc, 5);
            passed = completed_summary(argv[3], "transient");
            passed = check(summary_number(read_summary(argv[3]), "accepted_steps") == 2.0,
                         "B4.0 restart completes two remaining time steps") &&
                     passed;
            passed = compare_result_files(argv[2], argv[4], 0.0, false) && passed;
        } else if (mode == "mpi-equivalence" || mode == "restart") {
            require_argument_count(mode, argc, 5);
            passed = completed_summary(argv[3], "transient");
            const auto summary = read_summary(argv[3]);
            passed = check(summary_number(summary, "mpi_ranks") == (mode == "restart" ? 1.0 : 2.0),
                         "Production summary confirms the requested process count") &&
                     passed;
            if (mode == "restart")
                passed = check(summary_number(summary, "accepted_steps") == 5.0,
                             "Restart completes the five remaining time steps") &&
                         passed;
            passed =
                compare_result_files(argv[2], argv[4], mode == "restart" ? 0.0 : 1.0e-10, mode == "mpi-equivalence") &&
                passed;
        } else {
            throw std::invalid_argument("Unknown production result check '" + mode + "'");
        }
        if (passed) std::cout << "[PASS] production-executable result comparison " << mode << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] production result comparison raised: " << error.what() << '\n';
        return 1;
    }
}
