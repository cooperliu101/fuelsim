#include "support/exodus_result_reader.hpp"
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

struct FieldErrorMetrics final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_actual = 0.0;
    double maximum_reference = 0.0;
    double maximum_absolute_difference = 0.0;
    double maximum_pointwise_relative = 0.0;
    double maximum_zero_reference_difference = 0.0;
    std::size_t value_count = 0;
    std::size_t nonzero_reference_count = 0;
    std::size_t zero_reference_count = 0;
    std::size_t maximum_absolute_difference_index = 0;
    std::size_t maximum_pointwise_relative_index = 0;
    double maximum_pointwise_relative_actual = 0.0;
    double maximum_pointwise_relative_reference = 0.0;
    double maximum_absolute_difference_actual = 0.0;
    double maximum_absolute_difference_reference = 0.0;

    void add(double actual, double reference) {
        if (!std::isfinite(actual) || !std::isfinite(reference))
            throw std::invalid_argument("Production and reference field values must be finite");
        const double difference = actual - reference;
        difference_squared += difference * difference;
        reference_squared += reference * reference;
        maximum_actual = std::max(maximum_actual, std::abs(actual));
        maximum_reference = std::max(maximum_reference, std::abs(reference));
        if (std::abs(difference) > maximum_absolute_difference) {
            maximum_absolute_difference = std::abs(difference);
            maximum_absolute_difference_index = value_count;
            maximum_absolute_difference_actual = actual;
            maximum_absolute_difference_reference = reference;
        }
        if (reference != 0.0) {
            const double relative = std::abs(difference) / std::abs(reference);
            if (relative > maximum_pointwise_relative) {
                maximum_pointwise_relative = relative;
                maximum_pointwise_relative_index = value_count;
                maximum_pointwise_relative_actual = actual;
                maximum_pointwise_relative_reference = reference;
            }
            ++nonzero_reference_count;
        } else {
            maximum_zero_reference_difference = std::max(maximum_zero_reference_difference, std::abs(difference));
            ++zero_reference_count;
        }
        ++value_count;
    }

    double relative_l2() const {
        if (reference_squared == 0.0) throw std::domain_error("Relative L2 error is undefined");
        return std::sqrt(difference_squared / reference_squared);
    }

    double relative_absolute_peak() const {
        if (maximum_reference == 0.0) throw std::domain_error("Relative absolute-peak error is undefined");
        return std::abs(maximum_actual - maximum_reference) / maximum_reference;
    }

    double maximum_pointwise_relative_error() const {
        if (nonzero_reference_count == 0) throw std::domain_error("Pointwise relative error is undefined");
        return maximum_pointwise_relative;
    }
};

bool relative_metrics_below(const FieldErrorMetrics& metrics, double tolerance) {
    return metrics.relative_l2() < tolerance && metrics.relative_absolute_peak() < tolerance &&
           metrics.maximum_pointwise_relative_error() < tolerance;
}

void print_relative_metrics(const std::string& name, const FieldErrorMetrics& metrics) {
    std::cout << name << "_relative_l2=" << metrics.relative_l2() << '\n';
    std::cout << name << "_relative_absolute_peak=" << metrics.relative_absolute_peak() << '\n';
    std::cout << name << "_maximum_pointwise_relative=" << metrics.maximum_pointwise_relative_error() << '\n';
    std::cout << name << "_maximum_pointwise_relative_index=" << metrics.maximum_pointwise_relative_index << '\n';
    std::cout << name << "_maximum_pointwise_relative_actual=" << metrics.maximum_pointwise_relative_actual << '\n';
    std::cout << name << "_maximum_pointwise_relative_reference=" << metrics.maximum_pointwise_relative_reference
              << '\n';
    std::cout << name << "_maximum_absolute_difference_index=" << metrics.maximum_absolute_difference_index << '\n';
    std::cout << name << "_maximum_absolute_difference_actual=" << metrics.maximum_absolute_difference_actual << '\n';
    std::cout << name << "_maximum_absolute_difference_reference=" << metrics.maximum_absolute_difference_reference
              << '\n';
    std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n';
    std::cout << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
              << '\n';
}

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

std::vector<double> read_contact_pressure_reference(const std::string& path, std::vector<double>& coordinates) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE pressure reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE pressure reference is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::size_t pressure_column = column_index(header, "contact_pressure", path);
    const std::size_t coordinate_column = column_index(header, "y", path);
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
    const std::string& variable_name, const std::string& reference_path, double coordinate_tolerance) {
    std::vector<double> reference_coordinates;
    const std::vector<double> reference = read_contact_pressure_reference(reference_path, reference_coordinates);
    const std::vector<double>& pressure = results.nodal(variable_name);
    std::vector<std::pair<double, double>> actual;
    for (std::size_t node = 0; node < pressure.size(); ++node)
        if (std::isfinite(pressure[node])) actual.emplace_back(results.nodes[node][1], pressure[node]);
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
    bool passed = check(relative_metrics_below(temperature_error, tolerance),
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
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: fuelsim_production_result_tests "
                     "<m0|m1|m1-unstructured|m21|b3|b6|rz-fields|cartesian-fields|hex20-fields|hex20-transient> "
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
