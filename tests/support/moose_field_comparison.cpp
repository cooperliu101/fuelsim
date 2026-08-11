#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuelsim::test {
namespace {

struct CsvColumns final {
    std::size_t temperature;
    std::size_t radial_displacement;
    std::size_t axial_displacement;
    std::size_t id;
    std::size_t radius;
    std::size_t axial_coordinate;
};

struct ActualNodalField final {
    double radius = 0.0;
    double axial_coordinate = 0.0;
    double temperature = 0.0;
    double radial_displacement = 0.0;
    double axial_displacement = 0.0;
    bool present = false;
};

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> result;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        result.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos)
            return result;
        begin = separator + 1;
    }
}

std::size_t column_index(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end())
        throw std::invalid_argument("MOOSE CSV is missing column '" + name + "'");
    return static_cast<std::size_t>(found - header.begin());
}

double csv_double(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    if (column >= fields.size())
        throw std::invalid_argument("MOOSE CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("MOOSE CSV contains an invalid number: " + path);
    return value;
}

std::size_t csv_id(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    const double value = csv_double(fields, column, path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("MOOSE CSV contains an invalid node ID: " + path);
    return static_cast<std::size_t>(value);
}

NodalFieldComparison compare_values(const std::vector<ActualNodalField>& actual,
                                    const std::vector<NodalFieldReference>& reference) {
    if (actual.size() != reference.size())
        throw std::invalid_argument("MOOSE and fuelsim full-field node counts differ");
    NodalFieldComparison result;
    result.node_count = actual.size();
    for (std::size_t node = 0; node < actual.size(); ++node) {
        if (!actual[node].present)
            throw std::invalid_argument("fuelsim full-field mapping does not cover every source node");
        result.maximum_coordinate_difference =
            std::max({result.maximum_coordinate_difference, std::abs(actual[node].radius - reference[node].radius),
                      std::abs(actual[node].axial_coordinate - reference[node].axial_coordinate)});
        result.temperature.add(actual[node].temperature, reference[node].temperature);
        result.radial_displacement.add(actual[node].radial_displacement, reference[node].radial_displacement);
        result.axial_displacement.add(actual[node].axial_displacement, reference[node].axial_displacement);
    }
    return result;
}

std::vector<ActualNodalField> steady_values(const SteadyProblem& problem, const std::vector<double>& state,
                                            std::size_t source_node_count) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument("Steady full-field state size mismatch");
    std::vector<ActualNodalField> result(source_node_count);
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const RegionMesh& mesh = problem.region_mesh(region);
        const std::size_t offset = problem.region_node_offset(region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const std::size_t source = mesh.source_node_ids().at(local);
            if (source >= result.size() || result[source].present)
                throw std::invalid_argument("Invalid or duplicate fuelsim source-node mapping");
            const std::size_t global = offset + local;
            result[source] = {mesh.nodes()[local].r,
                              mesh.nodes()[local].z,
                              state[problem.dof_map().temperature(global)],
                              state[problem.dof_map().radial_displacement(global)],
                              state[problem.dof_map().axial_displacement(global)],
                              true};
        }
    }
    return result;
}

std::vector<ActualNodalField> transient_values(const TransientProblem& problem, const std::vector<double>& state,
                                               std::size_t source_node_count) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument("Transient full-field state size mismatch");
    std::vector<ActualNodalField> result(source_node_count);
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const RegionMesh& mesh = problem.region_mesh(region);
        const std::size_t offset = problem.region_node_offset(region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const std::size_t source = mesh.source_node_ids().at(local);
            if (source >= result.size() || result[source].present)
                throw std::invalid_argument("Invalid or duplicate fuelsim source-node mapping");
            const std::size_t global = offset + local;
            result[source] = {mesh.nodes()[local].r,
                              mesh.nodes()[local].z,
                              state[problem.dof_map().temperature(global)],
                              state[problem.dof_map().radial_displacement(global)],
                              state[problem.dof_map().axial_displacement(global)],
                              true};
        }
    }
    return result;
}

} // namespace

void FieldErrorMetrics::add(double actual, double reference) {
    if (!std::isfinite(actual) || !std::isfinite(reference))
        throw std::invalid_argument("Full-field values must be finite");
    const double difference = actual - reference;
    difference_squared += difference * difference;
    reference_squared += reference * reference;
    maximum_actual = std::max(maximum_actual, std::abs(actual));
    maximum_reference = std::max(maximum_reference, std::abs(reference));
    if (std::abs(difference) > maximum_absolute_difference) {
        maximum_absolute_difference = std::abs(difference);
        maximum_absolute_difference_index = value_count;
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

bool FieldErrorMetrics::has_relative_norm() const noexcept {
    return reference_squared > 0.0 && maximum_reference > 0.0 && nonzero_reference_count > 0;
}

double FieldErrorMetrics::relative_l2() const {
    if (!has_relative_norm())
        throw std::domain_error("Relative full-field norm is undefined");
    return std::sqrt(difference_squared / reference_squared);
}

double FieldErrorMetrics::relative_absolute_peak() const {
    if (!has_relative_norm())
        throw std::domain_error("Relative full-field peak is undefined");
    return std::abs(maximum_actual - maximum_reference) / maximum_reference;
}

double FieldErrorMetrics::maximum_pointwise_relative_error() const {
    if (!has_relative_norm())
        throw std::domain_error("Pointwise relative full-field error is undefined");
    return maximum_pointwise_relative;
}

double FieldErrorMetrics::absolute_l2() const noexcept {
    return std::sqrt(difference_squared);
}

double FieldErrorMetrics::absolute_peak() const noexcept {
    return std::abs(maximum_actual - maximum_reference);
}

std::vector<NodalFieldReference> read_moose_nodal_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE nodal reference is empty: " + path);
    const std::vector<std::string> header = split_csv_line(line);
    const CsvColumns columns = {column_index(header, "T"),      column_index(header, "disp_x"),
                                column_index(header, "disp_y"), column_index(header, "id"),
                                column_index(header, "x"),      column_index(header, "y")};
    std::vector<NodalFieldReference> result;
    std::vector<bool> present;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv_line(line);
        const std::size_t id = csv_id(fields, columns.id, path);
        if (id >= result.size()) {
            result.resize(id + 1);
            present.resize(id + 1, false);
        }
        if (present[id])
            throw std::invalid_argument("MOOSE nodal reference contains a duplicate node ID: " + path);
        result[id] = {csv_double(fields, columns.radius, path), csv_double(fields, columns.axial_coordinate, path),
                      csv_double(fields, columns.temperature, path),
                      csv_double(fields, columns.radial_displacement, path),
                      csv_double(fields, columns.axial_displacement, path)};
        present[id] = true;
    }
    if (result.empty() || std::any_of(present.begin(), present.end(), [](bool value) { return !value; }))
        throw std::invalid_argument("MOOSE nodal reference IDs must be contiguous: " + path);
    return result;
}

NodalFieldComparison compare_moose_nodal_fields(const SteadyProblem& problem, const std::vector<double>& state,
                                                const std::vector<NodalFieldReference>& reference) {
    return compare_values(steady_values(problem, state, reference.size()), reference);
}

NodalFieldComparison compare_moose_nodal_fields(const TransientProblem& problem, const std::vector<double>& state,
                                                const std::vector<NodalFieldReference>& reference) {
    return compare_values(transient_values(problem, state, reference.size()), reference);
}

std::vector<double> read_moose_contact_pressure_reference(const std::string& path, std::vector<double>& coordinates) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE pressure reference: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE pressure reference is empty: " + path);
    const std::vector<std::string> header = split_csv_line(line);
    const std::size_t pressure = column_index(header, "contact_pressure");
    const std::size_t coordinate = column_index(header, "y");
    std::vector<std::pair<double, double>> values;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv_line(line);
        values.emplace_back(csv_double(fields, coordinate, path), csv_double(fields, pressure, path));
    }
    std::sort(values.begin(), values.end());
    coordinates.clear();
    std::vector<double> result;
    coordinates.reserve(values.size());
    result.reserve(values.size());
    for (const auto& value : values) {
        coordinates.push_back(value.first);
        result.push_back(value.second);
    }
    if (result.empty())
        throw std::invalid_argument("MOOSE pressure reference has no values: " + path);
    return result;
}

FieldErrorMetrics compare_moose_contact_pressure(const std::vector<ContactNodeSummary>& actual,
                                                 const std::vector<double>& reference,
                                                 const std::vector<double>& reference_coordinates,
                                                 double coordinate_tolerance) {
    if (actual.size() != reference.size() || reference.size() != reference_coordinates.size())
        throw std::invalid_argument("MOOSE and fuelsim contact-pressure counts differ");
    FieldErrorMetrics result;
    for (std::size_t node = 0; node < actual.size(); ++node) {
        if (std::abs(actual[node].z - reference_coordinates[node]) > coordinate_tolerance)
            throw std::invalid_argument("MOOSE and fuelsim contact coordinates differ");
        result.add(actual[node].pressure, reference[node]);
    }
    return result;
}

bool relative_metrics_below(const FieldErrorMetrics& metrics, double tolerance) {
    return metrics.relative_l2() < tolerance && metrics.relative_absolute_peak() < tolerance &&
           metrics.maximum_pointwise_relative_error() < tolerance;
}

bool relative_metrics_below_with_pointwise_tolerance(const FieldErrorMetrics& metrics, double aggregate_tolerance,
                                                     double pointwise_tolerance) {
    return metrics.relative_l2() < aggregate_tolerance && metrics.relative_absolute_peak() < aggregate_tolerance &&
           metrics.maximum_pointwise_relative_error() < pointwise_tolerance;
}

bool absolute_metrics_below(const FieldErrorMetrics& metrics, double tolerance) {
    return metrics.absolute_l2() < tolerance && metrics.absolute_peak() < tolerance &&
           metrics.maximum_absolute_difference < tolerance;
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
    std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n';
    std::cout << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
              << '\n';
}

void print_absolute_metrics(const std::string& name, const FieldErrorMetrics& metrics) {
    std::cout << name << "_absolute_l2=" << metrics.absolute_l2() << '\n';
    std::cout << name << "_absolute_peak=" << metrics.absolute_peak() << '\n';
    std::cout << name << "_maximum_pointwise_absolute=" << metrics.maximum_absolute_difference << '\n';
}

} // namespace fuelsim::test
