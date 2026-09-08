#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
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
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max())
        || std::floor(value) != value)
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
        result.maximum_coordinate_difference = std::max({result.maximum_coordinate_difference,
            std::abs(actual[node].radius - reference[node].radius),
            std::abs(actual[node].axial_coordinate - reference[node].axial_coordinate)});
        result.temperature.add(actual[node].temperature, reference[node].temperature);
        result.radial_displacement.add(actual[node].radial_displacement, reference[node].radial_displacement);
        result.axial_displacement.add(actual[node].axial_displacement, reference[node].axial_displacement);
    }
    return result;
}

std::vector<ActualNodalField>
steady_values(const SteadyProblem& problem, const std::vector<double>& state, std::size_t source_node_count) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument("Steady full-field state size mismatch");
    std::vector<ActualNodalField> result(source_node_count);
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const std::size_t source = mesh.source_node_ids().at(local);
            if (source >= result.size())
                throw std::invalid_argument("Invalid fuelsim source-node mapping");
            const std::size_t global = fuelsim::rz::ProblemAccess::dof_map(problem).global_node(region, local);
            if (result[source].present) {
                const bool identical =
                    result[source].radius == mesh.nodes()[local].r
                    && result[source].axial_coordinate == mesh.nodes()[local].z
                    && result[source].temperature
                           == state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature,
                               global)]
                    && result[source].radial_displacement
                           == state[fuelsim::rz::ProblemAccess::dof_map(problem)
                                   .dof(fuelsim::Field::radial_displacement, global)]
                    && result[source].axial_displacement
                           == state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
                               global)];
                if (!identical)
                    throw std::invalid_argument("Shared source node does not map to one consistent RZ state");
                continue;
            }
            result[source] = {mesh.nodes()[local].r,
                mesh.nodes()[local].z,
                state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, global)],
                state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, global)],
                state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, global)],
                true};
        }
    }
    return result;
}

std::vector<ActualNodalField>
transient_values(const TransientProblem& problem, const std::vector<double>& state, std::size_t source_node_count) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument("Transient full-field state size mismatch");
    std::vector<ActualNodalField> result(source_node_count);
    for (std::size_t region = 0; region < fuelsim::rz::ProblemAccess::region_count(problem); ++region) {
        const RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const std::size_t source = mesh.source_node_ids().at(local);
            if (source >= result.size())
                throw std::invalid_argument("Invalid fuelsim source-node mapping");
            const std::size_t global = fuelsim::rz::ProblemAccess::dof_map(problem).global_node(region, local);
            if (result[source].present) {
                const bool identical =
                    result[source].radius == mesh.nodes()[local].r
                    && result[source].axial_coordinate == mesh.nodes()[local].z
                    && result[source].temperature
                           == state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature,
                               global)]
                    && result[source].radial_displacement
                           == state[fuelsim::rz::ProblemAccess::dof_map(problem)
                                   .dof(fuelsim::Field::radial_displacement, global)]
                    && result[source].axial_displacement
                           == state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
                               global)];
                if (!identical)
                    throw std::invalid_argument("Shared source node does not map to one consistent RZ state");
                continue;
            }
            result[source] = {mesh.nodes()[local].r,
                mesh.nodes()[local].z,
                state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, global)],
                state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, global)],
                state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, global)],
                true};
        }
    }
    return result;
}
} // namespace

std::vector<NodalFieldReference> read_moose_nodal_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE nodal reference is empty: " + path);
    const std::vector<std::string> header = split_csv_line(line);
    const CsvColumns columns = {column_index(header, "T"),
        column_index(header, "disp_x"),
        column_index(header, "disp_y"),
        column_index(header, "id"),
        column_index(header, "x"),
        column_index(header, "y")};
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
        const NodalFieldReference candidate = {csv_double(fields, columns.radius, path),
            csv_double(fields, columns.axial_coordinate, path),
            csv_double(fields, columns.temperature, path),
            csv_double(fields, columns.radial_displacement, path),
            csv_double(fields, columns.axial_displacement, path)};
        if (present[id]) {
            const NodalFieldReference& existing = result[id];
            if (std::abs(existing.radius - candidate.radius) > 1.0e-12
                || std::abs(existing.axial_coordinate - candidate.axial_coordinate) > 1.0e-12
                || std::abs(existing.temperature - candidate.temperature) > 1.0e-12
                || std::abs(existing.radial_displacement - candidate.radial_displacement) > 1.0e-12
                || std::abs(existing.axial_displacement - candidate.axial_displacement) > 1.0e-12)
                throw std::invalid_argument("MOOSE nodal reference has inconsistent duplicate node values: " + path);
            continue;
        }
        result[id] = candidate;
        present[id] = true;
    }
    if (result.empty() || std::any_of(present.begin(), present.end(), [](bool value) { return !value; }))
        throw std::invalid_argument("MOOSE nodal reference IDs must be contiguous: " + path);
    return result;
}

NodalFieldComparison compare_moose_nodal_fields(const SteadyProblem& problem,
    const std::vector<double>& state,
    const std::vector<NodalFieldReference>& reference) {
    return compare_values(steady_values(problem, state, reference.size()), reference);
}

NodalFieldComparison compare_moose_nodal_fields(const TransientProblem& problem,
    const std::vector<double>& state,
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

} // namespace fuelsim::test
