#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct DisplacementReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> displacement{};
};

struct PressureReference final {
    std::size_t id;
    std::array<double, 3> point;
    double pressure, normal_force;
};

struct MortarReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> displacement{};
    double pressure;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        values.push_back(value);
    return values;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete H20.24 reference row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacement(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.24 Abaqus displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.24 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 3, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path)}});
    }
    return result;
}

std::vector<PressureReference> read_pressure(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.24 Abaqus pressure reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "pressure,normal_force_x,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.24 pressure and force header in " + path);
    std::vector<PressureReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 2, path)),
            {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
            number(values, 0, path),
            number(values, 1, path)});
    }
    return result;
}

double read_normal_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.24 Abaqus reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected H20.24 reaction header in " + path);
    std::vector<std::string> final_values;
    while (std::getline(input, line))
        if (!line.empty())
            final_values = split_csv(line);
    if (final_values.empty() || std::abs(number(final_values, 0, path) - 1.0) > 1.0e-12)
        throw std::invalid_argument("H20.24 Abaqus reaction does not end at unit load");
    return std::abs(number(final_values, 1, path));
}

std::vector<MortarReference> read_mortar_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.24 MOOSE mortar nodal reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "contact_pressure,disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.24 MOOSE mortar nodal header in " + path);
    std::vector<MortarReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 4, path)),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
            number(values, 0, path)});
    }
    return result;
}

double read_mortar_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.24 MOOSE mortar reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time,reaction_x")
        throw std::invalid_argument("Unexpected H20.24 MOOSE mortar reaction header in " + path);
    std::vector<std::string> final_values;
    while (std::getline(input, line))
        if (!line.empty())
            final_values = split_csv(line);
    if (final_values.empty() || std::abs(number(final_values, 0, path) - 1.0) > 1.0e-12)
        throw std::invalid_argument("H20.24 MOOSE mortar reaction does not end at unit load");
    return std::abs(number(final_values, 1, path));
}

double coordinate_difference(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return std::max({std::abs(first[0] - second[0]), std::abs(first[1] - second[1]), std::abs(first[2] - second[2])});
}
} // namespace

bool fuelsim::test::check_hex20_nonmatching(const std::string& output,
    const std::string& displacement_path,
    const std::string& reaction_path,
    const std::string& pressure_path,
    const std::string& mortar_path,
    const std::string& mortar_reaction_path) {
    const auto result = read_final_exodus_results(output);
    bool passed = true;
    const auto displacement = read_displacement(displacement_path);
    const auto mortar = read_mortar_nodes(mortar_path);
    fuelsim::test::FieldErrorMetrics displacement_x;
    fuelsim::test::FieldErrorMetrics displacement_x_mortar;
    fuelsim::test::FieldErrorMetrics mortar_displacement_x_abaqus;
    std::vector<bool> present(result.nodes.size(), false);
    double maximum_coordinate_difference = 0.0;
    double maximum_mortar_coordinate_difference = 0.0;
    for (std::size_t source = 0; source < result.nodes.size(); ++source) {
        const auto found = std::find_if(displacement.begin(),
            displacement.end(),
            [&](const DisplacementReference& value) { return value.id == source; });
        if (found == displacement.end() || present[source])
            throw std::invalid_argument("H20.24 displacement source-node mapping is incomplete or repeated");
        const auto found_mortar = std::find_if(mortar.begin(), mortar.end(), [&](const MortarReference& value) {
            return value.id == source;
        });
        if (found_mortar == mortar.end())
            throw std::invalid_argument("H20.24 MOOSE mortar source-node mapping is incomplete");
        present[source] = true;
        maximum_coordinate_difference =
            std::max(maximum_coordinate_difference, coordinate_difference(result.nodes[source], found->point));
        maximum_mortar_coordinate_difference = std::max(maximum_mortar_coordinate_difference,
            coordinate_difference(result.nodes[source], found_mortar->point));
        const double actual = result.nodal("displacement_x")[source];
        displacement_x.add(actual, found->displacement[0]);
        displacement_x_mortar.add(actual, found_mortar->displacement[0]);
        mortar_displacement_x_abaqus.add(found_mortar->displacement[0], found->displacement[0]);
    }

    const auto pressure = read_pressure(pressure_path);
    std::vector<std::size_t> secondary_sources;
    const auto& projected = result.nodal("contact_projected_interface");
    const auto& actual_pressure = result.nodal("contact_pressure_interface");
    const auto& actual_force = result.nodal("contact_normal_force_interface");
    std::size_t active_nodes = 0, unprojected_nodes = 0;
    for (std::size_t node = 0; node < projected.size(); ++node) {
        if (std::isnan(projected[node]))
            continue;
        secondary_sources.push_back(node);
        if (projected[node] == 1.0 && actual_pressure[node] > 0.0)
            ++active_nodes;
        if (projected[node] != 1.0)
            ++unprojected_nodes;
    }
    fuelsim::test::FieldErrorMetrics contact_pressure;
    fuelsim::test::FieldErrorMetrics contact_force;
    fuelsim::test::FieldErrorMetrics contact_pressure_mortar;
    fuelsim::test::FieldErrorMetrics mortar_contact_pressure_abaqus;
    double maximum_pressure_coordinate_difference = 0.0;
    for (std::size_t node = 0; node < secondary_sources.size(); ++node) {
        const auto found = std::find_if(pressure.begin(), pressure.end(), [&](const PressureReference& value) {
            return value.id == secondary_sources.at(node);
        });
        if (found == pressure.end())
            throw std::invalid_argument("H20.24 pressure source-node mapping is incomplete");
        const auto found_mortar = std::find_if(mortar.begin(), mortar.end(), [&](const MortarReference& value) {
            return value.id == secondary_sources.at(node);
        });
        if (found_mortar == mortar.end())
            throw std::invalid_argument("H20.24 MOOSE mortar pressure source-node mapping is incomplete");
        maximum_pressure_coordinate_difference = std::max(maximum_pressure_coordinate_difference,
            coordinate_difference(result.nodes[found->id], found->point));
        contact_pressure.add(actual_pressure[secondary_sources[node]], found->pressure);
        contact_force.add(actual_force[secondary_sources[node]], found->normal_force);
        contact_pressure_mortar.add(actual_pressure[secondary_sources[node]], found_mortar->pressure);
        mortar_contact_pressure_abaqus.add(found_mortar->pressure, found->pressure);
    }
    const double total_force = result.global("contact_force_interface");
    const double reference_force = read_normal_reaction(reaction_path);
    const double mortar_force = read_mortar_reaction(mortar_reaction_path);
    const double force_error = std::abs(total_force - reference_force) / reference_force;
    const double mortar_force_error = std::abs(total_force - mortar_force) / mortar_force;
    const double mortar_abaqus_force_error = std::abs(mortar_force - reference_force) / reference_force;
    fuelsim::test::print_relative_metrics("h20_24_displacement_x", displacement_x);
    fuelsim::test::print_relative_metrics("h20_24_contact_pressure", contact_pressure);
    fuelsim::test::print_relative_metrics("h20_24_contact_force", contact_force);
    fuelsim::test::print_relative_metrics("h20_24_displacement_x_mortar", displacement_x_mortar);
    fuelsim::test::print_relative_metrics("h20_24_contact_pressure_mortar", contact_pressure_mortar);
    fuelsim::test::print_relative_metrics("h20_24_mortar_displacement_x_abaqus", mortar_displacement_x_abaqus);
    fuelsim::test::print_relative_metrics("h20_24_mortar_contact_pressure_abaqus", mortar_contact_pressure_abaqus);
    std::cout << "h20_24_normal_resultant=" << total_force << '\n'
              << "h20_24_abaqus_normal_resultant=" << reference_force << '\n'
              << "h20_24_mortar_normal_resultant=" << mortar_force << '\n'
              << "h20_24_normal_resultant_relative_error=" << force_error << '\n'
              << "h20_24_mortar_normal_resultant_relative_error=" << mortar_force_error << '\n'
              << "h20_24_mortar_abaqus_normal_resultant_relative_error=" << mortar_abaqus_force_error << '\n';

    constexpr double abaqus_displacement_tolerance = 1.0e-2;
    constexpr double field_tolerance = 6.0e-2;
    constexpr double abaqus_contact_force_tolerance = 1.0e-2;
    constexpr double mortar_displacement_tolerance = 1.0e-2;
    constexpr double mortar_pressure_tolerance = 9.0e-2;
    constexpr double mortar_abaqus_pressure_tolerance = 7.0e-2;
    constexpr double resultant_tolerance = 5.0e-3;
    passed =
        check(displacement.size() == result.nodes.size() && displacement_x.value_count == result.nodes.size(),
            "H20.24 compares all eighty-eight normal displacements")
        && check(pressure.size() == 13 && contact_pressure.value_count == 13,
            "H20.24 compares all thirteen secondary-face nodal pressures")
        && check(mortar.size() == result.nodes.size() && displacement_x_mortar.value_count == result.nodes.size(),
            "H20.24 compares all eighty-eight normal displacements against MOOSE mortar")
        && check(maximum_coordinate_difference < 1.0e-6 && maximum_pressure_coordinate_difference < 1.0e-6,
            "H20.24 Abaqus fields use the tracked nonmatching Exodus coordinates")
        && check(maximum_mortar_coordinate_difference < 1.0e-12,
            "H20.24 MOOSE mortar fields use the tracked nonmatching Exodus coordinates")
        && check(active_nodes == 13 && unprojected_nodes == 0,
            "H20.24 keeps all thirteen secondary contact nodes active and projected")
        && check(fuelsim::test::relative_metrics_below(displacement_x, abaqus_displacement_tolerance),
            "H20.24 normal-displacement three Abaqus errors are below 1 percent")
        && check(fuelsim::test::relative_metrics_below(contact_pressure, field_tolerance),
            "H20.24 contact-pressure three Abaqus errors are below 6 percent")
        && check(fuelsim::test::relative_metrics_below(contact_force, abaqus_contact_force_tolerance),
            "H20.24 nodal normal-contact-force three Abaqus errors are below 1 percent")
        && check(fuelsim::test::relative_metrics_below(displacement_x_mortar, mortar_displacement_tolerance),
            "H20.24 normal-displacement three MOOSE mortar errors are below 1 percent")
        && check(fuelsim::test::relative_metrics_below(contact_pressure_mortar, mortar_pressure_tolerance),
            "H20.24 Abaqus-style contact-pressure three MOOSE mortar diagnostic errors are below 9 percent")
        && check(fuelsim::test::relative_metrics_below(mortar_displacement_x_abaqus, mortar_displacement_tolerance),
            "H20.24 MOOSE mortar normal-displacement three Abaqus errors are below 1 percent")
        && check(
            fuelsim::test::relative_metrics_below(mortar_contact_pressure_abaqus, mortar_abaqus_pressure_tolerance),
            "H20.24 MOOSE mortar contact-pressure three Abaqus errors are below 7 percent")
        && check(force_error < resultant_tolerance, "H20.24 normal resultant agrees with Abaqus below 0.5 percent")
        && check(mortar_force_error < resultant_tolerance,
            "H20.24 normal resultant agrees with MOOSE mortar below 0.5 percent")
        && check(mortar_abaqus_force_error < resultant_tolerance,
            "H20.24 MOOSE mortar normal resultant agrees with Abaqus below 0.5 percent")
        && passed;
    return passed;
}
