#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct DisplacementReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> normal;
    double displacement;
};

struct ForceReference final {
    std::size_t id;
    std::array<double, 3> point;
    double force;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete H20.30 CSV row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacements(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.30 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_displacement,normal_x,normal_y,normal_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.30 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 4, path)),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
            number(values, 0, path)});
    }
    return result;
}

std::vector<ForceReference> read_forces(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.30 nodal-force reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_force,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.30 force header in " + path);
    std::vector<ForceReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            number(values, 0, path)});
    }
    return result;
}

std::array<double, 2> read_resultants(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.30 resultant reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "radial_contact_resultant,radial_primary_reaction")
        throw std::invalid_argument("Unexpected H20.30 resultant header in " + path);
    if (!std::getline(input, line))
        throw std::invalid_argument("Missing H20.30 resultant value in " + path);
    const std::vector<std::string> values = split_csv(line);
    return {std::abs(number(values, 0, path)), std::abs(number(values, 1, path))};
}

bool check(bool condition, const std::string& message) {
    if (!condition)
        std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

double coordinate_difference(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::max({std::abs(a[0] - b[0]), std::abs(a[1] - b[1]), std::abs(a[2] - b[2])});
}
} // namespace

namespace fuelsim::test {
bool check_hex20_curved(const std::string& output_path,
    const std::string& name,
    const std::string& displacement_path,
    const std::string& force_path,
    const std::string& resultant_path) {
    const auto output = read_final_exodus_results(output_path);
    const auto displacement = read_displacements(displacement_path);
    const auto force = read_forces(force_path);
    if (displacement.size() != output.nodes.size() || force.empty())
        throw std::runtime_error("H20.30 reference node counts differ");
    const std::array<std::string, 3> displacements = {"displacement_x", "displacement_y", "displacement_z"};
    const std::array<std::string, 3> forces = {"contact_normal_force_x_interface",
        "contact_normal_force_y_interface",
        "contact_normal_force_z_interface"};
    std::vector<bool> present(output.nodes.size(), false), force_present(output.nodes.size(), false);
    FieldErrorMetrics displacement_error, force_error;
    double coordinate_error = 0.0, force_coordinate_error = 0.0;
    for (const auto& row : displacement) {
        if (row.id >= present.size() || present[row.id])
            throw std::runtime_error("H20.30 repeated or invalid node ID");
        present[row.id] = true;
        coordinate_error = std::max(coordinate_error, coordinate_difference(output.nodes[row.id], row.point));
        double actual = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            actual += row.normal[component] * output.nodal(displacements[component]).at(row.id);
        displacement_error.add(actual, row.displacement);
    }
    for (const auto& row : force) {
        if (row.id >= force_present.size() || force_present[row.id])
            throw std::runtime_error("H20.30 repeated or invalid force node ID");
        force_present[row.id] = true;
        force_coordinate_error =
            std::max(force_coordinate_error, coordinate_difference(output.nodes[row.id], row.point));
        const double radius = std::hypot(row.point[0], row.point[1]);
        if (!(radius > 0.0))
            throw std::runtime_error("H20.30 radial force has zero radius");
        const double actual =
            -(row.point[0] * output.nodal(forces[0]).at(row.id) + row.point[1] * output.nodal(forces[1]).at(row.id))
            / radius;
        force_error.add(actual, row.force);
    }
    const auto& projected = output.nodal("contact_projected_interface");
    double resultant = 0.0;
    std::size_t count = 0;
    bool passed = true;
    for (std::size_t node = 0; node < output.nodes.size(); ++node) {
        if (std::isnan(projected[node]))
            continue;
        ++count;
        passed = check(projected[node] == 1.0 && output.nodal("contact_pressure_interface")[node] > 0.0
                           && force_present[node],
                     "H20.30 every contact constraint is active, projected and compared")
                 && passed;
        const auto& point = output.nodes[node];
        resultant -= (point[0] * output.nodal(forces[0])[node] + point[1] * output.nodal(forces[1])[node])
                     / std::hypot(point[0], point[1]);
    }
    const auto reference_resultants = read_resultants(resultant_path);
    if (!(reference_resultants[0] > 0.0))
        throw std::runtime_error("H20.30 reference resultant is not positive");
    const double resultant_error = std::abs(resultant - reference_resultants[0]) / reference_resultants[0];
    print_relative_metrics("h20_30_" + name + "_radial_displacement", displacement_error);
    print_relative_metrics("h20_30_" + name + "_nodal_radial_force", force_error);
    std::cout << "h20_30_" << name << "_radial_contact_resultant=" << resultant << '\n'
              << "scalar_contact_force_sum=" << output.global("contact_force_interface") << '\n'
              << "abaqus_radial_contact_resultant=" << reference_resultants[0] << '\n'
              << "abaqus_primary_radial_reaction=" << reference_resultants[1] << '\n'
              << "radial_resultant_relative_error=" << resultant_error << '\n';
    return check(count == force.size() && count > 0, "H20.30 all secondary nodal radial forces are compared")
           && check(coordinate_error < 1.0e-12 && force_coordinate_error < 1.0e-12,
               "H20.30 reference coordinates match the tracked mesh")
           && check(relative_metrics_below(displacement_error, 1.0e-2),
               "H20.30 radial displacement errors pass 1 percent")
           && check(relative_metrics_below(force_error, 1.0e-2), "H20.30 radial force errors pass 1 percent")
           && check(resultant_error < 1.0e-2, "H20.30 radial resultant error passes 1 percent") && passed;
}
} // namespace fuelsim::test
