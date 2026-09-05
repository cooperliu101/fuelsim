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
struct NodeReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 4> fields{};
};

struct ContactReference final {
    std::size_t id;
    std::array<double, 3> point;
    double pressure;
};

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
    if (index >= values.size()) throw std::invalid_argument("Incomplete three-dimensional contact row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_nodes(const std::string& thermal_path, const std::string& mechanical_path) {
    std::ifstream thermal(thermal_path);
    if (!thermal) throw std::runtime_error("Could not read three-dimensional thermal-contact nodes: " + thermal_path);
    std::string line;
    std::getline(thermal, line);
    if (line != "T,id,x,y,z") throw std::invalid_argument("Unexpected thermal-contact header in " + thermal_path);
    std::vector<NodeReference> result;
    while (std::getline(thermal, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, thermal_path)),
            {number(values, 2, thermal_path), number(values, 3, thermal_path), number(values, 4, thermal_path)},
            {number(values, 0, thermal_path), 0.0, 0.0, 0.0}});
    }
    std::ifstream mechanical(mechanical_path);
    if (!mechanical)
        throw std::runtime_error("Could not read three-dimensional friction-contact nodes: " + mechanical_path);
    std::getline(mechanical, line);
    if (line != "disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected friction-contact header in " + mechanical_path);
    std::size_t row = 0;
    while (std::getline(mechanical, line)) {
        const auto values = split_csv(line);
        if (row >= result.size() || result[row].id != static_cast<std::size_t>(number(values, 3, mechanical_path)))
            throw std::invalid_argument("Thermal and mechanical MOOSE node order differs");
        const std::array<double, 3> point = {
            number(values, 4, mechanical_path), number(values, 5, mechanical_path), number(values, 6, mechanical_path)};
        if (std::max({std::abs(point[0] - result[row].point[0]), std::abs(point[1] - result[row].point[1]),
                std::abs(point[2] - result[row].point[2])}) > 1.0e-12)
            throw std::invalid_argument("Thermal and mechanical MOOSE node coordinates differ");
        result[row].fields[1] = number(values, 0, mechanical_path);
        result[row].fields[2] = number(values, 1, mechanical_path);
        result[row].fields[3] = number(values, 2, mechanical_path);
        ++row;
    }
    if (row != result.size()) throw std::invalid_argument("Thermal and mechanical MOOSE node counts differ");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional contact pressure: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "contact_pressure,id,nodal_area,penetration,x,y,z")
        throw std::invalid_argument("Unexpected three-dimensional contact-pressure header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}, number(values, 0, path)});
    }
    return result;
}

} // namespace

namespace fuelsim::test {
bool check_hex8_sticking(const std::string& output_path, const std::string& thermal_path,
    const std::string& mechanical_path, const std::string& contact_path, const std::string& baseline_path) {
    const auto output = read_final_exodus_results(output_path);
    const auto baseline = read_final_exodus_results(baseline_path);
    const auto nodes = read_nodes(thermal_path, mechanical_path);
    if (nodes.size() != output.nodes.size() || baseline.nodes != output.nodes)
        throw std::runtime_error("B3.3 node mapping differs");
    std::array<FieldErrorMetrics, 4> fields;
    const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t n = 0; n < nodes.size(); ++n) {
        if (nodes[n].id != n) throw std::runtime_error("B3.3 node IDs are not contiguous");
        for (std::size_t c = 0; c < 3; ++c)
            if (std::abs(nodes[n].point[c] - output.nodes[n][c]) >= 1.0e-12)
                throw std::runtime_error("B3.3 nodal reference coordinates differ");
        for (std::size_t f = 0; f < 4; ++f) fields[f].add(output.nodal(names[f]).at(n), nodes[n].fields[f]);
    }
    bool passed = true;
    for (std::size_t f = 0; f < 4; ++f) {
        print_relative_metrics("b33_" + names[f], fields[f]);
        passed = check(relative_metrics_below(fields[f], 5.0e-3), "B3.3 nodal errors pass 0.5 percent") && passed;
    }
    const auto contact = read_contact(contact_path);
    const auto& projected = output.nodal("contact_projected_interface");
    const auto& pressures = output.nodal("contact_pressure_interface");
    std::vector<bool> seen(output.nodes.size(), false);
    FieldErrorMetrics pressure_error;
    for (const auto& row : contact) {
        if (row.id >= seen.size() || seen[row.id]) throw std::runtime_error("B3.3 contact reference ID is invalid");
        seen[row.id] = true;
        for (std::size_t c = 0; c < 3; ++c)
            if (std::abs(output.nodes[row.id][c] - row.point[c]) >= 1.0e-12)
                throw std::runtime_error("B3.3 contact coordinates differ");
        pressure_error.add(pressures.at(row.id), row.pressure);
    }
    print_relative_metrics("b33_contact_pressure", pressure_error);
    passed = check(relative_metrics_below(pressure_error, 5.0e-3), "B3.3 contact pressure errors pass") && passed;
    std::size_t count = 0, active = 0, sliding = 0;
    double excess = 0.0, maximum_pressure = 0.0, tangential_force = 0.0;
    for (std::size_t n = 0; n < projected.size(); ++n) {
        if (std::isnan(projected[n])) continue;
        ++count;
        passed = check(projected[n] == 1.0 && seen[n], "B3.3 contact nodes are projected and compared") && passed;
        if (pressures[n] > 0.0) {
            ++active;
            if (output.nodal("contact_sliding_interface")[n] == 1.0) ++sliding;
            excess = std::max(excess, output.nodal("contact_tangential_traction_interface")[n] - 0.2 * pressures[n]);
            maximum_pressure = std::max(maximum_pressure, pressures[n]);
        }
        tangential_force += output.nodal("contact_tangential_force_interface")[n];
    }
    passed = check(count == 4 && active == 4 && sliding == 0 && contact.size() == 4,
                 "B3.3 has four active nodes in the sticking branch") &&
             passed;
    passed = check(excess <= 1.0e-12 * maximum_pressure, "B3.3 respects the Coulomb cap") && passed;
    passed = check(output.global("contact_heat_rate_interface") > 0.0 &&
                       output.global("contact_force_interface") > 0.0 && tangential_force > 0.0,
                 "B3.3 transfers nonzero heat, normal force and friction force") &&
             passed;
    double difference_squared = 0.0, scale_squared = 0.0;
    for (std::size_t n = 0; n < nodes.size(); ++n) {
        const double actual = output.nodal("displacement_y")[n];
        const double difference = actual - baseline.nodal("displacement_y")[n];
        difference_squared += difference * difference;
        scale_squared += actual * actual;
    }
    const double effect = std::sqrt(difference_squared / scale_squared);
    std::cout << "b33_frictional_tangential_field_change=" << effect << '\n';
    return check(std::isfinite(effect) && effect > 1.0e-3, "B3.3 friction measurably changes the tangential field") &&
           passed;
}
} // namespace fuelsim::test
