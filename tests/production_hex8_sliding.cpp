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
struct NodeReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> displacement;
};

struct ContactReference final {
    std::size_t id;
    std::array<double, 3> point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2, opening, pressure;
};

struct ReactionReference final {
    std::array<double, 3> reaction;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("B3.4 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("B3.4 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("B3.4 CSV contains an invalid integer: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.4 Abaqus nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "id,x,y,z,u1,u2,u3")
        throw std::invalid_argument("Unexpected B3.4 nodal header: " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)}});
    }
    if (result.size() != 16) throw std::invalid_argument("B3.4 Abaqus reference must contain sixteen nodes");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.4 Abaqus contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected B3.4 contact header: " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back(
            {index_value(values, 0, path), {number(values, 1, path), number(values, 2, path), number(values, 3, path)},
                {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
                {number(values, 7, path), number(values, 8, path), number(values, 9, path)}, number(values, 10, path),
                number(values, 11, path), number(values, 12, path), number(values, 13, path)});
    }
    if (result.size() != 4) throw std::invalid_argument("B3.4 Abaqus reference must contain four contact nodes");
    return result;
}

ReactionReference read_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.4 Abaqus reaction reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "reaction_x,reaction_y,reaction_z" || !std::getline(input, line))
        throw std::invalid_argument("Unexpected B3.4 reaction reference: " + path);
    const std::vector<std::string> values = split_csv(line);
    return {{{number(values, 0, path), number(values, 1, path), number(values, 2, path)}}};
}

struct InterfaceValues final {
    double total_contact_force = 0.0, total_tangential_force = 0.0;
    std::size_t active_contact_nodes = 0;
};
} // namespace

namespace fuelsim::test {
bool check_hex8_sliding(const std::string& output_path, const std::string& nodes_path, const std::string& contact_path,
    const std::string& reaction_path) {
    const auto output = read_final_exodus_results(output_path);
    const auto nodes = read_nodes(nodes_path);
    if (output.nodes.size() != nodes.size()) throw std::runtime_error("B3.4 output node count differs");
    const std::array<std::string, 3> names = {"displacement_x", "displacement_y", "displacement_z"};
    std::vector<bool> present(nodes.size(), false);
    std::array<FieldErrorMetrics, 3> displacement;
    double maximum_coordinate_difference = 0.0;
    for (const auto& row : nodes) {
        if (row.id >= present.size() || present[row.id]) throw std::runtime_error("B3.4 invalid node ID");
        present[row.id] = true;
        for (std::size_t component = 0; component < 3; ++component) {
            maximum_coordinate_difference = std::max(
                maximum_coordinate_difference, std::abs(output.nodes[row.id][component] - row.point[component]));
            displacement[component].add(output.nodal(names[component]).at(row.id), row.displacement[component]);
        }
    }
    const auto reference = read_contact(contact_path);
    std::vector<std::size_t> source_nodes;
    InterfaceValues interface;
    const auto& projected = output.nodal("contact_projected_interface");
    for (std::size_t node = 0; node < projected.size(); ++node) {
        if (std::isnan(projected[node])) continue;
        if (projected[node] != 1.0) throw std::runtime_error("B3.4 lost a contact projection");
        source_nodes.push_back(node);
        if (output.nodal("contact_pressure_interface")[node] > 0.0) ++interface.active_contact_nodes;
        interface.total_contact_force += output.nodal("contact_normal_force_interface")[node];
        interface.total_tangential_force += output.nodal("contact_tangential_force_interface")[node];
    }
    if (source_nodes.size() != reference.size()) throw std::runtime_error("B3.4 contact node counts differ");
    fuelsim::test::FieldErrorMetrics normal_force, tangential_force_y, slip_y, opening, pressure;
    double maximum_theoretical_zero_force = 0.0;
    std::size_t active = 0, sliding = 0;
    for (std::size_t node = 0; node < source_nodes.size(); ++node) {
        const auto found = std::find_if(reference.begin(), reference.end(),
            [&](const ContactReference& value) { return value.id == source_nodes[node]; });
        if (found == reference.end()) throw std::invalid_argument("B3.4 Abaqus contact-node mapping is incomplete");
        maximum_coordinate_difference = std::max(
            maximum_coordinate_difference, std::max({std::abs(output.nodes[source_nodes[node]][0] - found->point[0]),
                                               std::abs(output.nodes[source_nodes[node]][1] - found->point[1]),
                                               std::abs(output.nodes[source_nodes[node]][2] - found->point[2])}));
        normal_force.add(-output.nodal("contact_normal_force_x_interface")[source_nodes[node]], found->normal_force[0]);
        tangential_force_y.add(
            -output.nodal("contact_tangential_force_y_interface")[source_nodes[node]], found->tangential_force[1]);
        slip_y.add(output.nodal("contact_total_slip_y_interface")[source_nodes[node]], -found->slip_2);
        opening.add(output.nodal("contact_gap_interface")[source_nodes[node]], found->opening);
        pressure.add(output.nodal("contact_pressure_interface")[source_nodes[node]], found->pressure);
        maximum_theoretical_zero_force = std::max(maximum_theoretical_zero_force,
            std::max({std::abs(output.nodal("contact_normal_force_y_interface")[source_nodes[node]]),
                std::abs(output.nodal("contact_normal_force_z_interface")[source_nodes[node]]),
                std::abs(output.nodal("contact_tangential_force_x_interface")[source_nodes[node]]),
                std::abs(output.nodal("contact_tangential_force_z_interface")[source_nodes[node]]),
                std::abs(found->normal_force[1]), std::abs(found->normal_force[2]),
                std::abs(found->tangential_force[0]), std::abs(found->tangential_force[2])}));
        if (output.nodal("contact_pressure_interface")[source_nodes[node]] > 0.0) ++active;
        if (output.nodal("contact_sliding_interface")[source_nodes[node]]) ++sliding;
    }

    const ReactionReference reaction = read_reaction(reaction_path);
    const double normal_resultant_error =
        std::abs(interface.total_contact_force - std::abs(reaction.reaction[0])) / std::abs(reaction.reaction[0]);
    const double tangential_resultant_error =
        std::abs(interface.total_tangential_force - std::abs(reaction.reaction[1])) / std::abs(reaction.reaction[1]);

    fuelsim::test::print_relative_metrics("b34_displacement_x", displacement[0]);
    fuelsim::test::print_relative_metrics("b34_displacement_y", displacement[1]);
    fuelsim::test::print_absolute_metrics("b34_displacement_z", displacement[2]);
    fuelsim::test::print_relative_metrics("b34_signed_normal_force_x", normal_force);
    fuelsim::test::print_relative_metrics("b34_signed_tangential_force_y", tangential_force_y);
    fuelsim::test::print_relative_metrics("b34_tangential_slip_y", slip_y);
    fuelsim::test::print_relative_metrics("b34_opening", opening);
    fuelsim::test::print_relative_metrics("b34_pressure", pressure);
    std::cout << "b34_displacement_z_zero_reference_count=" << displacement[2].zero_reference_count << '\n'
              << "b34_displacement_z_maximum_absolute_difference=" << displacement[2].maximum_absolute_difference
              << '\n'
              << "b34_theoretical_zero_contact_force_component_value_count=32\n"
              << "b34_theoretical_zero_contact_force_maximum_absolute_difference=" << maximum_theoretical_zero_force
              << '\n'
              << "b34_normal_resultant_relative_error=" << normal_resultant_error << '\n'
              << "b34_tangential_resultant_relative_error=" << tangential_resultant_error << '\n'
              << "b34_maximum_mesh_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << "b34_active_contact_nodes=" << active << '\n'
              << "b34_sliding_contact_nodes=" << sliding << '\n';

    constexpr double tolerance = 5.0e-3;
    return check(active == 4 && sliding > 0 && interface.active_contact_nodes == 4,
               "B3.4 keeps four active constraints and reaches the Coulomb sliding branch") &&
           check(fuelsim::test::relative_metrics_below(displacement[0], tolerance) &&
                     fuelsim::test::relative_metrics_below(displacement[1], tolerance),
               "B3.4 in-plane nodal displacements pass all three Abaqus metrics below 0.5 percent") &&
           check(displacement[2].maximum_absolute_difference < 5.0e-10,
               "B3.4 out-of-plane symmetry displacement agrees with Abaqus in absolute value") &&
           check(fuelsim::test::relative_metrics_below(normal_force, tolerance) &&
                     fuelsim::test::relative_metrics_below(tangential_force_y, tolerance),
               "B3.4 signed in-plane nodal contact-force fields pass all three Abaqus metrics below 0.5 percent") &&
           check(fuelsim::test::relative_metrics_below(slip_y, tolerance) &&
                     fuelsim::test::relative_metrics_below(opening, tolerance) &&
                     fuelsim::test::relative_metrics_below(pressure, tolerance),
               "B3.4 slip, opening, and pressure pass all three Abaqus metrics below 0.5 percent") &&
           check(maximum_theoretical_zero_force < 1.0e-4,
               "B3.4 theoretical-zero transverse contact-force components pass their absolute check") &&
           check(normal_resultant_error < tolerance && tangential_resultant_error < tolerance,
               "B3.4 normal and tangential resultants agree with Abaqus below 0.5 percent") &&
           check(maximum_coordinate_difference < 5.0e-10,
               "B3.4 Abaqus and Fuelsim use matching tracked HEX8 coordinates") &&
           check(std::abs(interface.total_tangential_force - 0.001 * interface.total_contact_force) <
                     1.0e-12 * interface.total_contact_force,
               "B3.4 Fuelsim sliding resultant lies on the Coulomb cap");
}

} // namespace fuelsim::test
