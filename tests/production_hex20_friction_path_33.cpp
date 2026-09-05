#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include "support/production_contact_history.hpp"
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
using StepState = fuelsim::test::OutputContactStep;

struct DisplacementReference final {
    std::size_t id;
    fuelsim::test::OutputPoint point;
    std::array<double, 3> displacement;
};

struct ContactReference final {
    std::size_t id;
    fuelsim::test::OutputPoint point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2;
};

struct ReactionReference final {
    std::array<double, 3> reaction;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) result.push_back(field);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("H20.33 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("H20.33 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("H20.33 CSV contains an invalid integer index: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<std::vector<DisplacementReference>> read_displacement_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.33 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,time,id,x,y,z,disp_x,disp_y,disp_z")
        throw std::invalid_argument("Unexpected H20.33 displacement header in " + path);
    std::vector<std::vector<DisplacementReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.33 displacement step is invalid: " + path);
        result[step - 1].push_back(
            {index_value(values, 2, path), {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 6, path), number(values, 7, path), number(values, 8, path)}});
    }
    return result;
}

std::vector<std::vector<ContactReference>> read_contact_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.33 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line !=
        "step,time,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected H20.33 contact header in " + path);
    std::vector<std::vector<ContactReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.33 contact step is invalid: " + path);
        result[step - 1].push_back(
            {index_value(values, 2, path), {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
                {number(values, 9, path), number(values, 10, path), number(values, 11, path)}, number(values, 12, path),
                number(values, 13, path)});
    }
    return result;
}

std::vector<ReactionReference> read_reaction_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.33 reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,time,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected H20.33 reaction header in " + path);
    std::vector<ReactionReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step != result.size() + 1 || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.33 reaction steps are not consecutive: " + path);
        result.push_back({{number(values, 2, path), number(values, 3, path), number(values, 4, path)}});
    }
    return result;
}

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool compare_abaqus_steps(const std::vector<StepState>& states, const std::string& displacement_path,
    const std::string& contact_path, const std::string& reaction_path) {
    const auto displacement = read_displacement_reference(displacement_path);
    const auto contact = read_contact_reference(contact_path);
    const auto reaction = read_reaction_reference(reaction_path);
    if (states.size() != displacement.size() || states.size() != contact.size() || states.size() != reaction.size())
        throw std::invalid_argument("H20.33 Fuelsim and Abaqus step counts differ");
    const auto& contact_nodes = states.front().source_nodes;
    constexpr double tolerance = 1.0e-2;
    constexpr double zero_tolerance = 1.2e-8;
    bool passed = true;
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_metrics;
        fuelsim::test::FieldErrorMetrics normal_force, tangential_force, slip_1, slip_2;
        fuelsim::test::GroupedFieldErrorMetrics displacement_vector, contact_force_vector, tangential_slip_vector;
        std::vector<bool> present(states[step].output.nodes.size(), false);
        std::vector<std::size_t> displacement_source_ids;
        double maximum_coordinate_difference = 0.0;
        if (displacement[step].size() != states[step].output.nodes.size())
            throw std::invalid_argument("Displacement reference must cover every production mesh node");
        for (std::size_t source = 0; source < states[step].output.nodes.size(); ++source) {
            const auto found = std::find_if(displacement[step].begin(), displacement[step].end(),
                [source](const DisplacementReference& value) { return value.id == source; });
            if (found == displacement[step].end() || present[source])
                throw std::invalid_argument("H20.33 displacement source-node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_difference = std::max(
                {maximum_coordinate_difference, std::abs(states[step].output.nodes[source][0] - found->point.x),
                    std::abs(states[step].output.nodes[source][1] - found->point.y),
                    std::abs(states[step].output.nodes[source][2] - found->point.z)});
            std::array<double, 3> actual_displacement{}, reference_displacement{};
            for (std::size_t component = 0; component < 3; ++component) {
                actual_displacement[component] =
                    states[step].output.nodal("displacement_" + std::string(1, "xyz"[component])).at(source);
                reference_displacement[component] = found->displacement[component];
                displacement_metrics[component].add(actual_displacement[component], reference_displacement[component]);
            }
            displacement_vector.add(actual_displacement.data(), reference_displacement.data(), 3);
            displacement_source_ids.push_back(source);
        }
        if (contact_nodes.size() != states[step].contact.size() || contact_nodes.size() != contact[step].size())
            throw std::invalid_argument("H20.33 contact-node counts differ");
        for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
            const auto found = std::find_if(contact[step].begin(), contact[step].end(),
                [&](const ContactReference& value) { return value.id == contact_nodes[node]; });
            if (found == contact[step].end()) throw std::invalid_argument("H20.33 contact-node mapping is incomplete");
            const fuelsim::test::OutputContact& actual = states[step].contact[node];
            maximum_coordinate_difference =
                std::max({maximum_coordinate_difference, std::abs(actual.x - found->point.x),
                    std::abs(actual.y - found->point.y), std::abs(actual.z - found->point.z)});
            // Abaqus contact-force output uses the opposite contact-pair action
            // convention from the Fuelsim secondary-side summary.
            normal_force.add(-actual.normal_contact_force[0], found->normal_force[0]);
            tangential_force.add(-actual.tangential_contact_force[1], found->tangential_force[1]);
            slip_1.add(actual.tangential_slip[2], found->slip_1);
            slip_2.add(actual.tangential_slip[1], -found->slip_2);
            const std::array<double, 2> actual_force = {-actual.normal_contact_force[0],
                                            -actual.tangential_contact_force[1]},
                                        reference_force = {found->normal_force[0], found->tangential_force[1]},
                                        actual_slip = {actual.tangential_slip[2], actual.tangential_slip[1]},
                                        reference_slip = {found->slip_1, -found->slip_2};
            contact_force_vector.add(actual_force.data(), reference_force.data(), 2);
            tangential_slip_vector.add(actual_slip.data(), reference_slip.data(), 2);
        }
        const std::string prefix = "h20_33_step_" + std::to_string(step + 1) + "_";
        for (std::size_t component = 0; component < 3; ++component) {
            const std::string name = prefix + "displacement_" + std::string(1, "xyz"[component]);
            if (displacement_metrics[component].has_relative_norm())
                fuelsim::test::print_relative_metrics(name, displacement_metrics[component]);
            else
                fuelsim::test::print_absolute_metrics(name, displacement_metrics[component]);
        }
        fuelsim::test::print_relative_metrics(prefix + "signed_normal_force_x", normal_force);
        fuelsim::test::print_relative_metrics(prefix + "signed_tangential_force_y", tangential_force);
        if (slip_1.has_relative_norm())
            fuelsim::test::print_relative_metrics(prefix + "tangential_slip_1", slip_1);
        else
            fuelsim::test::print_absolute_metrics(prefix + "tangential_slip_1", slip_1);
        fuelsim::test::print_relative_metrics(prefix + "tangential_slip_2", slip_2);
        fuelsim::test::print_grouped_relative_metrics(prefix + "displacement_vector", displacement_vector);
        fuelsim::test::print_grouped_relative_metrics(prefix + "contact_force_vector", contact_force_vector);
        fuelsim::test::print_grouped_relative_metrics(prefix + "tangential_slip_vector", tangential_slip_vector);
        std::cout << prefix << "displacement_vector_maximum_pointwise_source_node="
                  << displacement_source_ids.at(displacement_vector.maximum_pointwise_relative_index) << '\n'
                  << prefix << "contact_force_vector_maximum_pointwise_source_node="
                  << contact_nodes.at(contact_force_vector.maximum_pointwise_relative_index) << '\n'
                  << prefix << "tangential_slip_vector_maximum_pointwise_source_node="
                  << contact_nodes.at(tangential_slip_vector.maximum_pointwise_relative_index) << '\n';
        passed = check(maximum_coordinate_difference < 1.0e-7,
                     "H20.33 step " + std::to_string(step + 1) + " uses the tracked nonmatching Exodus mesh") &&
                 check(fuelsim::test::grouped_relative_metrics_below(displacement_vector, tolerance) &&
                           displacement_vector.maximum_zero_reference_difference < zero_tolerance,
                     "H20.33 step " + std::to_string(step + 1) +
                         " complete displacement-vector metrics are below 1 percent") &&
                 check(fuelsim::test::grouped_relative_metrics_below(contact_force_vector, tolerance) &&
                           contact_force_vector.maximum_zero_reference_difference < zero_tolerance,
                     "H20.33 step " + std::to_string(step + 1) +
                         " complete contact-force-vector metrics are below 1 percent") &&
                 check(fuelsim::test::grouped_relative_metrics_below(tangential_slip_vector, tolerance) &&
                           tangential_slip_vector.maximum_zero_reference_difference < zero_tolerance,
                     "H20.33 step " + std::to_string(step + 1) +
                         " complete tangential-slip-vector metrics are below 1 percent") &&
                 passed;
        double actual_normal = 0.0, actual_tangential = 0.0;
        for (const fuelsim::test::OutputContact& node : states[step].contact) {
            actual_normal += node.normal_contact_force[0];
            actual_tangential += node.tangential_contact_force[1];
        }
        const double normal_resultant_error =
                         std::abs(actual_normal - reaction[step].reaction[0]) / std::abs(reaction[step].reaction[0]),
                     tangential_resultant_error = std::abs(actual_tangential - reaction[step].reaction[1]) /
                                                  std::abs(reaction[step].reaction[1]);
        std::cout << prefix << "normal_resultant_relative_error=" << normal_resultant_error << '\n'
                  << prefix << "tangential_resultant_relative_error=" << tangential_resultant_error << '\n';
        passed = check(normal_resultant_error < tolerance && tangential_resultant_error < tolerance,
                     "H20.33 step " + std::to_string(step + 1) +
                         " signed normal and tangential resultants agree with Abaqus below 1 percent") &&
                 passed;
    }
    return passed;
}

std::array<std::size_t, 2> stick_slide_counts(const StepState& state) {
    std::array<std::size_t, 2> result{};
    for (const fuelsim::test::OutputContact& node : state.contact) {
        if (!(node.pressure > 0.0)) continue;
        ++result[node.sliding ? 1 : 0];
    }
    return result;
}

double tangential_resultant_y(const StepState& state) {
    double result = 0.0;
    for (const fuelsim::test::OutputContact& node : state.contact) result += node.tangential_contact_force[1];
    return result;
}

} // namespace

namespace fuelsim::test {
bool check_hex20_friction_path_33(const std::string& output_path, const std::string& displacement_path,
    const std::string& contact_path, const std::string& reaction_path) {
    const auto states = read_seven_contact_steps(output_path, 13);
    bool passed = true;
    if (states.size() == 7) {
        const auto first = stick_slide_counts(states[0]);
        const auto second = stick_slide_counts(states[1]);
        const auto forward_entry = stick_slide_counts(states[2]);
        const auto forward = stick_slide_counts(states[3]);
        const auto unload = stick_slide_counts(states[4]);
        const auto reverse = stick_slide_counts(states[5]);
        const auto restick = stick_slide_counts(states[6]);
        passed = check(first[0] == 13 && first[1] == 0 && second[0] == 13 && second[1] == 0,
                     "H20 friction path starts with two fully sticking states") &&
                 check(forward_entry[0] == 0 && forward_entry[1] == 13 && forward[0] == 0 && forward[1] == 13,
                     "H20 friction path enters and remains in committed forward sliding") &&
                 check(unload[0] == 0 && unload[1] == 13 && reverse[0] == 0 && reverse[1] == 13 &&
                           tangential_resultant_y(states[4]) < 0.0 && tangential_resultant_y(states[5]) < 0.0,
                     "H20 friction path covers committed unloading and reverse sliding with reversed signed force") &&
                 check(restick[0] == 13 && restick[1] == 0,
                     "H20 friction path finishes with all thirteen constraints restuck") &&
                 passed;
    }

    return compare_abaqus_steps(states, displacement_path, contact_path, reaction_path) && passed;
}
} // namespace fuelsim::test
