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
    std::array<double, 3> normal_force, tangential_force, tangent_first, tangent_second;
    double slip_1, slip_2, opening, pressure;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
        result.push_back(field);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size())
        throw std::invalid_argument("H20.36 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("H20.36 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max())
        || std::floor(result) != result)
        throw std::invalid_argument("H20.36 CSV contains an invalid integer index: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<std::vector<DisplacementReference>> read_displacement_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.36 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,time,id,x,y,z,disp_x,disp_y,disp_z")
        throw std::invalid_argument("Unexpected H20.36 displacement header in " + path);
    std::vector<std::vector<DisplacementReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.36 displacement step is invalid: " + path);
        result[step - 1].push_back({index_value(values, 2, path),
            {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
            {number(values, 6, path), number(values, 7, path), number(values, 8, path)}});
    }
    return result;
}

std::vector<std::vector<ContactReference>> read_contact_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.36 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,time,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,ctandir1_x,"
           "ctandir1_y,ctandir1_z,ctandir2_x,ctandir2_y,ctandir2_z,copen,cpress")
        throw std::invalid_argument("Unexpected H20.36 contact header in " + path);
    std::vector<std::vector<ContactReference>> result(7);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        const std::size_t step = index_value(values, 0, path);
        if (step == 0 || step > result.size() || std::abs(number(values, 1, path) - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("H20.36 contact step is invalid: " + path);
        result[step - 1].push_back({index_value(values, 2, path),
            {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
            {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
            {number(values, 9, path), number(values, 10, path), number(values, 11, path)},
            {number(values, 14, path), number(values, 15, path), number(values, 16, path)},
            {number(values, 17, path), number(values, 18, path), number(values, 19, path)},
            number(values, 12, path),
            number(values, 13, path),
            number(values, 20, path),
            number(values, 21, path)});
    }
    return result;
}

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double dot(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

bool compare_abaqus_steps(const std::vector<StepState>& states,
    const std::string& displacement_path,
    const std::string& contact_path) {
    const auto displacement = read_displacement_reference(displacement_path);
    const auto contact = read_contact_reference(contact_path);
    if (states.size() != displacement.size() || states.size() != contact.size())
        throw std::invalid_argument("H20.36 Fuelsim and Abaqus step counts differ");
    const auto& contact_nodes = states.front().source_nodes;
    constexpr double tolerance = 1.0e-2;
    constexpr double complete_slip_tolerance = 1.25e-2;
    constexpr double zero_tolerance = 1.2e-8;
    bool passed = true;
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<fuelsim::test::FieldErrorMetrics, 3> displacement_metrics;
        fuelsim::test::FieldErrorMetrics normal_force, circumferential_force, axial_force, slip_1, slip_2,
            slip_magnitude;
        fuelsim::test::GroupedFieldErrorMetrics displacement_vector, contact_force_vector, tangential_slip_vector;
        std::vector<bool> present(states[step].output.nodes.size(), false);
        std::vector<std::size_t> displacement_source_ids;
        double maximum_coordinate_difference = 0.0, maximum_tangent_basis_error = 0.0;
        if (displacement[step].size() != states[step].output.nodes.size())
            throw std::invalid_argument("Displacement reference must cover every production mesh node");
        for (std::size_t source = 0; source < states[step].output.nodes.size(); ++source) {
            const auto found = std::find_if(displacement[step].begin(),
                displacement[step].end(),
                [source](const DisplacementReference& value) { return value.id == source; });
            if (found == displacement[step].end() || present[source])
                throw std::invalid_argument("H20.36 displacement source-node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_difference = std::max({maximum_coordinate_difference,
                std::abs(states[step].output.nodes[source][0] - found->point.x),
                std::abs(states[step].output.nodes[source][1] - found->point.y),
                std::abs(states[step].output.nodes[source][2] - found->point.z)});
            const double radius = std::hypot(found->point.x, found->point.y);
            const double actual_x = states[step].output.nodal("displacement_x").at(source);
            const double actual_y = states[step].output.nodal("displacement_y").at(source);
            const double actual_z = states[step].output.nodal("displacement_z").at(source);
            const std::array<double, 3>
                actual_displacement = {(found->point.x * actual_x + found->point.y * actual_y) / radius,
                    (-found->point.y * actual_x + found->point.x * actual_y) / radius,
                    actual_z},
                reference_displacement = {
                    (found->point.x * found->displacement[0] + found->point.y * found->displacement[1]) / radius,
                    (-found->point.y * found->displacement[0] + found->point.x * found->displacement[1]) / radius,
                    found->displacement[2]};
            for (std::size_t component = 0; component < 3; ++component)
                displacement_metrics[component].add(actual_displacement[component], reference_displacement[component]);
            displacement_vector.add(actual_displacement.data(), reference_displacement.data(), 3);
            displacement_source_ids.push_back(source);
        }
        if (contact_nodes.size() != states[step].contact.size() || contact_nodes.size() != contact[step].size())
            throw std::invalid_argument("H20.36 contact-node counts differ");
        double actual_normal_resultant = 0.0, reference_normal_resultant = 0.0;
        double actual_circumferential_resultant = 0.0, reference_circumferential_resultant = 0.0;
        double actual_axial_resultant = 0.0, reference_axial_resultant = 0.0;
        bool abaqus_contact_closed = true;
        for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
            const auto found = std::find_if(contact[step].begin(),
                contact[step].end(),
                [&](const ContactReference& value) { return value.id == contact_nodes[node]; });
            if (found == contact[step].end())
                throw std::invalid_argument("H20.36 contact-node mapping is incomplete");
            const fuelsim::test::OutputContact& actual = states[step].contact[node];
            maximum_coordinate_difference = std::max({maximum_coordinate_difference,
                std::abs(actual.x - found->point.x),
                std::abs(actual.y - found->point.y),
                std::abs(actual.z - found->point.z)});
            const double radius = std::hypot(found->point.x, found->point.y);
            const std::array<double, 3> radial = {found->point.x / radius, found->point.y / radius, 0.0};
            const std::array<double, 3> circumferential = {-found->point.y / radius, found->point.x / radius, 0.0};
            const std::array<double, 3> actual_normal = {-actual.normal_contact_force[0],
                -actual.normal_contact_force[1],
                -actual.normal_contact_force[2]};
            const std::array<double, 3> actual_tangent = {-actual.tangential_contact_force[0],
                -actual.tangential_contact_force[1],
                -actual.tangential_contact_force[2]};
            const double actual_normal_value = dot(actual_normal, radial);
            const double reference_normal_value = dot(found->normal_force, radial);
            const double actual_circumferential = dot(actual_tangent, circumferential);
            const double reference_circumferential = dot(found->tangential_force, circumferential);
            normal_force.add(actual_normal_value, reference_normal_value);
            circumferential_force.add(actual_circumferential, reference_circumferential);
            axial_force.add(actual_tangent[2], found->tangential_force[2]);
            // For these C3D20 S6 faces, Abaqus local direction 1 opposes the
            // positive circumferential direction and direction 2 opposes the
            // positive global axial direction.
            slip_1.add(-dot(actual.tangential_slip, circumferential), found->slip_1);
            slip_2.add(-actual.tangential_slip[2], found->slip_2);
            const std::array<double, 3> actual_force = {actual_normal_value, actual_circumferential, actual_tangent[2]},
                                        reference_force = {reference_normal_value,
                                            reference_circumferential,
                                            found->tangential_force[2]};
            std::array<double, 3> reference_slip_vector{};
            for (std::size_t component = 0; component < 3; ++component)
                reference_slip_vector[component] =
                    found->slip_1 * found->tangent_first[component] + found->slip_2 * found->tangent_second[component];
            contact_force_vector.add(actual_force.data(), reference_force.data(), 3);
            tangential_slip_vector.add(actual.tangential_slip.data(), reference_slip_vector.data(), 3);
            slip_magnitude.add(
                std::hypot(actual.tangential_slip[0], actual.tangential_slip[1], actual.tangential_slip[2]),
                std::hypot(reference_slip_vector[0], reference_slip_vector[1], reference_slip_vector[2]));
            maximum_tangent_basis_error = std::max({maximum_tangent_basis_error,
                std::abs(dot(found->tangent_first, found->tangent_first) - 1.0),
                std::abs(dot(found->tangent_second, found->tangent_second) - 1.0),
                std::abs(dot(found->tangent_first, found->tangent_second))});
            actual_normal_resultant += actual_normal_value;
            reference_normal_resultant += reference_normal_value;
            actual_circumferential_resultant += actual_circumferential;
            reference_circumferential_resultant += reference_circumferential;
            actual_axial_resultant += actual_tangent[2];
            reference_axial_resultant += found->tangential_force[2];
            abaqus_contact_closed = abaqus_contact_closed && found->pressure > 0.0 && found->opening < 0.0;
        }
        const std::string prefix = "h20_36_step_" + std::to_string(step + 1) + "_";
        const double circumferential_displacement_difference_on_axial_scale =
            displacement_metrics[1].maximum_absolute_difference / displacement_metrics[2].maximum_reference;
        const double circumferential_force_difference_on_axial_scale =
            circumferential_force.maximum_absolute_difference / axial_force.maximum_reference;
        const double circumferential_slip_difference_on_axial_scale =
            slip_1.maximum_absolute_difference / slip_2.maximum_reference;
        fuelsim::test::print_relative_metrics(prefix + "displacement_radial", displacement_metrics[0]);
        if (displacement_metrics[1].has_relative_norm())
            fuelsim::test::print_relative_metrics(prefix + "displacement_circumferential", displacement_metrics[1]);
        fuelsim::test::print_absolute_metrics(prefix + "displacement_circumferential", displacement_metrics[1]);
        std::cout << prefix << "displacement_circumferential_maximum_difference_on_axial_scale="
                  << circumferential_displacement_difference_on_axial_scale << '\n';
        fuelsim::test::print_relative_metrics(prefix + "displacement_axial", displacement_metrics[2]);
        fuelsim::test::print_grouped_relative_metrics(prefix + "displacement_vector", displacement_vector);
        std::cout << prefix << "displacement_vector_maximum_pointwise_source_node="
                  << displacement_source_ids.at(displacement_vector.maximum_pointwise_relative_index) << '\n';
        passed =
            check(fuelsim::test::grouped_relative_metrics_below(displacement_vector, tolerance)
                      && displacement_vector.maximum_zero_reference_difference < zero_tolerance,
                "H20.36 step " + std::to_string(step + 1) + "complete displacement-vector metrics are below 1 percent")
            && passed;
        fuelsim::test::print_relative_metrics(prefix + "signed_normal_radial_force", normal_force);
        fuelsim::test::print_relative_metrics(prefix + "signed_tangential_circumferential_force",
            circumferential_force);
        fuelsim::test::print_relative_metrics(prefix + "signed_tangential_axial_force", axial_force);
        if (slip_1.has_relative_norm())
            fuelsim::test::print_relative_metrics(prefix + "tangential_slip_1", slip_1);
        fuelsim::test::print_absolute_metrics(prefix + "signed_tangential_circumferential_force",
            circumferential_force);
        fuelsim::test::print_absolute_metrics(prefix + "tangential_slip_1", slip_1);
        std::cout << prefix << "signed_tangential_circumferential_force_maximum_difference_on_axial_scale="
                  << circumferential_force_difference_on_axial_scale << '\n'
                  << prefix << "tangential_slip_1_maximum_difference_on_axial_scale="
                  << circumferential_slip_difference_on_axial_scale << '\n';
        fuelsim::test::print_relative_metrics(prefix + "tangential_slip_2", slip_2);
        fuelsim::test::print_grouped_relative_metrics(prefix + "contact_force_vector", contact_force_vector);
        fuelsim::test::print_grouped_relative_metrics(prefix + "tangential_slip_vector", tangential_slip_vector);
        fuelsim::test::print_relative_metrics(prefix + "tangential_slip_magnitude", slip_magnitude);
        std::cout << prefix << "contact_force_vector_maximum_pointwise_source_node="
                  << contact_nodes.at(contact_force_vector.maximum_pointwise_relative_index) << '\n'
                  << prefix << "tangential_slip_vector_maximum_pointwise_source_node="
                  << contact_nodes.at(tangential_slip_vector.maximum_pointwise_relative_index) << '\n';
        passed =
            check(maximum_coordinate_difference < 1.0e-12,
                "H20.36 step " + std::to_string(step + 1) + " uses the tracked quadratic Exodus mesh")
            && check(maximum_tangent_basis_error < 1.0e-12,
                "H20.36 step " + std::to_string(step + 1) + " uses an orthonormal Abaqus contact tangent basis")
            && check(abaqus_contact_closed,
                "H20.36 step " + std::to_string(step + 1) + " keeps all Abaqus contact nodes closed")
            && check(fuelsim::test::grouped_relative_metrics_below(contact_force_vector, tolerance)
                         && contact_force_vector.maximum_zero_reference_difference < zero_tolerance,
                "H20.36 step " + std::to_string(step + 1) + "complete contact-force-vector metrics are below 1 percent")
            && check(fuelsim::test::grouped_relative_metrics_below(tangential_slip_vector, complete_slip_tolerance)
                         && tangential_slip_vector.maximum_zero_reference_difference < zero_tolerance,
                "H20.36 step " + std::to_string(step + 1)
                    + "complete tangential-slip-vector metrics are below 1.25 percent")
            && passed;
        const double normal_resultant_error =
            std::abs(actual_normal_resultant - reference_normal_resultant) / std::abs(reference_normal_resultant);
        const double circumferential_resultant_difference =
            std::abs(actual_circumferential_resultant - reference_circumferential_resultant);
        const double circumferential_resultant_error =
            circumferential_resultant_difference / std::abs(reference_circumferential_resultant);
        const double axial_resultant_error =
            std::abs(actual_axial_resultant - reference_axial_resultant) / std::abs(reference_axial_resultant);
        const double circumferential_resultant_difference_on_axial_scale =
            circumferential_resultant_difference / std::abs(reference_axial_resultant);
        std::cout << prefix << "normal_resultant_relative_error=" << normal_resultant_error << '\n'
                  << prefix << "circumferential_resultant_relative_error=" << circumferential_resultant_error << '\n'
                  << prefix << "circumferential_resultant_absolute_difference=" << circumferential_resultant_difference
                  << '\n'
                  << prefix << "circumferential_resultant_difference_on_axial_scale="
                  << circumferential_resultant_difference_on_axial_scale << '\n'
                  << prefix << "axial_resultant_relative_error=" << axial_resultant_error << '\n';
        passed = check(normal_resultant_error < tolerance && axial_resultant_error < tolerance,
                     "H20.36 step " + std::to_string(step + 1)
                         + "signed normal and driven axial resultants agree with Abaqus below 1 percent")
                 && check(circumferential_resultant_difference_on_axial_scale < tolerance,
                     "H20.36 step " + std::to_string(step + 1)
                         + "undriven circumferential resultant difference is below 1 percent of the driven axial "
                           "resultant")
                 && passed;
    }
    return passed;
}

std::array<std::size_t, 2> stick_slide_counts(const StepState& state) {
    std::array<std::size_t, 2> result{};
    for (const fuelsim::test::OutputContact& node : state.contact) {
        if (!(node.pressure > 0.0))
            continue;
        ++result[node.sliding ? 1 : 0];
    }
    return result;
}

std::array<double, 2> tangential_resultant(const StepState& state) {
    std::array<double, 2> result{};
    for (const fuelsim::test::OutputContact& node : state.contact) {
        const double radius = std::hypot(node.x, node.y);
        result[0] += (-node.y * node.tangential_contact_force[0] + node.x * node.tangential_contact_force[1]) / radius;
        result[1] += node.tangential_contact_force[2];
    }
    return result;
}

} // namespace

namespace fuelsim::test {
bool check_hex20_friction_path_36(const std::string& output_path,
    const std::string& displacement_path,
    const std::string& contact_path) {
    const auto states = read_seven_contact_steps(output_path, 37);
    bool passed = true;
    if (states.size() == 7) {
        const auto first = stick_slide_counts(states[0]);
        const auto second = stick_slide_counts(states[1]);
        const auto mixed = stick_slide_counts(states[2]);
        const auto forward = stick_slide_counts(states[3]);
        const auto unload = stick_slide_counts(states[4]);
        const auto reverse = stick_slide_counts(states[5]);
        const auto restick = stick_slide_counts(states[6]);
        const std::array<double, 2> forward_resultant = tangential_resultant(states[3]);
        const std::array<double, 2> reverse_resultant = tangential_resultant(states[5]);
        passed = check(std::all_of(states.begin(),
                           states.end(),
                           [](const StepState& state) {
                               const auto count = stick_slide_counts(state);
                               return count[0] + count[1] == 37;
                           }),
                     "H20.36 keeps all 37 curved constraints closed throughout the path")
                 && check(first[0] == 37 && first[1] == 0 && second[0] == 37 && second[1] == 0,
                     "H20.36 starts with two fully sticking states")
                 && check(mixed[0] > 0 && mixed[1] > 0 && forward[1] > 0,
                     "H20.36 contains simultaneous sticking and sliding before forward sliding")
                 && check(unload[0] > 0 && unload[1] > 0 && reverse[0] == 0 && reverse[1] == 37
                              && forward_resultant[1] * reverse_resultant[1] < 0.0,
                     "H20.36 covers mixed unloading, committed full reverse sliding, and reverses the driven axial "
                     "tangential resultant")
                 && check(restick[0] == 37 && restick[1] == 0, "H20.36 finishes with all 37 constraints restuck")
                 && passed;
    }

    return compare_abaqus_steps(states, displacement_path, contact_path) && passed;
}
} // namespace fuelsim::test
