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
struct ContactReference final {
    std::size_t step, id;
    std::array<double, 3> point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2, opening, pressure;
};

struct ReactionReference final {
    std::size_t step;
    std::array<double, 3> reaction;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("B4.0 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("B4.0 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("B4.0 CSV contains an invalid integer: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.0 contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected B4.0 contact header: " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path), index_value(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}, number(values, 11, path),
            number(values, 12, path), number(values, 13, path), number(values, 14, path)});
    }
    if (result.size() != 16) throw std::invalid_argument("B4.0 contact reference must contain sixteen rows");
    return result;
}

std::vector<ReactionReference> read_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.0 reaction reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "step,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected B4.0 reaction header: " + path);
    std::vector<ReactionReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path),
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}});
    }
    if (result.size() != 4) throw std::invalid_argument("B4.0 reaction reference must contain four rows");
    return result;
}

bool compare(const std::string& output_path, const std::string& contact_path, const std::string& reaction_path) {
    const auto final = fuelsim::test::read_final_exodus_results(output_path);
    if (final.step_count != 5 || final.time != 4.0)
        throw std::invalid_argument("B4.0 output must contain the initial state and four physical time steps");
    const std::vector<ContactReference> reference = read_contact(contact_path);
    const std::vector<ReactionReference> reactions = read_reaction(reaction_path);
    fuelsim::test::FieldErrorMetrics normal_force, tangential_y, tangential_z, slip_y, slip_z, opening, pressure;
    double maximum_coordinate_error = 0.0, maximum_zero_force = 0.0, maximum_resultant_error = 0.0;
    bool all_active = true, first_sticks = true, later_path_reaches_sliding = true, biaxial_history = true,
         later_path_is_on_coulomb_circle = true;
    for (std::size_t step = 0; step < 4; ++step) {
        const auto frame = fuelsim::test::read_exodus_results(output_path, step + 2);
        if (frame.time != static_cast<double>(step + 1))
            throw std::invalid_argument("B4.0 output time does not match the tracked reference step");
        const auto& projected = frame.nodal("contact_projected_interface");
        std::vector<std::size_t> sources;
        for (std::size_t source = 0; source < projected.size(); ++source)
            if (std::isfinite(projected[source])) sources.push_back(source);
        if (sources.size() != 4) throw std::invalid_argument("B4.0 requires four contact constraints");
        std::array<double, 3> actual_resultant{};
        std::size_t sticking_count = 0, sliding_count = 0;
        for (std::size_t node = 0; node < sources.size(); ++node) {
            const auto found = std::find_if(reference.begin(), reference.end(),
                [&](const ContactReference& value) { return value.step == step + 1 && value.id == sources[node]; });
            if (found == reference.end()) throw std::invalid_argument("B4.0 contact-node mapping is incomplete");
            const auto source = sources[node];
            const auto scalar = [&](const std::string& field) {
                const double value = frame.nodal("contact_" + field + "_interface").at(source);
                if (!std::isfinite(value)) throw std::invalid_argument("B4.0 contact output is not finite");
                return value;
            };
            const std::array<double, 3> normal = {
                scalar("normal_force_x"), scalar("normal_force_y"), scalar("normal_force_z")};
            const std::array<double, 3> tangent = {
                scalar("tangential_force_x"), scalar("tangential_force_y"), scalar("tangential_force_z")};
            const std::array<double, 3> slip = {scalar("total_slip_x"), scalar("total_slip_y"), scalar("total_slip_z")};
            const bool sliding = scalar("sliding") == 1.0;
            if (step > 0)
                biaxial_history = biaxial_history && std::abs(scalar("elastic_slip_y")) > 0.0 &&
                                  std::abs(scalar("elastic_slip_z")) > 0.0;
            maximum_coordinate_error =
                std::max({maximum_coordinate_error, std::abs(frame.nodes.at(source)[0] - found->point[0]),
                    std::abs(frame.nodes.at(source)[1] - found->point[1]),
                    std::abs(frame.nodes.at(source)[2] - found->point[2])});
            normal_force.add(-normal[0], found->normal_force[0]);
            tangential_y.add(-tangent[1], found->tangential_force[1]);
            tangential_z.add(-tangent[2], found->tangential_force[2]);
            slip_y.add(slip[1], -found->slip_2);
            slip_z.add(slip[2], found->slip_1);
            opening.add(scalar("gap"), found->opening);
            pressure.add(scalar("pressure"), found->pressure);
            maximum_zero_force = std::max({maximum_zero_force, std::abs(normal[1]), std::abs(normal[2]),
                std::abs(tangent[0]), std::abs(found->normal_force[1]), std::abs(found->normal_force[2]),
                std::abs(found->tangential_force[0])});
            actual_resultant[0] -= normal[0];
            actual_resultant[1] -= tangent[1];
            actual_resultant[2] -= tangent[2];
            all_active = all_active && (scalar("projected") == 1.0) && scalar("pressure") > 0.0;
            first_sticks = first_sticks && (step != 0 || !sliding);
            if (step > 0) {
                const double tangential_magnitude = std::hypot(tangent[0], tangent[1], tangent[2]);
                later_path_is_on_coulomb_circle =
                    later_path_is_on_coulomb_circle &&
                    std::abs(tangential_magnitude - 0.3 * scalar("normal_force")) < 1.0e-9;
            }
            if (sliding)
                ++sliding_count;
            else
                ++sticking_count;
        }

        for (std::size_t component = 0; component < 3; ++component) {
            const double scale = std::max(std::abs(reactions[step].reaction[component]), 1.0);
            maximum_resultant_error = std::max(maximum_resultant_error,
                std::abs(actual_resultant[component] + reactions[step].reaction[component]) / scale);
        }
        std::cout << "b40_step_" << step + 1 << "_sticking_constraints=" << sticking_count << '\n'
                  << "b40_step_" << step + 1 << "_sliding_constraints=" << sliding_count << '\n';
        if (step > 0) later_path_reaches_sliding = later_path_reaches_sliding && sliding_count > 0;
    }
    fuelsim::test::print_relative_metrics("b40_signed_normal_force_x", normal_force);
    fuelsim::test::print_relative_metrics("b40_signed_tangential_force_y", tangential_y);
    fuelsim::test::print_relative_metrics("b40_signed_tangential_force_z", tangential_z);
    fuelsim::test::print_relative_metrics("b40_tangential_slip_y", slip_y);
    fuelsim::test::print_relative_metrics("b40_tangential_slip_z", slip_z);
    fuelsim::test::print_relative_metrics("b40_opening", opening);
    fuelsim::test::print_relative_metrics("b40_pressure", pressure);
    std::cout << "b40_zero_reference_force_count=48\n"
              << "b40_zero_reference_force_maximum_absolute_difference=" << maximum_zero_force << '\n'
              << "b40_maximum_resultant_relative_error=" << maximum_resultant_error << '\n'
              << "b40_maximum_mesh_coordinate_difference=" << maximum_coordinate_error << '\n';
    constexpr double tolerance = 1.0e-2;
    return check(maximum_coordinate_error < 1.0e-14, "B4.0 Abaqus and Fuelsim use the same tracked HEX8 mesh") &&
           check(all_active && first_sticks && later_path_reaches_sliding && later_path_is_on_coulomb_circle,
               "B4.0 keeps every constraint active, starts in sticking, and then reaches the Coulomb circle with both "
               "tangent components") &&
           check(biaxial_history, "B4.0 stores two nonzero committed elastic-slip components") &&
           check(fuelsim::test::relative_metrics_below(normal_force, tolerance) &&
                     fuelsim::test::relative_metrics_below(tangential_y, tolerance) &&
                     fuelsim::test::relative_metrics_below(tangential_z, tolerance),
               "B4.0 signed normal and two tangential nodal-force fields pass all three Abaqus metrics below 1 "
               "percent") &&
           check(fuelsim::test::relative_metrics_below(slip_y, tolerance) &&
                     fuelsim::test::relative_metrics_below(slip_z, tolerance),
               "B4.0 both total tangential-slip fields pass all three Abaqus metrics below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(opening, tolerance) &&
                     fuelsim::test::relative_metrics_below(pressure, tolerance),
               "B4.0 opening and pressure pass all three Abaqus metrics below 1 percent") &&
           check(maximum_zero_force < 1.0e-9,
               "B4.0 theoretical-zero transverse force components pass their absolute check") &&
           check(maximum_resultant_error < tolerance,
               "B4.0 signed normal and two tangential resultants agree with Abaqus below 1 percent");
}

} // namespace

namespace fuelsim::test {
bool check_hex8_biaxial_friction(
    const std::string& output_path, const std::string& contact_path, const std::string& reaction_path) {
    return compare(output_path, contact_path, reaction_path);
}
} // namespace fuelsim::test
