#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Reference final {
    std::size_t node;
    std::array<double, 3> point;
    double gap, pressure;
    std::array<double, 3> normal_force;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete H20.41 reference row: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid H20.41 reference value: " + path);
    return result;
}

std::vector<Reference> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.41 contact reference: " + path);
    std::string line;
    const std::string expected = "id,x,y,z,gap,pressure,normal_x,normal_y,normal_z";
    if (!std::getline(input, line) || line != expected)
        throw std::invalid_argument("Unexpected H20.41 reference header: " + path);
    std::vector<Reference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({static_cast<std::size_t>(number(values, 0, path)),
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}, number(values, 4, path),
            number(values, 5, path), {number(values, 6, path), number(values, 7, path), number(values, 8, path)}});
    }
    if (result.size() != 29) throw std::invalid_argument("H20.41 requires 29 secondary contact nodes");
    return result;
}

bool run(const std::string& output_path, const std::string& reference_path) {
    const std::vector<Reference> references = read_reference(reference_path);
    std::map<std::size_t, Reference> reference_by_node;
    for (const Reference& reference : references)
        if (!reference_by_node.emplace(reference.node, reference).second)
            throw std::invalid_argument("H20.41 contains a repeated reference node");

    const auto output = fuelsim::test::read_final_exodus_results(output_path);
    const auto& projected = output.nodal("contact_projected_interface");
    std::vector<std::size_t> source_nodes;
    for (std::size_t node = 0; node < projected.size(); ++node)
        if (std::isfinite(projected[node])) source_nodes.push_back(node);
    if (source_nodes.size() != references.size())
        throw std::invalid_argument("H20.41 contact output size differs from reference");
    const double penalty = 1e11; // Matched, tracked Abaqus penalty; not a new problem definition.
    bool passed = true;
    fuelsim::test::FieldErrorMetrics constraint_pressure, recovered_pressure;
    fuelsim::test::GroupedFieldErrorMetrics normal_force;
    std::size_t fuelsim_constraint_active = 0, fuelsim_recovered_positive = 0, abaqus_constraint_active = 0,
                abaqus_recovered_positive = 0, abaqus_recovery_extension = 0;
    double coordinate_error = 0.0, abaqus_resultant = 0.0;
    for (std::size_t node = 0; node < source_nodes.size(); ++node) {
        const Reference& reference = reference_by_node.at(source_nodes[node]);
        for (std::size_t component = 0; component < 3; ++component)
            coordinate_error = std::max(
                coordinate_error, std::abs(output.nodes.at(reference.node)[component] - reference.point[component]));
        const auto scalar = [&](const std::string& field) {
            const double value = output.nodal("contact_" + field + "_interface").at(source_nodes[node]);
            if (!std::isfinite(value)) throw std::invalid_argument("H20.41 contact output is not finite");
            return value;
        };
        passed = check(scalar("projected") == 1.0, "H20.41 projects every constraint") && passed;
        const double actual_constraint = scalar("constraint_pressure");
        const double reference_constraint = std::max(-penalty * reference.gap, 0.0);
        constraint_pressure.add(actual_constraint, reference_constraint);
        recovered_pressure.add(scalar("pressure"), reference.pressure);
        std::array<double, 3> actual_normal_force{};
        for (std::size_t component = 0; component < actual_normal_force.size(); ++component)
            actual_normal_force[component] = -scalar("normal_force_" + std::string(1, "xyz"[component]));
        normal_force.add(actual_normal_force.data(), reference.normal_force.data(), actual_normal_force.size());
        fuelsim_constraint_active += actual_constraint > 0.0 ? 1U : 0U;
        fuelsim_recovered_positive += scalar("pressure") > 0.0 ? 1U : 0U;
        abaqus_constraint_active += reference_constraint > 0.0 ? 1U : 0U;
        abaqus_recovered_positive += reference.pressure > 0.0 ? 1U : 0U;
        abaqus_recovery_extension += reference.gap > 0.0 && reference.pressure > 0.0 ? 1U : 0U;
        abaqus_resultant += reference.normal_force[0];
    }

    const auto side = std::find(output.side_set_names.begin(), output.side_set_names.end(), "secondary_contact");
    if (side == output.side_set_names.end()) throw std::invalid_argument("H20.41 output lacks secondary contact sides");
    const auto& faces = output.side_set_face_nodes.at(static_cast<std::size_t>(side - output.side_set_names.begin()));
    std::size_t partial_faces = 0;
    for (const auto& face : faces) {
        std::size_t active = 0;
        if (face.size() != 8) throw std::invalid_argument("H20.41 contact faces must contain eight nodes");
        for (const std::size_t source : face) active += reference_by_node.at(source).gap < 0.0 ? 1U : 0U;
        partial_faces += active > 0 && active < face.size() ? 1U : 0U;
    }

    const double total_contact_force = output.global("contact_force_interface");
    const double resultant_error =
        std::abs(total_contact_force - std::abs(abaqus_resultant)) / std::abs(abaqus_resultant);
    fuelsim::test::print_relative_metrics("h20_41_constraint_pressure", constraint_pressure);
    fuelsim::test::print_relative_metrics("h20_41_recovered_pressure", recovered_pressure);
    fuelsim::test::print_grouped_relative_metrics("h20_41_normal_force", normal_force);
    std::cout << "h20_41_fuelsim_constraint_active_nodes=" << fuelsim_constraint_active << '\n'
              << "h20_41_fuelsim_recovered_positive_nodes=" << fuelsim_recovered_positive << '\n'
              << "h20_41_abaqus_constraint_active_nodes=" << abaqus_constraint_active << '\n'
              << "h20_41_abaqus_recovered_positive_nodes=" << abaqus_recovered_positive << '\n'
              << "h20_41_abaqus_recovery_extension_nodes=" << abaqus_recovery_extension << '\n'
              << "h20_41_partial_faces=" << partial_faces << '\n'
              << "h20_41_fuelsim_resultant=" << total_contact_force << '\n'
              << "h20_41_abaqus_resultant=" << std::abs(abaqus_resultant) << '\n'
              << "h20_41_resultant_relative_error=" << resultant_error << '\n';

    constexpr double tolerance = 1.0e-2;
    passed = check(coordinate_error < 1.0e-14, "H20.41 uses identical Abaqus and Fuelsim coordinates") && passed;
    passed =
        check(fuelsim_constraint_active == 11 && fuelsim_recovered_positive == 16 && abaqus_constraint_active == 11 &&
                  abaqus_recovered_positive == 16 && abaqus_recovery_extension == 5 && partial_faces == 4,
            "H20.41 recovers the Abaqus positive-pressure support across partially active quadratic faces") &&
        passed;
    passed = check(fuelsim::test::relative_metrics_below(constraint_pressure, tolerance),
                 "H20.41 constraint-pressure metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(recovered_pressure, tolerance),
                 "H20.41 recovered-pressure metrics are below one percent") &&
             passed;
    passed = check(fuelsim::test::grouped_relative_metrics_below(normal_force, tolerance),
                 "H20.41 nodal normal-force metrics are below one percent") &&
             passed;
    passed = check(recovered_pressure.maximum_zero_reference_difference == 0.0,
                 "H20.41 reports zero pressure at every Abaqus zero-pressure node") &&
             passed;
    passed =
        check(resultant_error < tolerance, "H20.41 total normal force differs from Abaqus by less than one percent") &&
        passed;
    return passed;
}
} // namespace

namespace fuelsim::test {
bool check_hex20_partial_contact(const std::string& output_path, const std::string& reference_path) {
    return run(output_path, reference_path);
}
} // namespace fuelsim::test
