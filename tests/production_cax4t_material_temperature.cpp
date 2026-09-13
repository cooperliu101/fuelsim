#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Row = std::map<std::string, std::string>;
constexpr double relative_tolerance = 0.001;
constexpr double displacement_zero_tolerance = 1e-13;
constexpr double strain_zero_tolerance = 1e-13;
constexpr double stress_zero_tolerance = 1e-4;

std::vector<std::string> split(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::string value;
    while (std::getline(stream, value, ',')) {
        if (!value.empty() && value.back() == '\r')
            value.pop_back();
        result.push_back(value);
    }
    return result;
}

double number(const Row& row, const std::string& name) {
    const std::string& text = row.at(name);
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value))
        throw std::runtime_error("Invalid native field " + name);
    return value;
}

std::vector<Row> read_accepted_rows(const std::string& path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("Cannot read native reference " + path);
    std::string line;
    std::getline(stream, line);
    const auto header = split(line);
    std::vector<Row> result;
    while (std::getline(stream, line)) {
        if (line.empty())
            continue;
        const auto fields = split(line);
        if (fields.size() != header.size())
            throw std::runtime_error("Invalid native CSV width " + path);
        Row row;
        for (std::size_t i = 0; i < fields.size(); ++i)
            if (!row.emplace(header[i], fields[i]).second)
                throw std::runtime_error("Duplicate native CSV column");
        if (number(row, "frame") != 0)
            result.push_back(std::move(row));
    }
    return result;
}

void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

bool same_time(double a, double b) {
    // Native ODB frame times have binary32 precision; material values are binary64.
    return std::abs(a - b) <= 1.2e-7;
}

bool check(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics, double zero_tolerance) {
    using namespace fuelsim::test;
    if (metrics.has_relative_norm())
        print_relative_metrics(name, metrics);
    else
        print_absolute_metrics(name, metrics);
    return metrics.value_count > 0
           && (!metrics.has_relative_norm() || relative_metrics_below(metrics, relative_tolerance))
           && metrics.maximum_zero_reference_difference <= zero_tolerance;
}
} // namespace

namespace fuelsim::test {
bool check_cax4t_material_temperature(const std::string& output,
    const std::string& references,
    const std::string& case_name) {
    const bool elastic = case_name == "elastic_moduli" || case_name == "elastic_moduli_finite";
    const bool plastic = case_name == "plastic_yield" || case_name == "plastic_yield_finite"
                         || case_name == "plastic_hardening" || case_name == "plastic_hardening_finite";
    const bool creep = case_name == "creep_coefficient_finite" || case_name == "creep_exponent_finite";
    require(elastic || plastic || creep, "Unexpected CAX4T material temperature case");
    const std::size_t increments = elastic ? 8 : 10;
    const double dt = elastic ? 0.25 : plastic ? 0.2 : 0.1;
    const std::string prefix = references + "/" + (creep ? "production_" : "") + case_name;
    const auto nodes = read_accepted_rows(prefix + "_nodes.csv");
    const auto points = read_accepted_rows(prefix + "_points.csv");
    const auto frames = read_exodus_history(output);
    require(frames.size() == increments + 1 && nodes.size() == increments * 4 && points.size() == increments * 4,
        "Incomplete native or production increment coverage");
    require(same_time(frames.front().time, 0), "Production initial time differs");
    constexpr std::array<std::size_t, 4> point_map = {0, 1, 3, 2};
    const std::array<std::string, 4> components = {"rr", "zz", "hoop", "rz"};
    const std::array<std::string, 4> fields = {"stress_", "elastic_", "plastic_", "creep_"};
    const std::array<std::string, 4> native = {"s_", "ee_", "pe_", "ce_"};
    const std::array<std::string, 2> equivalent = {"equiv_plastic", "equiv_creep"};
    const std::array<std::string, 2> native_equivalent = {"peeq", "ceeq"};
    FieldErrorMetrics temperature, radial, axial, fixed_displacement, zero_stress, zero_strain;
    std::array<std::array<FieldErrorMetrics, 3>, 4> normal;
    std::array<FieldErrorMetrics, 2> scalars, scalar_increments;
    std::array<std::array<double, 4>, 2> previous_actual{}, previous_native{};
    for (std::size_t step = 0; step <= increments; ++step) {
        const auto& frame = frames[step];
        require(frame.nodes.size() == 4 && frame.block_element_counts == std::vector<std::size_t>{1}
                    && frame.element("material_point_count").at(0) == 4,
            "Expected one complete CAX4T element and four nodes");
        require(same_time(frame.time, dt * static_cast<double>(step)), "Production increment times differ");
        for (std::size_t q = 0; q < 4; ++q) {
            const std::string suffix = "_q" + std::to_string(point_map[q]);
            if (step == 0) {
                for (std::size_t t = 0; t < fields.size(); ++t)
                    for (const auto& component : components)
                        (t == 0 ? zero_stress : zero_strain)
                            .add(frame.element(fields[t] + component + suffix).at(0), 0);
                for (const auto& field : equivalent)
                    zero_strain.add(frame.element(field + suffix).at(0), 0);
                continue;
            }
            const Row& point = points[(step - 1) * 4 + q];
            const std::string expected_step = elastic   ? (step <= 4 ? "FIRST" : "SWAP")
                                              : plastic ? (step <= 5 ? "LOAD_FIXED_T" : "LOAD_SWAPPED_T")
                                                        : "LOAD_AND_HOLD";
            const std::size_t native_frame = elastic ? (step - 1) % 4 + 1 : plastic ? (step - 1) % 5 + 1 : step;
            require(number(point, "element") == 1 && number(point, "point") == static_cast<double>(q + 1)
                        && point.at("step") == expected_step
                        && number(point, "frame") == static_cast<double>(native_frame)
                        && same_time(number(point, "time"), frame.time),
                "Native material point order or time differs");
            const Row& node = nodes[(step - 1) * 4 + q];
            require(number(node, "node") == static_cast<double>(q + 1) && node.at("step") == expected_step
                        && number(node, "frame") == static_cast<double>(native_frame)
                        && same_time(number(node, "time"), frame.time)
                        && std::abs(number(node, "r") - frame.nodes[q][0]) < 1e-15
                        && std::abs(number(node, "z") - frame.nodes[q][1]) < 1e-15,
                "Native node identity, coordinates or time differs");
            temperature.add(frame.nodal("temperature").at(q), number(node, "temperature"));
            if (elastic)
                radial.add(frame.nodal("displacement_r").at(q), number(node, "ur"));
            else {
                fixed_displacement.add(frame.nodal("displacement_r").at(q), 0);
                fixed_displacement.add(number(node, "ur"), 0);
            }
            if (q >= 2)
                axial.add(frame.nodal("displacement_z").at(q), number(node, "uz"));
            else {
                fixed_displacement.add(frame.nodal("displacement_z").at(q), 0);
                fixed_displacement.add(number(node, "uz"), 0);
            }
            for (std::size_t t = 0; t < fields.size(); ++t) {
                const bool active = t < 2 || (t == 2 && plastic) || (t == 3 && creep);
                for (std::size_t c = 0; c < components.size(); ++c) {
                    const double actual = frame.element(fields[t] + components[c] + suffix).at(0);
                    if (active && c < 3)
                        normal[t][c].add(actual, number(point, native[t] + components[c]));
                    else {
                        // These diagonal loading paths have analytically zero shear.
                        // Native roundoff shear is checked absolutely, never used as a relative denominator.
                        (t == 0 ? zero_stress : zero_strain).add(actual, 0);
                        if (active)
                            (t == 0 ? zero_stress : zero_strain).add(number(point, native[t] + components[c]), 0);
                    }
                }
            }
            for (std::size_t t = 0; t < equivalent.size(); ++t) {
                const double actual = frame.element(equivalent[t] + suffix).at(0);
                if ((t == 0 && plastic) || (t == 1 && creep)) {
                    const double reference = number(point, native_equivalent[t]);
                    const double actual_increment = actual - previous_actual[t][q];
                    const double native_increment = reference - previous_native[t][q];
                    require(actual_increment > 1e-10 && native_increment > 1e-10,
                        "Expected an active inelastic increment at every material point");
                    scalars[t].add(actual, reference);
                    scalar_increments[t].add(actual_increment, native_increment);
                    previous_actual[t][q] = actual;
                    previous_native[t][q] = reference;
                } else
                    zero_strain.add(actual, 0);
            }
        }
    }
    bool passed = check("cax4t_temperature", temperature, 1e-10)
                  && check("cax4t_axial_displacement", axial, displacement_zero_tolerance)
                  && check("cax4t_fixed_displacement", fixed_displacement, displacement_zero_tolerance)
                  && check("cax4t_zero_stress", zero_stress, stress_zero_tolerance)
                  && check("cax4t_zero_strain", zero_strain, strain_zero_tolerance);
    if (elastic)
        passed = check("cax4t_radial_displacement", radial, displacement_zero_tolerance) && passed;
    for (std::size_t t = 0; t < fields.size(); ++t)
        if (t < 2 || (t == 2 && plastic) || (t == 3 && creep))
            for (std::size_t c = 0; c < 3; ++c)
                passed = check("cax4t_" + fields[t] + components[c],
                             normal[t][c],
                             t == 0 ? stress_zero_tolerance : strain_zero_tolerance)
                         && passed;
    for (std::size_t t = 0; t < equivalent.size(); ++t)
        if ((t == 0 && plastic) || (t == 1 && creep)) {
            passed = check("cax4t_" + equivalent[t], scalars[t], strain_zero_tolerance) && passed;
            passed =
                check("cax4t_" + equivalent[t] + "_increment", scalar_increments[t], strain_zero_tolerance) && passed;
        }
    std::cout << "cax4t_native_node_samples=" << nodes.size() << '\n'
              << "cax4t_native_material_samples=" << points.size() << '\n'
              << "cax4t_material_temperature_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace fuelsim::test
