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
using Row = std::map<std::string, double>;

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) {
        if (!field.empty() && field.back() == '\r') field.pop_back();
        result.push_back(field);
    }
    return result;
}

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("Cannot read Abaqus reference " + path);
    std::string line;
    std::getline(stream, line);
    const auto header = split(line);
    std::vector<Row> rows;
    while (std::getline(stream, line)) {
        if (line.empty()) continue;
        const auto fields = split(line);
        if (fields.size() != header.size()) throw std::runtime_error("Invalid CSV row in " + path);
        Row row;
        for (std::size_t i = 0; i < fields.size(); ++i) {
            std::size_t consumed = 0;
            const double value = std::stod(fields[i], &consumed);
            if (consumed != fields[i].size() || !std::isfinite(value))
                throw std::runtime_error("Invalid Abaqus reference number in " + path);
            row.emplace(header[i], value);
        }
        rows.push_back(std::move(row));
    }
    if (rows.empty()) throw std::runtime_error("Empty Abaqus reference " + path);
    return rows;
}

bool check_scalar(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics, double zero_tolerance) {
    if (metrics.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metrics);
    else {
        fuelsim::test::print_absolute_metrics(name, metrics);
        std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
                  << '\n';
    }
    return metrics.value_count > 0 &&
           (!metrics.has_relative_norm() || fuelsim::test::relative_metrics_below(metrics, 0.005)) &&
           metrics.maximum_zero_reference_difference <= zero_tolerance;
}

bool check_group(
    const std::string& name, const fuelsim::test::GroupedFieldErrorMetrics& metrics, double zero_tolerance) {
    if (metrics.has_relative_norm())
        fuelsim::test::print_grouped_relative_metrics(name, metrics);
    else
        std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
                  << '\n';
    return metrics.group_count > 0 &&
           (!metrics.has_relative_norm() || fuelsim::test::grouped_relative_metrics_below(metrics, 0.005)) &&
           metrics.maximum_zero_reference_difference <= zero_tolerance;
}

bool same_time(double reference, double actual) {
    // ODB times are float32; keep the tolerance relative even for the short preload.
    return std::abs(reference - actual) <= 8e-8 * std::max(std::abs(reference), std::abs(actual)) + 1e-16;
}
} // namespace

namespace fuelsim::test {
bool check_rz_abaqus(const std::string& output_path, const std::string& node_path, const std::string& point_path,
    const std::string& mechanisms) {
    const auto nodes = read_rows(node_path), points = read_rows(point_path);
    const auto frames = read_exodus_history(output_path);
    const bool contact = mechanisms == "contact";
    if (!contact && mechanisms != "plastic" && mechanisms != "creep" && mechanisms != "coupled" &&
        mechanisms != "thermal")
        throw std::runtime_error("Unknown RZ qualification mechanism");
    FieldErrorMetrics temperature, reaction, plastic, creep;
    FieldErrorMetrics pressure, gap, nodal_normal_force, normal_force, heat_rate;
    GroupedFieldErrorMetrics displacement;
    std::array<GroupedFieldErrorMetrics, 4> tensors;
    std::size_t node_row = 0, point_row = 0, accepted = 0;
    double max_plastic = 0.0, max_creep = 0.0;
    constexpr std::array<std::size_t, 4> point_map = {0, 1, 3, 2};
    const std::array<std::string, 4> prefixes = {"stress_", "elastic_", "plastic_", "creep_"};
    const std::array<std::string, 4> components = {"rr", "zz", "hoop", "rz"};
    for (const auto& frame : frames) {
        if (frame.time == 0.0) continue;
        ++accepted;
        double actual_reaction = 0.0, reference_reaction = 0.0;
        double reference_force = 0.0, reference_heat = 0.0, reference_opposite_force = 0.0;
        double bottom_z = frame.nodes.front()[1];
        for (const auto& node : frame.nodes) bottom_z = std::min(bottom_z, node[1]);
        for (std::size_t n = 0; n < frame.nodes.size(); ++n) {
            const auto& row = nodes.at(node_row++);
            // Abaqus stores frame times as single precision even with full field precision.
            if (!same_time(row.at("time"), frame.time) || row.at("node") != static_cast<double>(n + 1))
                throw std::runtime_error("Abaqus nodal history does not match production output");
            temperature.add(frame.nodal("temperature")[n], row.at("temperature"));
            const std::array<double, 2> actual = {frame.nodal("displacement_r")[n], frame.nodal("displacement_z")[n]};
            std::array<double, 2> reference = {row.at("ur"), row.at("uz")};
            // The axis/bottom intersection is analytically fixed in both components.
            // Audit its raw reference first, then compare to the exact prescribed zero.
            if ((frame.nodes[n][0] == 0.0 && (frame.nodes[n][1] == 0.0 || mechanisms == "thermal")) ||
                (contact && frame.nodes[n][1] == bottom_z)) {
                if (std::hypot(reference[0], reference[1]) > 1e-12)
                    throw std::runtime_error("Abaqus violates the fixed origin displacement");
                reference = {0.0, 0.0};
            }
            displacement.add(actual.data(), reference.data(), actual.size());
            if (frame.nodes[n][1] == bottom_z) {
                if (!contact) actual_reaction += frame.nodal("reaction_force_z")[n];
                reference_reaction += row.at("rf_z");
                if (contact) reference_heat += row.at("reaction_heat");
            }
            if (contact) {
                reference_opposite_force += row.at("contact_force_z");
                if (n == 2 || n == 3) {
                    pressure.add(frame.nodal("contact_pressure_interface")[n], row.at("contact_pressure"));
                    gap.add(frame.nodal("contact_gap_interface")[n], row.at("contact_gap"));
                    nodal_normal_force.add(
                        frame.nodal("contact_normal_force_interface")[n], -row.at("contact_force_z"));
                    reference_force -= row.at("contact_force_z");
                    if (row.at("contact_pressure") <= 1e6 || row.at("contact_gap") >= 0.0 ||
                        row.at("contact_heat_flux") <= 1e5 || frame.nodal("contact_projected_interface")[n] != 1.0)
                        throw std::runtime_error("RZ thermal and mechanical contact must both be active");
                }
            }
        }
        if (contact) {
            normal_force.add(frame.global("contact_force_interface"), reference_force);
            heat_rate.add(frame.global("contact_heat_rate_interface"), reference_heat);
            if (std::abs(reference_opposite_force) > 1e-7 || std::abs(reference_force - reference_reaction) > 1e-7)
                throw std::runtime_error("Abaqus contact reference violates force equilibrium");
        } else
            reaction.add(actual_reaction, reference_reaction);
        std::size_t elements = 0;
        for (const auto count : frame.block_element_counts) elements += count;
        for (std::size_t e = 0; e < elements; ++e) {
            for (std::size_t aq = 0; aq < 4; ++aq) {
                const auto& row = points.at(point_row++);
                if (!same_time(row.at("time"), frame.time) || row.at("element") != static_cast<double>(e + 1) ||
                    row.at("point") != static_cast<double>(aq + 1))
                    throw std::runtime_error("Abaqus material history does not match production output");
                const auto suffix = "_q" + std::to_string(point_map[aq]);
                for (std::size_t t = 0; t < (contact ? 1 : tensors.size()); ++t) {
                    std::array<double, 4> actual{}, reference{};
                    for (std::size_t c = 0; c < components.size(); ++c) {
                        const auto field = prefixes[t] + components[c];
                        const double weight = c == 3 ? std::sqrt(2.0) : 1.0;
                        actual[c] = weight * frame.element(field + suffix)[e];
                        reference[c] = weight * row.at(field);
                    }
                    tensors[t].add(actual.data(), reference.data(), actual.size());
                }
                if (!contact) {
                    plastic.add(frame.element("equiv_plastic" + suffix)[e], row.at("equiv_plastic"));
                    creep.add(frame.element("equiv_creep" + suffix)[e], row.at("equiv_creep"));
                }
                max_plastic = std::max(max_plastic, row.at("equiv_plastic"));
                max_creep = std::max(max_creep, row.at("equiv_creep"));
            }
        }
    }
    bool passed = accepted > 0 && node_row == nodes.size() && point_row == points.size();
    passed = check_scalar("rz_temperature", temperature, 1e-10) && passed;
    passed = check_group("rz_displacement", displacement, 1e-12) && passed;
    if (!contact) passed = check_scalar("rz_bottom_axial_reaction", reaction, 1e-7) && passed;
    for (std::size_t t = 0; t < (contact ? 1 : tensors.size()); ++t)
        passed = check_group("rz_" + prefixes[t] + "tensor", tensors[t], t == 0 ? 1e-3 : 1e-12) && passed;
    if (contact) {
        passed = check_scalar("rz_contact_pressure", pressure, 1e-3) && passed;
        passed = check_scalar("rz_contact_gap", gap, 1e-12) && passed;
        passed = check_scalar("rz_contact_nodal_normal_force", nodal_normal_force, 1e-7) && passed;
        passed = check_scalar("rz_contact_normal_force", normal_force, 1e-7) && passed;
        passed = check_scalar("rz_contact_heat_rate", heat_rate, 1e-10) && passed;
    } else {
        passed = check_scalar("rz_equivalent_plastic_strain", plastic, 1e-12) && passed;
        passed = check_scalar("rz_equivalent_creep_strain", creep, 1e-12) && passed;
    }
    if (mechanisms == "plastic" || mechanisms == "coupled") passed = max_plastic > 1e-5 && passed;
    if (mechanisms == "creep" || mechanisms == "coupled") passed = max_creep > 1e-6 && passed;
    std::cout << "rz_accepted_frames=" << accepted << "\nrz_abaqus_qualification=" << (passed ? "passed" : "failed")
              << '\n';
    return passed;
}
} // namespace fuelsim::test
