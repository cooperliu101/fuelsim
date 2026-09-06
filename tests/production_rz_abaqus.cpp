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
        mechanisms != "thermal" && mechanisms != "sliding" && mechanisms != "probe")
        throw std::runtime_error("Unknown RZ qualification mechanism");
    FieldErrorMetrics temperature, reaction, plastic, creep;
    FieldErrorMetrics pressure, gap, nodal_normal_force, normal_force, heat_rate;
    GroupedFieldErrorMetrics displacement, support_reaction;
    FieldErrorMetrics free_reaction;
    FieldErrorMetrics reaction_heat;
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
        double outer_r = frame.nodes.front()[0];
        for (const auto& node : frame.nodes) bottom_z = std::min(bottom_z, node[1]);
        for (const auto& node : frame.nodes) outer_r = std::max(outer_r, node[0]);
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
                (contact && frame.nodes[n][1] == bottom_z) ||
                (mechanisms == "sliding" && frame.nodes[n][0] == outer_r)) {
                if (std::hypot(reference[0], reference[1]) > 1e-12)
                    throw std::runtime_error("Abaqus violates the fixed origin displacement");
                reference = {0.0, 0.0};
            }
            displacement.add(actual.data(), reference.data(), actual.size());
            if (mechanisms == "probe") {
                const std::array<double, 2> force = {
                    frame.nodal("reaction_force_r")[n], frame.nodal("reaction_force_z")[n]};
                const std::array<double, 2> reference_force_vector = {row.at("rf_r"), row.at("rf_z")};
                support_reaction.add(force.data(), reference_force_vector.data(), 2);
                reaction_heat.add(frame.nodal("reaction_heat_flux")[n], row.at("reaction_heat"));
            }
            if (mechanisms == "sliding") {
                const std::array<double, 2> force = {
                    frame.nodal("reaction_force_r")[n], frame.nodal("reaction_force_z")[n]};
                const std::array<double, 2> reference_force_vector = {row.at("rf_r"), row.at("rf_z")};
                if ((n < 6 && n % 2 == 0) || (n >= 6 && n % 2 == 1))
                    support_reaction.add(force.data(), reference_force_vector.data(), force.size());
                else {
                    if (std::hypot(reference_force_vector[0], reference_force_vector[1]) > 1e-7)
                        throw std::runtime_error("Abaqus free node has a nonzero boundary reaction");
                    for (double component : force) free_reaction.add(component, 0.0);
                }
            }
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
            const double reported_count = frame.element("material_point_count")[e];
            if (reported_count != 1.0 && reported_count != 4.0 && reported_count != 9.0)
                throw std::runtime_error("Invalid RZ material point count");
            const auto point_count = static_cast<std::size_t>(reported_count);
            if (point_count == 1)
                for (std::size_t q = 1; q < 4; ++q)
                    for (const auto& component : components)
                        if (!std::isnan(frame.element("stress_" + component + "_q" + std::to_string(q))[e]))
                            throw std::runtime_error("Inactive CAX4RT output must not masquerade as a material point");
            for (std::size_t aq = 0; aq < point_count; ++aq) {
                const auto& row = points.at(point_row++);
                if (!same_time(row.at("time"), frame.time) || row.at("element") != static_cast<double>(e + 1) ||
                    row.at("point") != static_cast<double>(aq + 1))
                    throw std::runtime_error("Abaqus material history does not match production output");
                const auto suffix = "_q" + std::to_string(point_count == 9 ? aq : point_map[aq]);
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
    if (mechanisms == "sliding" || mechanisms == "probe") {
        passed = check_group("rz_support_reaction", support_reaction, 1e-7) && passed;
        if (mechanisms == "sliding") passed = check_scalar("rz_free_node_reaction", free_reaction, 1e-7) && passed;
        if (mechanisms == "probe") passed = check_scalar("rz_reaction_heat", reaction_heat, 1e-7) && passed;
    }
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

bool check_rz_sliding_abaqus(const std::string& output_path, const std::string& node_path,
    const std::string& point_path, const std::string& contact_path, bool require_segment_crossing) {
    bool passed = check_rz_abaqus(output_path, node_path, point_path, "sliding");
    const auto rows = read_rows(contact_path);
    const auto frames = read_exodus_history(output_path);
    FieldErrorMetrics pressure, gap, shear, slip, normal_force, tangent_force;
    GroupedFieldErrorMetrics normal_vector, tangent_vector;
    std::map<std::size_t, std::vector<double>> ownership;
    std::map<std::size_t, bool> previous_sliding;
    std::map<std::size_t, double> previous_slip, previous_shear;
    std::size_t row_index = 0, sticking = 0, sliding = 0;
    bool forward = false, reverse = false, reverse_sticking = false;
    for (const auto& frame : frames) {
        if (frame.time == 0.0) continue;
        double reference_normal = 0.0, reference_tangent = 0.0;
        std::array<double, 2> actual_normal_sum{}, reference_normal_sum{}, actual_tangent_sum{},
            reference_tangent_sum{};
        for (const std::size_t label : {2U, 4U, 6U}) {
            const auto& row = rows.at(row_index++);
            if (row.at("node") != static_cast<double>(label) || !same_time(row.at("time"), frame.time))
                throw std::runtime_error("RZ sliding contact reference does not match output");
            const auto n = label - 1;
            pressure.add(frame.nodal("contact_pressure_interface")[n], row.at("pressure"));
            gap.add(frame.nodal("contact_gap_interface")[n], row.at("gap"));
            shear.add(frame.nodal("contact_tangential_traction_interface")[n], row.at("shear"));
            slip.add(frame.nodal("contact_total_tangential_slip_interface")[n], row.at("slip"));
            normal_force.add(
                frame.nodal("contact_normal_force_interface")[n], std::hypot(row.at("normal_r"), row.at("normal_z")));
            tangent_force.add(frame.nodal("contact_tangential_force_interface")[n],
                std::copysign(std::hypot(row.at("tangential_r"), row.at("tangential_z")), row.at("shear")));
            reference_normal += std::hypot(row.at("normal_r"), row.at("normal_z"));
            reference_tangent +=
                std::copysign(std::hypot(row.at("tangential_r"), row.at("tangential_z")), row.at("shear"));
            if (frame.nodal("contact_projected_interface")[n] != 1.0)
                throw std::runtime_error("Sliding path loses its contact projection");
            ownership[label].push_back(frame.nodal("contact_primary_segment_interface")[n]);
            const double segment_value = frame.nodal("contact_primary_segment_interface")[n];
            if (segment_value < 0.0 || segment_value > 5.0 || std::floor(segment_value) != segment_value)
                throw std::runtime_error("Invalid primary segment in coaxial-ring output");
            const auto a = 6 + 2 * static_cast<std::size_t>(segment_value), b = a + 2;
            const double dr = frame.nodes[b][0] + frame.nodal("displacement_r")[b] - frame.nodes[a][0] -
                              frame.nodal("displacement_r")[a];
            const double dz = frame.nodes[b][1] + frame.nodal("displacement_z")[b] - frame.nodes[a][1] -
                              frame.nodal("displacement_z")[a];
            const double length = std::hypot(dr, dz);
            const double normal = frame.nodal("contact_normal_force_interface")[n];
            const double tangent = frame.nodal("contact_tangential_force_interface")[n];
            const std::array<double, 2> actual_n = {-normal * dz / length, normal * dr / length};
            const std::array<double, 2> actual_t = {-tangent * dr / length, -tangent * dz / length};
            const std::array<double, 2> reference_n = {row.at("normal_r"), row.at("normal_z")};
            const std::array<double, 2> reference_t = {row.at("tangential_r"), row.at("tangential_z")};
            normal_vector.add(actual_n.data(), reference_n.data(), 2);
            tangent_vector.add(actual_t.data(), reference_t.data(), 2);
            for (std::size_t c = 0; c < 2; ++c) {
                actual_normal_sum[c] += actual_n[c];
                reference_normal_sum[c] += reference_n[c];
                actual_tangent_sum[c] += actual_t[c];
                reference_tangent_sum[c] += reference_t[c];
            }
            const bool is_sliding = frame.nodal("contact_sliding_interface")[n] == 1.0;
            if (!(row.at("pressure") > 1e5) || !(row.at("gap") < 0.0))
                throw std::runtime_error("Friction qualification requires active compressive contact");
            const double friction_ratio = std::abs(row.at("shear")) / (0.2 * row.at("pressure"));
            const bool reference_sliding = std::abs(friction_ratio - 1.0) < 1e-7;
            if (friction_ratio > 1.0 + 1e-7 || reference_sliding != is_sliding)
                throw std::runtime_error("Coulomb sticking/sliding state differs from Abaqus");
            if (previous_sliding[label] && !reference_sliding &&
                (row.at("slip") - previous_slip[label]) * previous_shear[label] < 0.0)
                reverse_sticking = true;
            previous_sliding[label] = reference_sliding;
            previous_slip[label] = row.at("slip");
            previous_shear[label] = row.at("shear");
            sliding += is_sliding ? 1 : 0;
            sticking += is_sliding ? 0 : 1;
            forward = forward || (is_sliding && row.at("shear") > 0.0);
            reverse = reverse || (is_sliding && row.at("shear") < 0.0);
        }
        normal_force.add(frame.global("contact_force_interface"), reference_normal);
        tangent_force.add(frame.global("contact_tangential_force_interface"), reference_tangent);
        normal_vector.add(actual_normal_sum.data(), reference_normal_sum.data(), 2);
        tangent_vector.add(actual_tangent_sum.data(), reference_tangent_sum.data(), 2);
    }
    for (const auto& field :
        std::vector<std::pair<std::string, const FieldErrorMetrics*>>{{"pressure", &pressure}, {"gap", &gap},
            {"shear", &shear}, {"slip", &slip}, {"normal_force", &normal_force}, {"tangent_force", &tangent_force}})
        passed = check_scalar("rz_sliding_" + field.first, *field.second,
                     field.first == "gap" || field.first == "slip" ? 1e-12 : 1e-7) &&
                 passed;
    passed = check_group("rz_sliding_normal_force_vector", normal_vector, 1e-7) && passed;
    passed = check_group("rz_sliding_tangent_force_vector", tangent_vector, 1e-7) && passed;
    std::size_t maximum_crossed = 0;
    for (const auto& entry : ownership) {
        const auto limits = std::minmax_element(entry.second.begin(), entry.second.end());
        maximum_crossed = std::max(maximum_crossed, static_cast<std::size_t>(*limits.second - *limits.first));
    }
    std::cout << "rz_sliding_maximum_segments_crossed=" << maximum_crossed << '\n'
              << "rz_sliding_sticking_samples=" << sticking << '\n'
              << "rz_sliding_sliding_samples=" << sliding << '\n'
              << "rz_sliding_reverse_sticking=" << reverse_sticking << '\n';
    passed = passed && row_index == rows.size() && (!require_segment_crossing || maximum_crossed >= 2) &&
             sticking > 0 && forward && reverse && reverse_sticking;
    std::cout << "rz_friction_abaqus_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace fuelsim::test
