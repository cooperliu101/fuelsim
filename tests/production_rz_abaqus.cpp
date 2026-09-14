#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Row = std::map<std::string, double>;

constexpr double rz_relative_tolerance = 0.001;
constexpr double rz_temperature_zero_tolerance = 1e-11;
constexpr double rz_displacement_zero_tolerance = 1e-13;
constexpr double rz_force_zero_tolerance = 1e-8;
constexpr double rz_stress_zero_tolerance = 1e-4;
constexpr double rz_strain_zero_tolerance = 1e-13;
constexpr double rz_heat_rate_zero_tolerance = 1e-11;

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) {
        if (!field.empty() && field.back() == '\r')
            field.pop_back();
        result.push_back(field);
    }
    return result;
}

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("Cannot read Abaqus reference " + path);
    std::string line;
    std::getline(stream, line);
    const auto header = split(line);
    std::vector<Row> rows;
    while (std::getline(stream, line)) {
        if (line.empty())
            continue;
        const auto fields = split(line);
        if (fields.size() != header.size())
            throw std::runtime_error("Invalid CSV row in " + path);
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
    if (rows.empty())
        throw std::runtime_error("Empty Abaqus reference " + path);
    return rows;
}

bool check_scalar(const std::string& name,
    const fuelsim::test::FieldErrorMetrics& metrics,
    double zero_tolerance,
    bool absolute_only = false,
    double allowed_relative = rz_relative_tolerance) {
    if (absolute_only) {
        std::cout << name << "_maximum_absolute_difference=" << metrics.maximum_absolute_difference << '\n';
        return metrics.value_count > 0 && metrics.maximum_absolute_difference <= zero_tolerance;
    }
    if (metrics.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metrics);
    else {
        fuelsim::test::print_absolute_metrics(name, metrics);
        std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
                  << '\n';
    }
    return metrics.value_count > 0
           && (!metrics.has_relative_norm() || fuelsim::test::relative_metrics_below(metrics, allowed_relative))
           && metrics.maximum_zero_reference_difference <= zero_tolerance;
}

bool check_group(const std::string& name,
    const fuelsim::test::GroupedFieldErrorMetrics& metrics,
    double zero_tolerance,
    bool absolute_only = false) {
    if (absolute_only) {
        std::cout << name << "_maximum_absolute_difference=" << metrics.maximum_difference << '\n';
        return metrics.group_count > 0 && metrics.maximum_difference <= zero_tolerance;
    }
    if (metrics.has_relative_norm())
        fuelsim::test::print_grouped_relative_metrics(name, metrics);
    else
        std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
                  << '\n';
    return metrics.group_count > 0
           && (!metrics.has_relative_norm()
               || fuelsim::test::grouped_relative_metrics_below(metrics, rz_relative_tolerance))
           && metrics.maximum_zero_reference_difference <= zero_tolerance;
}

bool same_time(double reference, double actual) {
    // ODB times are float32; keep the tolerance relative even for the short preload.
    return std::abs(reference - actual) <= 8e-8 * std::max(std::abs(reference), std::abs(actual)) + 1e-16;
}

bool check_creep_integration(const std::string& output, const std::string& node_path, const std::string& point_path) {
    using namespace fuelsim::test;
    const auto frames = read_exodus_history(output);
    const auto nodes = read_rows(node_path), points = read_rows(point_path);
    if (frames.size() != 11)
        throw std::runtime_error("Creep ramp requires ten fixed increments");
    std::array<FieldErrorMetrics, 2> temperature, reaction, equivalent, plastic, heat, creep_increment,
        plastic_increment;
    std::array<GroupedFieldErrorMetrics, 2> displacement, stress, elastic, creep_tensor, plastic_tensor;
    FieldErrorMetrics cross_creep;
    GroupedFieldErrorMetrics cross_displacement;
    std::size_t ni = 0, pi = 0;
    double backward = 0.0, native_creep = 0.0;
    std::array<std::size_t, 2> simultaneous_samples{};
    std::map<std::size_t, std::array<double, 2>> old_native;
    constexpr std::array<std::size_t, 4> point_map = {0, 1, 3, 2};
    const std::array<std::string, 4> component = {"rr", "zz", "hoop", "rz"};
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        if (!same_time(frame.time, .2 * static_cast<double>(step)))
            throw std::runtime_error("Creep ramp time grid changed");
        const double load = std::min(.2 * static_cast<double>(step), 1.0);
        const double old_load = std::min(.2 * static_cast<double>(step - 1), 1.0);
        const double backward_increment = .2e-4 * load * load * load;
        const double native_increment = step <= 2 ? .2e-4 * old_load * old_load * old_load : backward_increment;
        backward += backward_increment;
        native_creep += native_increment;
        const double sigma = 1e8 * load;
        const double expected_plastic = std::max(0.0, (sigma - 5e7) / 2e11);
        const double old_plastic = std::max(0.0, (1e8 * old_load - 5e7) / 2e11);
        const std::array<double, 2> accumulated = {backward, native_creep};
        std::array<double, 2> bottom_force{};
        for (std::size_t n = 0; n < frame.nodes.size(); ++n) {
            const auto& row = nodes.at(ni++);
            if (!same_time(row.at("time"), frame.time) || row.at("node") != static_cast<double>(n + 1))
                throw std::runtime_error("Creep ramp node history mismatch");
            const std::array<std::array<double, 2>, 2> u = {
                {{frame.nodal("displacement_r")[n], frame.nodal("displacement_z")[n]}, {row.at("ur"), row.at("uz")}}};
            cross_displacement.add(u[0].data(), u[1].data(), 2);
            for (std::size_t solver = 0; solver < 2; ++solver) {
                const std::array<double, 2> expected = {
                    frame.nodes[n][0] * (-.3 * sigma / 2e11 - .5 * (accumulated[solver] + expected_plastic)),
                    frame.nodes[n][1] * (sigma / 2e11 + accumulated[solver] + expected_plastic)};
                displacement[solver].add(u[solver].data(), expected.data(), 2);
                temperature[solver].add(solver == 0 ? frame.nodal("temperature")[n] : row.at("temperature"), 600.0);
                heat[solver].add(solver == 0 ? frame.nodal("reaction_heat_flux")[n] : row.at("reaction_heat"), 0.0);
                if (frame.nodes[n][1] == 0.0)
                    bottom_force[solver] += solver == 0 ? frame.nodal("reaction_force_z")[n] : row.at("rf_z");
            }
        }
        for (std::size_t solver = 0; solver < 2; ++solver)
            reaction[solver].add(bottom_force[solver], -std::acos(-1.0) * 1e-6 * sigma);
        if (frame.element("material_point_count").size() != 1)
            throw std::runtime_error("Creep ramp requires the independent single-element specimen");
        const auto count = static_cast<std::size_t>(frame.element("material_point_count")[0]);
        if (count != 1 && count != 4 && count != 9)
            throw std::runtime_error("Invalid creep ramp material point count");
        for (std::size_t q = 0; q < count; ++q) {
            const auto& row = points.at(pi++);
            if (!same_time(row.at("time"), frame.time) || row.at("element") != 1
                || row.at("point") != static_cast<double>(q + 1))
                throw std::runtime_error("Creep ramp material point history mismatch");
            const std::string suffix = "_q" + std::to_string(count == 9 ? q : point_map[q]);
            cross_creep.add(frame.element("equiv_creep" + suffix)[0], row.at("equiv_creep"));
            for (std::size_t solver = 0; solver < 2; ++solver) {
                equivalent[solver].add(solver == 0 ? frame.element("equiv_creep" + suffix)[0] : row.at("equiv_creep"),
                    accumulated[solver]);
                const double c = solver == 0 ? frame.element("equiv_creep" + suffix)[0] : row.at("equiv_creep");
                const double p = solver == 0 ? frame.element("equiv_plastic" + suffix)[0] : row.at("equiv_plastic");
                const double dc =
                    c - (solver == 0 ? frames[step - 1].element("equiv_creep" + suffix)[0] : old_native[q][0]);
                const double dp =
                    p - (solver == 0 ? frames[step - 1].element("equiv_plastic" + suffix)[0] : old_native[q][1]);
                plastic[solver].add(p, expected_plastic);
                creep_increment[solver].add(dc, solver == 0 ? backward_increment : native_increment);
                plastic_increment[solver].add(dp, expected_plastic - old_plastic);
                if (dc > 1e-8 && dp > 1e-8)
                    ++simultaneous_samples[solver];
                if (solver == 1)
                    old_native[q] = {c, p};
                const std::array<std::string, 4> prefixes = {"stress", "elastic", "creep", "plastic"};
                const std::array<std::array<double, 4>, 4> expected = {{{0, sigma, 0, 0},
                    {-.3 * sigma / 2e11, sigma / 2e11, -.3 * sigma / 2e11, 0},
                    {-.5 * accumulated[solver], accumulated[solver], -.5 * accumulated[solver], 0},
                    {-.5 * expected_plastic, expected_plastic, -.5 * expected_plastic, 0}}};
                const std::array<GroupedFieldErrorMetrics*, 4> metrics = {&stress[solver],
                    &elastic[solver],
                    &creep_tensor[solver],
                    &plastic_tensor[solver]};
                for (std::size_t field = 0; field < prefixes.size(); ++field) {
                    std::array<double, 4> actual{};
                    for (std::size_t c = 0; c < 4; ++c) {
                        const std::string key = prefixes[field] + "_" + component[c];
                        actual[c] = solver == 0 ? frame.element(key + suffix)[0] : row.at(key);
                        if (c == 3)
                            actual[c] *= std::sqrt(2.0);
                    }
                    metrics[field]->add(actual.data(), expected[field].data(), 4);
                }
            }
        }
    }
    bool passed = ni == nodes.size() && pi == points.size();
    for (std::size_t solver = 0; solver < 2; ++solver) {
        const std::string prefix = solver == 0 ? "fuelsim_backward_euler_" : "abaqus_native_integration_";
        passed = check_scalar(prefix + "temperature", temperature[solver], rz_temperature_zero_tolerance) && passed;
        passed = check_scalar(prefix + "reaction", reaction[solver], rz_force_zero_tolerance) && passed;
        passed = check_scalar(prefix + "equiv_creep", equivalent[solver], rz_strain_zero_tolerance) && passed;
        passed = check_scalar(prefix + "equiv_plastic", plastic[solver], rz_strain_zero_tolerance) && passed;
        passed = check_scalar(prefix + "creep_increment", creep_increment[solver], rz_strain_zero_tolerance) && passed;
        passed =
            check_scalar(prefix + "plastic_increment", plastic_increment[solver], rz_strain_zero_tolerance) && passed;
        passed = check_scalar(prefix + "heat_reaction", heat[solver], rz_heat_rate_zero_tolerance) && passed;
        passed = check_group(prefix + "displacement", displacement[solver], rz_displacement_zero_tolerance) && passed;
        passed = check_group(prefix + "stress", stress[solver], rz_stress_zero_tolerance) && passed;
        passed = check_group(prefix + "elastic", elastic[solver], rz_strain_zero_tolerance) && passed;
        passed = check_group(prefix + "creep", creep_tensor[solver], rz_strain_zero_tolerance) && passed;
        passed = check_group(prefix + "plastic", plastic_tensor[solver], rz_strain_zero_tolerance) && passed;
        const auto count = static_cast<std::size_t>(frames.back().element("material_point_count")[0]);
        passed = simultaneous_samples[solver] == 3 * count && passed;
        std::cout << prefix << "simultaneous_plastic_creep_samples=" << simultaneous_samples[solver] << '\n';
    }
    print_relative_metrics("cross_solver_diagnostic_equiv_creep", cross_creep);
    print_grouped_relative_metrics("cross_solver_diagnostic_displacement", cross_displacement);
    std::cout << "expected_backward_final_creep=" << backward << "\nexpected_abaqus_final_creep=" << native_creep
              << "\ncreep_integration_difference_is_expected=true\n";
    // Distinguish this analytic-method qualification from strict solver equality.
    passed = cross_creep.maximum_zero_reference_difference > 1e-8 && cross_creep.maximum_absolute_difference > 1e-6
             && passed;
    std::cout << "rz_creep_integration_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace

namespace fuelsim::test {
bool check_rz8_recovery_abaqus(const std::string& output_path,
    const std::string& node_path,
    const std::string& point_path,
    const std::string& contact_path) {
    bool passed = check_rz_abaqus(output_path, node_path, point_path, "probe", true);
    const auto rows = read_rows(contact_path);
    const auto frames = read_exodus_history(output_path);
    FieldErrorMetrics pressure, shear, gap, normal, tangent;
    std::size_t index = 0, open_frames = 0, mixed_frames = 0, near_zero_count = 0;
    double near_zero_difference = 0;
    for (const auto& frame : frames) {
        if (frame.time == 0)
            continue;
        std::size_t count = 0, closed = 0;
        while (index < rows.size() && same_time(rows[index].at("time"), frame.time)) {
            const auto& row = rows[index++];
            const auto label = static_cast<std::size_t>(std::llround(row.at("node")));
            if (label == 0 || label > frame.nodes.size() || row.at("node") != static_cast<double>(label))
                throw std::runtime_error("CAX8T recovery reference has an invalid node label");
            const auto n = label - 1;
            const double raw_pressure = frame.nodal("contact_pressure_interface")[n];
            const double raw_shear = frame.nodal("contact_tangential_traction_interface")[n];
            const double area = frame.nodal("contact_tributary_area_interface")[n];
            const double fn = frame.nodal("contact_normal_force_interface")[n];
            const double ft = frame.nodal("contact_tangential_force_interface")[n];
            if (std::abs(raw_pressure * area - fn) > 1e-9 || std::abs(raw_shear * area - ft) > 1e-9)
                throw std::runtime_error("CAX8T recovery changed raw force/traction output");
            const double p = frame.nodal("contact_recovered_pressure_interface")[n];
            const double s = frame.nodal("contact_recovered_shear_interface")[n];
            for (const auto& item :
                std::array<std::pair<double, double>, 2>{{{p, row.at("pressure")}, {s, row.at("shear")}}}) {
                if (!std::isfinite(item.first))
                    throw std::runtime_error("Nonfinite recovered contact traction");
                if (std::abs(item.second) <= rz_stress_zero_tolerance) {
                    ++near_zero_count;
                    near_zero_difference = std::max(near_zero_difference, std::abs(item.first - item.second));
                }
            }
            if (std::abs(row.at("pressure")) > rz_stress_zero_tolerance)
                pressure.add(p, row.at("pressure"));
            if (std::abs(row.at("shear")) > rz_stress_zero_tolerance)
                shear.add(s, row.at("shear"));
            gap.add(frame.nodal("contact_gap_interface")[n], row.at("gap"));
            normal.add(fn, std::hypot(row.at("normal_r"), row.at("normal_z")));
            tangent.add(std::abs(ft), std::hypot(row.at("tangential_r"), row.at("tangential_z")));
            ++count;
            if (row.at("gap") < 0)
                ++closed;
        }
        if (count != 5)
            throw std::runtime_error("CAX8T two-edge recovery must compare five nodes each frame");
        open_frames += closed == 0 ? 1 : 0;
        mixed_frames += closed > 0 && closed < count ? 1 : 0;
    }
    passed = check_scalar("rz8_recovered_pressure", pressure, rz_stress_zero_tolerance) && passed;
    passed = check_scalar("rz8_recovered_shear", shear, rz_stress_zero_tolerance) && passed;
    passed = check_scalar("rz8_recovery_gap", gap, rz_displacement_zero_tolerance) && passed;
    passed = check_scalar("rz8_recovery_normal_force", normal, rz_force_zero_tolerance) && passed;
    // Vanishing reaction noise in Abaqus is checked in absolute units above for
    // recovered stresses; force magnitudes use the same absolute roundoff gate.
    std::cout << "rz8_recovery_tangent_force_maximum_absolute_difference=" << tangent.maximum_absolute_difference
              << '\n';
    passed = tangent.maximum_absolute_difference < rz_force_zero_tolerance && passed;
    std::cout << "rz8_recovery_near_zero_stress_count=" << near_zero_count << '\n'
              << "rz8_recovery_near_zero_stress_maximum_absolute_difference=" << near_zero_difference << '\n'
              << "rz8_recovery_open_frames=" << open_frames << '\n'
              << "rz8_recovery_partial_frames=" << mixed_frames << '\n';
    return passed && index == rows.size() && index == 70 && open_frames == 1 && mixed_frames == 5
           && near_zero_difference <= rz_stress_zero_tolerance;
}

bool check_rz_abaqus(const std::string& output_path,
    const std::string& node_path,
    const std::string& point_path,
    const std::string& mechanisms,
    bool prescribed_state_absolute_check) {
    if (mechanisms == "creep-integration")
        return check_creep_integration(output_path, node_path, point_path);
    const auto nodes = read_rows(node_path), points = read_rows(point_path);
    const auto frames = read_exodus_history(output_path);
    const bool contact = mechanisms == "contact";
    const std::string reference_name = std::filesystem::path(node_path).filename().string();
    const bool qualified_heat = reference_name == "b150_cax4t_finite_thermal_operator_nodes.csv"
                                || reference_name == "b91_cax4rt_finite_probe_nodes.csv"
                                || reference_name == "b917_cax4rt_finite_thermal_operators_nodes.csv";
    const double heat_tolerance = qualified_heat ? 0.02 : rz_relative_tolerance;
    const bool qualified_recovery_heat =
        prescribed_state_absolute_check && mechanisms == "probe" && reference_name == "b114_cax8t_recovery_nodes.csv";
    const double heat_absolute_tolerance = qualified_recovery_heat ? 2e-11 : rz_heat_rate_zero_tolerance;
    std::cout << "rz_reaction_heat_relative_tolerance=" << heat_tolerance << '\n'
              << "rz_reaction_heat_absolute_tolerance=" << heat_absolute_tolerance << '\n';
    if (!contact && mechanisms != "plastic" && mechanisms != "creep" && mechanisms != "coupled"
        && mechanisms != "thermal" && mechanisms != "sliding" && mechanisms != "probe")
        throw std::runtime_error("Unknown RZ qualification mechanism");
    FieldErrorMetrics temperature, reaction, plastic, creep;
    FieldErrorMetrics pressure, gap, nodal_normal_force, normal_force, heat_rate;
    FieldErrorMetrics recovered_pressure, recovered_shear;
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
        if (frame.time == 0.0)
            continue;
        ++accepted;
        double actual_reaction = 0.0, reference_reaction = 0.0;
        double reference_force = 0.0, reference_heat = 0.0, reference_opposite_force = 0.0;
        double bottom_z = frame.nodes.front()[1];
        double outer_r = frame.nodes.front()[0];
        for (const auto& node : frame.nodes)
            bottom_z = std::min(bottom_z, node[1]);
        for (const auto& node : frame.nodes)
            outer_r = std::max(outer_r, node[0]);
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
            if ((frame.nodes[n][0] == 0.0 && (frame.nodes[n][1] == 0.0 || mechanisms == "thermal"))
                || (contact && frame.nodes[n][1] == bottom_z)
                || (mechanisms == "sliding" && frame.nodes[n][0] == outer_r)) {
                if (std::hypot(reference[0], reference[1]) > rz_displacement_zero_tolerance)
                    throw std::runtime_error("Abaqus violates the fixed origin displacement");
                reference = {0.0, 0.0};
            }
            displacement.add(actual.data(), reference.data(), actual.size());
            if (mechanisms == "probe") {
                const std::array<double, 2> force = {frame.nodal("reaction_force_r")[n],
                    frame.nodal("reaction_force_z")[n]};
                const std::array<double, 2> reference_force_vector = {row.at("rf_r"), row.at("rf_z")};
                support_reaction.add(force.data(), reference_force_vector.data(), 2);
                reaction_heat.add(frame.nodal("reaction_heat_flux")[n], row.at("reaction_heat"));
            }
            if (mechanisms == "sliding") {
                const std::array<double, 2> force = {frame.nodal("reaction_force_r")[n],
                    frame.nodal("reaction_force_z")[n]};
                const std::array<double, 2> reference_force_vector = {row.at("rf_r"), row.at("rf_z")};
                if (std::hypot(reference_force_vector[0], reference_force_vector[1]) > rz_force_zero_tolerance)
                    support_reaction.add(force.data(), reference_force_vector.data(), force.size());
                else
                    for (double component : force)
                        free_reaction.add(component, 0.0);
            }
            if (frame.nodes[n][1] == bottom_z) {
                if (!contact)
                    actual_reaction += frame.nodal("reaction_force_z")[n];
                reference_reaction += row.at("rf_z");
                if (contact)
                    reference_heat += row.at("reaction_heat");
            }
            if (contact) {
                reference_opposite_force += row.at("contact_force_z");
                if (frame.nodal("contact_projected_interface")[n] == 1.0) {
                    pressure.add(frame.nodal("contact_pressure_interface")[n], row.at("contact_pressure"));
                    if (std::find(frame.nodal_variable_names.begin(),
                            frame.nodal_variable_names.end(),
                            "temperature_active")
                        != frame.nodal_variable_names.end()) {
                        recovered_pressure.add(frame.nodal("contact_recovered_pressure_interface")[n],
                            row.at("contact_pressure"));
                        recovered_shear.add(frame.nodal("contact_recovered_shear_interface")[n], 0.0);
                    }
                    gap.add(frame.nodal("contact_gap_interface")[n], row.at("contact_gap"));
                    nodal_normal_force.add(frame.nodal("contact_normal_force_interface")[n],
                        -row.at("contact_force_z"));
                    reference_force -= row.at("contact_force_z");
                    if (row.at("contact_pressure") <= 1e6 || row.at("contact_gap") >= 0.0
                        || frame.nodal("contact_projected_interface")[n] != 1.0)
                        throw std::runtime_error("RZ thermal and mechanical contact must both be active");
                }
            }
        }
        if (contact) {
            normal_force.add(frame.global("contact_force_interface"), reference_force);
            heat_rate.add(frame.global("contact_heat_rate_interface"), reference_heat);
            if (std::abs(reference_opposite_force) > rz_force_zero_tolerance
                || std::abs(reference_force - reference_reaction) > rz_force_zero_tolerance)
                throw std::runtime_error("Abaqus contact reference violates force equilibrium");
        } else
            reaction.add(actual_reaction, reference_reaction);
        std::size_t elements = 0;
        for (const auto count : frame.block_element_counts)
            elements += count;
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
            if (point_count == 4
                && std::find(frame.nodal_variable_names.begin(), frame.nodal_variable_names.end(), "temperature_active")
                       != frame.nodal_variable_names.end())
                for (std::size_t q = 4; q < 9; ++q)
                    for (const auto& name : frame.element_variable_names)
                        if (name.size() >= 3 && name.substr(name.size() - 3) == "_q" + std::to_string(q)
                            && !std::isnan(frame.element(name)[e]))
                            throw std::runtime_error("Inactive CAX8RT output must be NaN");
            for (std::size_t aq = 0; aq < point_count; ++aq) {
                const auto& row = points.at(point_row++);
                if (!same_time(row.at("time"), frame.time) || row.at("element") != static_cast<double>(e + 1)
                    || row.at("point") != static_cast<double>(aq + 1))
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
    passed = check_scalar("rz_temperature", temperature, rz_temperature_zero_tolerance, prescribed_state_absolute_check)
             && passed;
    passed =
        check_group("rz_displacement", displacement, rz_displacement_zero_tolerance, prescribed_state_absolute_check)
        && passed;
    if (mechanisms == "sliding" || mechanisms == "probe") {
        passed = check_group("rz_support_reaction",
                     support_reaction,
                     rz_force_zero_tolerance,
                     prescribed_state_absolute_check)
                 && passed;
        if (mechanisms == "sliding")
            passed = check_scalar("rz_free_node_reaction", free_reaction, rz_force_zero_tolerance) && passed;
        if (mechanisms == "probe")
            passed = check_scalar("rz_reaction_heat",
                         reaction_heat,
                         heat_absolute_tolerance,
                         prescribed_state_absolute_check,
                         heat_tolerance)
                     && passed;
    }
    if (!contact)
        passed =
            check_scalar("rz_bottom_axial_reaction", reaction, rz_force_zero_tolerance, prescribed_state_absolute_check)
            && passed;
    for (std::size_t t = 0; t < (contact ? 1 : tensors.size()); ++t)
        passed = check_group("rz_" + prefixes[t] + "tensor",
                     tensors[t],
                     t == 0 ? rz_stress_zero_tolerance : rz_strain_zero_tolerance,
                     prescribed_state_absolute_check)
                 && passed;
    if (contact) {
        passed = check_scalar("rz_contact_pressure", pressure, rz_stress_zero_tolerance) && passed;
        passed = check_scalar("rz_contact_gap", gap, rz_displacement_zero_tolerance) && passed;
        passed = check_scalar("rz_contact_nodal_normal_force", nodal_normal_force, rz_force_zero_tolerance) && passed;
        passed = check_scalar("rz_contact_normal_force", normal_force, rz_force_zero_tolerance) && passed;
        passed = check_scalar("rz_contact_heat_rate", heat_rate, rz_heat_rate_zero_tolerance) && passed;
        if (recovered_pressure.value_count > 0) {
            passed =
                check_scalar("rz_contact_recovered_pressure", recovered_pressure, rz_stress_zero_tolerance) && passed;
            passed = check_scalar("rz_contact_recovered_shear", recovered_shear, rz_stress_zero_tolerance) && passed;
        }
    } else {
        passed = check_scalar("rz_equivalent_plastic_strain",
                     plastic,
                     rz_strain_zero_tolerance,
                     prescribed_state_absolute_check)
                 && passed;
        passed =
            check_scalar("rz_equivalent_creep_strain", creep, rz_strain_zero_tolerance, prescribed_state_absolute_check)
            && passed;
    }
    if (mechanisms == "plastic" || mechanisms == "coupled")
        passed = max_plastic > 1e-5 && passed;
    if (mechanisms == "creep" || mechanisms == "coupled")
        passed = max_creep > 1e-6 && passed;
    std::cout << "rz_accepted_frames=" << accepted << "\nrz_abaqus_qualification=" << (passed ? "passed" : "failed")
              << '\n';
    return passed;
}

bool check_rz_sliding_abaqus(const std::string& output_path,
    const std::string& node_path,
    const std::string& point_path,
    const std::string& contact_path,
    bool require_segment_crossing,
    bool quadratic_contact) {
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
        if (frame.time == 0.0)
            continue;
        double reference_normal = 0.0, reference_tangent = 0.0;
        std::array<double, 2> actual_normal_sum{}, reference_normal_sum{}, actual_tangent_sum{},
            reference_tangent_sum{};
        std::size_t frame_contacts = 0;
        while (row_index < rows.size() && same_time(rows[row_index].at("time"), frame.time)) {
            const auto& row = rows.at(row_index++);
            const auto label = static_cast<std::size_t>(std::llround(row.at("node")));
            if (row.at("node") != static_cast<double>(label) || label == 0 || label > frame.nodes.size())
                throw std::runtime_error("RZ sliding contact reference does not match output");
            ++frame_contacts;
            const auto n = label - 1;
            pressure.add(frame.nodal(quadratic_contact ? "contact_recovered_pressure_interface"
                                                       : "contact_pressure_interface")[n],
                row.at("pressure"));
            shear.add(frame.nodal(quadratic_contact ? "contact_recovered_shear_interface"
                                                    : "contact_tangential_traction_interface")[n],
                row.at("shear"));
            gap.add(frame.nodal("contact_gap_interface")[n], row.at("gap"));
            if (!quadratic_contact)
                slip.add(frame.nodal("contact_total_tangential_slip_interface")[n], row.at("slip"));
            else
                slip.add(std::abs(frame.nodal("contact_total_tangential_slip_interface")[n]), std::abs(row.at("slip")));
            normal_force.add(frame.nodal("contact_normal_force_interface")[n],
                std::hypot(row.at("normal_r"), row.at("normal_z")));
            if (!quadratic_contact)
                tangent_force.add(frame.nodal("contact_tangential_force_interface")[n],
                    std::copysign(std::hypot(row.at("tangential_r"), row.at("tangential_z")), row.at("shear")));
            else
                tangent_force.add(std::abs(frame.nodal("contact_tangential_force_interface")[n]),
                    std::hypot(row.at("tangential_r"), row.at("tangential_z")));
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
            const double dr = frame.nodes[b][0] + frame.nodal("displacement_r")[b] - frame.nodes[a][0]
                              - frame.nodal("displacement_r")[a];
            const double dz = frame.nodes[b][1] + frame.nodal("displacement_z")[b] - frame.nodes[a][1]
                              - frame.nodal("displacement_z")[a];
            const double length = std::hypot(dr, dz);
            const double normal = frame.nodal("contact_normal_force_interface")[n];
            const double tangent = frame.nodal("contact_tangential_force_interface")[n];
            const std::array<double, 2> actual_n = {-normal * dz / length, normal * dr / length};
            const double tangent_orientation = quadratic_contact ? 1.0 : -1.0;
            const std::array<double, 2> actual_t = {tangent_orientation * tangent * dr / length,
                tangent_orientation * tangent * dz / length};
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
            const bool reference_sliding = std::abs(friction_ratio - 1.0) < 1e-8;
            if (friction_ratio > 1.0 + 1e-8 || reference_sliding != is_sliding)
                throw std::runtime_error("Coulomb sticking/sliding state differs from Abaqus");
            if (previous_sliding[label] && !reference_sliding
                && (row.at("slip") - previous_slip[label]) * previous_shear[label] < 0.0)
                reverse_sticking = true;
            previous_sliding[label] = reference_sliding;
            previous_slip[label] = row.at("slip");
            previous_shear[label] = row.at("shear");
            sliding += is_sliding ? 1 : 0;
            sticking += is_sliding ? 0 : 1;
            forward = forward || (is_sliding && row.at("shear") > 0.0);
            reverse = reverse || (is_sliding && row.at("shear") < 0.0);
        }
        if (frame_contacts == 0)
            throw std::runtime_error("RZ sliding frame has no contact reference values");
        normal_force.add(frame.global("contact_force_interface"), reference_normal);
        if (!quadratic_contact)
            tangent_force.add(frame.global("contact_tangential_force_interface"), reference_tangent);
        normal_vector.add(actual_normal_sum.data(), reference_normal_sum.data(), 2);
        tangent_vector.add(actual_tangent_sum.data(), reference_tangent_sum.data(), 2);
    }
    std::vector<std::pair<std::string, const FieldErrorMetrics*>> fields = {{"gap", &gap},
        {"slip", &slip},
        {"normal_force", &normal_force}};
    fields.push_back({"pressure", &pressure});
    fields.push_back({"shear", &shear});
    fields.push_back({"tangent_force", &tangent_force});
    for (const auto& field : fields)
        passed = check_scalar("rz_sliding_" + field.first,
                     *field.second,
                     field.first == "gap" || field.first == "slip" ? rz_displacement_zero_tolerance
                                                                   : rz_force_zero_tolerance)
                 && passed;
    if (!quadratic_contact) {
        passed = check_group("rz_sliding_normal_force_vector", normal_vector, rz_force_zero_tolerance) && passed;
        passed = check_group("rz_sliding_tangent_force_vector", tangent_vector, rz_force_zero_tolerance) && passed;
    }
    std::size_t maximum_crossed = 0;
    for (const auto& entry : ownership) {
        const auto limits = std::minmax_element(entry.second.begin(), entry.second.end());
        maximum_crossed = std::max(maximum_crossed, static_cast<std::size_t>(*limits.second - *limits.first));
    }
    std::cout << "rz_sliding_maximum_segments_crossed=" << maximum_crossed << '\n'
              << "rz_sliding_sticking_samples=" << sticking << '\n'
              << "rz_sliding_sliding_samples=" << sliding << '\n'
              << "rz_sliding_reverse_sticking=" << reverse_sticking << '\n';
    passed = passed && row_index == rows.size() && (!require_segment_crossing || maximum_crossed >= 2) && sticking > 0
             && forward && reverse && reverse_sticking;
    std::cout << "rz_friction_abaqus_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace fuelsim::test
