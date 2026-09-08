#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double moose_relative_tolerance = 1.0e-3;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double relative_error(double actual, double expected) {
    if (!std::isfinite(actual) || !std::isfinite(expected) || expected == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - expected) / std::abs(expected);
}

bool check_scalar_metrics(const std::string& name, double actual, double reference, double tolerance) {
    const double relative_l2 = relative_error(actual, reference);
    const double relative_absolute_peak = std::abs(std::abs(actual) - std::abs(reference)) / std::abs(reference);
    const double maximum_pointwise_relative = relative_l2;
    std::cout << name << "_relative_l2=" << relative_l2 << '\n';
    std::cout << name << "_relative_absolute_peak=" << relative_absolute_peak << '\n';
    std::cout << name << "_maximum_pointwise_relative=" << maximum_pointwise_relative << '\n';
    return check(relative_l2 < tolerance && relative_absolute_peak < tolerance
                     && maximum_pointwise_relative < tolerance,
        name + " three MOOSE error metrics pass");
}

struct J2HistoryValue final {
    double time = 0.0;
    double axial_displacement = 0.0;
    double axial_plastic = 0.0;
    double axial_stress = 0.0;
    double effective_plastic = 0.0;
    double hoop_plastic = 0.0;
    double radial_plastic = 0.0;
};

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        fields.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos)
            return fields;
        begin = separator + 1;
    }
}

std::size_t csv_column(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end())
        throw std::invalid_argument("MOOSE history CSV is missing column '" + name + "'");
    return static_cast<std::size_t>(found - header.begin());
}

double csv_value(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    if (column >= fields.size())
        throw std::invalid_argument("MOOSE history CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("MOOSE history CSV contains an invalid number: " + path);
    return value;
}

std::vector<J2HistoryValue> read_j2_history(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE J2 history: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE J2 history is empty: " + path);
    const std::vector<std::string> header = split_csv_line(line);
    const std::size_t time = csv_column(header, "time");
    const std::size_t axial_displacement = csv_column(header, "axial_displacement");
    const std::size_t axial_plastic = csv_column(header, "axial_plastic");
    const std::size_t axial_stress = csv_column(header, "axial_stress");
    const std::size_t effective_plastic = csv_column(header, "effective_plastic");
    const std::size_t hoop_plastic = csv_column(header, "hoop_plastic");
    const std::size_t radial_plastic = csv_column(header, "radial_plastic");
    std::vector<J2HistoryValue> values;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv_line(line);
        const J2HistoryValue value{csv_value(fields, time, path),
            csv_value(fields, axial_displacement, path),
            csv_value(fields, axial_plastic, path),
            csv_value(fields, axial_stress, path),
            csv_value(fields, effective_plastic, path),
            csv_value(fields, hoop_plastic, path),
            csv_value(fields, radial_plastic, path)};
        if (value.time > 0.0)
            values.push_back(value);
    }
    if (values.empty())
        throw std::invalid_argument("MOOSE J2 history contains no accepted time steps: " + path);
    return values;
}

bool check_history_metrics(const std::string& name,
    const fuelsim::test::FieldErrorMetrics& metrics,
    double tolerance,
    double zero_reference_absolute_tolerance) {
    fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below(metrics, tolerance)
                     && metrics.maximum_zero_reference_difference <= zero_reference_absolute_tolerance,
        name + " history three MOOSE error metrics and zero references pass");
}

double average(const fuelsim::test::ExodusResults& result, const std::string& name) {
    double value = 0.0;
    for (std::size_t q = 0; q < 4; ++q)
        value += result.element(name + "_q" + std::to_string(q)).at(0) / 4.0;
    return value;
}

double top_displacement(const fuelsim::test::ExodusResults& result) {
    double maximum_z = -std::numeric_limits<double>::infinity();
    for (const auto& point : result.nodes)
        maximum_z = std::max(maximum_z, point[1]);
    double value = 0.0;
    std::size_t count = 0;
    for (std::size_t node = 0; node < result.nodes.size(); ++node)
        if (result.nodes[node][1] == maximum_z) {
            ++count;
            value += result.nodal("displacement_z")[node];
        }
    if (count == 0)
        throw std::runtime_error("RZ single-element output has no top nodes");
    return value / static_cast<double>(count);
}

double maximum_trace(const fuelsim::test::ExodusResults& result, const std::string& prefix) {
    double value = 0.0;
    for (std::size_t q = 0; q < 4; ++q) {
        const auto suffix = "_q" + std::to_string(q);
        value = std::max(value,
            std::abs(result.element(prefix + "_rr" + suffix)[0] + result.element(prefix + "_zz" + suffix)[0]
                     + result.element(prefix + "_hoop" + suffix)[0]));
    }
    return value;
}

bool check_unload_history(const std::string& output_path, const std::string& reference_path) {
    const auto final = fuelsim::test::read_final_exodus_results(output_path);
    const auto reference = read_j2_history(reference_path);
    std::vector<J2HistoryValue> actual;
    for (std::size_t step = 1; step <= final.step_count; ++step) {
        const auto frame = fuelsim::test::read_exodus_results(output_path, step);
        if (frame.time <= 0.0)
            continue;
        actual.push_back({frame.time,
            top_displacement(frame),
            average(frame, "plastic_zz"),
            average(frame, "stress_zz"),
            average(frame, "equiv_plastic"),
            average(frame, "plastic_hoop"),
            average(frame, "plastic_rr")});
    }
    bool passed = true;
    passed = check(actual.size() == reference.size() && actual.size() == 30,
                 "J2 unload-reload compares all 30 accepted MOOSE steps")
             && passed;
    if (actual.size() != reference.size())
        return false;
    fuelsim::test::FieldErrorMetrics axial_displacement;
    fuelsim::test::FieldErrorMetrics axial_plastic;
    fuelsim::test::FieldErrorMetrics axial_stress;
    fuelsim::test::FieldErrorMetrics effective_plastic;
    fuelsim::test::FieldErrorMetrics hoop_plastic;
    fuelsim::test::FieldErrorMetrics radial_plastic;
    double maximum_time_difference = 0.0;
    for (std::size_t step = 0; step < actual.size(); ++step) {
        maximum_time_difference = std::max(maximum_time_difference, std::abs(actual[step].time - reference[step].time));
        axial_displacement.add(actual[step].axial_displacement, reference[step].axial_displacement);
        axial_plastic.add(actual[step].axial_plastic, reference[step].axial_plastic);
        axial_stress.add(actual[step].axial_stress, reference[step].axial_stress);
        effective_plastic.add(actual[step].effective_plastic, reference[step].effective_plastic);
        hoop_plastic.add(actual[step].hoop_plastic, reference[step].hoop_plastic);
        radial_plastic.add(actual[step].radial_plastic, reference[step].radial_plastic);
    }
    std::cout << "m22_j2_unload_reload_maximum_time_difference=" << maximum_time_difference << '\n';
    passed = check(maximum_time_difference < 1.0e-12, "J2 unload-reload time coordinates match MOOSE") && passed;
    passed = check_history_metrics("m22_j2_unload_reload_axial_displacement",
                 axial_displacement,
                 moose_relative_tolerance,
                 0.0)
             && passed;
    passed =
        check_history_metrics("m22_j2_unload_reload_axial_plastic", axial_plastic, moose_relative_tolerance, 1.0e-14)
        && passed;
    passed = check_history_metrics("m22_j2_unload_reload_axial_stress", axial_stress, moose_relative_tolerance, 1.0e-6)
             && passed;
    passed = check_history_metrics("m22_j2_unload_reload_effective_plastic",
                 effective_plastic,
                 moose_relative_tolerance,
                 1.0e-14)
             && passed;
    passed = check_history_metrics("m22_j2_unload_reload_hoop_plastic", hoop_plastic, moose_relative_tolerance, 1.0e-14)
             && passed;
    passed =
        check_history_metrics("m22_j2_unload_reload_radial_plastic", radial_plastic, moose_relative_tolerance, 1.0e-14)
        && passed;
    const J2HistoryValue& first_peak = actual.at(9);
    const J2HistoryValue& unloaded = actual.at(19);
    const J2HistoryValue& reloaded = actual.at(29);
    std::cout << "m22_j2_unload_reload_first_peak_stress=" << first_peak.axial_stress << '\n';
    std::cout << "m22_j2_unload_reload_unloaded_stress=" << unloaded.axial_stress << '\n';
    std::cout << "m22_j2_unload_reload_first_peak_effective_plastic=" << first_peak.effective_plastic << '\n';
    std::cout << "m22_j2_unload_reload_unloaded_effective_plastic=" << unloaded.effective_plastic << '\n';
    std::cout << "m22_j2_unload_reload_reloaded_effective_plastic=" << reloaded.effective_plastic << '\n';
    passed = check(first_peak.axial_stress > 2.0e8 && unloaded.axial_stress < 0.0,
                 "J2 path reaches tensile plasticity then reverses stress during unload")
             && passed;
    passed = check(std::abs(unloaded.effective_plastic - first_peak.effective_plastic) < 1.0e-14,
                 "J2 elastic unload preserves committed equivalent plastic strain")
             && passed;
    passed = check(reloaded.effective_plastic > first_peak.effective_plastic,
                 "J2 reload activates additional plastic strain")
             && passed;
    passed = check(maximum_trace(final, "plastic") < 1.0e-12, "J2 unload-reload plastic strain remains trace-free")
             && passed;
    return passed;
}

} // namespace

bool fuelsim::test::check_rz_inelastic(const std::string& output_path,
    const std::string& branch,
    const std::string& history_path) {
    const auto result = read_final_exodus_results(output_path);
    bool passed = check(result.block_element_counts.size() == 1 && result.block_element_counts[0] == 1
                            && result.nodes.size() == 4,
        "M2.2 has exactly one region and one Quad4");
    for (const auto temperature : result.nodal("temperature"))
        passed = check(std::abs(temperature - 600.0) <= 1.0e-12, "M2.2 temperature remains 600 K") && passed;
    const double stress = average(result, "stress_zz"), plastic = average(result, "equiv_plastic"),
                 creep = average(result, "equiv_creep"), displacement = top_displacement(result);
    if (branch == "j2_unload_reload")
        return check_unload_history(output_path, history_path) && passed;
    std::array<double, 4> expected{};
    if (branch == "j2")
        expected = {201980198.0198, 0.000990099009901, 0.0, 2.0e-6};
    else if (branch == "norton")
        expected = {99998007.620195, 0.0, 9.9991036503199e-5, 5.9997808603447e-7};
    else if (branch == "coupled_traction")
        expected = {200999992.08159, 4.9999406119027e-4, 2.4564701918764e-4, 1.7506410289082e-6};
    else if (branch == "coupled_displacement")
        expected = {200963368.63564, 4.8168428343307e-4, 5.1349887359505e-4, 2.0e-6};
    else
        throw std::invalid_argument("Unknown M2.2 branch: " + branch);
    const std::array<double, 4> values = {stress, plastic, creep, displacement};
    const std::array<std::string, 4> names = {"axial_stress",
        "equivalent_plastic",
        "equivalent_creep",
        "top_displacement"};
    for (std::size_t f = 0; f < 4; ++f)
        if (expected[f] != 0.0)
            passed =
                check_scalar_metrics("m22_" + branch + "_" + names[f], values[f], expected[f], moose_relative_tolerance)
                && passed;
    if (branch == "j2" || branch == "coupled_traction")
        passed = check(maximum_trace(result, "plastic") < 1.0e-12, "M2.2 plastic strain remains trace-free") && passed;
    if (branch == "norton" || branch == "coupled_traction")
        passed = check(maximum_trace(result, "creep") < 1.0e-12, "M2.2 creep strain remains trace-free") && passed;
    if (branch == "coupled_traction")
        passed = check(relative_error(plastic, 5.0e-4) < 1.0e-8 && relative_error(creep, 2.4564818025e-4) < 1.0e-8
                           && relative_error(displacement, 1.75064818025e-6) < 1.0e-8,
                     "Coupled traction retains its independent analytic history gate")
                 && passed;
    return passed;
}
