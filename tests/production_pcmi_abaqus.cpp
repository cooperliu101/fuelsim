#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Row = std::map<std::string, double>;

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Missing Abaqus reference: " + path);
    std::string line, item;
    std::getline(input, line);
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    std::istringstream header(line);
    std::vector<std::string> names;
    while (std::getline(header, item, ','))
        names.push_back(item);
    std::vector<Row> rows;
    while (std::getline(input, line)) {
        std::istringstream stream(line);
        Row row;
        std::size_t c = 0;
        while (std::getline(stream, item, ',')) {
            const double value = std::stod(item);
            if (!std::isfinite(value))
                throw std::runtime_error("Nonfinite Abaqus reference");
            row.emplace(names.at(c++), value);
        }
        if (c != names.size())
            throw std::runtime_error("Incomplete Abaqus row");
        rows.push_back(std::move(row));
    }
    return rows;
}

void require_time(const Row& row, double time) {
    if (std::abs(row.at("time") - time) > 2e-6)
        throw std::runtime_error("Abaqus and Fuelsim time grids differ");
}

std::vector<double> integrated_reference_times(const std::string& reference_path) {
    const auto path = std::filesystem::path(reference_path).parent_path() / "b13_integrated_reference_times.txt";
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Missing integrated reference time selection");
    std::vector<double> times;
    double time = 0;
    while (input >> time) {
        if (!std::isfinite(time) || time <= 0 || (!times.empty() && time <= times.back()))
            throw std::runtime_error("Reference times must be finite, positive and strictly increasing");
        times.push_back(time);
    }
    if (!input.eof() || times.empty() || times.front() != 0.0625 || times.back() != 6.0)
        throw std::runtime_error("Incomplete integrated reference time selection");
    return times;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6 && argc != 7)
        return 2;
    try {
        using namespace fuelsim::test;
        const auto frames = read_exodus_history(argv[1]);
        const auto nodes = read_rows(argv[2]), points = read_rows(argv[3]), contacts = read_rows(argv[4]);
        const std::string group(argv[5]);
        if (group != "small" && group != "finite" && group != "integrated")
            throw std::invalid_argument("Unknown PCMI validation group");
        const bool integrated = group == "integrated";
        const bool small = group == "small";
        const double relative_tolerance = small ? 0.005 : 0.001;
        std::map<std::string, FieldErrorMetrics> scalar;
        std::map<std::string, GroupedFieldErrorMetrics> tensor;
        const bool finite_cax4t =
            group == "finite" && std::filesystem::path(argv[2]).filename() == "b13_finite_cax4t_nodes.csv";
        // User-qualified small nonzero fields retain both aggregate relative gates.
        // Every sample must satisfy its relative OR absolute pointwise gate.
        const std::map<std::string, double> absolute_limits = finite_cax4t
                                                                  ? std::map<std::string, double>{{"stress", 20.0},
                                                                        {"elastic", 3e-10},
                                                                        {"displacement", 2e-12},
                                                                        {"creep", 1e-18},
                                                                        {"equiv_creep", 1e-18},
                                                                        {"gap", 1e-12},
                                                                        {"pressure", 5.0},
                                                                        {"normal_force", 1e-4},
                                                                        {"normal_resultant", 1e-4},
                                                                        {"support_force", 1e-4}}
                                                                  : std::map<std::string, double>{};
        std::map<std::string, std::size_t> pointwise_failures, absolute_acceptances;
        std::map<std::string, double> maximum_relative_failure_difference, relative_failure_reference;
        const auto check_point = [&](const std::string& name, double difference, double reference) {
            if (reference == 0.0 || difference < relative_tolerance * reference)
                return;
            if (difference > maximum_relative_failure_difference[name]) {
                maximum_relative_failure_difference[name] = difference;
                relative_failure_reference[name] = reference;
            }
            const auto limit = absolute_limits.find(name);
            if (limit != absolute_limits.end() && difference < limit->second)
                ++absolute_acceptances[name];
            else
                ++pointwise_failures[name];
        };
        const auto add_scalar = [&](const std::string& name, double actual, double reference) {
            scalar[name].add(actual, reference);
            check_point(name, std::abs(actual - reference), std::abs(reference));
        };
        const auto add_tensor =
            [&](const std::string& name, const double* actual, const double* reference, std::size_t count) {
                tensor[name].add(actual, reference, count);
                double difference = 0.0, norm = 0.0;
                for (std::size_t component = 0; component < count; ++component) {
                    difference = std::hypot(difference, actual[component] - reference[component]);
                    norm = std::hypot(norm, reference[component]);
                }
                check_point(name, difference, norm);
            };
        std::size_t ni = 0, pi = 0, ci = 0, active = 0, sliding = 0;
        double max_plastic = 0, max_creep = 0;
        std::map<std::pair<std::size_t, std::size_t>, std::array<double, 4>> previous_inelastic;
        std::size_t simultaneous_samples = 0;
        double maximum_backward_euler_error = 0;
        const std::size_t expected_frames =
            argc == 7 ? std::stoul(argv[6]) + 1 : (integrated ? 97U : (small ? 12U : 21U));
        if (frames.size() != expected_frames || frames.front().time != 0
            || std::abs(frames.back().time - (integrated ? 6.0 : 20.0)) > 1e-12)
            throw std::runtime_error("Incomplete PCMI time history");
        const auto selected_times = integrated ? integrated_reference_times(argv[2]) : std::vector<double>{};
        std::size_t compared_frames = 0;
        double previous_time = -1;
        for (const auto& frame : frames) {
            if (!std::isfinite(frame.time) || frame.time <= previous_time)
                throw std::runtime_error("Production time history must be strictly increasing");
            const double dt = frame.time - previous_time;
            previous_time = frame.time;
            if (frame.time == 0)
                continue;
            // Conservation uses every production frame, including times omitted from the external references.
            add_scalar("thermal_balance_heat", frame.global("conservation_global_thermal_balance"), 0.0);
            add_scalar("interface_balance_heat", frame.global("conservation_interface_heat_imbalance"), 0.0);
            if (integrated) {
                if (std::abs(frame.time - std::round(frame.time / 0.0625) * 0.0625) > 1e-12)
                    throw std::runtime_error("Integrated production time grid differs");
                if (compared_frames == selected_times.size() || frame.time < selected_times[compared_frames] - 2e-6)
                    continue;
                if (std::abs(frame.time - selected_times[compared_frames]) > 2e-6)
                    throw std::runtime_error("Selected reference time is missing from production history");
            }
            ++compared_frames;
            const bool quadratic =
                std::find(frame.nodal_variable_names.begin(), frame.nodal_variable_names.end(), "temperature_active")
                != frame.nodal_variable_names.end();
            std::vector<std::array<bool, 2>> support(frame.nodes.size(), {false, false});
            std::vector<bool> temperature_support(frame.nodes.size(), false);
            for (std::size_t side = 0; side < frame.side_set_names.size(); ++side) {
                const auto& name = frame.side_set_names[side];
                for (const auto& face : frame.side_set_face_nodes[side])
                    for (const auto n : face) {
                        support[n][0] = support[n][0] || name == "fuel_left";
                        support[n][1] = support[n][1] || name == "fuel_bottom" || name == "clad_bottom"
                                        || (small && name == "clad_top");
                        temperature_support[n] = temperature_support[n] || name == "clad_right";
                    }
            }
            for (std::size_t n = 0; n < frame.nodes.size(); ++n) {
                const auto& row = nodes.at(ni++);
                require_time(row, frame.time);
                if (row.at("node") != static_cast<double>(n + 1))
                    throw std::runtime_error("Node labels differ");
                add_scalar("temperature", frame.nodal("temperature")[n], row.at("temperature"));
                std::array<double, 2> a = {frame.nodal("displacement_r")[n], frame.nodal("displacement_z")[n]},
                                      b = {row.at("ur"), row.at("uz")};
                // Only prescribed zero components are normalized, after an
                // independent absolute check of both programs' constraint error.
                const std::array<bool, 2> fixed = {frame.nodes[n][0] == 0.0,
                    frame.nodes[n][1] == (integrated ? -0.001 : 0.0)};
                for (std::size_t component = 0; component < 2; ++component)
                    if (fixed[component]) {
                        add_scalar("prescribed_displacement_zero", a[component], 0.0);
                        add_scalar("prescribed_displacement_zero", b[component], 0.0);
                        a[component] = b[component] = 0.0;
                    }
                add_tensor("displacement", a.data(), b.data(), 2);
                for (std::size_t component = 0; component < 2; ++component)
                    if (support[n][component])
                        add_scalar("support_force",
                            frame.nodal(component == 0 ? "reaction_force_r" : "reaction_force_z")[n],
                            row.at(component == 0 ? "rf_r" : "rf_z"));
                if (!quadratic || frame.nodal("temperature_active")[n] == 1) {
                    const double actual_heat = frame.nodal("reaction_heat_flux")[n];
                    if (temperature_support[n])
                        add_scalar("boundary_heat", actual_heat, row.at("reaction_heat"));
                    else {
                        add_scalar("free_heat", actual_heat, 0.0);
                        add_scalar("free_heat", row.at("reaction_heat"), 0.0);
                    }
                }
            }
            for (std::size_t e = 0; e < frame.element("material_point_count").size(); ++e) {
                const auto count = static_cast<std::size_t>(frame.element("material_point_count")[e]);
                if ((quadratic && count != 4 && count != 9) || (!quadratic && count != 1 && count != 4))
                    throw std::runtime_error("Unexpected active material point count");
                constexpr std::array<std::size_t, 4> map = {0, 1, 3, 2};
                for (std::size_t q = 0; q < count; ++q) {
                    const auto& row = points.at(pi++);
                    require_time(row, frame.time);
                    if (row.at("element") != static_cast<double>(e + 1)
                        || row.at("point") != static_cast<double>(q + 1))
                        throw std::runtime_error("Material point labels differ");
                    const std::string suffix = "_q" + std::to_string(count == 9 ? q : map.at(q));
                    for (const auto& prefix :
                        {std::string("stress"), std::string("elastic"), std::string("plastic"), std::string("creep")}) {
                        std::array<double, 4> a{}, b{};
                        std::size_t c = 0;
                        for (const char* name : {"rr", "zz", "hoop", "rz"}) {
                            const std::string key = prefix + "_" + name;
                            a[c] = frame.element(key + suffix)[e];
                            b[c] = row.at(key);
                            if (c == 3) {
                                a[c] *= std::sqrt(2.0);
                                b[c] *= std::sqrt(2.0);
                            }
                            ++c;
                        }
                        add_tensor(prefix, a.data(), b.data(), 4);
                    }
                    for (const char* name : {"equiv_plastic", "equiv_creep"})
                        add_scalar(name, frame.element(std::string(name) + suffix)[e], row.at(name));
                    max_plastic = std::max(max_plastic, frame.element("equiv_plastic" + suffix)[e]);
                    max_creep = std::max(max_creep, frame.element("equiv_creep" + suffix)[e]);
                    // B13 small: every clad point must activate both mechanisms at every step.
                    // This scope matches the independent original history audit, not other B13 loads.
                    if (small && row.at("element") >= 25 && row.at("element") <= 34) {
                        auto& previous = previous_inelastic[{e, q}];
                        const std::array<double, 4> current = {frame.element("equiv_plastic" + suffix)[e],
                            frame.element("equiv_creep" + suffix)[e],
                            row.at("equiv_plastic"),
                            row.at("equiv_creep")};
                        for (std::size_t solver = 0; solver < 2; ++solver) {
                            const double dp = current[2 * solver] - previous[2 * solver];
                            const double dc = current[2 * solver + 1] - previous[2 * solver + 1];
                            if (!(dp > 1e-12 && dc > 1e-12))
                                throw std::runtime_error(
                                    "B13 small clad point did not activate plasticity and creep together");
                            std::array<double, 4> stress{};
                            std::size_t component = 0;
                            for (const char* name : {"rr", "zz", "hoop", "rz"}) {
                                const std::string key = "stress_" + std::string(name);
                                stress[component++] = solver == 0 ? frame.element(key + suffix)[e] : row.at(key);
                            }
                            const double mises = std::hypot(std::hypot(stress[0] - stress[1], stress[1] - stress[2]),
                                                     std::hypot(stress[2] - stress[0], std::sqrt(6.0) * stress[3]))
                                                 / std::sqrt(2.0);
                            const double expected = dt * 1e-5 * std::pow(mises / 5e6, 3);
                            const double error = std::abs(dc / expected - 1);
                            if (!std::isfinite(error) || error >= 1e-3)
                                throw std::runtime_error("B13 small creep increment violates backward Euler");
                            maximum_backward_euler_error = std::max(maximum_backward_euler_error, error);
                        }
                        previous = current;
                        ++simultaneous_samples;
                    }
                }
            }
            double total_normal = 0;
            const std::size_t contact_begin = ci;
            std::vector<bool> contact_nodes(frame.nodes.size(), false);
            while (ci < contacts.size() && std::abs(contacts[ci].at("time") - frame.time) < 2e-6) {
                const auto& row = contacts[ci++];
                const auto n = static_cast<std::size_t>(row.at("node")) - 1;
                if (contact_nodes.at(n))
                    throw std::runtime_error("Duplicate contact node in reference frame");
                contact_nodes[n] = true;
                const std::string pair = "_fuel_cladding";
                add_scalar("pressure",
                    frame.nodal("contact_" + std::string(quadratic ? "recovered_pressure" : "pressure") + pair).at(n),
                    row.at("pressure"));
                add_scalar("gap", frame.nodal("contact_gap" + pair).at(n), row.at("gap"));
                const double normal = std::hypot(row.at("normal_r"), row.at("normal_z"));
                add_scalar("normal_force", frame.nodal("contact_normal_force" + pair).at(n), normal);
                total_normal += normal;
                active += frame.nodal("contact_pressure" + pair).at(n) > 0 ? 1U : 0U;
                if (integrated) {
                    add_scalar("shear",
                        frame
                            .nodal(
                                "contact_" + std::string(quadratic ? "recovered_shear" : "tangential_traction") + pair)
                            .at(n),
                        row.at("shear"));
                    add_scalar("slip",
                        std::abs(frame.nodal("contact_total_tangential_slip" + pair).at(n)),
                        std::abs(row.at("slip")));
                    add_scalar("tangent_force",
                        std::abs(frame.nodal("contact_tangential_force" + pair).at(n)),
                        std::hypot(row.at("tangential_r"), row.at("tangential_z")));
                    sliding += frame.nodal("contact_sliding" + pair).at(n) == 1 ? 1U : 0U;
                }
            }
            if (ci == contact_begin)
                throw std::runtime_error("Missing contact frame");
            add_scalar("normal_resultant", frame.global("contact_force_fuel_cladding"), total_normal);
        }
        bool passed = ni == nodes.size() && pi == points.size() && ci == contacts.size() && active > 0
                      && max_plastic > 0 && max_creep > 0
                      && (!small
                          || (previous_inelastic.size()
                                  == 10 * static_cast<std::size_t>(frames.back().element("material_point_count").at(24))
                              && simultaneous_samples == previous_inelastic.size() * (frames.size() - 1)))
                      && (!integrated || (sliding > 0 && compared_frames == selected_times.size()));
        std::cout << std::scientific << std::setprecision(12);
        std::cout << "reference_time_scope=" << (integrated ? "selected_times" : "complete_history") << '\n'
                  << "production_time_steps=" << frames.size() - 1 << '\n'
                  << "compared_reference_time_steps=" << compared_frames << '\n';
        std::cout << "relative_tolerance=" << relative_tolerance << '\n';
        for (const auto& [name, m] : scalar) {
            const double absolute = name == "temperature"                                                   ? 1e-11
                                    : name == "pressure" || name == "shear"                                 ? 1e-4
                                    : name.find("force") != std::string::npos || name == "normal_resultant" ? 1e-8
                                    : name.find("heat") != std::string::npos                                ? 1e-8
                                                                                                            : 1e-13;
            if (m.has_relative_norm())
                print_relative_metrics(name, m);
            else
                print_absolute_metrics(name, m);
            passed = passed
                     && (!m.has_relative_norm()
                         || (m.relative_l2() < relative_tolerance && m.relative_absolute_peak() < relative_tolerance
                             && pointwise_failures[name] == 0))
                     && m.maximum_zero_reference_difference <= absolute;
        }
        for (const auto& [name, m] : tensor) {
            if (m.has_relative_norm())
                print_grouped_relative_metrics(name, m);
            std::cout << name << "_zero_reference_count=" << m.zero_reference_count << '\n'
                      << name << "_maximum_zero_reference_absolute_difference=" << m.maximum_zero_reference_difference
                      << '\n';
            passed = passed
                     && (!m.has_relative_norm()
                         || (m.relative_l2() < relative_tolerance && m.relative_absolute_peak() < relative_tolerance
                             && pointwise_failures[name] == 0))
                     && m.maximum_zero_reference_difference <= (name == "stress" ? 1e-4 : 1e-13);
        }
        for (const auto& [name, limit] : absolute_limits)
            std::cout << name << "_qualified_pointwise_absolute_tolerance=" << limit << '\n'
                      << name << "_absolute_qualified_samples=" << absolute_acceptances[name] << '\n'
                      << name
                      << "_maximum_relative_failure_absolute_difference=" << maximum_relative_failure_difference[name]
                      << '\n'
                      << name << "_relative_failure_reference_norm=" << relative_failure_reference[name] << '\n';
        for (const auto& [name, count] : pointwise_failures)
            if (count != 0)
                std::cerr << "[FAIL] " << name << " pointwise samples outside both tolerances=" << count << '\n';
        std::cout << "active_contact_samples=" << active << "\nsliding_contact_samples=" << sliding
                  << "\nmaximum_plastic=" << max_plastic << "\nmaximum_creep=" << max_creep
                  << "\nsimultaneous_clad_samples=" << simultaneous_samples
                  << "\nmaximum_backward_euler_relative_error=" << maximum_backward_euler_error << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
