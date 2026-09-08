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
        std::size_t ni = 0, pi = 0, ci = 0, active = 0, sliding = 0;
        double max_plastic = 0, max_creep = 0;
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
            previous_time = frame.time;
            if (frame.time == 0)
                continue;
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
            for (std::size_t n = 0; n < frame.nodes.size(); ++n) {
                const auto& row = nodes.at(ni++);
                require_time(row, frame.time);
                if (row.at("node") != static_cast<double>(n + 1))
                    throw std::runtime_error("Node labels differ");
                scalar["temperature"].add(frame.nodal("temperature")[n], row.at("temperature"));
                std::array<double, 2> a = {frame.nodal("displacement_r")[n], frame.nodal("displacement_z")[n]},
                                      b = {row.at("ur"), row.at("uz")};
                // Only prescribed zero components are normalized, after an
                // independent absolute check of both programs' constraint error.
                const std::array<bool, 2> fixed = {frame.nodes[n][0] == 0.0,
                    frame.nodes[n][1] == (integrated ? -0.001 : 0.0)};
                for (std::size_t component = 0; component < 2; ++component)
                    if (fixed[component]) {
                        scalar["prescribed_displacement_zero"].add(a[component], 0.0);
                        scalar["prescribed_displacement_zero"].add(b[component], 0.0);
                        a[component] = b[component] = 0.0;
                    }
                tensor["displacement"].add(a.data(), b.data(), 2);
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
                        tensor[prefix].add(a.data(), b.data(), 4);
                    }
                    for (const char* name : {"equiv_plastic", "equiv_creep"})
                        scalar[name].add(frame.element(std::string(name) + suffix)[e], row.at(name));
                    max_plastic = std::max(max_plastic, frame.element("equiv_plastic" + suffix)[e]);
                    max_creep = std::max(max_creep, frame.element("equiv_creep" + suffix)[e]);
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
                scalar["pressure"].add(
                    frame.nodal("contact_" + std::string(quadratic ? "recovered_pressure" : "pressure") + pair).at(n),
                    row.at("pressure"));
                scalar["gap"].add(frame.nodal("contact_gap" + pair).at(n), row.at("gap"));
                const double normal = std::hypot(row.at("normal_r"), row.at("normal_z"));
                scalar["normal_force"].add(frame.nodal("contact_normal_force" + pair).at(n), normal);
                total_normal += normal;
                active += frame.nodal("contact_pressure" + pair).at(n) > 0 ? 1U : 0U;
                if (integrated) {
                    scalar["shear"].add(
                        frame
                            .nodal(
                                "contact_" + std::string(quadratic ? "recovered_shear" : "tangential_traction") + pair)
                            .at(n),
                        row.at("shear"));
                    scalar["slip"].add(std::abs(frame.nodal("contact_total_tangential_slip" + pair).at(n)),
                        std::abs(row.at("slip")));
                    scalar["tangent_force"].add(std::abs(frame.nodal("contact_tangential_force" + pair).at(n)),
                        std::hypot(row.at("tangential_r"), row.at("tangential_z")));
                    sliding += frame.nodal("contact_sliding" + pair).at(n) == 1 ? 1U : 0U;
                }
            }
            if (ci == contact_begin)
                throw std::runtime_error("Missing contact frame");
            scalar["normal_resultant"].add(frame.global("contact_force_fuel_cladding"), total_normal);
        }
        bool passed = ni == nodes.size() && pi == points.size() && ci == contacts.size() && active > 0
                      && max_plastic > 0 && max_creep > 0
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
                                                                                                            : 1e-13;
            if (m.has_relative_norm())
                print_relative_metrics(name, m);
            else
                print_absolute_metrics(name, m);
            passed = passed && (!m.has_relative_norm() || relative_metrics_below(m, relative_tolerance))
                     && m.maximum_zero_reference_difference <= absolute;
        }
        for (const auto& [name, m] : tensor) {
            if (m.has_relative_norm())
                print_grouped_relative_metrics(name, m);
            std::cout << name << "_zero_reference_count=" << m.zero_reference_count << '\n'
                      << name << "_maximum_zero_reference_absolute_difference=" << m.maximum_zero_reference_difference
                      << '\n';
            passed = passed && (!m.has_relative_norm() || grouped_relative_metrics_below(m, relative_tolerance))
                     && m.maximum_zero_reference_difference <= (name == "stress" ? 1e-4 : 1e-13);
        }
        std::cout << "active_contact_samples=" << active << "\nsliding_contact_samples=" << sliding
                  << "\nmaximum_plastic=" << max_plastic << "\nmaximum_creep=" << max_creep << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
