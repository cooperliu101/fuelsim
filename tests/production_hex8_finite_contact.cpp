#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Point {
    double x, y, z;
};

struct NodeReference final {
    std::size_t step, source_node;
    Point reference;
    std::array<double, 3> displacement{}, reaction{};
};

struct ContactReference final {
    std::size_t step, source_node;
    Point current;
    std::array<double, 3> normal_force{}, tangential_force{};
    double slip_1, slip_2, gap, pressure;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete B4.3 CSV row: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B4.3 CSV number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double value = number(values, column, path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("Invalid B4.3 CSV index: " + path);
    return static_cast<std::size_t>(value);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.3 node reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,ux,uy,uz,rf_x,rf_y,rf_z")
        throw std::invalid_argument("Unexpected B4.3 node CSV header: " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path), index_value(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.3 contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,slip_1,"
                "slip_2,gap,pressure")
        throw std::invalid_argument("Unexpected B4.3 contact CSV header: " + path);
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
    return result;
}

std::array<double, 3> cross(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return {first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

double dot(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

} // namespace

namespace fuelsim::test {
bool check_hex8_finite_contact(
    const std::string& output_path, const std::string& node_path, const std::string& contact_path) {
    const auto frames = read_exodus_nodal_history(output_path);
    const auto nodes = read_nodes(node_path);
    const auto contacts = read_contact(contact_path);
    if (frames.size() != 21 || nodes.size() != 320 || contacts.size() != 80)
        throw std::invalid_argument("B4.3 requires twenty complete nodal/contact frames");
    const std::array<std::string, 3> axes = {"x", "y", "z"};
    std::array<FieldErrorMetrics, 3> displacement, normal, tangential, resultant, moment, center;
    FieldErrorMetrics gap, pressure;
    std::set<std::pair<std::size_t, std::size_t>> seen_nodes, seen_contacts;
    double max_coordinate = 0, max_current = 0, area_change = 0, normal_rotation = 0;
    bool active = true;
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        if (frame.nodes.size() != 16 || std::abs(frame.time - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("B4.3 production time or node count differs");
        for (const auto& ref : nodes) {
            if (ref.step != step) continue;
            if (ref.source_node >= 16 || !seen_nodes.emplace(step, ref.source_node).second)
                throw std::invalid_argument("B4.3 node association is not unique");
            const std::array<double, 3> reference_point = {ref.reference.x, ref.reference.y, ref.reference.z};
            for (std::size_t c = 0; c < 3; ++c) {
                max_coordinate =
                    std::max(max_coordinate, std::abs(frame.nodes[ref.source_node][c] - reference_point[c]));
                displacement[c].add(frame.nodal("displacement_" + axes[c]).at(ref.source_node), ref.displacement[c]);
            }
        }
        std::array<double, 3> force{}, ref_force{}, torque{}, ref_torque{}, weighted{}, ref_weighted{}, normal_sum{};
        double weight = 0, ref_weight = 0, area = 0;
        std::size_t count = 0, projected_count = 0;
        for (double projected : frame.nodal("contact_projected_interface"))
            if (projected == 1) ++projected_count;
        for (const auto& ref : contacts) {
            if (ref.step != step) continue;
            const auto n = ref.source_node;
            if ((n != 8 && n != 11 && n != 12 && n != 15) || !seen_contacts.emplace(step, n).second)
                throw std::invalid_argument("B4.3 contact association is not unique");
            ++count;
            active = active && frame.nodal("contact_projected_interface").at(n) == 1 &&
                     frame.nodal("contact_pressure_interface").at(n) > 0 &&
                     frame.nodal("contact_sliding_interface").at(n) == 0 && ref.pressure > 0;
            std::array<double, 3> point{}, expected_point = {ref.current.x, ref.current.y, ref.current.z}, f{}, ef{};
            for (std::size_t c = 0; c < 3; ++c) {
                point[c] = frame.nodes[n][c] + frame.nodal("displacement_" + axes[c]).at(n);
                max_current = std::max(max_current, std::abs(point[c] - expected_point[c]));
                const double nf = -frame.nodal("contact_normal_force_" + axes[c] + "_interface").at(n);
                const double tf = -frame.nodal("contact_tangential_force_" + axes[c] + "_interface").at(n);
                normal[c].add(nf, ref.normal_force[c]);
                tangential[c].add(tf, ref.tangential_force[c]);
                f[c] = nf + tf;
                ef[c] = ref.normal_force[c] + ref.tangential_force[c];
                force[c] += f[c];
                ref_force[c] += ef[c];
                normal_sum[c] += nf;
            }
            gap.add(frame.nodal("contact_gap_interface").at(n), ref.gap);
            pressure.add(frame.nodal("contact_pressure_interface").at(n), ref.pressure);
            area += frame.nodal("contact_tributary_area_interface").at(n);
            const auto m = cross(point, f), em = cross(expected_point, ef);
            const double w = std::sqrt(dot(f, f)), ew = std::sqrt(dot(ef, ef));
            weight += w;
            ref_weight += ew;
            for (std::size_t c = 0; c < 3; ++c) {
                torque[c] += m[c];
                ref_torque[c] += em[c];
                weighted[c] += w * point[c];
                ref_weighted[c] += ew * expected_point[c];
            }
        }
        if (count != 4 || projected_count != 4 || !(weight > 0 && ref_weight > 0))
            throw std::invalid_argument("B4.3 requires four uniquely projected active contact nodes");
        for (std::size_t c = 0; c < 3; ++c) {
            resultant[c].add(force[c], ref_force[c]);
            moment[c].add(torque[c], ref_torque[c]);
            center[c].add(weighted[c] / weight, ref_weighted[c] / ref_weight);
        }
        if (step == 20) {
            area_change = std::abs(area - 1);
            normal_rotation = std::hypot(normal_sum[1], normal_sum[2]) / std::abs(normal_sum[0]);
        }
    }
    bool passed = check(
        seen_nodes.size() == 320 && seen_contacts.size() == 80, "B4.3 covers all nodes and constraints at every step");
    const auto report = [&](const std::string& name, const FieldErrorMetrics& metric) {
        if (metric.has_relative_norm())
            print_relative_metrics("b43_" + name, metric);
        else
            print_absolute_metrics("b43_" + name, metric);
        return check((!metric.has_relative_norm() || relative_metrics_below(metric, 1e-2)) &&
                         metric.maximum_zero_reference_difference < 1e-7,
            "B4.3 " + name + " preserves its one-percent and zero-reference gates");
    };
    for (std::size_t c = 0; c < 3; ++c) {
        passed = report("displacement_" + axes[c], displacement[c]) && passed;
        passed = report("normal_force_" + axes[c], normal[c]) && passed;
        passed = report("tangential_force_" + axes[c], tangential[c]) && passed;
        passed = report("resultant_" + axes[c], resultant[c]) && passed;
        passed = report("moment_" + axes[c], moment[c]) && passed;
        passed = report("force_center_" + axes[c], center[c]) && passed;
    }
    passed = report("gap", gap) && passed;
    passed = report("pressure", pressure) && passed;
    std::cout << "b43_maximum_coordinate_difference=" << max_coordinate
              << "\nb43_maximum_current_coordinate_difference=" << max_current
              << "\nb43_final_secondary_area_change=" << area_change
              << "\nb43_final_normal_rotation_ratio=" << normal_rotation << '\n';
    return check(max_coordinate < 1e-14 && max_current < 1e-14, "B4.3 reference and current coordinates match") &&
           check(area_change > 0.02 && normal_rotation > 0.01, "B4.3 retains finite area change and normal rotation") &&
           check(active, "B4.3 retains four active frictionless constraints throughout the path") && passed;
}
} // namespace fuelsim::test
