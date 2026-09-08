#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
struct Point final {
    double x, y, z;
};

struct Tensor final {
    double xx, yy, zz, xy, yz, xz;
};

struct NodeReference final {
    std::size_t time, node;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t time, element;
    Point position;
    std::array<double, 3> heat_flux;
    Tensor stress;
    Tensor strain;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete Abaqus B5.8 row in " + path);
    return std::stod(values[index]);
}

std::size_t integer_time(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > 4.0)
        throw std::invalid_argument("Abaqus B5.8 time is not one of 1, 2, 3, or 4 in " + path);
    return static_cast<std::size_t>(rounded);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.8 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.8 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 10)
            throw std::invalid_argument("Unexpected Abaqus B5.8 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 2, path);
            if (std::abs(fields[field]) < 1.0e-20)
                fields[field] = 0.0;
        }
        result.push_back(
            {integer_time(number(values, 0, path), path), static_cast<std::size_t>(number(values, 1, path)), fields});
    }
    if (result.size() != 48)
        throw std::invalid_argument("Abaqus B5.8 nodal reference must contain 48 rows");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.8 integration points: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
           "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,e13_engineering,"
           "e23_engineering,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.8 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 23)
            throw std::invalid_argument("Unexpected Abaqus B5.8 integration-point column count in " + path);
        result.push_back({integer_time(number(values, 0, path), path),
            static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
            {number(values, 7, path), number(values, 8, path), number(values, 9, path)},
            {number(values, 10, path),
                number(values, 11, path),
                number(values, 12, path),
                number(values, 13, path),
                number(values, 15, path),
                number(values, 14, path)},
            {number(values, 16, path),
                number(values, 17, path),
                number(values, 18, path),
                0.5 * number(values, 19, path),
                0.5 * number(values, 21, path),
                0.5 * number(values, 20, path)}});
    }
    if (result.size() != 64)
        throw std::invalid_argument("Abaqus B5.8 integration-point reference must contain 64 rows");
    return result;
}

} // namespace

namespace fuelsim::test {
bool check_hex8_b58(const std::string& output_path,
    const std::string& nodal_path,
    const std::string& integration_path) {
    std::map<std::size_t, ExodusResults> frames;
    for (std::size_t time = 1; time <= 4; ++time) {
        auto output = read_exodus_results(output_path, time + 1);
        if (output.time != static_cast<double>(time) || output.step_count != 5 || output.nodes.size() != 12)
            throw std::invalid_argument("B5.8 requires twelve nodes and four one-second increments");
        frames.emplace(time, std::move(output));
    }
    bool passed = true;
    const std::array<std::string, 8> nodal_names = {"temperature",
        "displacement_x",
        "displacement_y",
        "displacement_z",
        "reaction_heat_flux",
        "reaction_force_x",
        "reaction_force_y",
        "reaction_force_z"};
    std::array<FieldErrorMetrics, 8> nodal;
    for (const auto& frame : frames)
        passed = check(frame.second.global("conservation_relative_thermal_balance") < 1e-11,
                     "B5.8 thermal conservation closes at every accepted time")
                 && passed;
    std::set<std::pair<std::size_t, std::size_t>> seen;
    for (const auto& reference : read_nodes(nodal_path)) {
        if (reference.node < 1 || reference.node > 12 || !seen.emplace(reference.time, reference.node).second)
            throw std::invalid_argument("B5.8 reference node mapping is invalid");
        const auto& output = frames.at(reference.time);
        for (std::size_t field = 0; field < 8; ++field)
            nodal[field].add(output.nodal(nodal_names[field]).at(reference.node - 1), reference.fields[field]);
    }
    const std::array<double, 8> nodal_zero_tolerances = {1e-12, 1e-14, 1e-14, 1e-14, 2e-4, 2e-2, 2e-2, 2e-2};
    for (std::size_t field = 0; field < 8; ++field) {
        print_relative_metrics("b58_" + nodal_names[field], nodal[field]);
        passed = check(relative_metrics_below(nodal[field], 1e-3)
                           && nodal[field].maximum_zero_reference_difference < nodal_zero_tolerances[field],
                     "B5.8 nodal field " + nodal_names[field] + " passes all three 0.1 percent metrics")
                 && passed;
    }
    std::array<FieldErrorMetrics, 15> integration;
    const std::array<std::string, 15> names = {"heat_flux_x",
        "heat_flux_y",
        "heat_flux_z",
        "stress_xx",
        "stress_yy",
        "stress_zz",
        "stress_xy",
        "stress_yz",
        "stress_xz",
        "infinitesimal_strain_xx",
        "infinitesimal_strain_yy",
        "infinitesimal_strain_zz",
        "infinitesimal_strain_xy",
        "infinitesimal_strain_yz",
        "infinitesimal_strain_xz"};
    std::set<std::tuple<std::size_t, std::size_t, std::size_t>> matched;
    double coordinate_error = 0;
    for (const auto& reference : read_integration(integration_path)) {
        if (reference.element < 1 || reference.element > 2)
            throw std::invalid_argument("B5.8 reference element is invalid");
        const auto& output = frames.at(reference.time);
        const auto element = reference.element - 1;
        std::size_t closest = 0;
        double distance = std::numeric_limits<double>::infinity();
        for (std::size_t q = 0; q < 8; ++q) {
            const auto suffix = "_q" + std::to_string(q);
            const double difference =
                std::hypot(output.element("reference_x" + suffix).at(element) - reference.position.x,
                    output.element("reference_y" + suffix).at(element) - reference.position.y,
                    output.element("reference_z" + suffix).at(element) - reference.position.z);
            if (difference < distance) {
                closest = q;
                distance = difference;
            }
        }
        if (!std::isfinite(distance) || !matched.emplace(reference.time, element, closest).second)
            throw std::invalid_argument("B5.8 integration point association is not unique");
        coordinate_error = std::max(coordinate_error, distance);
        const std::array<double, 15> expected = {reference.heat_flux[0],
            reference.heat_flux[1],
            reference.heat_flux[2],
            reference.stress.xx,
            reference.stress.yy,
            reference.stress.zz,
            reference.stress.xy,
            reference.stress.yz,
            reference.stress.xz,
            reference.strain.xx,
            reference.strain.yy,
            reference.strain.zz,
            reference.strain.xy,
            reference.strain.yz,
            reference.strain.xz};
        for (std::size_t field = 0; field < 15; ++field)
            integration[field].add(output.element(names[field] + "_q" + std::to_string(closest)).at(element),
                expected[field]);
    }
    for (std::size_t field = 0; field < 15; ++field) {
        print_relative_metrics("b58_" + names[field], integration[field]);
        passed = check(relative_metrics_below(integration[field], 1e-3)
                           && integration[field].maximum_zero_reference_difference < 1e-10,
                     "B5.8 integration field " + names[field] + " passes all three 0.1 percent metrics")
                 && passed;
    }
    return check(matched.size() == 64 && coordinate_error < 1e-12,
               "B5.8 matches all sixty-four material integration points across four times uniquely")
           && passed;
}
} // namespace fuelsim::test
