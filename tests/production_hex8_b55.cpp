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
struct Point final {
    double x, y, z;
};

struct Tensor final {
    double xx, yy, zz, xy, yz, xz;
};

struct NodeReference final {
    std::size_t node;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t element;
    Point position;
    std::array<double, 3> heat_flux;
    Tensor stress;
    Tensor strain;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.5 row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.5 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.5 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 9) throw std::invalid_argument("Unexpected Abaqus B5.5 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 1, path);
            if (std::abs(fields[field]) < 1.0e-20) fields[field] = 0.0;
        }
        result.push_back({static_cast<std::size_t>(number(values, 0, path)), fields});
    }
    if (result.size() != 12) throw std::invalid_argument("Abaqus B5.5 nodal reference must contain twelve nodes");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.5 integration points: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
                "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,e13_engineering,"
                "e23_engineering,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.5 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 22)
            throw std::invalid_argument("Unexpected Abaqus B5.5 integration-point column count in " + path);
        result.push_back({static_cast<std::size_t>(number(values, 0, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
            {number(values, 9, path), number(values, 10, path), number(values, 11, path), number(values, 12, path),
                number(values, 14, path), number(values, 13, path)},
            {number(values, 15, path), number(values, 16, path), number(values, 17, path),
                0.5 * number(values, 18, path), 0.5 * number(values, 20, path), 0.5 * number(values, 19, path)}});
    }
    if (result.size() != 16)
        throw std::invalid_argument("Abaqus B5.5 integration-point reference must contain sixteen rows");
    return result;
}

} // namespace

namespace fuelsim::test {
bool check_hex8_b55(
    const std::string& output_path, const std::string& nodal_path, const std::string& integration_path) {
    const auto output = read_final_exodus_results(output_path);
    if (output.time != 1.0 || output.step_count != 2 || output.nodes.size() != 12)
        throw std::invalid_argument("B5.5 requires twelve nodes, the initial frame, and one one-second increment");
    bool passed = true;
    const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y", "displacement_z",
        "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
    std::array<FieldErrorMetrics, 8> nodal;
    std::set<std::size_t> seen;
    for (const auto& reference : read_nodes(nodal_path)) {
        if (reference.node < 1 || reference.node > 12 || !seen.insert(reference.node).second)
            throw std::invalid_argument("B5.5 reference node mapping is invalid");
        for (std::size_t field = 0; field < 8; ++field)
            nodal[field].add(output.nodal(nodal_names[field]).at(reference.node - 1), reference.fields[field]);
    }
    const std::array<double, 8> nodal_zero_tolerances = {1e-12, 1e-14, 1e-14, 1e-14, 2e-8, 1e-2, 1e-2, 1e-2};
    for (std::size_t field = 0; field < 8; ++field) {
        print_relative_metrics("b55_" + nodal_names[field], nodal[field]);
        passed = check(relative_metrics_below(nodal[field], 1e-3) &&
                           nodal[field].maximum_zero_reference_difference < nodal_zero_tolerances[field],
                     "B5.5 nodal field " + nodal_names[field] + " passes all three 0.1 percent metrics") &&
                 passed;
    }
    std::array<FieldErrorMetrics, 15> integration;
    const std::array<std::string, 15> names = {"heat_flux_x", "heat_flux_y", "heat_flux_z", "stress_xx", "stress_yy",
        "stress_zz", "stress_xy", "stress_yz", "stress_xz", "infinitesimal_strain_xx", "infinitesimal_strain_yy",
        "infinitesimal_strain_zz", "infinitesimal_strain_xy", "infinitesimal_strain_yz", "infinitesimal_strain_xz"};
    std::set<std::pair<std::size_t, std::size_t>> matched;
    double coordinate_error = 0;
    for (const auto& reference : read_integration(integration_path)) {
        if (reference.element < 1 || reference.element > 2)
            throw std::invalid_argument("B5.5 reference element is invalid");
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
        if (!std::isfinite(distance) || !matched.emplace(element, closest).second)
            throw std::invalid_argument("B5.5 integration point association is not unique");
        coordinate_error = std::max(coordinate_error, distance);
        const std::array<double, 15> expected = {reference.heat_flux[0], reference.heat_flux[1], reference.heat_flux[2],
            reference.stress.xx, reference.stress.yy, reference.stress.zz, reference.stress.xy, reference.stress.yz,
            reference.stress.xz, reference.strain.xx, reference.strain.yy, reference.strain.zz, reference.strain.xy,
            reference.strain.yz, reference.strain.xz};
        for (std::size_t field = 0; field < 15; ++field)
            integration[field].add(
                output.element(names[field] + "_q" + std::to_string(closest)).at(element), expected[field]);
    }
    for (std::size_t field = 0; field < 15; ++field) {
        print_relative_metrics("b55_" + names[field], integration[field]);
        passed = check(relative_metrics_below(integration[field], 1e-3) &&
                           integration[field].maximum_zero_reference_difference <
                               (field < 3 ? 1e-8 : (field < 9 ? 1e-2 : 1e-14)),
                     "B5.5 integration field " + names[field] + " passes all three 0.1 percent metrics") &&
                 passed;
    }
    return check(matched.size() == 16 && coordinate_error < 1e-7,
               "B5.5 matches all sixteen material integration points uniquely") &&
           passed;
}
} // namespace fuelsim::test
