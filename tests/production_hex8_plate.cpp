#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    std::array<double, 3> point;
    std::array<double, 4> fields;
};

struct ElementReference final {
    std::array<double, 3> values;
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

std::size_t column(const std::vector<std::string>& header, const std::string& name, const std::string& path) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE CSV is missing '" + name + "': " + path);
    return static_cast<std::size_t>(found - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("MOOSE CSV row is incomplete: " + path);
    const double value = std::stod(values[index]);
    if (!std::isfinite(value)) throw std::invalid_argument("MOOSE CSV contains a non-finite value: " + path);
    return value;
}

std::vector<NodeReference> read_nodes(const std::string& path, std::size_t node_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE nodal reference is empty: " + path);
    const auto header = split_csv(line);
    const std::array<std::size_t, 8> columns = {column(header, "id", path), column(header, "x", path),
        column(header, "y", path), column(header, "z", path), column(header, "T", path), column(header, "disp_x", path),
        column(header, "disp_y", path), column(header, "disp_z", path)};
    std::vector<NodeReference> result(node_count);
    std::vector<bool> present(node_count, false);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        const std::size_t id = static_cast<std::size_t>(number(values, columns[0], path));
        if (id >= node_count) throw std::invalid_argument("MOOSE node ID is outside the Exodus mesh: " + path);
        const NodeReference candidate = {
            {number(values, columns[1], path), number(values, columns[2], path), number(values, columns[3], path)},
            {number(values, columns[4], path), number(values, columns[5], path), number(values, columns[6], path),
                number(values, columns[7], path)}};
        if (present[id]) {
            const NodeReference& existing = result[id];
            if (std::abs(existing.point[0] - candidate.point[0]) > 1.0e-12 ||
                std::abs(existing.point[1] - candidate.point[1]) > 1.0e-12 ||
                std::abs(existing.point[2] - candidate.point[2]) > 1.0e-12 ||
                !std::equal(existing.fields.begin(), existing.fields.end(), candidate.fields.begin(),
                    [](double left, double right) { return std::abs(left - right) < 1.0e-12; }))
                throw std::invalid_argument("MOOSE shared-node rows disagree: " + path);
            continue;
        }
        result[id] = candidate;
        present[id] = true;
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("MOOSE nodal reference does not cover every Exodus node: " + path);
    return result;
}

std::vector<ElementReference> read_elements(const std::string& path, std::size_t element_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE element reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE element reference is empty: " + path);
    const auto header = split_csv(line);
    const std::array<std::size_t, 4> columns = {column(header, "id", path), column(header, "stress_xx", path),
        column(header, "effective_plastic_strain", path), column(header, "effective_creep_strain", path)};
    std::vector<ElementReference> result(element_count);
    std::vector<bool> present(element_count, false);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        const std::size_t id = static_cast<std::size_t>(number(values, columns[0], path));
        if (id >= element_count || present[id]) throw std::invalid_argument("Invalid MOOSE element ID: " + path);
        for (std::size_t value = 0; value < result[id].values.size(); ++value)
            result[id].values[value] = number(values, columns[value + 1], path);
        present[id] = true;
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("MOOSE element reference does not cover every Exodus element: " + path);
    return result;
}

bool check_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below(metrics, 5.0e-3), name + " three errors are below 0.5 percent");
}

} // namespace

namespace fuelsim::test {
bool check_hex8_shared_plate(
    const std::string& output_path, const std::string& node_path, const std::string& element_path) {
    const auto result = read_final_exodus_results(output_path);
    const auto nodes = read_nodes(node_path, result.nodes.size());
    const bool transient = !element_path.empty();
    std::array<FieldErrorMetrics, 4> nodal;
    const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        for (std::size_t component = 0; component < 3; ++component)
            if (std::abs(result.nodes[node][component] - nodes[node].point[component]) > 1.0e-12)
                throw std::runtime_error("Shared-node reference coordinates differ");
        for (std::size_t field = 0; field < 4; ++field)
            nodal[field].add(result.nodal(names[field]).at(node), nodes[node].fields[field]);
    }
    bool passed = true;
    for (std::size_t field = 0; field < 4; ++field) {
        if (transient)
            passed = check_metrics("b36_" + names[field], nodal[field]) && passed;
        else {
            print_relative_metrics("b35_" + names[field], nodal[field]);
            if (field == 0)
                passed = check(relative_metrics_below(nodal[field], 1.0e-3), "B3.5 temperature errors pass") && passed;
            else
                passed = check(nodal[field].maximum_zero_reference_difference < 1.0e-10,
                             "B3.5 zero displacement references pass; nonzero relative errors remain diagnostic") &&
                         passed;
        }
    }
    if (!transient) {
        std::array<double, 3> sum{};
        std::array<std::size_t, 3> count{};
        for (std::size_t node = 0; node < nodes.size(); ++node)
            for (std::size_t plane = 0; plane < 3; ++plane)
                if (std::abs(result.nodes[node][0] - static_cast<double>(plane)) < 1.0e-12) {
                    sum[plane] += result.nodal("temperature").at(node);
                    ++count[plane];
                }
        for (std::size_t plane = 0; plane < 3; ++plane) {
            if (count[plane] == 0) throw std::runtime_error("B3.5 temperature plane is empty");
            sum[plane] /= static_cast<double>(count[plane]);
        }
        const double meat_flux = 10.0 * (sum[1] - sum[0]), clad_flux = 20.0 * (sum[2] - sum[1]);
        passed = check(std::abs(sum[1] - 500.0) < 1.0e-9, "B3.5 analytic interface temperature is 500 K") && passed;
        return check(std::abs(meat_flux - clad_flux) / std::max({std::abs(meat_flux), std::abs(clad_flux), 1.0}) <
                         1.0e-10,
                   "B3.5 interface heat flux is continuous") &&
               passed;
    }
    const std::size_t element_count = result.element("stress_xx_q0").size();
    const auto elements = read_elements(element_path, element_count);
    const std::array<std::string, 3> material_names = {"stress_xx", "equiv_plastic", "equiv_creep"};
    std::array<FieldErrorMetrics, 3> material;
    for (std::size_t element = 0; element < element_count; ++element)
        for (std::size_t field = 0; field < 3; ++field) {
            double average = 0.0;
            for (std::size_t q = 0; q < 8; ++q)
                average += result.element(material_names[field] + "_q" + std::to_string(q)).at(element) / 8.0;
            material[field].add(average, elements[element].values[field]);
        }
    for (std::size_t field = 0; field < 3; ++field)
        passed = check_metrics("b36_" + material_names[field], material[field]) && passed;
    return check(material[1].maximum_reference > 0.0 && material[2].maximum_reference > 0.0,
               "Both plasticity and creep are active in the meat-clad plate") &&
           passed;
}
} // namespace fuelsim::test
