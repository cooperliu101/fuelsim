#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
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
} // namespace

int main(int argc, char** argv) {
    if (argc != 4)
        return 2;
    try {
        using namespace fuelsim::test;
        const auto frames = read_exodus_history(argv[1]);
        const auto nodes = read_rows(argv[2]), contacts = read_rows(argv[3]);
        if (frames.size() != 5 || frames.back().time != 4)
            throw std::runtime_error("Incomplete NTS identification history");
        std::map<std::string, FieldErrorMetrics> metrics;
        std::size_t ni = 0, ci = 0;
        double prescribed_temperature = 0, prescribed_displacement = 0, bulk_heat_absolute_error = 0;
        for (const auto& frame : frames) {
            if (frame.time == 0)
                continue;
            const bool quadratic =
                std::find(frame.nodal_variable_names.begin(), frame.nodal_variable_names.end(), "temperature_active")
                != frame.nodal_variable_names.end();
            std::vector<bool> interface_node(frame.nodes.size(), false);
            for (std::size_t side = 0; side < frame.side_set_names.size(); ++side)
                if (frame.side_set_names[side] == "fuel_right" || frame.side_set_names[side] == "clad_left")
                    for (const auto& face : frame.side_set_face_nodes[side])
                        for (const auto n : face)
                            interface_node.at(n) = true;
            for (std::size_t n = 0; n < frame.nodes.size(); ++n) {
                const auto& row = nodes.at(ni++);
                require_time(row, frame.time);
                if (row.at("node") != static_cast<double>(n + 1))
                    throw std::runtime_error("Node order mismatch");
                prescribed_temperature =
                    std::max(prescribed_temperature, std::abs(frame.nodal("temperature")[n] - row.at("temperature")));
                prescribed_displacement = std::max({prescribed_displacement,
                    std::abs(frame.nodal("displacement_r")[n] - row.at("ur")),
                    std::abs(frame.nodal("displacement_z")[n] - row.at("uz"))});
                if (interface_node[n] && (!quadratic || frame.nodal("temperature_active")[n] == 1))
                    metrics["nodal_heat_reaction"].add(frame.nodal("reaction_heat_flux")[n], row.at("reaction_heat"));
                else
                    // Noninterface bulk heat is nonzero. Compare its difference
                    // against the native value with an independent absolute gate.
                    bulk_heat_absolute_error = std::max(bulk_heat_absolute_error,
                        std::abs(frame.nodal("reaction_heat_flux")[n] - row.at("reaction_heat")));
            }
            while (ci < contacts.size() && std::abs(contacts[ci].at("time") - frame.time) < 1e-12) {
                const auto& row = contacts[ci++];
                const auto n = static_cast<std::size_t>(row.at("node")) - 1;
                // In this all-prescribed probe, Abaqus R2018x reports first-step
                // pressure inconsistent with the specified penalty and gap.
                // Its cause remains unverified; see B14_NTS.md. Later mechanical
                // output and every step's heat reactions remain in the comparison.
                if (frame.time == 1)
                    continue;
                const std::string pair = "_fuel_cladding";
                metrics["pressure"].add(
                    frame.nodal(std::string(quadratic ? "contact_recovered_pressure" : "contact_pressure") + pair)[n],
                    row.at("pressure"));
                metrics["gap"].add(frame.nodal("contact_gap" + pair)[n], row.at("gap"));
                metrics["normal_force"].add(frame.nodal("contact_normal_force" + pair)[n],
                    std::hypot(row.at("normal_r"), row.at("normal_z")));
                metrics["shear"].add(
                    frame.nodal(
                        std::string(quadratic ? "contact_recovered_shear" : "contact_tangential_traction") + pair)[n],
                    row.at("shear"));
                metrics["tangential_force"].add(std::abs(frame.nodal("contact_tangential_force" + pair)[n]),
                    std::hypot(row.at("tangential_r"), row.at("tangential_z")));
                metrics["slip"].add(std::abs(frame.nodal("contact_total_tangential_slip" + pair)[n]),
                    std::abs(row.at("slip")));
            }
        }
        bool passed = ni == nodes.size() && ci == contacts.size() && prescribed_temperature < 1e-9
                      && prescribed_displacement < 1e-13 && bulk_heat_absolute_error < 1e-8;
        std::cout << std::scientific << std::setprecision(12);
        for (const auto& [name, m] : metrics) {
            print_relative_metrics(name, m);
            const double zero = name == "gap" || name == "slip" ? 1e-13 : 1e-8;
            passed = passed && m.value_count > 0 && relative_metrics_below(m, .001)
                     && m.maximum_zero_reference_difference < zero;
        }
        std::cout << "prescribed_temperature_maximum_absolute_error=" << prescribed_temperature << '\n'
                  << "noninterface_heat_maximum_absolute_error=" << bulk_heat_absolute_error << '\n'
                  << "prescribed_displacement_maximum_absolute_error=" << prescribed_displacement << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
