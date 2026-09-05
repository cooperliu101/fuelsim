#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        result.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos) return result;
        begin = separator + 1;
    }
}

double final_csv_value(const std::string& path, const std::string& name) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE scalar reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE scalar reference is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE scalar reference is missing column: " + name);
    const std::size_t column = static_cast<std::size_t>(found - header.begin());
    std::vector<std::string> final_row;
    while (std::getline(input, line))
        if (!line.empty()) final_row = split_csv(line);
    if (column >= final_row.size())
        throw std::invalid_argument("MOOSE scalar reference has no final value for: " + name);
    return std::stod(final_row[column]);
}

struct CsvTable final {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};

CsvTable read_csv(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE integration-point reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE integration-point reference is empty: " + path);
    CsvTable result;
    result.header = split_csv(line);
    while (std::getline(input, line))
        if (!line.empty()) result.rows.push_back(split_csv(line));
    return result;
}

double csv_value(const CsvTable& table, const std::vector<std::string>& row, const std::string& name) {
    const auto found = std::find(table.header.begin(), table.header.end(), name);
    if (found == table.header.end())
        throw std::invalid_argument("MOOSE integration-point reference is missing column: " + name);
    const std::size_t column = static_cast<std::size_t>(found - table.header.begin());
    if (column >= row.size()) throw std::invalid_argument("MOOSE integration-point reference row is incomplete");
    return std::stod(row[column]);
}

struct RegionAverages final {
    double plastic = 0.0;
    double creep = 0.0;
    double stress = 0.0;
    double maximum_plastic = 0.0;
    double maximum_creep = 0.0;

    struct Point final {
        double radius;
        double axial_coordinate;
        double stress;
        double plastic;
        double creep;
    };

    std::vector<Point> points;
};

std::vector<RegionAverages::Point> read_cladding_points(
    const std::string& coordinate_path, const std::string& value_path) {
    const CsvTable coordinates = read_csv(coordinate_path);
    const CsvTable values = read_csv(value_path);
    if (values.rows.empty() || coordinates.rows.size() != values.rows.size() * 4)
        throw std::invalid_argument("MOOSE integration-point coordinate and value counts differ");
    std::vector<RegionAverages::Point> result;
    result.reserve(coordinates.rows.size());
    constexpr std::array<std::size_t, 4> moose_qp_from_fuelsim = {0, 1, 3, 2};
    for (const std::vector<std::string>& value_row : values.rows) {
        const std::size_t element = static_cast<std::size_t>(csv_value(values, value_row, "id"));
        for (const std::size_t qp : moose_qp_from_fuelsim) {
            const auto coordinate = std::find_if(
                coordinates.rows.begin(), coordinates.rows.end(), [&coordinates, element, qp](const auto& row) {
                    return static_cast<std::size_t>(csv_value(coordinates, row, "elem_id")) == element &&
                           static_cast<std::size_t>(csv_value(coordinates, row, "qp_id")) == qp;
                });
            if (coordinate == coordinates.rows.end())
                throw std::invalid_argument("MOOSE integration-point coordinate is missing");
            const std::string suffix = std::to_string(qp);
            result.push_back({csv_value(coordinates, *coordinate, "x"), csv_value(coordinates, *coordinate, "y"),
                csv_value(values, value_row, "stress_q" + suffix), csv_value(values, value_row, "plastic_q" + suffix),
                csv_value(values, value_row, "creep_q" + suffix)});
        }
    }
    return result;
}

RegionAverages region_averages(const fuelsim::test::ExodusResults& result) {
    const auto found = std::find(result.block_names.begin(), result.block_names.end(), "clad");
    if (found == result.block_names.end()) throw std::invalid_argument("M5.7 output has no clad block");
    const auto block = static_cast<std::size_t>(found - result.block_names.begin());
    std::size_t offset = 0;
    for (std::size_t b = 0; b < block; ++b) offset += result.block_element_counts[b];
    RegionAverages averages;
    double measure = 0.0;
    for (std::size_t e = offset; e < offset + result.block_element_counts[block]; ++e) {
        for (std::size_t q = 0; q < 4; ++q) {
            const auto suffix = "_q" + std::to_string(q);
            const double weight = result.element("reference_measure" + suffix)[e];
            const double rr = result.element("stress_rr" + suffix)[e], zz = result.element("stress_zz" + suffix)[e],
                         hoop = result.element("stress_hoop" + suffix)[e],
                         shear = result.element("stress_rz" + suffix)[e];
            const double mean = (rr + zz + hoop) / 3.0;
            const double stress = std::sqrt(1.5 * ((rr - mean) * (rr - mean) + (zz - mean) * (zz - mean) +
                                                      (hoop - mean) * (hoop - mean) + 2.0 * shear * shear));
            const double plastic = result.element("equiv_plastic" + suffix)[e];
            const double creep = result.element("equiv_creep" + suffix)[e];
            if (!(weight > 0.0) || !std::isfinite(weight))
                throw std::runtime_error("M5.7 reference integration measure is invalid");
            measure += weight;
            averages.stress += weight * stress;
            averages.plastic += weight * plastic;
            averages.creep += weight * creep;
            averages.maximum_plastic = std::max(averages.maximum_plastic, plastic);
            averages.maximum_creep = std::max(averages.maximum_creep, creep);
            averages.points.push_back({result.element("reference_r" + suffix)[e],
                result.element("reference_z" + suffix)[e], stress, plastic, creep});
        }
    }
    if (!(measure > 0.0)) throw std::runtime_error("M5.7 clad integration measure is empty");
    averages.stress /= measure;
    averages.plastic /= measure;
    averages.creep /= measure;
    return averages;
}

bool check_cladding_fields(const fuelsim::test::ExodusResults& result, const std::string& coordinate_path,
    const std::string& value_path, const std::string& scalar_path, double tolerance) {
    const auto cladding = region_averages(result);
    const auto reference = read_cladding_points(coordinate_path, value_path);
    if (cladding.points.size() != reference.size() || reference.empty())
        throw std::runtime_error("RZ cladding integration point count differs");
    fuelsim::test::FieldErrorMetrics stress, plastic, creep;
    double coordinate_error = 0.0;
    for (std::size_t i = 0; i < reference.size(); ++i) {
        const auto& a = cladding.points[i];
        const auto& b = reference[i];
        coordinate_error = std::max(
            {coordinate_error, std::abs(a.radius - b.radius), std::abs(a.axial_coordinate - b.axial_coordinate)});
        stress.add(a.stress, b.stress);
        plastic.add(a.plastic, b.plastic);
        creep.add(a.creep, b.creep);
    }
    fuelsim::test::print_relative_metrics("rz_clad_qp_equivalent_stress", stress);
    fuelsim::test::print_relative_metrics("rz_clad_qp_equivalent_plastic_strain", plastic);
    fuelsim::test::print_relative_metrics("rz_clad_qp_equivalent_creep_strain", creep);
    bool passed = check(coordinate_error < 1.0e-12 && fuelsim::test::relative_metrics_below(stress, tolerance) &&
                            fuelsim::test::relative_metrics_below(plastic, tolerance) &&
                            fuelsim::test::relative_metrics_below(creep, tolerance) && cladding.maximum_plastic > 0.0 &&
                            cladding.maximum_creep > 0.0,
        "RZ cladding all integration-point fields retain three configured error gates and active inelastic branches");
    for (const auto& entry : {std::pair<const char*, double>{"average_effective_plastic", cladding.plastic},
             {"average_effective_creep", cladding.creep}, {"average_vonmises_stress", cladding.stress}}) {
        const double expected = final_csv_value(scalar_path, entry.first);
        const double error = std::abs(entry.second - expected) / std::abs(expected);
        passed =
            check(expected != 0.0 && error < tolerance, std::string("RZ cladding weighted ") + entry.first) && passed;
        std::cout << "rz_clad_" << entry.first << "_relative_error=" << error << '\n';
    }
    return passed;
}
} // namespace

bool fuelsim::test::check_rz_integrated(const std::string& output_path, const std::string& coordinate_path,
    const std::string& value_path, const std::string& scalar_path) {
    const auto result = read_final_exodus_results(output_path);
    const auto initial = read_exodus_results(output_path, 1);
    if (initial.time != 0.0 || result.nodes != initial.nodes)
        throw std::runtime_error("M5.7 output must preserve the initial frame");
    bool passed = check_cladding_fields(result, coordinate_path, value_path, scalar_path, 5.0e-3);
    const std::string pair = "fuel_cladding";
    const auto& projected = result.nodal("contact_projected_" + pair);
    const auto& pressure = result.nodal("contact_pressure_" + pair);
    const auto& segment = result.nodal("contact_primary_segment_" + pair);
    const auto& initial_segment = initial.nodal("contact_primary_segment_" + pair);
    std::size_t nodes = 0, sliding = 0;
    double maximum_pressure = 0.0, maximum_excess = 0.0;
    for (std::size_t node = 0; node < projected.size(); ++node) {
        if (std::isnan(projected[node])) continue;
        ++nodes;
        passed = check(projected[node] == 1.0 && pressure[node] > 0.0 && segment[node] >= initial_segment[node] + 2.0,
                     "M5.7 every active contact node crosses two segments") &&
                 passed;
        sliding += result.nodal("contact_sliding_" + pair)[node] == 1.0 ? 1U : 0U;
        maximum_pressure = std::max(maximum_pressure, pressure[node]);
        maximum_excess = std::max(maximum_excess,
            std::abs(result.nodal("contact_tangential_traction_" + pair)[node]) - 0.002 * pressure[node]);
    }
    return check(nodes > 0 && sliding >= nodes - 1 && maximum_excess <= 1.0e-10 * maximum_pressure &&
                     std::abs(result.global("contact_heat_rate_" + pair)) > 0.0 &&
                     std::abs(result.global("contact_tangential_force_" + pair)) > 0.0,
               "M5.7 contact retains sliding, Coulomb cap and nonzero thermal/mechanical transfer") &&
           passed;
}

bool fuelsim::test::check_rz_pcmi(const std::string& output_path, const std::string& coordinate_path,
    const std::string& value_path, const std::string& scalar_path, const std::string& contact_path, double tolerance) {
    const auto result = read_final_exodus_results(output_path);
    bool passed = check_cladding_fields(result, coordinate_path, value_path, scalar_path, tolerance);
    const auto side_count = [&](const std::string& name) {
        const auto found = std::find(result.side_set_names.begin(), result.side_set_names.end(), name);
        if (found == result.side_set_names.end()) throw std::runtime_error("PCMI output is missing side set " + name);
        return result.side_set_sizes.at(static_cast<std::size_t>(found - result.side_set_names.begin()));
    };
    std::size_t elements = 0, offset = 0;
    for (std::size_t b = 0; b < result.block_names.size(); ++b) {
        elements += result.block_element_counts[b];
        if (result.block_names[b] == "fuel")
            for (std::size_t e = offset; e < offset + result.block_element_counts[b]; ++e)
                for (std::size_t q = 0; q < 4; ++q)
                    passed = check(result.element("equiv_plastic_q" + std::to_string(q))[e] == 0.0 &&
                                       result.element("equiv_creep_q" + std::to_string(q))[e] == 0.0,
                                 "PCMI fuel remains elastic at every material point") &&
                             passed;
        offset += result.block_element_counts[b];
    }
    passed = check(result.nodes.size() == 53 && elements == 34 && side_count("fuel_right") == 4 &&
                       side_count("clad_left") == 5,
                 "PCMI retains the tracked four-to-five nonmatching interface") &&
             passed;
    const auto& projected = result.nodal("contact_projected_fuel_cladding");
    std::size_t active = 0;
    double minimum_gap = 0.0;
    for (std::size_t node = 0; node < projected.size(); ++node) {
        if (std::isnan(projected[node])) continue;
        passed = check(projected[node] == 1.0, "PCMI taller cladding contains every fuel-node projection") && passed;
        active += result.nodal("contact_pressure_fuel_cladding")[node] > 0.0 ? 1U : 0U;
        minimum_gap = std::min(minimum_gap, result.nodal("contact_gap_fuel_cladding")[node]);
    }
    const auto reference = read_csv(contact_path);
    double force = 0.0;
    for (const auto& row : reference.rows)
        force += csv_value(reference, row, "contact_pressure") * csv_value(reference, row, "nodal_area");
    const double actual_force = result.global("contact_force_fuel_cladding");
    FieldErrorMetrics force_error;
    force_error.add(actual_force, force);
    print_relative_metrics("pcmi_total_contact_force", force_error);
    return check(active > 0 && minimum_gap < 0.0 && force > 0.0 && relative_metrics_below(force_error, tolerance),
               "PCMI develops contact pressure and retains its force-resultant gate") &&
           passed;
}
