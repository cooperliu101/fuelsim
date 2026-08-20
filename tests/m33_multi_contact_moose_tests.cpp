#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        fields.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos) return fields;
        begin = separator + 1;
    }
}

std::size_t column(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE contact CSV is missing column: " + name);
    return static_cast<std::size_t>(found - header.begin());
}

double value(const std::vector<std::string>& fields, std::size_t index, const std::string& path) {
    if (index >= fields.size()) throw std::invalid_argument("MOOSE contact CSV row is incomplete: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(fields[index], &parsed);
    if (parsed != fields[index].size() || !std::isfinite(result))
        throw std::invalid_argument("MOOSE contact CSV contains an invalid number: " + path);
    return result;
}

struct ContactReference final {
    std::vector<double> coordinates;
    std::vector<double> pressures;
    std::vector<double> radii;
    std::vector<double> axial_coordinates;
    double total_force = 0.0;
};

ContactReference read_contact_reference(const std::string& path) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE contact CSV: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE contact CSV is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const std::size_t pressure = column(header, "contact_pressure");
    const std::size_t coordinate = column(header, "y");
    const std::size_t radius = column(header, "x");
    const std::size_t radial_displacement = column(header, "disp_x");
    const std::size_t axial_displacement = column(header, "disp_y");
    std::vector<std::tuple<double, double, double, double>> values;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv(line);
        const double y = value(fields, coordinate, path);
        const double p = value(fields, pressure, path);
        values.push_back({y, value(fields, radius, path) + value(fields, radial_displacement, path),
            y + value(fields, axial_displacement, path), p});
    }
    if (values.empty()) throw std::invalid_argument("MOOSE contact CSV has no values: " + path);
    std::sort(values.begin(), values.end());
    ContactReference result;
    result.coordinates.reserve(values.size());
    result.pressures.reserve(values.size());
    result.radii.reserve(values.size());
    result.axial_coordinates.reserve(values.size());
    for (const auto& entry : values) {
        result.coordinates.push_back(std::get<0>(entry));
        result.radii.push_back(std::get<1>(entry));
        result.axial_coordinates.push_back(std::get<2>(entry));
        result.pressures.push_back(std::get<3>(entry));
    }
    for (std::size_t node = 0; node + 1 < result.coordinates.size(); ++node) {
        const double dr = result.radii[node + 1] - result.radii[node];
        const double dz = result.axial_coordinates[node + 1] - result.axial_coordinates[node];
        const double edge_length = std::hypot(dr, dz);
        const double area_node =
            2.0 * pi * 0.5 * edge_length * (2.0 * result.radii[node] + result.radii[node + 1]) / 3.0;
        const double area_next =
            2.0 * pi * 0.5 * edge_length * (2.0 * result.radii[node + 1] + result.radii[node]) / 3.0;
        result.total_force += result.pressures[node] * area_node + result.pressures[node + 1] * area_next;
    }
    return result;
}

bool run_comparison(const std::string& input_path, const std::string& nodal_reference_path,
    const std::string& first_contact_reference_path, const std::string& second_contact_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::steady)
        throw std::invalid_argument("M3.3 multi-contact comparison requires a steady input card");
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::read_exodus_quad4(definition.mesh_file);
    bool passed =
        check(source.element_blocks().size() == 3 && source.element_block_ids().size() == source.elements().size(),
            "M3.3 multi-contact reads the tracked three-block MOOSE mesh") &&
        check(definition.spatial.regions.size() == 3 && definition.spatial.contacts.size() == 2,
            "M3.3 multi-contact input contains three regions and two contact pairs");
    fuelsim::SteadyProblem problem(definition.spatial, source);
    const fuelsim::SolverOptions solver = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(problem,
        {definition.steady_execution.load_steps, definition.steady_execution.cutback_factor,
            definition.steady_execution.maximum_cutbacks_per_step, definition.steady_execution.minimum_load_increment},
        solver);
    passed = check(result.completed && result.solve.converged, "M3.3 multi-contact solve converges") && passed;
    const std::vector<fuelsim::test::NodalFieldReference> nodal_reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, result.solve.state, nodal_reference);
    constexpr double tolerance = 1.0e-2;
    passed = check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                 "M3.3 multi-contact compares every MOOSE node at matching coordinates") &&
             check(fuelsim::test::relative_metrics_below(fields.temperature, tolerance),
                 "M3.3 multi-contact temperature three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(fields.radial_displacement, tolerance),
                 "M3.3 multi-contact radial-displacement three full-field errors pass") &&
             check(fuelsim::test::relative_metrics_below(fields.axial_displacement, tolerance),
                 "M3.3 multi-contact axial-displacement three full-field errors pass") &&
             passed;
    const std::array<std::string, 2> contact_names = {"pellet_to_inner_clad", "inner_to_outer_clad"};
    const std::array<std::string, 2> contact_reference_paths = {
        first_contact_reference_path, second_contact_reference_path};
    for (std::size_t contact = 0; contact < contact_names.size(); ++contact) {
        const std::vector<fuelsim::ContactNodeSummary> actual =
            fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, contact, result.solve.state);
        const ContactReference reference = read_contact_reference(contact_reference_paths[contact]);
        const fuelsim::test::FieldErrorMetrics pressure =
            fuelsim::test::compare_moose_contact_pressure(actual, reference.pressures, reference.coordinates, 1.0e-12);
        const fuelsim::InterfaceSummary interface =
            fuelsim::rz::ProblemAccess::summarize_interface(problem, contact, result.solve.state);
        const double force_difference = std::abs(interface.total_contact_force - reference.total_force);
        const double force_scale = std::max({1.0, std::abs(interface.total_contact_force), reference.total_force});
        passed = check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                     "M3.3 multi-contact " + contact_names[contact] + " pressure errors pass") &&
                 check(interface.projected_contact_nodes == actual.size() &&
                           interface.active_contact_nodes == actual.size() && interface.unprojected_contact_nodes == 0,
                     "M3.3 multi-contact " + contact_names[contact] + " remains projected and active") &&
                 check(force_difference < tolerance * force_scale,
                     "M3.3 multi-contact " + contact_names[contact] + " total force agrees") &&
                 passed;
        fuelsim::test::print_relative_metrics("m33_" + contact_names[contact] + "_pressure", pressure);
        std::cout << "m33_" << contact_names[contact] << "_fuelsim_total_force=" << interface.total_contact_force
                  << '\n'
                  << "m33_" << contact_names[contact] << "_moose_total_force=" << reference.total_force << '\n';
    }
    fuelsim::test::print_relative_metrics("m33_multi_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m33_multi_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m33_multi_axial_displacement", fields.axial_displacement);
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_m33_multi_contact_moose_tests <case.fsi> <all-nodes.csv> "
                     "<first-contact.csv> <second-contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M3.3 multi-contact MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3], argv[4])) return 1;
        std::cout << "[PASS] M3.3 multi-contact MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.3 multi-contact comparison raised: " << error.what() << '\n';
        return 1;
    }
}
