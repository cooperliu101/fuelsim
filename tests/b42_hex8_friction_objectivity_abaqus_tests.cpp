#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct PathStep final {
    std::string name;
    double angle;
    std::array<double, 3> relative_displacement;
    bool sliding;
};

struct ReferenceNode final {
    std::size_t step, source_node;
    fuelsim::CartesianPoint3 current_point;
    std::array<double, 3> normal_force{}, tangential_force{};
    double slip_1, slip_2, gap, pressure;
};

const std::vector<PathStep> path = [] {
    std::vector<PathStep> result = {
        {"STICK", 0.0, {-0.01, 3.0e-5, 4.0e-5}, false}, {"SLIDE", 0.0, {-0.01, 1.2e-3, 1.6e-3}, true}};
    for (std::size_t increment = 1; increment <= 20; ++increment) {
        std::ostringstream name;
        name << "ROTATE_" << std::setw(2) << std::setfill('0') << increment;
        result.push_back({name.str(), 0.35 * static_cast<double>(increment) / 20.0, {-0.01, 1.2e-3, 1.6e-3}, true});
    }
    result.push_back({"RESTICK", 0.35, {-0.01, 1.17e-3, 1.56e-3}, false});
    return result;
}();

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::array<double, 3> rotate(double angle, const std::array<double, 3>& value) {
    const double cosine = std::cos(angle), sine = std::sin(angle);
    return {cosine * value[0] - sine * value[1], sine * value[0] + cosine * value[1], value[2]};
}

fuelsim::CartesianPoint3 current_point(
    const fuelsim::CartesianPoint3& reference, const PathStep& step, bool secondary) {
    const std::array<double, 3> shifted = {reference.x + (secondary ? step.relative_displacement[0] : 0.0),
        reference.y + (secondary ? step.relative_displacement[1] : 0.0),
        reference.z + (secondary ? step.relative_displacement[2] : 0.0)};
    const std::array<double, 3> rotated = rotate(step.angle, shifted);
    return {rotated[0], rotated[1], rotated[2]};
}

void write_labels(std::ofstream& output, const std::vector<std::size_t>& labels) {
    for (std::size_t index = 0; index < labels.size(); ++index) {
        output << labels[index] + 1;
        if ((index + 1) % 16 == 0 || index + 1 == labels.size())
            output << '\n';
        else
            output << ", ";
    }
}

void write_abaqus_input(const std::string& output_path, const fuelsim::UnstructuredHex8Mesh& mesh) {
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("Could not write B4.2 Abaqus input: " + output_path);
    output << std::setprecision(16) << "*Heading\n"
           << "** B4.2 C3D8 finite-sliding friction objectivity and resticking.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    output << "*Element, type=C3D8, elset=PRIMARY\n1";
    for (const std::size_t node : mesh.elements().at(0).nodes) output << ", " << node + 1;
    output << "\n*Element, type=C3D8, elset=SECONDARY\n2";
    for (const std::size_t node : mesh.elements().at(1).nodes) output << ", " << node + 1;
    output << "\n*Nset, nset=PRIMARY_ALL\n";
    write_labels(output, mesh.node_sets().at(0).nodes);
    output << "*Nset, nset=SECONDARY_ALL\n";
    write_labels(output, mesh.node_sets().at(1).nodes);
    output << "*Nset, nset=SECONDARY_CONTACT_NODES\n";
    write_labels(output, mesh.node_sets().at(2).nodes);
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.25\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e5,\n"
           << "*Friction, slip tolerance=1.e-4\n0.5,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n";
    for (const PathStep& step : path) {
        output << "*Step, name=" << step.name << ", nlgeom=YES, inc=1\n"
               << "*Static\n1., 1., 1., 1.\n"
               << "*Boundary, op=NEW\n";
        for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
            const bool secondary = node >= mesh.elements().at(1).nodes.front();
            const fuelsim::CartesianPoint3& reference = mesh.nodes()[node];
            const fuelsim::CartesianPoint3 current = current_point(reference, step, secondary);
            output << node + 1 << ", 1, 1, " << current.x - reference.x << '\n'
                   << node + 1 << ", 2, 2, " << current.y - reference.y << '\n'
                   << node + 1 << ", 3, 3, " << current.z - reference.z << '\n';
        }
        output << "*Output, field, frequency=1\n"
               << "*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n"
               << "*End Step\n";
    }
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete B4.2 CSV row: " + input_path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B4.2 CSV number: " + input_path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    const double value = number(values, column, input_path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("Invalid B4.2 CSV index: " + input_path);
    return static_cast<std::size_t>(value);
}

std::vector<std::vector<ReferenceNode>> read_reference(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Could not read B4.2 Abaqus reference: " + input_path);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,"
                                              "tangential_z,slip_1,slip_2,gap,pressure")
        throw std::invalid_argument("Unexpected B4.2 CSV header: " + input_path);
    std::vector<std::vector<ReferenceNode>> result(path.size());
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        const std::size_t step = index_value(values, 0, input_path);
        if (step == 0 || step > path.size()) throw std::invalid_argument("Invalid B4.2 CSV step: " + input_path);
        result[step - 1].push_back({step - 1, index_value(values, 1, input_path),
            {number(values, 2, input_path), number(values, 3, input_path), number(values, 4, input_path)},
            {number(values, 5, input_path), number(values, 6, input_path), number(values, 7, input_path)},
            {number(values, 8, input_path), number(values, 9, input_path), number(values, 10, input_path)},
            number(values, 11, input_path), number(values, 12, input_path), number(values, 13, input_path),
            number(values, 14, input_path)});
    }
    for (const auto& values : result)
        if (values.size() != 4) throw std::invalid_argument("Each B4.2 step must contain four contact nodes");
    return result;
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
}

fuelsim::SpatialDefinition definition(fuelsim::StrainFormulation strain) {
    fuelsim::SpatialDefinition result;
    result.regions = {{"primary", "primary", material(), 0.0, 300.0, -1, "", strain},
        {"secondary", "secondary", material(), 0.0, 300.0, -1, "", strain}};
    fuelsim::ContactDefinition contact;
    contact.name = "objective_interface";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.mechanical = true;
    contact.penalty = 1.0e5;
    contact.friction_coefficient = 0.5;
    contact.friction_elastic_slip = 1.0e-4;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    result.contacts.push_back(contact);
    return result;
}

std::vector<double> prescribed_state(const fuelsim::SteadyProblem& problem, const PathStep& step) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> result = problem.initial_state();
    for (std::size_t region = 0; region < 2; ++region) {
        const fuelsim::Hex8RegionMesh& mesh = spatial.region_mesh(region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const fuelsim::CartesianPoint3& reference = mesh.nodes()[local];
            const fuelsim::CartesianPoint3 current = current_point(reference, step, region == 1);
            const std::size_t global = spatial.global_node(region, local);
            result[spatial.dof(fuelsim::Field::displacement_x, global)] = current.x - reference.x;
            result[spatial.dof(fuelsim::Field::displacement_y, global)] = current.y - reference.y;
            result[spatial.dof(fuelsim::Field::displacement_z, global)] = current.z - reference.z;
        }
    }
    return result;
}

bool compare_case(const fuelsim::UnstructuredHex8Mesh& mesh, const std::vector<std::vector<ReferenceNode>>& reference,
    fuelsim::StrainFormulation strain, const std::string& prefix) {
    fuelsim::SteadyProblem problem(definition(strain), mesh);
    const std::vector<std::size_t> source_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    fuelsim::test::FieldErrorMetrics normal, tangent_first, tangent_second, gap, pressure, resultant, slip;
    double maximum_coordinate_difference = 0.0, maximum_objective_history_error = 0.0;
    std::array<double, 3> sliding_history_before_rotation{};
    bool physical_state_sequence_matches = true, reported_end_states_match = true, biaxial_history = true;
    for (std::size_t step = 0; step < path.size(); ++step) {
        const std::vector<double> state = prescribed_state(problem, path[step]);
        problem.validate_state(state);
        const std::vector<fuelsim::CartesianContactNodeSummary> actual =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        if (actual.size() != source_nodes.size()) throw std::invalid_argument("B4.2 contact-node count differs");
        std::array<double, 3> actual_resultant{}, reference_resultant{};
        const std::array<double, 3> normal_direction = rotate(path[step].angle, {1.0, 0.0, 0.0}),
                                    first_direction = rotate(path[step].angle, {0.0, 1.0, 0.0}),
                                    second_direction = {0.0, 0.0, 1.0};
        for (std::size_t node = 0; node < actual.size(); ++node) {
            const auto found = std::find_if(reference[step].begin(), reference[step].end(),
                [&](const ReferenceNode& value) { return value.source_node == source_nodes[node]; });
            if (found == reference[step].end()) throw std::invalid_argument("B4.2 source-node mapping is incomplete");
            const fuelsim::CartesianPoint3 expected_point =
                current_point(mesh.nodes().at(source_nodes[node]), path[step], true);
            maximum_coordinate_difference =
                std::max({maximum_coordinate_difference, std::abs(expected_point.x - found->current_point.x),
                    std::abs(expected_point.y - found->current_point.y),
                    std::abs(expected_point.z - found->current_point.z)});
            reported_end_states_match = reported_end_states_match && actual[node].projected &&
                                        (step != 0 || !actual[node].sliding) &&
                                        (step + 1 != path.size() || !actual[node].sliding);
            double actual_normal_component = 0.0, reference_normal_component = 0.0, actual_tangent_first = 0.0,
                   reference_tangent_first = 0.0, actual_tangent_second = 0.0, reference_tangent_second = 0.0;
            for (std::size_t component = 0; component < 3; ++component) {
                const double actual_normal = -actual[node].normal_contact_force[component],
                             actual_tangent = -actual[node].tangential_contact_force[component];
                actual_normal_component += normal_direction[component] * actual_normal;
                reference_normal_component += normal_direction[component] * found->normal_force[component];
                actual_tangent_first += first_direction[component] * actual_tangent;
                reference_tangent_first += first_direction[component] * found->tangential_force[component];
                actual_tangent_second += second_direction[component] * actual_tangent;
                reference_tangent_second += second_direction[component] * found->tangential_force[component];
                actual_resultant[component] += actual_normal + actual_tangent;
                reference_resultant[component] += found->normal_force[component] + found->tangential_force[component];
            }
            normal.add(actual_normal_component, reference_normal_component);
            tangent_first.add(actual_tangent_first, reference_tangent_first);
            tangent_second.add(actual_tangent_second, reference_tangent_second);
            const double actual_ratio =
                std::hypot(actual_tangent_first, actual_tangent_second) / (0.5 * std::abs(actual_normal_component));
            const double reference_ratio = std::hypot(reference_tangent_first, reference_tangent_second) /
                                           (0.5 * std::abs(reference_normal_component));
            physical_state_sequence_matches =
                physical_state_sequence_matches &&
                (path[step].sliding ? std::abs(actual_ratio - 1.0) < 2.0e-3 && std::abs(reference_ratio - 1.0) < 2.0e-3
                                    : actual_ratio < 0.75 && reference_ratio < 0.75);
            gap.add(actual[node].gap, found->gap);
            pressure.add(actual[node].pressure, found->pressure);
            slip.add(std::abs(path[step].relative_displacement[2]), std::abs(found->slip_1));
            slip.add(std::abs(path[step].relative_displacement[1]), std::abs(found->slip_2));
        }
        for (std::size_t component = 0; component < 3; ++component)
            resultant.add(actual_resultant[component], reference_resultant[component]);
        problem.commit_internal_state(state);
        const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        for (const auto& history : histories) {
            const auto& value = history.cartesian_elastic_tangential_slip;
            biaxial_history =
                biaxial_history && std::abs(value[0]) + std::abs(value[1]) > 0.0 && std::abs(value[2]) > 0.0;
            if (step == 1) sliding_history_before_rotation = value;
            if (step >= 2 && step < path.size() - 1) {
                const std::array<double, 3> expected = rotate(path[step].angle, sliding_history_before_rotation);
                for (std::size_t component = 0; component < 3; ++component)
                    maximum_objective_history_error =
                        std::max(maximum_objective_history_error, std::abs(value[component] - expected[component]));
            }
        }
    }
    fuelsim::test::print_relative_metrics(prefix + "normal_nodal_force", normal);
    fuelsim::test::print_relative_metrics(prefix + "tangential_nodal_force_first", tangent_first);
    fuelsim::test::print_relative_metrics(prefix + "tangential_nodal_force_second", tangent_second);
    fuelsim::test::print_relative_metrics(prefix + "contact_resultant", resultant);
    fuelsim::test::print_relative_metrics(prefix + "gap", gap);
    fuelsim::test::print_relative_metrics(prefix + "pressure", pressure);
    fuelsim::test::print_relative_metrics(prefix + "total_slip", slip);
    std::cout << prefix << "maximum_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << prefix << "maximum_objective_history_error=" << maximum_objective_history_error << '\n'
              << prefix << "physical_state_sequence_matches=" << physical_state_sequence_matches << '\n'
              << prefix << "reported_end_states_match=" << reported_end_states_match << '\n'
              << prefix << "biaxial_history=" << biaxial_history << '\n';
    constexpr double tolerance = 1.0e-2, zero_tolerance = 1.0e-8;
    return check(maximum_coordinate_difference < 3.0e-8,
               prefix + "uses the exact current coordinates from the Abaqus path") &&
           check(physical_state_sequence_matches && reported_end_states_match && biaxial_history,
               prefix + "matches the sticking-sliding-rotation-resticking sequence with biaxial history") &&
           check(maximum_objective_history_error < 1.0e-12,
               prefix + "rotates the nonzero two-component elastic-slip history objectively") &&
           check(fuelsim::test::relative_metrics_below(normal, tolerance) &&
                     normal.maximum_zero_reference_difference < zero_tolerance,
               prefix + "normal nodal-force metrics and zero references agree with Abaqus below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(tangent_first, tolerance) &&
                     tangent_first.maximum_zero_reference_difference < zero_tolerance &&
                     fuelsim::test::relative_metrics_below(tangent_second, tolerance) &&
                     tangent_second.maximum_zero_reference_difference < zero_tolerance,
               prefix +
                   "both tangent-plane nodal-force metrics and zero references agree with Abaqus below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(resultant, tolerance) &&
                     fuelsim::test::relative_metrics_below(gap, tolerance) &&
                     fuelsim::test::relative_metrics_below(pressure, tolerance) &&
                     fuelsim::test::relative_metrics_below(slip, tolerance),
               prefix + "resultant, gap, pressure, and total-slip metrics agree with Abaqus below 1 percent");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 4 && std::string(argv[1]) == "--generate") {
            const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(argv[2]);
            write_abaqus_input(std::string(argv[3]) + "/b42_hex8_sts_friction_objectivity.inp", mesh);
            return 0;
        }
        if (argc != 3) {
            std::cerr << "Usage: fuelsim_b42_hex8_friction_objectivity_abaqus_tests <mesh.e> <reference.csv>\n"
                         "   or: fuelsim_b42_hex8_friction_objectivity_abaqus_tests --generate <mesh.e> "
                         "<abaqus-directory>\n";
            return 2;
        }
        const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(argv[1]);
        const std::vector<std::vector<ReferenceNode>> reference = read_reference(argv[2]);
        const bool small = compare_case(mesh, reference, fuelsim::StrainFormulation::small, "b42_small_");
        const bool finite = compare_case(mesh, reference, fuelsim::StrainFormulation::finite, "b42_finite_");
        if (small && finite) std::cout << "[PASS] B4.2 HEX8 friction-objectivity Abaqus comparison\n";
        return small && finite ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.2 raised: " << error.what() << '\n';
        return 1;
    }
}
