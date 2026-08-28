#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

struct NodeReference final {
    double temperature = 0.0, reaction_heat_flux = 0.0;
    std::array<double, 3> displacement{}, reaction_force{};
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

std::vector<NodeReference> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B5.22 Abaqus reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,temperature_k,reaction_heat_flux_w,u1_m,u2_m,u3_m,rf1_n,rf2_n,rf3_n")
        throw std::invalid_argument("Unexpected B5.22 Abaqus header");
    std::map<std::size_t, NodeReference> rows;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 9) throw std::invalid_argument("Unexpected B5.22 Abaqus column count");
        const std::size_t node = std::stoul(values[0]);
        if (node < 1 || rows.count(node) != 0) throw std::invalid_argument("Invalid or duplicate B5.22 Abaqus node");
        rows[node] = {std::stod(values[1]), std::stod(values[2]),
            {std::stod(values[3]), std::stod(values[4]), std::stod(values[5])},
            {std::stod(values[6]), std::stod(values[7]), std::stod(values[8])}};
    }
    std::vector<NodeReference> result;
    result.reserve(rows.size());
    for (std::size_t node = 1; node <= rows.size(); ++node) {
        const auto found = rows.find(node);
        if (found == rows.end()) throw std::invalid_argument("B5.22 Abaqus reference is incomplete");
        result.push_back(found->second);
    }
    return result;
}

fuelsim::CartesianPoint3 cylindrical(double radius, double angle, double z) {
    return {radius * std::cos(angle), radius * std::sin(angle), z};
}

fuelsim::Hex8Element append_annular(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map, double r0, double r1, double a0, double a1) {
    const std::array<std::array<double, 3>, 8> logical = {{{r0, a0, 0.0}, {r1, a0, 0.0}, {r1, a1, 0.0}, {r0, a1, 0.0},
        {r0, a0, 1.0}, {r1, a0, 1.0}, {r1, a1, 1.0}, {r0, a1, 1.0}}};
    fuelsim::Hex8Element element{};
    for (std::size_t local = 0; local < logical.size(); ++local) {
        const auto inserted = node_map.emplace(logical[local], nodes.size());
        if (inserted.second) nodes.push_back(cylindrical(logical[local][0], logical[local][1], logical[local][2]));
        element.nodes[local] = inserted.first->second;
    }
    return element;
}

fuelsim::UnstructuredHex8Mesh mesh(std::size_t facets, double primary_radius = 1.0, double secondary_radius = 1.005) {
    std::vector<double> angles(facets + 1);
    for (std::size_t angle = 0; angle <= facets; ++angle)
        angles[angle] = static_cast<double>(angle) * pi / (2.0 * static_cast<double>(facets));
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_map, secondary_map;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, secondary_contact;
    for (std::size_t angle = 0; angle + 1 < angles.size(); ++angle) {
        elements.push_back(
            append_annular(nodes, primary_map, primary_radius - 0.2, primary_radius, angles[angle], angles[angle + 1]));
        blocks.push_back(1);
        primary_contact.push_back({elements.size() - 1, 1});
    }
    for (std::size_t angle = 0; angle + 1 < angles.size(); ++angle) {
        elements.push_back(append_annular(
            nodes, secondary_map, secondary_radius, secondary_radius + 0.195, angles[angle], angles[angle + 1]));
        blocks.push_back(2);
        secondary_contact.push_back({elements.size() - 1, 3});
    }
    std::vector<std::size_t> primary_all, secondary_all;
    for (const auto& entry : primary_map) primary_all.push_back(entry.second);
    for (const auto& entry : secondary_map) secondary_all.push_back(entry.second);
    std::sort(primary_all.begin(), primary_all.end());
    std::sort(secondary_all.begin(), secondary_all.end());
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "primary"}, {2, "secondary"}}, {{10, "primary_all", primary_all}, {20, "secondary_all", secondary_all}},
        {{50, "primary_contact", primary_contact}, {60, "secondary_contact", secondary_contact}});
}

fuelsim::SpatialDefinition definition(bool swap = false) {
    fuelsim::SpatialDefinition result;
    const fuelsim::ThermoelasticProperties material =
        fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    result.regions = {{"primary", "primary", material, 0.0, 300.0, -1, "", fuelsim::StrainFormulation::finite},
        {"secondary", "secondary", material, 0.0, 400.0, -1, "", fuelsim::StrainFormulation::finite}};
    fuelsim::ContactDefinition contact;
    contact.name = "faceted_thermal_contact";
    contact.primary = swap ? "secondary_contact" : "primary_contact";
    contact.secondary = swap ? "primary_contact" : "secondary_contact";
    contact.thermal = true;
    contact.mechanical = true;
    contact.gap_conductivity = 1.0;
    contact.minimum_gap = 1.0;
    contact.penalty = 1.0e7;
    contact.friction_coefficient = 0.1;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    contact.gap_heat_conductance_law = fuelsim::GapHeatConductanceLaw::affine;
    contact.gap_conductance = 0.0;
    contact.gap_conductance_pressure_derivative = 1.0e-3;
    result.contacts.push_back(contact);
    return result;
}

std::vector<double> state(const fuelsim::SteadyProblem& problem, const fuelsim::UnstructuredHex8Mesh& source_mesh) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> result = problem.initial_state();
    for (std::size_t local = 0; local < spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t source = spatial.region_mesh(1).source_node_ids()[local],
                          global = spatial.global_node(1, local);
        const fuelsim::CartesianPoint3& point = source_mesh.nodes()[source];
        const double radius = std::hypot(point.x, point.y);
        result[spatial.dof(fuelsim::Field::temperature, global)] = 400.0;
        result[spatial.dof(fuelsim::Field::displacement_x, global)] = -0.015 * point.x / radius;
        result[spatial.dof(fuelsim::Field::displacement_y, global)] = -0.015 * point.y / radius;
    }
    return result;
}

std::vector<double> contact_residual(const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state, fuelsim::SpatialContributionType type, double& conservation_error) {
    std::vector<double> global(state.size());
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != type) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, nullptr);
        for (std::size_t row = 0; row < dofs.size(); ++row) global[dofs[row]] += residual[row];
        if (type == fuelsim::SpatialContributionType::thermal_contact) {
            double balance = 0.0;
            for (std::size_t row = 0; row < 8; ++row) balance += residual[row];
            conservation_error = std::max(conservation_error, std::abs(balance));
        }
    }
    return global;
}

struct Comparison final {
    bool passed = false;
    double resultant_force_error = 0.0, pressure_weighted_force_error = 0.0, heat_rate_error = 0.0,
           normal_integration_error = 0.0;
};

Comparison compare(const std::vector<NodeReference>& reference, const std::string& name, std::size_t facets,
    double primary_radius = 1.0, double secondary_radius = 1.005, bool swap = false) {
    const std::size_t primary_node_count = 4 * (facets + 1);
    const fuelsim::UnstructuredHex8Mesh source_mesh = mesh(facets, primary_radius, secondary_radius);
    if (reference.size() != source_mesh.nodes().size())
        throw std::invalid_argument("B5.22 Abaqus reference and Fuelsim mesh node counts differ");
    fuelsim::SteadyProblem problem(definition(swap), source_mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> current = state(problem, source_mesh);
    problem.validate_state(current);
    double conservation_error = 0.0;
    const std::vector<double> thermal = contact_residual(spatial, current,
                                  fuelsim::SpatialContributionType::thermal_contact, conservation_error),
                              mechanical = contact_residual(spatial, current,
                                  fuelsim::SpatialContributionType::mechanical_contact, conservation_error);
    std::map<std::size_t, std::size_t> source_to_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
            source_to_global[spatial.region_mesh(region).source_node_ids()[local]] = spatial.global_node(region, local);
    fuelsim::test::FieldErrorMetrics nodal_heat;
    double state_difference = 0.0, actual_secondary_heat = 0.0, reference_secondary_heat = 0.0;
    std::array<double, 3> actual_secondary_force{}, reference_secondary_force{};
    const std::array<fuelsim::Field, 3> fields = {
        fuelsim::Field::displacement_x, fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
    for (std::size_t source = 0; source < reference.size(); ++source) {
        const std::size_t global = source_to_global.at(source);
        nodal_heat.add(thermal[spatial.dof(fuelsim::Field::temperature, global)], reference[source].reaction_heat_flux);
        state_difference = std::max(state_difference,
            std::abs(current[spatial.dof(fuelsim::Field::temperature, global)] - reference[source].temperature));
        for (std::size_t component = 0; component < 3; ++component)
            state_difference = std::max(state_difference,
                std::abs(current[spatial.dof(fields[component], global)] - reference[source].displacement[component]));
        if (source >= primary_node_count) {
            actual_secondary_heat += thermal[spatial.dof(fuelsim::Field::temperature, global)];
            reference_secondary_heat += reference[source].reaction_heat_flux;
            for (std::size_t component = 0; component < 3; ++component) {
                actual_secondary_force[component] -= mechanical[spatial.dof(fields[component], global)];
                reference_secondary_force[component] -= reference[source].reaction_force[component];
            }
        }
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, current);
    const std::vector<fuelsim::CartesianContactNodeSummary> contact =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, current);
    std::set<std::size_t> primary_faces;
    for (const auto& point : contact) primary_faces.insert(point.primary_face);
    fuelsim::test::print_relative_metrics("b522_" + name + "_nodal_reaction_heat_flux", nodal_heat);
    const double heat_rate_error =
        std::abs(actual_secondary_heat - reference_secondary_heat) / std::abs(reference_secondary_heat);
    double force_difference_squared = 0.0, reference_force_squared = 0.0;
    for (std::size_t component = 0; component < 3; ++component) {
        force_difference_squared +=
            std::pow(actual_secondary_force[component] - reference_secondary_force[component], 2);
        reference_force_squared += reference_secondary_force[component] * reference_secondary_force[component];
    }
    const double resultant_force_error = std::sqrt(force_difference_squared / reference_force_squared);
    const double pressure_weighted_force_reference = std::abs(reference_secondary_heat) / (1.0e-3 * 100.0);
    const double pressure_weighted_force_error =
        std::abs(interface.total_contact_force - pressure_weighted_force_reference) /
        std::abs(pressure_weighted_force_reference);
    const double angle = pi / (2.0 * static_cast<double>(facets));
    const double penetration = (primary_radius - (secondary_radius - 0.015)) * std::cos(0.5 * angle);
    const double integration_radius = swap ? primary_radius : secondary_radius - 0.015;
    const double facet_force = 1.0e7 * penetration * 2.0 * integration_radius * std::sin(0.5 * angle);
    std::array<double, 3> analytical_force{};
    for (std::size_t facet = 0; facet < facets; ++facet) {
        const double normal_angle = (static_cast<double>(facet) + 0.5) * angle;
        analytical_force[0] += facet_force * std::cos(normal_angle);
        analytical_force[1] += facet_force * std::sin(normal_angle);
    }
    double normal_difference_squared = 0.0, analytical_squared = 0.0;
    for (std::size_t component = 0; component < 3; ++component) {
        normal_difference_squared += std::pow(actual_secondary_force[component] - analytical_force[component], 2);
        analytical_squared += analytical_force[component] * analytical_force[component];
    }
    const double normal_integration_error = std::sqrt(normal_difference_squared / analytical_squared);
    std::cout << "b522_" << name << "_nodal_heat_zero_reference_count=" << nodal_heat.zero_reference_count << '\n'
              << "b522_" << name << "_nodal_heat_zero_reference_maximum_absolute_difference="
              << nodal_heat.maximum_zero_reference_difference << '\n'
              << "b522_" << name << "_total_heat_rate_relative_error=" << heat_rate_error << '\n'
              << "b522_" << name << "_resultant_contact_force_relative_l2=" << resultant_force_error << '\n'
              << "b522_" << name << "_pressure_weighted_contact_force_relative_error=" << pressure_weighted_force_error
              << '\n'
              << "b522_" << name << "_normal_integration_relative_error=" << normal_integration_error << '\n'
              << "b522_" << name << "_actual_secondary_force=" << actual_secondary_force[0] << ','
              << actual_secondary_force[1] << ',' << actual_secondary_force[2] << '\n'
              << "b522_" << name << "_reference_secondary_force=" << reference_secondary_force[0] << ','
              << reference_secondary_force[1] << ',' << reference_secondary_force[2] << '\n'
              << "b522_" << name << "_analytical_secondary_force=" << analytical_force[0] << ',' << analytical_force[1]
              << ',' << analytical_force[2] << '\n'
              << "b522_" << name << "_state_maximum_absolute_difference=" << state_difference << '\n'
              << "b522_" << name << "_thermal_conservation_maximum_absolute=" << conservation_error << '\n'
              << "b522_" << name << "_active_primary_face_count=" << primary_faces.size() << '\n';
    const double thermal_tolerance = swap ? 7.0e-3 : 1.0e-6;
    const bool passed =
        check(state_difference < 3.0e-8, "B5.22 " + name + " Fuelsim and Abaqus use the same faceted-cylinder state") &&
        check(contact.size() == 2 * (facets + 1) && primary_faces.size() == facets &&
                  interface.active_contact_nodes == contact.size(),
            "B5.22 " + name + " activates every curved-surface contact node and primary facet") &&
        check(
            fuelsim::test::relative_metrics_below(nodal_heat, thermal_tolerance) && heat_rate_error < thermal_tolerance,
            "B5.22 " + name +
                " nodal and total pressure-dependent thermal-contact rates satisfy the designation-specific limit") &&
        check(pressure_weighted_force_error < thermal_tolerance,
            "B5.22 " + name +
                " pressure-weighted scalar contact force satisfies the designation-specific Abaqus heat-rate limit") &&
        check(normal_integration_error < 1.0e-12,
            "B5.22 " + name + " mechanical resultant matches analytical constant-facet normal integration") &&
        check(nodal_heat.maximum_zero_reference_difference < 1.0e-10 && conservation_error < 1.0e-10 &&
                  std::abs(std::abs(interface.total_heat_rate) - std::abs(actual_secondary_heat)) < 1.0e-10,
            "B5.22 " + name + " curved thermal contact remains exactly conservative and has exact zero support");
    return {passed, resultant_force_error, pressure_weighted_force_error, heat_rate_error, normal_integration_error};
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests "
                     "<facets3.csv> <facets6.csv> <facets12.csv> <swapped.csv> <tight.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const Comparison facets3 = compare(read_reference(argv[1]), "facets3", 3);
        const Comparison facets6 = compare(read_reference(argv[2]), "facets6", 6);
        const Comparison facets12 = compare(read_reference(argv[3]), "facets12", 12);
        const Comparison swapped = compare(read_reference(argv[4]), "swapped", 3, 1.0, 1.005, true);
        const Comparison tight = compare(read_reference(argv[5]), "tight", 3, 0.5, 0.505);
        const double first_contraction = facets6.resultant_force_error / facets3.resultant_force_error,
                     second_contraction = facets12.resultant_force_error / facets6.resultant_force_error;
        std::cout << "b522_facets3_to_facets6_force_error_contraction=" << first_contraction << '\n'
                  << "b522_facets6_to_facets12_force_error_contraction=" << second_contraction << '\n';
        const bool passed =
            facets3.passed && facets6.passed && facets12.passed && swapped.passed && tight.passed &&
            check(first_contraction < 0.5 && second_contraction < 0.5,
                "B5.22 Abaqus mechanical resultant difference contracts under each curved-facet refinement") &&
            check(swapped.resultant_force_error < 2.0e-2 && tight.resultant_force_error < 2.0e-2,
                "B5.22 designation exchange and tighter curvature keep the mechanical resultant difference below two "
                "percent");
        if (passed) std::cout << "[PASS] B5.22 Abaqus C3D8T faceted thermal contact\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.22 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
