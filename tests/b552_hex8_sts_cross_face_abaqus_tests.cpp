#include "fuelsim/core/steady_problem.hpp"
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
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::size_t facets = 3;

struct ContactReference final {
    double pressure = 0.0;
    std::array<double, 3> normal_force{};
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

std::map<std::size_t, ContactReference> read_contact_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B5.52 Abaqus contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,shear_stress1_pa,shear_stress2_pa,"
                "normal_force1_n,normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,"
                "contact_heat_flux_w")
        throw std::invalid_argument("Unexpected B5.52 Abaqus contact header");
    std::map<std::size_t, ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 17) throw std::invalid_argument("Unexpected B5.52 Abaqus contact column count");
        const std::size_t node = std::stoul(values[0]);
        if (node == 0 || result.count(node - 1) != 0)
            throw std::invalid_argument("Invalid or duplicate B5.52 Abaqus contact node");
        result[node - 1] = {
            std::stod(values[5]), {std::stod(values[10]), std::stod(values[11]), std::stod(values[12])}};
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

fuelsim::UnstructuredHex8Mesh mesh() {
    std::array<double, facets + 1> angles{};
    for (std::size_t angle = 0; angle <= facets; ++angle)
        angles[angle] = static_cast<double>(angle) * pi / (2.0 * static_cast<double>(facets));
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_map, secondary_map;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, secondary_contact;
    for (std::size_t angle = 0; angle < facets; ++angle) {
        elements.push_back(append_annular(nodes, primary_map, 0.8, 1.0, angles[angle], angles[angle + 1]));
        blocks.push_back(1);
        primary_contact.push_back({elements.size() - 1, 1});
    }
    for (std::size_t angle = 0; angle < facets; ++angle) {
        elements.push_back(append_annular(nodes, secondary_map, 1.005, 1.2, angles[angle], angles[angle + 1]));
        blocks.push_back(2);
        secondary_contact.push_back({elements.size() - 1, 3});
    }
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "primary"}, {2, "secondary"}}, {},
        {{50, "primary_contact", primary_contact}, {60, "secondary_contact", secondary_contact}});
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    const fuelsim::ThermoelasticProperties material =
        fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    result.regions = {{"primary", "primary", material, 0.0, 300.0, -1, "", fuelsim::StrainFormulation::finite},
        {"secondary", "secondary", material, 0.0, 300.0, -1, "", fuelsim::StrainFormulation::finite}};
    fuelsim::ContactDefinition contact;
    contact.name = "nonuniform_curved_contact";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.thermal = false;
    contact.mechanical = true;
    contact.penalty = 1.0e7;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    result.contacts.push_back(contact);
    return result;
}

double inward_displacement(double angle, double z) { return 0.014 + 0.003 * std::cos(2.0 * angle) + 0.002 * (z - 0.5); }

std::vector<double> prescribed_state(
    const fuelsim::SteadyProblem& problem, const fuelsim::UnstructuredHex8Mesh& source_mesh) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> result = problem.initial_state();
    for (std::size_t local = 0; local < spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t source = spatial.region_mesh(1).source_node_ids()[local],
                          global = spatial.global_node(1, local);
        const fuelsim::CartesianPoint3& point = source_mesh.nodes()[source];
        const double radius = std::hypot(point.x, point.y), angle = std::atan2(point.y, point.x),
                     inward = inward_displacement(angle, point.z);
        result[spatial.dof(fuelsim::Field::displacement_x, global)] = -inward * point.x / radius;
        result[spatial.dof(fuelsim::Field::displacement_y, global)] = -inward * point.y / radius;
    }
    return result;
}

double mechanical_contact_directional_error(
    fuelsim::SteadyProblem& problem, const std::vector<double>& state, double perturbation) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> direction(state.size()), analytic(state.size()), plus_residual(state.size()),
        minus_residual(state.size());
    for (std::size_t dof = 0; dof < state.size(); ++dof) direction[dof] = std::sin(static_cast<double>(dof + 1));
    const auto assemble = [&](const std::vector<double>& current, std::vector<double>& residual,
                              std::vector<double>* jacobian_action) {
        problem.validate_state(current);
        for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
            if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
                continue;
            std::vector<std::size_t> dofs;
            problem.contribution_dofs(contribution, dofs);
            std::vector<double> local_state(dofs.size());
            for (std::size_t local = 0; local < dofs.size(); ++local) local_state[local] = current[dofs[local]];
            std::vector<double> local_residual, local_jacobian;
            spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, local_residual,
                jacobian_action == nullptr ? nullptr : &local_jacobian);
            for (std::size_t row = 0; row < dofs.size(); ++row) {
                residual[dofs[row]] += local_residual[row];
                if (jacobian_action != nullptr)
                    for (std::size_t column = 0; column < dofs.size(); ++column)
                        (*jacobian_action)[dofs[row]] +=
                            local_jacobian[row * dofs.size() + column] * direction[dofs[column]];
            }
        }
    };
    std::vector<double> unused(state.size()), plus = state, minus = state;
    assemble(state, unused, &analytic);
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += perturbation * direction[dof];
        minus[dof] -= perturbation * direction[dof];
    }
    assemble(plus, plus_residual, nullptr);
    assemble(minus, minus_residual, nullptr);
    problem.validate_state(state);
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
        difference_squared += std::pow(analytic[row] - finite_difference, 2);
        reference_squared += finite_difference * finite_difference;
    }
    return std::sqrt(difference_squared / reference_squared);
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests "
                         "<Abaqus reference directory>\n";
            return 2;
        }
        const fuelsim::UnstructuredHex8Mesh source_mesh = mesh();
        fuelsim::SteadyProblem problem(definition(), source_mesh);
        std::vector<double> state = prescribed_state(problem, source_mesh);
        problem.validate_state(state);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        const std::vector<fuelsim::CartesianContactNodeSummary> contact =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        const std::vector<std::size_t> sources =
            fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        const std::map<std::size_t, ContactReference> reference =
            read_contact_reference(std::string(argv[1]) + "/b552_hex8_sts_cross_face_contact.csv");
        if (contact.size() != sources.size() || reference.size() != sources.size())
            throw std::invalid_argument("B5.52 Abaqus and fuelsim contact-node counts differ");
        fuelsim::test::FieldErrorMetrics pressure;
        std::array<fuelsim::test::FieldErrorMetrics, 3> normal_force;
        for (std::size_t constraint = 0; constraint < contact.size(); ++constraint) {
            const auto found = reference.find(sources[constraint]);
            if (found == reference.end()) throw std::invalid_argument("B5.52 Abaqus contact node is missing");
            pressure.add(contact[constraint].pressure, found->second.pressure);
            for (std::size_t component = 0; component < 3; ++component)
                normal_force[component].add(
                    -contact[constraint].normal_contact_force[component], found->second.normal_force[component]);
        }
        fuelsim::test::print_relative_metrics("b552_contact_pressure", pressure);
        for (std::size_t component = 0; component < 3; ++component)
            if (normal_force[component].has_relative_norm())
                fuelsim::test::print_relative_metrics(
                    "b552_contact_normal_force_" + std::to_string(component + 1), normal_force[component]);
        const fuelsim::InterfaceSummary interface =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, state);
        const auto partition = spatial.finite_region_partition_summary(0, state);
        const double jacobian_error = mechanical_contact_directional_error(problem, state, 1.0e-8);
        std::cout << std::scientific << std::setprecision(12) << "b552_dofs=" << problem.dof_count() << '\n'
                  << "b552_elements=" << source_mesh.elements().size() << '\n'
                  << "b552_active_contact_nodes=" << interface.active_contact_nodes << '\n'
                  << "b552_cross_face_constraints=" << partition.cross_face_constraint_count << '\n'
                  << "b552_mechanical_contact_jacobian_directional_relative_error=" << jacobian_error << '\n';
        bool passed = check(problem.dof_count() == 128 && source_mesh.elements().size() == 6,
            "B5.52 remains a lightweight six-element, 128-degree-of-freedom operator regression");
        passed = check(problem.jacobian_sparsity_is_state_dependent(),
                     "B5.52 finite-sliding averaged contact declares its changing active Jacobian support") &&
                 passed;
        passed = check(interface.active_contact_nodes == 8 && contact.size() == 8 && partition.constraint_count == 8 &&
                           partition.active_primary_face_count == 3 && partition.cross_face_constraint_count == 4 &&
                           partition.maximum_owners_per_integration_point == 1 && partition.all_projected,
                     "B5.52 keeps all eight constraints active and uniquely partitions four cross-face regions") &&
                 passed;
        passed = check(fuelsim::test::relative_metrics_below(pressure, 1.0e-2) && pressure.zero_reference_count == 0 &&
                           pressure.maximum_zero_reference_difference == 0.0,
                     "B5.52 contact pressure passes relative L2, relative absolute peak, and maximum pointwise "
                     "relative error below one percent without a denominator floor") &&
                 passed;
        passed = check(fuelsim::test::relative_metrics_below(normal_force[0], 1.0e-2) &&
                           fuelsim::test::relative_metrics_below(normal_force[1], 1.0e-2) &&
                           fuelsim::test::relative_metrics_below(normal_force[2], 1.0e-2),
                     "B5.52 Cartesian normal nodal forces agree with Abaqus below one percent") &&
                 passed;
        passed = check(jacobian_error < 2.0e-5,
                     "B5.52 finite-region contact Jacobian matches a centered directional difference") &&
                 passed;
        if (passed) std::cout << "[PASS] B5.52 lightweight cross-primary-face STS Abaqus verification\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
