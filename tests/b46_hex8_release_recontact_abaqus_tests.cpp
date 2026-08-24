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
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::array<fuelsim::Field, 3> displacement_fields = {
    fuelsim::Field::displacement_x, fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};

struct PathStep final {
    const char* name;
    std::array<double, 3> displacement;
    bool projected, active;
};

constexpr std::array<PathStep, 6> path = {
    {{"CLOSE", {-0.01, 0.05, 0.04}, true, true}, {"SLIDE", {-0.01, 0.95, 0.20}, true, true},
        {"OPEN", {0.02, 0.95, 0.20}, true, false}, {"OPEN_CROSS", {0.02, 0.05, 0.40}, true, false},
        {"RECONTACT", {-0.01, 0.05, 0.40}, true, true}, {"SLIDE_OUT", {-0.01, 2.20, 0.40}, false, false}}};

struct NodeReference final {
    std::size_t step, id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> displacement{}, reaction{};
};

struct ContactReference final {
    std::size_t step, id;
    fuelsim::CartesianPoint3 current;
    std::array<double, 3> normal_force{};
    double gap, pressure;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::Hex8Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map, double x0, double x1, double y0, double y1, double z0,
    double z1) {
    const std::array<fuelsim::CartesianPoint3, 8> points = {{{x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
        {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}}};
    fuelsim::Hex8Element element{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const std::array<double, 3> key = {points[local].x, points[local].y, points[local].z};
        const auto inserted = node_map.emplace(key, nodes.size());
        if (inserted.second) nodes.push_back(points[local]);
        element.nodes[local] = inserted.first->second;
    }
    return element;
}

fuelsim::UnstructuredHex8Mesh generate_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    const fuelsim::Hex8Element primary_lower = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 0.0, 1.0, -1.0, 2.0),
                               primary_upper = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 1.0, 2.0, -1.0, 2.0),
                               secondary = append_cuboid(nodes, secondary_nodes, 1.0, 2.0, 0.1, 0.9, 0.1, 0.9);
    std::vector<std::size_t> primary_all, secondary_all;
    for (const auto& entry : primary_nodes) primary_all.push_back(entry.second);
    for (const auto& entry : secondary_nodes) secondary_all.push_back(entry.second);
    std::sort(primary_all.begin(), primary_all.end());
    std::sort(secondary_all.begin(), secondary_all.end());
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), {primary_lower, primary_upper, secondary}, {1, 1, 2},
        {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all", primary_all}, {20, "secondary_all", secondary_all},
            {30, "secondary_contact_nodes", {12, 15, 16, 19}}},
        {{40, "primary_contact", {{0, 1}, {1, 1}}}, {50, "secondary_contact", {{2, 3}}}});
}

void write_labels(std::ofstream& output, const std::vector<std::size_t>& nodes) {
    for (std::size_t index = 0; index < nodes.size(); ++index)
        output << nodes[index] + 1 << (index + 1 == nodes.size() ? "\n" : ", ");
}

void write_input(const std::string& path_name, const fuelsim::UnstructuredHex8Mesh& mesh) {
    std::ofstream output(path_name);
    if (!output) throw std::runtime_error("Could not write B4.6 Abaqus input: " + path_name);
    output << std::setprecision(16) << "*Heading\n** B4.6 C3D8 release, open crossing, recontact, and slide-out.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node)
        output << node + 1 << ", " << mesh.nodes()[node].x << ", " << mesh.nodes()[node].y << ", "
               << mesh.nodes()[node].z << '\n';
    for (std::int64_t block = 1; block <= 2; ++block) {
        output << "*Element, type=C3D8, elset=" << (block == 1 ? "PRIMARY" : "SECONDARY") << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != block) continue;
            output << element + 1;
            for (const std::size_t node : mesh.elements()[element].nodes) output << ", " << node + 1;
            output << '\n';
        }
    }
    for (const auto& set : mesh.node_sets()) {
        output << "*Nset, nset=" << set.name << '\n';
        write_labels(output, set.nodes);
    }
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.25\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n*Surface Behavior, pressure-overclosure=LINEAR\n1.e5,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n";
    for (std::size_t step = 0; step < path.size(); ++step) {
        output << "*Step, name=" << path[step].name << ", nlgeom=NO, inc=100\n"
               << "*Static\n0.1, 1., 1.e-8, 0.1\n"
               << (step == 0 ? "*Boundary\nprimary_all, 1, 3, 0.\n" : "*Boundary, op=MOD\n");
        for (std::size_t component = 0; component < 3; ++component)
            output << "secondary_all, " << component + 1 << ", " << component + 1 << ", "
                   << path[step].displacement[component] << '\n';
        output << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
    }
}

void generate(const std::string& directory) {
    const auto mesh = generate_mesh();
    fuelsim::write_exodus_hex8(directory + "/b46_hex8_release_recontact_mesh.e", mesh);
    write_input(directory + "/b46_hex8_release_recontact.inp", mesh);
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column) { return std::stod(values.at(column)); }

std::size_t index_value(const std::vector<std::string>& values, std::size_t column) {
    const double value = number(values, column);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        value != std::floor(value))
        throw std::invalid_argument("Invalid B4.6 CSV index");
    return static_cast<std::size_t>(value);
}

std::vector<std::vector<NodeReference>> read_nodes(const std::string& path_name) {
    std::ifstream input(path_name);
    if (!input) throw std::runtime_error("Could not read B4.6 node reference: " + path_name);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,ux,uy,uz,rfx,rfy,rfz")
        throw std::invalid_argument("Invalid B4.6 node CSV");
    std::vector<std::vector<NodeReference>> result(path.size());
    while (std::getline(input, line)) {
        const auto values = split(line);
        const std::size_t step = index_value(values, 0);
        result.at(step - 1).push_back(
            {step - 1, index_value(values, 1), {number(values, 2), number(values, 3), number(values, 4)},
                {number(values, 5), number(values, 6), number(values, 7)},
                {number(values, 8), number(values, 9), number(values, 10)}});
    }
    return result;
}

std::vector<std::vector<ContactReference>> read_contact(const std::string& path_name) {
    std::ifstream input(path_name);
    if (!input) throw std::runtime_error("Could not read B4.6 contact reference: " + path_name);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,normal_x,normal_y,normal_z,gap,pressure")
        throw std::invalid_argument("Invalid B4.6 contact CSV");
    std::vector<std::vector<ContactReference>> result(path.size());
    while (std::getline(input, line)) {
        const auto values = split(line);
        const std::size_t step = index_value(values, 0);
        result.at(step - 1).push_back(
            {step - 1, index_value(values, 1), {number(values, 2), number(values, 3), number(values, 4)},
                {number(values, 5), number(values, 6), number(values, 7)}, number(values, 8), number(values, 9)});
    }
    return result;
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    result.regions = {{"primary", "primary", material(), 0.0, 300.0, -1, "", fuelsim::StrainFormulation::small},
        {"secondary", "secondary", material(), 0.0, 300.0, -1, "", fuelsim::StrainFormulation::small}};
    fuelsim::ContactDefinition contact;
    contact.name = "release_recontact";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.mechanical = true;
    contact.penalty = 1.0e5;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    result.contacts.push_back(contact);
    return result;
}

std::vector<double> state_for_step(const fuelsim::SteadyProblem& problem, std::size_t step) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> state = problem.initial_state();
    for (std::size_t local = 0; local < spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = spatial.global_node(1, local);
        for (std::size_t component = 0; component < 3; ++component)
            state[spatial.dof(displacement_fields[component], global)] = path[step].displacement[component];
    }
    return state;
}

std::vector<double> contact_residual(
    const fuelsim::cartesian::SpatialAssembly& spatial, const std::vector<double>& state, double* jacobian_error) {
    std::vector<double> global(state.size());
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(
            contribution, local, nullptr, nullptr, 0.0, residual, jacobian_error == nullptr ? nullptr : &jacobian);
        for (std::size_t row = 0; row < dofs.size(); ++row) global[dofs[row]] += residual[row];
        if (jacobian_error == nullptr) continue;
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double perturbation = 1.0e-8;
        for (std::size_t column = 0; column < local.size(); ++column) {
            direction[column] = std::sin(static_cast<double>(column + 1));
            plus[column] += perturbation * direction[column];
            minus[column] -= perturbation * direction[column];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                analytic += jacobian[row * local.size() + column] * direction[column];
            const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += (analytic - reference) * (analytic - reference);
            reference_squared += reference * reference;
        }
        if (reference_squared > 0.0)
            *jacobian_error = std::max(*jacobian_error, std::sqrt(difference_squared / reference_squared));
    }
    return global;
}

void print_metric(const std::string& name, const fuelsim::test::FieldErrorMetrics& metric) {
    if (metric.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metric);
    else {
        fuelsim::test::print_absolute_metrics(name, metric);
        std::cout << name << "_zero_reference_count=" << metric.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << metric.maximum_zero_reference_difference
                  << '\n';
    }
}

bool compare(const fuelsim::UnstructuredHex8Mesh& mesh, const std::string& node_path, const std::string& contact_path) {
    const auto nodes = read_nodes(node_path);
    const auto contacts = read_contact(contact_path);
    fuelsim::SteadyProblem problem(definition(), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto sources = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    std::map<std::size_t, std::size_t> source_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
            source_global[spatial.region_mesh(region).source_node_ids()[local]] = spatial.global_node(region, local);
    fuelsim::test::FieldErrorMetrics displacement[3], reaction[3], normal[3], gap, pressure, resultant[3];
    double coordinate_error = 0.0, current_coordinate_error = 0.0, jacobian_error = 0.0, maximum_action_reaction = 0.0;
    bool states_match = true;
    for (std::size_t step = 0; step + 1 < path.size(); ++step) {
        std::vector<double> state = state_for_step(problem, step);
        problem.validate_state(state);
        const auto actual = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        const auto residual = contact_residual(spatial, state, path[step].active ? &jacobian_error : nullptr);
        std::array<double, 3> balance{}, actual_resultant{}, reference_resultant{};
        for (const auto& reference : nodes[step]) {
            const auto& point = mesh.nodes().at(reference.id);
            coordinate_error = std::max({coordinate_error, std::abs(point.x - reference.point.x),
                std::abs(point.y - reference.point.y), std::abs(point.z - reference.point.z)});
            const std::size_t global = source_global.at(reference.id);
            for (std::size_t component = 0; component < 3; ++component) {
                const double expected_displacement = reference.id >= 12 ? path[step].displacement[component] : 0.0;
                displacement[component].add(expected_displacement, reference.displacement[component]);
                reaction[component].add(
                    residual[spatial.dof(displacement_fields[component], global)], reference.reaction[component]);
                balance[component] += residual[spatial.dof(displacement_fields[component], global)];
            }
        }
        for (std::size_t node = 0; node < actual.size(); ++node) {
            const auto found = std::find_if(contacts[step].begin(), contacts[step].end(),
                [&](const ContactReference& value) { return value.id == sources[node]; });
            if (found == contacts[step].end()) throw std::invalid_argument("B4.6 contact mapping is incomplete");
            const auto& initial = mesh.nodes().at(sources[node]);
            current_coordinate_error =
                std::max({current_coordinate_error, std::abs(initial.x + path[step].displacement[0] - found->current.x),
                    std::abs(initial.y + path[step].displacement[1] - found->current.y),
                    std::abs(initial.z + path[step].displacement[2] - found->current.z)});
            const std::size_t expected_face = step == 0 || step >= 3 ? 0 : 1;
            states_match = states_match && actual[node].projected == path[step].projected &&
                           (actual[node].pressure > 0.0) == path[step].active &&
                           actual[node].primary_face == expected_face;
            for (std::size_t component = 0; component < 3; ++component) {
                const double actual_force = -actual[node].normal_contact_force[component];
                normal[component].add(actual_force, found->normal_force[component]);
                actual_resultant[component] += actual_force;
                reference_resultant[component] += found->normal_force[component];
            }
            gap.add(actual[node].gap, found->gap);
            pressure.add(actual[node].pressure, found->pressure);
        }
        for (std::size_t component = 0; component < 3; ++component) {
            maximum_action_reaction = std::max(maximum_action_reaction, std::abs(balance[component]));
            resultant[component].add(actual_resultant[component], reference_resultant[component]);
        }
        problem.commit_internal_state(state);
    }
    const auto& final_nodes = nodes.back();
    const auto& final_contact = contacts.back();
    bool abaqus_released = true, abaqus_no_projection_sentinel = true;
    for (const auto& reference : final_nodes) {
        for (std::size_t component = 0; component < 3; ++component) {
            const double expected_displacement = reference.id >= 12 ? path.back().displacement[component] : 0.0;
            displacement[component].add(expected_displacement, reference.displacement[component]);
            abaqus_released = abaqus_released && std::abs(reference.reaction[component]) < 1.0e-7;
        }
    }
    for (const auto& reference : final_contact) {
        abaqus_released = abaqus_released && std::abs(reference.normal_force[0]) < 1.0e-8 &&
                          std::abs(reference.normal_force[1]) < 1.0e-8 &&
                          std::abs(reference.normal_force[2]) < 1.0e-8 && std::abs(reference.pressure) < 1.0e-8;
        abaqus_no_projection_sentinel = abaqus_no_projection_sentinel && reference.gap < -1.0e30;
    }
    bool fuelsim_rejected = false;
    try {
        problem.validate_state(state_for_step(problem, path.size() - 1));
    } catch (const std::domain_error&) { fuelsim_rejected = true; }
    for (std::size_t component = 0; component < 3; ++component) {
        print_metric("b46_displacement_" + std::to_string(component), displacement[component]);
        print_metric("b46_reaction_" + std::to_string(component), reaction[component]);
        print_metric("b46_normal_force_" + std::to_string(component), normal[component]);
        print_metric("b46_resultant_" + std::to_string(component), resultant[component]);
    }
    print_metric("b46_gap", gap);
    print_metric("b46_pressure", pressure);
    std::cout << "b46_coordinate_error=" << coordinate_error << '\n'
              << "b46_current_coordinate_error=" << current_coordinate_error << '\n'
              << "b46_jacobian_directional_error=" << jacobian_error << '\n'
              << "b46_action_reaction_maximum_absolute=" << maximum_action_reaction << '\n'
              << "b46_abaqus_final_natural_release=" << abaqus_released << '\n'
              << "b46_fuelsim_final_lost_projection_rejected=" << fuelsim_rejected << '\n';
    const auto passes = [](const fuelsim::test::FieldErrorMetrics& metric) {
        return !metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, 1.0e-2);
    };
    bool metrics = passes(gap) && passes(pressure);
    for (std::size_t component = 0; component < 3; ++component)
        metrics = metrics && passes(displacement[component]) && passes(reaction[component]) &&
                  passes(normal[component]) && passes(resultant[component]);
    const bool complete_reference =
        nodes.size() == path.size() && contacts.size() == path.size() &&
        std::all_of(nodes.begin(), nodes.end(), [](const auto& step) { return step.size() == 20; }) &&
        std::all_of(contacts.begin(), contacts.end(), [](const auto& step) { return step.size() == 4; });
    return check(complete_reference, "B4.6 reads all six tracked Abaqus states") &&
           check(coordinate_error < 3.0e-8 && current_coordinate_error < 3.0e-8,
               "B4.6 uses identical reference and current coordinates") &&
           check(states_match, "B4.6 matches contact states and uniquely transfers all points across the internal face "
                               "while closed and open") &&
           check(metrics, "B4.6 matched states pass all three field metrics below 1 percent") &&
           check(jacobian_error < 2.0e-5, "B4.6 active contact Jacobians match centered directional differences") &&
           check(maximum_action_reaction < 1.0e-8, "B4.6 contact residual preserves three-component action-reaction") &&
           check(abaqus_released && abaqus_no_projection_sentinel,
               "B4.6 Abaqus naturally releases the secondary face after it leaves the complete primary surface") &&
           check(
               fuelsim_rejected, "B4.6 Fuelsim deliberately rejects the same activated-contact lost-projection state");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--generate") {
            generate(argv[2]);
            return 0;
        }
        if (argc != 4) {
            std::cerr << "Usage: fuelsim_b46_hex8_release_recontact_abaqus_tests <mesh.e> <nodes.csv> <contact.csv>\n"
                         "   or: fuelsim_b46_hex8_release_recontact_abaqus_tests --generate <abaqus-directory>\n";
            return 2;
        }
        std::cout << std::scientific << std::setprecision(12);
        const bool passed = compare(fuelsim::read_exodus_hex8(argv[1]), argv[2], argv[3]);
        if (passed) std::cout << "[PASS] B4.6 HEX8 release, recontact, and slide-out comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.6 raised: " << error.what() << '\n';
        return 1;
    }
}
