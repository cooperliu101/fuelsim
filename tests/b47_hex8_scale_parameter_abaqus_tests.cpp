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
#include <utility>
#include <vector>

namespace {
constexpr std::array<fuelsim::Field, 3> displacement_fields = {fuelsim::Field::displacement_x,
    fuelsim::Field::displacement_y,
    fuelsim::Field::displacement_z};
constexpr std::array<double, 3> slip_fractions = {0.25, 0.75, 1.5};
constexpr std::array<const char*, 3> step_names = {"STICK_LOW", "STICK_HIGH", "SLIDE"};

struct CaseSpec final {
    std::string name;
    double y_size, z_size, penalty, friction, slip_tolerance, primary_modulus, secondary_modulus;
};

struct NodeReference final {
    std::size_t step, id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> displacement{}, reaction{};
};

struct ContactReference final {
    std::size_t step, id;
    std::array<double, 3> normal_force{}, tangential_force{};
    double slip_1, slip_2, gap, pressure;
};

std::vector<CaseSpec> specs() {
    return {{"unit", 1.0, 1.0, 1.0e8, 0.30, 1.0e-5, 1.0e9, 1.0e9},
        {"aspect16", 4.0, 0.25, 2.5e7, 0.55, 2.0e-5, 1.0e8, 1.0e9},
        {"area4", 4.0, 1.0, 5.0e8, 0.20, 5.0e-6, 1.0e9, 1.0e8},
        {"area_quarter", 0.25, 1.0, 1.0e9, 0.40, 4.0e-5, 1.0e10, 1.0e8}};
}

double characteristic_length(const CaseSpec& spec) {
    return std::sqrt(spec.y_size * spec.z_size);
}

double elastic_slip(const CaseSpec& spec) {
    return spec.slip_tolerance * characteristic_length(spec);
}

double normal_closure(const CaseSpec& spec) {
    return -1.0e-4 * characteristic_length(spec);
}

std::array<double, 3> displacement(const CaseSpec& spec, std::size_t step) {
    const double magnitude = slip_fractions[step] * elastic_slip(spec);
    return {normal_closure(spec), 0.6 * magnitude, 0.8 * magnitude};
}

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::UnstructuredHex8Mesh mesh(const CaseSpec& spec) {
    const std::vector<fuelsim::CartesianPoint3> nodes = {{0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, spec.y_size, 0.0},
        {0.0, spec.y_size, 0.0},
        {0.0, 0.0, spec.z_size},
        {1.0, 0.0, spec.z_size},
        {1.0, spec.y_size, spec.z_size},
        {0.0, spec.y_size, spec.z_size},
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {2.0, spec.y_size, 0.0},
        {1.0, spec.y_size, 0.0},
        {1.0, 0.0, spec.z_size},
        {2.0, 0.0, spec.z_size},
        {2.0, spec.y_size, spec.z_size},
        {1.0, spec.y_size, spec.z_size}};
    const std::vector<fuelsim::Hex8Element> elements = {{{{0, 1, 2, 3, 4, 5, 6, 7}}},
        {{{8, 9, 10, 11, 12, 13, 14, 15}}}};
    return fuelsim::UnstructuredHex8Mesh(nodes,
        elements,
        {1, 2},
        {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all", {0, 1, 2, 3, 4, 5, 6, 7}},
            {20, "secondary_all", {8, 9, 10, 11, 12, 13, 14, 15}},
            {30, "secondary_contact_nodes", {8, 11, 12, 15}}},
        {{40, "primary_contact", {{0, 1}}}, {50, "secondary_contact", {{1, 3}}}});
}

void labels(std::ofstream& output, const std::vector<std::size_t>& values) {
    for (std::size_t index = 0; index < values.size(); ++index)
        output << values[index] + 1 << (index + 1 == values.size() ? "\n" : ", ");
}

void write_input(const std::string& path_name, const CaseSpec& spec, const fuelsim::UnstructuredHex8Mesh& generated) {
    std::ofstream output(path_name);
    if (!output)
        throw std::runtime_error("Could not write B4.7 Abaqus input: " + path_name);
    output << std::setprecision(16) << "*Heading\n** B4.7 " << spec.name
           << " C3D8 scale and friction-parameter probe.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n*Node\n";
    for (std::size_t node = 0; node < generated.nodes().size(); ++node)
        output << node + 1 << ", " << generated.nodes()[node].x << ", " << generated.nodes()[node].y << ", "
               << generated.nodes()[node].z << '\n';
    output << "*Element, type=C3D8, elset=PRIMARY\n1";
    for (const auto node : generated.elements()[0].nodes)
        output << ", " << node + 1;
    output << "\n*Element, type=C3D8, elset=SECONDARY\n2";
    for (const auto node : generated.elements()[1].nodes)
        output << ", " << node + 1;
    output << "\n*Nset, nset=PRIMARY_ALL\n";
    labels(output, generated.node_sets()[0].nodes);
    output << "*Nset, nset=SECONDARY_ALL\n";
    labels(output, generated.node_sets()[1].nodes);
    output << "*Nset, nset=SECONDARY_CONTACT_NODES\n";
    labels(output, generated.node_sets()[2].nodes);
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=PRIMARY_MATERIAL\n*Elastic\n"
           << spec.primary_modulus << ", 0.0\n*Material, name=SECONDARY_MATERIAL\n*Elastic\n"
           << spec.secondary_modulus << ", 0.0\n"
           << "*Solid Section, elset=PRIMARY, material=PRIMARY_MATERIAL\n,\n"
           << "*Solid Section, elset=SECONDARY, material=SECONDARY_MATERIAL\n,\n"
           << "*Surface Interaction, name=CONTACT\n*Surface Behavior, penalty=LINEAR\n"
           << spec.penalty << ",\n*Friction, slip tolerance=" << spec.slip_tolerance << '\n'
           << spec.friction
           << ",\n*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n";
    for (std::size_t step = 0; step < step_names.size(); ++step) {
        const auto value = displacement(spec, step);
        output << "*Step, name=" << step_names[step] << ", nlgeom=NO, inc=100\n"
               << "*Static\n1., 1., 1.e-8, 1.\n*Boundary\nPRIMARY_ALL, 1, 3, 0.\n";
        for (std::size_t component = 0; component < 3; ++component)
            output << "SECONDARY_ALL, " << component + 1 << ", " << component + 1 << ", " << value[component] << '\n';
        output << "*Output, field, frequency=1\n*Node Output\nRF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
    }
}

void generate(const std::string& directory) {
    for (const auto& spec : specs()) {
        const auto generated = mesh(spec);
        const std::string base = directory + "/b47_hex8_scale_" + spec.name;
        fuelsim::write_exodus_hex8(base + "_mesh.e", generated);
        write_input(base + ".inp", spec, generated);
    }
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column) {
    return std::stod(values.at(column));
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column) {
    const double value = number(values, column);
    if (value < 0.0 || value != std::floor(value)
        || value > static_cast<double>(std::numeric_limits<std::size_t>::max()))
        throw std::invalid_argument("Invalid B4.7 CSV index");
    return static_cast<std::size_t>(value);
}

std::vector<std::vector<NodeReference>> read_nodes(const std::string& path_name) {
    std::ifstream input(path_name);
    if (!input)
        throw std::runtime_error("Could not read B4.7 node reference: " + path_name);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,ux,uy,uz,rfx,rfy,rfz")
        throw std::invalid_argument("Invalid B4.7 node CSV");
    std::vector<std::vector<NodeReference>> result(step_names.size());
    while (std::getline(input, line)) {
        const auto values = split(line);
        const std::size_t step = index_value(values, 0);
        result.at(step - 1).push_back({step - 1,
            index_value(values, 1),
            {number(values, 2), number(values, 3), number(values, 4)},
            {number(values, 5), number(values, 6), number(values, 7)},
            {number(values, 8), number(values, 9), number(values, 10)}});
    }
    return result;
}

std::vector<std::vector<ContactReference>> read_contact(const std::string& path_name) {
    std::ifstream input(path_name);
    if (!input)
        throw std::runtime_error("Could not read B4.7 contact reference: " + path_name);
    std::string line;
    if (!std::getline(input, line)
        || line != "step,id,normal_x,normal_y,normal_z,tangent_x,tangent_y,tangent_z,slip1,slip2,gap,pressure")
        throw std::invalid_argument("Invalid B4.7 contact CSV");
    std::vector<std::vector<ContactReference>> result(step_names.size());
    while (std::getline(input, line)) {
        const auto values = split(line);
        const std::size_t step = index_value(values, 0);
        result.at(step - 1).push_back({step - 1,
            index_value(values, 1),
            {number(values, 2), number(values, 3), number(values, 4)},
            {number(values, 5), number(values, 6), number(values, 7)},
            number(values, 8),
            number(values, 9),
            number(values, 10),
            number(values, 11)});
    }
    return result;
}

fuelsim::ThermoelasticProperties material(double modulus) {
    return fuelsim::test::thermoelastic(0.0, 10.0, modulus, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
}

fuelsim::SpatialDefinition definition(const CaseSpec& spec) {
    fuelsim::SpatialDefinition result;
    result.regions = {
        {"primary", "primary", material(spec.primary_modulus), 0.0, 300.0, -1, "", fuelsim::StrainFormulation::small},
        {"secondary",
            "secondary",
            material(spec.secondary_modulus),
            0.0,
            300.0,
            -1,
            "",
            fuelsim::StrainFormulation::small}};
    fuelsim::ContactDefinition contact;
    contact.name = "scale_parameter";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.mechanical = true;
    contact.penalty = spec.penalty;
    contact.friction_coefficient = spec.friction;
    contact.friction_slip_tolerance = spec.slip_tolerance;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::small;
    result.contacts.push_back(contact);
    return result;
}

std::vector<double> state_for_step(const fuelsim::SteadyProblem& problem, const CaseSpec& spec, std::size_t step) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> state = problem.initial_state();
    const auto value = displacement(spec, step);
    for (std::size_t local = 0; local < spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = spatial.global_node(1, local);
        for (std::size_t component = 0; component < 3; ++component)
            state[spatial.dof(displacement_fields[component], global)] = value[component];
    }
    return state;
}

std::vector<double> contact_residual(const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state,
    double* jacobian_error) {
    std::vector<double> global(state.size());
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        for (std::size_t row = 0; row < dofs.size(); ++row)
            global[dofs[row]] += residual[row];
        if (jacobian_error == nullptr)
            continue;
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

bool compare_case(const CaseSpec& spec,
    const std::string& mesh_path,
    const std::string& node_path,
    const std::string& contact_path) {
    const auto generated = fuelsim::read_exodus_hex8(mesh_path);
    const auto nodes = read_nodes(node_path);
    const auto contacts = read_contact(contact_path);
    fuelsim::SteadyProblem problem(definition(spec), generated);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto sources = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    std::map<std::size_t, std::size_t> source_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
            source_global[spatial.region_mesh(region).source_node_ids()[local]] = spatial.global_node(region, local);
    fuelsim::test::FieldErrorMetrics displacement_error[3], reaction_error[3], normal_error, tangential_y, tangential_z,
        slip_y, slip_z, gap_error, pressure_error, resultant_error[3], identified_elastic_slip;
    double coordinate_error = 0.0, jacobian_error = 0.0, balance_error = 0.0;
    bool states_match = true;
    for (std::size_t step = 0; step < step_names.size(); ++step) {
        const auto expected_displacement = displacement(spec, step);
        std::vector<double> state = state_for_step(problem, spec, step);
        problem.validate_state(state);
        const auto actual = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        const auto residual = contact_residual(spatial, state, &jacobian_error);
        std::array<double, 3> actual_resultant{}, reference_resultant{}, balance{};
        for (const auto& reference : nodes[step]) {
            const auto& point = generated.nodes().at(reference.id);
            coordinate_error = std::max({coordinate_error,
                std::abs(point.x - reference.point.x),
                std::abs(point.y - reference.point.y),
                std::abs(point.z - reference.point.z)});
            const std::size_t global = source_global.at(reference.id);
            for (std::size_t component = 0; component < 3; ++component) {
                displacement_error[component].add(reference.id >= 8 ? expected_displacement[component] : 0.0,
                    reference.displacement[component]);
                reaction_error[component].add(residual[spatial.dof(displacement_fields[component], global)],
                    reference.reaction[component]);
                balance[component] += residual[spatial.dof(displacement_fields[component], global)];
            }
        }
        for (std::size_t node = 0; node < actual.size(); ++node) {
            const auto found = std::find_if(contacts[step].begin(),
                contacts[step].end(),
                [&](const ContactReference& value) { return value.id == sources[node]; });
            if (found == contacts[step].end())
                throw std::invalid_argument("B4.7 contact mapping is incomplete");
            normal_error.add(-actual[node].normal_contact_force[0], found->normal_force[0]);
            tangential_y.add(-actual[node].tangential_contact_force[1], found->tangential_force[1]);
            tangential_z.add(-actual[node].tangential_contact_force[2], found->tangential_force[2]);
            slip_y.add(actual[node].tangential_slip[1], -found->slip_2);
            slip_z.add(actual[node].tangential_slip[2], found->slip_1);
            gap_error.add(actual[node].gap, found->gap);
            pressure_error.add(actual[node].pressure, found->pressure);
            for (std::size_t component = 0; component < 3; ++component) {
                const double actual_force =
                    -actual[node].normal_contact_force[component] - actual[node].tangential_contact_force[component];
                const double reference_force = found->normal_force[component] + found->tangential_force[component];
                actual_resultant[component] += actual_force;
                reference_resultant[component] += reference_force;
            }
            states_match = states_match && actual[node].projected && actual[node].pressure > 0.0
                           && actual[node].sliding == (step == 2);
            if (step == 0) {
                const double normal_magnitude = std::hypot(found->normal_force[0],
                                 std::hypot(found->normal_force[1], found->normal_force[2])),
                             tangent_magnitude = std::hypot(found->tangential_force[0],
                                 std::hypot(found->tangential_force[1], found->tangential_force[2])),
                             slip_magnitude = std::hypot(found->slip_1, found->slip_2);
                identified_elastic_slip.add(spec.friction * normal_magnitude * slip_magnitude / tangent_magnitude,
                    elastic_slip(spec));
            }
        }
        for (std::size_t component = 0; component < 3; ++component) {
            resultant_error[component].add(actual_resultant[component], reference_resultant[component]);
            balance_error = std::max(balance_error, std::abs(balance[component]));
        }
        problem.commit_internal_state(state);
    }
    const std::string prefix = "b47_" + spec.name + "_";
    for (std::size_t component = 0; component < 3; ++component) {
        print_metric(prefix + "displacement_" + std::to_string(component), displacement_error[component]);
        print_metric(prefix + "reaction_" + std::to_string(component), reaction_error[component]);
        print_metric(prefix + "resultant_" + std::to_string(component), resultant_error[component]);
    }
    print_metric(prefix + "normal_force", normal_error);
    print_metric(prefix + "tangential_force_y", tangential_y);
    print_metric(prefix + "tangential_force_z", tangential_z);
    print_metric(prefix + "slip_y", slip_y);
    print_metric(prefix + "slip_z", slip_z);
    print_metric(prefix + "gap", gap_error);
    print_metric(prefix + "pressure", pressure_error);
    print_metric(prefix + "identified_elastic_slip", identified_elastic_slip);
    std::cout << prefix << "characteristic_length=" << characteristic_length(spec) << '\n'
              << prefix << "absolute_elastic_slip=" << elastic_slip(spec) << '\n'
              << prefix << "jacobian_directional_error=" << jacobian_error << '\n'
              << prefix << "action_reaction_maximum_absolute=" << balance_error << '\n';
    const auto passes = [](const fuelsim::test::FieldErrorMetrics& metric) {
        return !metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, 1.0e-2);
    };
    bool metrics = passes(normal_error) && passes(tangential_y) && passes(tangential_z) && passes(slip_y)
                   && passes(slip_z) && passes(gap_error) && passes(pressure_error) && passes(identified_elastic_slip);
    for (std::size_t component = 0; component < 3; ++component)
        metrics = metrics && passes(displacement_error[component]) && passes(reaction_error[component])
                  && passes(resultant_error[component]);
    return check(nodes.size() == 3 && contacts.size() == 3 && nodes[0].size() == 16 && contacts[0].size() == 4,
               "B4.7 " + spec.name + " reads every tracked node and contact constraint")
           && check(coordinate_error < 1.0e-12, "B4.7 " + spec.name + " uses the tracked Exodus geometry")
           && check(states_match, "B4.7 " + spec.name + " has two sticking states followed by full sliding")
           && check(metrics, "B4.7 " + spec.name + " fields and identified elastic slip pass 1 percent")
           && check(jacobian_error < 1.0e-6,
               "B4.7 " + spec.name + " sticking and sliding Jacobians match centered differences")
           && check(balance_error < 1.0e-8, "B4.7 " + spec.name + " preserves three-component action-reaction");
}

bool check_default_slip_tolerance() {
    CaseSpec spec = specs().front();
    spec.slip_tolerance = 5.0e-3;
    const fuelsim::UnstructuredHex8Mesh generated = mesh(spec);
    fuelsim::SpatialDefinition explicit_definition = definition(spec);
    fuelsim::SpatialDefinition default_definition = explicit_definition;
    default_definition.contacts[0].friction_slip_tolerance = 0.0;
    fuelsim::SteadyProblem explicit_problem(std::move(explicit_definition), generated);
    fuelsim::SteadyProblem default_problem(std::move(default_definition), generated);
    const std::vector<double> explicit_state = state_for_step(explicit_problem, spec, 1);
    const std::vector<double> default_state = state_for_step(default_problem, spec, 1);
    const std::vector<double> explicit_residual =
        contact_residual(fuelsim::cartesian::ProblemAccess::view(explicit_problem), explicit_state, nullptr);
    const std::vector<double> default_residual =
        contact_residual(fuelsim::cartesian::ProblemAccess::view(default_problem), default_state, nullptr);
    return check(explicit_state == default_state && explicit_residual == default_residual,
        "B4.7 omitted surface-friction slip_tolerance equals the Abaqus default 0.005");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--generate") {
            generate(argv[2]);
            return 0;
        }
        constexpr int arguments_per_case = 4;
        const auto cases = specs();
        if (argc != 1 + arguments_per_case * static_cast<int>(cases.size())) {
            std::cerr << "Usage: fuelsim_b47_hex8_scale_parameter_abaqus_tests "
                         "<name> <mesh.e> <nodes.csv> <contact.csv> [... ]\n";
            return 2;
        }
        std::cout << std::scientific << std::setprecision(12);
        bool passed = true;
        for (std::size_t index = 0; index < cases.size(); ++index) {
            const int argument = 1 + static_cast<int>(index) * arguments_per_case;
            if (argv[argument] != cases[index].name)
                throw std::invalid_argument("B4.7 case order differs");
            passed = compare_case(cases[index], argv[argument + 1], argv[argument + 2], argv[argument + 3]) && passed;
        }
        const double minimum_ratio = 0.1, maximum_ratio = 100.0;
        passed = check(cases[2].secondary_modulus / cases[2].primary_modulus == minimum_ratio
                           && cases[3].primary_modulus / cases[3].secondary_modulus == maximum_ratio,
                     "B4.7 decks span material-stiffness ratios from 0.1 to 100")
                 && passed;
        passed = check_default_slip_tolerance() && passed;
        if (passed)
            std::cout << "[PASS] B4.7 HEX8 scale and contact-parameter matrix\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.7 raised: " << error.what() << '\n';
        return 1;
    }
}
