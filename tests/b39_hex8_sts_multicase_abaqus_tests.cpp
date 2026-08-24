#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
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
struct CaseSpec final {
    std::string name;
    std::vector<double> primary_y, primary_z, secondary_y, secondary_z;
    double angle, gap, penalty, primary_modulus, secondary_modulus, closure;
};

struct GeneratedCase final {
    fuelsim::UnstructuredHex8Mesh mesh;
    std::vector<std::size_t> primary_outer_nodes, secondary_outer_nodes, secondary_contact_nodes;
};

struct DisplacementReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> normal;
    double displacement;
};

struct ForceReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double force;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<CaseSpec> case_specs() {
    return {{"matching", {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, 0.0, 0.0, 1.0e11, 1.0e9, 1.0e9, -1.0e-5},
        {"nonmatching", {0.0, 0.5, 1.0}, {0.0, 0.5, 1.0}, {0.0, 0.28, 0.67, 1.0}, {0.0, 0.43, 1.0}, 0.0, 0.0, 1.0e11,
            1.0e9, 1.0e9, -1.0e-5},
        {"tilted_gap", {0.0, 0.33, 0.7, 1.0}, {0.0, 0.38, 1.0}, {0.0, 0.48, 1.0}, {0.0, 0.22, 0.63, 1.0},
            0.4363323129985824, 2.0e-6, 5.0e10, 1.4e9, 8.0e8, -1.2e-5}};
}

fuelsim::CartesianPoint3 transform(double normal_coordinate, double tangent_coordinate, double z, double angle) {
    const double cosine = std::cos(angle), sine = std::sin(angle);
    return {normal_coordinate * cosine - tangent_coordinate * sine,
        normal_coordinate * sine + tangent_coordinate * cosine, z};
}

fuelsim::Hex8Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map, double normal_lower, double normal_upper,
    double tangent_lower, double tangent_upper, double z_lower, double z_upper, double angle) {
    const std::array<std::array<double, 3>, 8> local_points = {
        {{normal_lower, tangent_lower, z_lower}, {normal_upper, tangent_lower, z_lower},
            {normal_upper, tangent_upper, z_lower}, {normal_lower, tangent_upper, z_lower},
            {normal_lower, tangent_lower, z_upper}, {normal_upper, tangent_lower, z_upper},
            {normal_upper, tangent_upper, z_upper}, {normal_lower, tangent_upper, z_upper}}};
    fuelsim::Hex8Element element{};
    for (std::size_t local = 0; local < local_points.size(); ++local) {
        const auto inserted = node_map.emplace(local_points[local], nodes.size());
        if (inserted.second)
            nodes.push_back(transform(local_points[local][0], local_points[local][1], local_points[local][2], angle));
        element.nodes[local] = inserted.first->second;
    }
    return element;
}

GeneratedCase generate_case(const CaseSpec& spec) {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, primary_outer, secondary_contact, secondary_outer;
    for (std::size_t z = 0; z + 1 < spec.primary_z.size(); ++z)
        for (std::size_t y = 0; y + 1 < spec.primary_y.size(); ++y) {
            const std::size_t element = elements.size();
            elements.push_back(append_cuboid(nodes, primary_nodes, 0.0, 1.0, spec.primary_y[y], spec.primary_y[y + 1],
                spec.primary_z[z], spec.primary_z[z + 1], spec.angle));
            blocks.push_back(1);
            primary_contact.push_back({element, 1});
            primary_outer.push_back({element, 3});
        }
    for (std::size_t z = 0; z + 1 < spec.secondary_z.size(); ++z)
        for (std::size_t y = 0; y + 1 < spec.secondary_y.size(); ++y) {
            const std::size_t element = elements.size();
            elements.push_back(
                append_cuboid(nodes, secondary_nodes, 1.0 + spec.gap, 2.0 + spec.gap, spec.secondary_y[y],
                    spec.secondary_y[y + 1], spec.secondary_z[z], spec.secondary_z[z + 1], spec.angle));
            blocks.push_back(2);
            secondary_contact.push_back({element, 3});
            secondary_outer.push_back({element, 1});
        }

    std::vector<std::size_t> primary_outer_nodes, secondary_outer_nodes, secondary_contact_nodes;
    for (const auto& entry : primary_nodes)
        if (entry.first[0] == 0.0) primary_outer_nodes.push_back(entry.second);
    for (const auto& entry : secondary_nodes) {
        if (entry.first[0] == 2.0 + spec.gap) secondary_outer_nodes.push_back(entry.second);
        if (entry.first[0] == 1.0 + spec.gap) secondary_contact_nodes.push_back(entry.second);
    }
    for (std::vector<std::size_t>* values : {&primary_outer_nodes, &secondary_outer_nodes, &secondary_contact_nodes})
        std::sort(values->begin(), values->end());
    std::vector<fuelsim::NodeSet> node_sets = {{10, "primary_outer_nodes", primary_outer_nodes},
        {11, "secondary_outer_nodes", secondary_outer_nodes}, {12, "secondary_contact_nodes", secondary_contact_nodes}};
    GeneratedCase result{
        fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(blocks),
            {{1, "primary"}, {2, "secondary"}}, std::move(node_sets),
            {{20, "primary_contact", std::move(primary_contact)}, {21, "primary_outer", std::move(primary_outer)},
                {30, "secondary_contact", std::move(secondary_contact)},
                {31, "secondary_outer", std::move(secondary_outer)}}),
        std::move(primary_outer_nodes), std::move(secondary_outer_nodes), std::move(secondary_contact_nodes)};
    return result;
}

void write_label_set(std::ofstream& output, const std::vector<std::size_t>& labels) {
    for (std::size_t index = 0; index < labels.size(); ++index) {
        output << labels[index] + 1;
        if ((index + 1) % 16 == 0 || index + 1 == labels.size())
            output << '\n';
        else
            output << ", ";
    }
}

void write_abaqus_input(const std::string& path, const CaseSpec& spec, const GeneratedCase& generated) {
    const fuelsim::UnstructuredHex8Mesh& mesh = generated.mesh;
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B3.9 Abaqus input: " + path);
    output << std::setprecision(16) << "*Heading\n"
           << "** B3.9: " << spec.name << " C3D8 small-sliding surface-to-surface contact.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    const auto write_block = [&](const std::string& name, std::int64_t id) {
        output << "*Element, type=C3D8, elset=" << name << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != id) continue;
            output << element + 1;
            for (const std::size_t node : mesh.elements()[element].nodes) output << ", " << node + 1;
            output << '\n';
        }
    };
    write_block("PRIMARY", 1);
    write_block("SECONDARY", 2);
    output << "*Nset, nset=PRIMARY_OUTER\n";
    write_label_set(output, generated.primary_outer_nodes);
    output << "*Nset, nset=SECONDARY_OUTER\n";
    write_label_set(output, generated.secondary_outer_nodes);
    output << "*Nset, nset=SECONDARY_CONTACT_NODES\n";
    write_label_set(output, generated.secondary_contact_nodes);
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\n"
           << "PRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\n"
           << "SECONDARY, S6\n"
           << "*Material, name=PRIMARY_MATERIAL\n"
           << "*Elastic\n"
           << spec.primary_modulus << ", 0.0\n"
           << "*Material, name=SECONDARY_MATERIAL\n"
           << "*Elastic\n"
           << spec.secondary_modulus << ", 0.0\n"
           << "*Solid Section, elset=PRIMARY, material=PRIMARY_MATERIAL\n,\n"
           << "*Solid Section, elset=SECONDARY, material=SECONDARY_MATERIAL\n,\n"
           << "*Surface Interaction, name=PENALTY_CONTACT\n"
           << "*Surface Behavior, penalty=LINEAR\n"
           << spec.penalty << ",\n"
           << "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
           << "*Step, name=LOAD, nlgeom=NO, inc=100\n"
           << "*Static\n0.1, 1., 1.e-8, 0.1\n"
           << "*Boundary\n"
           << "PRIMARY_OUTER, 1, 3, 0.\n";
    const double cosine = std::cos(spec.angle), sine = std::sin(spec.angle);
    output << "SECONDARY_OUTER, 1, 1, " << spec.closure * cosine << '\n'
           << "SECONDARY_OUTER, 2, 2, " << spec.closure * sine << '\n'
           << "SECONDARY_OUTER, 3, 3, 0.\n"
           << "*Output, field, frequency=1\n"
           << "*Node Output\nCOORD, RF, U\n"
           << "*Contact Output\nCSTRESS, CDISP, CFORCE\n"
           << "*End Step\n";
}

void write_material(std::ofstream& output, const std::string& name, double modulus) {
    output << "  [" << name << "]\n"
           << "    [thermal]\n      function = constant_thermophysical\n      conductivity = 1\n"
           << "      density = 1\n      specific_heat = 1\n    []\n"
           << "    [elasticity]\n      function = constant_isotropic\n      young_modulus = " << modulus
           << "\n      poisson_ratio = 0.0\n    []\n  []\n";
}

void write_fsi(const std::string& path, const std::string& mesh_relative_path, const CaseSpec& spec) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B3.9 Fuelsim input: " + path);
    output << std::setprecision(16) << "[Case]\n  version = 3\n  problem = steady\n"
           << "  geometry = cartesian_3d\n[]\n\n[Mesh]\n  type = exodus\n  file = " << mesh_relative_path
           << "\n[]\n\n[Materials]\n";
    write_material(output, "primary", spec.primary_modulus);
    write_material(output, "secondary", spec.secondary_modulus);
    output << "[]\n\n[Regions]\n"
           << "  [primary]\n    block = primary\n    material = primary\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n"
           << "  [secondary]\n    block = secondary\n    material = secondary\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n[]\n\n"
           << "[Contact]\n  [interface]\n    primary = primary_contact\n    secondary = secondary_contact\n"
           << "    [mechanical]\n      formulation = penalty\n      discretization = surface_to_surface\n"
           << "      sliding = small\n      penalty = " << spec.penalty << "\n      mu = 0\n    []\n  []\n[]\n\n"
           << "[BoundaryConditions]\n"
           << "  [primary_temperature]\n    type = dirichlet\n    boundary = primary_outer\n"
           << "    field = temperature\n    value = 300\n  []\n"
           << "  [secondary_temperature]\n    type = dirichlet\n    boundary = secondary_outer\n"
           << "    field = temperature\n    value = 300\n  []\n";
    for (std::size_t component = 0; component < 3; ++component)
        output << "  [primary_displacement_" << component << "]\n    type = dirichlet\n"
               << "    boundary = primary_outer\n    field = displacement_"
               << (component == 0 ? "x" : (component == 1 ? "y" : "z")) << "\n    value = 0\n  []\n";
    const double cosine = std::cos(spec.angle), sine = std::sin(spec.angle);
    const std::array<double, 3> displacement = {spec.closure * cosine, spec.closure * sine, 0.0};
    for (std::size_t component = 0; component < 3; ++component)
        output << "  [secondary_displacement_" << component << "]\n    type = dirichlet\n"
               << "    boundary = secondary_outer\n    field = displacement_"
               << (component == 0 ? "x" : (component == 1 ? "y" : "z")) << "\n    value = " << displacement[component]
               << "\n    scale_with_load = true\n  []\n";
    output << "[]\n\n[Executioner]\n  type = steady\n  load_steps = 4\n[]\n\n"
           << "[Solver]\n  absolute_tolerance = 1e-8\n  relative_tolerance = 1e-11\n"
           << "  step_tolerance = 1e-12\n  maximum_iterations = 30\n  linear_solver = direct\n"
           << "  direct_factorization = mumps\n  field_residual_scaling = true\n"
           << "  linear_relative_tolerance = 1e-11\n  maximum_linear_iterations = 400\n[]\n\n"
           << "[Outputs]\n  console = true\n[]\n";
}

void generate_references(const std::string& abaqus_directory, const std::string& fuelsim_directory) {
    for (const CaseSpec& spec : case_specs()) {
        const GeneratedCase generated = generate_case(spec);
        const std::string base = "b39_hex8_sts_" + spec.name;
        fuelsim::write_exodus_hex8(abaqus_directory + '/' + base + "_mesh.e", generated.mesh);
        write_abaqus_input(abaqus_directory + '/' + base + ".inp", spec, generated);
        write_fsi(
            fuelsim_directory + "/steady_hex8_sts_" + spec.name + "_abaqus.fsi", "../abaqus/" + base + "_mesh.e", spec);
    }
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete B3.9 CSV row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.9 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_displacement,normal_x,normal_y,normal_z,id,x,y,z")
        throw std::invalid_argument("Unexpected B3.9 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 4, path)),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}, number(values, 0, path)});
    }
    return result;
}

std::vector<ForceReference> read_forces(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.9 nodal-force reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_force,id,x,y,z") throw std::invalid_argument("Unexpected B3.9 force header in " + path);
    std::vector<ForceReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)}, number(values, 0, path)});
    }
    return result;
}

std::array<double, 2> read_resultants(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B3.9 resultant reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_contact_resultant,normal_outer_reaction")
        throw std::invalid_argument("Unexpected B3.9 resultant header in " + path);
    if (!std::getline(input, line)) throw std::invalid_argument("Missing B3.9 resultant row in " + path);
    const std::vector<std::string> values = split_csv(line);
    return {std::abs(number(values, 0, path)), std::abs(number(values, 1, path))};
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = definition.solver.absolute_tolerance;
    result.relative_tolerance = definition.solver.relative_tolerance;
    result.step_tolerance = definition.solver.step_tolerance;
    result.maximum_iterations = definition.solver.maximum_iterations;
    result.linear_solver = definition.solver.linear_solver;
    result.direct_factorization = definition.solver.direct_factorization;
    result.field_residual_scaling = definition.solver.field_residual_scaling;
    result.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    return result;
}

double coordinate_difference(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return std::max({std::abs(first.x - second.x), std::abs(first.y - second.y), std::abs(first.z - second.z)});
}

bool compare_case(const std::string& name, const std::string& case_path, const std::string& displacement_path,
    const std::string& force_path, const std::string& reaction_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(
        solve.completed && solve.solve.converged && solve.completed_steps == definition.steady_execution.load_steps,
        "B3.9 " + name + " Fuelsim solve completes every load step");
    const std::vector<DisplacementReference> displacements = read_displacements(displacement_path);
    const std::vector<ForceReference> forces = read_forces(force_path);
    if (displacements.empty() || forces.empty()) throw std::invalid_argument("B3.9 reference fields are empty");
    const std::array<double, 3> normal = displacements.front().normal;
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    fuelsim::test::FieldErrorMetrics displacement_error, force_error;
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_error = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto reference = std::find_if(displacements.begin(), displacements.end(),
                [source](const DisplacementReference& value) { return value.id == source; });
            if (reference == displacements.end() || present[source])
                throw std::invalid_argument("B3.9 displacement node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_error =
                std::max(maximum_coordinate_error, coordinate_difference(mesh.nodes()[source], reference->point));
            const std::size_t global = spatial.global_node(region, local);
            double actual = 0.0;
            for (std::size_t component = 0; component < 3; ++component)
                actual += normal[component] * solve.solve.state[fields[component + 1].begin + global];
            displacement_error.add(actual, reference->displacement);
        }
    }
    const std::vector<fuelsim::CartesianContactNodeSummary> summaries =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const std::vector<std::size_t> sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    double maximum_force_coordinate_error = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const auto reference = std::find_if(
            forces.begin(), forces.end(), [&](const ForceReference& value) { return value.id == sources[node]; });
        if (reference == forces.end()) throw std::invalid_argument("B3.9 contact-force node mapping is incomplete");
        maximum_force_coordinate_error = std::max(
            maximum_force_coordinate_error, coordinate_difference(mesh.nodes()[reference->id], reference->point));
        force_error.add(summaries[node].contact_force, reference->force);
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    const std::array<double, 2> reference_resultants = read_resultants(reaction_path);
    const double resultant_error =
        std::abs(interface.total_contact_force - reference_resultants[0]) / reference_resultants[0];
    std::size_t mechanical_contributions = 0;
    std::array<double, 3> residual_balance{}, normal_resultant{};
    double maximum_jacobian_directional_error = 0.0, maximum_tangential_force = 0.0;
    for (const fuelsim::CartesianContactNodeSummary& summary : summaries) {
        for (std::size_t component = 0; component < 3; ++component)
            normal_resultant[component] -= summary.normal_contact_force[component];
        maximum_tangential_force = std::max(maximum_tangential_force, summary.tangential_force);
    }
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        ++mechanical_contributions;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = solve.solve.state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double perturbation = 1.0e-10;
        for (std::size_t index = 0; index < local.size(); ++index) {
            direction[index] = std::sin(static_cast<double>(index + 1));
            plus[index] += perturbation * direction[index];
            minus[index] -= perturbation * direction[index];
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
            difference_squared += std::pow(analytic - reference, 2);
            reference_squared += reference * reference;
            const std::size_t global = dofs[row];
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (global >= field.begin && global < field.end) residual_balance[component] += residual[row];
            }
        }
        maximum_jacobian_directional_error =
            std::max(maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
    }
    double normal_direction_error = 0.0;
    const double normal_resultant_magnitude = std::hypot(normal_resultant[0], normal_resultant[1], normal_resultant[2]);
    for (std::size_t component = 0; component < 3; ++component)
        normal_direction_error = std::max(normal_direction_error,
            std::abs(normal_resultant[component] / normal_resultant_magnitude - normal[component]));
    fuelsim::test::print_relative_metrics("b39_" + name + "_normal_displacement", displacement_error);
    fuelsim::test::print_relative_metrics("b39_" + name + "_nodal_normal_force", force_error);
    std::cout << "b39_" << name << "_normal_resultant=" << interface.total_contact_force << '\n'
              << "b39_" << name << "_abaqus_normal_resultant=" << reference_resultants[0] << '\n'
              << "b39_" << name << "_abaqus_outer_normal_reaction=" << reference_resultants[1] << '\n'
              << "b39_" << name << "_normal_resultant_relative_error=" << resultant_error << '\n'
              << "b39_" << name << "_mechanical_constraint_count=" << mechanical_contributions << '\n'
              << "b39_" << name << "_contact_jacobian_directional_error=" << maximum_jacobian_directional_error << '\n'
              << "b39_" << name << "_contact_residual_balance=" << residual_balance[0] << ',' << residual_balance[1]
              << ',' << residual_balance[2] << '\n'
              << "b39_" << name << "_normal_direction_error=" << normal_direction_error << '\n';
    constexpr double tolerance = 1.0e-2;
    return check(displacement_error.value_count == mesh.nodes().size() &&
                     std::all_of(present.begin(), present.end(), [](bool value) { return value; }),
               "B3.9 " + name + " compares every mesh-node normal displacement") &&
           check(force_error.value_count == summaries.size() && forces.size() == summaries.size(),
               "B3.9 " + name + " compares every secondary nodal normal force") &&
           check(maximum_coordinate_error < 1.0e-12 && maximum_force_coordinate_error < 1.0e-12,
               "B3.9 " + name + " Abaqus references use the tracked Exodus coordinates") &&
           check(interface.active_contact_nodes == summaries.size() && interface.unprojected_contact_nodes == 0,
               "B3.9 " + name + " keeps every secondary constraint active and projected") &&
           check(mechanical_contributions == summaries.size(),
               "B3.9 " + name + " stores one averaged constraint for every unique secondary contact node") &&
           check(maximum_jacobian_directional_error < 1.0e-7,
               "B3.9 " + name + " averaged-contact Jacobians match centered directional differences") &&
           check(std::abs(residual_balance[0]) < 1.0e-8 && std::abs(residual_balance[1]) < 1.0e-8 &&
                     std::abs(residual_balance[2]) < 1.0e-8,
               "B3.9 " + name + " contact residuals preserve three-component action-reaction balance") &&
           check(normal_direction_error < 1.0e-12 && maximum_tangential_force < 1.0e-10,
               "B3.9 " + name + " frictionless resultant follows the common reference-surface normal") &&
           check(fuelsim::test::relative_metrics_below(displacement_error, tolerance),
               "B3.9 " + name + " normal-displacement three Abaqus errors are below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(force_error, tolerance),
               "B3.9 " + name + " nodal-normal-force three Abaqus errors are below 1 percent") &&
           check(
               resultant_error < tolerance, "B3.9 " + name + " normal resultant agrees with Abaqus below 1 percent") &&
           passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--generate") {
        try {
            generate_references(argv[2], argv[3]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] B3.9 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    constexpr int arguments_per_case = 5;
    if (argc < 1 + arguments_per_case || (argc - 1) % arguments_per_case != 0) {
        std::cerr << "Usage: fuelsim_b39_hex8_sts_multicase_abaqus_tests "
                     "<name> <case.fsi> <displacement.csv> <force.csv> <reaction.csv> [...]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B3.9 HEX8 Abaqus surface contact comparison\n");
        bool passed = true;
        for (int argument = 1; argument < argc; argument += arguments_per_case)
            passed = compare_case(argv[argument], argv[argument + 1], argv[argument + 2], argv[argument + 3],
                         argv[argument + 4]) &&
                     passed;
        if (passed && session.rank() == 0) std::cout << "[PASS] B3.9 HEX8 multi-aspect Abaqus STS comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B3.9 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
