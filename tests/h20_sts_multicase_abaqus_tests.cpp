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
    double angle, gap, penalty, primary_modulus, secondary_modulus, closure, side_traction;
};

struct GeneratedCase final {
    fuelsim::UnstructuredHex20Mesh mesh;
    std::vector<std::pair<std::string, double>> secondary_normal_displacements;
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
    return {
        {"biaxial", {0.0, 0.5, 1.0}, {0.0, 0.5, 1.0}, {0.0, 0.28, 0.67, 1.0}, {0.0, 0.43, 1.0}, 0.0, 0.0, 1.0e11, 1.0e9,
            1.0e9, -1.0e-5, 0.0},
        {"tilted", {0.0, 0.22, 0.58, 1.0}, {0.0, 0.45, 1.0}, {0.0, 0.4, 1.0}, {0.0, 0.25, 0.72, 1.0},
            0.4363323129985824, 0.0, 1.0e11, 1.0e9, 1.0e9, -1.0e-5, 0.0},
        {"graded", {0.0, 0.33, 0.7, 1.0}, {0.0, 0.38, 1.0}, {0.0, 0.48, 1.0}, {0.0, 0.22, 0.63, 1.0}, 0.0, 0.0, 1.0e11,
            2.0e9, 6.5e8, -1.0e-5, -1.0e3},
        {"soft_penalty", {0.0, 0.5, 1.0}, {0.0, 0.5, 1.0}, {0.0, 0.28, 0.67, 1.0}, {0.0, 0.43, 1.0}, 0.0, 0.0, 1.0e10,
            1.0e9, 1.0e9, -1.0e-5, 0.0},
        {"initial_gap", {0.0, 0.3, 0.64, 1.0}, {0.0, 0.41, 1.0}, {0.0, 0.46, 1.0}, {0.0, 0.2, 0.61, 1.0}, 0.0, 2.0e-6,
            5.0e10, 1.4e9, 8.0e8, -1.2e-5, 0.0},
    };
}

fuelsim::CartesianPoint3 transform(double normal_coordinate, double tangent_coordinate, double z, double angle) {
    const double cosine = std::cos(angle), sine = std::sin(angle);
    return {normal_coordinate * cosine - tangent_coordinate * sine,
        normal_coordinate * sine + tangent_coordinate * cosine, z};
}

fuelsim::Hex20Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map, double normal_lower, double normal_upper,
    double tangent_lower, double tangent_upper, double z_lower, double z_upper, double angle) {
    const std::array<std::array<double, 3>, 8> corners = {
        {{normal_lower, tangent_lower, z_lower}, {normal_upper, tangent_lower, z_lower},
            {normal_upper, tangent_upper, z_lower}, {normal_lower, tangent_upper, z_lower},
            {normal_lower, tangent_lower, z_upper}, {normal_upper, tangent_lower, z_upper},
            {normal_upper, tangent_upper, z_upper}, {normal_lower, tangent_upper, z_upper}}};
    const std::array<std::pair<std::size_t, std::size_t>, 12> edges = {
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}}};
    std::array<std::array<double, 3>, 20> local_points{};
    std::copy(corners.begin(), corners.end(), local_points.begin());
    for (std::size_t edge = 0; edge < edges.size(); ++edge)
        for (std::size_t component = 0; component < 3; ++component)
            local_points[8 + edge][component] =
                0.5 * (corners[edges[edge].first][component] + corners[edges[edge].second][component]);
    fuelsim::Hex20Element element{};
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
    std::vector<fuelsim::Hex20Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, primary_outer, secondary_contact, secondary_outer, secondary_top;
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
            if (z + 2 == spec.secondary_z.size()) secondary_top.push_back({element, 5});
        }

    std::vector<std::size_t> primary_outer_nodes, secondary_outer_nodes, secondary_contact_nodes;
    for (const auto& entry : primary_nodes)
        if (entry.first[0] == 0.0) primary_outer_nodes.push_back(entry.second);
    for (const auto& entry : secondary_nodes) {
        if (entry.first[0] == 2.0 + spec.gap) secondary_outer_nodes.push_back(entry.second);
        if (entry.first[0] == 1.0 + spec.gap) secondary_contact_nodes.push_back(entry.second);
    }
    std::sort(primary_outer_nodes.begin(), primary_outer_nodes.end());
    std::sort(secondary_outer_nodes.begin(), secondary_outer_nodes.end());
    std::sort(secondary_contact_nodes.begin(), secondary_contact_nodes.end());
    std::vector<fuelsim::NodeSet> node_sets = {{10, "primary_outer", primary_outer_nodes},
        {11, "secondary_outer", secondary_outer_nodes}, {12, "secondary_contact_nodes", secondary_contact_nodes}};
    std::vector<std::pair<std::string, double>> load_sets = {{"secondary_outer", spec.closure}};

    return {
        fuelsim::UnstructuredHex20Mesh(std::move(nodes), std::move(elements), std::move(blocks),
            {{1, "primary"}, {2, "secondary"}}, std::move(node_sets),
            {{20, "primary_contact", std::move(primary_contact)}, {21, "primary_outer", std::move(primary_outer)},
                {30, "secondary_contact", std::move(secondary_contact)},
                {31, "secondary_outer", std::move(secondary_outer)}, {32, "secondary_top", std::move(secondary_top)}}),
        std::move(load_sets)};
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
    const fuelsim::UnstructuredHex20Mesh& mesh = generated.mesh;
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write H20.29 Abaqus input: " + path);
    output << std::setprecision(16) << "*Heading\n"
           << "** H20.29: " << spec.name << " C3D20 small-sliding surface-to-surface contact.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    const std::array<std::size_t, 20> abaqus_order = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
    const auto write_block = [&](const std::string& name, std::int64_t id) {
        output << "*Element, type=C3D20, elset=" << name << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != id) continue;
            output << element + 1;
            for (const std::size_t local : abaqus_order) output << ", " << mesh.elements()[element].nodes[local] + 1;
            output << '\n';
        }
    };
    write_block("PRIMARY", 1);
    write_block("SECONDARY", 2);
    for (const fuelsim::NodeSet& set : mesh.node_sets()) {
        output << "*Nset, nset=" << set.name << '\n';
        write_label_set(output, set.nodes);
    }
    std::vector<std::size_t> secondary_top_elements;
    for (const fuelsim::ElementSide& side : mesh.side_set("secondary_top").sides)
        secondary_top_elements.push_back(side.element);
    output << "*Elset, elset=SECONDARY_TOP_ELEMENTS\n";
    write_label_set(output, secondary_top_elements);
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\n"
           << "PRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\n"
           << "SECONDARY, S6\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_TOP\n"
           << "SECONDARY_TOP_ELEMENTS, S2\n"
           << "*Material, name=PRIMARY_MATERIAL\n"
           << "*Elastic\n"
           << spec.primary_modulus << ", 0.25\n"
           << "*Material, name=SECONDARY_MATERIAL\n"
           << "*Elastic\n"
           << spec.secondary_modulus << ", 0.25\n"
           << "*Solid Section, elset=PRIMARY, material=PRIMARY_MATERIAL\n,\n"
           << "*Solid Section, elset=SECONDARY, material=SECONDARY_MATERIAL\n,\n"
           << "*Surface Interaction, name=PENALTY_CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n"
           << spec.penalty << ",\n"
           << "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
           << "*Step, name=LOAD, nlgeom=NO, inc=100\n"
           << "*Static\n0.1, 1., 1.e-8, 0.1\n"
           << "*Boundary\n"
           << "primary_outer, 1, 3, 0.\n";
    const double cosine = std::cos(spec.angle), sine = std::sin(spec.angle);
    for (const auto& load : generated.secondary_normal_displacements) {
        output << load.first << ", 1, 1, " << load.second * cosine << '\n'
               << load.first << ", 2, 2, " << load.second * sine << '\n'
               << load.first << ", 3, 3, 0.\n";
    }
    if (spec.side_traction != 0.0)
        output << "*Dsload\nSECONDARY_TOP, TRVEC, " << spec.side_traction << ", 1., 0., 0.\n";
    output << "*Output, field, frequency=1\n"
           << "*Node Output\nCOORD, RF, U\n"
           << "*Contact Output\nCSTRESS, CDISP, CFORCE\n"
           << "*End Step\n";
}

void write_material(std::ofstream& output, const std::string& name, double modulus) {
    output << "  [" << name << "]\n"
           << "    [thermal]\n"
           << "      function = constant_thermophysical\n"
           << "      conductivity = 10\n"
           << "      density = 1\n"
           << "      specific_heat = 1\n"
           << "    []\n"
           << "    [elasticity]\n"
           << "      function = constant_isotropic\n"
           << "      young_modulus = " << modulus << '\n'
           << "      poisson_ratio = 0.25\n"
           << "    []\n"
           << "  []\n";
}

void write_fsi(const std::string& path, const std::string& mesh_relative_path, const CaseSpec& spec,
    const GeneratedCase& generated) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write H20.29 Fuelsim input: " + path);
    output << std::setprecision(16) << "[Case]\n  version = 3\n  problem = steady\n  geometry = cartesian_3d\n[]\n\n"
           << "[Mesh]\n  type = exodus\n  file = " << mesh_relative_path << "\n[]\n\n"
           << "[Materials]\n";
    write_material(output, "primary", spec.primary_modulus);
    write_material(output, "secondary", spec.secondary_modulus);
    output << "[]\n\n"
           << "[Regions]\n"
           << "  [primary]\n    block = primary\n    material = primary\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n"
           << "  [secondary]\n    block = secondary\n    material = secondary\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n[]\n\n"
           << "[Contact]\n  [interface]\n    primary = primary_contact\n    secondary = secondary_contact\n"
           << "    [mechanical]\n      formulation = penalty\n      discretization = surface_to_surface\n"
           << "      penalty = " << spec.penalty << "\n      mu = 0\n    []\n  []\n[]\n\n"
           << "[BoundaryConditions]\n"
           << "  [primary_temperature]\n    type = dirichlet\n    boundary = primary_outer\n"
           << "    field = temperature\n    value = 300\n  []\n"
           << "  [secondary_temperature]\n    type = dirichlet\n    boundary = secondary_outer\n"
           << "    field = temperature\n    value = 300\n  []\n";
    for (std::size_t component = 0; component < 3; ++component) {
        output << "  [primary_displacement_" << component << "]\n    type = dirichlet\n"
               << "    boundary = primary_outer\n    field = displacement_"
               << (component == 0 ? "x" : (component == 1 ? "y" : "z")) << "\n    value = 0\n  []\n";
    }
    const double cosine = std::cos(spec.angle), sine = std::sin(spec.angle);
    for (std::size_t load = 0; load < generated.secondary_normal_displacements.size(); ++load) {
        const std::string& boundary = generated.secondary_normal_displacements[load].first;
        const double normal_displacement = generated.secondary_normal_displacements[load].second;
        const std::array<double, 3> values = {normal_displacement * cosine, normal_displacement * sine, 0.0};
        for (std::size_t component = 0; component < 3; ++component) {
            output << "  [secondary_displacement_" << load << '_' << component
                   << "]\n    type = dirichlet\n    boundary = " << boundary << "\n    field = displacement_"
                   << (component == 0 ? "x" : (component == 1 ? "y" : "z")) << "\n    value = " << values[component]
                   << "\n    scale_with_load = true\n  []\n";
        }
    }
    if (spec.side_traction != 0.0)
        output << "  [secondary_top_traction]\n    type = traction\n    boundary = secondary_top\n"
               << "    field = displacement_x\n    value = " << spec.side_traction
               << "\n    configuration = reference\n    scale_with_load = true\n  []\n";
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
        const std::string base = "h20_29_hex20_sts_" + spec.name;
        fuelsim::write_exodus_hex20(abaqus_directory + '/' + base + "_mesh.e", generated.mesh);
        write_abaqus_input(abaqus_directory + '/' + base + ".inp", spec, generated);
        write_fsi(fuelsim_directory + "/steady_hex20_sts_" + spec.name + "_abaqus.fsi", "../abaqus/" + base + "_mesh.e",
            spec, generated);
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
    if (index >= values.size()) throw std::invalid_argument("Incomplete H20.29 CSV row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.29 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_displacement,normal_x,normal_y,normal_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.29 displacement header in " + path);
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
    if (!input) throw std::runtime_error("Could not read H20.29 nodal-force reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_force,id,x,y,z") throw std::invalid_argument("Unexpected H20.29 force header in " + path);
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
    if (!input) throw std::runtime_error("Could not read H20.29 resultant reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_contact_resultant,normal_outer_reaction")
        throw std::invalid_argument("Unexpected H20.29 resultant header in " + path);
    if (!std::getline(input, line)) throw std::invalid_argument("Missing H20.29 reaction value in " + path);
    const std::vector<std::string> values = split_csv(line);
    return {std::abs(number(values, 0, path)), std::abs(number(values, 1, path))};
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.direct_factorization = definition.solver.direct_factorization;
    options.field_residual_scaling = definition.solver.field_residual_scaling;
    options.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    options.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    return options;
}

double coordinate_difference(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return std::max({std::abs(first.x - second.x), std::abs(first.y - second.y), std::abs(first.z - second.z)});
}

bool compare_case(const std::string& name, const std::string& case_path, const std::string& displacement_path,
    const std::string& force_path, const std::string& reaction_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(
        solve.completed && solve.solve.converged && solve.completed_steps == definition.steady_execution.load_steps,
        "H20.29 " + name + " Fuelsim solve completes every load step");
    const std::vector<DisplacementReference> displacement = read_displacements(displacement_path);
    const std::vector<ForceReference> force = read_forces(force_path);
    if (displacement.empty() || force.empty()) throw std::invalid_argument("H20.29 reference fields are empty");
    const std::array<double, 3> normal = displacement.front().normal;
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    fuelsim::test::FieldErrorMetrics displacement_error, force_error;
    std::vector<bool> present(mesh.nodes().size(), false);
    double maximum_coordinate_error = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto reference = std::find_if(displacement.begin(), displacement.end(),
                [source](const DisplacementReference& value) { return value.id == source; });
            if (reference == displacement.end() || present[source])
                throw std::invalid_argument("H20.29 displacement node mapping is incomplete or repeated");
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
    const auto summaries = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    double maximum_force_coordinate_error = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const auto reference = std::find_if(
            force.begin(), force.end(), [&](const ForceReference& value) { return value.id == source_nodes[node]; });
        if (reference == force.end()) throw std::invalid_argument("H20.29 contact-force node mapping is incomplete");
        maximum_force_coordinate_error = std::max(
            maximum_force_coordinate_error, coordinate_difference(mesh.nodes()[reference->id], reference->point));
        force_error.add(summaries[node].contact_force, reference->force);
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    const std::array<double, 2> reference_resultants = read_resultants(reaction_path);
    const double resultant_error =
        std::abs(interface.total_contact_force - reference_resultants[0]) / reference_resultants[0];
    fuelsim::test::print_relative_metrics("h20_29_" + name + "_normal_displacement", displacement_error);
    fuelsim::test::print_relative_metrics("h20_29_" + name + "_nodal_normal_force", force_error);
    std::cout << "h20_29_" << name << "_normal_resultant=" << interface.total_contact_force << '\n'
              << "h20_29_" << name << "_abaqus_normal_resultant=" << reference_resultants[0] << '\n'
              << "h20_29_" << name << "_abaqus_outer_normal_reaction=" << reference_resultants[1] << '\n'
              << "h20_29_" << name << "_normal_resultant_relative_error=" << resultant_error << '\n';
    constexpr double tolerance = 1.0e-2;
    return check(displacement_error.value_count == mesh.nodes().size() &&
                     std::all_of(present.begin(), present.end(), [](bool value) { return value; }),
               "H20.29 " + name + " compares every mesh-node normal displacement") &&
           check(force_error.value_count == summaries.size() && force.size() == summaries.size(),
               "H20.29 " + name + " compares every secondary nodal normal force") &&
           check(maximum_coordinate_error < 1.0e-12 && maximum_force_coordinate_error < 1.0e-12,
               "H20.29 " + name + " Abaqus references use the tracked Exodus coordinates") &&
           check(interface.active_contact_nodes == summaries.size() && interface.unprojected_contact_nodes == 0,
               "H20.29 " + name + " keeps every secondary contact constraint active and projected") &&
           check(fuelsim::test::relative_metrics_below(displacement_error, tolerance),
               "H20.29 " + name + " normal-displacement three Abaqus errors are below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(force_error, tolerance),
               "H20.29 " + name + " nodal-normal-force three Abaqus errors are below 1 percent") &&
           check(resultant_error < tolerance,
               "H20.29 " + name + " normal resultant agrees with Abaqus below 1 percent") &&
           passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--generate") {
        try {
            generate_references(argv[2], argv[3]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] H20.29 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    constexpr std::size_t arguments_per_case = 5;
    const std::size_t argument_count = argc > 0 ? static_cast<std::size_t>(argc - 1) : 0;
    if (argument_count < arguments_per_case || argument_count % arguments_per_case != 0) {
        std::cerr << "Usage: fuelsim_h20_29_hex20_sts_multicase_abaqus_tests "
                     "<name> <case.fsi> <displacement.csv> <force.csv> <reaction.csv> [...]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.29 Abaqus comparison\n");
        bool passed = true;
        for (int argument = 1; argument < argc; argument += static_cast<int>(arguments_per_case))
            passed = compare_case(argv[argument], argv[argument + 1], argv[argument + 2], argv[argument + 3],
                         argv[argument + 4]) &&
                     passed;
        if (passed && session.rank() == 0) std::cout << "[PASS] H20.29 multi-aspect Abaqus STS comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.29 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
