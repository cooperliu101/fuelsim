#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
constexpr std::size_t step_count = 4;

struct PathStep final {
    const char* name;
    std::array<double, 3> pair_a;
    std::array<double, 3> pair_b;
};

constexpr std::array<PathStep, step_count> path = {{{"STICK", {-0.01, 0.012, 0.014}, {-0.01, -0.012, 0.014}},
    {"CROSS", {-0.01, 0.36, 0.42}, {-0.01, -0.36, 0.42}}, {"FORWARD", {-0.01, 0.54, 0.63}, {-0.01, -0.54, 0.63}},
    {"FAR_SLIDE", {-0.01, 0.72, 0.84}, {-0.01, -0.72, 0.84}}}};

constexpr std::array<double, 2> friction = {0.3, 0.5};
constexpr double elastic_slip = 0.02;

struct ContactReference final {
    std::size_t step, id, state;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2, gap, pressure;
};

struct ReactionReference final {
    std::size_t step;
    std::array<std::array<double, 3>, 2> pair;
    std::array<double, 3> global;
};

struct StepState final {
    double time;
    std::vector<double> solution;
    std::array<std::vector<fuelsim::CartesianContactNodeSummary>, 2> contact;
    std::array<std::vector<fuelsim::ContactPointHistory>, 2> histories;
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

std::vector<std::size_t> sorted_nodes(const std::map<std::array<double, 3>, std::size_t>& node_map) {
    std::vector<std::size_t> result;
    for (const auto& entry : node_map) result.push_back(entry.second);
    std::sort(result.begin(), result.end());
    return result;
}

fuelsim::UnstructuredHex8Mesh generate_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::array<std::map<std::array<double, 3>, std::size_t>, 4> node_maps;
    const fuelsim::Hex8Element primary_a_lower = append_cuboid(nodes, node_maps[0], 0.0, 1.0, 0.0, 1.0, -1.0, 2.0),
                               primary_a_upper = append_cuboid(nodes, node_maps[0], 0.0, 1.0, 1.0, 2.0, -1.0, 2.0),
                               secondary_a = append_cuboid(nodes, node_maps[1], 1.0, 2.0, 0.1, 0.9, 0.1, 0.9),
                               primary_b_lower = append_cuboid(nodes, node_maps[2], 3.0, 4.0, -1.0, 2.0, 0.0, 1.0),
                               primary_b_upper = append_cuboid(nodes, node_maps[2], 3.0, 4.0, -1.0, 2.0, 1.0, 2.0),
                               secondary_b = append_cuboid(nodes, node_maps[3], 4.0, 5.0, 0.1, 0.9, 0.1, 0.9);
    const std::vector<std::size_t> primary_a_all = sorted_nodes(node_maps[0]),
                                   secondary_a_all = sorted_nodes(node_maps[1]),
                                   primary_b_all = sorted_nodes(node_maps[2]),
                                   secondary_b_all = sorted_nodes(node_maps[3]);
    std::vector<fuelsim::ElementSide> primary_a_faces, secondary_a_faces, primary_b_faces, secondary_b_faces;
    for (std::size_t element : {std::size_t(0), std::size_t(1)})
        for (std::size_t side = 0; side < 6; ++side) primary_a_faces.push_back({element, side});
    for (std::size_t side = 0; side < 6; ++side) secondary_a_faces.push_back({2, side});
    for (std::size_t element : {std::size_t(3), std::size_t(4)})
        for (std::size_t side = 0; side < 6; ++side) primary_b_faces.push_back({element, side});
    for (std::size_t side = 0; side < 6; ++side) secondary_b_faces.push_back({5, side});
    const std::vector<std::size_t> secondary_a_contact = {
        secondary_a.nodes[0], secondary_a.nodes[3], secondary_a.nodes[4], secondary_a.nodes[7]};
    const std::vector<std::size_t> secondary_b_contact = {
        secondary_b.nodes[0], secondary_b.nodes[3], secondary_b.nodes[4], secondary_b.nodes[7]};
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        {primary_a_lower, primary_a_upper, secondary_a, primary_b_lower, primary_b_upper, secondary_b},
        {1, 1, 2, 3, 3, 4}, {{1, "primary_a"}, {2, "secondary_a"}, {3, "primary_b"}, {4, "secondary_b"}},
        {{10, "primary_a_all", primary_a_all}, {11, "secondary_a_all", secondary_a_all},
            {12, "secondary_a_contact_nodes", secondary_a_contact}, {20, "primary_b_all", primary_b_all},
            {21, "secondary_b_all", secondary_b_all}, {22, "secondary_b_contact_nodes", secondary_b_contact}},
        {{30, "primary_a_contact", {{0, 1}, {1, 1}}}, {31, "secondary_a_contact", {{2, 3}}},
            {40, "primary_b_contact", {{3, 1}, {4, 1}}}, {41, "secondary_b_contact", {{5, 3}}},
            {50, "primary_a_all", primary_a_faces}, {51, "secondary_a_all", secondary_a_faces},
            {60, "primary_b_all", primary_b_faces}, {61, "secondary_b_all", secondary_b_faces}});
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
    if (!output) throw std::runtime_error("Could not write B4.8 Abaqus input: " + output_path);
    output << std::setprecision(16) << "*Heading\n"
           << "** B4.8 C3D8 two-pair finite-sliding surface-to-surface friction path.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const auto& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    const std::array<std::string, 4> block_names = {"PRIMARY_A", "SECONDARY_A", "PRIMARY_B", "SECONDARY_B"};
    for (std::size_t block = 0; block < block_names.size(); ++block) {
        output << "*Element, type=C3D8, elset=" << block_names[block] << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != static_cast<std::int64_t>(block + 1)) continue;
            output << element + 1;
            for (const std::size_t node : mesh.elements()[element].nodes) output << ", " << node + 1;
            output << '\n';
        }
    }
    for (const auto& set : mesh.node_sets()) {
        output << "*Nset, nset=" << set.name << '\n';
        write_labels(output, set.nodes);
    }
    output << "*Surface, type=ELEMENT, name=PRIMARY_A_CONTACT\nPRIMARY_A, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_A_CONTACT\nSECONDARY_A, S6\n"
           << "*Surface, type=ELEMENT, name=PRIMARY_B_CONTACT\nPRIMARY_B, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_B_CONTACT\nSECONDARY_B, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.25\n";
    for (const auto& block : block_names) output << "*Solid Section, elset=" << block << ", material=ELASTIC\n,\n";
    output << "*Surface Interaction, name=PAIR_A_INTERACTION\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e5,\n"
           << "*Friction, slip tolerance=0.025\n0.3,\n"
           << "*Surface Interaction, name=PAIR_B_INTERACTION\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e5,\n"
           << "*Friction, slip tolerance=0.025\n0.5,\n"
           << "*Contact Pair, interaction=PAIR_A_INTERACTION, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_A_CONTACT, PRIMARY_A_CONTACT\n"
           << "*Contact Pair, interaction=PAIR_B_INTERACTION, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_B_CONTACT, PRIMARY_B_CONTACT\n";
    for (std::size_t step = 0; step < path.size(); ++step) {
        output << "*Step, name=" << path[step].name << ", nlgeom=YES, inc=100\n"
               << "*Static\n0.1, 1., 1.e-8, 0.1\n"
               << (step == 0 ? "*Boundary\nprimary_a_all, 1, 3, 0.\nprimary_b_all, 1, 3, 0.\n" : "*Boundary, op=MOD\n");
        for (std::size_t component = 0; component < 3; ++component) {
            output << "secondary_a_all, " << component + 1 << ", " << component + 1 << ", "
                   << path[step].pair_a[component] << '\n';
            output << "secondary_b_all, " << component + 1 << ", " << component + 1 << ", "
                   << path[step].pair_b[component] << '\n';
        }
        output << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n"
               << "*Output, history, frequency=1\n*Energy Output\nALLFD\n*End Step\n";
    }
}

void write_fsi(const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("Could not write B4.8 Fuelsim input: " + output_path);
    output << "[Case]\n  version = 3\n  problem = transient\n  geometry = cartesian_3d\n[]\n\n"
           << "[Mesh]\n  type = exodus\n  file = ../abaqus/b48_hex8_multi_contact_mesh.e\n[]\n\n"
           << "[Materials]\n  [elastic]\n    [thermal]\n      function = constant_thermophysical\n"
           << "      conductivity = 1\n      density = 1\n      specific_heat = 1\n    []\n"
           << "    [elasticity]\n      function = constant_isotropic\n      young_modulus = 1e9\n"
           << "      poisson_ratio = 0.25\n    []\n  []\n[]\n\n[Regions]\n";
    for (const std::string& region :
        {std::string("primary_a"), std::string("secondary_a"), std::string("primary_b"), std::string("secondary_b")})
        output << "  [" << region << "]\n    block = " << region
               << "\n    material = elastic\n    strain = finite\n    initial_temperature = 300\n"
               << "    volumetric_heat_source = 0\n  []\n";
    output << "[]\n\n[TimeFunctions]\n";
    const std::array<std::pair<std::string, std::array<double, step_count>>, 6> functions = {
        {{"pair_a_x", {-0.01, -0.01, -0.01, -0.01}}, {"pair_a_y", {0.012, 0.36, 0.54, 0.72}},
            {"pair_a_z", {0.014, 0.42, 0.63, 0.84}}, {"pair_b_x", {-0.01, -0.01, -0.01, -0.01}},
            {"pair_b_y", {-0.012, -0.36, -0.54, -0.72}}, {"pair_b_z", {0.014, 0.42, 0.63, 0.84}}}};
    for (const auto& function : functions) {
        output << "  [" << function.first << "]\n    type = piecewise_linear\n    times = 0 1 2 3 4\n"
               << "    values = 0";
        for (double value : function.second) output << ' ' << value;
        output << "\n  []\n";
    }
    output << "[]\n\n[Contact]\n";
    for (std::size_t pair = 0; pair < 2; ++pair) {
        const char letter = static_cast<char>('a' + pair);
        output << "  [pair_" << letter << "]\n    primary = primary_" << letter
               << "_contact\n    secondary = secondary_" << letter << "_contact\n"
               << "    [mechanical]\n      formulation = penalty\n      discretization = surface_to_surface\n"
               << "      sliding = finite\n      penalty = 1e5\n      mu = " << friction[pair]
               << "\n      slip_tolerance = 0.025\n    []\n  []\n";
    }
    output << "[]\n\n[BoundaryConditions]\n";
    for (const std::string& region :
        {std::string("primary_a"), std::string("secondary_a"), std::string("primary_b"), std::string("secondary_b")})
        output << "  [" << region << "_temperature]\n    type = dirichlet\n    boundary = " << region
               << "_all\n    field = temperature\n    value = 300\n  []\n";
    for (const std::string& primary : {std::string("primary_a"), std::string("primary_b")})
        for (const std::string& component : {std::string("x"), std::string("y"), std::string("z")})
            output << "  [" << primary << '_' << component << "]\n    type = dirichlet\n    boundary = " << primary
                   << "_all\n    field = displacement_" << component << "\n    value = 0\n  []\n";
    for (const std::string& secondary : {std::string("secondary_a"), std::string("secondary_b")})
        for (const std::string& component : {std::string("x"), std::string("y"), std::string("z")})
            output << "  [" << secondary << '_' << component << "]\n    type = dirichlet\n    boundary = " << secondary
                   << "_all\n    field = displacement_" << component << "\n    value = 1\n    function = pair_"
                   << secondary.back() << '_' << component << "\n  []\n";
    output << "[]\n\n[Executioner]\n  type = transient\n  end_time = 4\n  initial_time_step = 1\n"
           << "  minimum_time_step = 0.125\n  maximum_time_step = 1\n  growth_factor = 1\n"
           << "  cutback_factor = 0.5\n  maximum_cutbacks = 3\n  load_ramp_time = 0\n[]\n\n"
           << "[Solver]\n  absolute_tolerance = 1e-8\n  relative_tolerance = 1e-11\n"
           << "  step_tolerance = 1e-12\n  maximum_iterations = 40\n  linear_solver = direct\n"
           << "  direct_factorization = mumps\n  field_residual_scaling = true\n"
           << "  linear_relative_tolerance = 1e-11\n  maximum_linear_iterations = 400\n[]\n\n"
           << "[Outputs]\n  console = true\n[]\n";
}

void generate(const std::string& abaqus_directory, const std::string& fuelsim_directory) {
    const fuelsim::UnstructuredHex8Mesh mesh = generate_mesh();
    fuelsim::write_exodus_hex8(abaqus_directory + "/b48_hex8_multi_contact_mesh.e", mesh);
    write_abaqus_input(abaqus_directory + "/b48_hex8_multi_contact.inp", mesh);
    write_fsi(fuelsim_directory + "/transient_hex8_multi_contact_path_abaqus.fsi");
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete B4.8 CSV row: " + input_path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B4.8 CSV number: " + input_path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    const double value = number(values, column, input_path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("Invalid B4.8 CSV index: " + input_path);
    return static_cast<std::size_t>(value);
}

std::vector<ContactReference> read_contact(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Could not read B4.8 contact reference: " + input_path);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,"
                                              "tangential_z,slip1,slip2,gap,pressure,state")
        throw std::invalid_argument("Unexpected B4.8 contact CSV header: " + input_path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, input_path), index_value(values, 1, input_path),
            index_value(values, 15, input_path),
            {number(values, 2, input_path), number(values, 3, input_path), number(values, 4, input_path)},
            {number(values, 5, input_path), number(values, 6, input_path), number(values, 7, input_path)},
            {number(values, 8, input_path), number(values, 9, input_path), number(values, 10, input_path)},
            number(values, 11, input_path), number(values, 12, input_path), number(values, 13, input_path),
            number(values, 14, input_path)});
    }
    if (result.size() != step_count * 4) throw std::invalid_argument("B4.8 contact reference must have sixteen rows");
    return result;
}

std::vector<ReactionReference> read_reactions(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Could not read B4.8 reaction reference: " + input_path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,pair_a_x,pair_a_y,pair_a_z,pair_b_x,pair_b_y,pair_b_z,global_x,global_y,global_z")
        throw std::invalid_argument("Unexpected B4.8 reaction CSV header: " + input_path);
    std::vector<ReactionReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, input_path),
            {{{number(values, 1, input_path), number(values, 2, input_path), number(values, 3, input_path)},
                {number(values, 4, input_path), number(values, 5, input_path), number(values, 6, input_path)}}},
            {number(values, 7, input_path), number(values, 8, input_path), number(values, 9, input_path)}});
    }
    if (result.size() != step_count) throw std::invalid_argument("B4.8 reaction reference must have four rows");
    return result;
}

std::vector<double> read_energy(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Could not read B4.8 energy reference: " + input_path);
    std::string line;
    if (!std::getline(input, line) || line != "step,allfd")
        throw std::invalid_argument("Unexpected B4.8 energy CSV header: " + input_path);
    std::vector<double> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        if (index_value(values, 0, input_path) != result.size() + 1)
            throw std::invalid_argument("B4.8 energy steps are not consecutive");
        result.push_back(number(values, 1, input_path));
    }
    if (result.size() != step_count) throw std::invalid_argument("B4.8 energy reference must have four rows");
    return result;
}

class Recorder final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        StepState state;
        state.time = step.time;
        state.solution = problem.committed_solution();
        for (std::size_t pair = 0; pair < 2; ++pair) {
            state.contact[pair] =
                fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, pair, state.solution);
            state.histories[pair] = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(pair);
        }
        states.push_back(std::move(state));
    }

    std::vector<StepState> states;
};

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& input) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = input.solver.absolute_tolerance;
    result.relative_tolerance = input.solver.relative_tolerance;
    result.step_tolerance = input.solver.step_tolerance;
    result.maximum_iterations = input.solver.maximum_iterations;
    result.linear_solver = input.solver.linear_solver;
    result.direct_factorization = input.solver.direct_factorization;
    result.field_residual_scaling = input.solver.field_residual_scaling;
    result.linear_relative_tolerance = input.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = input.solver.maximum_linear_iterations;
    return result;
}

fuelsim::TransientTimeOptions time_options(
    const fuelsim::FuelSimCaseDefinition& input, double end_time, double initial_time_step) {
    fuelsim::TransientTimeOptions result = input.transient_execution;
    result.end_time = end_time;
    result.initial_time_step = initial_time_step;
    return result;
}

bool histories_equal(const std::vector<std::vector<fuelsim::ContactPointHistory>>& first,
    const std::vector<std::vector<fuelsim::ContactPointHistory>>& second) {
    if (first.size() != second.size()) return false;
    for (std::size_t pair = 0; pair < first.size(); ++pair) {
        if (first[pair].size() != second[pair].size()) return false;
        for (std::size_t point = 0; point < first[pair].size(); ++point) {
            const auto& left = first[pair][point];
            const auto& right = second[pair][point];
            if (left.elastic_tangential_slip != right.elastic_tangential_slip || left.sliding != right.sliding ||
                left.normal_multiplier != right.normal_multiplier ||
                left.cartesian_elastic_tangential_slip != right.cartesian_elastic_tangential_slip ||
                left.cartesian_total_tangential_slip != right.cartesian_total_tangential_slip ||
                left.cartesian_tangent_basis_initialized != right.cartesian_tangent_basis_initialized ||
                left.cartesian_contact_normal != right.cartesian_contact_normal ||
                left.cartesian_contact_tangent_first != right.cartesian_contact_tangent_first)
                return false;
        }
    }
    return true;
}

double jacobian_error(const fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double perturbation = 1.0e-7;
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
        }
        if (reference_squared > 0.0) maximum = std::max(maximum, std::sqrt(difference_squared / reference_squared));
    }
    return maximum;
}

double contact_action_reaction_error(const fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, nullptr);
        std::array<double, 3> resultant{};
        const auto& fields = spatial.field_layout();
        for (std::size_t row = 0; row < dofs.size(); ++row)
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = fields.at(component + 1);
                if (dofs[row] >= field.begin && dofs[row] < field.end) resultant[component] += residual[row];
            }
        for (double value : resultant) maximum = std::max(maximum, std::abs(value));
    }
    return maximum;
}

std::array<double, 3> reference_slip(const ContactReference& value) { return {0.0, -value.slip_2, value.slip_1}; }

std::array<double, 3> reference_elastic_slip(const ContactReference& value, double mu) {
    std::array<double, 3> result{};
    const double normal =
        std::sqrt(value.normal_force[0] * value.normal_force[0] + value.normal_force[1] * value.normal_force[1] +
                  value.normal_force[2] * value.normal_force[2]);
    if (normal == 0.0 || mu == 0.0) return result;
    for (std::size_t component = 0; component < 3; ++component)
        result[component] = -value.tangential_force[component] * elastic_slip / (mu * normal);
    return result;
}

bool compare(const fuelsim::TransientProblem& problem, const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::vector<StepState>& states, const std::array<std::vector<ContactReference>, 2>& references,
    const std::vector<ReactionReference>& reactions, const std::vector<double>& allfd) {
    if (states.size() != step_count) throw std::invalid_argument("B4.8 Fuelsim path must contain four states");
    std::array<std::vector<std::size_t>, 2> sources;
    for (std::size_t pair = 0; pair < 2; ++pair)
        sources[pair] = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, pair);
    std::array<fuelsim::test::FieldErrorMetrics, 2> normal, tangent_y, tangent_z, slip_y, slip_z, gap, pressure,
        resultant, force_center, dissipation;
    fuelsim::test::FieldErrorMetrics global_reaction, global_dissipation;
    std::array<std::vector<std::array<double, 3>>, 2> previous_actual_plastic, previous_reference_plastic;
    for (std::size_t pair = 0; pair < 2; ++pair) {
        previous_actual_plastic[pair].resize(sources[pair].size());
        previous_reference_plastic[pair].resize(sources[pair].size());
    }
    std::array<double, 2> actual_cumulative_dissipation{}, reference_cumulative_dissipation{};
    double maximum_coordinate_difference = 0.0, maximum_zero_reference_difference = 0.0,
           maximum_action_reaction_difference = 0.0, maximum_abaqus_shear_normal_component = 0.0,
           maximum_abaqus_shear_normal_resultant = 0.0;
    bool all_active = true, first_sticks = true, later_slides = true, crossed_different_primary_faces = false,
         states_match = true;
    std::array<std::size_t, 2> first_primary_face = {
        std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()};
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<std::array<double, 3>, 2> actual_resultant{};
        std::array<double, 3> actual_global{};
        for (std::size_t pair = 0; pair < 2; ++pair) {
            std::array<double, 3> actual_center_sum{}, reference_center_sum{};
            double actual_center_weight = 0.0, reference_center_weight = 0.0;
            double actual_increment_dissipation = 0.0, reference_increment_dissipation = 0.0;
            double abaqus_shear_normal_resultant = 0.0;
            std::size_t sliding_nodes = 0;
            for (std::size_t node = 0; node < states[step].contact[pair].size(); ++node) {
                const auto found =
                    std::find_if(references[pair].begin(), references[pair].end(), [&](const ContactReference& value) {
                        return value.step == step + 1 && value.id == sources[pair][node];
                    });
                if (found == references[pair].end()) throw std::invalid_argument("B4.8 contact mapping is incomplete");
                const auto& actual = states[step].contact[pair][node];
                const auto& point = mesh.nodes().at(sources[pair][node]);
                maximum_coordinate_difference =
                    std::max({maximum_coordinate_difference, std::abs(point.x - found->point.x),
                        std::abs(point.y - found->point.y), std::abs(point.z - found->point.z)});
                normal[pair].add(-actual.normal_contact_force[0], found->normal_force[0]);
                tangent_y[pair].add(-actual.tangential_contact_force[1], found->tangential_force[1]);
                tangent_z[pair].add(-actual.tangential_contact_force[2], found->tangential_force[2]);
                const auto& displacement = pair == 0 ? path[step].pair_a : path[step].pair_b;
                slip_y[pair].add(displacement[1], -found->slip_2);
                slip_z[pair].add(displacement[2], found->slip_1);
                gap[pair].add(actual.gap, found->gap);
                pressure[pair].add(actual.pressure, found->pressure);
                maximum_zero_reference_difference =
                    std::max({maximum_zero_reference_difference, std::abs(actual.normal_contact_force[1]),
                        std::abs(actual.normal_contact_force[2]), std::abs(actual.tangential_contact_force[0]),
                        std::abs(found->normal_force[1]), std::abs(found->normal_force[2])});
                maximum_abaqus_shear_normal_component =
                    std::max(maximum_abaqus_shear_normal_component, std::abs(found->tangential_force[0]));
                abaqus_shear_normal_resultant += found->tangential_force[0];
                const std::array<double, 3> actual_force = {
                    -actual.normal_contact_force[0] - actual.tangential_contact_force[0],
                    -actual.normal_contact_force[1] - actual.tangential_contact_force[1],
                    -actual.normal_contact_force[2] - actual.tangential_contact_force[2]};
                const std::array<double, 3> actual_current = {
                    point.x + displacement[0], point.y + displacement[1], point.z + displacement[2]};
                const std::array<double, 3> reference_current = actual_current;
                for (std::size_t component = 0; component < 3; ++component) {
                    actual_resultant[pair][component] += actual_force[component];
                    actual_global[component] += actual_force[component];
                    actual_center_sum[component] += actual.contact_force * actual_current[component];
                    reference_center_sum[component] += std::abs(found->normal_force[0]) * reference_current[component];
                }
                actual_center_weight += actual.contact_force;
                reference_center_weight += std::abs(found->normal_force[0]);
                const std::array<double, 3> reference_total_slip = reference_slip(*found),
                                            reference_elastic = reference_elastic_slip(*found, friction[pair]);
                std::array<double, 3> actual_plastic{}, reference_plastic{};
                for (std::size_t component = 0; component < 3; ++component) {
                    const double actual_total_slip = component == 0 ? 0.0 : displacement[component];
                    actual_plastic[component] =
                        actual_total_slip -
                        states[step].histories[pair][node].cartesian_elastic_tangential_slip[component];
                    reference_plastic[component] = reference_total_slip[component] - reference_elastic[component];
                    if (component > 0) {
                        actual_increment_dissipation -=
                            (-actual.tangential_contact_force[component]) *
                            (actual_plastic[component] - previous_actual_plastic[pair][node][component]);
                        reference_increment_dissipation -=
                            found->tangential_force[component] *
                            (reference_plastic[component] - previous_reference_plastic[pair][node][component]);
                    }
                }
                previous_actual_plastic[pair][node] = actual_plastic;
                previous_reference_plastic[pair][node] = reference_plastic;
                const std::size_t actual_state = actual.pressure <= 0.0 ? 0 : (actual.sliding ? 2 : 1);
                states_match = states_match && actual_state == found->state;
                all_active = all_active && actual.projected && actual.pressure > 0.0 && found->pressure > 0.0;
                first_sticks = first_sticks && (step != 0 || (!actual.sliding && found->state == 1));
                if (actual.sliding) ++sliding_nodes;
                if (step == 0 && node == 0) first_primary_face[pair] = actual.primary_face;
                if (step > 0 && actual.primary_face != first_primary_face[pair]) crossed_different_primary_faces = true;
            }
            later_slides = later_slides && (step == 0 || sliding_nodes > 0);
            maximum_abaqus_shear_normal_resultant =
                std::max(maximum_abaqus_shear_normal_resultant, std::abs(abaqus_shear_normal_resultant));
            actual_cumulative_dissipation[pair] += actual_increment_dissipation;
            reference_cumulative_dissipation[pair] += reference_increment_dissipation;
            dissipation[pair].add(
                actual_cumulative_dissipation[pair], step == 0 ? 0.0 : reference_cumulative_dissipation[pair]);
            for (std::size_t component = 0; component < 3; ++component) {
                resultant[pair].add(actual_resultant[pair][component], -reactions[step].pair[pair][component]);
                force_center[pair].add(actual_center_sum[component] / actual_center_weight,
                    reference_center_sum[component] / reference_center_weight);
                maximum_action_reaction_difference = std::max(maximum_action_reaction_difference,
                    std::abs(actual_resultant[pair][component] + reactions[step].pair[pair][component]));
            }
        }
        const double actual_total_dissipation = actual_cumulative_dissipation[0] + actual_cumulative_dissipation[1];
        global_dissipation.add(actual_total_dissipation, allfd[step]);
        for (std::size_t component = 0; component < 3; ++component)
            global_reaction.add(actual_global[component], -reactions[step].global[component]);
    }
    const auto print = [](const std::string& name, const fuelsim::test::FieldErrorMetrics& metric) {
        if (metric.has_relative_norm())
            fuelsim::test::print_relative_metrics(name, metric);
        else
            fuelsim::test::print_absolute_metrics(name, metric);
    };
    for (std::size_t pair = 0; pair < 2; ++pair) {
        const std::string prefix = std::string("b48_pair_") + static_cast<char>('a' + pair) + '_';
        print(prefix + "normal_force_x", normal[pair]);
        print(prefix + "tangential_force_y", tangent_y[pair]);
        print(prefix + "tangential_force_z", tangent_z[pair]);
        print(prefix + "slip_y", slip_y[pair]);
        print(prefix + "slip_z", slip_z[pair]);
        print(prefix + "gap", gap[pair]);
        print(prefix + "pressure", pressure[pair]);
        print(prefix + "resultant", resultant[pair]);
        print(prefix + "force_center", force_center[pair]);
        print(prefix + "friction_dissipation", dissipation[pair]);
    }
    print("b48_global_reaction", global_reaction);
    print("b48_global_friction_dissipation", global_dissipation);
    std::cout << "b48_maximum_mesh_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << "b48_maximum_zero_reference_force_difference=" << maximum_zero_reference_difference << '\n'
              << "b48_maximum_action_reaction_difference=" << maximum_action_reaction_difference << '\n'
              << "b48_abaqus_maximum_nodal_shear_normal_component=" << maximum_abaqus_shear_normal_component << '\n'
              << "b48_abaqus_maximum_shear_normal_resultant=" << maximum_abaqus_shear_normal_resultant << '\n';
    constexpr double tolerance = 1.0e-2, zero_tolerance = 1.0e-7;
    const auto passes = [&](const fuelsim::test::FieldErrorMetrics& metric) {
        return (!metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, tolerance)) &&
               metric.maximum_zero_reference_difference < zero_tolerance;
    };
    bool fields_pass = passes(global_reaction) && passes(global_dissipation);
    for (std::size_t pair = 0; pair < 2; ++pair)
        fields_pass = fields_pass && passes(normal[pair]) && passes(tangent_y[pair]) && passes(tangent_z[pair]) &&
                      passes(slip_y[pair]) && passes(slip_z[pair]) && passes(gap[pair]) && passes(pressure[pair]) &&
                      passes(resultant[pair]) && passes(force_center[pair]) && passes(dissipation[pair]);
    return check(maximum_coordinate_difference < 3.0e-8, "B4.8 Abaqus and Fuelsim use the same tracked mesh") &&
           check(all_active && first_sticks && later_slides && crossed_different_primary_faces,
               "B4.8 keeps both pairs active, starts in sticking, slides both pairs, and transfers primary-face "
               "ownership") &&
           check(
               states_match, "B4.8 open, sticking, and sliding states agree with Abaqus at every accepted increment") &&
           check(maximum_zero_reference_difference < zero_tolerance,
               "B4.8 theoretical-zero contact-force components pass the separate absolute check") &&
           check(maximum_abaqus_shear_normal_resultant < 1.0e-10,
               "B4.8 Abaqus nodal shear leakage normal to the plane cancels in each pair resultant") &&
           check(fields_pass, "B4.8 all per-pair nodal fields, resultants, force centers, friction dissipation, and "
                              "global reactions pass the one-percent metrics");
}

bool run(const std::string& input_path, const std::array<std::string, 2>& contact_paths,
    const std::string& reaction_path, const std::string& energy_path, const std::string& checkpoint_path,
    const fuelsim::PetscSession& session) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    Recorder recorder;
    fuelsim::TransientProblem full(input.spatial, mesh);
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(full, time_options(input, 4.0, 1.0), solver, &recorder);
    bool passed = check(result.completed && recorder.states.size() == step_count && result.rejected_steps.empty(),
        "B4.8 completes the same four accepted increments for both contact pairs without rejected steps");
    const std::array<std::vector<ContactReference>, 2> contacts = {
        read_contact(contact_paths[0]), read_contact(contact_paths[1])};
    passed = compare(full, mesh, recorder.states, contacts, read_reactions(reaction_path), read_energy(energy_path)) &&
             passed;

    Recorder split_recorder;
    fuelsim::TransientProblem split(input.spatial, mesh);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(input, 2.0, 1.0), solver, &split_recorder);
    passed = check(first.completed && split_recorder.states.size() == 2,
                 "B4.8 reaches the two-pair finite-sliding checkpoint") &&
             passed;
    session.collective_root_action(
        [&]() { fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step); });
    fuelsim::TransientProblem restarted(input.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    Recorder restart_recorder;
    const fuelsim::TransientResult second =
        fuelsim::solve_transient(restarted, time_options(input, 4.0, restored_step), solver, &restart_recorder);
    double maximum_restart_difference = 0.0;
    if (second.completed && restarted.committed_solution().size() == full.committed_solution().size())
        for (std::size_t dof = 0; dof < full.committed_solution().size(); ++dof)
            maximum_restart_difference = std::max(maximum_restart_difference,
                std::abs(restarted.committed_solution()[dof] - full.committed_solution()[dof]));
    const bool history_match =
        histories_equal(fuelsim::cartesian::ProblemAccess::committed_contact_histories(restarted),
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(full));
    passed = check(second.completed && restart_recorder.states.size() == 2 && maximum_restart_difference == 0.0 &&
                       history_match,
                 "B4.8 checkpoint restart exactly reproduces both contact pairs and their vector friction histories") &&
             passed;
    std::vector<double> sliding_trial = full.committed_solution();
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(full);
    for (std::size_t pair = 0; pair < 2; ++pair) {
        const std::size_t region = pair == 0 ? 1 : 3;
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local) {
            const std::size_t global = spatial.global_node(region, local);
            sliding_trial[spatial.dof(fuelsim::Field::displacement_y, global)] += pair == 0 ? 1.2e-3 : -1.2e-3;
            sliding_trial[spatial.dof(fuelsim::Field::displacement_z, global)] += 1.4e-3;
        }
    }
    spatial.validate_state(sliding_trial);
    const double directional_error = jacobian_error(full, sliding_trial);
    const double action_reaction_error = contact_action_reaction_error(full, sliding_trial);
    std::cout << "b48_contact_jacobian_directional_error=" << directional_error << '\n'
              << "b48_contact_action_reaction_maximum_absolute=" << action_reaction_error << '\n'
              << "b48_restart_maximum_absolute_difference=" << maximum_restart_difference << '\n';
    passed = check(directional_error < 1.0e-5,
                 "B4.8 sliding contact Jacobians match centered directional differences away from transitions") &&
             check(action_reaction_error < 1.0e-10,
                 "B4.8 every mechanical contact constraint has equal and opposite three-component force") &&
             passed;
    session.collective_root_action([&]() {
        if (std::remove(checkpoint_path.c_str()) != 0 && errno != ENOENT)
            throw std::runtime_error("B4.8 could not remove its checkpoint artifact");
    });
    passed = check(true, "B4.8 removes its checkpoint artifact") && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--generate") {
        try {
            generate(argv[2], argv[3]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] B4.8 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc != 7) {
        std::cerr << "Usage: fuelsim_b48_hex8_multi_contact_path_abaqus_tests <case.fsi> <pair-a.csv> "
                     "<pair-b.csv> <reaction.csv> <energy.csv> <checkpoint>\n"
                     "   or: fuelsim_b48_hex8_multi_contact_path_abaqus_tests --generate <abaqus-directory> "
                     "<fuelsim-directory>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B4.8 HEX8 two-pair finite-sliding Abaqus comparison\n");
        const bool passed = run(argv[1], {argv[2], argv[3]}, argv[4], argv[5], argv[6], session);
        if (passed && session.rank() == 0) std::cout << "[PASS] B4.8 HEX8 two-pair finite-sliding Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.8 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
