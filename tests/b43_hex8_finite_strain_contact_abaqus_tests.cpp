#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/contact.hpp"
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
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    std::size_t step, source_node;
    fuelsim::CartesianPoint3 reference;
    std::array<double, 3> displacement{}, reaction{};
};

struct ContactReference final {
    std::size_t step, source_node;
    fuelsim::CartesianPoint3 current;
    std::array<double, 3> normal_force{}, tangential_force{};
    double slip_1, slip_2, gap, pressure;
};

struct StepState final {
    double time;
    std::vector<double> solution;
    std::vector<fuelsim::CartesianContactNodeSummary> contact;
};

constexpr double end_time = 20.0;

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::UnstructuredHex8Mesh generate_mesh() {
    const std::vector<fuelsim::CartesianPoint3> nodes = {{0.0, 0.0, -0.25}, {1.0, 0.0, -0.25}, {1.0, 1.0, -0.25},
        {0.0, 1.0, -0.25}, {0.0, 0.0, 1.25}, {1.0, 0.0, 1.25}, {1.0, 1.0, 1.25}, {0.0, 1.0, 1.25}, {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {2.0, 1.0, 1.0},
        {1.0, 1.0, 1.0}};
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{0, 1, 2, 3, 4, 5, 6, 7}}}, {{{8, 9, 10, 11, 12, 13, 14, 15}}}};
    std::vector<fuelsim::NodeSet> node_sets = {{10, "primary_all_nodes", {0, 1, 2, 3, 4, 5, 6, 7}},
        {20, "secondary_all_nodes", {8, 9, 10, 11, 12, 13, 14, 15}}, {30, "secondary_contact_nodes", {8, 11, 12, 15}}};
    const std::vector<fuelsim::ElementSide> primary_all = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}},
                                            secondary_all = {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
    return fuelsim::UnstructuredHex8Mesh(nodes, elements, {1, 2}, {{1, "primary"}, {2, "secondary"}},
        std::move(node_sets),
        {{40, "primary_all", primary_all}, {41, "primary_y0", {{0, 0}}}, {42, "primary_y1", {{0, 2}}},
            {43, "primary_z0", {{0, 4}}}, {44, "primary_z1", {{0, 5}}}, {50, "primary_contact", {{0, 1}}},
            {60, "secondary_all", secondary_all}, {61, "secondary_y0", {{1, 0}}}, {62, "secondary_y1", {{1, 2}}},
            {63, "secondary_z0", {{1, 4}}}, {64, "secondary_z1", {{1, 5}}}, {70, "secondary_contact", {{1, 3}}}});
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
    if (!output) throw std::runtime_error("Could not write B4.3 Abaqus input: " + output_path);
    output << std::setprecision(16) << "*Heading\n"
           << "** B4.3 deformable C3D8 finite-strain finite-sliding STS contact.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    output << "*Element, type=C3D8, elset=PRIMARY\n1";
    for (const std::size_t node : mesh.elements()[0].nodes) output << ", " << node + 1;
    output << "\n*Element, type=C3D8, elset=SECONDARY\n2";
    for (const std::size_t node : mesh.elements()[1].nodes) output << ", " << node + 1;
    output << "\n*Nset, nset=PRIMARY_ALL_NODES\n";
    write_labels(output, mesh.node_sets()[0].nodes);
    output << "*Nset, nset=SECONDARY_ALL_NODES\n";
    write_labels(output, mesh.node_sets()[1].nodes);
    output << "*Nset, nset=SECONDARY_CONTACT_NODES\n";
    write_labels(output, mesh.node_sets()[2].nodes);
    output << "*Nset, nset=PRIMARY_Y0\n1, 2, 5, 6\n"
           << "*Nset, nset=PRIMARY_Y1\n3, 4, 7, 8\n"
           << "*Nset, nset=PRIMARY_Z0\n1, 2, 3, 4\n"
           << "*Nset, nset=PRIMARY_Z1\n5, 6, 7, 8\n"
           << "*Nset, nset=SECONDARY_Y0\n9, 10, 13, 14\n"
           << "*Nset, nset=SECONDARY_Y1\n11, 12, 15, 16\n"
           << "*Nset, nset=SECONDARY_Z0\n9, 10, 11, 12\n"
           << "*Nset, nset=SECONDARY_Z1\n13, 14, 15, 16\n";
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e7, 0.\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e8,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
           << "*Step, name=DEFORM, nlgeom=YES, inc=20\n"
           << "*Static\n1., 20., 1., 1.\n"
           << "*Boundary\n"
           << "PRIMARY_Y0, 1, 1, -0.015\nPRIMARY_Y1, 1, 1, 0.015\n"
           << "PRIMARY_Y0, 2, 2, -0.05\nPRIMARY_Y1, 2, 2, 0.05\n"
           << "PRIMARY_Z0, 3, 3, -0.06\nPRIMARY_Z1, 3, 3, 0.06\n"
           << "SECONDARY_Y0, 1, 1, -0.095\nSECONDARY_Y1, 1, 1, -0.065\n"
           << "SECONDARY_Y0, 2, 2, -0.05\nSECONDARY_Y1, 2, 2, 0.05\n"
           << "SECONDARY_Z0, 3, 3, -0.02\nSECONDARY_Z1, 3, 3, 0.06\n";
    output << "*Output, field, frequency=1\n"
           << "*Node Output\nCOORD, RF, U\n"
           << "*Contact Output\nCSTRESS, CDISP, CFORCE\n"
           << "*End Step\n";
}

void write_fsi(const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("Could not write B4.3 Fuelsim input: " + output_path);
    output << std::setprecision(16) << "[Case]\n  version = 3\n  problem = transient\n"
           << "  geometry = cartesian_3d\n[]\n\n"
           << "[Mesh]\n  type = exodus\n  file = ../abaqus/b43_hex8_finite_strain_contact_mesh.e\n[]\n\n"
           << "[Materials]\n  [elastic]\n"
           << "    [thermal]\n      function = constant_thermophysical\n      conductivity = 1\n"
           << "      density = 1\n      specific_heat = 1\n    []\n"
           << "    [elasticity]\n      function = constant_isotropic\n      young_modulus = 1e7\n"
           << "      poisson_ratio = 0\n    []\n  []\n[]\n\n"
           << "[Regions]\n"
           << "  [primary]\n    block = primary\n    material = elastic\n    strain = finite\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n"
           << "  [secondary]\n    block = secondary\n    material = elastic\n    strain = finite\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n[]\n\n"
           << "[TimeFunctions]\n";
    const std::array<std::pair<const char*, double>, 12> drives = {
        {{"primary_y0_x", -0.015}, {"primary_y1_x", 0.015}, {"primary_y0_y", -0.05}, {"primary_y1_y", 0.05},
            {"primary_z0_z", -0.06}, {"primary_z1_z", 0.06}, {"secondary_y0_x", -0.095}, {"secondary_y1_x", -0.065},
            {"secondary_y0_y", -0.05}, {"secondary_y1_y", 0.05}, {"secondary_z0_z", -0.02}, {"secondary_z1_z", 0.06}}};
    for (const auto& drive : drives)
        output << "  [" << drive.first << "]\n    type = piecewise_linear\n    times = 0 20\n    values = 0 "
               << drive.second << "\n  []\n";
    output << "[]\n\n[Contact]\n  [interface]\n    primary = primary_contact\n"
           << "    secondary = secondary_contact\n    [mechanical]\n      formulation = penalty\n"
           << "      discretization = surface_to_surface\n      sliding = finite\n      penalty = 1e8\n"
           << "      mu = 0\n    []\n  []\n[]\n\n"
           << "[BoundaryConditions]\n"
           << "  [primary_temperature]\n    type = dirichlet\n    boundary = primary_all\n"
           << "    field = temperature\n    value = 300\n  []\n"
           << "  [secondary_temperature]\n    type = dirichlet\n    boundary = secondary_all\n"
           << "    field = temperature\n    value = 300\n  []\n";
    for (const auto& drive : drives) {
        const std::string name = drive.first;
        const std::size_t final_separator = name.rfind('_');
        output << "  [" << name << "]\n    type = dirichlet\n    boundary = " << name.substr(0, final_separator)
               << "\n    field = displacement_" << name.substr(final_separator + 1)
               << "\n    value = 1\n    function = " << name << "\n  []\n";
    }
    output << "[]\n\n[Executioner]\n  type = transient\n  end_time = 20\n  initial_time_step = 1\n"
           << "  minimum_time_step = 1\n  maximum_time_step = 1\n  growth_factor = 1\n"
           << "  cutback_factor = 0.5\n  maximum_cutbacks = 0\n  load_ramp_time = 0\n"
           << "  include_thermal_time_term = false\n[]\n\n"
           << "[Solver]\n  absolute_tolerance = 1e-7\n  relative_tolerance = 1e-11\n"
           << "  step_tolerance = 1e-12\n  maximum_iterations = 40\n  linear_solver = direct\n"
           << "  direct_factorization = mumps\n  field_residual_scaling = true\n"
           << "  linear_relative_tolerance = 1e-11\n  maximum_linear_iterations = 400\n[]\n\n"
           << "[Outputs]\n  console = true\n[]\n";
}

void generate(const std::string& abaqus_directory, const std::string& fuelsim_directory) {
    const fuelsim::UnstructuredHex8Mesh mesh = generate_mesh();
    fuelsim::write_exodus_hex8(abaqus_directory + "/b43_hex8_finite_strain_contact_mesh.e", mesh);
    write_abaqus_input(abaqus_directory + "/b43_hex8_finite_strain_contact.inp", mesh);
    write_fsi(fuelsim_directory + "/transient_hex8_finite_strain_contact_abaqus.fsi");
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete B4.3 CSV row: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B4.3 CSV number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double value = number(values, column, path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("Invalid B4.3 CSV index: " + path);
    return static_cast<std::size_t>(value);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.3 node reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,ux,uy,uz,rf_x,rf_y,rf_z")
        throw std::invalid_argument("Unexpected B4.3 node CSV header: " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path), index_value(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.3 contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,slip_1,"
                "slip_2,gap,pressure")
        throw std::invalid_argument("Unexpected B4.3 contact CSV header: " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path), index_value(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}, number(values, 11, path),
            number(values, 12, path), number(values, 13, path), number(values, 14, path)});
    }
    return result;
}

class Recorder final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        states.push_back({step.time, problem.committed_solution(),
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, problem.committed_solution())});
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

std::size_t local_node(const fuelsim::Hex8RegionMesh& region, std::size_t source) {
    const auto found = std::find(region.source_node_ids().begin(), region.source_node_ids().end(), source);
    if (found == region.source_node_ids().end()) throw std::logic_error("B4.3 local contact node is missing");
    return static_cast<std::size_t>(found - region.source_node_ids().begin());
}

fuelsim::Quad4SurfaceContactLocalValues contact_state(const fuelsim::spatial_detail::SpatialLayout& spatial,
    const fuelsim::Hex8RegionMesh& primary, const fuelsim::Hex8RegionMesh& secondary,
    const std::vector<double>& state) {
    constexpr std::array<std::size_t, 4> secondary_sources = {8, 12, 15, 11}, primary_sources = {1, 2, 6, 5};
    fuelsim::Quad4SurfaceContactLocalValues result{};
    for (std::size_t node = 0; node < 4; ++node) {
        const std::size_t secondary_global = spatial.global_node(1, local_node(secondary, secondary_sources[node])),
                          primary_global = spatial.global_node(0, local_node(primary, primary_sources[node]));
        for (std::size_t component = 0; component < 3; ++component) {
            const fuelsim::Field field = component == 0   ? fuelsim::Field::displacement_x
                                         : component == 1 ? fuelsim::Field::displacement_y
                                                          : fuelsim::Field::displacement_z;
            const std::size_t offset = 8 * (component + 1);
            result[offset + node] = state[spatial.dof(field, secondary_global)];
            result[offset + 4 + node] = state[spatial.dof(field, primary_global)];
        }
    }
    return result;
}

std::array<double, 2> constraint_parent_coordinate(std::size_t source) {
    if (source == 8) return {-0.5, -0.5};
    if (source == 12) return {0.5, -0.5};
    if (source == 15) return {0.5, 0.5};
    if (source == 11) return {-0.5, 0.5};
    throw std::logic_error("B4.3 contact constraint source node is unexpected");
}

std::array<double, 3> cross(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return {first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

double dot(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

bool compare(const fuelsim::TransientProblem& problem, const fuelsim::UnstructuredHex8Mesh& mesh,
    const std::vector<StepState>& states, const std::vector<NodeReference>& nodes,
    const std::vector<ContactReference>& contacts) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    fuelsim::test::FieldErrorMetrics displacement_x, displacement_y, displacement_z, normal_x, normal_y, normal_z, gap,
        pressure, tangent_first_force, tangent_second_force, total_slip, resultant_x, resultant_y, resultant_z,
        moment_x, moment_y, moment_z, force_center_x, force_center_y, force_center_z;
    double maximum_coordinate_difference = 0.0, maximum_current_coordinate_difference = 0.0, final_area_change = 0.0,
           final_normal_rotation = 0.0;
    const std::vector<std::size_t> contact_sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    const fuelsim::Hex8RegionMesh& primary_region = spatial.region_mesh(0);
    const fuelsim::Hex8RegionMesh& secondary_region = spatial.region_mesh(1);
    constexpr std::array<std::size_t, 4> secondary_face_sources = {8, 12, 15, 11}, primary_face_sources = {1, 2, 6, 5};
    fuelsim::Quad4FaceCoordinates secondary_face{}, primary_face{};
    for (std::size_t node = 0; node < 4; ++node) {
        secondary_face[node] = mesh.nodes().at(secondary_face_sources[node]);
        primary_face[node] = mesh.nodes().at(primary_face_sources[node]);
    }
    std::vector<fuelsim::ContactPointHistory> slip_histories(contact_sources.size());
    std::vector<double> committed_state(spatial.dof_count(), 0.0);
    bool all_projected_active_and_sticking = true;
    for (std::size_t step = 0; step < states.size(); ++step) {
        const fuelsim::Quad4SurfaceContactLocalValues current_local = contact_state(spatial, primary_region,
                                                          secondary_region, states[step].solution),
                                                      committed_local = contact_state(
                                                          spatial, primary_region, secondary_region, committed_state);
        std::array<double, 3> actual_resultant{}, reference_resultant{}, actual_moment{}, reference_moment{},
            actual_center_sum{}, reference_center_sum{};
        double actual_center_weight = 0.0, reference_center_weight = 0.0;
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const fuelsim::Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
            for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
                const std::size_t source = region_mesh.source_node_ids()[local];
                const auto found = std::find_if(nodes.begin(), nodes.end(),
                    [&](const NodeReference& value) { return value.step == step + 1 && value.source_node == source; });
                if (found == nodes.end()) throw std::invalid_argument("B4.3 node mapping is incomplete");
                const fuelsim::CartesianPoint3& point = mesh.nodes()[source];
                maximum_coordinate_difference =
                    std::max({maximum_coordinate_difference, std::abs(point.x - found->reference.x),
                        std::abs(point.y - found->reference.y), std::abs(point.z - found->reference.z)});
                const std::size_t global = spatial.global_node(region, local);
                const std::array<std::size_t, 3> dofs = {spatial.dof(fuelsim::Field::displacement_x, global),
                    spatial.dof(fuelsim::Field::displacement_y, global),
                    spatial.dof(fuelsim::Field::displacement_z, global)};
                displacement_x.add(states[step].solution[dofs[0]], found->displacement[0]);
                displacement_y.add(states[step].solution[dofs[1]], found->displacement[1]);
                displacement_z.add(states[step].solution[dofs[2]], found->displacement[2]);
            }
        }
        for (std::size_t node = 0; node < states[step].contact.size(); ++node) {
            const auto found = std::find_if(contacts.begin(), contacts.end(), [&](const ContactReference& value) {
                return value.step == step + 1 && value.source_node == contact_sources[node];
            });
            if (found == contacts.end()) throw std::invalid_argument("B4.3 contact mapping is incomplete");
            const auto& actual = states[step].contact[node];
            const std::array<double, 2> parent = constraint_parent_coordinate(contact_sources[node]);
            const fuelsim::Quad4FaceQuadraturePoint point =
                fuelsim::make_quad4_face_quadrature_point(secondary_face, parent[0], parent[1], 1.0);
            const fuelsim::Quad4ToQuad4MechanicalGeometry geometry{
                secondary_face, primary_face, point.shape, point.derivative_xi, point.derivative_eta, 1.0, -1.0};
            const fuelsim::CartesianContactPointValue local_value = fuelsim::compute_quad4_to_quad4_contact_value(
                {1.0e8, 1.0e6, false, 0.0}, geometry, current_local, committed_local, slip_histories[node]);
            all_projected_active_and_sticking = all_projected_active_and_sticking && actual.projected &&
                                                actual.pressure > 0.0 && !actual.sliding && found->pressure > 0.0 &&
                                                local_value.projected && !local_value.sliding;
            const std::size_t secondary_local = local_node(secondary_region, contact_sources[node]),
                              secondary_global = spatial.global_node(1, secondary_local);
            const std::array<double, 3> actual_point = {
                mesh.nodes()[contact_sources[node]].x +
                    states[step].solution[spatial.dof(fuelsim::Field::displacement_x, secondary_global)],
                mesh.nodes()[contact_sources[node]].y +
                    states[step].solution[spatial.dof(fuelsim::Field::displacement_y, secondary_global)],
                mesh.nodes()[contact_sources[node]].z +
                    states[step].solution[spatial.dof(fuelsim::Field::displacement_z, secondary_global)]};
            const std::array<double, 3> reference_point = {found->current.x, found->current.y, found->current.z};
            for (std::size_t component = 0; component < 3; ++component)
                maximum_current_coordinate_difference = std::max(maximum_current_coordinate_difference,
                    std::abs(actual_point[component] - reference_point[component]));
            std::array<double, 3> actual_total_force{}, reference_total_force{};
            normal_x.add(-actual.normal_contact_force[0], found->normal_force[0]);
            normal_y.add(-actual.normal_contact_force[1], found->normal_force[1]);
            normal_z.add(-actual.normal_contact_force[2], found->normal_force[2]);
            gap.add(actual.gap, found->gap);
            pressure.add(actual.pressure, found->pressure);
            for (std::size_t component = 0; component < 3; ++component) {
                actual_total_force[component] =
                    -actual.normal_contact_force[component] - actual.tangential_contact_force[component];
                reference_total_force[component] = found->normal_force[component] + found->tangential_force[component];
                actual_resultant[component] += actual_total_force[component];
                reference_resultant[component] += reference_total_force[component];
            }
            const std::array<double, 3> tangent_second = cross(local_value.normal, local_value.tangent_first);
            std::array<double, 3> actual_tangential_force{}, reference_tangential_force{};
            for (std::size_t component = 0; component < 3; ++component) {
                actual_tangential_force[component] = -actual.tangential_contact_force[component];
                reference_tangential_force[component] = found->tangential_force[component];
            }
            tangent_first_force.add(dot(actual_tangential_force, local_value.tangent_first),
                dot(reference_tangential_force, local_value.tangent_first));
            tangent_second_force.add(
                dot(actual_tangential_force, tangent_second), dot(reference_tangential_force, tangent_second));
            total_slip.add(std::sqrt(dot(local_value.elastic_tangential_slip, local_value.elastic_tangential_slip)),
                std::hypot(found->slip_1, found->slip_2));
            const std::array<double, 3> actual_node_moment = cross(actual_point, actual_total_force),
                                        reference_node_moment = cross(reference_point, reference_total_force);
            const double actual_weight = std::sqrt(dot(actual_total_force, actual_total_force)),
                         reference_weight = std::sqrt(dot(reference_total_force, reference_total_force));
            actual_center_weight += actual_weight;
            reference_center_weight += reference_weight;
            for (std::size_t component = 0; component < 3; ++component) {
                actual_moment[component] += actual_node_moment[component];
                reference_moment[component] += reference_node_moment[component];
                actual_center_sum[component] += actual_weight * actual_point[component];
                reference_center_sum[component] += reference_weight * reference_point[component];
            }
            slip_histories[node].sliding = local_value.sliding;
            slip_histories[node].cartesian_elastic_tangential_slip = local_value.elastic_tangential_slip;
            slip_histories[node].cartesian_tangent_basis_initialized = true;
            slip_histories[node].cartesian_contact_normal = local_value.normal;
            slip_histories[node].cartesian_contact_tangent_first = local_value.tangent_first;
        }
        resultant_x.add(actual_resultant[0], reference_resultant[0]);
        resultant_y.add(actual_resultant[1], reference_resultant[1]);
        resultant_z.add(actual_resultant[2], reference_resultant[2]);
        moment_x.add(actual_moment[0], reference_moment[0]);
        moment_y.add(actual_moment[1], reference_moment[1]);
        moment_z.add(actual_moment[2], reference_moment[2]);
        force_center_x.add(
            actual_center_sum[0] / actual_center_weight, reference_center_sum[0] / reference_center_weight);
        force_center_y.add(
            actual_center_sum[1] / actual_center_weight, reference_center_sum[1] / reference_center_weight);
        force_center_z.add(
            actual_center_sum[2] / actual_center_weight, reference_center_sum[2] / reference_center_weight);
        if (step + 1 == states.size()) {
            double area = 0.0;
            std::array<double, 3> normal_sum{};
            for (const auto& contact : states[step].contact) {
                area += contact.tributary_area;
                for (std::size_t component = 0; component < 3; ++component)
                    normal_sum[component] -= contact.normal_contact_force[component];
            }
            final_area_change = std::abs(area - 1.0);
            final_normal_rotation = std::hypot(normal_sum[1], normal_sum[2]) / std::abs(normal_sum[0]);
        }
        committed_state = states[step].solution;
    }
    const auto print_metric = [](const std::string& name, const fuelsim::test::FieldErrorMetrics& value) {
        if (value.has_relative_norm())
            fuelsim::test::print_relative_metrics(name, value);
        else {
            fuelsim::test::print_absolute_metrics(name, value);
            std::cout << name << "_zero_reference_count=" << value.zero_reference_count << '\n'
                      << name
                      << "_maximum_zero_reference_absolute_difference=" << value.maximum_zero_reference_difference
                      << '\n';
        }
    };
    print_metric("b43_displacement_x", displacement_x);
    print_metric("b43_displacement_y", displacement_y);
    print_metric("b43_displacement_z", displacement_z);
    print_metric("b43_normal_force_x", normal_x);
    print_metric("b43_normal_force_y", normal_y);
    print_metric("b43_normal_force_z", normal_z);
    print_metric("b43_local_tangent_first_force", tangent_first_force);
    print_metric("b43_local_tangent_second_force", tangent_second_force);
    print_metric("b43_total_slip", total_slip);
    print_metric("b43_gap", gap);
    print_metric("b43_pressure", pressure);
    print_metric("b43_resultant_x", resultant_x);
    print_metric("b43_resultant_y", resultant_y);
    print_metric("b43_resultant_z", resultant_z);
    print_metric("b43_moment_x", moment_x);
    print_metric("b43_moment_y", moment_y);
    print_metric("b43_moment_z", moment_z);
    print_metric("b43_force_center_x", force_center_x);
    print_metric("b43_force_center_y", force_center_y);
    print_metric("b43_force_center_z", force_center_z);
    std::cout << "b43_maximum_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << "b43_maximum_current_coordinate_difference=" << maximum_current_coordinate_difference << '\n'
              << "b43_final_secondary_area_change=" << final_area_change << '\n'
              << "b43_final_normal_rotation_ratio=" << final_normal_rotation << '\n';
    constexpr double tolerance = 1.0e-2, zero_tolerance = 1.0e-7;
    const auto field_passes = [&](const fuelsim::test::FieldErrorMetrics& value) {
        return (!value.has_relative_norm() || fuelsim::test::relative_metrics_below(value, tolerance)) &&
               value.maximum_zero_reference_difference < zero_tolerance;
    };
    return check(maximum_coordinate_difference < 1.0e-14, "B4.3 Abaqus and Fuelsim use the same tracked mesh") &&
           check(maximum_current_coordinate_difference < 1.0e-14,
               "B4.3 Abaqus and Fuelsim use the same current contact coordinates at every increment") &&
           check(final_area_change > 0.02 && final_normal_rotation > 0.01,
               "B4.3 produces material current-area change and contact-normal rotation") &&
           check(all_projected_active_and_sticking,
               "B4.3 keeps all four constraints projected, active, and frictionlessly nonsliding") &&
           check(field_passes(displacement_x) && field_passes(displacement_y) && field_passes(displacement_z),
               "B4.3 all three full-mesh displacement fields agree with Abaqus below 1 percent") &&
           check(field_passes(normal_x) && field_passes(normal_y) && field_passes(normal_z),
               "B4.3 all three secondary nodal normal-force fields agree with Abaqus below 1 percent") &&
           check(field_passes(tangent_first_force) && field_passes(tangent_second_force),
               "B4.3 both local tangent-plane force fields pass their separate theoretical-zero checks") &&
           check(field_passes(total_slip), "B4.3 cumulative total slip agrees with Abaqus below 1 percent") &&
           check(field_passes(gap) && field_passes(pressure),
               "B4.3 gap and pressure fields agree with Abaqus below 1 percent") &&
           check(field_passes(resultant_x) && field_passes(resultant_y) && field_passes(resultant_z),
               "B4.3 complete contact resultants agree with Abaqus below 1 percent") &&
           check(field_passes(moment_x) && field_passes(moment_y) && field_passes(moment_z),
               "B4.3 complete contact moments agree with Abaqus below 1 percent") &&
           check(field_passes(force_center_x) && field_passes(force_center_y) && field_passes(force_center_z),
               "B4.3 contact force centers agree with Abaqus below 1 percent");
}

bool run(const std::string& input_path, const std::string& node_path, const std::string& contact_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(input.mesh_file);
    fuelsim::TransientProblem problem(input.spatial, mesh);
    Recorder recorder;
    const fuelsim::TransientResult result =
        fuelsim::solve_transient(problem, input.transient_execution, solver_options(input), &recorder);
    bool passed = check(result.completed && recorder.states.size() == 20 && result.rejected_steps.empty(),
        "B4.3 completes the same twenty one-unit increments without rejected steps");
    passed = compare(problem, mesh, recorder.states, read_nodes(node_path), read_contact(contact_path)) && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--generate") {
        try {
            generate(argv[2], argv[3]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] B4.3 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_b43_hex8_finite_strain_contact_abaqus_tests <case.fsi> <nodes.csv> "
                     "<contact.csv>\n"
                     "   or: fuelsim_b43_hex8_finite_strain_contact_abaqus_tests --generate <abaqus-directory> "
                     "<fuelsim-directory>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B4.3 HEX8 finite-strain Abaqus contact comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3]);
        if (passed && session.rank() == 0) std::cout << "[PASS] B4.3 HEX8 finite-strain contact Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.3 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
