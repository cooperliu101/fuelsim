#include "fuelsim/core/cartesian3d_hex20.hpp"
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
    bool quadratic_surface;
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

struct GeometryEvidence final {
    double minimum_contact_radius = std::numeric_limits<double>::max();
    double maximum_contact_radius_error = 0.0;
    double maximum_face_planarity_error = 0.0;
    double minimum_face_normal_dot = 1.0;
    double minimum_within_face_normal_dot = 1.0;
};

struct ContactNumericalEvidence final {
    std::array<double, 3> residual_balance{};
    double maximum_jacobian_directional_error = 0.0;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<CaseSpec> case_specs() { return {{"faceted_cylinder", false}, {"quadratic_cylinder", true}}; }

fuelsim::CartesianPoint3 cylindrical_point(double radius, double angle, double z) {
    return {radius * std::cos(angle), radius * std::sin(angle), z};
}

fuelsim::CartesianPoint3 subtract(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

fuelsim::CartesianPoint3 cross(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y, first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

double dot(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

fuelsim::CartesianPoint3 unit(const fuelsim::CartesianPoint3& value) {
    const double measure = std::sqrt(dot(value, value));
    if (!(measure > 0.0)) throw std::invalid_argument("H20.30 geometry has an undefined surface normal");
    return {value.x / measure, value.y / measure, value.z / measure};
}

fuelsim::Hex20Element append_annular_element(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map, double radius_lower, double radius_upper,
    double angle_lower, double angle_upper, double z_lower, double z_upper, bool quadratic_surface) {
    const std::array<std::array<double, 3>, 8> logical_corners = {
        {{radius_lower, angle_lower, z_lower}, {radius_upper, angle_lower, z_lower},
            {radius_upper, angle_upper, z_lower}, {radius_lower, angle_upper, z_lower},
            {radius_lower, angle_lower, z_upper}, {radius_upper, angle_lower, z_upper},
            {radius_upper, angle_upper, z_upper}, {radius_lower, angle_upper, z_upper}}};
    const std::array<std::pair<std::size_t, std::size_t>, 12> edges = {
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}}};
    std::array<std::array<double, 3>, 20> logical_points{};
    std::copy(logical_corners.begin(), logical_corners.end(), logical_points.begin());
    for (std::size_t edge = 0; edge < edges.size(); ++edge)
        for (std::size_t component = 0; component < 3; ++component)
            logical_points[8 + edge][component] =
                0.5 * (logical_corners[edges[edge].first][component] + logical_corners[edges[edge].second][component]);

    fuelsim::Hex20Element element{};
    for (std::size_t local = 0; local < logical_points.size(); ++local) {
        const auto inserted = node_map.emplace(logical_points[local], nodes.size());
        if (inserted.second) {
            fuelsim::CartesianPoint3 point{};
            if (local < logical_corners.size() || quadratic_surface) {
                point = cylindrical_point(logical_points[local][0], logical_points[local][1], logical_points[local][2]);
            } else {
                const auto edge = edges[local - logical_corners.size()];
                const fuelsim::CartesianPoint3 first = cylindrical_point(
                    logical_corners[edge.first][0], logical_corners[edge.first][1], logical_corners[edge.first][2]);
                const fuelsim::CartesianPoint3 second = cylindrical_point(
                    logical_corners[edge.second][0], logical_corners[edge.second][1], logical_corners[edge.second][2]);
                point = {0.5 * (first.x + second.x), 0.5 * (first.y + second.y), 0.5 * (first.z + second.z)};
            }
            nodes.push_back(point);
        }
        element.nodes[local] = inserted.first->second;
    }
    return element;
}

std::vector<std::size_t> side_nodes(
    const std::vector<fuelsim::Hex20Element>& elements, const std::vector<fuelsim::ElementSide>& sides) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {
        {{{0, 1, 5, 4, 8, 13, 16, 12}}, {{1, 2, 6, 5, 9, 14, 17, 13}}, {{2, 3, 7, 6, 10, 15, 18, 14}},
            {{3, 0, 4, 7, 11, 12, 19, 15}}, {{0, 3, 2, 1, 11, 10, 9, 8}}, {{4, 5, 6, 7, 16, 17, 18, 19}}}};
    std::vector<std::size_t> result;
    for (const fuelsim::ElementSide& side : sides)
        for (std::size_t local : face_nodes.at(side.local_side))
            result.push_back(elements.at(side.element).nodes[local]);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

fuelsim::UnstructuredHex20Mesh generate_mesh(const CaseSpec& spec) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const std::array<double, 5> angles = {0.0, pi / 8.0, pi / 4.0, 3.0 * pi / 8.0, pi / 2.0};
    const std::array<double, 3> axial = {0.0, 0.5, 1.0};
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    std::vector<fuelsim::Hex20Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, primary_inner, primary_back, secondary_contact, secondary_outer,
        secondary_back, secondary_theta_lower, secondary_theta_upper;
    for (std::size_t z = 0; z + 1 < axial.size(); ++z)
        for (std::size_t angle = 0; angle + 1 < angles.size(); ++angle) {
            const std::size_t element = elements.size();
            elements.push_back(append_annular_element(nodes, primary_nodes, 0.9, 1.0, angles[angle], angles[angle + 1],
                axial[z], axial[z + 1], spec.quadratic_surface));
            blocks.push_back(1);
            primary_contact.push_back({element, 1});
            primary_inner.push_back({element, 3});
            if (z == 0) primary_back.push_back({element, 4});
        }
    for (std::size_t z = 0; z + 1 < axial.size(); ++z)
        for (std::size_t angle = 0; angle + 1 < angles.size(); ++angle) {
            const std::size_t element = elements.size();
            elements.push_back(append_annular_element(nodes, secondary_nodes, 1.0, 1.1, angles[angle],
                angles[angle + 1], axial[z], axial[z + 1], spec.quadratic_surface));
            blocks.push_back(2);
            secondary_contact.push_back({element, 3});
            secondary_outer.push_back({element, 1});
            if (z == 0) secondary_back.push_back({element, 4});
            if (angle == 0) secondary_theta_lower.push_back({element, 0});
            if (angle + 2 == angles.size()) secondary_theta_upper.push_back({element, 2});
        }

    const std::vector<fuelsim::NodeSet> node_sets = {{10, "primary_inner", side_nodes(elements, primary_inner)},
        {11, "primary_back", side_nodes(elements, primary_back)},
        {12, "secondary_contact_nodes", side_nodes(elements, secondary_contact)},
        {13, "secondary_outer", side_nodes(elements, secondary_outer)},
        {14, "secondary_back", side_nodes(elements, secondary_back)},
        {15, "secondary_theta_lower", side_nodes(elements, secondary_theta_lower)},
        {16, "secondary_theta_upper", side_nodes(elements, secondary_theta_upper)}};
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "primary"}, {2, "secondary"}}, node_sets,
        {{20, "primary_contact", std::move(primary_contact)}, {21, "primary_inner", std::move(primary_inner)},
            {22, "primary_back", std::move(primary_back)}, {30, "secondary_contact", std::move(secondary_contact)},
            {31, "secondary_outer", std::move(secondary_outer)}, {32, "secondary_back", std::move(secondary_back)},
            {33, "secondary_theta_lower", std::move(secondary_theta_lower)},
            {34, "secondary_theta_upper", std::move(secondary_theta_upper)}});
}

GeometryEvidence inspect_contact_geometry(const fuelsim::UnstructuredHex20Mesh& mesh) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {
        {{{0, 1, 5, 4, 8, 13, 16, 12}}, {{1, 2, 6, 5, 9, 14, 17, 13}}, {{2, 3, 7, 6, 10, 15, 18, 14}},
            {{3, 0, 4, 7, 11, 12, 19, 15}}, {{0, 3, 2, 1, 11, 10, 9, 8}}, {{4, 5, 6, 7, 16, 17, 18, 19}}}};
    GeometryEvidence result;
    std::vector<fuelsim::CartesianPoint3> face_normals;
    for (const fuelsim::ElementSide& side : mesh.side_set("secondary_contact").sides) {
        fuelsim::Quad8FaceCoordinates coordinates{};
        for (std::size_t local = 0; local < coordinates.size(); ++local) {
            coordinates[local] = mesh.nodes()[mesh.elements()[side.element].nodes[face_nodes[side.local_side][local]]];
            const double radius = std::hypot(coordinates[local].x, coordinates[local].y);
            result.minimum_contact_radius = std::min(result.minimum_contact_radius, radius);
            result.maximum_contact_radius_error = std::max(result.maximum_contact_radius_error, std::abs(radius - 1.0));
        }
        const fuelsim::CartesianPoint3 plane_normal =
            unit(cross(subtract(coordinates[1], coordinates[0]), subtract(coordinates[3], coordinates[0])));
        for (const fuelsim::CartesianPoint3& point : coordinates)
            result.maximum_face_planarity_error = std::max(
                result.maximum_face_planarity_error, std::abs(dot(subtract(point, coordinates[0]), plane_normal)));
        const fuelsim::Quad8FaceMechanicalQuadraturePoint center =
            fuelsim::make_quad8_face_mechanical_point(coordinates, 0.0, 0.0, 1.0);
        face_normals.push_back(unit(cross(center.tangent_xi, center.tangent_eta)));
        const fuelsim::Quad8FaceMechanicalQuadraturePoint first =
            fuelsim::make_quad8_face_mechanical_point(coordinates, -0.75, -0.75, 1.0);
        const fuelsim::Quad8FaceMechanicalQuadraturePoint second =
            fuelsim::make_quad8_face_mechanical_point(coordinates, 0.75, -0.75, 1.0);
        result.minimum_within_face_normal_dot = std::min(
            result.minimum_within_face_normal_dot, std::abs(dot(unit(cross(first.tangent_xi, first.tangent_eta)),
                                                       unit(cross(second.tangent_xi, second.tangent_eta)))));
    }
    for (std::size_t first = 0; first < face_normals.size(); ++first)
        for (std::size_t second = first + 1; second < face_normals.size(); ++second)
            result.minimum_face_normal_dot =
                std::min(result.minimum_face_normal_dot, std::abs(dot(face_normals[first], face_normals[second])));
    return result;
}

ContactNumericalEvidence inspect_contact_numerics(const fuelsim::SteadyProblem& problem,
    const fuelsim::cartesian::SpatialAssembly& spatial, const std::vector<double>& state) {
    ContactNumericalEvidence result;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local) local_state[local] = state[dofs[local]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local_state.size()), plus = local_state, minus = local_state;
        constexpr double perturbation = 1.0e-10;
        for (std::size_t local = 0; local < local_state.size(); ++local) {
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] += perturbation * direction[local];
            minus[local] -= perturbation * direction[local];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local_state.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[row * local_state.size() + column] * direction[column];
            const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - finite_difference, 2);
            reference_squared += finite_difference * finite_difference;
        }
        if (!(reference_squared > 0.0)) throw std::logic_error("H20.30 active contact has a zero Jacobian direction");
        result.maximum_jacobian_directional_error =
            std::max(result.maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
        for (std::size_t local = 0; local < dofs.size(); ++local)
            for (std::size_t component = 0; component < result.residual_balance.size(); ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (dofs[local] >= field.begin && dofs[local] < field.end)
                    result.residual_balance[component] += residual[local];
            }
    }
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

void write_abaqus_input(const std::string& path, const CaseSpec& spec, const fuelsim::UnstructuredHex20Mesh& mesh) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write H20.30 Abaqus input: " + path);
    output << std::setprecision(16) << "*Heading\n"
           << "** H20.30: " << spec.name << " curved C3D20 small-sliding surface-to-surface contact.\n"
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
            for (std::size_t local : abaqus_order) output << ", " << mesh.elements()[element].nodes[local] + 1;
            output << '\n';
        }
    };
    write_block("PRIMARY", 1);
    write_block("SECONDARY", 2);
    for (const fuelsim::NodeSet& set : mesh.node_sets()) {
        output << "*Nset, nset=" << set.name << '\n';
        write_label_set(output, set.nodes);
    }
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_OUTER\nSECONDARY, S4\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.25\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=PENALTY_CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e11,\n"
           << "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
           << "*Step, name=LOAD, nlgeom=NO, inc=100\n*Static\n0.1, 1., 1.e-8, 0.1\n"
           << "*Boundary\n"
           << "primary_inner, 1, 2, 0.\nprimary_back, 3, 3, 0.\nsecondary_back, 3, 3, 0.\n"
           << "secondary_theta_lower, 2, 2, 0.\nsecondary_theta_upper, 1, 1, 0.\n"
           << "*Dsload\nSECONDARY_OUTER, P, 10000.\n"
           << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
           << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
}

void write_material(std::ofstream& output, const std::string& name) {
    output << "  [" << name << "]\n"
           << "    [thermal]\n      function = constant_thermophysical\n      conductivity = 10\n"
           << "      density = 1\n      specific_heat = 1\n    []\n"
           << "    [elasticity]\n      function = constant_isotropic\n      young_modulus = 1e9\n"
           << "      poisson_ratio = 0.25\n    []\n  []\n";
}

void write_fsi(const std::string& path, const std::string& mesh_relative_path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write H20.30 Fuelsim input: " + path);
    output << "[Case]\n  version = 3\n  problem = steady\n  geometry = cartesian_3d\n[]\n\n"
           << "[Mesh]\n  type = exodus\n  file = " << mesh_relative_path << "\n[]\n\n[Materials]\n";
    write_material(output, "primary");
    write_material(output, "secondary");
    output << "[]\n\n[Regions]\n"
           << "  [primary]\n    block = primary\n    material = primary\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n"
           << "  [secondary]\n    block = secondary\n    material = secondary\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n[]\n\n"
           << "[Contact]\n  [interface]\n    primary = primary_contact\n    secondary = secondary_contact\n"
           << "    [mechanical]\n      formulation = penalty\n      discretization = surface_to_surface\n"
           << "      penalty = 1e11\n      mu = 0\n    []\n  []\n[]\n\n"
           << "[BoundaryConditions]\n"
           << "  [primary_temperature]\n    type = dirichlet\n    boundary = primary_inner\n"
           << "    field = temperature\n    value = 300\n  []\n"
           << "  [secondary_temperature]\n    type = dirichlet\n    boundary = secondary_outer\n"
           << "    field = temperature\n    value = 300\n  []\n"
           << "  [primary_fix_x]\n    type = dirichlet\n    boundary = primary_inner\n"
           << "    field = displacement_x\n    value = 0\n  []\n"
           << "  [primary_fix_y]\n    type = dirichlet\n    boundary = primary_inner\n"
           << "    field = displacement_y\n    value = 0\n  []\n"
           << "  [primary_fix_z]\n    type = dirichlet\n    boundary = primary_back\n"
           << "    field = displacement_z\n    value = 0\n  []\n"
           << "  [secondary_fix_z]\n    type = dirichlet\n    boundary = secondary_back\n"
           << "    field = displacement_z\n    value = 0\n  []\n"
           << "  [secondary_symmetry_lower]\n    type = dirichlet\n    boundary = secondary_theta_lower\n"
           << "    field = displacement_y\n    value = 0\n  []\n"
           << "  [secondary_symmetry_upper]\n    type = dirichlet\n    boundary = secondary_theta_upper\n"
           << "    field = displacement_x\n    value = 0\n  []\n"
           << "  [outer_pressure]\n    type = pressure\n    boundary = secondary_outer\n"
           << "    value = 10000\n    configuration = reference\n    scale_with_load = true\n  []\n"
           << "[]\n\n[Executioner]\n  type = steady\n  load_steps = 4\n[]\n\n"
           << "[Solver]\n  absolute_tolerance = 1e-8\n  relative_tolerance = 1e-11\n"
           << "  step_tolerance = 1e-12\n  maximum_iterations = 30\n  linear_solver = direct\n"
           << "  direct_factorization = mumps\n  field_residual_scaling = true\n"
           << "  linear_relative_tolerance = 1e-11\n  maximum_linear_iterations = 400\n[]\n\n"
           << "[Outputs]\n  console = true\n[]\n";
}

void generate_references(const std::string& abaqus_directory, const std::string& fuelsim_directory) {
    for (const CaseSpec& spec : case_specs()) {
        const fuelsim::UnstructuredHex20Mesh mesh = generate_mesh(spec);
        const std::string base = "h20_30_hex20_sts_" + spec.name;
        fuelsim::write_exodus_hex20(abaqus_directory + '/' + base + "_mesh.e", mesh);
        write_abaqus_input(abaqus_directory + '/' + base + ".inp", spec, mesh);
        write_fsi(
            fuelsim_directory + "/steady_hex20_sts_" + spec.name + "_abaqus.fsi", "../abaqus/" + base + "_mesh.e");
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
    if (index >= values.size()) throw std::invalid_argument("Incomplete H20.30 CSV row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.30 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_displacement,normal_x,normal_y,normal_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.30 displacement header in " + path);
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
    if (!input) throw std::runtime_error("Could not read H20.30 nodal-force reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "normal_force,id,x,y,z") throw std::invalid_argument("Unexpected H20.30 force header in " + path);
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
    if (!input) throw std::runtime_error("Could not read H20.30 resultant reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "radial_contact_resultant,radial_primary_reaction")
        throw std::invalid_argument("Unexpected H20.30 resultant header in " + path);
    if (!std::getline(input, line)) throw std::invalid_argument("Missing H20.30 resultant value in " + path);
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
    const std::string& force_path, const std::string& resultant_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(case_path);
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    const GeometryEvidence geometry = inspect_contact_geometry(mesh);
    fuelsim::SteadyProblem problem(definition.spatial, mesh);
    const fuelsim::SteadyResult solve =
        fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
    bool passed = check(
        solve.completed && solve.solve.converged && solve.completed_steps == definition.steady_execution.load_steps,
        "H20.30 " + name + " Fuelsim solve completes every load step");
    const std::vector<DisplacementReference> displacement = read_displacements(displacement_path);
    const std::vector<ForceReference> force = read_forces(force_path);
    if (displacement.empty() || force.empty()) throw std::invalid_argument("H20.30 reference fields are empty");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const ContactNumericalEvidence numerical = inspect_contact_numerics(problem, spatial, solve.solve.state);
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
                throw std::invalid_argument("H20.30 displacement node mapping is incomplete or repeated");
            present[source] = true;
            maximum_coordinate_error =
                std::max(maximum_coordinate_error, coordinate_difference(mesh.nodes()[source], reference->point));
            const std::size_t global = spatial.global_node(region, local);
            double actual = 0.0;
            for (std::size_t component = 0; component < 3; ++component)
                actual += reference->normal[component] * solve.solve.state[fields[component + 1].begin + global];
            displacement_error.add(actual, reference->displacement);
        }
    }
    const auto summaries = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    double maximum_force_coordinate_error = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const auto reference = std::find_if(
            force.begin(), force.end(), [&](const ForceReference& value) { return value.id == source_nodes[node]; });
        if (reference == force.end()) throw std::invalid_argument("H20.30 contact-force node mapping is incomplete");
        maximum_force_coordinate_error = std::max(
            maximum_force_coordinate_error, coordinate_difference(mesh.nodes()[reference->id], reference->point));
        const double radius = std::hypot(reference->point.x, reference->point.y);
        const std::array<double, 3> radial = {reference->point.x / radius, reference->point.y / radius, 0.0};
        double actual_force = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            actual_force += radial[component] * summaries[node].normal_contact_force[component];
        force_error.add(-actual_force, reference->force);
    }
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
    const std::array<double, 2> reference_resultants = read_resultants(resultant_path);
    double radial_resultant = 0.0;
    for (std::size_t node = 0; node < summaries.size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[source_nodes[node]];
        const double radius = std::hypot(point.x, point.y);
        radial_resultant -=
            (point.x * summaries[node].normal_contact_force[0] + point.y * summaries[node].normal_contact_force[1]) /
            radius;
    }
    const double resultant_error = std::abs(radial_resultant - reference_resultants[0]) / reference_resultants[0];
    fuelsim::test::print_relative_metrics("h20_30_" + name + "_radial_displacement", displacement_error);
    fuelsim::test::print_relative_metrics("h20_30_" + name + "_nodal_radial_force", force_error);
    std::cout << "h20_30_" << name << "_radial_contact_resultant=" << radial_resultant << '\n'
              << "h20_30_" << name << "_scalar_contact_force_sum=" << interface.total_contact_force << '\n'
              << "h20_30_" << name << "_abaqus_radial_contact_resultant=" << reference_resultants[0] << '\n'
              << "h20_30_" << name << "_abaqus_primary_radial_reaction=" << reference_resultants[1] << '\n'
              << "h20_30_" << name << "_radial_resultant_relative_error=" << resultant_error << '\n'
              << "h20_30_" << name << "_minimum_contact_node_radius=" << geometry.minimum_contact_radius << '\n'
              << "h20_30_" << name << "_maximum_contact_node_radius_error=" << geometry.maximum_contact_radius_error
              << '\n'
              << "h20_30_" << name << "_maximum_face_planarity_error=" << geometry.maximum_face_planarity_error << '\n'
              << "h20_30_" << name << "_minimum_between_face_normal_dot=" << geometry.minimum_face_normal_dot << '\n'
              << "h20_30_" << name << "_minimum_within_face_normal_dot=" << geometry.minimum_within_face_normal_dot
              << '\n'
              << "h20_30_" << name << "_contact_residual_balance=" << numerical.residual_balance[0] << ','
              << numerical.residual_balance[1] << ',' << numerical.residual_balance[2] << '\n'
              << "h20_30_" << name
              << "_contact_jacobian_directional_error=" << numerical.maximum_jacobian_directional_error << '\n';
    constexpr double tolerance = 1.0e-2;
    const bool faceted_geometry = geometry.minimum_contact_radius < 0.99 &&
                                  geometry.maximum_face_planarity_error < 1.0e-12 &&
                                  geometry.minimum_within_face_normal_dot > 1.0 - 1.0e-12;
    const bool quadratic_geometry = geometry.maximum_contact_radius_error < 1.0e-12 &&
                                    geometry.maximum_face_planarity_error > 1.0e-3 &&
                                    geometry.minimum_within_face_normal_dot < 0.99;
    return check(displacement_error.value_count == mesh.nodes().size() &&
                     std::all_of(present.begin(), present.end(), [](bool value) { return value; }),
               "H20.30 " + name + " compares every mesh-node radial displacement") &&
           check(force_error.value_count == summaries.size() && force.size() == summaries.size(),
               "H20.30 " + name + " compares every secondary nodal radial force") &&
           check(maximum_coordinate_error < 1.0e-12 && maximum_force_coordinate_error < 1.0e-12,
               "H20.30 " + name + " Abaqus references use the tracked Exodus coordinates") &&
           check(interface.active_contact_nodes == summaries.size() && interface.unprojected_contact_nodes == 0,
               "H20.30 " + name + " keeps every secondary contact constraint active and projected") &&
           check(geometry.minimum_face_normal_dot < 0.95,
               "H20.30 " + name + " changes the contact normal across the assembled cylinder") &&
           check(name == "faceted_cylinder" ? faceted_geometry : quadratic_geometry,
               "H20.30 " + name + " has the requested faceted or genuine quadratic contact geometry") &&
           check(std::abs(numerical.residual_balance[0]) < 1.0e-8 && std::abs(numerical.residual_balance[1]) < 1.0e-8 &&
                     std::abs(numerical.residual_balance[2]) < 1.0e-8,
               "H20.30 " + name + " contact residual is action-reaction conservative") &&
           check(numerical.maximum_jacobian_directional_error < 1.0e-7,
               "H20.30 " + name + " contact Jacobian matches a centered directional difference") &&
           check(fuelsim::test::relative_metrics_below(displacement_error, tolerance),
               "H20.30 " + name + " radial-displacement three Abaqus errors are below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(force_error, tolerance),
               "H20.30 " + name + " nodal-radial-force three Abaqus errors are below 1 percent") &&
           check(resultant_error < tolerance, "H20.30 " + name + " radial resultant agrees below 1 percent") && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--generate") {
        try {
            generate_references(argv[2], argv[3]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] H20.30 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    constexpr std::size_t arguments_per_case = 5;
    const std::size_t argument_count = argc > 0 ? static_cast<std::size_t>(argc - 1) : 0;
    if (argument_count < arguments_per_case || argument_count % arguments_per_case != 0) {
        std::cerr << "Usage: fuelsim_h20_30_hex20_curved_sts_abaqus_tests "
                     "<name> <case.fsi> <displacement.csv> <force.csv> <resultant.csv> [...]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.30 curved Abaqus comparison\n");
        bool passed = true;
        for (int argument = 1; argument < argc; argument += static_cast<int>(arguments_per_case))
            passed = compare_case(argv[argument], argv[argument + 1], argv[argument + 2], argv[argument + 3],
                         argv[argument + 4]) &&
                     passed;
        if (passed && session.rank() == 0) std::cout << "[PASS] H20.30 curved Abaqus STS comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.30 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
