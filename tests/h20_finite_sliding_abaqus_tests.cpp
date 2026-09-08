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
struct ReferenceNode final {
    std::size_t step, source_node;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> normal_force{}, tangential_force{};
    double gap, pressure;
};

struct PathStep final {
    const char* name;
    std::array<double, 3> displacement;
};

constexpr std::array<PathStep, 4> path = {{{"CLOSE", {-0.01, 0.05, 0.04}},
    {"CROSS", {-0.01, 1.0, 0.4}},
    {"STRADDLE", {-0.01, 0.2, -0.3}},
    {"RETURN", {-0.01, 0.0, 0.1}}}};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::Hex20Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map,
    double x0,
    double x1,
    double y0,
    double y1,
    double z0,
    double z1) {
    const std::array<fuelsim::CartesianPoint3, 8> corners = {{{x0, y0, z0},
        {x1, y0, z0},
        {x1, y1, z0},
        {x0, y1, z0},
        {x0, y0, z1},
        {x1, y0, z1},
        {x1, y1, z1},
        {x0, y1, z1}}};
    const std::array<std::pair<std::size_t, std::size_t>, 12> edges = {
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}}};
    std::array<fuelsim::CartesianPoint3, 20> points{};
    std::copy(corners.begin(), corners.end(), points.begin());
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const fuelsim::CartesianPoint3& first = corners[edges[edge].first];
        const fuelsim::CartesianPoint3& second = corners[edges[edge].second];
        points[8 + edge] = {0.5 * (first.x + second.x), 0.5 * (first.y + second.y), 0.5 * (first.z + second.z)};
    }
    fuelsim::Hex20Element element{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const std::array<double, 3> key = {points[local].x, points[local].y, points[local].z};
        const auto inserted = node_map.emplace(key, nodes.size());
        if (inserted.second)
            nodes.push_back(points[local]);
        element.nodes[local] = inserted.first->second;
    }
    return element;
}

fuelsim::UnstructuredHex20Mesh generate_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    const fuelsim::Hex20Element primary_lower = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 0.0, 1.0, -1.0, 2.0),
                                primary_upper = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 1.0, 2.0, -1.0, 2.0),
                                secondary = append_cuboid(nodes, secondary_nodes, 1.0, 2.0, 0.1, 0.9, 0.1, 0.9);
    std::vector<std::size_t> primary_all, secondary_all;
    for (const auto& entry : primary_nodes)
        primary_all.push_back(entry.second);
    for (const auto& entry : secondary_nodes)
        secondary_all.push_back(entry.second);
    std::sort(primary_all.begin(), primary_all.end());
    std::sort(secondary_all.begin(), secondary_all.end());
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes),
        {primary_lower, primary_upper, secondary},
        {1, 1, 2},
        {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all", std::move(primary_all)}, {20, "secondary_all", std::move(secondary_all)}},
        {{30, "primary_contact", {{0, 1}, {1, 1}}}, {40, "secondary_contact", {{2, 3}}}});
}

void write_labels(std::ofstream& output, const std::vector<std::size_t>& nodes) {
    for (std::size_t index = 0; index < nodes.size(); ++index) {
        output << nodes[index] + 1;
        if ((index + 1) % 16 == 0 || index + 1 == nodes.size())
            output << '\n';
        else
            output << ", ";
    }
}

void write_abaqus_input(const std::string& output_path, const fuelsim::UnstructuredHex20Mesh& mesh, bool nlgeom) {
    std::ofstream output(output_path);
    if (!output)
        throw std::runtime_error("Could not write H20.39 Abaqus input: " + output_path);
    output << std::setprecision(16) << "*Heading\n"
           << "** H20.39 HEX20 finite-sliding surface-to-surface cross-face friction path.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    constexpr std::array<std::size_t, 20> abaqus_order =
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
    const auto write_block = [&](const std::string& name, std::int64_t block) {
        output << "*Element, type=C3D20, elset=" << name << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != block)
                continue;
            output << element + 1;
            for (const std::size_t local : abaqus_order)
                output << ", " << mesh.elements()[element].nodes[local] + 1;
            output << '\n';
        }
    };
    write_block("PRIMARY", 1);
    write_block("SECONDARY", 2);
    for (const fuelsim::NodeSet& set : mesh.node_sets()) {
        output << "*Nset, nset=" << set.name << '\n';
        write_labels(output, set.nodes);
    }
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\n"
           << "PRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\n"
           << "SECONDARY, S6\n"
           << "*Material, name=ELASTIC\n"
           << "*Elastic\n1.e9, 0.25\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e5,\n"
           << "*Friction, slip tolerance=1.e-4\n0.5,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n";
    for (std::size_t step = 0; step < path.size(); ++step) {
        output << "*Step, name=" << path[step].name << ", nlgeom=" << (nlgeom ? "YES" : "NO") << ", inc=100\n"
               << "*Static\n0.1, 1., 1.e-8, 0.1\n"
               << (step == 0 ? "*Boundary\nprimary_all, 1, 3, 0.\n" : "*Boundary, op=MOD\n");
        for (std::size_t component = 0; component < 3; ++component)
            output << "secondary_all, " << component + 1 << ", " << component + 1 << ", "
                   << path[step].displacement[component] << '\n';
        output << "*Output, field, frequency=1\n"
               << "*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n"
               << "*End Step\n";
    }
}

void generate(const std::string& directory) {
    const fuelsim::UnstructuredHex20Mesh mesh = generate_mesh();
    fuelsim::write_exodus_hex20(directory + "/h20_39_hex20_finite_sliding_mesh.e", mesh);
    write_abaqus_input(directory + "/h20_39_hex20_finite_sliding_small.inp", mesh, false);
    write_abaqus_input(directory + "/h20_39_hex20_finite_sliding_finite.inp", mesh, true);
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    if (column >= values.size())
        throw std::invalid_argument("Incomplete H20.39 CSV row: " + input_path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid H20.39 CSV number: " + input_path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    const double value = number(values, column, input_path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max())
        || std::floor(value) != value)
        throw std::invalid_argument("Invalid H20.39 CSV index: " + input_path);
    return static_cast<std::size_t>(value);
}

std::vector<std::vector<ReferenceNode>> read_reference(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input)
        throw std::runtime_error("Could not read H20.39 Abaqus reference: " + input_path);
    std::string line;
    if (!std::getline(input, line)
        || line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,gap,pressure")
        throw std::invalid_argument("Unexpected H20.39 CSV header: " + input_path);
    std::vector<std::vector<ReferenceNode>> result(path.size());
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split(line);
        const std::size_t step = index_value(values, 0, input_path);
        if (step == 0 || step > path.size())
            throw std::invalid_argument("Invalid H20.39 CSV step: " + input_path);
        result[step - 1].push_back({step - 1,
            index_value(values, 1, input_path),
            {number(values, 2, input_path), number(values, 3, input_path), number(values, 4, input_path)},
            {number(values, 5, input_path), number(values, 6, input_path), number(values, 7, input_path)},
            {number(values, 8, input_path), number(values, 9, input_path), number(values, 10, input_path)},
            number(values, 11, input_path),
            number(values, 12, input_path)});
    }
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
    contact.name = "finite_sliding_interface";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.mechanical = true;
    contact.penalty = 1.0e5;
    contact.friction_coefficient = 0.5;
    contact.friction_slip_tolerance = 1.0e-4;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    result.contacts.push_back(contact);
    return result;
}

bool compare_case(const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::string& reference_path,
    fuelsim::StrainFormulation strain,
    const std::string& prefix) {
    const std::vector<std::vector<ReferenceNode>> reference = read_reference(reference_path);
    fuelsim::SteadyProblem problem(definition(strain), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::vector<std::size_t> source_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    fuelsim::test::FieldErrorMetrics normal, tangential, resultant, gap, pressure, force_center;
    double maximum_coordinate_difference = 0.0, maximum_balance = 0.0;
    std::vector<double> state = problem.initial_state();
    bool all_projected_and_sliding = true, all_basis_initialized = true;
    for (std::size_t step = 0; step < path.size(); ++step) {
        state = problem.initial_state();
        const fuelsim::Hex20RegionMesh& secondary = spatial.hex20_region_mesh(1);
        for (std::size_t local = 0; local < secondary.nodes().size(); ++local) {
            const std::size_t global = spatial.global_node(1, local);
            state[spatial.dof(fuelsim::Field::displacement_x, global)] = path[step].displacement[0];
            state[spatial.dof(fuelsim::Field::displacement_y, global)] = path[step].displacement[1];
            state[spatial.dof(fuelsim::Field::displacement_z, global)] = path[step].displacement[2];
        }
        problem.validate_state(state);
        const std::vector<fuelsim::CartesianContactNodeSummary> actual =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        if (actual.size() != source_nodes.size() || actual.size() != reference[step].size())
            throw std::invalid_argument("H20.39 contact-node counts differ: " + reference_path);
        std::array<double, 3> actual_resultant{}, reference_resultant{};
        double actual_weight = 0.0, reference_weight = 0.0, actual_center_y = 0.0, actual_center_z = 0.0,
               reference_center_y = 0.0, reference_center_z = 0.0;
        for (std::size_t node = 0; node < actual.size(); ++node) {
            const auto found = std::find_if(reference[step].begin(),
                reference[step].end(),
                [&](const ReferenceNode& value) { return value.source_node == source_nodes[node]; });
            if (found == reference[step].end())
                throw std::invalid_argument("H20.39 source-node mapping is incomplete: " + reference_path);
            const fuelsim::CartesianPoint3& point = mesh.nodes().at(source_nodes[node]);
            maximum_coordinate_difference = std::max({maximum_coordinate_difference,
                std::abs(point.x - found->point.x),
                std::abs(point.y - found->point.y),
                std::abs(point.z - found->point.z)});
            all_projected_and_sliding = all_projected_and_sliding && actual[node].projected && actual[node].sliding;
            for (std::size_t component = 0; component < 3; ++component) {
                const double actual_normal = -actual[node].normal_contact_force[component],
                             actual_tangent = -actual[node].tangential_contact_force[component];
                normal.add(actual_normal, found->normal_force[component]);
                if (component != 0)
                    tangential.add(actual_tangent, found->tangential_force[component]);
                actual_resultant[component] += actual_normal;
                reference_resultant[component] += found->normal_force[component];
            }
            gap.add(actual[node].gap, found->gap);
            pressure.add(actual[node].pressure, found->pressure);
            const double actual_force = -actual[node].normal_contact_force[0], reference_force = found->normal_force[0];
            actual_weight += actual_force;
            reference_weight += reference_force;
            actual_center_y += actual_force * (point.y + path[step].displacement[1]);
            actual_center_z += actual_force * (point.z + path[step].displacement[2]);
            reference_center_y += reference_force * (found->point.y + path[step].displacement[1]);
            reference_center_z += reference_force * (found->point.z + path[step].displacement[2]);
        }
        for (std::size_t component = 0; component < 3; ++component) {
            double actual_tangential_resultant = 0.0, reference_tangential_resultant = 0.0;
            for (std::size_t node = 0; node < actual.size(); ++node) {
                const auto found = std::find_if(reference[step].begin(),
                    reference[step].end(),
                    [&](const ReferenceNode& value) { return value.source_node == source_nodes[node]; });
                actual_tangential_resultant -= actual[node].tangential_contact_force[component];
                reference_tangential_resultant += found->tangential_force[component];
            }
            if (std::abs(reference_tangential_resultant) < 1.0e-10)
                reference_tangential_resultant = 0.0;
            resultant.add(actual_tangential_resultant, reference_tangential_resultant);
            maximum_balance =
                std::max(maximum_balance, std::abs(actual_resultant[component] - reference_resultant[component]));
        }
        force_center.add(actual_center_y / actual_weight, reference_center_y / reference_weight);
        force_center.add(actual_center_z / actual_weight, reference_center_z / reference_weight);
        problem.commit_internal_state(state);
        const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        all_basis_initialized = all_basis_initialized && histories.size() == actual.size()
                                && std::all_of(histories.begin(), histories.end(), [](const auto& history) {
                                       return history.sliding && history.cartesian_tangent_basis_initialized;
                                   });
    }
    fuelsim::test::print_relative_metrics(prefix + "normal_nodal_force", normal);
    fuelsim::test::print_relative_metrics(prefix + "tangential_nodal_force", tangential);
    fuelsim::test::print_relative_metrics(prefix + "tangential_resultant", resultant);
    fuelsim::test::print_relative_metrics(prefix + "gap", gap);
    fuelsim::test::print_relative_metrics(prefix + "pressure", pressure);
    fuelsim::test::print_relative_metrics(prefix + "force_center", force_center);
    std::cout << prefix << "maximum_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << prefix << "maximum_resultant_absolute_difference=" << maximum_balance << '\n';
    constexpr double tolerance = 1.0e-2, zero_tolerance = 1.0e-8;
    return check(maximum_coordinate_difference < 3.0e-8, prefix + "uses the exact tracked Exodus coordinates in Abaqus")
           && check(all_projected_and_sliding && all_basis_initialized,
               prefix + "keeps all contact nodes projected and all node-centered constraint histories sliding")
           && check(fuelsim::test::relative_metrics_below(normal, tolerance)
                        && normal.maximum_zero_reference_difference < zero_tolerance,
               prefix + "normal nodal-force metrics and zero-reference values agree with Abaqus below 1 percent")
           && check(fuelsim::test::relative_metrics_below(tangential, tolerance)
                        && tangential.maximum_zero_reference_difference < zero_tolerance,
               prefix + "two tangent-plane nodal-force components agree with Abaqus below 1 percent")
           && check(fuelsim::test::relative_metrics_below(resultant, tolerance)
                        && resultant.maximum_zero_reference_difference < zero_tolerance,
               prefix + "complete three-dimensional tangential resultants agree with Abaqus below 1 percent")
           && check(fuelsim::test::relative_metrics_below(gap, tolerance)
                        && fuelsim::test::relative_metrics_below(pressure, tolerance),
               prefix + "gap and pressure metrics agree with Abaqus below 1 percent")
           && check(fuelsim::test::relative_metrics_below(force_center, tolerance),
               prefix + "normal-force center follows the finite-sliding path below 1 percent");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--generate") {
            generate(argv[2]);
            return 0;
        }
        if (argc != 4) {
            std::cerr << "Usage: fuelsim_h20_39_finite_sliding_abaqus_tests <mesh.e> <small.csv> <finite.csv>\n"
                         "   or: fuelsim_h20_39_finite_sliding_abaqus_tests --generate <abaqus-directory>\n";
            return 2;
        }
        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(argv[1]);
        const bool small = compare_case(mesh, argv[2], fuelsim::StrainFormulation::small, "h20_39_small_");
        const bool finite = compare_case(mesh, argv[3], fuelsim::StrainFormulation::finite, "h20_39_finite_");
        if (small && finite)
            std::cout << "[PASS] H20.39 HEX20 finite-sliding Abaqus comparison\n";
        return small && finite ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.39 raised: " << error.what() << '\n';
        return 1;
    }
}
