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
constexpr double pi = 3.141592653589793238462643383279502884;

struct PathStep final {
    const char* name;
    double rotation, axial_shift;
};

constexpr std::array<PathStep, 4> path = {
    {{"CLOSE", 0.03, 0.04}, {"CROSS", pi / 4.0, 0.3}, {"STRADDLE", pi / 32.0, -0.2}, {"RETURN", 0.0, 0.1}}};

struct ReferenceNode final {
    std::size_t source_node;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> normal_force{}, tangential_force{};
    double slip_first, slip_second, gap, pressure;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::CartesianPoint3 cylindrical(double radius, double angle, double z) {
    return {radius * std::cos(angle), radius * std::sin(angle), z};
}

using CoordinateKey = std::array<std::int64_t, 3>;

CoordinateKey coordinate_key(const fuelsim::CartesianPoint3& point) {
    constexpr double scale = 1.0e12;
    return {static_cast<std::int64_t>(std::llround(scale * point.x)),
        static_cast<std::int64_t>(std::llround(scale * point.y)),
        static_cast<std::int64_t>(std::llround(scale * point.z))};
}

fuelsim::Hex20Element append_annular_element(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<CoordinateKey, std::size_t>& node_map, double radius_lower, double radius_upper, double angle_lower,
    double angle_upper, double z_lower, double z_upper) {
    const std::array<std::array<double, 3>, 8> corners = {
        {{radius_lower, angle_lower, z_lower}, {radius_upper, angle_lower, z_lower},
            {radius_upper, angle_upper, z_lower}, {radius_lower, angle_upper, z_lower},
            {radius_lower, angle_lower, z_upper}, {radius_upper, angle_lower, z_upper},
            {radius_upper, angle_upper, z_upper}, {radius_lower, angle_upper, z_upper}}};
    const std::array<std::pair<std::size_t, std::size_t>, 12> edges = {
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}}};
    std::array<std::array<double, 3>, 20> points{};
    std::copy(corners.begin(), corners.end(), points.begin());
    for (std::size_t edge = 0; edge < edges.size(); ++edge)
        for (std::size_t component = 0; component < 3; ++component)
            points[8 + edge][component] =
                0.5 * (corners[edges[edge].first][component] + corners[edges[edge].second][component]);
    fuelsim::Hex20Element element{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const fuelsim::CartesianPoint3 point = cylindrical(points[local][0], points[local][1], points[local][2]);
        const auto inserted = node_map.emplace(coordinate_key(point), nodes.size());
        if (inserted.second) nodes.push_back(point);
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
        for (const std::size_t local : face_nodes.at(side.local_side))
            result.push_back(elements.at(side.element).nodes[local]);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

fuelsim::UnstructuredHex20Mesh generate_mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<CoordinateKey, std::size_t> primary_nodes, secondary_nodes;
    std::vector<fuelsim::Hex20Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, secondary_contact;
    for (std::size_t angle = 0; angle < 16; ++angle) {
        const std::size_t element = elements.size();
        const double lower = -pi / 4.0 + static_cast<double>(angle) * pi / 16.0;
        elements.push_back(append_annular_element(nodes, primary_nodes, 0.9, 1.0, lower, lower + pi / 16.0, -0.5, 1.5));
        blocks.push_back(1);
        primary_contact.push_back({element, 1});
    }
    for (std::size_t angle = 0; angle < 8; ++angle) {
        const std::size_t element = elements.size();
        const double lower = static_cast<double>(angle) * pi / 16.0;
        elements.push_back(
            append_annular_element(nodes, secondary_nodes, 1.0, 1.1, lower, lower + pi / 16.0, 0.0, 1.0));
        blocks.push_back(2);
        secondary_contact.push_back({element, 3});
    }
    std::vector<std::size_t> primary_all, secondary_all;
    for (const auto& entry : primary_nodes) primary_all.push_back(entry.second);
    for (const auto& entry : secondary_nodes) secondary_all.push_back(entry.second);
    std::sort(primary_all.begin(), primary_all.end());
    std::sort(secondary_all.begin(), secondary_all.end());
    const std::vector<std::size_t> secondary_contact_nodes = side_nodes(elements, secondary_contact);
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all", std::move(primary_all)}, {20, "secondary_all", std::move(secondary_all)},
            {21, "secondary_contact_nodes", secondary_contact_nodes}},
        {{30, "primary_contact", std::move(primary_contact)}, {40, "secondary_contact", std::move(secondary_contact)}});
}

std::array<double, 3> secondary_displacement(const fuelsim::CartesianPoint3& point, const PathStep& step) {
    const double radius = std::hypot(point.x, point.y), angle = std::atan2(point.y, point.x) + step.rotation;
    return {
        (radius - 0.001) * std::cos(angle) - point.x, (radius - 0.001) * std::sin(angle) - point.y, step.axial_shift};
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
    if (!output) throw std::runtime_error("Could not write H20.40 Abaqus input: " + output_path);
    output << std::setprecision(16) << "*Heading\n"
           << "** H20.40 true-quadratic curved HEX20 finite-sliding friction path.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    constexpr std::array<std::size_t, 20> abaqus_order = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
    const auto write_block = [&](const std::string& name, std::int64_t block) {
        output << "*Element, type=C3D20, elset=" << name << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != block) continue;
            output << element + 1;
            for (const std::size_t local : abaqus_order) output << ", " << mesh.elements()[element].nodes[local] + 1;
            output << '\n';
        }
    };
    write_block("PRIMARY", 1);
    write_block("SECONDARY", 2);
    for (const fuelsim::NodeSet& set : mesh.node_sets()) {
        output << "*Nset, nset=" << set.name << '\n';
        write_labels(output, set.nodes);
    }
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.25\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e6,\n"
           << "*Friction, slip tolerance=1.e-4\n0.5,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n";
    const auto secondary_all_entry = std::find_if(
        mesh.node_sets().begin(), mesh.node_sets().end(), [](const auto& set) { return set.name == "secondary_all"; });
    if (secondary_all_entry == mesh.node_sets().end()) throw std::logic_error("H20.40 secondary_all set is missing");
    for (std::size_t step = 0; step < path.size(); ++step) {
        output << "*Step, name=" << path[step].name << ", nlgeom=" << (nlgeom ? "YES" : "NO") << ", inc=100\n"
               << "*Static\n0.025, 1., 1.e-8, 0.025\n"
               << (step == 0 ? "*Boundary\nprimary_all, 1, 3, 0.\n" : "*Boundary, op=MOD\n");
        for (const std::size_t node : secondary_all_entry->nodes) {
            const std::array<double, 3> displacement = secondary_displacement(mesh.nodes()[node], path[step]);
            for (std::size_t component = 0; component < 3; ++component)
                output << node + 1 << ", " << component + 1 << ", " << component + 1 << ", " << displacement[component]
                       << '\n';
        }
        output << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
    }
}

void generate(const std::string& directory) {
    const fuelsim::UnstructuredHex20Mesh mesh = generate_mesh();
    fuelsim::write_exodus_hex20(directory + "/h20_40_hex20_curved_finite_sliding_mesh.e", mesh);
    write_abaqus_input(directory + "/h20_40_hex20_curved_finite_sliding_small.inp", mesh, false);
    write_abaqus_input(directory + "/h20_40_hex20_curved_finite_sliding_finite.inp", mesh, true);
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete H20.40 CSV row: " + input_path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid H20.40 CSV number: " + input_path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    const double value = number(values, column, input_path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("Invalid H20.40 CSV index: " + input_path);
    return static_cast<std::size_t>(value);
}

std::vector<std::vector<ReferenceNode>> read_reference(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Could not read H20.40 Abaqus reference: " + input_path);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,"
                                              "tangential_z,slip_1,slip_2,gap,pressure")
        throw std::invalid_argument("Unexpected H20.40 CSV header: " + input_path);
    std::vector<std::vector<ReferenceNode>> result(path.size());
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        const std::size_t step = index_value(values, 0, input_path);
        if (step == 0 || step > path.size()) throw std::invalid_argument("Invalid H20.40 CSV step: " + input_path);
        result[step - 1].push_back({index_value(values, 1, input_path),
            {number(values, 2, input_path), number(values, 3, input_path), number(values, 4, input_path)},
            {number(values, 5, input_path), number(values, 6, input_path), number(values, 7, input_path)},
            {number(values, 8, input_path), number(values, 9, input_path), number(values, 10, input_path)},
            number(values, 11, input_path), number(values, 12, input_path), number(values, 13, input_path),
            number(values, 14, input_path)});
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
    contact.name = "curved_finite_sliding_interface";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.mechanical = true;
    contact.penalty = 1.0e6;
    contact.friction_coefficient = 0.5;
    contact.friction_slip_tolerance = 1.0e-4;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    result.contacts.push_back(contact);
    return result;
}

bool compare_case(const fuelsim::UnstructuredHex20Mesh& mesh, const std::string& reference_path,
    fuelsim::StrainFormulation strain, const std::string& prefix) {
    const std::vector<std::vector<ReferenceNode>> reference = read_reference(reference_path);
    fuelsim::SteadyProblem problem(definition(strain), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::vector<std::size_t> source_nodes =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    fuelsim::test::FieldErrorMetrics normal_radial, tangential_circumferential, tangential_axial, resultant, gap,
        pressure;
    double maximum_coordinate_difference = 0.0;
    bool all_projected_and_sliding = true, all_basis_initialized = true;
    for (std::size_t step = 0; step < path.size(); ++step) {
        const fuelsim::Hex20RegionMesh& secondary = spatial.hex20_region_mesh(1);
        std::vector<double> state;
        constexpr std::size_t increments_per_step = 40;
        for (std::size_t increment = 1; increment <= increments_per_step; ++increment) {
            state = problem.initial_state();
            const double fraction = static_cast<double>(increment) / static_cast<double>(increments_per_step);
            for (std::size_t local = 0; local < secondary.nodes().size(); ++local) {
                const std::size_t global = spatial.global_node(1, local);
                const std::array<double, 3> current = secondary_displacement(secondary.nodes()[local], path[step]);
                const std::array<double, 3> previous =
                    step == 0 ? std::array<double, 3>{}
                              : secondary_displacement(secondary.nodes()[local], path[step - 1]);
                for (std::size_t component = 0; component < 3; ++component) {
                    const double displacement =
                        previous[component] + fraction * (current[component] - previous[component]);
                    const fuelsim::Field field = component == 0   ? fuelsim::Field::displacement_x
                                                 : component == 1 ? fuelsim::Field::displacement_y
                                                                  : fuelsim::Field::displacement_z;
                    state[spatial.dof(field, global)] = displacement;
                }
            }
            problem.validate_state(state);
            if (increment != increments_per_step) problem.commit_internal_state(state);
        }
        const std::vector<fuelsim::CartesianContactNodeSummary> actual =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        if (actual.size() != source_nodes.size() || actual.size() != reference[step].size())
            throw std::invalid_argument("H20.40 contact-node counts differ: " + reference_path);
        std::array<double, 3> actual_resultant{}, reference_resultant{};
        for (std::size_t node = 0; node < actual.size(); ++node) {
            const auto found = std::find_if(reference[step].begin(), reference[step].end(),
                [&](const ReferenceNode& value) { return value.source_node == source_nodes[node]; });
            if (found == reference[step].end())
                throw std::invalid_argument("H20.40 source-node mapping is incomplete: " + reference_path);
            const fuelsim::CartesianPoint3& point = mesh.nodes().at(source_nodes[node]);
            maximum_coordinate_difference = std::max({maximum_coordinate_difference, std::abs(point.x - found->point.x),
                std::abs(point.y - found->point.y), std::abs(point.z - found->point.z)});
            const double current_angle = std::atan2(point.y, point.x) + path[step].rotation,
                         cosine = std::cos(current_angle), sine = std::sin(current_angle);
            std::array<double, 3> actual_normal{}, actual_tangent{};
            for (std::size_t component = 0; component < 3; ++component) {
                actual_normal[component] = -actual[node].normal_contact_force[component];
                actual_tangent[component] = -actual[node].tangential_contact_force[component];
                actual_resultant[component] += actual_normal[component] + actual_tangent[component];
                reference_resultant[component] += found->normal_force[component] + found->tangential_force[component];
            }
            normal_radial.add(cosine * actual_normal[0] + sine * actual_normal[1],
                cosine * found->normal_force[0] + sine * found->normal_force[1]);
            tangential_circumferential.add(-sine * actual_tangent[0] + cosine * actual_tangent[1],
                -sine * found->tangential_force[0] + cosine * found->tangential_force[1]);
            tangential_axial.add(actual_tangent[2], found->tangential_force[2]);
            gap.add(actual[node].gap, found->gap);
            pressure.add(actual[node].pressure, found->pressure);
            all_projected_and_sliding = all_projected_and_sliding && actual[node].projected && actual[node].sliding;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            if (std::abs(reference_resultant[component]) < 1.0e-10) reference_resultant[component] = 0.0;
            resultant.add(actual_resultant[component], reference_resultant[component]);
        }
        problem.commit_internal_state(state);
        const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        all_basis_initialized = all_basis_initialized && histories.size() == 72 &&
                                std::all_of(histories.begin(), histories.end(), [](const auto& history) {
                                    return history.sliding && history.cartesian_tangent_basis_initialized;
                                });
    }
    fuelsim::test::print_relative_metrics(prefix + "normal_radial_nodal_force", normal_radial);
    fuelsim::test::print_relative_metrics(
        prefix + "tangential_circumferential_nodal_force", tangential_circumferential);
    fuelsim::test::print_relative_metrics(prefix + "tangential_axial_nodal_force", tangential_axial);
    fuelsim::test::print_relative_metrics(prefix + "complete_contact_resultant", resultant);
    fuelsim::test::print_relative_metrics(prefix + "gap", gap);
    fuelsim::test::print_relative_metrics(prefix + "pressure", pressure);
    std::cout << prefix << "maximum_coordinate_difference=" << maximum_coordinate_difference << '\n';
    constexpr double tolerance = 1.0e-2, zero_tolerance = 1.0e-8;
    return check(maximum_coordinate_difference < 3.1e-8,
               prefix + "uses the tracked Exodus coordinates within Abaqus output-database precision") &&
           check(all_projected_and_sliding && all_basis_initialized,
               prefix + "keeps every curved contact node projected and all 72 histories sliding") &&
           check(fuelsim::test::relative_metrics_below(normal_radial, tolerance) &&
                     normal_radial.maximum_zero_reference_difference < zero_tolerance,
               prefix + "radial normal nodal-force metrics agree with Abaqus below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(tangential_circumferential, tolerance) &&
                     tangential_circumferential.maximum_zero_reference_difference < zero_tolerance &&
                     fuelsim::test::relative_metrics_below(tangential_axial, tolerance) &&
                     tangential_axial.maximum_zero_reference_difference < zero_tolerance,
               prefix + "two curved tangent-plane nodal-force metrics agree with Abaqus below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(resultant, tolerance) &&
                     resultant.maximum_zero_reference_difference < zero_tolerance,
               prefix + "complete contact-resultant metrics agree with Abaqus below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(gap, tolerance) &&
                     fuelsim::test::relative_metrics_below(pressure, tolerance),
               prefix + "curved gap and pressure metrics agree with Abaqus below 1 percent");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--generate") {
            generate(argv[2]);
            return 0;
        }
        if (argc != 4) {
            std::cerr << "Usage: fuelsim_h20_40_curved_finite_sliding_abaqus_tests <mesh.e> <small.csv> <finite.csv>\n"
                         "   or: fuelsim_h20_40_curved_finite_sliding_abaqus_tests --generate <abaqus-directory>\n";
            return 2;
        }
        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(argv[1]);
        const bool small = compare_case(mesh, argv[2], fuelsim::StrainFormulation::small, "h20_40_small_");
        const bool finite = compare_case(mesh, argv[3], fuelsim::StrainFormulation::finite, "h20_40_finite_");
        if (small && finite) std::cout << "[PASS] H20.40 curved HEX20 finite-sliding Abaqus comparison\n";
        return small && finite ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.40 raised: " << error.what() << '\n';
        return 1;
    }
}
