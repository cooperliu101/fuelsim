#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/contact.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::array<std::array<std::size_t, 4>, 6> face_nodes = {
    {{{0, 1, 5, 4}}, {{1, 2, 6, 5}}, {{2, 3, 7, 6}}, {{0, 4, 7, 3}}, {{0, 3, 2, 1}}, {{4, 5, 6, 7}}}};

struct GeneratedCase final {
    std::string name;
    fuelsim::UnstructuredHex8Mesh mesh;
    std::vector<std::array<double, 3>> displacement;
};

struct NodeReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> displacement{};
};

struct ContactReference final {
    bool secondary;
    std::size_t id;
    fuelsim::CartesianPoint3 current;
    std::array<double, 3> force{};
    double gap, pressure;
};

struct ReconstructedContact final {
    std::map<std::size_t, std::array<double, 3>> secondary_force, primary_force;
    std::map<std::size_t, double> gap, pressure, area;
    std::array<double, 3> balance{};
    double jacobian_error = 0.0;
    std::size_t points = 0;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
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

std::vector<std::size_t> side_nodes(
    const std::vector<fuelsim::Hex8Element>& elements, const std::vector<fuelsim::ElementSide>& sides) {
    std::vector<std::size_t> result;
    for (const auto& side : sides)
        for (const std::size_t local : face_nodes[side.local_side])
            result.push_back(elements[side.element].nodes[local]);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

GeneratedCase faceted_case() {
    const std::array<double, 4> angles = {0.0, pi / 6.0, pi / 3.0, pi / 2.0};
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_map, secondary_map;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, secondary_contact;
    for (std::size_t angle = 0; angle + 1 < angles.size(); ++angle) {
        elements.push_back(append_annular(nodes, primary_map, 0.8, 1.0, angles[angle], angles[angle + 1]));
        blocks.push_back(1);
        primary_contact.push_back({elements.size() - 1, 1});
    }
    for (std::size_t angle = 0; angle + 1 < angles.size(); ++angle) {
        elements.push_back(append_annular(nodes, secondary_map, 1.005, 1.2, angles[angle], angles[angle + 1]));
        blocks.push_back(2);
        secondary_contact.push_back({elements.size() - 1, 3});
    }
    std::vector<std::size_t> primary_all, secondary_all;
    for (const auto& entry : primary_map) primary_all.push_back(entry.second);
    for (const auto& entry : secondary_map) secondary_all.push_back(entry.second);
    std::sort(primary_all.begin(), primary_all.end());
    std::sort(secondary_all.begin(), secondary_all.end());
    const std::vector<std::size_t> primary_contact_nodes = side_nodes(elements, primary_contact),
                                   secondary_contact_nodes = side_nodes(elements, secondary_contact);
    std::vector<std::array<double, 3>> displacement(nodes.size());
    for (const std::size_t node : secondary_all) {
        const double radius = std::hypot(nodes[node].x, nodes[node].y);
        displacement[node] = {-0.015 * nodes[node].x / radius, -0.015 * nodes[node].y / radius, 0.0};
    }
    return {"faceted_cylinder",
        fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(blocks),
            {{1, "primary"}, {2, "secondary"}},
            {{10, "primary_all", primary_all}, {20, "secondary_all", secondary_all},
                {30, "primary_contact_nodes", primary_contact_nodes},
                {40, "secondary_contact_nodes", secondary_contact_nodes}},
            {{50, "primary_contact", primary_contact}, {60, "secondary_contact", secondary_contact}}),
        std::move(displacement)};
}

GeneratedCase warped_case() {
    constexpr double warp = 0.125;
    constexpr double initial_gap = 1.0 / 256.0;
    constexpr double imposed_x = -3.0 / 256.0;
    std::vector<fuelsim::CartesianPoint3> nodes = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0 + warp, 1.0, 1.0}, {0.0, 1.0, 1.0}, {1.0 + initial_gap, 0.0, 0.0},
        {2.0 + initial_gap, 0.0, 0.0}, {2.0 + initial_gap, 1.0, 0.0}, {1.0 + initial_gap, 1.0, 0.0},
        {1.0 + initial_gap, 0.0, 1.0}, {2.0 + initial_gap, 0.0, 1.0}, {2.0 + warp + initial_gap, 1.0, 1.0},
        {1.0 + warp + initial_gap, 1.0, 1.0}};
    std::vector<fuelsim::Hex8Element> elements(2);
    elements[0].nodes = {0, 1, 2, 3, 4, 5, 6, 7};
    elements[1].nodes = {8, 9, 10, 11, 12, 13, 14, 15};
    std::vector<std::array<double, 3>> displacement(nodes.size());
    for (std::size_t node = 8; node < nodes.size(); ++node) displacement[node] = {imposed_x, 0.0, 0.0};
    return {"warped_bilinear",
        fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), {1, 2}, {{1, "primary"}, {2, "secondary"}},
            {{10, "primary_all", {0, 1, 2, 3, 4, 5, 6, 7}}, {20, "secondary_all", {8, 9, 10, 11, 12, 13, 14, 15}},
                {30, "primary_contact_nodes", {1, 2, 5, 6}}, {40, "secondary_contact_nodes", {8, 11, 12, 15}}},
            {{50, "primary_contact", {{0, 1}}}, {60, "secondary_contact", {{1, 3}}}}),
        std::move(displacement)};
}

std::vector<GeneratedCase> cases() {
    std::vector<GeneratedCase> result;
    result.push_back(faceted_case());
    result.push_back(warped_case());
    return result;
}

void write_labels(std::ofstream& output, const std::vector<std::size_t>& nodes) {
    for (std::size_t index = 0; index < nodes.size(); ++index)
        output << nodes[index] + 1 << (index + 1 == nodes.size() ? "\n" : ", ");
}

void write_input(const std::string& path, const GeneratedCase& generated, const std::string& smoothing) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B4.5 Abaqus input: " + path);
    const auto& mesh = generated.mesh;
    output << std::setprecision(16) << "*Heading\n** B4.5 " << generated.name
           << " C3D8 nonplanar finite-sliding contact probe.\n*Preprint, echo=NO, model=NO, history=NO, "
              "contact=YES\n*Node\n";
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
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.0\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n*Surface Behavior, pressure-overclosure=LINEAR\n1.e7,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.";
    if (!smoothing.empty()) output << ", sliding transition=" << smoothing << " smoothing";
    output << "\nSECONDARY_CONTACT, PRIMARY_CONTACT\n*Step, name=LOAD, nlgeom=YES, inc=1\n"
           << "*Static\n1., 1., 1., 1.\n*Boundary, op=NEW\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node)
        for (std::size_t component = 0; component < 3; ++component)
            output << node + 1 << ", " << component + 1 << ", " << component + 1 << ", "
                   << generated.displacement[node][component] << '\n';
    output << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
           << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
}

void generate(const std::string& directory) {
    for (const auto& generated : cases()) {
        const std::string base = directory + "/b45_hex8_" + generated.name;
        fuelsim::write_exodus_hex8(base + "_mesh.e", generated.mesh);
        write_input(base + ".inp", generated, "");
        if (generated.name == "warped_bilinear") {
            write_input(base + "_linear.inp", generated, "LINEAR");
            write_input(base + "_quadratic.inp", generated, "QUADRATIC");
        }
    }
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column) { return std::stod(values.at(column)); }

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.5 node reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "id,x,y,z,ux,uy,uz") throw std::invalid_argument("Invalid B4.5 node CSV");
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const auto values = split(line);
        result.push_back(
            {static_cast<std::size_t>(number(values, 0)), {number(values, 1), number(values, 2), number(values, 3)},
                {number(values, 4), number(values, 5), number(values, 6)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.5 contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "side,id,x,y,z,force_x,force_y,force_z,gap,pressure")
        throw std::invalid_argument("Invalid B4.5 contact CSV");
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        const auto values = split(line);
        result.push_back({values[0] == "secondary", static_cast<std::size_t>(number(values, 1)),
            {number(values, 2), number(values, 3), number(values, 4)},
            {number(values, 5), number(values, 6), number(values, 7)}, number(values, 8), number(values, 9)});
    }
    return result;
}

fuelsim::Quad4FaceCoordinates face_coordinates(
    const fuelsim::UnstructuredHex8Mesh& mesh, const fuelsim::ElementSide& side) {
    fuelsim::Quad4FaceCoordinates result{};
    for (std::size_t node = 0; node < 4; ++node)
        result[node] = mesh.nodes()[mesh.elements()[side.element].nodes[face_nodes[side.local_side][node]]];
    return result;
}

ReconstructedContact reconstruct(const GeneratedCase& generated) {
    const auto& mesh = generated.mesh;
    const auto& primary_sides = mesh.side_set("primary_contact").sides;
    const auto& secondary_sides = mesh.side_set("secondary_contact").sides;
    if (primary_sides.size() != secondary_sides.size()) throw std::logic_error("B4.5 paired face count differs");
    ReconstructedContact result;
    for (std::size_t face = 0; face < secondary_sides.size(); ++face) {
        const auto secondary = face_coordinates(mesh, secondary_sides[face]);
        const auto primary = face_coordinates(mesh, primary_sides[face]);
        fuelsim::Quad4SurfaceContactLocalValues state{}, committed{};
        for (std::size_t node = 0; node < 4; ++node) {
            const std::size_t secondary_source = mesh.elements()[secondary_sides[face].element]
                                                     .nodes[face_nodes[secondary_sides[face].local_side][node]];
            for (std::size_t component = 0; component < 3; ++component)
                state[8 * (component + 1) + node] = generated.displacement[secondary_source][component];
        }
        for (std::size_t point = 0; point < 4; ++point) {
            const double xi = point == 0 || point == 3 ? -0.5 : 0.5, eta = point < 2 ? -0.5 : 0.5;
            const auto quadrature = fuelsim::make_quad4_face_quadrature_point(secondary, xi, eta, 1.0);
            const fuelsim::Quad4ToQuad4MechanicalGeometry geometry{
                secondary, primary, quadrature.shape, quadrature.derivative_xi, quadrature.derivative_eta, 1.0, -1.0};
            fuelsim::Quad4SurfaceContactLocalJacobian jacobian{};
            const auto residual = fuelsim::compute_quad4_to_quad4_contact(
                {1.0e7, 0.0, false, 0.0}, geometry, state, committed, {}, &jacobian);
            const auto value =
                fuelsim::compute_quad4_to_quad4_contact_value({1.0e7, 0.0, false, 0.0}, geometry, state, committed, {});
            if (!value.projected || !(value.pressure > 0.0))
                throw std::logic_error("B4.5 reconstructed point is inactive");
            fuelsim::Quad4SurfaceContactLocalValues direction{}, plus = state, minus = state;
            constexpr double perturbation = 1.0e-8;
            for (std::size_t column = 0; column < direction.size(); ++column) {
                direction[column] = std::sin(static_cast<double>(column + 1));
                plus[column] += perturbation * direction[column];
                minus[column] -= perturbation * direction[column];
            }
            const auto plus_residual = fuelsim::compute_quad4_to_quad4_contact(
                {1.0e7, 0.0, false, 0.0}, geometry, plus, committed, {}, nullptr);
            const auto minus_residual = fuelsim::compute_quad4_to_quad4_contact(
                {1.0e7, 0.0, false, 0.0}, geometry, minus, committed, {}, nullptr);
            double difference_squared = 0.0, reference_squared = 0.0;
            for (std::size_t row = 0; row < residual.size(); ++row) {
                double analytic = 0.0;
                for (std::size_t column = 0; column < direction.size(); ++column)
                    analytic += jacobian[row * direction.size() + column] * direction[column];
                const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
                difference_squared += (analytic - reference) * (analytic - reference);
                reference_squared += reference * reference;
            }
            if (reference_squared > 0.0)
                result.jacobian_error =
                    std::max(result.jacobian_error, std::sqrt(difference_squared / reference_squared));
            const std::size_t secondary_source = mesh.elements()[secondary_sides[face].element]
                                                     .nodes[face_nodes[secondary_sides[face].local_side][point]];
            result.gap[secondary_source] = value.gap;
            result.pressure[secondary_source] = value.pressure;
            result.area[secondary_source] += value.tributary_area;
            for (std::size_t component = 0; component < 3; ++component) {
                for (std::size_t local = 0; local < 4; ++local) {
                    const std::size_t secondary_node = mesh.elements()[secondary_sides[face].element]
                                                           .nodes[face_nodes[secondary_sides[face].local_side][local]],
                                      primary_node = mesh.elements()[primary_sides[face].element]
                                                         .nodes[face_nodes[primary_sides[face].local_side][local]];
                    result.secondary_force[secondary_node][component] -= residual[8 * (component + 1) + local];
                    result.primary_force[primary_node][component] += residual[8 * (component + 1) + 4 + local];
                    result.balance[component] +=
                        residual[8 * (component + 1) + local] + residual[8 * (component + 1) + 4 + local];
                }
            }
            ++result.points;
        }
    }
    return result;
}

double norm(const std::array<double, 3>& value) { return std::hypot(value[0], std::hypot(value[1], value[2])); }

double primary_face_warp(const GeneratedCase& generated) {
    const auto side = generated.mesh.side_set("primary_contact").sides.front();
    const auto face = face_coordinates(generated.mesh, side);
    const std::array<double, 3> first = {face[1].x - face[0].x, face[1].y - face[0].y, face[1].z - face[0].z},
                                second = {face[3].x - face[0].x, face[3].y - face[0].y, face[3].z - face[0].z},
                                diagonal = {face[2].x - face[0].x, face[2].y - face[0].y, face[2].z - face[0].z};
    const std::array<double, 3> normal = {first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2], first[0] * second[1] - first[1] * second[0]};
    return std::abs(normal[0] * diagonal[0] + normal[1] * diagonal[1] + normal[2] * diagonal[2]) / norm(normal);
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

bool compare(const GeneratedCase& generated, const std::string& mesh_path, const std::string& node_path,
    const std::string& contact_path, const std::string& variant, bool require_transverse_force) {
    fuelsim::UnstructuredHex8Mesh tracked_mesh = fuelsim::read_exodus_hex8(mesh_path);
    if (tracked_mesh.nodes().size() != generated.mesh.nodes().size() ||
        tracked_mesh.elements().size() != generated.mesh.elements().size())
        throw std::invalid_argument("B4.5 tracked mesh size differs from its generator");
    double tracked_coordinate_error = 0.0;
    for (std::size_t node = 0; node < tracked_mesh.nodes().size(); ++node)
        tracked_coordinate_error =
            std::max({tracked_coordinate_error, std::abs(tracked_mesh.nodes()[node].x - generated.mesh.nodes()[node].x),
                std::abs(tracked_mesh.nodes()[node].y - generated.mesh.nodes()[node].y),
                std::abs(tracked_mesh.nodes()[node].z - generated.mesh.nodes()[node].z)});
    for (std::size_t element = 0; element < tracked_mesh.elements().size(); ++element)
        if (tracked_mesh.elements()[element].nodes != generated.mesh.elements()[element].nodes)
            throw std::invalid_argument("B4.5 tracked mesh connectivity differs from its generator");
    GeneratedCase tracked{generated.name + variant, std::move(tracked_mesh), generated.displacement};
    const auto nodes = read_nodes(node_path);
    const auto contacts = read_contact(contact_path);
    const auto actual = reconstruct(tracked);
    const double face_warp = primary_face_warp(tracked);
    fuelsim::test::FieldErrorMetrics displacement[3], secondary_force[3], primary_force[3], gap, pressure, area;
    double coordinate_error = 0.0;
    for (const auto& reference : nodes) {
        const auto& point = tracked.mesh.nodes().at(reference.id);
        coordinate_error = std::max({coordinate_error, std::abs(point.x - reference.point.x),
            std::abs(point.y - reference.point.y), std::abs(point.z - reference.point.z)});
        for (std::size_t component = 0; component < 3; ++component)
            displacement[component].add(
                tracked.displacement[reference.id][component], reference.displacement[component]);
    }
    std::array<double, 3> actual_resultant{}, reference_resultant{};
    double maximum_normal_angle = 0.0;
    for (const auto& reference : contacts) {
        const auto& forces = reference.secondary ? actual.secondary_force : actual.primary_force;
        const auto found = forces.find(reference.id);
        if (found == forces.end()) throw std::logic_error("B4.5 contact node mapping is incomplete");
        for (std::size_t component = 0; component < 3; ++component) {
            (reference.secondary ? secondary_force[component] : primary_force[component])
                .add(found->second[component], reference.force[component]);
            if (reference.secondary) {
                actual_resultant[component] += found->second[component];
                reference_resultant[component] += reference.force[component];
            }
        }
        if (reference.secondary) {
            gap.add(actual.gap.at(reference.id), reference.gap);
            pressure.add(actual.pressure.at(reference.id), reference.pressure);
            if (generated.name == "warped_bilinear" && reference.pressure > 0.0) {
                area.add(actual.area.at(reference.id), norm(reference.force) / reference.pressure);
                const double actual_norm = norm(found->second), reference_norm = norm(reference.force);
                if (actual_norm > 0.0 && reference_norm > 0.0) {
                    double cosine = 0.0;
                    for (std::size_t component = 0; component < 3; ++component)
                        cosine +=
                            found->second[component] * reference.force[component] / (actual_norm * reference_norm);
                    cosine = std::max(-1.0, std::min(1.0, cosine));
                    maximum_normal_angle = std::max(maximum_normal_angle, std::acos(cosine) * 180.0 / pi);
                }
            }
        }
    }
    fuelsim::test::FieldErrorMetrics resultant;
    resultant.add(norm(actual_resultant), norm(reference_resultant));
    for (std::size_t component = 0; component < 3; ++component) {
        print_metric("b45_" + tracked.name + "_displacement_" + std::to_string(component), displacement[component]);
        print_metric(
            "b45_" + tracked.name + "_secondary_force_" + std::to_string(component), secondary_force[component]);
        print_metric("b45_" + tracked.name + "_primary_force_" + std::to_string(component), primary_force[component]);
    }
    print_metric("b45_" + tracked.name + "_gap", gap);
    print_metric("b45_" + tracked.name + "_pressure", pressure);
    print_metric("b45_" + tracked.name + "_area", area);
    print_metric("b45_" + tracked.name + "_resultant", resultant);
    const auto passes = [](const fuelsim::test::FieldErrorMetrics& metric) {
        return !metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, 1.0e-2);
    };
    bool metrics = passes(gap) && passes(pressure) && passes(area) && passes(resultant) && passes(secondary_force[0]) &&
                   passes(primary_force[0]);
    for (std::size_t component = 0; component < 3; ++component) metrics = metrics && passes(displacement[component]);
    if (require_transverse_force)
        for (std::size_t component = 1; component < 3; ++component)
            metrics = metrics && passes(secondary_force[component]) && passes(primary_force[component]);
    const bool qualified_transverse =
        generated.name != "warped_bilinear" || ((!passes(secondary_force[1]) || !passes(secondary_force[2])) &&
                                                   fuelsim::test::relative_metrics_below(secondary_force[1], 0.15) &&
                                                   fuelsim::test::relative_metrics_below(secondary_force[2], 0.15) &&
                                                   fuelsim::test::relative_metrics_below(primary_force[1], 0.15) &&
                                                   fuelsim::test::relative_metrics_below(primary_force[2], 0.15) &&
                                                   maximum_normal_angle > 0.4 && maximum_normal_angle < 0.5);
    std::cout << "b45_" << tracked.name << "_tracked_coordinate_error=" << tracked_coordinate_error << '\n'
              << "b45_" << tracked.name << "_abaqus_coordinate_error=" << coordinate_error << '\n'
              << "b45_" << tracked.name << "_maximum_normal_angle_degrees=" << maximum_normal_angle << '\n'
              << "b45_" << tracked.name << "_primary_face_warp=" << face_warp << '\n'
              << "b45_" << tracked.name << "_jacobian_directional_error=" << actual.jacobian_error << '\n'
              << "b45_" << tracked.name << "_point_count=" << actual.points << '\n'
              << "b45_" << tracked.name << "_balance=" << actual.balance[0] << ',' << actual.balance[1] << ','
              << actual.balance[2] << '\n';
    return check(
               nodes.size() == tracked.mesh.nodes().size(), "B4.5 compares every prescribed mesh-node displacement") &&
           check(tracked_coordinate_error < 1.0e-7, "B4.5 reads the tracked Exodus geometry and connectivity") &&
           check(coordinate_error < 1.0e-7, "B4.5 uses the tracked Abaqus geometry") &&
           check(generated.name != "warped_bilinear" || face_warp > 0.1,
               "B4.5 warped case has a genuinely noncoplanar primary face") &&
           check(actual.jacobian_error < 1.0e-6,
               "B4.5 active nonplanar contact Jacobians match centered directional differences") &&
           check(qualified_transverse,
               "B4.5 retains the bounded warped-face transverse-force mismatch as a qualified boundary") &&
           check(std::abs(actual.balance[0]) < 1.0e-8 && std::abs(actual.balance[1]) < 1.0e-8 &&
                     std::abs(actual.balance[2]) < 1.0e-8,
               "B4.5 reconstructed contact preserves action-reaction") &&
           check(metrics,
               require_transverse_force
                   ? "B4.5 displacement, all force components, gap, pressure, area, and resultant pass 1 percent"
                   : "B4.5 displacement, normal-dominant force, gap, pressure, area, and resultant pass 1 percent");
}

double contact_reference_difference(const std::string& first_path, const std::string& second_path) {
    const auto first = read_contact(first_path);
    const auto second = read_contact(second_path);
    if (first.size() != second.size()) throw std::invalid_argument("B4.5 smoothing references have different sizes");
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t row = 0; row < first.size(); ++row) {
        if (first[row].secondary != second[row].secondary || first[row].id != second[row].id)
            throw std::invalid_argument("B4.5 smoothing references have different node ordering");
        for (std::size_t component = 0; component < 3; ++component) {
            const double difference = first[row].force[component] - second[row].force[component];
            difference_squared += difference * difference;
            reference_squared += second[row].force[component] * second[row].force[component];
        }
    }
    return std::sqrt(difference_squared / reference_squared);
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--generate") {
        try {
            generate(argv[2]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] B4.5 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc != 11) {
        std::cerr << "Usage: fuelsim_b45_hex8_nonplanar_sts_abaqus_tests <faceted_mesh.e> <faceted_nodes.csv> "
                     "<faceted_contact.csv> <warped_mesh.e> <warped_default_nodes.csv> <warped_default_contact.csv> "
                     "<warped_linear_nodes.csv> <warped_linear_contact.csv> <warped_quadratic_nodes.csv> "
                     "<warped_quadratic_contact.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const std::vector<GeneratedCase> generated = cases();
        bool passed = compare(generated[0], argv[1], argv[2], argv[3], "", true);
        passed = compare(generated[1], argv[4], argv[5], argv[6], "_default", false) && passed;
        passed = compare(generated[1], argv[4], argv[7], argv[8], "_linear", false) && passed;
        passed = compare(generated[1], argv[4], argv[9], argv[10], "_quadratic", false) && passed;
        const double default_linear = contact_reference_difference(argv[6], argv[8]);
        const double default_quadratic = contact_reference_difference(argv[6], argv[10]);
        std::cout << "b45_warped_default_linear_force_relative_difference=" << default_linear << '\n'
                  << "b45_warped_default_quadratic_force_relative_difference=" << default_quadratic << '\n';
        passed = check(default_linear < 1.0e-14,
                     "B4.5 Abaqus default C3D8 finite-sliding transition is exactly linear smoothing") &&
                 check(default_quadratic > 1.0e-4,
                     "B4.5 Abaqus quadratic smoothing measurably changes the warped-face nodal forces") &&
                 passed;
        if (passed) std::cout << "[PASS] B4.5 HEX8 nonplanar surface-to-surface comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.5 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
