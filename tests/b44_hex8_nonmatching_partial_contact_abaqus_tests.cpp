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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct PathStep final {
    const char* name;
    double normal_offset, shift_y, shift_z;
    bool partial;
};

constexpr double gap_quantum = 1.0 / 1024.0;
constexpr std::array<PathStep, 6> path = {{{"FULL_CLOSE", -14.0 * gap_quantum, 0.0, 0.0, false},
    {"PARTIAL_ORIGIN", -6.0 * gap_quantum, 0.0, 0.0, true}, {"PARTIAL_LOW", -6.0 * gap_quantum, 0.125, 0.125, true},
    {"CROSS_PARTIAL", -6.0 * gap_quantum, 0.25, 0.25, true}, {"CROSS_VERTEX", -6.0 * gap_quantum, 0.375, 0.375, true},
    {"FULL_RETURN", -14.0 * gap_quantum, 0.375, 0.375, false}}};

struct NodeReference final {
    std::size_t step, source_node;
    fuelsim::CartesianPoint3 reference;
    std::array<double, 3> displacement{};
};

struct ContactReference final {
    std::size_t step, source_node;
    fuelsim::CartesianPoint3 current;
    std::array<double, 3> normal_force{}, tangential_force{};
    double slip_1, slip_2, gap, pressure;
};

struct ContactCheck final {
    double jacobian_error = 0.0;
    std::array<double, 3> action_reaction{};
    std::size_t differentiated_contributions = 0;
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
    fuelsim::Hex8Element result{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const std::array<double, 3> key = {points[local].x, points[local].y, points[local].z};
        const auto inserted = node_map.emplace(key, nodes.size());
        if (inserted.second) nodes.push_back(points[local]);
        result.nodes[local] = inserted.first->second;
    }
    return result;
}

fuelsim::UnstructuredHex8Mesh generate_mesh() {
    constexpr std::array<double, 5> primary_y = {-0.5, 0.0, 0.5, 1.0, 1.5}, primary_z = {-0.5, 0.0, 0.5, 1.0, 1.5};
    constexpr std::array<double, 4> secondary_y = {0.0, 0.375, 0.625, 1.0};
    constexpr std::array<double, 3> secondary_z = {0.0, 0.375, 1.0};
    constexpr std::array<double, 6> secondary_gaps = {
        gap_quantum, 7.0 * gap_quantum, 3.0 * gap_quantum, 9.0 * gap_quantum, 5.0 * gap_quantum, 11.0 * gap_quantum};
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::vector<fuelsim::ElementSide> primary_contact, primary_contact_patch, secondary_contact;
    std::vector<std::size_t> primary_patch_nodes, secondary_all, secondary_contact_nodes;
    for (std::size_t z = 0; z + 1 < primary_z.size(); ++z)
        for (std::size_t y = 0; y + 1 < primary_y.size(); ++y) {
            const std::size_t element = elements.size();
            elements.push_back(append_cuboid(
                nodes, primary_nodes, 0.0, 1.0, primary_y[y], primary_y[y + 1], primary_z[z], primary_z[z + 1]));
            blocks.push_back(1);
            primary_contact.push_back({element, 1});
            if (y >= 1 && y <= 2 && z >= 1 && z <= 2) {
                primary_contact_patch.push_back({element, 1});
                for (const std::size_t local : std::array<std::size_t, 4>{1, 2, 5, 6})
                    primary_patch_nodes.push_back(elements.back().nodes[local]);
            }
        }
    std::size_t tile = 0;
    for (std::size_t z = 0; z + 1 < secondary_z.size(); ++z)
        for (std::size_t y = 0; y + 1 < secondary_y.size(); ++y) {
            std::map<std::array<double, 3>, std::size_t> tile_nodes;
            const std::size_t element = elements.size();
            const double contact_x = 1.0 + secondary_gaps[tile++];
            elements.push_back(append_cuboid(nodes, tile_nodes, contact_x, contact_x + 1.0, secondary_y[y],
                secondary_y[y + 1], secondary_z[z], secondary_z[z + 1]));
            blocks.push_back(2);
            secondary_contact.push_back({element, 3});
            for (const auto& entry : tile_nodes) secondary_all.push_back(entry.second);
            for (const std::size_t local : std::array<std::size_t, 4>{0, 3, 4, 7})
                secondary_contact_nodes.push_back(elements.back().nodes[local]);
        }
    std::vector<std::size_t> primary_all;
    for (const auto& entry : primary_nodes) primary_all.push_back(entry.second);
    for (std::vector<std::size_t>* values : {&primary_all, &secondary_all, &secondary_contact_nodes})
        std::sort(values->begin(), values->end());
    std::sort(primary_patch_nodes.begin(), primary_patch_nodes.end());
    primary_patch_nodes.erase(
        std::unique(primary_patch_nodes.begin(), primary_patch_nodes.end()), primary_patch_nodes.end());
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all", std::move(primary_all)}, {20, "secondary_all", std::move(secondary_all)},
            {30, "secondary_contact_nodes", std::move(secondary_contact_nodes)},
            {40, "primary_patch_nodes", std::move(primary_patch_nodes)}},
        {{40, "primary_contact", std::move(primary_contact)}, {50, "secondary_contact", std::move(secondary_contact)},
            {60, "primary_contact_patch", std::move(primary_contact_patch)}});
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

std::array<double, 3> secondary_displacement(const fuelsim::CartesianPoint3&, const PathStep& step) {
    return {step.normal_offset, step.shift_y, step.shift_z};
}

void write_abaqus_input(const std::string& path_name, const fuelsim::UnstructuredHex8Mesh& mesh, bool swapped) {
    std::ofstream output(path_name);
    if (!output) throw std::runtime_error("Could not write B4.4 Abaqus input: " + path_name);
    output << std::setprecision(16) << "*Heading\n"
           << "** B4.4 nonmatching C3D8 finite-sliding partial-contact path"
           << (swapped ? " with exchanged primary and secondary.\n" : ".\n")
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n*Node\n";
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        const auto& point = mesh.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    const auto write_block = [&](const char* name, std::int64_t block) {
        output << "*Element, type=C3D8, elset=" << name << '\n';
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            if (mesh.element_block_ids()[element] != block) continue;
            output << element + 1;
            for (const std::size_t node : mesh.elements()[element].nodes) output << ", " << node + 1;
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
           << "*Elset, elset=PRIMARY_PATCH\n";
    for (const fuelsim::ElementSide& side : mesh.side_sets()[2].sides) output << side.element + 1 << ",\n";
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT_PATCH\nPRIMARY_PATCH, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.0\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=CONTACT\n"
           << "*Surface Behavior, pressure-overclosure=LINEAR\n1.e7,\n"
           << "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.\n";
    output << (swapped ? "PRIMARY_CONTACT_PATCH, SECONDARY_CONTACT\n" : "SECONDARY_CONTACT, PRIMARY_CONTACT\n");
    const std::size_t first_step = swapped ? 1 : 0, end_step = swapped ? 2 : path.size();
    for (std::size_t step_index = first_step; step_index < end_step; ++step_index) {
        const PathStep& step = path[step_index];
        output << "*Step, name=" << step.name << ", nlgeom=YES, inc=1\n"
               << "*Static\n1., 1., 1., 1.\n*Boundary, op=NEW\nPRIMARY_ALL, 1, 3, 0.\n";
        for (const std::size_t node : mesh.node_sets()[1].nodes) {
            const std::array<double, 3> displacement = secondary_displacement(mesh.nodes()[node], step);
            for (std::size_t component = 0; component < 3; ++component)
                output << node + 1 << ", " << component + 1 << ", " << component + 1 << ", " << displacement[component]
                       << '\n';
        }
        output << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
    }
}

void generate(const std::string& directory) {
    const fuelsim::UnstructuredHex8Mesh mesh = generate_mesh();
    fuelsim::write_exodus_hex8(directory + "/b44_hex8_nonmatching_partial_mesh.e", mesh);
    write_abaqus_input(directory + "/b44_hex8_nonmatching_partial_contact.inp", mesh, false);
    write_abaqus_input(directory + "/b44_hex8_nonmatching_partial_contact_swapped.inp", mesh, true);
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& file) {
    if (column >= values.size()) throw std::invalid_argument("Incomplete B4.4 CSV row: " + file);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B4.4 CSV number: " + file);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& file) {
    const double value = number(values, column, file);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("Invalid B4.4 CSV index: " + file);
    return static_cast<std::size_t>(value);
}

std::vector<NodeReference> read_nodes(const std::string& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("Could not read B4.4 node reference: " + file);
    std::string line;
    if (!std::getline(input, line) || line != "step,id,x,y,z,ux,uy,uz")
        throw std::invalid_argument("Unexpected B4.4 node CSV header: " + file);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, file), index_value(values, 1, file),
            {number(values, 2, file), number(values, 3, file), number(values, 4, file)},
            {number(values, 5, file), number(values, 6, file), number(values, 7, file)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("Could not read B4.4 contact reference: " + file);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,slip_1,"
                "slip_2,gap,pressure")
        throw std::invalid_argument("Unexpected B4.4 contact CSV header: " + file);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, file), index_value(values, 1, file),
            {number(values, 2, file), number(values, 3, file), number(values, 4, file)},
            {number(values, 5, file), number(values, 6, file), number(values, 7, file)},
            {number(values, 8, file), number(values, 9, file), number(values, 10, file)}, number(values, 11, file),
            number(values, 12, file), number(values, 13, file), number(values, 14, file)});
    }
    return result;
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
}

fuelsim::SpatialDefinition definition(bool swapped = false) {
    fuelsim::SpatialDefinition result;
    result.regions = {{"primary", "primary", material(), 0.0, 300.0, -1, "", fuelsim::StrainFormulation::finite},
        {"secondary", "secondary", material(), 0.0, 300.0, -1, "", fuelsim::StrainFormulation::finite}};
    fuelsim::ContactDefinition contact;
    contact.name = swapped ? "nonmatching_partial_swapped" : "nonmatching_partial";
    contact.primary = swapped ? "secondary_contact" : "primary_contact";
    contact.secondary = swapped ? "primary_contact_patch" : "secondary_contact";
    contact.mechanical = true;
    contact.penalty = 1.0e7;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    result.contacts.push_back(contact);
    return result;
}

std::vector<double> state_for_step(
    const fuelsim::SteadyProblem& problem, const fuelsim::UnstructuredHex8Mesh& mesh, const PathStep& step) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> state = problem.initial_state();
    const fuelsim::Hex8RegionMesh& secondary = spatial.region_mesh(1);
    for (std::size_t local = 0; local < secondary.nodes().size(); ++local) {
        const std::size_t source = secondary.source_node_ids()[local], global = spatial.global_node(1, local);
        const std::array<double, 3> displacement = secondary_displacement(mesh.nodes()[source], step);
        state[spatial.dof(fuelsim::Field::displacement_x, global)] = displacement[0];
        state[spatial.dof(fuelsim::Field::displacement_y, global)] = displacement[1];
        state[spatial.dof(fuelsim::Field::displacement_z, global)] = displacement[2];
    }
    return state;
}

ContactCheck contact_checks(const fuelsim::cartesian::SpatialAssembly& spatial, const std::vector<double>& state) {
    ContactCheck result;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        for (std::size_t row = 0; row < dofs.size(); ++row) {
            const std::size_t global = dofs[row];
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (global >= field.begin && global < field.end) result.action_reaction[component] += residual[row];
            }
        }
        double residual_scale = 0.0;
        for (const double value : residual) residual_scale = std::max(residual_scale, std::abs(value));
        if (residual_scale == 0.0) continue;
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double step = 1.0e-8;
        for (std::size_t index = 0; index < local.size(); ++index) {
            direction[index] = std::sin(static_cast<double>(index + 1));
            plus[index] += step * direction[index];
            minus[index] -= step * direction[index];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                analytic += jacobian[row * local.size() + column] * direction[column];
            const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
            difference_squared += (analytic - reference) * (analytic - reference);
            reference_squared += reference * reference;
        }
        if (reference_squared > 0.0) {
            ++result.differentiated_contributions;
            result.jacobian_error = std::max(result.jacobian_error, std::sqrt(difference_squared / reference_squared));
        }
    }
    return result;
}

std::array<double, 3> cross(const std::array<double, 3>& first, const std::array<double, 3>& second) {
    return {first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0]};
}

double magnitude(const std::array<double, 3>& value) { return std::hypot(value[0], std::hypot(value[1], value[2])); }

void print_metric(const std::string& name, const fuelsim::test::FieldErrorMetrics& value) {
    if (value.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, value);
    else {
        fuelsim::test::print_absolute_metrics(name, value);
        std::cout << name << "_zero_reference_count=" << value.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << value.maximum_zero_reference_difference
                  << '\n';
    }
}

bool run(const std::string& mesh_path, const std::string& node_path, const std::string& contact_path,
    double& partial_resultant, double& partial_reference_resultant) {
    const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(mesh_path);
    const std::vector<NodeReference> node_reference = read_nodes(node_path);
    const std::vector<ContactReference> contact_reference = read_contact(contact_path);
    fuelsim::SteadyProblem problem(definition(), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::vector<std::size_t> contact_sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    fuelsim::test::FieldErrorMetrics displacement_x, displacement_y, displacement_z, normal_x, normal_y, normal_z,
        tangent_y, tangent_z, gap, pressure, total_slip, resultant_x, resultant_y, resultant_z, moment_x, moment_y,
        moment_z, center_x, center_y, center_z;
    double maximum_reference_coordinate_difference = 0.0, maximum_current_coordinate_difference = 0.0;
    bool state_matches = true, all_projected = true, partial_states = true, full_states = true,
         four_face_vertex_visited = false;
    std::set<std::size_t> visited_primary_faces;
    std::vector<std::size_t> first_ownership;
    std::vector<bool> uninterrupted_contact(contact_sources.size(), true);
    bool ownership_changed = false;
    ContactCheck derivative_check;
    for (std::size_t step = 0; step < path.size(); ++step) {
        std::vector<double> state = state_for_step(problem, mesh, path[step]);
        problem.validate_state(state);
        const std::vector<fuelsim::CartesianContactNodeSummary> actual =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        if (actual.size() != contact_sources.size()) throw std::logic_error("B4.4 contact-node count changed");
        std::size_t active = 0;
        std::array<double, 3> actual_resultant{}, reference_resultant{}, actual_moment{}, reference_moment{},
            actual_center_sum{}, reference_center_sum{};
        double actual_weight = 0.0, reference_weight = 0.0;
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const fuelsim::Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
            for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
                const std::size_t source = region_mesh.source_node_ids()[local];
                const auto found = std::find_if(node_reference.begin(), node_reference.end(),
                    [&](const auto& value) { return value.step == step + 1 && value.source_node == source; });
                if (found == node_reference.end()) throw std::logic_error("B4.4 node mapping is incomplete");
                const fuelsim::CartesianPoint3& reference = mesh.nodes()[source];
                maximum_reference_coordinate_difference =
                    std::max({maximum_reference_coordinate_difference, std::abs(reference.x - found->reference.x),
                        std::abs(reference.y - found->reference.y), std::abs(reference.z - found->reference.z)});
                const std::size_t global = spatial.global_node(region, local);
                displacement_x.add(state[spatial.dof(fuelsim::Field::displacement_x, global)], found->displacement[0]);
                displacement_y.add(state[spatial.dof(fuelsim::Field::displacement_y, global)], found->displacement[1]);
                displacement_z.add(state[spatial.dof(fuelsim::Field::displacement_z, global)], found->displacement[2]);
            }
        }
        std::vector<std::size_t> ownership;
        for (std::size_t node = 0; node < actual.size(); ++node) {
            const auto found = std::find_if(contact_reference.begin(), contact_reference.end(), [&](const auto& value) {
                return value.step == step + 1 && value.source_node == contact_sources[node];
            });
            if (found == contact_reference.end()) throw std::logic_error("B4.4 contact mapping is incomplete");
            const auto& value = actual[node];
            all_projected = all_projected && value.projected;
            visited_primary_faces.insert(value.primary_face);
            ownership.push_back(value.primary_face);
            const bool actual_active = value.pressure > 1.0e-8, reference_active = found->pressure > 1.0e-8;
            state_matches = state_matches && actual_active == reference_active;
            if (!actual_active || !reference_active) uninterrupted_contact[node] = false;
            if (actual_active) ++active;
            const std::size_t secondary_local =
                static_cast<std::size_t>(std::find(spatial.region_mesh(1).source_node_ids().begin(),
                                             spatial.region_mesh(1).source_node_ids().end(), contact_sources[node]) -
                                         spatial.region_mesh(1).source_node_ids().begin());
            const std::size_t global = spatial.global_node(1, secondary_local);
            const std::array<double, 3> actual_point = {mesh.nodes()[contact_sources[node]].x +
                                                            state[spatial.dof(fuelsim::Field::displacement_x, global)],
                                            mesh.nodes()[contact_sources[node]].y +
                                                state[spatial.dof(fuelsim::Field::displacement_y, global)],
                                            mesh.nodes()[contact_sources[node]].z +
                                                state[spatial.dof(fuelsim::Field::displacement_z, global)]},
                                        reference_point = {found->current.x, found->current.y, found->current.z};
            for (std::size_t component = 0; component < 3; ++component)
                maximum_current_coordinate_difference = std::max(maximum_current_coordinate_difference,
                    std::abs(actual_point[component] - reference_point[component]));
            four_face_vertex_visited = four_face_vertex_visited || (std::abs(actual_point[1] - 0.5) < 1.0e-14 &&
                                                                       std::abs(actual_point[2] - 0.5) < 1.0e-14);
            std::array<double, 3> actual_force{}, reference_force{};
            for (std::size_t component = 0; component < 3; ++component) {
                actual_force[component] =
                    -value.normal_contact_force[component] - value.tangential_contact_force[component];
                reference_force[component] = found->normal_force[component] + found->tangential_force[component];
                actual_resultant[component] += actual_force[component];
                reference_resultant[component] += reference_force[component];
            }
            normal_x.add(-value.normal_contact_force[0], found->normal_force[0]);
            normal_y.add(-value.normal_contact_force[1], found->normal_force[1]);
            normal_z.add(-value.normal_contact_force[2], found->normal_force[2]);
            tangent_y.add(-value.tangential_contact_force[1], found->tangential_force[1]);
            tangent_z.add(-value.tangential_contact_force[2], found->tangential_force[2]);
            gap.add(value.gap, found->gap);
            pressure.add(value.pressure, found->pressure);
            const double contact_episode_slip =
                actual_active && uninterrupted_contact[node] ? std::hypot(path[step].shift_y, path[step].shift_z) : 0.0;
            total_slip.add(contact_episode_slip, std::hypot(found->slip_1, found->slip_2));
            const std::array<double, 3> actual_node_moment = cross(actual_point, actual_force),
                                        reference_node_moment = cross(reference_point, reference_force);
            const double node_actual_weight = magnitude(actual_force),
                         node_reference_weight = magnitude(reference_force);
            actual_weight += node_actual_weight;
            reference_weight += node_reference_weight;
            for (std::size_t component = 0; component < 3; ++component) {
                actual_moment[component] += actual_node_moment[component];
                reference_moment[component] += reference_node_moment[component];
                actual_center_sum[component] += node_actual_weight * actual_point[component];
                reference_center_sum[component] += node_reference_weight * reference_point[component];
            }
        }
        if (step == 0)
            first_ownership = ownership;
        else
            ownership_changed = ownership_changed || ownership != first_ownership;
        if (path[step].partial)
            partial_states = partial_states && active > 0 && active < actual.size();
        else
            full_states = full_states && active == actual.size();
        resultant_x.add(actual_resultant[0], reference_resultant[0]);
        if (step == 1) {
            partial_resultant = std::abs(actual_resultant[0]);
            partial_reference_resultant = std::abs(reference_resultant[0]);
        }
        resultant_y.add(actual_resultant[1], reference_resultant[1]);
        resultant_z.add(actual_resultant[2], reference_resultant[2]);
        moment_x.add(actual_moment[0], reference_moment[0]);
        moment_y.add(actual_moment[1], reference_moment[1]);
        moment_z.add(actual_moment[2], reference_moment[2]);
        center_x.add(actual_center_sum[0] / actual_weight, reference_center_sum[0] / reference_weight);
        center_y.add(actual_center_sum[1] / actual_weight, reference_center_sum[1] / reference_weight);
        center_z.add(actual_center_sum[2] / actual_weight, reference_center_sum[2] / reference_weight);
        if (step == 1) derivative_check = contact_checks(spatial, state);
        problem.commit_internal_state(state);
    }
    for (const auto& entry : std::array<std::pair<const char*, const fuelsim::test::FieldErrorMetrics*>, 19>{
             {{"displacement_x", &displacement_x}, {"displacement_y", &displacement_y},
                 {"displacement_z", &displacement_z}, {"normal_force_x", &normal_x}, {"normal_force_y", &normal_y},
                 {"normal_force_z", &normal_z}, {"tangent_force_y", &tangent_y}, {"tangent_force_z", &tangent_z},
                 {"gap", &gap}, {"pressure", &pressure}, {"total_slip", &total_slip}, {"resultant_x", &resultant_x},
                 {"resultant_y", &resultant_y}, {"resultant_z", &resultant_z}, {"moment_x", &moment_x},
                 {"moment_y", &moment_y}, {"moment_z", &moment_z}, {"force_center_x", &center_x},
                 {"force_center_y", &center_y}}})
        print_metric(std::string("b44_") + entry.first, *entry.second);
    print_metric("b44_force_center_z", center_z);
    std::cout << "b44_maximum_reference_coordinate_difference=" << maximum_reference_coordinate_difference << '\n'
              << "b44_maximum_current_coordinate_difference=" << maximum_current_coordinate_difference << '\n'
              << "b44_visited_primary_faces=" << visited_primary_faces.size() << '\n'
              << "b44_four_face_common_vertex_visited=" << four_face_vertex_visited << '\n'
              << "b44_jacobian_directional_error=" << derivative_check.jacobian_error << '\n'
              << "b44_differentiated_contributions=" << derivative_check.differentiated_contributions << '\n'
              << "b44_action_reaction=" << derivative_check.action_reaction[0] << ','
              << derivative_check.action_reaction[1] << ',' << derivative_check.action_reaction[2] << '\n';
    constexpr double relative_tolerance = 1.0e-2, force_zero_tolerance = 1.0e-6, displacement_zero_tolerance = 1.0e-12;
    const auto passes = [&](const fuelsim::test::FieldErrorMetrics& metric, double zero_tolerance) {
        return (!metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, relative_tolerance)) &&
               metric.maximum_zero_reference_difference < zero_tolerance;
    };
    return check(node_reference.size() == path.size() * mesh.nodes().size(),
               "B4.4 compares every mesh node at every accepted increment") &&
           check(contact_reference.size() == path.size() * contact_sources.size(),
               "B4.4 compares every secondary contact node at every accepted increment") &&
           check(maximum_reference_coordinate_difference < 1.0e-14 && maximum_current_coordinate_difference < 1.0e-14,
               "B4.4 uses the exact tracked mesh and current coordinates") &&
           check(all_projected && partial_states && full_states && state_matches,
               "B4.4 distinguishes full and partial physical contact while every constraint remains projected") &&
           check(visited_primary_faces.size() >= 4 && ownership_changed && four_face_vertex_visited,
               "B4.4 finite sliding transfers ownership through a four-face common primary vertex") &&
           check(derivative_check.differentiated_contributions > 0 && derivative_check.jacobian_error < 1.0e-5,
               "B4.4 active partial-contact Jacobians match centered directional differences away from transitions") &&
           check(std::abs(derivative_check.action_reaction[0]) < force_zero_tolerance &&
                     std::abs(derivative_check.action_reaction[1]) < force_zero_tolerance &&
                     std::abs(derivative_check.action_reaction[2]) < force_zero_tolerance,
               "B4.4 contact residuals preserve three-component action-reaction balance") &&
           check(passes(displacement_x, displacement_zero_tolerance) &&
                     passes(displacement_y, displacement_zero_tolerance) &&
                     passes(displacement_z, displacement_zero_tolerance),
               "B4.4 all displacement fields pass three metrics and separate zero-reference checks") &&
           check(passes(normal_x, force_zero_tolerance) && passes(normal_y, force_zero_tolerance) &&
                     passes(normal_z, force_zero_tolerance) && passes(tangent_y, force_zero_tolerance) &&
                     passes(tangent_z, force_zero_tolerance),
               "B4.4 global normal and both local tangent force fields agree with Abaqus below 1 percent") &&
           check(passes(gap, displacement_zero_tolerance) && passes(pressure, force_zero_tolerance) &&
                     passes(total_slip, displacement_zero_tolerance),
               "B4.4 gap, pressure, and total-slip fields agree with Abaqus below 1 percent") &&
           check(passes(resultant_x, force_zero_tolerance) && passes(resultant_y, force_zero_tolerance) &&
                     passes(resultant_z, force_zero_tolerance) && passes(moment_x, force_zero_tolerance) &&
                     passes(moment_y, force_zero_tolerance) && passes(moment_z, force_zero_tolerance),
               "B4.4 resultants and moments agree with Abaqus below 1 percent") &&
           check(passes(center_x, displacement_zero_tolerance) && passes(center_y, displacement_zero_tolerance) &&
                     passes(center_z, displacement_zero_tolerance),
               "B4.4 contact force centers agree with Abaqus below 1 percent");
}

bool run_swapped(const fuelsim::UnstructuredHex8Mesh& mesh, const std::string& node_path,
    const std::string& contact_path, double base_resultant, double base_reference_resultant) {
    const std::vector<NodeReference> node_reference = read_nodes(node_path);
    const std::vector<ContactReference> contact_reference = read_contact(contact_path);
    fuelsim::SteadyProblem problem(definition(true), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::vector<std::size_t> contact_sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    std::vector<double> state = state_for_step(problem, mesh, path[1]);
    problem.validate_state(state);
    const std::vector<fuelsim::CartesianContactNodeSummary> actual =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
    fuelsim::test::FieldErrorMetrics displacement_x, displacement_y, displacement_z, normal_x, normal_y, normal_z,
        tangent_y, tangent_z, gap, pressure, total_slip, resultant_x, resultant_y, resultant_z, moment_x, moment_y,
        moment_z, center_x, center_y, center_z;
    double maximum_reference_coordinate_difference = 0.0, maximum_current_coordinate_difference = 0.0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const fuelsim::Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const auto found = std::find_if(node_reference.begin(), node_reference.end(),
                [&](const auto& value) { return value.step == 1 && value.source_node == source; });
            if (found == node_reference.end()) throw std::logic_error("B4.4 swapped node mapping is incomplete");
            const fuelsim::CartesianPoint3& reference = mesh.nodes()[source];
            maximum_reference_coordinate_difference =
                std::max({maximum_reference_coordinate_difference, std::abs(reference.x - found->reference.x),
                    std::abs(reference.y - found->reference.y), std::abs(reference.z - found->reference.z)});
            const std::size_t global = spatial.global_node(region, local);
            displacement_x.add(state[spatial.dof(fuelsim::Field::displacement_x, global)], found->displacement[0]);
            displacement_y.add(state[spatial.dof(fuelsim::Field::displacement_y, global)], found->displacement[1]);
            displacement_z.add(state[spatial.dof(fuelsim::Field::displacement_z, global)], found->displacement[2]);
        }
    }
    std::array<double, 3> actual_resultant{}, reference_resultant{}, actual_moment{}, reference_moment{},
        actual_center_sum{}, reference_center_sum{};
    double actual_weight = 0.0, reference_weight = 0.0;
    std::size_t active = 0, reference_active_count = 0;
    std::size_t fuelsim_positive_gap_pressure = 0, abaqus_positive_gap_pressure = 0;
    double fuelsim_penalty_law_difference = 0.0, abaqus_penalty_law_difference = 0.0;
    std::set<std::size_t> visited_primary_faces;
    bool state_matches = true, all_projected = true;
    for (std::size_t node = 0; node < actual.size(); ++node) {
        const auto found = std::find_if(contact_reference.begin(), contact_reference.end(),
            [&](const auto& value) { return value.step == 1 && value.source_node == contact_sources[node]; });
        if (found == contact_reference.end()) throw std::logic_error("B4.4 swapped contact mapping is incomplete");
        const auto& value = actual[node];
        all_projected = all_projected && value.projected;
        visited_primary_faces.insert(value.primary_face);
        const bool actual_active = value.pressure > 1.0e-8, reference_active = found->pressure > 1.0e-8;
        state_matches = state_matches && actual_active == reference_active;
        if (actual_active) ++active;
        if (reference_active) ++reference_active_count;
        const std::size_t secondary_local =
            static_cast<std::size_t>(std::find(spatial.region_mesh(0).source_node_ids().begin(),
                                         spatial.region_mesh(0).source_node_ids().end(), contact_sources[node]) -
                                     spatial.region_mesh(0).source_node_ids().begin());
        const std::size_t global = spatial.global_node(0, secondary_local);
        const std::array<double, 3> actual_point = {
            mesh.nodes()[contact_sources[node]].x + state[spatial.dof(fuelsim::Field::displacement_x, global)],
            mesh.nodes()[contact_sources[node]].y + state[spatial.dof(fuelsim::Field::displacement_y, global)],
            mesh.nodes()[contact_sources[node]].z + state[spatial.dof(fuelsim::Field::displacement_z, global)]};
        const std::array<double, 3> reference_point = {found->current.x, found->current.y, found->current.z};
        for (std::size_t component = 0; component < 3; ++component)
            maximum_current_coordinate_difference = std::max(
                maximum_current_coordinate_difference, std::abs(actual_point[component] - reference_point[component]));
        std::array<double, 3> actual_force{}, reference_force{};
        for (std::size_t component = 0; component < 3; ++component) {
            actual_force[component] =
                -value.normal_contact_force[component] - value.tangential_contact_force[component];
            reference_force[component] = found->normal_force[component] + found->tangential_force[component];
            actual_resultant[component] += actual_force[component];
            reference_resultant[component] += reference_force[component];
        }
        normal_x.add(-value.normal_contact_force[0], found->normal_force[0]);
        normal_y.add(-value.normal_contact_force[1], found->normal_force[1]);
        normal_z.add(-value.normal_contact_force[2], found->normal_force[2]);
        tangent_y.add(-value.tangential_contact_force[1], found->tangential_force[1]);
        tangent_z.add(-value.tangential_contact_force[2], found->tangential_force[2]);
        gap.add(value.gap, found->gap);
        pressure.add(value.pressure, found->pressure);
        fuelsim_penalty_law_difference =
            std::max(fuelsim_penalty_law_difference, std::abs(value.pressure - std::max(-1.0e7 * value.gap, 0.0)));
        abaqus_penalty_law_difference =
            std::max(abaqus_penalty_law_difference, std::abs(found->pressure - std::max(-1.0e7 * found->gap, 0.0)));
        if (value.gap > 0.0 && value.pressure > 0.0) ++fuelsim_positive_gap_pressure;
        if (found->gap > 0.0 && found->pressure > 0.0) ++abaqus_positive_gap_pressure;
        total_slip.add(0.0, std::hypot(found->slip_1, found->slip_2));
        const std::array<double, 3> actual_node_moment = cross(actual_point, actual_force),
                                    reference_node_moment = cross(reference_point, reference_force);
        const double node_actual_weight = magnitude(actual_force), node_reference_weight = magnitude(reference_force);
        actual_weight += node_actual_weight;
        reference_weight += node_reference_weight;
        for (std::size_t component = 0; component < 3; ++component) {
            actual_moment[component] += actual_node_moment[component];
            reference_moment[component] += reference_node_moment[component];
            actual_center_sum[component] += node_actual_weight * actual_point[component];
            reference_center_sum[component] += node_reference_weight * reference_point[component];
        }
    }
    resultant_x.add(actual_resultant[0], reference_resultant[0]);
    resultant_y.add(actual_resultant[1], reference_resultant[1]);
    resultant_z.add(actual_resultant[2], reference_resultant[2]);
    moment_x.add(actual_moment[0], reference_moment[0]);
    moment_y.add(actual_moment[1], reference_moment[1]);
    moment_z.add(actual_moment[2], reference_moment[2]);
    center_x.add(actual_center_sum[0] / actual_weight, reference_center_sum[0] / reference_weight);
    center_y.add(actual_center_sum[1] / actual_weight, reference_center_sum[1] / reference_weight);
    center_z.add(actual_center_sum[2] / actual_weight, reference_center_sum[2] / reference_weight);
    for (const auto& entry : std::array<std::pair<const char*, const fuelsim::test::FieldErrorMetrics*>, 20>{
             {{"displacement_x", &displacement_x}, {"displacement_y", &displacement_y},
                 {"displacement_z", &displacement_z}, {"normal_force_x", &normal_x}, {"normal_force_y", &normal_y},
                 {"normal_force_z", &normal_z}, {"tangent_force_y", &tangent_y}, {"tangent_force_z", &tangent_z},
                 {"gap", &gap}, {"pressure", &pressure}, {"total_slip", &total_slip}, {"resultant_x", &resultant_x},
                 {"resultant_y", &resultant_y}, {"resultant_z", &resultant_z}, {"moment_x", &moment_x},
                 {"moment_y", &moment_y}, {"moment_z", &moment_z}, {"force_center_x", &center_x},
                 {"force_center_y", &center_y}, {"force_center_z", &center_z}}})
        print_metric(std::string("b44_swapped_") + entry.first, *entry.second);
    const double swapped_resultant = std::abs(actual_resultant[0]),
                 swapped_reference_resultant = std::abs(reference_resultant[0]),
                 fuelsim_bias = std::abs(swapped_resultant - base_resultant) / base_resultant,
                 abaqus_bias =
                     std::abs(swapped_reference_resultant - base_reference_resultant) / base_reference_resultant;
    std::cout << "b44_swapped_maximum_reference_coordinate_difference=" << maximum_reference_coordinate_difference
              << '\n'
              << "b44_swapped_maximum_current_coordinate_difference=" << maximum_current_coordinate_difference << '\n'
              << "b44_swapped_visited_primary_faces=" << visited_primary_faces.size() << '\n'
              << "b44_swapped_fuelsim_penalty_law_maximum_absolute_difference=" << fuelsim_penalty_law_difference
              << '\n'
              << "b44_swapped_abaqus_penalty_law_maximum_absolute_difference=" << abaqus_penalty_law_difference << '\n'
              << "b44_swapped_fuelsim_positive_gap_pressure_count=" << fuelsim_positive_gap_pressure << '\n'
              << "b44_swapped_abaqus_positive_gap_pressure_count=" << abaqus_positive_gap_pressure << '\n'
              << "b44_fuelsim_one_sided_resultant_bias=" << fuelsim_bias << '\n'
              << "b44_abaqus_one_sided_resultant_bias=" << abaqus_bias << '\n';
    constexpr double relative_tolerance = 1.0e-2, force_zero_tolerance = 1.0e-6, displacement_zero_tolerance = 1.0e-12;
    const auto passes = [&](const fuelsim::test::FieldErrorMetrics& metric, double zero_tolerance) {
        return (!metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, relative_tolerance)) &&
               metric.maximum_zero_reference_difference < zero_tolerance;
    };
    return check(node_reference.size() == mesh.nodes().size() && contact_reference.size() == contact_sources.size(),
               "B4.4 swapped comparison covers every mesh and secondary contact node") &&
           check(maximum_reference_coordinate_difference < 1.0e-14 && maximum_current_coordinate_difference < 1.0e-14,
               "B4.4 swapped Abaqus and Fuelsim coordinates agree") &&
           check(all_projected && active == reference_active_count && active > 0 && active < actual.size() &&
                     visited_primary_faces.size() >= 2 && state_matches,
               "B4.4 swapped designation reproduces the partially active recovered nodal field") &&
           check(passes(displacement_x, displacement_zero_tolerance) &&
                     passes(displacement_y, displacement_zero_tolerance) &&
                     passes(displacement_z, displacement_zero_tolerance),
               "B4.4 swapped full-mesh displacement fields agree below 1 percent") &&
           check(passes(normal_y, force_zero_tolerance) && passes(normal_z, force_zero_tolerance) &&
                     passes(tangent_y, force_zero_tolerance) && passes(tangent_z, force_zero_tolerance) &&
                     passes(total_slip, displacement_zero_tolerance),
               "B4.4 swapped theoretical-zero transverse forces and slip remain bounded") &&
           check(fuelsim_penalty_law_difference > 1.0e4 && abaqus_penalty_law_difference > 1.0e4 &&
                     fuelsim_positive_gap_pressure == abaqus_positive_gap_pressure && abaqus_positive_gap_pressure > 0,
               "B4.4 distinguishes recovered nodal pressure from the internal constraint penalty law") &&
           check(passes(normal_x, force_zero_tolerance) && passes(gap, displacement_zero_tolerance) &&
                     passes(pressure, force_zero_tolerance) && passes(resultant_x, force_zero_tolerance) &&
                     passes(moment_y, force_zero_tolerance) && passes(moment_z, force_zero_tolerance) &&
                     passes(center_x, displacement_zero_tolerance) && passes(center_y, displacement_zero_tolerance) &&
                     passes(center_z, displacement_zero_tolerance),
               "B4.4 swapped nodal force, gap, recovered pressure, resultants, moments, and force center pass 1 "
               "percent") &&
           check(fuelsim_bias > 0.4 && abaqus_bias > 0.4 && std::abs(fuelsim_bias - abaqus_bias) < 1.0e-12,
               "B4.4 reproduces the exchanged-side discretization bias instead of qualifying it");
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--generate") {
        try {
            generate(argv[2]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] B4.4 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc != 6) {
        std::cerr << "Usage: fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests <mesh.e> <nodes.csv> "
                     "<contact.csv> <swapped_nodes.csv> <swapped_contact.csv>\n"
                     "   or: fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests --generate "
                     "<abaqus-directory>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        double partial_resultant = 0.0, partial_reference_resultant = 0.0;
        const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(argv[1]);
        const bool base_passed = run(argv[1], argv[2], argv[3], partial_resultant, partial_reference_resultant);
        const bool passed =
            base_passed && run_swapped(mesh, argv[4], argv[5], partial_resultant, partial_reference_resultant);
        if (passed) std::cout << "[PASS] B4.4 HEX8 nonmatching finite-sliding partial-contact comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.4 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
