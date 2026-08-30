#include "fuelsim/solver/solve_workflows.hpp"
#include "support/abaqus_hex8_full_field.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double time_step = 0.02;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::size_t increment_count = 20;

enum class ScanPath { monotonic, contact_cycle, friction_reversal, nonmatching_contact_cycle };

struct ScanParameters final {
    double penalty = 1.0e9, friction_coefficient = 0.05, slip_tolerance = 0.005, conductance = 50.0,
           pressure_conductance = 0.001, primary_poisson = 0.28, secondary_poisson = 0.3, bending_traction = 0.0;
    double initial_gap = 0.0;
    double normal_displacement_scale = 1.0, tangential_displacement_scale = 1.0;
    bool traction_controlled = false, elastic_only = false, anchor_bending = false, reduced_integration = false,
         nonmatching_mesh = false;
    ScanPath path = ScanPath::monotonic;
};

struct StepSnapshot final {
    double time = 0.0, load_factor = 0.0;
    std::vector<double> state;
    std::vector<std::vector<fuelsim::CartesianMaterialPointState>> material;
    std::vector<fuelsim::CartesianContactNodeSummary> contact;
    std::vector<fuelsim::ContactPointHistory> contact_history;
    fuelsim::TransientConservationSummary conservation;
};

struct ContactReference final {
    std::size_t increment = 0, node = 0, state = 0;
    double time = 0.0;
    fuelsim::CartesianPoint3 position{};
    double opening = 0.0, pressure = 0.0, slip_first = 0.0, slip_second = 0.0;
    std::array<double, 3> normal_force{}, shear_force{};
    double heat_flux = 0.0;
};

struct EnergyReference final {
    std::size_t increment = 0;
    double time = 0.0, internal = 0.0, elastic = 0.0, plastic = 0.0, creep = 0.0, friction = 0.0, external_work = 0.0,
           boundary_heat_rate = 0.0;
};

struct NodeReference final {
    std::size_t increment = 0, node = 0;
    double time = 0.0;
    std::array<double, 8> fields{};
};

struct IntegrationReference final {
    std::size_t increment = 0, element = 0, point = 0;
    double time = 0.0;
    fuelsim::CartesianPoint3 position{};
    double temperature = 0.0;
    std::array<double, 3> heat_flux{};
    fuelsim::SymmetricTensor3Values stress{}, logarithmic_strain{}, elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0, integration_volume = 0.0;
};

class SnapshotObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        StepSnapshot snapshot;
        snapshot.time = step.time;
        snapshot.load_factor = step.load_factor;
        snapshot.state = problem.committed_solution();
        snapshot.conservation = step.conservation;
        for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
            const std::size_t elements =
                fuelsim::cartesian::ProblemAccess::region_mesh(problem, region).elements().size();
            for (std::size_t element = 0; element < elements; ++element) {
                const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, region, element);
                snapshot.material.emplace_back(history.begin(), history.end());
            }
        }
        snapshot.contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, snapshot.state);
        snapshot.contact_history = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
        _snapshots.push_back(std::move(snapshot));
    }

    const std::vector<StepSnapshot>& snapshots() const noexcept { return _snapshots; }

  private:
    std::vector<StepSnapshot> _snapshots;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.23 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.23 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.23 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "increment,time_s,node,temperature_k,u1_m,u2_m,u3_m,reaction_heat_flux_w,rf1_n,rf2_n,rf3_n")
        throw std::invalid_argument("Unexpected Abaqus B5.23 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11) throw std::invalid_argument("Unexpected Abaqus B5.23 nodal column count in " + path);
        NodeReference reference;
        reference.increment = positive_integer(number(values, 0, path), path);
        reference.time = number(values, 1, path);
        reference.node = positive_integer(number(values, 2, path), path);
        for (std::size_t field = 0; field < reference.fields.size(); ++field) {
            reference.fields[field] = number(values, field + 3, path);
            if (std::abs(reference.fields[field]) < 1.0e-20) reference.fields[field] = 0.0;
        }
        result.push_back(reference);
    }
    if (result.size() != increment_count * 24)
        throw std::invalid_argument("Abaqus B5.23 nodal reference must contain 480 rows");
    return result;
}

fuelsim::SymmetricTensor3Values tensor(
    const std::vector<std::string>& values, std::size_t start, const std::string& path, bool engineering_shear) {
    const double shear_scale = engineering_shear ? 0.5 : 1.0;
    return {number(values, start, path), number(values, start + 1, path), number(values, start + 2, path),
        shear_scale * number(values, start + 3, path), shear_scale * number(values, start + 5, path),
        shear_scale * number(values, start + 4, path)};
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.23 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected =
        "increment,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl1_w_m2,hfl2_w_m2,hfl3_w_m2,"
        "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,le11,le22,le33,le12_engineering,le13_engineering,"
        "le23_engineering,ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,pe11,pe22,pe33,"
        "pe12_engineering,pe13_engineering,pe23_engineering,peeq,ce11,ce22,ce33,ce12_engineering,"
        "ce13_engineering,ce23_engineering,ceeq,ivol_m3";
    if (line != expected) throw std::invalid_argument("Unexpected Abaqus B5.23 integration header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 44)
            throw std::invalid_argument("Unexpected Abaqus B5.23 integration column count in " + path);
        IntegrationReference reference;
        reference.increment = positive_integer(number(values, 0, path), path);
        reference.time = number(values, 1, path);
        reference.element = positive_integer(number(values, 2, path), path);
        reference.point = positive_integer(number(values, 3, path), path);
        reference.position = {number(values, 4, path), number(values, 5, path), number(values, 6, path)};
        reference.temperature = number(values, 7, path);
        reference.heat_flux = {number(values, 8, path), number(values, 9, path), number(values, 10, path)};
        reference.stress = tensor(values, 11, path, false);
        reference.logarithmic_strain = tensor(values, 17, path, true);
        reference.elastic_strain = tensor(values, 23, path, true);
        reference.plastic_strain = tensor(values, 29, path, true);
        reference.equivalent_plastic_strain = number(values, 35, path);
        reference.creep_strain = tensor(values, 36, path, true);
        reference.equivalent_creep_strain = number(values, 42, path);
        reference.integration_volume = number(values, 43, path);
        result.push_back(reference);
    }
    if (result.size() != increment_count * 4 * 8)
        throw std::invalid_argument("Abaqus B5.23 integration reference must contain 640 rows");
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.23 contact: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected =
        "increment,time_s,node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,normal_force1_n,"
        "normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,contact_heat_flux_w,state";
    if (line != expected) throw std::invalid_argument("Unexpected Abaqus B5.23 contact header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 18) throw std::invalid_argument("Unexpected Abaqus B5.23 contact column count in " + path);
        ContactReference value;
        value.increment = positive_integer(number(values, 0, path), path);
        value.time = number(values, 1, path);
        value.node = positive_integer(number(values, 2, path), path);
        value.position = {number(values, 3, path), number(values, 4, path), number(values, 5, path)};
        value.opening = number(values, 6, path);
        value.pressure = number(values, 7, path);
        value.slip_first = number(values, 8, path);
        value.slip_second = number(values, 9, path);
        for (std::size_t component = 0; component < 3; ++component) {
            value.normal_force[component] = number(values, 10 + component, path);
            value.shear_force[component] = number(values, 13 + component, path);
        }
        value.heat_flux = number(values, 16, path);
        value.state = positive_integer(number(values, 17, path), path);
        result.push_back(value);
    }
    if (result.size() != increment_count * 4)
        throw std::invalid_argument("Abaqus B5.23 contact reference must contain 80 rows");
    return result;
}

std::vector<EnergyReference> read_energy(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.23 energy: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "increment,time_s,allie_j,allse_j,allpd_j,allcd_j,allfd_j,allwk_j,boundary_heat_rate_w")
        throw std::invalid_argument("Unexpected Abaqus B5.23 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 9) throw std::invalid_argument("Unexpected Abaqus B5.23 energy column count in " + path);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            number(values, 2, path), number(values, 3, path), number(values, 4, path), number(values, 5, path),
            number(values, 6, path), number(values, 7, path), number(values, 8, path)});
    }
    if (result.size() != increment_count)
        throw std::invalid_argument("Abaqus B5.23 energy reference must contain 20 rows");
    return result;
}

std::array<double, 6> components(const fuelsim::SymmetricTensor3Values& value) {
    return {value.xx, value.yy, value.zz, value.xy, value.yz, value.xz};
}

using Matrix3 = std::array<std::array<double, 3>, 3>;

fuelsim::SymmetricTensor3Values logarithmic_strain(
    const fuelsim::Hex8QuadraturePoint& point, const fuelsim::Hex8LocalValues& state) {
    Matrix3 deformation = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                deformation[component][direction] +=
                    state[8 * (component + 1) + node] * point.gradient[node][direction];
    Matrix3 left{};
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t column = 0; column < 3; ++column)
            for (std::size_t inner = 0; inner < 3; ++inner)
                left[row][column] += deformation[row][inner] * deformation[column][inner];
    Matrix3 vectors = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (std::size_t sweep = 0; sweep < 32; ++sweep) {
        std::size_t p = 0, q = 1;
        double largest = std::abs(left[p][q]);
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t column = row + 1; column < 3; ++column)
                if (std::abs(left[row][column]) > largest) {
                    p = row;
                    q = column;
                    largest = std::abs(left[row][column]);
                }
        const double scale = std::max({1.0, std::abs(left[0][0]), std::abs(left[1][1]), std::abs(left[2][2])});
        if (largest <= 1.0e-15 * scale) break;
        const double tau = (left[q][q] - left[p][p]) / (2.0 * left[p][q]);
        const double tangent = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent), sine = tangent * cosine;
        const double app = left[p][p], aqq = left[q][q], apq = left[p][q];
        left[p][p] = app - tangent * apq;
        left[q][q] = aqq + tangent * apq;
        left[p][q] = left[q][p] = 0.0;
        for (std::size_t row = 0; row < 3; ++row) {
            if (row == p || row == q) continue;
            const double arp = left[row][p], arq = left[row][q];
            left[row][p] = left[p][row] = cosine * arp - sine * arq;
            left[row][q] = left[q][row] = sine * arp + cosine * arq;
        }
        for (std::size_t row = 0; row < 3; ++row) {
            const double vrp = vectors[row][p], vrq = vectors[row][q];
            vectors[row][p] = cosine * vrp - sine * vrq;
            vectors[row][q] = sine * vrp + cosine * vrq;
        }
    }
    Matrix3 result{};
    for (std::size_t mode = 0; mode < 3; ++mode) {
        if (!(left[mode][mode] > 0.0))
            throw std::domain_error("B5.23 logarithmic strain requires a positive left stretch tensor");
        const double value = 0.5 * std::log(left[mode][mode]);
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t column = 0; column < 3; ++column)
                result[row][column] += value * vectors[row][mode] * vectors[column][mode];
    }
    return {result[0][0], result[1][1], result[2][2], result[0][1], result[1][2], result[0][2]};
}

fuelsim::UnstructuredHex8Mesh mesh(std::size_t through_thickness_elements = 1, std::size_t tangential_elements = 2,
    double distortion = 0.0, double initial_gap = 0.0, bool nonmatching = false) {
    if (through_thickness_elements == 0 || tangential_elements == 0)
        throw std::invalid_argument("B5.23 mesh divisions must be positive");
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> element_blocks;
    std::array<std::vector<fuelsim::ElementSide>, 10> faces;
    const std::array<std::size_t, 2> tangential_divisions =
        nonmatching ? std::array<std::size_t, 2>{2, 1}
                    : std::array<std::size_t, 2>{tangential_elements, tangential_elements};
    const std::array<double, 2> y_lower =
                                    nonmatching ? std::array<double, 2>{0.0, 0.1} : std::array<double, 2>{0.0, 0.0},
                                y_upper =
                                    nonmatching ? std::array<double, 2>{2.0, 0.9} : std::array<double, 2>{1.0, 1.0},
                                z_lower =
                                    nonmatching ? std::array<double, 2>{-1.0, 0.1} : std::array<double, 2>{0.0, 0.0},
                                z_upper =
                                    nonmatching ? std::array<double, 2>{2.0, 0.9} : std::array<double, 2>{1.0, 1.0};
    const std::array<std::size_t, 2> block_offsets = {
        0, 2 * (tangential_divisions[0] + 1) * (through_thickness_elements + 1)};
    const auto node = [&](std::size_t block, std::size_t z, std::size_t y, std::size_t x) {
        return block_offsets[block] + z * (tangential_divisions[block] + 1) * (through_thickness_elements + 1) +
               y * (through_thickness_elements + 1) + x;
    };
    for (std::size_t block = 0; block < 2; ++block)
        for (std::size_t z = 0; z < 2; ++z)
            for (std::size_t y = 0; y <= tangential_divisions[block]; ++y)
                for (std::size_t x = 0; x <= through_thickness_elements; ++x)
                    nodes.push_back({static_cast<double>(block) + (block == 1 ? initial_gap : 0.0) +
                                         static_cast<double>(x) / static_cast<double>(through_thickness_elements) +
                                         distortion * std::sin(pi * static_cast<double>(y) /
                                                               static_cast<double>(tangential_divisions[block])),
                        y_lower[block] + (y_upper[block] - y_lower[block]) * static_cast<double>(y) /
                                             static_cast<double>(tangential_divisions[block]),
                        z == 0 ? z_lower[block] : z_upper[block]});
    for (std::size_t block = 0; block < 2; ++block)
        for (std::size_t y = 0; y < tangential_divisions[block]; ++y)
            for (std::size_t x = 0; x < through_thickness_elements; ++x) {
                const std::size_t element_index = elements.size();
                elements.push_back({{{node(block, 0, y, x), node(block, 0, y, x + 1), node(block, 0, y + 1, x + 1),
                    node(block, 0, y + 1, x), node(block, 1, y, x), node(block, 1, y, x + 1),
                    node(block, 1, y + 1, x + 1), node(block, 1, y + 1, x)}}});
                element_blocks.push_back(block == 0 ? 1 : 2);
                if (block == 0 && x == 0) faces[0].push_back({element_index, 3});
                if (block == 0 && x + 1 == through_thickness_elements) faces[1].push_back({element_index, 1});
                if (block == 1 && x == 0) faces[2].push_back({element_index, 3});
                if (block == 1 && x + 1 == through_thickness_elements) faces[3].push_back({element_index, 1});
                if (block == 1 && x + 1 == through_thickness_elements && 2 * y < tangential_divisions[block])
                    faces[6].push_back({element_index, 1});
                if (block == 1 && x + 1 == through_thickness_elements && 2 * y >= tangential_divisions[block])
                    faces[7].push_back({element_index, 1});
                if (block == 1 && y == 0) faces[4].push_back({element_index, 0});
                if (block == 1 && y + 1 == tangential_divisions[block]) faces[8].push_back({element_index, 2});
                if (block == 1) faces[5].push_back({element_index, 4});
                if (block == 1) faces[9].push_back({element_index, 5});
            }

    std::vector<std::size_t> primary_outer, secondary_outer, secondary_y0, secondary_z0;
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t y = 0; y <= tangential_divisions[0]; ++y) primary_outer.push_back(node(0, z, y, 0));
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t y = 0; y <= tangential_divisions[1]; ++y)
            secondary_outer.push_back(node(1, z, y, through_thickness_elements));
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t x = 0; x <= through_thickness_elements; ++x) secondary_y0.push_back(node(1, z, 0, x));
    for (std::size_t y = 0; y <= tangential_divisions[1]; ++y)
        for (std::size_t x = 0; x <= through_thickness_elements; ++x) secondary_z0.push_back(node(1, 0, y, x));
    std::vector<fuelsim::ElementSide> secondary_all_surface;
    for (const std::size_t face : {2U, 3U, 4U, 5U, 8U, 9U})
        secondary_all_surface.insert(secondary_all_surface.end(), faces[face].begin(), faces[face].end());
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(element_blocks),
        {{1, "primary"}, {2, "secondary"}},
        {{11, "primary_outer_nodes", primary_outer}, {12, "secondary_outer_nodes", secondary_outer},
            {13, "secondary_y0", secondary_y0}, {14, "secondary_z0", secondary_z0}},
        {{21, "primary_outer", faces[0]}, {22, "primary_contact", faces[1]}, {23, "secondary_contact", faces[2]},
            {24, "secondary_outer", faces[3]}, {25, "secondary_y0_surface", faces[4]},
            {26, "secondary_z0_surface", faces[5]}, {27, "secondary_outer_lower", faces[6]},
            {28, "secondary_outer_upper", faces[7]}, {29, "secondary_all_surface", secondary_all_surface}});
}

fuelsim::ThermoelasticProperties primary_material(double poisson = 0.28) {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(0.0, 15.0, 1.2e8, poisson, 8.0e-6, 300.0, -1.0e5, 0.0, 0.0, 1.0, 1.0);
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>(*result.functions);
    functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
        {{"conductivity", 15.0}, {"density", 100.0}, {"specific_heat", 1.0}, {"reference_temperature", 300.0},
            {"conductivity_temperature_coefficient", 0.015}, {"density_temperature_coefficient", 0.0},
            {"specific_heat_temperature_coefficient", 0.001}});
    result.functions = std::move(functions);
    return result;
}

fuelsim::ThermoelasticProperties secondary_material(double poisson = 0.3, bool elastic_only = false) {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(0.0, 10.0, 1.0e8, poisson, 1.0e-5, 300.0, -5.0e4, 0.0, 0.0, 1.0, 1.0);
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>(*result.functions);
    functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
        {{"conductivity", 10.0}, {"density", 100.0}, {"specific_heat", 1.0}, {"reference_temperature", 300.0},
            {"conductivity_temperature_coefficient", 0.01}, {"density_temperature_coefficient", 0.0},
            {"specific_heat_temperature_coefficient", 0.001}});
    result.functions = std::move(functions);
    if (elastic_only) return result;
    result = fuelsim::test::with_plasticity(std::move(result), 2.0e5, 1.0e7, 300.0, -100.0, -5000.0);
    return fuelsim::test::with_norton(std::move(result), 1.0e-4, 2.0e5, 3.0, 300.0, 2.0e-7);
}

fuelsim::BoundaryConditionDefinition configured_boundary(const std::string& name, fuelsim::BoundaryConditionType type,
    const std::string& boundary, fuelsim::Field field, const std::string& function, bool use_displaced_geometry) {
    fuelsim::BoundaryConditionDefinition result{name, type, boundary, field, 1.0, false, function};
    result.use_displaced_geometry = use_displaced_geometry;
    result.configuration_explicit = true;
    return result;
}

fuelsim::SpatialDefinition definition(const ScanParameters& parameters = {}) {
    fuelsim::SpatialDefinition result;
    result.regions = {{"primary", "primary", primary_material(parameters.primary_poisson), 0.0, 300.0, -1, "",
                          fuelsim::StrainFormulation::finite},
        {"secondary", "secondary", secondary_material(parameters.secondary_poisson, parameters.elastic_only), 0.0,
            300.0, -1, "", fuelsim::StrainFormulation::finite}};
    if (parameters.reduced_integration)
        for (fuelsim::RegionDefinition& region : result.regions)
            region.hex8_element_formulation = fuelsim::Hex8ElementFormulation::c3d8rt;
    const std::vector<double> times =
        parameters.path == ScanPath::monotonic ? std::vector<double>{0.0, 0.4}
        : parameters.path == ScanPath::nonmatching_contact_cycle
            ? std::vector<double>{0.0, 0.08, 0.10, 0.18, 0.20, 0.40, 0.48, 0.50, 0.58, 0.60, 0.80, 0.88, 0.90}
            : std::vector<double>{0.0, 0.1, 0.2, 0.3, 0.4};
    if (parameters.path == ScanPath::monotonic) {
        result.time_tables.emplace_back("secondary_temperature", times, std::vector<double>{300.0, 500.0});
        result.time_tables.emplace_back("pressure", times, std::vector<double>{0.0, 3.5e5});
        result.time_tables.emplace_back("normal_x", times, std::vector<double>{0.0, 0.0});
        result.time_tables.emplace_back("tangential_y", times, std::vector<double>{0.0, 2.0e-2});
        result.time_tables.emplace_back("tangential_z", times, std::vector<double>{0.0, 1.0e-2});
        result.time_tables.emplace_back(
            "bending_traction", times, std::vector<double>{0.0, parameters.bending_traction});
    } else if (parameters.path == ScanPath::contact_cycle) {
        result.time_tables.emplace_back(
            "secondary_temperature", times, std::vector<double>{300.0, 301.0, 301.0, 301.0, 301.0});
        result.time_tables.emplace_back("pressure", times, std::vector<double>(times.size(), 0.0));
        result.time_tables.emplace_back("normal_x", times, std::vector<double>{0.0, -2.0e-4, 0.0, -2.0e-4, -2.0e-4});
        result.time_tables.emplace_back("tangential_y", times, std::vector<double>(times.size(), 0.0));
        result.time_tables.emplace_back("tangential_z", times, std::vector<double>(times.size(), 0.0));
        result.time_tables.emplace_back("bending_traction", times, std::vector<double>(times.size(), 0.0));
    } else if (parameters.path == ScanPath::nonmatching_contact_cycle) {
        result.time_tables.emplace_back("secondary_temperature", times,
            std::vector<double>{
                300.0, 308.0, 310.0, 310.0, 310.0, 310.0, 310.0, 310.0, 310.0, 310.0, 310.0, 310.0, 310.0});
        result.time_tables.emplace_back("pressure", times, std::vector<double>(times.size(), 0.0));
        result.time_tables.emplace_back("normal_x", times,
            std::vector<double>{0.0, 0.0, -5.0e-4, -5.0e-4, 0.0, 0.0, 0.0, -5.0e-4, -5.0e-4, 0.0, 0.0, 0.0, -5.0e-4});
        result.time_tables.emplace_back("tangential_y", times,
            std::vector<double>{0.0, 0.04, 0.05, 0.05, 0.05, 0.95, 0.95, 0.95, 0.95, 0.95, 0.05, 0.05, 0.05});
        result.time_tables.emplace_back("tangential_z", times,
            std::vector<double>{0.0, 0.032, 0.04, 0.04, 0.04, 0.20, 0.20, 0.20, 0.20, 0.20, 0.40, 0.40, 0.40});
        result.time_tables.emplace_back("bending_traction", times, std::vector<double>(times.size(), 0.0));
    } else {
        result.time_tables.emplace_back(
            "secondary_temperature", times, std::vector<double>{300.0, 400.0, 400.0, 400.0, 400.0});
        result.time_tables.emplace_back("pressure", times, std::vector<double>(times.size(), 0.0));
        result.time_tables.emplace_back("normal_x", times,
            std::vector<double>{0.0, -1.0e-3 * parameters.normal_displacement_scale,
                -1.0e-3 * parameters.normal_displacement_scale, -1.0e-3 * parameters.normal_displacement_scale,
                -1.0e-3 * parameters.normal_displacement_scale});
        result.time_tables.emplace_back("tangential_y", times,
            std::vector<double>{0.0, 2.0e-5 * parameters.tangential_displacement_scale,
                1.2e-2 * parameters.tangential_displacement_scale, -4.0e-3 * parameters.tangential_displacement_scale,
                -3.98e-3 * parameters.tangential_displacement_scale});
        result.time_tables.emplace_back("tangential_z", times, std::vector<double>(times.size(), 0.0));
        result.time_tables.emplace_back("bending_traction", times, std::vector<double>(times.size(), 0.0));
    }
    result.boundary_conditions = {
        {"primary_temperature", fuelsim::BoundaryConditionType::dirichlet, "primary_outer", fuelsim::Field::temperature,
            300.0},
        {"primary_fix_x", fuelsim::BoundaryConditionType::dirichlet, "primary_outer", fuelsim::Field::displacement_x,
            0.0},
        {"primary_fix_y", fuelsim::BoundaryConditionType::dirichlet, "primary_outer", fuelsim::Field::displacement_y,
            0.0},
        {"primary_fix_z", fuelsim::BoundaryConditionType::dirichlet, "primary_outer", fuelsim::Field::displacement_z,
            0.0},
        {"secondary_temperature", fuelsim::BoundaryConditionType::dirichlet, "secondary_outer",
            fuelsim::Field::temperature, 1.0, false, "secondary_temperature"},
    };
    if (parameters.path != ScanPath::monotonic) {
        const std::string tangential_boundary =
            parameters.path == ScanPath::nonmatching_contact_cycle ? "secondary_all_surface" : "secondary_outer";
        result.boundary_conditions.push_back({"secondary_normal_x", fuelsim::BoundaryConditionType::dirichlet,
            "secondary_outer", fuelsim::Field::displacement_x, 1.0, false, "normal_x"});
        result.boundary_conditions.push_back({"secondary_tangential_y", fuelsim::BoundaryConditionType::dirichlet,
            tangential_boundary, fuelsim::Field::displacement_y, 1.0, false, "tangential_y"});
        result.boundary_conditions.push_back({"secondary_tangential_z", fuelsim::BoundaryConditionType::dirichlet,
            tangential_boundary, fuelsim::Field::displacement_z, 1.0, false, "tangential_z"});
    } else if (parameters.traction_controlled) {
        result.boundary_conditions.push_back(
            configured_boundary("lower_bending_traction", fuelsim::BoundaryConditionType::traction,
                "secondary_outer_lower", fuelsim::Field::displacement_z, "bending_traction", true));
        result.boundary_conditions.push_back(
            configured_boundary("upper_bending_traction", fuelsim::BoundaryConditionType::traction,
                "secondary_outer_upper", fuelsim::Field::displacement_z, "bending_traction", true));
        result.boundary_conditions.back().value = -1.0;
        if (parameters.anchor_bending) {
            result.boundary_conditions.push_back({"secondary_anchor_y", fuelsim::BoundaryConditionType::dirichlet,
                "secondary_y0_surface", fuelsim::Field::displacement_y, 0.0});
            result.boundary_conditions.push_back({"secondary_anchor_z", fuelsim::BoundaryConditionType::dirichlet,
                "secondary_y0_surface", fuelsim::Field::displacement_z, 0.0});
        }
    } else {
        result.boundary_conditions.push_back({"secondary_tangential_y", fuelsim::BoundaryConditionType::dirichlet,
            "secondary_outer", fuelsim::Field::displacement_y, 1.0, false, "tangential_y"});
        result.boundary_conditions.push_back({"secondary_tangential_z", fuelsim::BoundaryConditionType::dirichlet,
            "secondary_outer", fuelsim::Field::displacement_z, 1.0, false, "tangential_z"});
    }
    if (parameters.path == ScanPath::monotonic)
        result.boundary_conditions.push_back(
            configured_boundary("outer_pressure", fuelsim::BoundaryConditionType::pressure, "secondary_outer",
                fuelsim::Field::displacement_x, "pressure", true));
    fuelsim::ContactDefinition contact;
    contact.name = "coupled_contact";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.thermal = true;
    contact.mechanical = true;
    contact.gap_conductivity = 1.0;
    contact.minimum_gap = 1.0;
    contact.penalty = parameters.penalty;
    contact.friction_coefficient = parameters.friction_coefficient;
    contact.friction_slip_tolerance = parameters.slip_tolerance;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    contact.gap_heat_conductance_law = fuelsim::GapHeatConductanceLaw::affine;
    contact.gap_conductance = parameters.conductance;
    contact.gap_conductance_pressure_derivative = parameters.pressure_conductance;
    result.contacts.push_back(contact);
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = 1.0e-8;
    result.relative_tolerance = 1.0e-12;
    result.maximum_iterations = 60;
    result.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    result.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    result.field_residual_scaling = true;
    result.residual_reduction_tolerance = 1.0e-10;
    result.temperature_residual_absolute_tolerance = 1.0e-6;
    result.mechanical_residual_absolute_tolerance = 1.0e-4;
    return result;
}

std::vector<double> raw_residual(fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    std::vector<double> result(problem.dof_count(), 0.0);
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            result[workspace.dofs[local]] += workspace.residual[local];
    }
    return result;
}

std::vector<double> assembled_contact_residual(const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state, fuelsim::SpatialContributionType type, double* conservation_error = nullptr) {
    std::vector<double> result(state.size(), 0.0);
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != type) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, nullptr);
        double balance = 0.0;
        for (std::size_t row = 0; row < dofs.size(); ++row) {
            result[dofs[row]] += residual[row];
            if (type == fuelsim::SpatialContributionType::thermal_contact && row < 8) balance += residual[row];
        }
        if (conservation_error != nullptr) *conservation_error = std::max(*conservation_error, std::abs(balance));
    }
    return result;
}

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance,
    double qualified_pointwise_absolute_tolerance = 0.0) {
    if (metrics.has_relative_norm()) {
        const bool aggregate_passed =
            metrics.relative_l2() < relative_tolerance && metrics.relative_absolute_peak() < relative_tolerance;
        const double maximum_pointwise_absolute_difference =
            std::abs(metrics.maximum_pointwise_relative_actual - metrics.maximum_pointwise_relative_reference);
        const bool pointwise_passed =
            metrics.maximum_pointwise_relative_error() < relative_tolerance ||
            (qualified_pointwise_absolute_tolerance > 0.0 &&
                maximum_pointwise_absolute_difference < qualified_pointwise_absolute_tolerance);
        if (!aggregate_passed || !pointwise_passed) return false;
    }
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

using GroupedFieldErrorMetrics = fuelsim::test::GroupedFieldErrorMetrics;

void print_grouped_metrics(const std::string& name, const GroupedFieldErrorMetrics& metrics) {
    fuelsim::test::print_grouped_relative_metrics(name, metrics);
}

bool grouped_metrics_pass(const GroupedFieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance) {
    if (metrics.has_relative_norm() && !fuelsim::test::grouped_relative_metrics_below(metrics, relative_tolerance))
        return false;
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

struct ScanResponse final {
    bool completed = false, full_field_passed = false, transition_verified = false;
    std::size_t accepted_steps = 0, rejected_steps = 0, active_contact_nodes = 0, open_contact_steps = 0,
                active_contact_steps = 0, sticking_nodes_seen = 0, sliding_nodes_seen = 0;
    int nonlinear_iterations = 0;
    double contact_force = 0.0, tangential_force = 0.0, contact_heat_rate = 0.0, maximum_penetration = 0.0,
           average_contact_temperature = 0.0, maximum_plastic_strain = 0.0, maximum_creep_strain = 0.0,
           maximum_equivalent_stress = 0.0, maximum_absolute_z_displacement = 0.0, friction_dissipation = 0.0,
           external_work = 0.0, boundary_heat_rate = 0.0,
           minimum_tangential_y_resultant = std::numeric_limits<double>::infinity(),
           maximum_tangential_y_resultant = -std::numeric_limits<double>::infinity();
    std::string failure_message;
};

ScanResponse run_scan(std::size_t through_thickness_elements, std::size_t tangential_elements, double step,
    const ScanParameters& parameters, double distortion, const std::string& case_name,
    const std::string& reference_directory) {
    const fuelsim::SpatialDefinition case_definition = definition(parameters);
    const fuelsim::UnstructuredHex8Mesh case_mesh = mesh(through_thickness_elements, tangential_elements, distortion,
        parameters.initial_gap, parameters.nonmatching_mesh);
    fuelsim::TransientProblem problem(case_definition, case_mesh);
    fuelsim::test::AbaqusHex8SnapshotObserver observer;
    const bool transition_case = case_name.rfind("b526_", 0) == 0 || case_name.rfind("b53", 0) == 0 ||
                                 parameters.path == ScanPath::nonmatching_contact_cycle;
    const double end_time = parameters.path == ScanPath::nonmatching_contact_cycle ? 0.9
                            : transition_case ? (parameters.path == ScanPath::contact_cycle ? 0.3 : 0.4)
                                              : 0.2;
    const fuelsim::TransientResult solve =
        fuelsim::solve_transient(problem, {end_time, step, step, step, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
    ScanResponse response;
    response.completed = solve.completed;
    response.accepted_steps = solve.accepted_steps.size();
    response.rejected_steps = solve.rejected_steps.size();
    response.nonlinear_iterations = solve.total_nonlinear_iterations;
    if (!solve.rejected_steps.empty())
        response.failure_message =
            std::string(fuelsim::solve_failure_category_name(solve.rejected_steps.back().failure_category)) + ": " +
            solve.rejected_steps.back().failure_message;
    if (!solve.completed || solve.accepted_steps.empty()) return response;
    fuelsim::test::AbaqusHex8FullFieldOptions comparison;
    comparison.case_name = case_name;
    comparison.reference_prefix = reference_directory + "/" + case_name;
    comparison.expected_steps = static_cast<std::size_t>(std::llround(end_time / step));
    comparison.time_step = step;
    comparison.reduced_integration = parameters.reduced_integration;
    comparison.use_contact_summary_total_slip = true;
    if (case_name.rfind("b524_", 0) == 0) {
        comparison.reaction_heat_flux_pointwise_absolute_tolerance = 5.0e-2;
        comparison.inelastic_pointwise_relative_tolerance = 4.0e-2;
        comparison.energy_pointwise_relative_tolerance = 3.0e-2;
        comparison.contact_replayed_heat_rate_relative_tolerance = 4.0e-2;
        comparison.contact_replayed_heat_rate_pointwise_relative_tolerance = 4.0e-2;
        comparison.contact_replayed_heat_rate_pointwise_absolute_tolerance = 1.0;
        comparison.contact_total_heat_rate_relative_tolerance = 4.0e-2;
        comparison.gate_contact_slip = false;
        comparison.gate_contact_pressure = false;
        comparison.gate_contact_state = false;
    }
    if (case_name.rfind("b525_", 0) == 0) {
        comparison.reaction_heat_flux_pointwise_absolute_tolerance = 1.0;
        comparison.contact_pointwise_relative_tolerance = 3.0e-2;
        comparison.contact_slip_pointwise_relative_tolerance = 2.5e-1;
        comparison.contact_replayed_heat_rate_relative_tolerance = 3.0e-2;
        comparison.contact_replayed_heat_rate_pointwise_absolute_tolerance = 1.0;
        comparison.tangent_basis_tolerance = 1.0e-5;
        if (case_name == "b525_poisson_0499_refined") {
            comparison.contact_relative_tolerance = 1.0e-2;
            comparison.contact_pointwise_relative_tolerance = 1.0e-2;
        }
    }
    if (case_name.rfind("b526_", 0) == 0 || case_name.rfind("b53", 0) == 0 || case_name.rfind("b540_", 0) == 0) {
        comparison.bulk_relative_tolerance = 1.0e-2;
        comparison.contact_relative_tolerance = 5.0e-3;
        comparison.contact_pointwise_relative_tolerance = 1.25e-2;
        comparison.energy_relative_tolerance = 1.0e-2;
        comparison.reaction_heat_flux_pointwise_absolute_tolerance = 2.0e-1;
        comparison.contact_replayed_heat_rate_relative_tolerance = 4.0e-2;
        comparison.contact_replayed_heat_rate_pointwise_relative_tolerance = 4.0e-2;
        comparison.contact_replayed_heat_rate_pointwise_absolute_tolerance = 1.0;
        comparison.contact_total_heat_rate_relative_tolerance = 4.0e-2;
        comparison.contact_slip_pointwise_absolute_tolerance = 5.0e-6;
        if (case_name == "b526_contact_cycle" || case_name == "b538_contact_cycle") {
            comparison.displacement_pointwise_absolute_tolerance = 1.0e-12;
            comparison.logarithmic_strain_pointwise_absolute_tolerance = 1.0e-11;
            comparison.external_work_pointwise_absolute_tolerance = 1.0e-12;
            if (case_name == "b526_contact_cycle") comparison.gate_contact_slip = false;
        }
    }
    if (case_name == "b540_nonmatching_contact_cycle") {
        comparison.displacement_pointwise_absolute_tolerance = 5.0e-8;
        comparison.reaction_pointwise_absolute_tolerance = 5.0e-2;
        comparison.stress_pointwise_absolute_tolerance = 1.0e-1;
        comparison.logarithmic_strain_pointwise_absolute_tolerance = 2.0e-8;
        comparison.elastic_strain_pointwise_absolute_tolerance = 2.0e-10;
        comparison.energy_pointwise_relative_tolerance = 7.5e-2;
        comparison.hourglass_energy_pointwise_absolute_tolerance = 1.0e-8;
        comparison.external_work_pointwise_absolute_tolerance = 1.0e-12;
        comparison.contact_total_heat_rate_pointwise_relative_tolerance = 2.0e-1;
        comparison.gate_contact_state = false;
    }
    response.full_field_passed = fuelsim::test::compare_abaqus_hex8_full_field(
        problem, case_definition, case_mesh, observer.snapshots(), comparison);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.committed_state);
    response.active_contact_nodes = interface.active_contact_nodes;
    response.contact_force = interface.total_contact_force;
    response.tangential_force = interface.total_tangential_force;
    response.contact_heat_rate = interface.total_heat_rate;
    response.maximum_penetration = std::max(0.0, -interface.minimum_contact_gap);
    response.boundary_heat_rate = solve.accepted_steps.back().conservation.dirichlet_heat_input_rate;
    bool active_seen = false, reopened_seen = false, recontact_seen = false, crossed_primary_face = false;
    bool final_active = false, final_sticking = false;
    std::size_t first_active_primary_face = std::numeric_limits<std::size_t>::max();
    for (const fuelsim::test::AbaqusHex8StepSnapshot& snapshot : observer.snapshots()) {
        std::size_t active_nodes = 0, sticking_nodes = 0, sliding_nodes = 0;
        double tangential_y_resultant = 0.0;
        for (const fuelsim::CartesianContactNodeSummary& point : snapshot.contact) {
            if (!(point.pressure > 0.0)) continue;
            if (first_active_primary_face == std::numeric_limits<std::size_t>::max())
                first_active_primary_face = point.primary_face;
            else if (point.primary_face != first_active_primary_face)
                crossed_primary_face = true;
            ++active_nodes;
            ++(point.sliding ? sliding_nodes : sticking_nodes);
            tangential_y_resultant += point.tangential_contact_force[1];
        }
        response.sticking_nodes_seen += sticking_nodes;
        response.sliding_nodes_seen += sliding_nodes;
        if (active_nodes == 0) {
            ++response.open_contact_steps;
            if (active_seen) reopened_seen = true;
        } else {
            ++response.active_contact_steps;
            if (reopened_seen) recontact_seen = true;
            active_seen = true;
        }
        response.minimum_tangential_y_resultant =
            std::min(response.minimum_tangential_y_resultant, tangential_y_resultant);
        response.maximum_tangential_y_resultant =
            std::max(response.maximum_tangential_y_resultant, tangential_y_resultant);
        final_active = active_nodes > 0;
        final_sticking = final_active && sliding_nodes == 0;
    }
    if (parameters.path == ScanPath::contact_cycle)
        response.transition_verified = active_seen && reopened_seen && recontact_seen && final_active;
    else if (parameters.path == ScanPath::nonmatching_contact_cycle)
        response.transition_verified =
            active_seen && crossed_primary_face && reopened_seen && recontact_seen && final_active;
    else if (parameters.path == ScanPath::friction_reversal)
        response.transition_verified = response.sticking_nodes_seen > 0 && response.sliding_nodes_seen > 0 &&
                                       response.minimum_tangential_y_resultant < 0.0 &&
                                       response.maximum_tangential_y_resultant > 0.0 && final_sticking;
    else
        response.transition_verified = true;
    for (const fuelsim::TransientAcceptedStep& accepted : solve.accepted_steps) {
        response.friction_dissipation += accepted.conservation.friction_dissipation_increment;
        response.external_work += accepted.conservation.trapezoidal_pressure_traction_work_increment +
                                  accepted.conservation.trapezoidal_dirichlet_reaction_work_increment;
    }
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::map<std::size_t, std::size_t> source_to_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
            source_to_global.emplace(
                spatial.region_mesh(region).source_node_ids()[local], spatial.global_node(region, local));
    const std::vector<std::size_t> contact_sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    for (const std::size_t source : contact_sources)
        response.average_contact_temperature +=
            solve.committed_state[spatial.dof(fuelsim::Field::temperature, source_to_global.at(source))];
    response.average_contact_temperature /= static_cast<double>(contact_sources.size());
    for (std::size_t node = 0; node < spatial.node_count(); ++node)
        response.maximum_absolute_z_displacement = std::max(response.maximum_absolute_z_displacement,
            std::abs(solve.committed_state[spatial.dof(fuelsim::Field::displacement_z, node)]));
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t element = 0; element < spatial.region_mesh(region).elements().size(); ++element)
            for (const fuelsim::CartesianMaterialPointState& point :
                fuelsim::cartesian::ProblemAccess::material_history(problem, region, element)) {
                response.maximum_plastic_strain =
                    std::max(response.maximum_plastic_strain, point.equivalent_plastic_strain);
                response.maximum_creep_strain = std::max(response.maximum_creep_strain, point.equivalent_creep_strain);
                const double mean = (point.stress.xx + point.stress.yy + point.stress.zz) / 3.0;
                const double deviator_x = point.stress.xx - mean, deviator_y = point.stress.yy - mean,
                             deviator_z = point.stress.zz - mean;
                response.maximum_equivalent_stress = std::max(response.maximum_equivalent_stress,
                    std::sqrt(1.5 * (deviator_x * deviator_x + deviator_y * deviator_y + deviator_z * deviator_z +
                                        2.0 * (point.stress.xy * point.stress.xy + point.stress.yz * point.stress.yz +
                                                  point.stress.xz * point.stress.xz))));
            }
    return response;
}

void print_scan(const std::string& name, const ScanResponse& response) {
    const std::string prefix =
        name.rfind("b5", 0) == 0 ? name : (name.rfind("poisson_", 0) == 0 ? "b525_" + name : "b524_" + name);
    std::cout << prefix << "_completed=" << response.completed << '\n'
              << prefix << "_accepted_steps=" << response.accepted_steps << '\n'
              << prefix << "_rejected_steps=" << response.rejected_steps << '\n'
              << prefix << "_full_field_passed=" << response.full_field_passed << '\n'
              << prefix << "_transition_verified=" << response.transition_verified << '\n'
              << prefix << "_nonlinear_iterations=" << response.nonlinear_iterations << '\n'
              << prefix << "_active_contact_nodes=" << response.active_contact_nodes << '\n'
              << prefix << "_open_contact_steps=" << response.open_contact_steps << '\n'
              << prefix << "_active_contact_steps=" << response.active_contact_steps << '\n'
              << prefix << "_sticking_nodes_seen=" << response.sticking_nodes_seen << '\n'
              << prefix << "_sliding_nodes_seen=" << response.sliding_nodes_seen << '\n'
              << prefix << "_minimum_tangential_y_resultant=" << response.minimum_tangential_y_resultant << '\n'
              << prefix << "_maximum_tangential_y_resultant=" << response.maximum_tangential_y_resultant << '\n'
              << prefix << "_contact_force=" << response.contact_force << '\n'
              << prefix << "_tangential_force=" << response.tangential_force << '\n'
              << prefix << "_contact_heat_rate=" << response.contact_heat_rate << '\n'
              << prefix << "_maximum_penetration=" << response.maximum_penetration << '\n'
              << prefix << "_average_contact_temperature=" << response.average_contact_temperature << '\n'
              << prefix << "_maximum_plastic_strain=" << response.maximum_plastic_strain << '\n'
              << prefix << "_maximum_creep_strain=" << response.maximum_creep_strain << '\n'
              << prefix << "_maximum_equivalent_stress=" << response.maximum_equivalent_stress << '\n'
              << prefix << "_maximum_absolute_z_displacement=" << response.maximum_absolute_z_displacement << '\n'
              << prefix << "_friction_dissipation=" << response.friction_dissipation << '\n'
              << prefix << "_external_work=" << response.external_work << '\n'
              << prefix << "_boundary_heat_rate=" << response.boundary_heat_rate << '\n';
    if (!response.failure_message.empty())
        std::cout << prefix << "_failure_message=" << response.failure_message << '\n';
}

void write_scan_response(
    const std::filesystem::path& path, const std::string& case_name, const ScanResponse& response) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not open B5.24/B5.25 response output: " + path.string());
    output << std::scientific << std::setprecision(17) << case_name << '\n'
           << response.completed << ' ' << response.full_field_passed << ' ' << response.transition_verified << ' '
           << response.accepted_steps << ' ' << response.rejected_steps << ' ' << response.active_contact_nodes << ' '
           << response.open_contact_steps << ' ' << response.active_contact_steps << ' ' << response.sticking_nodes_seen
           << ' ' << response.sliding_nodes_seen << ' ' << response.nonlinear_iterations << '\n'
           << response.contact_force << ' ' << response.tangential_force << ' ' << response.contact_heat_rate << ' '
           << response.maximum_penetration << ' ' << response.average_contact_temperature << ' '
           << response.maximum_plastic_strain << ' ' << response.maximum_creep_strain << ' '
           << response.maximum_equivalent_stress << ' ' << response.maximum_absolute_z_displacement << ' '
           << response.friction_dissipation << ' ' << response.external_work << ' ' << response.boundary_heat_rate
           << ' ' << response.minimum_tangential_y_resultant << ' ' << response.maximum_tangential_y_resultant << '\n';
    if (!output) throw std::runtime_error("Could not write B5.24/B5.25 response output: " + path.string());
}

ScanResponse read_scan_response(const std::filesystem::path& directory, const std::string& case_name) {
    const std::filesystem::path path = directory / (case_name + ".txt");
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not open B5.24/B5.25 response input: " + path.string());
    std::string stored_case;
    ScanResponse response;
    input >> stored_case >> response.completed >> response.full_field_passed >> response.transition_verified >>
        response.accepted_steps >> response.rejected_steps >> response.active_contact_nodes >>
        response.open_contact_steps >> response.active_contact_steps >> response.sticking_nodes_seen >>
        response.sliding_nodes_seen >> response.nonlinear_iterations >> response.contact_force >>
        response.tangential_force >> response.contact_heat_rate >> response.maximum_penetration >>
        response.average_contact_temperature >> response.maximum_plastic_strain >> response.maximum_creep_strain >>
        response.maximum_equivalent_stress >> response.maximum_absolute_z_displacement >>
        response.friction_dissipation >> response.external_work >> response.boundary_heat_rate >>
        response.minimum_tangential_y_resultant >> response.maximum_tangential_y_resultant;
    if (!input || stored_case != case_name)
        throw std::runtime_error("Invalid B5.24/B5.25 response input: " + path.string());
    return response;
}

ScanResponse run_named_scan(const std::string& name, const std::string& reference_directory) {
    ScanParameters parameters;
    std::size_t through_thickness_elements = 2, tangential_elements = 1;
    double step = 0.02, distortion = 0.0;
    if (name == "b524_mesh_medium")
        through_thickness_elements = 3;
    else if (name == "b524_mesh_fine")
        through_thickness_elements = 4;
    else if (name == "b524_time_coarse")
        step = 0.04;
    else if (name == "b524_time_fine")
        step = 0.01;
    else if (name == "b524_penalty_low")
        parameters.penalty = 5.0e8;
    else if (name == "b524_penalty_high")
        parameters.penalty = 2.0e9;
    else if (name == "b524_friction_low")
        parameters.friction_coefficient = 0.01;
    else if (name == "b524_friction_high")
        parameters.friction_coefficient = 0.10;
    else if (name == "b524_slip_low")
        parameters.slip_tolerance = 0.0025;
    else if (name == "b524_slip_high")
        parameters.slip_tolerance = 0.0100;
    else if (name == "b524_thermal_low")
        parameters.pressure_conductance = 0.0005;
    else if (name == "b524_thermal_high")
        parameters.pressure_conductance = 0.0020;
    else if (name == "b526_contact_cycle") {
        parameters.elastic_only = true;
        parameters.initial_gap = 1.0e-4;
        parameters.path = ScanPath::contact_cycle;
    } else if (name == "b526_friction_reversal") {
        parameters.elastic_only = true;
        parameters.path = ScanPath::friction_reversal;
    } else if (name == "b538_contact_cycle") {
        parameters.elastic_only = true;
        parameters.initial_gap = 1.0e-4;
        parameters.path = ScanPath::contact_cycle;
        parameters.reduced_integration = true;
    } else if (name == "b539_friction_reversal") {
        parameters.elastic_only = true;
        parameters.path = ScanPath::friction_reversal;
        parameters.reduced_integration = true;
        parameters.normal_displacement_scale = 0.125;
    } else if (name == "b540_nonmatching_contact_cycle") {
        through_thickness_elements = 1;
        parameters.elastic_only = true;
        parameters.initial_gap = 5.0e-4;
        parameters.path = ScanPath::nonmatching_contact_cycle;
        parameters.reduced_integration = true;
        parameters.nonmatching_mesh = true;
    } else if (name != "b524_mesh_coarse" && name.rfind("b525_", 0) != 0)
        throw std::invalid_argument("Unknown B5.24, B5.25, B5.26, B5.38, B5.39, or B5.40 selected case: " + name);
    if (name.rfind("b525_", 0) == 0) {
        through_thickness_elements = name == "b525_poisson_0499_refined" ? 2 : 1;
        tangential_elements = 16;
        distortion = 0.08;
        parameters.traction_controlled = true;
        parameters.elastic_only = true;
        parameters.anchor_bending = true;
        parameters.bending_traction = 2.0e4;
        parameters.friction_coefficient = 0.2;
        if (name == "b525_poisson_030")
            parameters.primary_poisson = parameters.secondary_poisson = 0.30;
        else if (name == "b525_poisson_045")
            parameters.primary_poisson = parameters.secondary_poisson = 0.45;
        else if (name == "b525_poisson_049")
            parameters.primary_poisson = parameters.secondary_poisson = 0.49;
        else if (name == "b525_poisson_0499" || name == "b525_poisson_0499_refined")
            parameters.primary_poisson = parameters.secondary_poisson = 0.499;
        else
            throw std::invalid_argument("Unknown B5.25 selected case: " + name);
    }
    return run_scan(
        through_thickness_elements, tangential_elements, step, parameters, distortion, name, reference_directory);
}

bool validate_aggregate_responses(const std::filesystem::path& directory) {
    const ScanResponse scan_base = read_scan_response(directory, "b524_mesh_coarse");
    const ScanResponse mesh_medium = read_scan_response(directory, "b524_mesh_medium");
    const ScanResponse mesh_fine = read_scan_response(directory, "b524_mesh_fine");
    const ScanResponse time_coarse = read_scan_response(directory, "b524_time_coarse");
    const ScanResponse time_fine = read_scan_response(directory, "b524_time_fine");
    const ScanResponse penalty_low = read_scan_response(directory, "b524_penalty_low");
    const ScanResponse penalty_high = read_scan_response(directory, "b524_penalty_high");
    const ScanResponse friction_low = read_scan_response(directory, "b524_friction_low");
    const ScanResponse friction_high = read_scan_response(directory, "b524_friction_high");
    const ScanResponse slip_low = read_scan_response(directory, "b524_slip_low");
    const ScanResponse slip_high = read_scan_response(directory, "b524_slip_high");
    const ScanResponse thermal_low = read_scan_response(directory, "b524_thermal_low");
    const ScanResponse thermal_high = read_scan_response(directory, "b524_thermal_high");
    const std::array<std::pair<const char*, const ScanResponse*>, 13> scans = {{{"b524_mesh_coarse", &scan_base},
        {"b524_mesh_medium", &mesh_medium}, {"b524_mesh_fine", &mesh_fine}, {"b524_time_coarse", &time_coarse},
        {"b524_time_fine", &time_fine}, {"b524_penalty_low", &penalty_low}, {"b524_penalty_high", &penalty_high},
        {"b524_friction_low", &friction_low}, {"b524_friction_high", &friction_high}, {"b524_slip_low", &slip_low},
        {"b524_slip_high", &slip_high}, {"b524_thermal_low", &thermal_low}, {"b524_thermal_high", &thermal_high}}};
    bool passed = true, scans_completed = true;
    for (const auto& scan : scans) {
        print_scan(scan.first, *scan.second);
        scans_completed = scans_completed && scan.second->completed && scan.second->rejected_steps == 0 &&
                          scan.second->active_contact_nodes > 0 && scan.second->full_field_passed;
    }
    passed = check(scans_completed, "B5.24 split pressure-controlled scans all converge, retain active contact, and "
                                    "pass their independent Abaqus full-field comparisons") &&
             passed;
    const auto contraction = [](double coarse, double medium, double fine) {
        const double first = std::abs(medium - coarse), second = std::abs(fine - medium);
        return first == 0.0 ? (second == 0.0 ? 0.0 : std::numeric_limits<double>::infinity()) : second / first;
    };
    const double mesh_force_contraction =
                     contraction(scan_base.tangential_force, mesh_medium.tangential_force, mesh_fine.tangential_force),
                 mesh_heat_contraction = contraction(
                     scan_base.contact_heat_rate, mesh_medium.contact_heat_rate, mesh_fine.contact_heat_rate),
                 mesh_temperature_contraction = contraction(scan_base.average_contact_temperature,
                     mesh_medium.average_contact_temperature, mesh_fine.average_contact_temperature),
                 mesh_plastic_contraction = contraction(scan_base.maximum_plastic_strain,
                     mesh_medium.maximum_plastic_strain, mesh_fine.maximum_plastic_strain),
                 mesh_work_contraction =
                     contraction(scan_base.external_work, mesh_medium.external_work, mesh_fine.external_work),
                 mesh_work_relative_change =
                     std::abs(mesh_fine.external_work - scan_base.external_work) / scan_base.external_work,
                 time_temperature_contraction = contraction(time_coarse.average_contact_temperature,
                     scan_base.average_contact_temperature, time_fine.average_contact_temperature),
                 time_creep_contraction = contraction(
                     time_coarse.maximum_creep_strain, scan_base.maximum_creep_strain, time_fine.maximum_creep_strain),
                 time_tangential_force_contraction =
                     contraction(time_coarse.tangential_force, scan_base.tangential_force, time_fine.tangential_force),
                 time_boundary_heat_contraction = contraction(
                     time_coarse.boundary_heat_rate, scan_base.boundary_heat_rate, time_fine.boundary_heat_rate);
    std::cout << "b524_mesh_tangential_force_contraction=" << mesh_force_contraction << '\n'
              << "b524_mesh_contact_heat_rate_contraction=" << mesh_heat_contraction << '\n'
              << "b524_mesh_contact_temperature_contraction=" << mesh_temperature_contraction << '\n'
              << "b524_mesh_maximum_plastic_strain_contraction=" << mesh_plastic_contraction << '\n'
              << "b524_mesh_external_work_contraction=" << mesh_work_contraction << '\n'
              << "b524_mesh_external_work_coarse_fine_relative_change=" << mesh_work_relative_change << '\n'
              << "b524_time_contact_temperature_contraction=" << time_temperature_contraction << '\n'
              << "b524_time_maximum_creep_strain_contraction=" << time_creep_contraction << '\n'
              << "b524_time_tangential_force_contraction=" << time_tangential_force_contraction << '\n'
              << "b524_time_boundary_heat_rate_contraction=" << time_boundary_heat_contraction << '\n';
    passed = check(mesh_force_contraction < 0.8 && mesh_heat_contraction < 0.8 && mesh_temperature_contraction < 0.8 &&
                       mesh_plastic_contraction < 0.8 && mesh_work_relative_change < 5.0e-2,
                 "B5.24 split mesh scans retain contraction of four tracked responses and keep the shortened-path "
                 "external-work change below five percent") &&
             passed;
    passed = check(time_temperature_contraction < 0.75 && time_creep_contraction < 0.75 &&
                       time_tangential_force_contraction < 0.85 && time_boundary_heat_contraction < 0.75,
                 "B5.24 split time-step scans retain contraction of the four tracked responses") &&
             passed;
    passed =
        check(penalty_low.maximum_penetration > scan_base.maximum_penetration &&
                  scan_base.maximum_penetration > penalty_high.maximum_penetration &&
                  std::abs(penalty_high.contact_force - penalty_low.contact_force) / scan_base.contact_force < 1.0e-4,
            "B5.24 split penalty scans reduce penetration while preserving the pressure-controlled resultant") &&
        passed;
    passed = check(friction_low.tangential_force < scan_base.tangential_force &&
                       scan_base.tangential_force < friction_high.tangential_force &&
                       friction_low.friction_dissipation > scan_base.friction_dissipation &&
                       scan_base.friction_dissipation >= friction_high.friction_dissipation,
                 "B5.24 split friction scans retain the resistance and dissipation response") &&
             passed;
    passed = check(slip_low.friction_dissipation > scan_base.friction_dissipation &&
                       scan_base.friction_dissipation > slip_high.friction_dissipation &&
                       std::abs(slip_high.contact_force - slip_low.contact_force) / scan_base.contact_force < 1.0e-3,
                 "B5.24 split slip-tolerance scans retain regularized dissipation and the normal resultant") &&
             passed;
    passed =
        check(thermal_low.contact_heat_rate < scan_base.contact_heat_rate &&
                  scan_base.contact_heat_rate < thermal_high.contact_heat_rate &&
                  thermal_low.average_contact_temperature > scan_base.average_contact_temperature &&
                  scan_base.average_contact_temperature > thermal_high.average_contact_temperature &&
                  std::abs(thermal_high.contact_force - thermal_low.contact_force) / scan_base.contact_force < 1.0e-4,
            "B5.24 split thermal-contact scans retain the heat-transfer and temperature response") &&
        passed;

    const ScanResponse poisson_030 = read_scan_response(directory, "b525_poisson_030");
    const ScanResponse poisson_045 = read_scan_response(directory, "b525_poisson_045");
    const ScanResponse poisson_049 = read_scan_response(directory, "b525_poisson_049");
    const ScanResponse poisson_0499 = read_scan_response(directory, "b525_poisson_0499");
    const ScanResponse poisson_0499_refined = read_scan_response(directory, "b525_poisson_0499_refined");
    const std::array<std::pair<const char*, const ScanResponse*>, 5> poisson_scans = {
        {{"b525_poisson_030", &poisson_030}, {"b525_poisson_045", &poisson_045}, {"b525_poisson_049", &poisson_049},
            {"b525_poisson_0499", &poisson_0499}, {"b525_poisson_0499_refined", &poisson_0499_refined}}};
    bool poisson_scans_completed = true;
    for (const auto& scan : poisson_scans) {
        print_scan(scan.first, *scan.second);
        poisson_scans_completed = poisson_scans_completed && scan.second->completed &&
                                  scan.second->rejected_steps == 0 && scan.second->active_contact_nodes > 0 &&
                                  scan.second->full_field_passed && scan.second->maximum_equivalent_stress > 0.0 &&
                                  scan.second->maximum_absolute_z_displacement > 0.0;
    }
    passed = check(poisson_scans_completed,
                 "B5.25 split distorted bending scans remain finite and pass their Abaqus full-field comparisons") &&
             passed;
    const double near_incompressible_displacement_change = std::abs(poisson_0499.maximum_absolute_z_displacement -
                                                                    poisson_049.maximum_absolute_z_displacement) /
                                                           poisson_049.maximum_absolute_z_displacement,
                 near_incompressible_stress_change =
                     std::abs(poisson_0499.maximum_equivalent_stress - poisson_049.maximum_equivalent_stress) /
                     poisson_049.maximum_equivalent_stress,
                 near_incompressible_reaction_change =
                     std::abs(poisson_0499.contact_force - poisson_049.contact_force) / poisson_049.contact_force,
                 refined_displacement_change = std::abs(poisson_0499_refined.maximum_absolute_z_displacement -
                                                        poisson_0499.maximum_absolute_z_displacement) /
                                               poisson_0499.maximum_absolute_z_displacement,
                 refined_stress_change =
                     std::abs(poisson_0499_refined.maximum_equivalent_stress - poisson_0499.maximum_equivalent_stress) /
                     poisson_0499.maximum_equivalent_stress,
                 refined_reaction_change = std::abs(poisson_0499_refined.contact_force - poisson_0499.contact_force) /
                                           poisson_0499.contact_force;
    std::cout << "b525_poisson_049_to_0499_displacement_relative_change=" << near_incompressible_displacement_change
              << '\n'
              << "b525_poisson_049_to_0499_stress_relative_change=" << near_incompressible_stress_change << '\n'
              << "b525_poisson_049_to_0499_reaction_relative_change=" << near_incompressible_reaction_change << '\n'
              << "b525_poisson_0499_refined_displacement_relative_change=" << refined_displacement_change << '\n'
              << "b525_poisson_0499_refined_stress_relative_change=" << refined_stress_change << '\n'
              << "b525_poisson_0499_refined_reaction_relative_change=" << refined_reaction_change << '\n';
    passed = check(poisson_0499.maximum_absolute_z_displacement > poisson_049.maximum_absolute_z_displacement &&
                       near_incompressible_displacement_change < 1.0e-2 && near_incompressible_stress_change < 1.0e-2 &&
                       near_incompressible_reaction_change < 1.0e-2,
                 "B5.25 split scans retain the smooth Poisson-ratio 0.49 to 0.499 limit") &&
             passed;
    passed = check(refined_displacement_change < 2.0e-1 && refined_stress_change < 1.0e-1 &&
                       refined_reaction_change < 1.0e-2,
                 "B5.25 split shortened-path scans retain the declared thickness-refinement limits") &&
             passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 6 || argc > 8) {
        std::cerr << "Usage: fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests "
                     "<nodal.csv> <integration.csv> <contact.csv> <energy.csv> <Abaqus reference directory> "
                     "[selected B5.24, B5.25, B5.26, B5.38, B5.39, or B5.40 case [response output] | "
                     "--aggregate-responses <response directory>]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        if (argc == 8 && std::string(argv[6]) == "--aggregate-responses") {
            const bool passed = validate_aggregate_responses(argv[7]);
            if (passed) std::cout << "[PASS] B5.24/B5.25 split response aggregation\n";
            return passed ? 0 : 1;
        }
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.23 integrated Abaqus comparison\n");
        if (argc >= 7) {
            if (std::string(argv[6]) == "--aggregate-responses")
                throw std::invalid_argument("The response aggregation mode requires a response directory");
            if (argc == 8) std::filesystem::remove(argv[7]);
            const ScanResponse response = run_named_scan(argv[6], argv[5]);
            print_scan(argv[6], response);
            if (argc == 8) write_scan_response(argv[7], argv[6], response);
            return response.completed && response.rejected_steps == 0 && response.active_contact_nodes > 0 &&
                           response.full_field_passed && response.transition_verified
                       ? 0
                       : 1;
        }
        const std::vector<NodeReference> node_reference = read_nodes(argv[1]);
        const std::vector<IntegrationReference> integration_reference = read_integration(argv[2]);
        const std::vector<ContactReference> contact_reference = read_contact(argv[3]);
        const std::vector<EnergyReference> energy_reference = read_energy(argv[4]);
        const fuelsim::UnstructuredHex8Mesh input_mesh = mesh(2, 1);
        fuelsim::TransientProblem problem(definition(), input_mesh);
        SnapshotObserver observer;
        const fuelsim::TransientResult solve = fuelsim::solve_transient(
            problem, {0.4, time_step, time_step, time_step, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
        bool passed = check(solve.completed && observer.snapshots().size() == increment_count &&
                                solve.accepted_steps.size() == increment_count && solve.rejected_steps.empty(),
            "B5.23 completes all twenty fixed increments without reducing the time step");
        std::cout << "b523_completed=" << solve.completed << '\n'
                  << "b523_committed_time=" << solve.committed_time << '\n'
                  << "b523_rejected_steps=" << solve.rejected_steps.size() << '\n';
        for (const fuelsim::TransientRejectedStep& rejected : solve.rejected_steps)
            std::cout << "b523_rejected_time=" << rejected.attempted_end_time
                      << " category=" << fuelsim::solve_failure_category_name(rejected.failure_category)
                      << " message=" << rejected.failure_message << '\n';
        if (observer.snapshots().empty()) return 1;
        double maximum_plastic = 0.0, maximum_creep = 0.0;
        for (const StepSnapshot& snapshot : observer.snapshots())
            for (const auto& element : snapshot.material)
                for (const fuelsim::CartesianMaterialPointState& point : element) {
                    maximum_plastic = std::max(maximum_plastic, point.equivalent_plastic_strain);
                    maximum_creep = std::max(maximum_creep, point.equivalent_creep_strain);
                }
        const StepSnapshot& final = observer.snapshots().back();
        const std::size_t active = static_cast<std::size_t>(std::count_if(final.contact.begin(), final.contact.end(),
            [](const fuelsim::CartesianContactNodeSummary& point) { return point.pressure > 0.0; }));
        std::cout << "b523_accepted_steps=" << observer.snapshots().size() << '\n'
                  << "b523_maximum_equivalent_plastic_strain=" << maximum_plastic << '\n'
                  << "b523_maximum_equivalent_creep_strain=" << maximum_creep << '\n'
                  << "b523_final_active_contact_nodes=" << active << '\n'
                  << "b523_final_stored_heat_rate=" << final.conservation.stored_heat_rate << '\n'
                  << "b523_final_boundary_heat_rate=" << final.conservation.dirichlet_heat_input_rate << '\n';
        passed =
            check(maximum_plastic > 0.0 && maximum_creep > 0.0, "B5.23 activates both plasticity and Norton creep") &&
            passed;
        passed =
            check(active == 4, "B5.23 keeps all four secondary contact nodes active at the final increment") && passed;

        std::map<std::size_t, std::size_t> source_to_global;
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        for (std::size_t region = 0; region < spatial.region_count(); ++region)
            for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
                source_to_global.emplace(
                    spatial.region_mesh(region).source_node_ids()[local], spatial.global_node(region, local));
        if (source_to_global.size() != 24)
            throw std::invalid_argument("Fuelsim B5.23 source-node map does not contain 24 independent nodes");

        std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
        GroupedFieldErrorMetrics displacement_vector_metrics, reaction_force_vector_metrics;
        fuelsim::TransientProblem reaction_problem(definition(), input_mesh);
        const auto& reaction_dofs = fuelsim::cartesian::ProblemAccess::dof_map(reaction_problem);
        std::vector<unsigned char> constrained_reaction_dofs(reaction_problem.dof_count(), 0U);
        for (const fuelsim::DirichletCondition& condition : reaction_problem.dirichlet_conditions())
            constrained_reaction_dofs.at(condition.dof) = 1U;
        const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
            fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
        for (std::size_t increment = 1; increment <= increment_count; ++increment) {
            const StepSnapshot& snapshot = observer.snapshots().at(increment - 1);
            reaction_problem.begin_time_step({snapshot.time, snapshot.load_factor, true});
            const std::vector<double> reaction = raw_residual(reaction_problem, snapshot.state);
            for (const NodeReference& reference : node_reference) {
                if (reference.increment != increment) continue;
                if (reference.node > 24 || std::abs(reference.time - snapshot.time) > 1.0e-7)
                    throw std::invalid_argument("Abaqus B5.23 nodal label or time lies outside the path");
                const std::size_t global = source_to_global.at(reference.node - 1);
                std::array<double, 3> actual_displacement{}, reference_displacement{}, actual_reaction_force{},
                    reference_reaction_force{};
                for (std::size_t field = 0; field < fields.size(); ++field) {
                    const std::size_t dof = reaction_dofs.dof(fields[field], global);
                    nodal_metrics[field].add(snapshot.state[dof], reference.fields[field]);
                    const double reaction_value = constrained_reaction_dofs[dof] != 0U ? reaction[dof] : 0.0;
                    nodal_metrics[4 + field].add(reaction_value, reference.fields[4 + field]);
                    if (field > 0) {
                        actual_displacement[field - 1] = snapshot.state[dof];
                        reference_displacement[field - 1] = reference.fields[field];
                        actual_reaction_force[field - 1] = reaction_value;
                        reference_reaction_force[field - 1] = reference.fields[4 + field];
                    }
                }
                displacement_vector_metrics.add(
                    actual_displacement.data(), reference_displacement.data(), actual_displacement.size());
                reaction_force_vector_metrics.add(
                    actual_reaction_force.data(), reference_reaction_force.data(), actual_reaction_force.size());
            }
            reaction_problem.commit_time_step(snapshot.state);
        }
        const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y",
            "displacement_z", "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
        const std::array<double, 8> zero_tolerances = {1.0e-8, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-2, 1.0, 1.0, 1.0};
        for (std::size_t field = 0; field < nodal_metrics.size(); ++field)
            fuelsim::test::print_relative_metrics("b523_" + nodal_names[field], nodal_metrics[field]);
        print_grouped_metrics("b523_displacement_vector", displacement_vector_metrics);
        print_grouped_metrics("b523_reaction_force_vector", reaction_force_vector_metrics);
        passed = check(metrics_pass(nodal_metrics[0], 5.0e-3, zero_tolerances[0]),
                     "B5.23 temperature field metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(nodal_metrics[4], 5.0e-3, zero_tolerances[4], 5.0e-2),
                     "B5.23 constrained thermal-reaction aggregate metrics are below 0.5 percent and the near-zero "
                     "pointwise value satisfies the explicit 0.05 W absolute qualification") &&
                 passed;
        passed = check(grouped_metrics_pass(displacement_vector_metrics, 5.0e-3, 1.0e-10),
                     "B5.23 displacement-vector field metrics are below 0.5 percent") &&
                 passed;
        passed = check(grouped_metrics_pass(reaction_force_vector_metrics, 5.0e-3, 1.0),
                     "B5.23 reaction-force-vector field metrics are below 0.5 percent") &&
                 passed;

        std::array<fuelsim::test::FieldErrorMetrics, 37> integration_metrics;
        GroupedFieldErrorMetrics heat_flux_vector_metrics, stress_tensor_metrics, logarithmic_strain_tensor_metrics,
            elastic_strain_tensor_metrics, plastic_strain_tensor_metrics, creep_strain_tensor_metrics;
        double maximum_coordinate_difference = 0.0;
        constexpr std::array<std::size_t, 8> gauss_to_material_node = {0, 1, 3, 2, 4, 5, 7, 6};
        const fuelsim::IsotropicThermoelasticMaterial primary_constitutive(primary_material()),
            secondary_constitutive(secondary_material());
        for (const IntegrationReference& reference : integration_reference) {
            if (reference.increment > observer.snapshots().size() || reference.element > 4 || reference.point > 8)
                throw std::invalid_argument("Abaqus B5.23 integration index lies outside the path");
            const StepSnapshot& snapshot = observer.snapshots().at(reference.increment - 1);
            if (std::abs(reference.time - snapshot.time) > 1.0e-7)
                throw std::invalid_argument("Abaqus B5.23 integration time lies outside the path");
            const std::size_t region = reference.element <= 2 ? 0 : 1;
            const std::size_t local_element = (reference.element - 1) % 2;
            const fuelsim::Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
            const fuelsim::Hex8Element& mesh_element = region_mesh.elements().at(local_element);
            const fuelsim::Hex8Geometry& geometry = spatial.region_element_geometry(region, local_element);
            fuelsim::Hex8LocalValues local_state{};
            for (std::size_t local = 0; local < 8; ++local) {
                const std::size_t global = spatial.global_node(region, mesh_element.nodes[local]);
                local_state[local] = snapshot.state[spatial.dof(fuelsim::Field::temperature, global)];
                local_state[8 + local] = snapshot.state[spatial.dof(fuelsim::Field::displacement_x, global)];
                local_state[16 + local] = snapshot.state[spatial.dof(fuelsim::Field::displacement_y, global)];
                local_state[24 + local] = snapshot.state[spatial.dof(fuelsim::Field::displacement_z, global)];
            }
            std::size_t closest = 0;
            double closest_squared = std::numeric_limits<double>::max();
            for (std::size_t q = 0; q < 8; ++q) {
                fuelsim::CartesianPoint3 current = geometry.points[q].position;
                for (std::size_t node = 0; node < 8; ++node) {
                    current.x += geometry.points[q].shape[node] * local_state[8 + node];
                    current.y += geometry.points[q].shape[node] * local_state[16 + node];
                    current.z += geometry.points[q].shape[node] * local_state[24 + node];
                }
                const double distance_squared = std::pow(current.x - reference.position.x, 2) +
                                                std::pow(current.y - reference.position.y, 2) +
                                                std::pow(current.z - reference.position.z, 2);
                if (distance_squared < closest_squared) {
                    closest = q;
                    closest_squared = distance_squared;
                }
            }
            maximum_coordinate_difference = std::max(maximum_coordinate_difference, std::sqrt(closest_squared));
            const fuelsim::Hex8QuadraturePoint& point = geometry.points[closest];
            fuelsim::Hex8LocalAdValues passive{};
            for (std::size_t local = 0; local < local_state.size(); ++local) passive[local] = local_state[local];
            const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::finite);
            double current_volume = 0.0;
            for (const fuelsim::Hex8QuadraturePoint& volume_point : geometry.points)
                current_volume += fuelsim::evaluate_cartesian_incremental_kinematics(
                    volume_point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::finite)
                                      .current_weighted_measure.value();
            const double material_temperature = local_state[gauss_to_material_node[closest]];
            const fuelsim::IsotropicThermoelasticMaterial& constitutive =
                region == 0 ? primary_constitutive : secondary_constitutive;
            const double conductivity = constitutive
                                            .conductivity(adlite::Scalar(material_temperature),
                                                {snapshot.time, point.position.x, point.position.y, point.position.z})
                                            .value();
            std::array<double, 3> heat_flux{};
            for (std::size_t node = 0; node < 8; ++node)
                for (std::size_t component = 0; component < 3; ++component)
                    heat_flux[component] -=
                        conductivity * kinematics.current_gradient[node][component].value() * local_state[node];
            const fuelsim::CartesianMaterialPointState& material =
                snapshot.material.at(reference.element - 1).at(closest);
            const std::array<double, 6> actual_stress = components(material.stress);
            const std::array<double, 6> expected_stress = components(reference.stress);
            std::array<double, 6> actual_logarithmic = components(logarithmic_strain(point, local_state));
            const double local_trace = actual_logarithmic[0] + actual_logarithmic[1] + actual_logarithmic[2];
            const double selective_trace = std::log(current_volume / geometry.reference_volume);
            for (std::size_t component = 0; component < 3; ++component)
                actual_logarithmic[component] += (selective_trace - local_trace) / 3.0;
            const std::array<double, 6> expected_logarithmic = components(reference.logarithmic_strain);
            const std::array<double, 6> expected_elastic = components(reference.elastic_strain);
            const std::array<double, 6> expected_plastic = components(reference.plastic_strain);
            const std::array<double, 6> expected_creep = components(reference.creep_strain);
            for (std::size_t component = 0; component < 3; ++component)
                integration_metrics[component].add(heat_flux[component], reference.heat_flux[component]);
            for (std::size_t component = 0; component < 6; ++component) {
                integration_metrics[3 + component].add(actual_stress[component], expected_stress[component]);
                integration_metrics[9 + component].add(actual_logarithmic[component], expected_logarithmic[component]);
                integration_metrics[15 + component].add(
                    material.elastic_strain[component], expected_elastic[component]);
                integration_metrics[21 + component].add(
                    material.plastic_strain[component], expected_plastic[component]);
                integration_metrics[28 + component].add(material.creep_strain[component], expected_creep[component]);
            }
            heat_flux_vector_metrics.add(heat_flux.data(), reference.heat_flux.data(), heat_flux.size());
            stress_tensor_metrics.add(actual_stress.data(), expected_stress.data(), actual_stress.size());
            logarithmic_strain_tensor_metrics.add(
                actual_logarithmic.data(), expected_logarithmic.data(), actual_logarithmic.size());
            elastic_strain_tensor_metrics.add(
                material.elastic_strain.data(), expected_elastic.data(), material.elastic_strain.size());
            plastic_strain_tensor_metrics.add(
                material.plastic_strain.data(), expected_plastic.data(), material.plastic_strain.size());
            creep_strain_tensor_metrics.add(
                material.creep_strain.data(), expected_creep.data(), material.creep_strain.size());
            integration_metrics[27].add(material.equivalent_plastic_strain, reference.equivalent_plastic_strain);
            integration_metrics[34].add(material.equivalent_creep_strain, reference.equivalent_creep_strain);
            integration_metrics[35].add(material_temperature, reference.temperature);
            const double selective_volume = point.weighted_measure / geometry.reference_volume * current_volume;
            integration_metrics[36].add(selective_volume, reference.integration_volume);
        }
        const std::array<std::string, 37> integration_names = {"heat_flux_x", "heat_flux_y", "heat_flux_z", "stress_xx",
            "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "log_strain_xx", "log_strain_yy",
            "log_strain_zz", "log_strain_xy", "log_strain_yz", "log_strain_xz", "elastic_strain_xx",
            "elastic_strain_yy", "elastic_strain_zz", "elastic_strain_xy", "elastic_strain_yz", "elastic_strain_xz",
            "plastic_strain_xx", "plastic_strain_yy", "plastic_strain_zz", "plastic_strain_xy", "plastic_strain_yz",
            "plastic_strain_xz", "equivalent_plastic_strain", "creep_strain_xx", "creep_strain_yy", "creep_strain_zz",
            "creep_strain_xy", "creep_strain_yz", "creep_strain_xz", "equivalent_creep_strain", "material_temperature",
            "integration_volume"};
        for (std::size_t field = 0; field < integration_metrics.size(); ++field)
            fuelsim::test::print_relative_metrics("b523_" + integration_names[field], integration_metrics[field]);
        const std::array<std::pair<const char*, const GroupedFieldErrorMetrics*>, 6> grouped_integration = {
            {{"heat_flux_vector", &heat_flux_vector_metrics}, {"stress_tensor", &stress_tensor_metrics},
                {"logarithmic_strain_tensor", &logarithmic_strain_tensor_metrics},
                {"elastic_strain_tensor", &elastic_strain_tensor_metrics},
                {"plastic_strain_tensor", &plastic_strain_tensor_metrics},
                {"creep_strain_tensor", &creep_strain_tensor_metrics}}};
        for (const auto& field : grouped_integration)
            print_grouped_metrics("b523_" + std::string(field.first), *field.second);
        passed = check(grouped_metrics_pass(stress_tensor_metrics, 5.0e-3, 1.0e-12),
                     "B5.23 stress-tensor field metrics are below 0.5 percent") &&
                 passed;
        passed = check(grouped_metrics_pass(logarithmic_strain_tensor_metrics, 5.0e-3, 1.0e-12),
                     "B5.23 logarithmic-strain-tensor field metrics are below 0.5 percent") &&
                 passed;
        passed = check(grouped_metrics_pass(elastic_strain_tensor_metrics, 5.0e-3, 1.0e-12),
                     "B5.23 elastic-strain-tensor field metrics are below 0.5 percent") &&
                 passed;
        passed = check(grouped_metrics_pass(plastic_strain_tensor_metrics, 5.0e-3, 1.0e-12),
                     "B5.23 plastic-strain-tensor field metrics are below 0.5 percent") &&
                 passed;
        passed = check(grouped_metrics_pass(creep_strain_tensor_metrics, 5.0e-3, 1.0e-12),
                     "B5.23 creep-strain-tensor field metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(integration_metrics[27], 5.0e-3, 1.0e-12),
                     "B5.23 equivalent-plastic-strain metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(integration_metrics[34], 5.0e-3, 1.0e-12),
                     "B5.23 equivalent-creep-strain metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(integration_metrics[35], 5.0e-3, 1.0e-12),
                     "B5.23 material-temperature metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(integration_metrics[36], 5.0e-3, 1.0e-12),
                     "B5.23 integration-volume metrics are below 0.5 percent") &&
                 passed;
        std::cout << "b523_maximum_integration_coordinate_difference=" << maximum_coordinate_difference << '\n';
        passed = check(maximum_coordinate_difference < 1.0e-4,
                     "B5.23 maps all 640 current integration-point coordinates uniquely") &&
                 passed;

        const std::vector<std::size_t> contact_sources =
            fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        if (contact_sources.size() != 4)
            throw std::invalid_argument("Fuelsim B5.23 contact surface does not contain four output nodes");
        std::array<fuelsim::test::FieldErrorMetrics, 10> contact_metrics;
        fuelsim::test::FieldErrorMetrics total_contact_heat_rate_metrics;
        std::vector<std::array<double, 3>> previous_actual_plastic(contact_sources.size()),
            previous_reference_plastic(contact_sources.size());
        double actual_friction_dissipation = 0.0, reconstructed_actual_friction_dissipation = 0.0,
               reconstructed_reference_friction_dissipation = 0.0, maximum_contact_heat_conservation_error = 0.0;
        bool contact_states_match = true;
        std::array<fuelsim::test::FieldErrorMetrics, 9> energy_metrics;
        double cumulative_elastic = 0.0, cumulative_plastic = 0.0, cumulative_creep = 0.0,
               cumulative_external_work = 0.0, cumulative_stored_heat = 0.0, cumulative_reference_heat = 0.0;
        constexpr double friction_coefficient = 0.05;
        const double maximum_elastic_slip = 0.005 * std::sqrt(0.5);
        for (std::size_t increment = 1; increment <= increment_count; ++increment) {
            const StepSnapshot& snapshot = observer.snapshots().at(increment - 1);
            const std::vector<double> thermal_contact = assembled_contact_residual(spatial, snapshot.state,
                fuelsim::SpatialContributionType::thermal_contact, &maximum_contact_heat_conservation_error);
            std::array<double, 3> actual_resultant{}, reference_resultant{}, actual_moment{}, reference_moment{},
                actual_center_sum{}, reference_center_sum{};
            double actual_center_weight = 0.0, reference_center_weight = 0.0, actual_contact_heat_rate = 0.0,
                   reference_contact_heat_rate = 0.0;
            for (std::size_t node = 0; node < contact_sources.size(); ++node) {
                const std::size_t source = contact_sources[node], source_label = source + 1;
                const auto found = std::find_if(
                    contact_reference.begin(), contact_reference.end(), [&](const ContactReference& value) {
                        return value.increment == increment && value.node == source_label;
                    });
                if (found == contact_reference.end())
                    throw std::invalid_argument("Abaqus B5.23 contact-node mapping is incomplete");
                if (std::abs(found->time - snapshot.time) > 1.0e-7)
                    throw std::invalid_argument("Abaqus B5.23 contact time lies outside the path");
                const fuelsim::CartesianContactNodeSummary& actual = snapshot.contact.at(node);
                const std::size_t global = source_to_global.at(source);
                const std::array<fuelsim::Field, 3> displacement_fields = {
                    fuelsim::Field::displacement_x, fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
                std::array<double, 3> actual_position = {
                    input_mesh.nodes().at(source).x, input_mesh.nodes().at(source).y, input_mesh.nodes().at(source).z};
                for (std::size_t component = 0; component < 3; ++component)
                    actual_position[component] += snapshot.state[spatial.dof(displacement_fields[component], global)];
                const std::array<double, 3> reference_position = {
                    found->position.x, found->position.y, found->position.z};
                for (std::size_t component = 0; component < 3; ++component)
                    contact_metrics[0].add(actual_position[component], reference_position[component]);
                contact_metrics[1].add(actual.gap, found->opening);
                contact_metrics[2].add(actual.pressure, found->pressure);

                const auto primary_found = std::find_if(input_mesh.nodes().begin(), input_mesh.nodes().begin() + 12,
                    [&](const fuelsim::CartesianPoint3& point) {
                        return point.x == 1.0 && point.y == input_mesh.nodes().at(source).y &&
                               point.z == input_mesh.nodes().at(source).z;
                    });
                if (primary_found == input_mesh.nodes().begin() + 12)
                    throw std::logic_error("B5.23 matching primary contact node was not found");
                const std::size_t primary_source = static_cast<std::size_t>(primary_found - input_mesh.nodes().begin());
                const std::size_t primary_global = source_to_global.at(primary_source);
                std::array<double, 3> relative{};
                for (std::size_t component = 0; component < 3; ++component) {
                    const double primary_position =
                        (component == 0      ? input_mesh.nodes().at(primary_source).x
                            : component == 1 ? input_mesh.nodes().at(primary_source).y
                                             : input_mesh.nodes().at(primary_source).z) +
                        snapshot.state[spatial.dof(displacement_fields[component], primary_global)];
                    relative[component] = actual_position[component] - primary_position;
                }
                std::array<double, 3> normal{};
                if (actual.contact_force > 0.0)
                    for (std::size_t component = 0; component < 3; ++component)
                        normal[component] = actual.normal_contact_force[component] / actual.contact_force;
                double normal_relative = 0.0;
                for (std::size_t component = 0; component < 3; ++component)
                    normal_relative += relative[component] * normal[component];
                std::array<double, 3> actual_total_slip{},
                    reference_total_slip = {0.0, -found->slip_second, found->slip_first};
                for (std::size_t component = 0; component < 3; ++component)
                    actual_total_slip[component] = relative[component] - normal_relative * normal[component];
                const double actual_slip_norm =
                                 std::hypot(actual_total_slip[0], actual_total_slip[1], actual_total_slip[2]),
                             reference_slip_norm =
                                 std::hypot(reference_total_slip[0], reference_total_slip[1], reference_total_slip[2]);
                contact_metrics[3].add(actual_slip_norm, reference_slip_norm);

                std::array<double, 3> actual_force{}, reference_force{}, actual_elastic{}, reference_elastic{};
                double reference_normal_norm = 0.0;
                for (std::size_t component = 0; component < 3; ++component) {
                    actual_force[component] =
                        -actual.normal_contact_force[component] - actual.tangential_contact_force[component];
                    reference_force[component] = found->normal_force[component] + found->shear_force[component];
                    contact_metrics[4].add(-actual.normal_contact_force[component], found->normal_force[component]);
                    contact_metrics[5].add(-actual.tangential_contact_force[component], found->shear_force[component]);
                    reference_normal_norm += found->normal_force[component] * found->normal_force[component];
                }
                reference_normal_norm = std::sqrt(reference_normal_norm);
                if (actual.contact_force > 0.0)
                    for (std::size_t component = 0; component < 3; ++component)
                        actual_elastic[component] = actual.tangential_contact_force[component] * maximum_elastic_slip /
                                                    (friction_coefficient * actual.contact_force);
                if (reference_normal_norm > 0.0)
                    for (std::size_t component = 0; component < 3; ++component)
                        reference_elastic[component] = -found->shear_force[component] * maximum_elastic_slip /
                                                       (friction_coefficient * reference_normal_norm);
                for (std::size_t component = 0; component < 3; ++component) {
                    const double actual_plastic = actual_total_slip[component] - actual_elastic[component];
                    const double reference_plastic = reference_total_slip[component] - reference_elastic[component];
                    reconstructed_actual_friction_dissipation +=
                        actual.tangential_contact_force[component] *
                        (actual_plastic - previous_actual_plastic[node][component]);
                    reconstructed_reference_friction_dissipation -=
                        found->shear_force[component] *
                        (reference_plastic - previous_reference_plastic[node][component]);
                    previous_actual_plastic[node][component] = actual_plastic;
                    previous_reference_plastic[node][component] = reference_plastic;
                    actual_resultant[component] += actual_force[component];
                    reference_resultant[component] += reference_force[component];
                    actual_center_sum[component] += actual.contact_force * actual_position[component];
                    reference_center_sum[component] += reference_normal_norm * reference_position[component];
                }
                actual_center_weight += actual.contact_force;
                reference_center_weight += reference_normal_norm;
                const auto add_moment = [](std::array<double, 3>& moment, const std::array<double, 3>& point,
                                            const std::array<double, 3>& force) {
                    moment[0] += point[1] * force[2] - point[2] * force[1];
                    moment[1] += point[2] * force[0] - point[0] * force[2];
                    moment[2] += point[0] * force[1] - point[1] * force[0];
                };
                add_moment(actual_moment, actual_position, actual_force);
                add_moment(reference_moment, reference_position, reference_force);
                const std::size_t actual_state = actual.pressure <= 0.0 ? 0 : actual.sliding ? 2 : 1;
                contact_states_match = contact_states_match && actual_state == found->state;
                const std::size_t temperature_dof = spatial.dof(fuelsim::Field::temperature, global);
                contact_metrics[9].add(thermal_contact[temperature_dof], found->heat_flux);
                actual_contact_heat_rate += thermal_contact[temperature_dof];
                reference_contact_heat_rate += found->heat_flux;
            }
            total_contact_heat_rate_metrics.add(actual_contact_heat_rate, reference_contact_heat_rate);
            for (std::size_t component = 0; component < 3; ++component) {
                contact_metrics[6].add(actual_resultant[component], reference_resultant[component]);
                contact_metrics[7].add(actual_moment[component], reference_moment[component]);
                contact_metrics[8].add(actual_center_sum[component] / actual_center_weight,
                    reference_center_sum[component] / reference_center_weight);
            }

            const EnergyReference& energy = energy_reference.at(increment - 1);
            if (energy.increment != increment || std::abs(energy.time - snapshot.time) > 1.0e-7)
                throw std::invalid_argument("Abaqus B5.23 energy time lies outside the path");
            cumulative_elastic += snapshot.conservation.elastic_energy_change;
            cumulative_plastic += snapshot.conservation.plastic_dissipation_increment;
            cumulative_creep += snapshot.conservation.creep_dissipation_increment;
            actual_friction_dissipation += snapshot.conservation.friction_dissipation_increment;
            cumulative_external_work += snapshot.conservation.trapezoidal_pressure_traction_work_increment +
                                        snapshot.conservation.trapezoidal_dirichlet_reaction_work_increment;
            cumulative_stored_heat += snapshot.conservation.stored_heat_rate * time_step;
            cumulative_reference_heat += energy.boundary_heat_rate * time_step;
            energy_metrics[0].add(cumulative_elastic + cumulative_plastic + cumulative_creep, energy.internal);
            energy_metrics[1].add(cumulative_elastic, energy.elastic);
            energy_metrics[2].add(cumulative_plastic, energy.plastic);
            energy_metrics[3].add(cumulative_creep, energy.creep);
            energy_metrics[4].add(actual_friction_dissipation, energy.friction);
            energy_metrics[5].add(reconstructed_reference_friction_dissipation, energy.friction);
            energy_metrics[6].add(cumulative_external_work, energy.external_work);
            energy_metrics[7].add(snapshot.conservation.dirichlet_heat_input_rate, energy.boundary_heat_rate);
            energy_metrics[8].add(cumulative_stored_heat, cumulative_reference_heat);
        }
        std::cout << "b523_reconstructed_actual_friction_dissipation=" << reconstructed_actual_friction_dissipation
                  << '\n';
        const std::array<std::string, 10> contact_names = {"current_coordinate", "opening", "pressure",
            "slip_magnitude", "normal_force_vector", "shear_force_vector", "resultant_force_vector",
            "resultant_moment_vector", "normal_force_center", "contact_heat_flux"};
        for (std::size_t field = 0; field < contact_metrics.size(); ++field)
            fuelsim::test::print_relative_metrics("b523_contact_" + contact_names[field], contact_metrics[field]);
        fuelsim::test::print_relative_metrics("b523_contact_total_heat_rate", total_contact_heat_rate_metrics);
        passed = check(metrics_pass(contact_metrics[6], 5.0e-3, 1.0),
                     "B5.23 contact resultant-force metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(contact_metrics[7], 5.0e-3, 1.0),
                     "B5.23 contact resultant-moment metrics are below 0.5 percent") &&
                 passed;
        passed = check(metrics_pass(contact_metrics[8], 5.0e-3, 1.0e-8),
                     "B5.23 contact normal-force-center metrics are below 0.5 percent") &&
                 passed;
        std::cout << "b523_contact_history_point_count=" << final.contact_history.size() << '\n'
                  << "b523_contact_states_match=" << contact_states_match << '\n'
                  << "b523_maximum_contact_heat_conservation_error=" << maximum_contact_heat_conservation_error << '\n';
        passed = check(maximum_contact_heat_conservation_error < 1.0e-10,
                     "B5.23 thermal contact is discretely conservative at every increment") &&
                 passed;

        const std::array<std::string, 9> energy_names = {"internal_energy", "elastic_energy", "plastic_dissipation",
            "creep_dissipation", "friction_dissipation", "abaqus_reconstructed_friction_dissipation", "external_work",
            "boundary_heat_rate", "stored_heat"};
        for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
            fuelsim::test::print_relative_metrics("b523_" + energy_names[field], energy_metrics[field]);
            if (field != 4 && field != 5)
                passed = check(metrics_pass(energy_metrics[field], 1.0e-2, 1.0e-8),
                             "B5.23 " + energy_names[field] + " metrics are below one percent") &&
                         passed;
        }
        passed = check(actual_friction_dissipation > 0.0 && energy_reference.back().friction > 0.0,
                     "B5.23 activates nonzero friction dissipation in both Fuelsim and Abaqus") &&
                 passed;

        if (passed && session.rank() == 0) std::cout << "[PASS] B5.23 integrated Fuelsim path solve\n";
        return passed ? 0 : 1;

    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.23 integrated comparison raised: " << error.what() << '\n';
        return 1;
    }
}
