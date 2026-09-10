#include "solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/field_error_metrics.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t stage_count = 20;
constexpr double time_step = 1.0e-3;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::array<fuelsim::CartesianPoint3, 8> nodes = {{{0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 0.0, 1.0},
    {1.0, 1.0, 1.0},
    {0.0, 1.0, 1.0}}};

struct Deformation final {
    double stretch, shear, rotation;
};

struct NodeReference final {
    std::size_t stage, node;
    double time;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t stage, point;
    double time;
    fuelsim::CartesianPoint3 position;
    double temperature;
    fuelsim::SymmetricTensor3Values stress, integrated_strain, elastic_strain, logarithmic_strain;
    double integration_volume;
};

struct EnergyReference final {
    std::size_t stage;
    double time, internal_energy, strain_energy, external_work;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete Abaqus B5.14 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.14 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

double zero_noise(double value, double threshold) {
    return std::abs(value) < threshold ? 0.0 : value;
}

fuelsim::SymmetricTensor3Values
tensor(const std::vector<std::string>& values, std::size_t first, const std::string& path, double zero_threshold) {
    fuelsim::SymmetricTensor3Values result = {number(values, first, path),
        number(values, first + 1, path),
        number(values, first + 2, path),
        0.5 * number(values, first + 3, path),
        0.5 * number(values, first + 5, path),
        0.5 * number(values, first + 4, path)};
    for (double* value : {&result.xx, &result.yy, &result.zz, &result.xy, &result.yz, &result.xz})
        *value = zero_noise(*value, zero_threshold);
    return result;
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.14 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.14 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11)
            throw std::invalid_argument("Unexpected Abaqus B5.14 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 3, path);
            if (field != 0)
                fields[field] = zero_noise(fields[field], field >= 5 ? 1.0e-3 : (field == 4 ? 1.0e-10 : 1.0e-16));
        }
        result.push_back({positive_integer(number(values, 0, path), path),
            positive_integer(number(values, 2, path), path),
            number(values, 1, path),
            fields});
    }
    if (result.size() != 8 * stage_count)
        throw std::invalid_argument("Abaqus B5.14 nodal reference has an unexpected row count");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.14 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected = "stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
                                 "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
                                 "e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,"
                                 "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
                                 "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,ivol_m3";
    if (line != expected)
        throw std::invalid_argument("Unexpected Abaqus B5.14 integration header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 33)
            throw std::invalid_argument("Unexpected Abaqus B5.14 integration column count in " + path);
        if (positive_integer(number(values, 2, path), path) != 1)
            throw std::invalid_argument("Abaqus B5.14 reference must contain one element");
        fuelsim::SymmetricTensor3Values stress = {number(values, 8, path),
            number(values, 9, path),
            number(values, 10, path),
            number(values, 11, path),
            number(values, 13, path),
            number(values, 12, path)};
        for (double* component : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *component = zero_noise(*component, 1.0);
        result.push_back({positive_integer(number(values, 0, path), path),
            positive_integer(number(values, 3, path), path),
            number(values, 1, path),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            number(values, 7, path),
            stress,
            tensor(values, 14, path, 1.0e-15),
            tensor(values, 20, path, 1.0e-15),
            tensor(values, 26, path, 1.0e-15),
            number(values, 32, path)});
    }
    if (result.size() != 8 * stage_count)
        throw std::invalid_argument("Abaqus B5.14 integration reference has an unexpected row count");
    return result;
}

std::vector<EnergyReference> read_energy(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.14 energy: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,internal_energy_j,strain_energy_j,external_work_j")
        throw std::invalid_argument("Unexpected Abaqus B5.14 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5)
            throw std::invalid_argument("Unexpected Abaqus B5.14 energy column count in " + path);
        result.push_back({positive_integer(number(values, 0, path), path),
            number(values, 1, path),
            number(values, 2, path),
            number(values, 3, path),
            number(values, 4, path)});
    }
    if (result.size() != stage_count)
        throw std::invalid_argument("Abaqus B5.14 energy reference has an unexpected row count");
    return result;
}

Deformation deformation(std::size_t stage) {
    if (stage < 1 || stage > stage_count)
        throw std::out_of_range("B5.14 stage lies outside 1 through 20");
    if (stage <= 5)
        return {1.0 + 0.006 * static_cast<double>(stage), 0.0, 0.0};
    if (stage <= 10)
        return {1.03, 0.008 * static_cast<double>(stage - 5), 0.0};
    if (stage <= 15)
        return {1.03, 0.04, 6.0 * static_cast<double>(stage - 10) * pi / 180.0};
    return {1.03, 0.04, (30.0 - 6.0 * static_cast<double>(stage - 15)) * pi / 180.0};
}

fuelsim::CartesianPoint3 current_position(const fuelsim::CartesianPoint3& point, const Deformation& value) {
    const double base_x = value.stretch * point.x + value.shear * point.y;
    const double cosine = std::cos(value.rotation), sine = std::sin(value.rotation);
    return {cosine * base_x - sine * point.y, sine * base_x + cosine * point.y, point.z};
}

fuelsim::SymmetricTensor3Values logarithmic_strain(const Deformation& value) {
    const double cosine = std::cos(value.rotation), sine = std::sin(value.rotation);
    const double f00 = cosine * value.stretch;
    const double f01 = cosine * value.shear - sine;
    const double f10 = sine * value.stretch;
    const double f11 = sine * value.shear + cosine;
    const double b00 = f00 * f00 + f01 * f01;
    const double b01 = f00 * f10 + f01 * f11;
    const double b11 = f10 * f10 + f11 * f11;
    const double mean = 0.5 * (b00 + b11);
    const double radius = std::hypot(0.5 * (b00 - b11), b01);
    const double first = mean + radius, second = mean - radius;
    if (!(first > 0.0) || !(second > 0.0) || !(radius > 0.0))
        throw std::domain_error("B5.14 logarithmic strain requires two positive distinct in-plane stretches");
    const double beta = (std::log(first) - std::log(second)) / (first - second);
    const double alpha = (first * std::log(second) - second * std::log(first)) / (first - second);
    return {0.5 * (alpha + beta * b00), 0.5 * (alpha + beta * b11), 0.0, 0.5 * beta * b01, 0.0, 0.0};
}

fuelsim::UnstructuredHex8Mesh mesh() {
    return fuelsim::UnstructuredHex8Mesh(std::vector<fuelsim::CartesianPoint3>(nodes.begin(), nodes.end()),
        {{{{0, 1, 2, 3, 4, 5, 6, 7}}}},
        {1},
        {{1, "solid"}},
        {},
        {});
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    const fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 1.0, 2.0e11, 0.3, 0.0, 600.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    result.regions.push_back({"solid", "solid", properties, 0.0, 600.0, -1, "", fuelsim::StrainFormulation::finite});
    return result;
}

std::vector<double>
prescribed_state(const fuelsim::TransientProblem& problem, const Deformation& value, double temperature) {
    std::vector<double> result = problem.committed_solution();
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        const fuelsim::CartesianPoint3 current = current_position(nodes[node], value);
        result[dofs.dof(fuelsim::Field::temperature, node)] = temperature;
        result[dofs.dof(fuelsim::Field::displacement_x, node)] = current.x - nodes[node].x;
        result[dofs.dof(fuelsim::Field::displacement_y, node)] = current.y - nodes[node].y;
        result[dofs.dof(fuelsim::Field::displacement_z, node)] = current.z - nodes[node].z;
    }
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

std::array<double, 6> components(const fuelsim::SymmetricTensor3Values& value) {
    return {value.xx, value.yy, value.zz, value.xy, value.yz, value.xz};
}

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance) {
    if (metrics.has_relative_norm() && !fuelsim::test::relative_metrics_below(metrics, relative_tolerance))
        return false;
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

void print_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    if (metrics.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metrics);
    else
        fuelsim::test::print_absolute_metrics(name, metrics);
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests "
                     "<nodes.csv> <integration.csv> <energy.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim Abaqus B5.14 finite-strain elastic comparison\n");
        const std::vector<NodeReference> node_reference = read_nodes(argv[1]);
        const std::vector<IntegrationReference> integration_reference = read_integration(argv[2]);
        const std::vector<EnergyReference> energy_reference = read_energy(argv[3]);
        fuelsim::TransientProblem problem(definition(), mesh());
        const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
        const fuelsim::Hex8Geometry& geometry =
            fuelsim::cartesian::ProblemAccess::region_element_geometry(problem, 0, 0);
        std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
        std::array<fuelsim::test::FieldErrorMetrics, 26> integration_metrics;
        std::array<fuelsim::test::FieldErrorMetrics, 2> energy_metrics;
        fuelsim::test::FieldErrorMetrics logarithmic_reference_metric;
        double cumulative_elastic_energy = 0.0;
        bool passed = true;
        for (std::size_t stage = 1; stage <= stage_count; ++stage) {
            const Deformation value = deformation(stage);
            const double time = time_step * static_cast<double>(stage);
            std::vector<double> state = prescribed_state(problem, value, 600.0);
            problem.begin_time_step({time, 1.0, true});
            const std::vector<double> reaction = raw_residual(problem, state);
            for (const NodeReference& reference : node_reference) {
                if (reference.stage != stage)
                    continue;
                if (reference.node < 1 || reference.node > 8 || std::abs(reference.time - time) > 1.0e-12)
                    throw std::invalid_argument("Abaqus B5.14 nodal label or time lies outside the path");
                const std::size_t node = reference.node - 1;
                const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature,
                    fuelsim::Field::displacement_x,
                    fuelsim::Field::displacement_y,
                    fuelsim::Field::displacement_z};
                for (std::size_t field = 0; field < fields.size(); ++field) {
                    const std::size_t dof = dofs.dof(fields[field], node);
                    nodal_metrics[field].add(state[dof], reference.fields[field]);
                    nodal_metrics[4 + field].add(reaction[dof], reference.fields[4 + field]);
                }
            }
            problem.commit_time_step(state);
            cumulative_elastic_energy += problem.last_conservation_summary().elastic_energy_change;
            const fuelsim::CartesianMaterialHistory& material =
                fuelsim::cartesian::ProblemAccess::material_history(problem, 0, 0);
            const fuelsim::SymmetricTensor3Values exact_logarithmic = logarithmic_strain(value);
            for (const IntegrationReference& reference : integration_reference) {
                if (reference.stage != stage)
                    continue;
                std::size_t closest = 0;
                double closest_squared = std::numeric_limits<double>::max();
                for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                    const fuelsim::CartesianPoint3 point = current_position(geometry.points[q].position, value);
                    const double distance_squared = std::pow(point.x - reference.position.x, 2)
                                                    + std::pow(point.y - reference.position.y, 2)
                                                    + std::pow(point.z - reference.position.z, 2);
                    if (distance_squared < closest_squared) {
                        closest = q;
                        closest_squared = distance_squared;
                    }
                }
                integration_metrics[0].add(std::sqrt(closest_squared), 0.0);
                const std::array<double, 6> actual_stress = components(material[closest].stress);
                const std::array<double, 6> expected_stress = components(reference.stress);
                const std::array<double, 6> actual_elastic = material[closest].elastic_strain;
                const std::array<double, 6> expected_integrated = components(reference.integrated_strain);
                const std::array<double, 6> expected_elastic = components(reference.elastic_strain);
                const std::array<double, 6> actual_logarithmic = components(exact_logarithmic);
                const std::array<double, 6> expected_logarithmic = components(reference.logarithmic_strain);
                for (std::size_t component = 0; component < 6; ++component) {
                    integration_metrics[1 + component].add(actual_stress[component], expected_stress[component]);
                    integration_metrics[7 + component].add(actual_elastic[component], expected_integrated[component]);
                    integration_metrics[13 + component].add(actual_elastic[component], expected_elastic[component]);
                    integration_metrics[19 + component].add(actual_logarithmic[component],
                        expected_logarithmic[component]);
                    logarithmic_reference_metric.add(expected_integrated[component], expected_logarithmic[component]);
                }
                integration_metrics[25].add(geometry.points[closest].weighted_measure * value.stretch,
                    reference.integration_volume);
            }
            const EnergyReference& energy = energy_reference.at(stage - 1);
            if (energy.stage != stage || std::abs(energy.time - time) > 1.0e-12)
                throw std::invalid_argument("Abaqus B5.14 energy stage or time lies outside the path");
            energy_metrics[0].add(cumulative_elastic_energy, energy.internal_energy);
            energy_metrics[1].add(cumulative_elastic_energy, energy.strain_energy);
        }

        const std::array<std::string, 8> nodal_names = {"temperature",
            "displacement_x",
            "displacement_y",
            "displacement_z",
            "reaction_heat_flux",
            "reaction_force_x",
            "reaction_force_y",
            "reaction_force_z"};
        const std::array<double, 8> nodal_zero_tolerances =
            {1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-10, 1.0e-2, 1.0e-2, 1.0e-2};
        for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
            print_metrics("b514_" + nodal_names[field], nodal_metrics[field]);
            passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, nodal_zero_tolerances[field]),
                         "B5.14 " + nodal_names[field] + " metrics are below 0.1 percent")
                     && passed;
        }
        const std::array<std::string, 26> integration_names = {"coordinate",
            "stress_xx",
            "stress_yy",
            "stress_zz",
            "stress_xy",
            "stress_yz",
            "stress_xz",
            "integrated_strain_xx",
            "integrated_strain_yy",
            "integrated_strain_zz",
            "integrated_strain_xy",
            "integrated_strain_yz",
            "integrated_strain_xz",
            "elastic_strain_xx",
            "elastic_strain_yy",
            "elastic_strain_zz",
            "elastic_strain_xy",
            "elastic_strain_yz",
            "elastic_strain_xz",
            "logarithmic_strain_xx",
            "logarithmic_strain_yy",
            "logarithmic_strain_zz",
            "logarithmic_strain_xy",
            "logarithmic_strain_yz",
            "logarithmic_strain_xz",
            "integration_volume"};
        for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
            print_metrics("b514_" + integration_names[field], integration_metrics[field]);
            const double zero_tolerance = field == 0 ? 1.0e-12 : (field <= 6 ? 1.0e-2 : 1.0e-14);
            passed = check(metrics_pass(integration_metrics[field], 1.0e-3, zero_tolerance),
                         "B5.14 " + integration_names[field] + " metrics are below 0.1 percent")
                     && passed;
        }
        fuelsim::test::print_relative_metrics("b514_integrated_vs_exact_logarithmic_strain",
            logarithmic_reference_metric);
        passed = check(logarithmic_reference_metric.maximum_absolute_difference > 1.0e-8,
                     "B5.14 distinguishes Abaqus integrated strain E from exact logarithmic strain LE")
                 && passed;
        const std::array<std::string, 2> energy_names = {"internal_energy", "strain_energy"};
        for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
            print_metrics("b514_" + energy_names[field], energy_metrics[field]);
            passed = check(metrics_pass(energy_metrics[field], 1.0e-3, 1.0e-6),
                         "B5.14 " + energy_names[field] + " metrics are below 0.1 percent")
                     && passed;
        }
        if (!passed)
            return 1;
        std::cout << "[PASS] Abaqus B5.14 finite-strain elastic kinematics comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] Abaqus B5.14 finite-strain elastic test raised: " << error.what() << '\n';
        return 1;
    }
}
