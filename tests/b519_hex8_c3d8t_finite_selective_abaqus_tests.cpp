#include "c3d8_types.hpp"
#include "c3d8t.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "support/c3d_recovery.hpp"
#include "support/field_error_metrics.hpp"
#include "support/test_support.hpp"
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
constexpr std::size_t local_size = fuelsim::hex8_local_dof_count;
constexpr std::array<std::size_t, 8> gauss_to_material_node = {0, 1, 3, 2, 4, 5, 7, 6};

struct NodalReference final {
    fuelsim::CartesianPoint3 position;
    double temperature;
    fuelsim::CartesianPoint3 displacement;
    double reaction_flux;
    fuelsim::CartesianPoint3 reaction;
};

struct IntegrationReference final {
    fuelsim::CartesianPoint3 position;
    double temperature;
    fuelsim::SymmetricTensor3Values stress;
    double volume;
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

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size())
        throw std::invalid_argument("Incomplete Abaqus B5.19 row in " + path);
    return std::stod(values[column]);
}

std::array<NodalReference, 8> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.19 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,x_m,y_m,z_m,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n")
        throw std::invalid_argument("Unexpected Abaqus B5.19 nodal header in " + path);
    std::array<NodalReference, 8> result{};
    std::array<bool, 8> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 12)
            throw std::invalid_argument("Unexpected Abaqus B5.19 nodal column count in " + path);
        const double label_value = number(values, 0, path);
        const double rounded = std::round(label_value);
        if (std::abs(label_value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > 8.0)
            throw std::invalid_argument("Abaqus B5.19 node label lies outside 1 through 8");
        const std::size_t node = static_cast<std::size_t>(rounded) - 1;
        if (present[node])
            throw std::invalid_argument("Duplicate Abaqus B5.19 nodal row");
        present[node] = true;
        result[node] = {{number(values, 1, path), number(values, 2, path), number(values, 3, path)},
            number(values, 4, path),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            number(values, 8, path),
            {number(values, 9, path), number(values, 10, path), number(values, 11, path)}};
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Abaqus B5.19 nodal reference does not contain all eight nodes");
    return result;
}

std::array<IntegrationReference, 8> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus B5.19 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected = "integration_point,x_m,y_m,z_m,temperature_k,"
                                 "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
                                 "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,ivol_m3";
    if (line != expected)
        throw std::invalid_argument("Unexpected Abaqus B5.19 integration header in " + path);
    std::array<IntegrationReference, 8> result{};
    std::array<bool, 8> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 18)
            throw std::invalid_argument("Unexpected Abaqus B5.19 integration column count in " + path);
        const double point_value = number(values, 0, path);
        const double rounded = std::round(point_value);
        if (std::abs(point_value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > 8.0)
            throw std::invalid_argument("Abaqus B5.19 integration point lies outside 1 through 8");
        const std::size_t point = static_cast<std::size_t>(rounded) - 1;
        if (present[point])
            throw std::invalid_argument("Duplicate Abaqus B5.19 integration row");
        present[point] = true;
        result[point] = {{number(values, 1, path), number(values, 2, path), number(values, 3, path)},
            number(values, 4, path),
            {number(values, 5, path),
                number(values, 6, path),
                number(values, 7, path),
                number(values, 8, path),
                number(values, 10, path),
                number(values, 9, path)},
            number(values, 17, path)};
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Abaqus B5.19 integration reference does not contain all eight points");
    return result;
}

fuelsim::Hex8Coordinates unit_cube() {
    return {{{0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 1.0, 1.0}}};
}

fuelsim::Hex8LocalValues prescribed_state() {
    fuelsim::Hex8LocalValues result{};
    result.fill(600.0);
    const std::array<double, 8> ux = {0.0, 0.024, 0.034, -0.008, 0.006, 0.020, 0.046, -0.016};
    const std::array<double, 8> uy = {0.0, -0.004, 0.018, 0.022, -0.010, 0.008, 0.032, 0.014};
    const std::array<double, 8> uz = {0.0, 0.006, -0.008, 0.004, 0.028, 0.020, 0.036, 0.022};
    for (std::size_t node = 0; node < 8; ++node) {
        result[8 + node] = ux[node];
        result[16 + node] = uy[node];
        result[24 + node] = uz[node];
    }
    return result;
}

fuelsim::CartesianPoint3 current_position(const fuelsim::Hex8QuadraturePoint& point,
    const fuelsim::Hex8Coordinates& coordinates,
    const fuelsim::Hex8LocalValues& state) {
    fuelsim::CartesianPoint3 result{};
    for (std::size_t node = 0; node < 8; ++node) {
        result.x += point.shape[node] * (coordinates[node].x + state[8 + node]);
        result.y += point.shape[node] * (coordinates[node].y + state[16 + node]);
        result.z += point.shape[node] * (coordinates[node].z + state[24 + node]);
    }
    return result;
}

double squared_distance(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    const double dx = first.x - second.x, dy = first.y - second.y, dz = first.z - second.z;
    return dx * dx + dy * dy + dz * dz;
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
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests <nodes.csv> <integration.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const std::array<NodalReference, 8> nodal = read_nodes(argv[1]);
        const std::array<IntegrationReference, 8> integration = read_integration(argv[2]);
        const fuelsim::Hex8Coordinates coordinates = unit_cube();
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
        const fuelsim::Hex8LocalValues state = prescribed_state();
        const fuelsim::ThermoelasticProperties properties =
            fuelsim::test::thermoelastic(0.0, 1.0, 2.0e11, 0.25, 0.0, 600.0, 0.0, 0.0, 0.0, 1.0, 1.0);
        const fuelsim::CartesianRegionData data{fuelsim::IsotropicThermoelasticMaterial(properties),
            0.0,
            1.0,
            fuelsim::StrainFormulation::finite};

        fuelsim::Hex8LocalJacobian jacobian{};
        const fuelsim::Hex8LocalResidual residual =
            fuelsim::compute_c3d8_thermoelastic(data, geometry, state, nullptr, 0.0, &jacobian);
        const std::array<fuelsim::SymmetricTensor3Values, 8> stresses =
            fuelsim::compute_c3d8_stress(data, geometry, state);

        std::array<fuelsim::test::FieldErrorMetrics, 3> reaction_metrics, coordinate_metrics;
        fuelsim::test::FieldErrorMetrics reaction_flux_metrics;
        for (std::size_t node = 0; node < 8; ++node) {
            const std::array<double, 3> actual_reaction = {residual[8 + node],
                residual[16 + node],
                residual[24 + node]};
            const std::array<double, 3> reference_reaction = {nodal[node].reaction.x,
                nodal[node].reaction.y,
                nodal[node].reaction.z};
            const std::array<double, 3> actual_position = {coordinates[node].x + state[8 + node],
                coordinates[node].y + state[16 + node],
                coordinates[node].z + state[24 + node]};
            const std::array<double, 3> reference_position = {
                std::abs(nodal[node].position.x) < 1.0e-20 ? 0.0 : nodal[node].position.x,
                std::abs(nodal[node].position.y) < 1.0e-20 ? 0.0 : nodal[node].position.y,
                std::abs(nodal[node].position.z) < 1.0e-20 ? 0.0 : nodal[node].position.z};
            for (std::size_t component = 0; component < 3; ++component) {
                reaction_metrics[component].add(actual_reaction[component], reference_reaction[component]);
                coordinate_metrics[component].add(actual_position[component], reference_position[component]);
            }
            reaction_flux_metrics.add(residual[node], nodal[node].reaction_flux);
            if (std::abs(state[node] - nodal[node].temperature) > 1.0e-12)
                throw std::invalid_argument("Abaqus B5.19 nodal temperature does not match the prescribed state");
        }

        std::array<fuelsim::test::FieldErrorMetrics, 6> stress_metrics;
        std::array<fuelsim::test::FieldErrorMetrics, 3> point_coordinate_metrics;
        fuelsim::test::FieldErrorMetrics volume_metrics, point_temperature_metrics;
        const double current_volume =
            fuelsim::test::recover_c3d8t(geometry, state, {}, fuelsim::StrainFormulation::finite).current_volume;
        std::array<bool, 8> reference_used{};
        for (std::size_t q = 0; q < 8; ++q) {
            const fuelsim::CartesianPoint3 point_position = current_position(geometry.points[q], coordinates, state);
            std::size_t matched = 8;
            double minimum_distance = std::numeric_limits<double>::infinity();
            for (std::size_t candidate = 0; candidate < 8; ++candidate) {
                if (reference_used[candidate])
                    continue;
                const double distance = squared_distance(point_position, integration[candidate].position);
                if (distance < minimum_distance) {
                    minimum_distance = distance;
                    matched = candidate;
                }
            }
            if (matched == 8 || minimum_distance > 1.0e-20)
                throw std::invalid_argument("Abaqus B5.19 integration-point coordinates cannot be uniquely matched");
            reference_used[matched] = true;
            const std::array<double, 3> actual_position = {point_position.x, point_position.y, point_position.z};
            const std::array<double, 3> reference_position = {integration[matched].position.x,
                integration[matched].position.y,
                integration[matched].position.z};
            for (std::size_t component = 0; component < 3; ++component)
                point_coordinate_metrics[component].add(actual_position[component], reference_position[component]);
            const std::array<double, 6> actual_stress =
                {stresses[q].xx, stresses[q].yy, stresses[q].zz, stresses[q].xy, stresses[q].yz, stresses[q].xz};
            const std::array<double, 6> reference_stress = {integration[matched].stress.xx,
                integration[matched].stress.yy,
                integration[matched].stress.zz,
                integration[matched].stress.xy,
                integration[matched].stress.yz,
                integration[matched].stress.xz};
            for (std::size_t component = 0; component < 6; ++component)
                stress_metrics[component].add(actual_stress[component], reference_stress[component]);
            const double selective_volume =
                geometry.points[q].weighted_measure / geometry.reference_volume * current_volume;
            volume_metrics.add(selective_volume, integration[matched].volume);
            point_temperature_metrics.add(state[gauss_to_material_node[q]], integration[matched].temperature);
        }

        const std::array<double, 24> direction = {0.37,
            -0.21,
            0.14,
            -0.43,
            0.28,
            0.09,
            -0.16,
            0.31,
            -0.11,
            0.29,
            -0.35,
            0.18,
            0.42,
            -0.08,
            0.24,
            -0.32,
            0.19,
            0.33,
            -0.27,
            0.12,
            -0.38,
            0.23,
            0.07,
            -0.17};
        constexpr double epsilon = 1.0e-7;
        fuelsim::Hex8LocalValues plus = state, minus = state;
        for (std::size_t column = 0; column < 24; ++column) {
            plus[8 + column] += epsilon * direction[column];
            minus[8 + column] -= epsilon * direction[column];
        }
        const fuelsim::Hex8LocalResidual plus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, plus);
        const fuelsim::Hex8LocalResidual minus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, minus);
        double directional_difference_squared = 0.0, directional_reference_squared = 0.0;
        for (std::size_t row = 8; row < 32; ++row) {
            double exact = 0.0;
            for (std::size_t column = 0; column < 24; ++column)
                exact += jacobian[row * local_size + 8 + column] * direction[column];
            const double centered = (plus_residual[row] - minus_residual[row]) / (2.0 * epsilon);
            directional_difference_squared += (exact - centered) * (exact - centered);
            directional_reference_squared += centered * centered;
        }
        const double directional_relative_error =
            std::sqrt(directional_difference_squared / directional_reference_squared);

        const std::array<std::string, 3> axes = {"x", "y", "z"};
        const std::array<std::string, 6> tensor_names = {"s11", "s22", "s33", "s12", "s23", "s13"};
        for (std::size_t component = 0; component < 3; ++component) {
            print_metrics("b519_nodal_reaction_" + axes[component], reaction_metrics[component]);
            print_metrics("b519_nodal_coordinate_" + axes[component], coordinate_metrics[component]);
            print_metrics("b519_integration_coordinate_" + axes[component], point_coordinate_metrics[component]);
        }
        for (std::size_t component = 0; component < 6; ++component)
            print_metrics("b519_stress_" + tensor_names[component], stress_metrics[component]);
        print_metrics("b519_integration_volume", volume_metrics);
        print_metrics("b519_integration_temperature", point_temperature_metrics);
        print_metrics("b519_nodal_reaction_flux", reaction_flux_metrics);
        std::cout << "b519_mechanical_jacobian_directional_relative_error=" << directional_relative_error << '\n';

        bool passed = true;
        for (std::size_t component = 0; component < 3; ++component) {
            passed = check(metrics_pass(reaction_metrics[component], 1.0e-10, 1.0e-6),
                         "fuelsim finite selective-integration nodal reactions match Abaqus C3D8T")
                     && passed;
            passed = check(metrics_pass(coordinate_metrics[component], 1.0e-12, 1.0e-14),
                         "fuelsim nodal current coordinates match Abaqus C3D8T")
                     && passed;
            passed = check(metrics_pass(point_coordinate_metrics[component], 1.0e-12, 1.0e-14),
                         "fuelsim integration-point current coordinates match Abaqus C3D8T")
                     && passed;
        }
        for (std::size_t component = 0; component < 6; ++component)
            passed = check(metrics_pass(stress_metrics[component], 1.0e-10, 1.0e-6),
                         "fuelsim finite selective-integration stresses match Abaqus C3D8T")
                     && passed;
        passed = check(metrics_pass(volume_metrics, 1.0e-12, 1.0e-14),
                     "fuelsim current integration volumes match Abaqus C3D8T")
                 && passed;
        passed = check(metrics_pass(point_temperature_metrics, 1.0e-12, 1.0e-14),
                     "fuelsim integration-point temperatures match Abaqus C3D8T")
                 && passed;
        passed = check(reaction_flux_metrics.absolute_l2() < 1.0e-12,
                     "uniform temperature produces zero thermal reaction flux")
                 && passed;
        passed = check(directional_relative_error < 1.0e-7,
                     "the narrow automatic-differentiation tangent includes the finite selective-volume chain")
                 && passed;
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
