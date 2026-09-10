#include "c3d8_types.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t node_count = fuelsim::hex8_node_count;
using Matrix8 = std::array<double, node_count * node_count>;

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

Matrix8 read_abaqus_total_matrix(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8T capacity reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "input_local_node,output_local_node,node,temperature_k,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus C3D8T capacity header in " + path);
    Matrix8 result{};
    std::array<bool, node_count * node_count> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5)
            throw std::invalid_argument("Unexpected Abaqus C3D8T capacity column count in " + path);
        const std::size_t column = static_cast<std::size_t>(std::stoul(values[0]));
        const std::size_t row = static_cast<std::size_t>(std::stoul(values[1]));
        const std::size_t node = static_cast<std::size_t>(std::stoul(values[2]));
        if (column < 1 || column > node_count || row < 1 || row > node_count || node != (column - 1) * node_count + row)
            throw std::invalid_argument("Abaqus C3D8T capacity row has inconsistent labels");
        const double temperature = std::stod(values[3]);
        const double expected_temperature = row == column ? 301.0 : 300.0;
        if (std::abs(temperature - expected_temperature) > 1.0e-12)
            throw std::invalid_argument("Abaqus C3D8T capacity temperature does not match the declared probe");
        const std::size_t index = (row - 1) * node_count + column - 1;
        if (present[index])
            throw std::invalid_argument("Duplicate Abaqus C3D8T capacity matrix entry");
        present[index] = true;
        result[index] = std::stod(values[4]);
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Abaqus C3D8T capacity reference must contain all 64 matrix entries");
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

fuelsim::ThermoelasticProperties properties() {
    return fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
}

double relative_frobenius_error(const Matrix8& actual, const Matrix8& expected) {
    double difference_squared = 0.0, expected_squared = 0.0;
    for (std::size_t entry = 0; entry < actual.size(); ++entry) {
        const double difference = actual[entry] - expected[entry];
        difference_squared += difference * difference;
        expected_squared += expected[entry] * expected[entry];
    }
    return std::sqrt(difference_squared / expected_squared);
}

double maximum_row_sum_error(const Matrix8& matrix, double expected) {
    double maximum = 0.0;
    for (std::size_t row = 0; row < node_count; ++row) {
        double sum = 0.0;
        for (std::size_t column = 0; column < node_count; ++column)
            sum += matrix[row * node_count + column];
        maximum = std::max(maximum, std::abs(sum - expected));
    }
    return maximum;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests <capacity.csv>\n";
        return 2;
    }
    try {
        const Matrix8 abaqus_total = read_abaqus_total_matrix(argv[1]);
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
        const fuelsim::CartesianRegionData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 0.0, 1.0};
        fuelsim::Hex8LocalValues committed_state{};
        for (std::size_t node = 0; node < node_count; ++node)
            committed_state[node] = 300.0;
        const fuelsim::CartesianMaterialHistory committed_material(geometry.points.size());
        Matrix8 conduction{}, fuelsim_total{}, fuelsim_capacity{}, analytic_consistent{}, analytic_lumped{};
        for (std::size_t column = 0; column < node_count; ++column) {
            fuelsim::Hex8LocalValues state = committed_state;
            state[column] += 1.0;
            const fuelsim::Hex8LocalResidual steady = fuelsim::compute_c3d8_thermoelastic(data, geometry, state);
            const fuelsim::Hex8LocalResidual transient =
                fuelsim::compute_c3d8_transient(data, geometry, state, committed_state, committed_material, 1.0);
            for (std::size_t row = 0; row < node_count; ++row) {
                const std::size_t entry = row * node_count + column;
                conduction[entry] = steady[row];
                fuelsim_total[entry] = transient[row];
                fuelsim_capacity[entry] = transient[row] - steady[row];
            }
        }
        constexpr double volumetric_capacity = 2000.0 * 3000.0;
        for (const fuelsim::Hex8QuadraturePoint& point : geometry.points)
            for (std::size_t row = 0; row < node_count; ++row)
                for (std::size_t column = 0; column < node_count; ++column)
                    analytic_consistent[row * node_count + column] +=
                        volumetric_capacity * point.weighted_measure * point.shape[row] * point.shape[column];
        constexpr double nodal_capacity = volumetric_capacity / 8.0;
        for (std::size_t node = 0; node < node_count; ++node)
            analytic_lumped[node * node_count + node] = nodal_capacity;

        Matrix8 abaqus_capacity{};
        for (std::size_t entry = 0; entry < abaqus_total.size(); ++entry)
            abaqus_capacity[entry] = abaqus_total[entry] - conduction[entry];
        const double abaqus_lumped_error = relative_frobenius_error(abaqus_capacity, analytic_lumped);
        const double fuelsim_lumped_error = relative_frobenius_error(fuelsim_capacity, analytic_lumped);
        const double fuelsim_total_error = relative_frobenius_error(fuelsim_total, abaqus_total);
        const double consistent_capacity_difference = relative_frobenius_error(analytic_consistent, abaqus_capacity);
        const double abaqus_row_sum_error = maximum_row_sum_error(abaqus_capacity, nodal_capacity);
        const double fuelsim_row_sum_error = maximum_row_sum_error(fuelsim_capacity, nodal_capacity);
        double maximum_abaqus_off_diagonal = 0.0, maximum_fuelsim_off_diagonal = 0.0;
        for (std::size_t row = 0; row < node_count; ++row)
            for (std::size_t column = 0; column < node_count; ++column)
                if (row != column) {
                    maximum_abaqus_off_diagonal =
                        std::max(maximum_abaqus_off_diagonal, std::abs(abaqus_capacity[row * node_count + column]));
                    maximum_fuelsim_off_diagonal =
                        std::max(maximum_fuelsim_off_diagonal, std::abs(fuelsim_capacity[row * node_count + column]));
                }

        std::cout << "b50_abaqus_lumped_capacity_relative_frobenius_error=" << abaqus_lumped_error << '\n'
                  << "b50_fuelsim_lumped_capacity_relative_frobenius_error=" << fuelsim_lumped_error << '\n'
                  << "b50_fuelsim_total_relative_frobenius_error=" << fuelsim_total_error << '\n'
                  << "b50_consistent_capacity_relative_frobenius_difference=" << consistent_capacity_difference << '\n'
                  << "b50_abaqus_capacity_row_sum_maximum_absolute_error=" << abaqus_row_sum_error << '\n'
                  << "b50_fuelsim_capacity_row_sum_maximum_absolute_error=" << fuelsim_row_sum_error << '\n'
                  << "b50_abaqus_capacity_maximum_off_diagonal=" << maximum_abaqus_off_diagonal << '\n'
                  << "b50_fuelsim_capacity_maximum_off_diagonal=" << maximum_fuelsim_off_diagonal << '\n';

        bool passed = true;
        passed = check(abaqus_lumped_error < 2.0e-12,
                     "Abaqus C3D8T transient response identifies the diagonal lumped heat-capacity matrix")
                 && passed;
        passed = check(fuelsim_lumped_error < 2.0e-12,
                     "fuelsim production HEX8 uses the diagonal lumped heat-capacity matrix")
                 && passed;
        passed = check(fuelsim_total_error < 2.0e-12,
                     "fuelsim conduction plus lumped heat capacity matches the Abaqus total reaction matrix")
                 && passed;
        passed = check(consistent_capacity_difference > 0.5,
                     "the nonuniform transient probe distinguishes the replaced consistent capacity matrix")
                 && passed;
        passed = check(abaqus_row_sum_error < 1.0e-6 && fuelsim_row_sum_error < 1.0e-6,
                     "both capacity matrices preserve the same uniform-heating nodal row sum")
                 && passed;
        passed = check(maximum_abaqus_off_diagonal < 1.0e-6 && maximum_fuelsim_off_diagonal < 1.0e-6,
                     "Abaqus and fuelsim lumped capacity matrices have zero off-diagonal entries")
                 && passed;
        if (passed)
            std::cout << "[PASS] B5.0 Abaqus C3D8T heat-capacity identification\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
