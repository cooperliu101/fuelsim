#include "c3d8_types.hpp"
#include "c3d8rt.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "support/test_support.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Matrix8 = std::array<double, 64>;

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

Matrix8 read_matrix(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8RT capacity reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "input_local_node,output_local_node,node,temperature_k,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus C3D8RT capacity header");
    Matrix8 result{};
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = split(line);
        const std::size_t column = std::stoul(fields[0]) - 1, row = std::stoul(fields[1]) - 1;
        result[row * 8 + column] = std::stod(fields[4]);
    }
    return result;
}

double relative_error(const Matrix8& actual, const Matrix8& expected) {
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t entry = 0; entry < actual.size(); ++entry) {
        numerator += (actual[entry] - expected[entry]) * (actual[entry] - expected[entry]);
        denominator += expected[entry] * expected[entry];
    }
    return std::sqrt(numerator / denominator);
}

bool check_case(const std::string& path, const fuelsim::Hex8Coordinates& coordinates) {
    const Matrix8 abaqus_total = read_matrix(path);
    const fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
    const fuelsim::CartesianRegionData data{fuelsim::IsotropicThermoelasticMaterial(properties),
        0.0,
        1.0,
        fuelsim::StrainFormulation::small,
        fuelsim::Hex8ElementFormulation::c3d8rt,
        300.0};
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8rt_geometry(coordinates);
    fuelsim::Hex8LocalValues committed{};
    for (std::size_t node = 0; node < 8; ++node)
        committed[node] = 300.0;
    const fuelsim::CartesianMaterialHistory history(1);
    Matrix8 fuelsim_total{}, fuelsim_capacity{}, expected_capacity{};
    for (std::size_t column = 0; column < 8; ++column) {
        fuelsim::Hex8LocalValues state = committed;
        state[column] += 1.0;
        const fuelsim::Hex8LocalResidual steady = fuelsim::compute_c3d8_thermoelastic(data, geometry, state);
        const fuelsim::Hex8LocalResidual transient =
            fuelsim::compute_c3d8_transient(data, geometry, state, committed, history, 1.0);
        for (std::size_t row = 0; row < 8; ++row) {
            fuelsim_total[row * 8 + column] = transient[row];
            fuelsim_capacity[row * 8 + column] = transient[row] - steady[row];
        }
        expected_capacity[column * 8 + column] = 6.0e6 * geometry.reduced_capacity_points[column].weighted_measure;
    }
    Matrix8 abaqus_capacity{};
    for (std::size_t entry = 0; entry < abaqus_total.size(); ++entry)
        abaqus_capacity[entry] = abaqus_total[entry] - (fuelsim_total[entry] - fuelsim_capacity[entry]);
    const double total_error = relative_error(fuelsim_total, abaqus_total);
    const double fuelsim_capacity_error = relative_error(fuelsim_capacity, expected_capacity);
    const double abaqus_capacity_error = relative_error(abaqus_capacity, expected_capacity);
    std::cout << "C3D8RT capacity " << path << ": total=" << total_error * 100.0
              << "%, fuelsim capacity=" << fuelsim_capacity_error * 100.0
              << "%, Abaqus capacity=" << abaqus_capacity_error * 100.0 << "%\n";
    return total_error < 1.0e-9 && fuelsim_capacity_error < 1.0e-12 && abaqus_capacity_error < 1.0e-9;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3)
        return 2;
    const fuelsim::Hex8Coordinates regular = {{{0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 1.0, 1.0}}};
    const fuelsim::Hex8Coordinates warped = {{{0.00, 0.00, 0.00},
        {1.20, 0.10, -0.05},
        {1.10, 1.00, 0.10},
        {-0.10, 0.90, 0.00},
        {0.05, -0.05, 1.00},
        {1.15, 0.00, 1.20},
        {1.00, 1.10, 1.10},
        {-0.05, 1.00, 0.90}}};
    return check_case(argv[1], regular) && check_case(argv[2], warped) ? 0 : 1;
}
