#include "c3d8_types.hpp"
#include "c3d8t.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "support/test_support.hpp"
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
struct NodalReference final {
    fuelsim::CartesianPoint3 current;
    double temperature;
    fuelsim::CartesianPoint3 displacement;
    double reaction_heat_flux;
};

struct IntegrationPointReference final {
    fuelsim::CartesianPoint3 current, heat_flux;
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

std::array<NodalReference, 8> read_nodal(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus finite-heat nodal reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,x_current_m,y_current_m,z_current_m,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus finite-heat nodal header in " + path);
    std::array<NodalReference, 8> result{};
    std::array<bool, 8> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 9)
            throw std::invalid_argument("Unexpected Abaqus finite-heat nodal column count in " + path);
        const std::size_t node = std::stoul(values[0]);
        if (node < 1 || node > 8 || present[node - 1])
            throw std::invalid_argument("Invalid or duplicate Abaqus finite-heat node label");
        present[node - 1] = true;
        result[node - 1] = {{std::stod(values[1]), std::stod(values[2]), std::stod(values[3])},
            std::stod(values[4]),
            {std::stod(values[5]), std::stod(values[6]), std::stod(values[7])},
            std::stod(values[8])};
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Abaqus finite-heat nodal reference must contain eight nodes");
    return result;
}

std::array<IntegrationPointReference, 8> read_integration_points(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus finite-heat integration-point reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "element,integration_point,x_current_m,y_current_m,z_current_m,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2")
        throw std::invalid_argument("Unexpected Abaqus finite-heat integration-point header in " + path);
    std::array<IntegrationPointReference, 8> result{};
    std::array<bool, 8> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 8 || std::stoul(values[0]) != 1)
            throw std::invalid_argument("Unexpected Abaqus finite-heat integration-point row in " + path);
        const std::size_t point = std::stoul(values[1]);
        if (point < 1 || point > 8 || present[point - 1])
            throw std::invalid_argument("Invalid or duplicate Abaqus finite-heat integration point");
        present[point - 1] = true;
        result[point - 1] = {{std::stod(values[2]), std::stod(values[3]), std::stod(values[4])},
            {std::stod(values[5]), std::stod(values[6]), std::stod(values[7])}};
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Abaqus finite-heat reference must contain eight integration points");
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

double relative_error(double actual, double expected) {
    return std::abs(actual - expected) / std::abs(expected);
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests <nodal.csv> <integration_points.csv>\n";
        return 2;
    }
    try {
        const std::array<NodalReference, 8> nodal = read_nodal(argv[1]);
        const std::array<IntegrationPointReference, 8> integration_points = read_integration_points(argv[2]);
        const fuelsim::Hex8Coordinates coordinates = unit_cube();
        const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates);
        fuelsim::Hex8LocalValues state{};
        double coordinate_maximum_difference = 0.0, displacement_maximum_difference = 0.0;
        for (std::size_t node = 0; node < 8; ++node) {
            const fuelsim::CartesianPoint3 expected_current{1.5 * coordinates[node].x,
                1.25 * coordinates[node].y,
                0.8 * coordinates[node].z};
            const fuelsim::CartesianPoint3 expected_displacement{expected_current.x - coordinates[node].x,
                expected_current.y - coordinates[node].y,
                expected_current.z - coordinates[node].z};
            coordinate_maximum_difference = std::max({coordinate_maximum_difference,
                std::abs(nodal[node].current.x - expected_current.x),
                std::abs(nodal[node].current.y - expected_current.y),
                std::abs(nodal[node].current.z - expected_current.z)});
            displacement_maximum_difference = std::max({displacement_maximum_difference,
                std::abs(nodal[node].displacement.x - expected_displacement.x),
                std::abs(nodal[node].displacement.y - expected_displacement.y),
                std::abs(nodal[node].displacement.z - expected_displacement.z)});
            state[node] = nodal[node].temperature;
            state[8 + node] = nodal[node].displacement.x;
            state[16 + node] = nodal[node].displacement.y;
            state[24 + node] = nodal[node].displacement.z;
        }
        const fuelsim::ThermoelasticProperties properties =
            fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
        const fuelsim::CartesianRegionData data{fuelsim::IsotropicThermoelasticMaterial(properties),
            0.0,
            1.0,
            fuelsim::StrainFormulation::finite};
        const fuelsim::Hex8LocalResidual fuelsim_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, state);
        double abaqus_hot_reaction = 0.0, fuelsim_hot_reaction = 0.0;
        for (std::size_t node : {std::size_t{1}, std::size_t{2}, std::size_t{5}, std::size_t{6}}) {
            abaqus_hot_reaction += nodal[node].reaction_heat_flux;
            fuelsim_hot_reaction += fuelsim_residual[node];
        }
        constexpr double current_heat_rate = 4.0 * (1.25 * 0.8) / 1.5 * 100.0;
        constexpr double current_heat_flux = -4.0 * 100.0 / 1.5;
        double heat_flux_maximum_difference = 0.0, transverse_heat_flux_maximum = 0.0;
        for (const IntegrationPointReference& point : integration_points) {
            heat_flux_maximum_difference =
                std::max(heat_flux_maximum_difference, std::abs(point.heat_flux.x - current_heat_flux));
            transverse_heat_flux_maximum =
                std::max({transverse_heat_flux_maximum, std::abs(point.heat_flux.y), std::abs(point.heat_flux.z)});
        }
        const double abaqus_current_error = relative_error(abaqus_hot_reaction, current_heat_rate);
        const double fuelsim_current_error = relative_error(fuelsim_hot_reaction, current_heat_rate);
        const double production_error = relative_error(fuelsim_hot_reaction, abaqus_hot_reaction);
        std::cout << "b51_current_coordinate_maximum_difference=" << coordinate_maximum_difference << '\n'
                  << "b51_displacement_maximum_difference=" << displacement_maximum_difference << '\n'
                  << "b51_abaqus_current_configuration_heat_rate_relative_error=" << abaqus_current_error << '\n'
                  << "b51_fuelsim_current_configuration_heat_rate_relative_error=" << fuelsim_current_error << '\n'
                  << "b51_production_heat_rate_relative_error=" << production_error << '\n'
                  << "b51_abaqus_heat_flux_maximum_difference=" << heat_flux_maximum_difference << '\n'
                  << "b51_abaqus_transverse_heat_flux_maximum=" << transverse_heat_flux_maximum << '\n';

        bool passed = true;
        passed = check(coordinate_maximum_difference < 2.0e-7 && displacement_maximum_difference < 2.0e-8,
                     "Abaqus reaches the prescribed finite deformation")
                 && passed;
        passed = check(abaqus_current_error < 1.0e-7 && heat_flux_maximum_difference < 2.0e-5
                           && transverse_heat_flux_maximum < 1.0e-10,
                     "Abaqus C3D8T evaluates conductivity on the current element dimensions")
                 && passed;
        passed = check(fuelsim_current_error < 1.0e-7,
                     "fuelsim finite-strain HEX8 evaluates conductivity on the current element dimensions")
                 && passed;
        passed = check(production_error < 1.0e-7, "fuelsim production finite-strain heat rate matches Abaqus C3D8T")
                 && passed;
        if (passed)
            std::cout << "[PASS] B5.1 Abaqus C3D8T finite-deformation heat identification\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
