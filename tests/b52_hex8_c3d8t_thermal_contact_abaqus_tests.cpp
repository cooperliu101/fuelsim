#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/contact.hpp"
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
struct NodalReference final {
    double temperature, reaction_heat_flux;
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

std::array<NodalReference, 16> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus thermal-contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,temperature_k,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus thermal-contact header in " + path);
    std::array<NodalReference, 16> result{};
    std::array<bool, 16> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 3)
            throw std::invalid_argument("Unexpected Abaqus thermal-contact column count in " + path);
        const std::size_t node = std::stoul(values[0]);
        if (node < 1 || node > 16 || present[node - 1])
            throw std::invalid_argument("Invalid or duplicate Abaqus thermal-contact node label");
        present[node - 1] = true;
        result[node - 1] = {std::stod(values[1]), std::stod(values[2])};
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Abaqus thermal-contact reference must contain 16 nodes");
    return result;
}

fuelsim::Hex8Coordinates cube(double x_shift) {
    return {{{x_shift, 0.0, 0.0}, {x_shift + 1.0, 0.0, 0.0}, {x_shift + 1.0, 1.0, 0.0}, {x_shift, 1.0, 0.0},
        {x_shift, 0.0, 1.0}, {x_shift + 1.0, 0.0, 1.0}, {x_shift + 1.0, 1.0, 1.0}, {x_shift, 1.0, 1.0}}};
}

double relative_error(double actual, double expected) { return std::abs(actual - expected) / std::abs(expected); }
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests <nodal.csv>\n";
        return 2;
    }
    try {
        const std::array<NodalReference, 16> reference = read_reference(argv[1]);
        constexpr double conductivity = 4.0, conductance = 40.0, hot_temperature = 400.0, cold_temperature = 300.0;
        constexpr double expected_heat_rate =
            (hot_temperature - cold_temperature) / (1.0 / conductivity + 1.0 / conductance + 1.0 / conductivity);
        constexpr double expected_secondary_interface = hot_temperature - expected_heat_rate / conductivity;
        constexpr double expected_primary_interface = cold_temperature + expected_heat_rate / conductivity;
        const std::array<std::size_t, 4> hot_nodes{{0, 3, 4, 7}}, cold_nodes{{9, 10, 13, 14}},
            secondary_nodes{{1, 2, 6, 5}}, primary_nodes{{8, 11, 15, 12}};
        double abaqus_hot_rate = 0.0, abaqus_cold_rate = 0.0, temperature_maximum_difference = 0.0;
        for (std::size_t node : hot_nodes) abaqus_hot_rate += reference[node].reaction_heat_flux;
        for (std::size_t node : cold_nodes) abaqus_cold_rate += reference[node].reaction_heat_flux;
        for (std::size_t node : hot_nodes)
            temperature_maximum_difference =
                std::max(temperature_maximum_difference, std::abs(reference[node].temperature - hot_temperature));
        for (std::size_t node : cold_nodes)
            temperature_maximum_difference =
                std::max(temperature_maximum_difference, std::abs(reference[node].temperature - cold_temperature));
        for (std::size_t node : secondary_nodes)
            temperature_maximum_difference = std::max(
                temperature_maximum_difference, std::abs(reference[node].temperature - expected_secondary_interface));
        for (std::size_t node : primary_nodes)
            temperature_maximum_difference = std::max(
                temperature_maximum_difference, std::abs(reference[node].temperature - expected_primary_interface));

        const fuelsim::ThermoelasticProperties properties =
            fuelsim::test::thermoelastic(0.0, conductivity, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
        const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties), 0.0, 1.0};
        std::array<fuelsim::Hex8LocalValues, 2> states{};
        for (std::size_t node = 0; node < 8; ++node) {
            states[0][node] = reference[node].temperature;
            states[1][node] = reference[8 + node].temperature;
        }
        const fuelsim::Hex8Coordinates secondary_cube = cube(0.0), primary_cube = cube(1.1);
        const fuelsim::Hex8LocalResidual secondary_volume =
            fuelsim::compute_hex8_thermoelastic(data, fuelsim::make_hex8_geometry(secondary_cube), states[0]);
        const fuelsim::Hex8LocalResidual primary_volume =
            fuelsim::compute_hex8_thermoelastic(data, fuelsim::make_hex8_geometry(primary_cube), states[1]);
        const fuelsim::Quad4FaceCoordinates secondary_face = {{{1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0},
                                                {1.0, 0.0, 1.0}}},
                                            primary_face = {
                                                {{1.1, 0.0, 0.0}, {1.1, 1.0, 0.0}, {1.1, 1.0, 1.0}, {1.1, 0.0, 1.0}}};
        const fuelsim::Quad4FaceGeometry face_geometry = fuelsim::make_quad4_face_geometry(secondary_face);
        fuelsim::Quad4SurfaceContactLocalValues contact_state{};
        for (std::size_t node = 0; node < 4; ++node) {
            contact_state[node] = reference[secondary_nodes[node]].temperature;
            contact_state[4 + node] = reference[primary_nodes[node]].temperature;
        }
        fuelsim::Quad4SurfaceContactLocalResidual contact_residual{};
        double value_heat_rate = 0.0, value_measure = 0.0, gap_maximum_difference = 0.0;
        for (const fuelsim::Quad4FaceQuadraturePoint& point : face_geometry.points) {
            const fuelsim::Quad4ToQuad4HeatGeometry contact_geometry{
                secondary_face, primary_face, point.shape, point.derivative_xi, point.derivative_eta, 1.0};
            const fuelsim::Quad4SurfaceContactLocalResidual contribution =
                fuelsim::compute_quad4_to_quad4_gap_heat({conductance * 0.1, 1.0e-6}, contact_geometry, contact_state);
            for (std::size_t dof = 0; dof < contact_residual.size(); ++dof) contact_residual[dof] += contribution[dof];
            const fuelsim::CartesianHeatQuadratureValue value = fuelsim::compute_quad4_to_quad4_gap_heat_value(
                {conductance * 0.1, 1.0e-6}, contact_geometry, contact_state);
            value_heat_rate += value.heat_flux * value.weighted_measure;
            value_measure += value.weighted_measure;
            gap_maximum_difference = std::max(gap_maximum_difference, std::abs(value.gap - 0.1));
        }
        double secondary_contact_rate = 0.0, primary_contact_rate = 0.0, free_residual_maximum = 0.0;
        for (std::size_t node = 0; node < 4; ++node) {
            secondary_contact_rate += contact_residual[node];
            primary_contact_rate += contact_residual[4 + node];
            free_residual_maximum = std::max(
                free_residual_maximum, std::abs(secondary_volume[secondary_nodes[node]] + contact_residual[node]));
            free_residual_maximum = std::max(
                free_residual_maximum, std::abs(primary_volume[primary_nodes[node] - 8] + contact_residual[4 + node]));
        }
        double fuelsim_hot_rate = 0.0, fuelsim_cold_rate = 0.0;
        for (std::size_t node : hot_nodes) fuelsim_hot_rate += secondary_volume[node];
        for (std::size_t node : cold_nodes) fuelsim_cold_rate += primary_volume[node - 8];

        const double abaqus_heat_rate_error = relative_error(abaqus_hot_rate, expected_heat_rate);
        const double fuelsim_heat_rate_error = relative_error(fuelsim_hot_rate, expected_heat_rate);
        const double code_to_code_heat_rate_error = relative_error(fuelsim_hot_rate, abaqus_hot_rate);
        const double abaqus_balance = std::abs(abaqus_hot_rate + abaqus_cold_rate);
        const double fuelsim_balance = std::max(
            std::abs(fuelsim_hot_rate + fuelsim_cold_rate), std::abs(secondary_contact_rate + primary_contact_rate));
        std::cout << "b52_temperature_maximum_absolute_difference=" << temperature_maximum_difference << '\n'
                  << "b52_abaqus_heat_rate_relative_error=" << abaqus_heat_rate_error << '\n'
                  << "b52_fuelsim_heat_rate_relative_error=" << fuelsim_heat_rate_error << '\n'
                  << "b52_code_to_code_heat_rate_relative_error=" << code_to_code_heat_rate_error << '\n'
                  << "b52_abaqus_heat_balance_absolute=" << abaqus_balance << '\n'
                  << "b52_fuelsim_heat_balance_absolute=" << fuelsim_balance << '\n'
                  << "b52_fuelsim_free_node_residual_maximum=" << free_residual_maximum << '\n'
                  << "b52_fuelsim_free_node_residual_relative=" << free_residual_maximum / expected_heat_rate << '\n'
                  << "b52_fuelsim_contact_value_heat_rate=" << value_heat_rate << '\n'
                  << "b52_fuelsim_contact_value_heat_rate_relative_error="
                  << relative_error(value_heat_rate, expected_heat_rate) << '\n'
                  << "b52_fuelsim_contact_measure=" << value_measure << '\n'
                  << "b52_fuelsim_gap_maximum_difference=" << gap_maximum_difference << '\n';

        bool passed = true;
        passed = check(temperature_maximum_difference < 2.0e-5,
                     "Abaqus interface and boundary temperatures match the one-dimensional resistance solution") &&
                 passed;
        passed = check(abaqus_heat_rate_error < 1.0e-7 && fuelsim_heat_rate_error < 1.0e-7 &&
                           code_to_code_heat_rate_error < 1.0e-7,
                     "Abaqus and fuelsim match the fixed-gap contact heat rate") &&
                 passed;
        passed = check(abaqus_balance < 1.0e-12 && fuelsim_balance < 1.0e-11,
                     "both codes preserve strict two-sided thermal balance") &&
                 passed;
        passed = check(free_residual_maximum / expected_heat_rate < 2.0e-7 &&
                           relative_error(value_heat_rate, expected_heat_rate) < 1.0e-6 &&
                           std::abs(value_measure - 1.0) < 1.0e-14 && gap_maximum_difference < 1.0e-14,
                     "fuelsim volume and surface-contact kernels satisfy the matched free-node equilibrium") &&
                 passed;
        if (passed) std::cout << "[PASS] B5.2 Abaqus C3D8T fixed-gap thermal contact\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
