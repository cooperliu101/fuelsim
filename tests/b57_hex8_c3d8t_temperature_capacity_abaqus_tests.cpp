#include "cartesian3d_hex8.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t node_count = fuelsim::hex8_node_count;
using NodalValues = std::array<double, node_count>;
using Matrix8 = std::array<double, node_count * node_count>;

struct StateReference final {
    NodalValues temperature{};
    std::array<NodalValues, 2> reaction{};
    std::array<std::array<bool, node_count>, 2> present{};
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

std::map<std::string, StateReference> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus temperature-capacity reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "state,copy,local_node,node,temperature_k,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus temperature-capacity header in " + path);
    std::map<std::string, StateReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 6)
            throw std::invalid_argument("Unexpected Abaqus temperature-capacity column count in " + path);
        const std::size_t copy = values[1] == "steady" ? 0 : values[1] == "transient" ? 1 : 2;
        const std::size_t node = static_cast<std::size_t>(std::stoul(values[2]));
        if (copy > 1 || node < 1 || node > node_count)
            throw std::invalid_argument("Invalid Abaqus temperature-capacity copy or node label");
        StateReference& state = result[values[0]];
        if (state.present[copy][node - 1])
            throw std::invalid_argument("Duplicate temperature-capacity row");
        state.present[copy][node - 1] = true;
        const double temperature = std::stod(values[4]);
        if (copy == 0)
            state.temperature[node - 1] = temperature;
        else if (std::abs(state.temperature[node - 1] - temperature) > 1.0e-12)
            throw std::invalid_argument("Steady and transient temperature-capacity copies differ in temperature");
        state.reaction[copy][node - 1] = std::stod(values[5]);
    }
    if (result.size() != 17)
        throw std::invalid_argument("Temperature-capacity reference must contain 17 states");
    for (const auto& [name, state] : result)
        for (const auto& copy : state.present)
            if (std::find(copy.begin(), copy.end(), false) != copy.end())
                throw std::invalid_argument("Incomplete temperature-capacity state " + name);
    return result;
}

fuelsim::Hex8Coordinates distorted_coordinates() {
    return {{{0.00, 0.00, 0.00},
        {1.20, 0.10, -0.05},
        {1.10, 1.00, 0.10},
        {-0.10, 0.90, 0.00},
        {0.05, -0.10, 1.00},
        {1.30, 0.00, 1.10},
        {1.00, 1.20, 0.90},
        {-0.20, 1.00, 1.20}}};
}

fuelsim::ThermoelasticProperties properties() {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>(*result.functions);
    functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
        {{"conductivity", 4.0},
            {"density", 2000.0},
            {"specific_heat", 3000.0},
            {"reference_temperature", 300.0},
            {"conductivity_temperature_coefficient", 0.0},
            {"density_temperature_coefficient", -1.0},
            {"specific_heat_temperature_coefficient", 4.0}});
    result.functions = std::move(functions);
    return result;
}

NodalValues capacity_reference(const StateReference& state) {
    NodalValues result{};
    for (std::size_t node = 0; node < node_count; ++node)
        result[node] = state.reaction[1][node] - state.reaction[0][node];
    return result;
}

double relative_error(const NodalValues& actual, const NodalValues& expected) {
    double difference_squared = 0.0, expected_squared = 0.0;
    for (std::size_t node = 0; node < node_count; ++node) {
        const double difference = actual[node] - expected[node];
        difference_squared += difference * difference;
        expected_squared += expected[node] * expected[node];
    }
    return std::sqrt(difference_squared / expected_squared);
}

double relative_error(const Matrix8& actual, const Matrix8& expected) {
    double difference_squared = 0.0, expected_squared = 0.0;
    for (std::size_t entry = 0; entry < actual.size(); ++entry) {
        const double difference = actual[entry] - expected[entry];
        difference_squared += difference * difference;
        expected_squared += expected[entry] * expected[entry];
    }
    return std::sqrt(difference_squared / expected_squared);
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests <capacity.csv>\n";
        return 2;
    }
    try {
        constexpr std::array<double, node_count> old_temperature =
            {310.0, 330.0, 350.0, 370.0, 390.0, 410.0, 430.0, 450.0};
        constexpr double perturbation = 1.0e-3;
        const std::map<std::string, StateReference> reference = read_reference(argv[1]);
        const StateReference& base = reference.at("BASE");
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(distorted_coordinates());
        const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 0.0, 1.0};
        fuelsim::Hex8LocalValues committed_state{}, state{};
        for (std::size_t node = 0; node < node_count; ++node) {
            committed_state[node] = old_temperature[node];
            state[node] = base.temperature[node];
        }
        const fuelsim::CartesianMaterialHistory history(geometry.points.size());
        fuelsim::Hex8LocalJacobian transient_jacobian{}, steady_jacobian{};
        const fuelsim::Hex8LocalResidual transient =
            fuelsim::compute_hex8_transient(data, geometry, state, committed_state, history, 1.0, &transient_jacobian);
        const fuelsim::Hex8LocalResidual steady = fuelsim::compute_hex8_transient(data,
            geometry,
            state,
            committed_state,
            history,
            1.0,
            &steady_jacobian,
            false);
        NodalValues fuelsim_capacity{}, analytic_capacity{}, enthalpy_capacity{};
        for (std::size_t node = 0; node < node_count; ++node) {
            fuelsim_capacity[node] = transient[node] - steady[node];
            const double temperature = state[node], increment = temperature - committed_state[node];
            const double density = 2000.0 - (temperature - 300.0);
            const double specific_heat = 3000.0 + 4.0 * (temperature - 300.0);
            analytic_capacity[node] =
                geometry.capacity_points[node].weighted_measure * density * specific_heat * increment;
            const double integrated_specific_heat =
                3000.0 * increment
                + 2.0
                      * ((temperature - 300.0) * (temperature - 300.0)
                          - (committed_state[node] - 300.0) * (committed_state[node] - 300.0));
            enthalpy_capacity[node] =
                geometry.capacity_points[node].weighted_measure * density * integrated_specific_heat;
        }

        Matrix8 abaqus_jacobian{}, fuelsim_jacobian{};
        for (std::size_t column = 0; column < node_count; ++column) {
            std::ostringstream plus_name, minus_name;
            plus_name << 'D' << (column < 10 ? "0" : "") << column << "_PLUS";
            minus_name << 'D' << (column < 10 ? "0" : "") << column << "_MINUS";
            const NodalValues plus = capacity_reference(reference.at(plus_name.str()));
            const NodalValues minus = capacity_reference(reference.at(minus_name.str()));
            for (std::size_t row = 0; row < node_count; ++row)
                abaqus_jacobian[row * node_count + column] = (plus[row] - minus[row]) / (2.0 * perturbation);
        }
        for (std::size_t row = 0; row < node_count; ++row)
            for (std::size_t column = 0; column < node_count; ++column)
                fuelsim_jacobian[row * node_count + column] =
                    transient_jacobian[row * 32 + column] - steady_jacobian[row * 32 + column];

        const NodalValues abaqus_capacity = capacity_reference(base);
        const double residual_error = relative_error(fuelsim_capacity, abaqus_capacity);
        const double analytic_error = relative_error(analytic_capacity, abaqus_capacity);
        const double enthalpy_difference = relative_error(enthalpy_capacity, abaqus_capacity);
        const double jacobian_error = relative_error(fuelsim_jacobian, abaqus_jacobian);
        double abaqus_off_diagonal = 0.0, fuelsim_off_diagonal = 0.0;
        for (std::size_t row = 0; row < node_count; ++row)
            for (std::size_t column = 0; column < node_count; ++column)
                if (row != column) {
                    abaqus_off_diagonal =
                        std::max(abaqus_off_diagonal, std::abs(abaqus_jacobian[row * node_count + column]));
                    fuelsim_off_diagonal =
                        std::max(fuelsim_off_diagonal, std::abs(fuelsim_jacobian[row * node_count + column]));
                }

        std::cout << std::scientific << "b57_temperature_capacity_residual_relative_error=" << residual_error << '\n'
                  << "b57_temperature_capacity_analytic_relative_error=" << analytic_error << '\n'
                  << "b57_temperature_capacity_jacobian_relative_frobenius_error=" << jacobian_error << '\n'
                  << "b57_enthalpy_increment_relative_difference=" << enthalpy_difference << '\n'
                  << "b57_abaqus_capacity_maximum_off_diagonal=" << abaqus_off_diagonal << '\n'
                  << "b57_fuelsim_capacity_maximum_off_diagonal=" << fuelsim_off_diagonal << '\n';
        const bool passed =
            check(residual_error < 2.0e-12, "Fuelsim temperature-dependent lumped-capacity residual matches Abaqus")
            && check(analytic_error < 2.0e-12,
                "Abaqus uses current rho(T) times cp(T) times the backward-Euler temperature increment")
            && check(jacobian_error < 2.0e-9,
                "Fuelsim temperature-dependent lumped-capacity Jacobian matches the Abaqus centered tangent")
            && check(enthalpy_difference > 1.0e-3,
                "the finite temperature increments distinguish the identified rule from an integrated-enthalpy rule")
            && check(abaqus_off_diagonal < 2.0e-4 && fuelsim_off_diagonal == 0.0,
                "Abaqus and Fuelsim temperature-dependent capacity tangents remain diagonal");
        if (passed)
            std::cout << "[PASS] B5.7 Abaqus C3D8T temperature-dependent heat capacity\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
