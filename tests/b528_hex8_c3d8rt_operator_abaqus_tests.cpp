#include "c3d8_types.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "support/test_support.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t local_size = fuelsim::hex8_local_dof_count;

struct Step final {
    fuelsim::Hex8LocalValues state{};
    fuelsim::Hex8LocalResidual reaction{};
};

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

std::map<std::string, Step> read_steps(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8RT operator reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,node,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n"
        && line != "case,node,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n")
        throw std::invalid_argument("Unexpected Abaqus C3D8RT operator header");
    std::map<std::string, Step> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = split(line);
        if (fields.size() != 10)
            throw std::invalid_argument("Unexpected Abaqus C3D8RT operator row");
        const std::size_t node = static_cast<std::size_t>(std::stoul(fields[1])) - 1;
        Step& step = result[fields[0]];
        step.state[node] = std::stod(fields[2]);
        step.state[8 + node] = std::stod(fields[3]);
        step.state[16 + node] = std::stod(fields[4]);
        step.state[24 + node] = std::stod(fields[5]);
        step.reaction[node] = std::stod(fields[6]);
        step.reaction[8 + node] = std::stod(fields[7]);
        step.reaction[16 + node] = std::stod(fields[8]);
        step.reaction[24 + node] = std::stod(fields[9]);
    }
    if (result.size() != 65)
        throw std::invalid_argument("Abaqus C3D8RT operator must contain 65 steps");
    return result;
}

fuelsim::Hex8Coordinates regular_coordinates() {
    return {{{0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 1.0, 1.0}}};
}

fuelsim::Hex8Coordinates warped_coordinates() {
    return {{{0.00, 0.00, 0.00},
        {1.20, 0.10, -0.05},
        {1.10, 1.00, 0.10},
        {-0.10, 0.90, 0.00},
        {0.05, -0.05, 1.00},
        {1.15, 0.00, 1.20},
        {1.00, 1.10, 1.10},
        {-0.05, 1.00, 0.90}}};
}

fuelsim::Hex8Coordinates holdout_coordinates() {
    return {{{0.00, 0.00, 0.00},
        {1.40, -0.10, 0.12},
        {1.18, 1.15, -0.08},
        {-0.15, 0.92, 0.18},
        {0.12, -0.08, 0.95},
        {1.30, 0.18, 1.28},
        {0.93, 1.26, 1.05},
        {-0.12, 1.08, 1.22}}};
}

double block_relative_error(const fuelsim::Hex8LocalJacobian& actual,
    const fuelsim::Hex8LocalJacobian& expected,
    std::size_t row_begin,
    std::size_t row_end,
    std::size_t column_begin,
    std::size_t column_end) {
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t row = row_begin; row < row_end; ++row)
        for (std::size_t column = column_begin; column < column_end; ++column) {
            const double difference = actual[row * local_size + column] - expected[row * local_size + column];
            numerator += difference * difference;
            denominator += expected[row * local_size + column] * expected[row * local_size + column];
        }
    return std::sqrt(numerator / denominator);
}

double vector_relative_error(const fuelsim::Hex8LocalResidual& actual,
    const fuelsim::Hex8LocalResidual& expected,
    std::size_t begin,
    std::size_t end) {
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t row = begin; row < end; ++row) {
        const double difference = actual[row] - expected[row];
        numerator += difference * difference;
        denominator += expected[row] * expected[row];
    }
    return std::sqrt(numerator / denominator);
}

bool check_case(const std::string& path,
    const fuelsim::Hex8Coordinates& coordinates,
    double thermal_limit,
    fuelsim::StrainFormulation strain_formulation = fuelsim::StrainFormulation::small) {
    const std::map<std::string, Step> steps = read_steps(path);
    const Step& base = steps.at("BASE");
    const fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
    const fuelsim::CartesianRegionData data{fuelsim::IsotropicThermoelasticMaterial(properties),
        0.0,
        0.0,
        strain_formulation,
        fuelsim::Hex8ElementFormulation::c3d8rt,
        300.0};
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
    fuelsim::Hex8LocalJacobian actual{};
    const fuelsim::Hex8LocalResidual actual_residual =
        fuelsim::compute_c3d8_thermoelastic(data, geometry, base.state, nullptr, 0.0, &actual);
    fuelsim::Hex8LocalJacobian expected{};
    for (std::size_t column = 0; column < local_size; ++column) {
        std::ostringstream plus_name, minus_name;
        plus_name << 'D' << (column < 10 ? "0" : "") << column << "_PLUS";
        minus_name << 'D' << (column < 10 ? "0" : "") << column << "_MINUS";
        const double perturbation = column < 8 ? 1.0e-3 : 1.0e-7;
        for (std::size_t row = 0; row < local_size; ++row)
            expected[row * local_size + column] =
                (steps.at(plus_name.str()).reaction[row] - steps.at(minus_name.str()).reaction[row])
                / (2.0 * perturbation);
    }
    const double thermal_error = block_relative_error(actual, expected, 0, 8, 0, 8);
    const double mechanical_error = block_relative_error(actual, expected, 8, 32, 8, 32);
    const double coupling_error = block_relative_error(actual, expected, 8, 32, 0, 8);
    const double thermal_residual_error = vector_relative_error(actual_residual, base.reaction, 0, 8);
    const double mechanical_residual_error = vector_relative_error(actual_residual, base.reaction, 8, 32);

    fuelsim::Hex8LocalJacobian finite_difference{};
    for (std::size_t column = 0; column < local_size; ++column) {
        const double perturbation = column < 8 ? 1.0e-4 : 1.0e-8;
        fuelsim::Hex8LocalValues plus = base.state, minus = base.state;
        plus[column] += perturbation;
        minus[column] -= perturbation;
        const fuelsim::Hex8LocalResidual plus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, plus);
        const fuelsim::Hex8LocalResidual minus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, minus);
        for (std::size_t row = 0; row < local_size; ++row)
            finite_difference[row * local_size + column] =
                (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
    }
    const double local_thermal_error = block_relative_error(actual, finite_difference, 0, 8, 0, 8);
    const double local_mechanical_error = block_relative_error(actual, finite_difference, 8, 32, 8, 32);
    const double local_coupling_error = block_relative_error(actual, finite_difference, 8, 32, 0, 8);
    std::cout << "C3D8RT operator " << path << ": Abaqus residual thermal=" << thermal_residual_error * 100.0
              << "%, residual mechanical=" << mechanical_residual_error * 100.0
              << "%, tangent thermal=" << thermal_error * 100.0 << "%, mechanical=" << mechanical_error * 100.0
              << "%, coupling=" << coupling_error * 100.0
              << "%; centered-difference thermal=" << local_thermal_error * 100.0
              << "%, mechanical=" << local_mechanical_error * 100.0 << "%, coupling=" << local_coupling_error * 100.0
              << "%\n";
    const bool finite = strain_formulation == fuelsim::StrainFormulation::finite;
    return thermal_residual_error < (finite ? 1.0e-3 : thermal_limit)
           && mechanical_residual_error < (finite ? 1.0e-2 : 1.0e-6) && thermal_error < thermal_limit
           && mechanical_error < (finite ? 1.0e-2 : 1.0e-6) && coupling_error < (finite ? 1.0e-4 : 1.0e-6)
           && local_thermal_error < 1.0e-9 && local_mechanical_error < 2.0e-7 && local_coupling_error < 2.0e-7;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests regular.csv warped.csv holdout.csv "
                     "finite_regular.csv finite_warped.csv\n";
        return 2;
    }
    bool ok = check_case(argv[1], regular_coordinates(), 1.0e-8);
    ok = check_case(argv[2], warped_coordinates(), 0.01) && ok;
    ok = check_case(argv[3], holdout_coordinates(), 0.01) && ok;
    ok = check_case(argv[4], regular_coordinates(), 1.0e-3, fuelsim::StrainFormulation::finite) && ok;
    ok = check_case(argv[5], warped_coordinates(), 1.0e-3, fuelsim::StrainFormulation::finite) && ok;
    return ok ? 0 : 1;
}
