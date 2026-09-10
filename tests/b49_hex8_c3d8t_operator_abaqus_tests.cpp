#include "cartesian3d_hex8.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t local_size = fuelsim::hex8_local_dof_count;
constexpr std::array<std::size_t, 8> gauss_to_material_node = {0, 1, 3, 2, 4, 5, 7, 6};

struct NodalStep final {
    fuelsim::Hex8LocalValues state{};
    fuelsim::Hex8LocalResidual reaction{};
    std::array<bool, 8> present{};
};

struct IntegrationPointReference final {
    fuelsim::CartesianPoint3 point;
    double output_temperature;
    std::array<double, 3> heat_flux;
    fuelsim::SymmetricTensor3Values strain, stress;
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
        throw std::invalid_argument("Incomplete Abaqus C3D8T row in " + path);
    return std::stod(values[column]);
}

std::map<std::string, NodalStep> read_nodal_steps(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8T nodal reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,node,temperature_k,ux_m,uy_m,uz_m,rfl_w,rf_x_n,rf_y_n,rf_z_n")
        throw std::invalid_argument("Unexpected Abaqus C3D8T nodal header in " + path);
    std::map<std::string, NodalStep> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 10)
            throw std::invalid_argument("Unexpected Abaqus C3D8T nodal column count in " + path);
        const std::size_t node = static_cast<std::size_t>(number(values, 1, path));
        if (node < 1 || node > 8)
            throw std::invalid_argument("Abaqus C3D8T node label lies outside 1 through 8");
        NodalStep& step = result[values[0]];
        if (step.present[node - 1])
            throw std::invalid_argument("Duplicate Abaqus C3D8T nodal row");
        step.present[node - 1] = true;
        step.state[node - 1] = number(values, 2, path);
        step.state[8 + node - 1] = number(values, 3, path);
        step.state[16 + node - 1] = number(values, 4, path);
        step.state[24 + node - 1] = number(values, 5, path);
        step.reaction[node - 1] = number(values, 6, path);
        step.reaction[8 + node - 1] = number(values, 7, path);
        step.reaction[16 + node - 1] = number(values, 8, path);
        step.reaction[24 + node - 1] = number(values, 9, path);
    }
    if (result.size() != 65)
        throw std::invalid_argument("Abaqus C3D8T nodal reference must contain 65 steps");
    for (const auto& entry : result)
        if (std::find(entry.second.present.begin(), entry.second.present.end(), false) != entry.second.present.end())
            throw std::invalid_argument("Abaqus C3D8T step does not contain all eight nodes");
    return result;
}

std::vector<IntegrationPointReference> read_integration_points(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8T integration-point reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
           "e11,e22,e33,e12,e13,e23,s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa")
        throw std::invalid_argument("Unexpected Abaqus C3D8T integration-point header in " + path);
    std::vector<IntegrationPointReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 21)
            throw std::invalid_argument("Unexpected Abaqus C3D8T integration-point column count in " + path);
        if (static_cast<std::size_t>(number(values, 0, path)) != 1)
            throw std::invalid_argument("Abaqus C3D8T reference contains an unexpected element label");
        result.push_back({{number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            number(values, 5, path),
            {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
            {number(values, 9, path),
                number(values, 10, path),
                number(values, 11, path),
                0.5 * number(values, 12, path),
                0.5 * number(values, 14, path),
                0.5 * number(values, 13, path)},
            {number(values, 15, path),
                number(values, 16, path),
                number(values, 17, path),
                number(values, 18, path),
                number(values, 20, path),
                number(values, 19, path)}});
    }
    if (result.size() != 8)
        throw std::invalid_argument("Abaqus C3D8T reference must contain eight integration points");
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

fuelsim::ThermoelasticProperties properties(bool temperature_dependent) {
    fuelsim::ThermoelasticProperties result = fuelsim::test::thermoelastic(0.0,
        4.0,
        2.0e11,
        0.25,
        1.2e-5,
        300.0,
        temperature_dependent ? -1.0e8 : 0.0,
        temperature_dependent ? 1.0e-4 : 0.0,
        temperature_dependent ? 2.0e-8 : 0.0,
        2000.0,
        3000.0);
    if (temperature_dependent) {
        fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
        auto functions = std::make_shared<fuelsim::MaterialFunctionSet>(*result.functions);
        functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
            {{"conductivity", 4.0},
                {"density", 2000.0},
                {"specific_heat", 3000.0},
                {"reference_temperature", 300.0},
                {"conductivity_temperature_coefficient", 0.01},
                {"density_temperature_coefficient", 0.0},
                {"specific_heat_temperature_coefficient", 0.0}});
        result.functions = std::move(functions);
    }
    return result;
}

double relative_difference(double numerator_squared, double denominator_squared) {
    return std::sqrt(numerator_squared / denominator_squared);
}

double relative_vector_error(const fuelsim::Hex8LocalResidual& actual,
    const fuelsim::Hex8LocalResidual& expected,
    std::size_t begin) {
    double difference_squared = 0.0, expected_squared = 0.0;
    for (std::size_t row = begin; row < local_size; ++row) {
        difference_squared += (actual[row] - expected[row]) * (actual[row] - expected[row]);
        expected_squared += expected[row] * expected[row];
    }
    return relative_difference(difference_squared, expected_squared);
}

double block_error(const fuelsim::Hex8LocalJacobian& actual,
    const fuelsim::Hex8LocalJacobian& expected,
    std::size_t row_begin,
    std::size_t row_end,
    std::size_t column_begin,
    std::size_t column_end) {
    double difference_squared = 0.0, expected_squared = 0.0;
    for (std::size_t row = row_begin; row < row_end; ++row)
        for (std::size_t column = column_begin; column < column_end; ++column) {
            const double difference = actual[row * local_size + column] - expected[row * local_size + column];
            difference_squared += difference * difference;
            expected_squared += expected[row * local_size + column] * expected[row * local_size + column];
        }
    return relative_difference(difference_squared, expected_squared);
}

fuelsim::Hex8LocalResidual abaqus_c3d8t_residual(const fuelsim::Hex8Geometry& geometry,
    const fuelsim::Hex8LocalValues& state,
    const fuelsim::IsotropicThermoelasticMaterial& material) {
    const fuelsim::CartesianThermoelasticData data{material, 0.0, 0.0};
    fuelsim::Hex8LocalResidual result = fuelsim::compute_hex8_thermoelastic(data, geometry, state);
    for (std::size_t row = 8; row < local_size; ++row)
        result[row] = 0.0;
    fuelsim::Hex8LocalAdValues passive{};
    for (std::size_t dof = 0; dof < local_size; ++dof)
        passive[dof] = state[dof];
    double volume = 0.0, average_trace = 0.0, average_temperature = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        average_temperature += state[node] / 8.0;
    std::array<fuelsim::SymmetricTensor3Values, 8> point_stresses{};
    double element_pressure = 0.0;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const fuelsim::Hex8QuadraturePoint& point = geometry.points[q];
        const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(point,
            passive,
            fuelsim::Hex8LocalValues{},
            fuelsim::StrainFormulation::small);
        volume += point.weighted_measure;
        average_trace += point.weighted_measure
                         * (kinematics.strain_increment.xx.value() + kinematics.strain_increment.yy.value()
                             + kinematics.strain_increment.zz.value());
    }
    average_trace /= volume;
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const fuelsim::Hex8QuadraturePoint& point = geometry.points[q];
        const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(point,
            passive,
            fuelsim::Hex8LocalValues{},
            fuelsim::StrainFormulation::small);
        fuelsim::SymmetricTensor3 selective_strain = kinematics.strain_increment;
        const adlite::Scalar correction =
            (average_trace - selective_strain.xx - selective_strain.yy - selective_strain.zz) / 3.0;
        selective_strain.xx += correction;
        selective_strain.yy += correction;
        selective_strain.zz += correction;
        const adlite::Scalar point_temperature = state[gauss_to_material_node[q]];
        const fuelsim::SymmetricTensor3 point_imposed = material.eigenstrain(point_temperature);
        const fuelsim::SymmetricTensor3 average_imposed = material.eigenstrain(average_temperature);
        const fuelsim::SymmetricTensor3 adjusted_strain{selective_strain.xx + point_imposed.xx - average_imposed.xx,
            selective_strain.yy + point_imposed.yy - average_imposed.yy,
            selective_strain.zz + point_imposed.zz - average_imposed.zz,
            selective_strain.xy + point_imposed.xy - average_imposed.xy,
            selective_strain.yz + point_imposed.yz - average_imposed.yz,
            selective_strain.xz + point_imposed.xz - average_imposed.xz};
        const fuelsim::SymmetricTensor3 stress = material.stress(adjusted_strain, point_temperature);
        point_stresses[q] = {stress.xx.value(),
            stress.yy.value(),
            stress.zz.value(),
            stress.xy.value(),
            stress.yz.value(),
            stress.xz.value()};
        element_pressure +=
            point.weighted_measure / volume * (stress.xx.value() + stress.yy.value() + stress.zz.value()) / 3.0;
    }
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const fuelsim::Hex8QuadraturePoint& point = geometry.points[q];
        fuelsim::SymmetricTensor3Values stress = point_stresses[q];
        const double point_pressure = (stress.xx + stress.yy + stress.zz) / 3.0;
        stress.xx += element_pressure - point_pressure;
        stress.yy += element_pressure - point_pressure;
        stress.zz += element_pressure - point_pressure;
        for (std::size_t node = 0; node < 8; ++node) {
            const double gx = point.gradient[node][0], gy = point.gradient[node][1], gz = point.gradient[node][2];
            result[8 + node] += point.weighted_measure * (stress.xx * gx + stress.xy * gy + stress.xz * gz);
            result[16 + node] += point.weighted_measure * (stress.xy * gx + stress.yy * gy + stress.yz * gz);
            result[24 + node] += point.weighted_measure * (stress.xz * gx + stress.yz * gy + stress.zz * gz);
        }
    }
    return result;
}

bool compare_operator(const std::map<std::string, NodalStep>& steps, bool temperature_dependent) {
    const NodalStep& base = steps.at("BASE");
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
    const fuelsim::IsotropicThermoelasticMaterial material(properties(temperature_dependent));
    const fuelsim::CartesianThermoelasticData data{material, 0.0, 0.0};
    fuelsim::Hex8LocalJacobian fuelsim_jacobian{};
    const fuelsim::Hex8LocalResidual fuelsim_residual =
        fuelsim::compute_hex8_thermoelastic(data, geometry, base.state, nullptr, 0.0, &fuelsim_jacobian);
    fuelsim::Hex8LocalJacobian abaqus_jacobian{};
    for (std::size_t column = 0; column < local_size; ++column) {
        std::ostringstream plus_name, minus_name;
        plus_name << 'D' << (column < 10 ? "0" : "") << column << "_PLUS";
        minus_name << 'D' << (column < 10 ? "0" : "") << column << "_MINUS";
        const double perturbation = column < 8 ? 1.0e-3 : 1.0e-7;
        const NodalStep& plus = steps.at(plus_name.str());
        const NodalStep& minus = steps.at(minus_name.str());
        for (std::size_t dof = 0; dof < local_size; ++dof) {
            const double expected_plus = base.state[dof] + (dof == column ? perturbation : 0.0);
            const double expected_minus = base.state[dof] - (dof == column ? perturbation : 0.0);
            if (std::abs(plus.state[dof] - expected_plus) > 1.0e-12
                || std::abs(minus.state[dof] - expected_minus) > 1.0e-12)
                throw std::invalid_argument("Abaqus C3D8T perturbation state does not match the declared probe");
        }
        for (std::size_t row = 0; row < local_size; ++row)
            abaqus_jacobian[row * local_size + column] =
                (plus.reaction[row] - minus.reaction[row]) / (2.0 * perturbation);
    }

    const double thermal_error = block_error(fuelsim_jacobian, abaqus_jacobian, 0, 8, 0, 8);
    const double mechanical_error = block_error(fuelsim_jacobian, abaqus_jacobian, 8, 32, 8, 32);
    const double coupling_error = block_error(fuelsim_jacobian, abaqus_jacobian, 8, 32, 0, 8);
    fuelsim::Hex8LocalJacobian compatible_jacobian{};
    for (std::size_t column = 0; column < local_size; ++column) {
        const double perturbation = column < 8 ? 1.0e-3 : 1.0e-7;
        fuelsim::Hex8LocalValues plus = base.state, minus = base.state;
        plus[column] += perturbation;
        minus[column] -= perturbation;
        const fuelsim::Hex8LocalResidual plus_residual = abaqus_c3d8t_residual(geometry, plus, material);
        const fuelsim::Hex8LocalResidual minus_residual = abaqus_c3d8t_residual(geometry, minus, material);
        for (std::size_t row = 0; row < local_size; ++row)
            compatible_jacobian[row * local_size + column] =
                (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
    }
    const double compatible_error = block_error(compatible_jacobian, abaqus_jacobian, 0, 32, 0, 32);
    double maximum_thermal_displacement = 0.0;
    for (std::size_t row = 0; row < 8; ++row)
        for (std::size_t column = 8; column < 32; ++column)
            maximum_thermal_displacement = std::max(maximum_thermal_displacement,
                std::max(std::abs(fuelsim_jacobian[row * local_size + column]),
                    std::abs(abaqus_jacobian[row * local_size + column])));
    double coupling_column_difference = 0.0, coupling_scale = 0.0, coupling_sum_difference = 0.0,
           coupling_sum_scale = 0.0;
    for (std::size_t row = 8; row < 32; ++row) {
        double fuelsim_sum = 0.0, abaqus_sum = 0.0;
        for (std::size_t column = 0; column < 8; ++column) {
            const double abaqus = abaqus_jacobian[row * local_size + column];
            coupling_column_difference =
                std::max(coupling_column_difference, std::abs(abaqus - abaqus_jacobian[row * local_size]));
            coupling_scale = std::max(coupling_scale, std::abs(abaqus));
            fuelsim_sum += fuelsim_jacobian[row * local_size + column];
            abaqus_sum += abaqus;
        }
        coupling_sum_difference = std::max(coupling_sum_difference, std::abs(fuelsim_sum - abaqus_sum));
        coupling_sum_scale = std::max(coupling_sum_scale, std::abs(abaqus_sum));
    }
    const fuelsim::Hex8LocalResidual compatible_residual = abaqus_c3d8t_residual(geometry, base.state, material);
    const double thermal_residual_error = relative_vector_error(fuelsim_residual, base.reaction, 0);
    double thermal_difference_squared = 0.0, thermal_scale_squared = 0.0;
    for (std::size_t row = 0; row < 8; ++row) {
        thermal_difference_squared +=
            (fuelsim_residual[row] - base.reaction[row]) * (fuelsim_residual[row] - base.reaction[row]);
        thermal_scale_squared += base.reaction[row] * base.reaction[row];
    }
    const double isolated_thermal_residual_error =
        relative_difference(thermal_difference_squared, thermal_scale_squared);
    const double production_mechanical_residual_error = relative_vector_error(fuelsim_residual, base.reaction, 8);
    const double compatible_mechanical_error = relative_vector_error(compatible_residual, base.reaction, 8);

    std::cout << "b49_ktt_relative_frobenius_error=" << thermal_error << '\n';
    std::cout << "b49_kuu_relative_frobenius_error=" << mechanical_error << '\n';
    std::cout << "b49_kut_relative_frobenius_error=" << coupling_error << '\n';
    std::cout << "b49_c3d8t_compatible_jacobian_relative_frobenius_error=" << compatible_error << '\n';
    std::cout << "b49_ktu_maximum_absolute=" << maximum_thermal_displacement << '\n';
    std::cout << "b49_abaqus_kut_column_uniformity=" << coupling_column_difference / coupling_scale << '\n';
    std::cout << "b49_uniform_temperature_coupling_error=" << coupling_sum_difference / coupling_sum_scale << '\n';
    std::cout << "b49_thermal_residual_relative_error=" << isolated_thermal_residual_error << '\n';
    std::cout << "b49_production_mechanical_residual_relative_error=" << production_mechanical_residual_error << '\n';
    std::cout << "b49_c3d8t_compatible_mechanical_residual_relative_error=" << compatible_mechanical_error << '\n';
    (void)thermal_residual_error;
    bool passed = check(thermal_error < 1.0e-9, "C3D8T and fuelsim thermal-conduction tangent blocks agree")
                  && check(compatible_error < 1.0e-8,
                      "the identified selective-integration C3D8T operator matches the Abaqus tangent")
                  && check(maximum_thermal_displacement < 1.0e-8,
                      "C3D8T and fuelsim thermal residuals have zero displacement derivative")
                  && check(temperature_dependent || coupling_column_difference / coupling_scale < 1.0e-9,
                      "constant-property Abaqus first-order coupled element uses one constant expansion temperature")
                  && check(coupling_sum_difference / coupling_sum_scale < 1.0e-9,
                      "uniform temperature perturbations give the same mechanical coupling")
                  && check(isolated_thermal_residual_error < 1.0e-12, "C3D8T and fuelsim base thermal residuals agree")
                  && check(compatible_mechanical_error < 1.0e-12,
                      "C3D8T reactions match selective mechanical integration and constant expansion temperature")
                  && check(coupling_error < 1.0e-8,
                      "fuelsim and C3D8T temperature-to-mechanics tangent blocks use the same element temperature")
                  && check(mechanical_error < 1.0e-8 && production_mechanical_residual_error < 1.0e-12,
                      "the production selective mechanical operator matches C3D8T")
                  && true;
    return passed;
}

double tensor_maximum_difference(const fuelsim::SymmetricTensor3Values& first,
    const fuelsim::SymmetricTensor3Values& second) {
    return std::max({std::abs(first.xx - second.xx),
        std::abs(first.yy - second.yy),
        std::abs(first.zz - second.zz),
        std::abs(first.xy - second.xy),
        std::abs(first.yz - second.yz),
        std::abs(first.xz - second.xz)});
}

bool compare_integration_points(const NodalStep& base,
    const std::vector<IntegrationPointReference>& references,
    bool temperature_dependent) {
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
    const fuelsim::IsotropicThermoelasticMaterial material(properties(temperature_dependent));
    const fuelsim::CartesianThermoelasticData data{material, 0.0, 0.0};
    const std::array<fuelsim::SymmetricTensor3Values, 8> production_stresses =
        fuelsim::compute_hex8_stress(data, geometry, base.state);
    fuelsim::Hex8LocalAdValues passive{};
    for (std::size_t dof = 0; dof < local_size; ++dof)
        passive[dof] = base.state[dof];
    double average_temperature = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        average_temperature += base.state[node] / 8.0;
    double volume = 0.0, average_trace = 0.0;
    for (const fuelsim::Hex8QuadraturePoint& point : geometry.points) {
        const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(point,
            passive,
            fuelsim::Hex8LocalValues{},
            fuelsim::StrainFormulation::small);
        volume += point.weighted_measure;
        average_trace += point.weighted_measure
                         * (kinematics.strain_increment.xx.value() + kinematics.strain_increment.yy.value()
                             + kinematics.strain_increment.zz.value());
    }
    average_trace /= volume;
    double maximum_coordinate_difference = 0.0, maximum_strain_difference = 0.0, maximum_heat_flux_difference = 0.0,
           maximum_constant_temperature_stress_difference = 0.0, maximum_production_stress_difference = 0.0;
    std::array<bool, 8> used{};
    for (const IntegrationPointReference& reference : references) {
        std::size_t closest = 0;
        double closest_squared = std::numeric_limits<double>::max();
        for (std::size_t q = 0; q < geometry.points.size(); ++q) {
            const fuelsim::CartesianPoint3& point = geometry.points[q].position;
            const double distance_squared = (point.x - reference.point.x) * (point.x - reference.point.x)
                                            + (point.y - reference.point.y) * (point.y - reference.point.y)
                                            + (point.z - reference.point.z) * (point.z - reference.point.z);
            if (distance_squared < closest_squared) {
                closest = q;
                closest_squared = distance_squared;
            }
        }
        if (used[closest])
            throw std::invalid_argument("Abaqus integration-point coordinate mapping is not unique");
        used[closest] = true;
        maximum_coordinate_difference = std::max(maximum_coordinate_difference, std::sqrt(closest_squared));
        const fuelsim::Hex8QuadraturePoint& point = geometry.points[closest];
        const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(point,
            passive,
            fuelsim::Hex8LocalValues{},
            fuelsim::StrainFormulation::small);
        fuelsim::SymmetricTensor3 selective_strain = kinematics.strain_increment;
        const adlite::Scalar correction =
            (average_trace - selective_strain.xx - selective_strain.yy - selective_strain.zz) / 3.0;
        selective_strain.xx += correction;
        selective_strain.yy += correction;
        selective_strain.zz += correction;
        const fuelsim::SymmetricTensor3Values strain = {selective_strain.xx.value(),
            selective_strain.yy.value(),
            selective_strain.zz.value(),
            selective_strain.xy.value(),
            selective_strain.yz.value(),
            selective_strain.xz.value()};
        maximum_strain_difference =
            std::max(maximum_strain_difference, tensor_maximum_difference(strain, reference.strain));
        std::array<double, 3> gradient_temperature{};
        const adlite::Scalar point_temperature = base.state[gauss_to_material_node[closest]];
        for (std::size_t node = 0; node < 8; ++node)
            for (std::size_t direction = 0; direction < 3; ++direction)
                gradient_temperature[direction] += point.gradient[node][direction] * base.state[node];
        const double conductivity = material.conductivity(point_temperature).value();
        for (std::size_t direction = 0; direction < 3; ++direction)
            maximum_heat_flux_difference = std::max(maximum_heat_flux_difference,
                std::abs(-conductivity * gradient_temperature[direction] - reference.heat_flux[direction]));
        const fuelsim::SymmetricTensor3 point_imposed = material.eigenstrain(point_temperature);
        const fuelsim::SymmetricTensor3 average_imposed = material.eigenstrain(average_temperature);
        const fuelsim::SymmetricTensor3 adjusted_strain{selective_strain.xx + point_imposed.xx - average_imposed.xx,
            selective_strain.yy + point_imposed.yy - average_imposed.yy,
            selective_strain.zz + point_imposed.zz - average_imposed.zz,
            selective_strain.xy + point_imposed.xy - average_imposed.xy,
            selective_strain.yz + point_imposed.yz - average_imposed.yz,
            selective_strain.xz + point_imposed.xz - average_imposed.xz};
        const fuelsim::SymmetricTensor3 abaqus_temperature_stress = material.stress(adjusted_strain, point_temperature);
        const fuelsim::SymmetricTensor3Values abaqus_temperature_values = {abaqus_temperature_stress.xx.value(),
            abaqus_temperature_stress.yy.value(),
            abaqus_temperature_stress.zz.value(),
            abaqus_temperature_stress.xy.value(),
            abaqus_temperature_stress.yz.value(),
            abaqus_temperature_stress.xz.value()};
        maximum_constant_temperature_stress_difference = std::max(maximum_constant_temperature_stress_difference,
            tensor_maximum_difference(abaqus_temperature_values, reference.stress));
        maximum_production_stress_difference = std::max(maximum_production_stress_difference,
            tensor_maximum_difference(production_stresses[closest], reference.stress));
        (void)reference.output_temperature;
    }
    std::cout << "b49_integration_point_coordinate_maximum_difference=" << maximum_coordinate_difference << '\n';
    std::cout << "b49_integration_point_strain_maximum_difference=" << maximum_strain_difference << '\n';
    std::cout << "b49_integration_point_heat_flux_maximum_difference=" << maximum_heat_flux_difference << '\n';
    std::cout << "b49_constant_temperature_stress_maximum_difference=" << maximum_constant_temperature_stress_difference
              << '\n';
    std::cout << "b49_production_stress_maximum_difference=" << maximum_production_stress_difference << '\n';
    return check(maximum_coordinate_difference < 1.0e-12, "C3D8T and fuelsim integration-point coordinates agree")
           && check(maximum_strain_difference < 1.0e-15,
               "C3D8T and fuelsim non-affine strains agree at all eight points")
           && check(maximum_heat_flux_difference < 1.0e-10, "C3D8T and fuelsim heat fluxes agree at all eight points")
           && check(maximum_constant_temperature_stress_difference < 1.0e-3,
               "C3D8T stresses use the element-average expansion temperature")
           && check(maximum_production_stress_difference < 1.0e-3,
               "fuelsim production selective-integration stresses match C3D8T")
           && true;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: fuelsim_b49_hex8_c3d8t_operator_abaqus_tests <nodal.csv> <integration_points.csv> "
                     "[temperature_dependent]\n";
        return 2;
    }
    try {
        std::cout << std::scientific;
        const bool temperature_dependent = argc == 4 && std::string(argv[3]) == "temperature_dependent";
        if (argc == 4 && !temperature_dependent)
            throw std::invalid_argument("Unknown C3D8T operator material mode");
        const std::map<std::string, NodalStep> steps = read_nodal_steps(argv[1]);
        const std::vector<IntegrationPointReference> points = read_integration_points(argv[2]);
        bool passed = compare_operator(steps, temperature_dependent);
        passed = compare_integration_points(steps.at("BASE"), points, temperature_dependent) && passed;
        if (!passed)
            return 1;
        std::cout << "[PASS] Abaqus C3D8T operator identification\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] Abaqus C3D8T operator identification raised: " << error.what() << '\n';
        return 1;
    }
}
