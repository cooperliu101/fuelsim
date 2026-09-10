#include "boundary_types.hpp"
#include "c3d20_types.hpp"
#include "contact_types.hpp"
#include "quad8_face_boundary.hpp"
#include "quad8_face_contact.hpp"
#include "support/element_evaluation.hpp"

#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool near(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance * std::max({1.0, std::abs(actual), std::abs(expected)});
}

fuelsim::Hex20Coordinates unit_cube() {
    return {{{0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 1.0, 1.0},
        {0.5, 0.0, 0.0},
        {1.0, 0.5, 0.0},
        {0.5, 1.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.0, 0.0, 0.5},
        {1.0, 0.0, 0.5},
        {1.0, 1.0, 0.5},
        {0.0, 1.0, 0.5},
        {0.5, 0.0, 1.0},
        {1.0, 0.5, 1.0},
        {0.5, 1.0, 1.0},
        {0.0, 0.5, 1.0}}};
}

fuelsim::ThermoelasticProperties material(bool inelastic = false) {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e5, 0.25, 1.2e-5, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
    if (inelastic) {
        result = fuelsim::test::with_norton(std::move(result), 1.0e-6, 10.0, 3.0, 300.0);
        result = fuelsim::test::with_plasticity(std::move(result), 20.0, 10.0, 300.0);
    }
    return result;
}

double directional_jacobian_error(const fuelsim::CartesianThermoelasticData& data,
    const fuelsim::Hex20Geometry& geometry,
    const fuelsim::Hex20LocalValues& state,
    const fuelsim::Hex20LocalValues& old,
    const fuelsim::CartesianMaterialHistory& history,
    std::size_t first_row = 0,
    std::size_t row_count = fuelsim::hex20_local_dof_count) {
    fuelsim::Hex20LocalJacobian jacobian{};
    (void)fuelsim::compute_c3d20_transient(data, geometry, state, old, history, 0.5, &jacobian);
    fuelsim::Hex20LocalValues direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.29 * static_cast<double>(dof + 1));
    constexpr double step = 2.0e-7;
    auto plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const auto plus_residual = fuelsim::compute_c3d20_transient(data, geometry, plus, old, history, 0.5);
    const auto minus_residual = fuelsim::compute_c3d20_transient(data, geometry, minus, old, history, 0.5);
    double error = 0.0, scale = 0.0;
    for (std::size_t row = first_row; row < first_row + row_count; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            analytic += jacobian[row * state.size() + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        error = std::max(error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    return error / scale;
}

double residual_path_error(const fuelsim::CartesianThermoelasticData& data,
    const fuelsim::Hex20Geometry& geometry,
    const fuelsim::Hex20LocalValues& state,
    const fuelsim::Hex20LocalValues& old,
    const fuelsim::CartesianMaterialHistory& history) {
    fuelsim::Hex20LocalJacobian jacobian{};
    const auto system_residual = fuelsim::compute_c3d20_transient(data, geometry, state, old, history, 0.5, &jacobian);
    const auto residual = fuelsim::compute_c3d20_transient(data, geometry, state, old, history, 0.5);
    double error = 0.0, scale = 0.0;
    for (std::size_t row = 0; row < residual.size(); ++row) {
        error = std::max(error, std::abs(residual[row] - system_residual[row]));
        scale = std::max({scale, 1.0, std::abs(residual[row]), std::abs(system_residual[row])});
    }
    return error / scale;
}

double material_state_error(const fuelsim::CartesianMaterialPointState& first,
    const fuelsim::CartesianMaterialPointState& second) {
    double error = 0.0, scale = 1.0;
    const std::array<const std::array<double, 6>*, 3> first_histories = {&first.elastic_strain,
        &first.plastic_strain,
        &first.creep_strain};
    const std::array<const std::array<double, 6>*, 3> second_histories = {&second.elastic_strain,
        &second.plastic_strain,
        &second.creep_strain};
    for (std::size_t history = 0; history < first_histories.size(); ++history)
        for (std::size_t component = 0; component < 6; ++component) {
            const double first_value = (*first_histories[history])[component];
            const double second_value = (*second_histories[history])[component];
            error = std::max(error, std::abs(first_value - second_value));
            scale = std::max({scale, std::abs(first_value), std::abs(second_value)});
        }
    const std::array<double, 8> first_scalars = {first.stress.xx,
        first.stress.yy,
        first.stress.zz,
        first.stress.xy,
        first.stress.yz,
        first.stress.xz,
        first.equivalent_plastic_strain,
        first.equivalent_creep_strain};
    const std::array<double, 8> second_scalars = {second.stress.xx,
        second.stress.yy,
        second.stress.zz,
        second.stress.xy,
        second.stress.yz,
        second.stress.xz,
        second.equivalent_plastic_strain,
        second.equivalent_creep_strain};
    for (std::size_t component = 0; component < first_scalars.size(); ++component) {
        error = std::max(error, std::abs(first_scalars[component] - second_scalars[component]));
        scale = std::max({scale, std::abs(first_scalars[component]), std::abs(second_scalars[component])});
    }
    return error / scale;
}

bool test_material_value_paths() {
    fuelsim::ThermoelasticProperties elastic =
        fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e5, 0.25, 1.2e-5, 300.0, -10.0, 1.0e-5, 2.0e-8);
    std::array<fuelsim::ThermoelasticProperties, 4> properties = {elastic,
        fuelsim::test::with_plasticity(elastic, 20.0, 10.0, 300.0, -0.01, 0.02),
        fuelsim::test::with_norton(elastic, 1.0e-6, 10.0, 3.0, 300.0, 1.0e-8, 0.01, 1.0e-3),
        fuelsim::test::with_plasticity(
            fuelsim::test::with_norton(elastic, 1.0e-6, 10.0, 3.0, 300.0, 1.0e-8, 0.01, 1.0e-3),
            20.0,
            10.0,
            300.0,
            -0.01,
            0.02)};
    const fuelsim::SymmetricTensor3Values strain{0.02, -0.004, 0.002, 0.003, -0.001, 0.002};
    const fuelsim::SymmetricTensor3 active_strain{strain.xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz};
    fuelsim::CartesianRotation rotation;
    constexpr double angle = 0.31;
    rotation.xx = std::cos(angle);
    rotation.xy = -std::sin(angle);
    rotation.yx = std::sin(angle);
    rotation.yy = std::cos(angle);
    double maximum_error = 0.0;
    for (const fuelsim::ThermoelasticProperties& property : properties) {
        const fuelsim::IsotropicThermoelasticMaterial model(property);
        const fuelsim::CartesianMaterialPointState committed{};
        const fuelsim::CartesianMaterialPointState active =
            model.response(active_strain, adlite::Scalar(315.0), 0.5, committed).trial_state;
        const fuelsim::CartesianMaterialPointState values = model.response_values(strain, 315.0, 0.5, committed);
        maximum_error = std::max(maximum_error, material_state_error(active, values));
        const fuelsim::CartesianMaterialPointState active_incremental =
            model.incremental_response(active_strain, rotation, adlite::Scalar(325.0), 315.0, 0.5, active).trial_state;
        const fuelsim::CartesianMaterialPointState values_incremental =
            model.incremental_response_values(strain, rotation, 325.0, 315.0, 0.5, active);
        maximum_error = std::max(maximum_error, material_state_error(active_incremental, values_incremental));
    }
    std::cout << "hex20_material_value_path_relative_error=" << maximum_error << '\n';
    return check(maximum_error < 2.0e-14,
        "ordinary-double Cartesian material values match passive automatic differentiation for all built-in branches");
}

bool test_geometry_and_constant_strain() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::make_c3d20_geometry(coordinates);
    double thermal_volume = 0.0, mechanical_volume = 0.0;
    for (const fuelsim::Hex20ThermalQuadraturePoint& point : geometry.thermal_points) {
        thermal_volume += point.weighted_measure;
        double temperature_sum = 0.0;
        std::array<double, 3> temperature_gradient{};
        for (std::size_t node = 0; node < 8; ++node) {
            temperature_sum += point.temperature_shape[node];
            for (std::size_t direction = 0; direction < 3; ++direction)
                temperature_gradient[direction] += point.temperature_gradient[node][direction];
        }
        if (!check(near(temperature_sum, 1.0, 2.0e-14), "HEX20 thermal shape functions form a partition of unity")
            || !check(std::max({std::abs(temperature_gradient[0]),
                          std::abs(temperature_gradient[1]),
                          std::abs(temperature_gradient[2])})
                          < 2.0e-14,
                "HEX20 thermal shape gradients sum to zero"))
            return false;
    }
    for (const fuelsim::Hex20MechanicalQuadraturePoint& point : geometry.mechanical_points) {
        mechanical_volume += point.weighted_measure;
        double temperature_sum = 0.0, displacement_sum = 0.0;
        std::array<double, 3> displacement_gradient{};
        for (std::size_t node = 0; node < 8; ++node)
            temperature_sum += point.temperature_shape[node];
        for (std::size_t node = 0; node < 20; ++node) {
            displacement_sum += point.displacement_shape[node];
            for (std::size_t direction = 0; direction < 3; ++direction)
                displacement_gradient[direction] += point.displacement_gradient[node][direction];
        }
        if (!check(near(temperature_sum, 1.0, 2.0e-14) && near(displacement_sum, 1.0, 2.0e-14),
                "HEX20 mechanical points carry separate U2 and T1 partitions of unity")
            || !check(std::max({std::abs(displacement_gradient[0]),
                          std::abs(displacement_gradient[1]),
                          std::abs(displacement_gradient[2])})
                          < 2.0e-14,
                "HEX20 mechanical displacement shape gradients sum to zero"))
            return false;
    }
    if (!check(near(thermal_volume, 1.0, 2.0e-14), "HEX20 thermal 8-point integration recovers unit volume")
        || !check(near(mechanical_volume, 1.0, 2.0e-14), "HEX20 mechanical 27-point integration recovers unit volume"))
        return false;
    fuelsim::Hex20LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node)
        state[node] = 300.0;
    for (std::size_t node = 0; node < 20; ++node) {
        const auto& point = coordinates[node];
        state[8 + node] = 0.01 * point.x + 0.004 * point.y + 0.008 * point.z;
        state[28 + node] = 0.004 * point.x - 0.02 * point.y - 0.006 * point.z;
        state[48 + node] = 0.008 * point.x - 0.006 * point.y + 0.03 * point.z;
    }
    const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(material()), 0.0, 0.0};
    const auto stresses = fuelsim::compute_c3d20_stress(data, geometry, state);
    const double lambda = 2.0e5 * 0.25 / (1.25 * 0.5), shear = 2.0e5 / 2.5, trace = 0.02;
    for (const auto& stress : stresses)
        if (!check(near(stress.xx, lambda * trace + 2.0 * shear * 0.01, 3.0e-13)
                       && near(stress.yy, lambda * trace - 2.0 * shear * 0.02, 3.0e-13)
                       && near(stress.zz, lambda * trace + 2.0 * shear * 0.03, 3.0e-13)
                       && near(stress.xy, 2.0 * shear * 0.004, 3.0e-13)
                       && near(stress.yz, -2.0 * shear * 0.006, 3.0e-13)
                       && near(stress.xz, 2.0 * shear * 0.008, 3.0e-13),
                "quadratic displacement interpolation reproduces all constant-strain stresses"))
            return false;
    return true;
}

bool test_finite_thermal_operators() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::make_c3d20_geometry(coordinates);
    const fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e5, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2.0, 3.0);
    const fuelsim::CartesianThermoelasticData conduction_data{fuelsim::IsotropicThermoelasticMaterial(properties),
        0.0,
        1.0,
        fuelsim::StrainFormulation::finite};
    fuelsim::Hex20LocalValues state{};
    constexpr double gradient_x = 7.0, gradient_y = -5.0, gradient_z = 3.0;
    for (std::size_t node = 0; node < 8; ++node)
        state[node] = 300.0 + gradient_x * coordinates[node].x + gradient_y * coordinates[node].y
                      + gradient_z * coordinates[node].z;
    for (std::size_t node = 0; node < 20; ++node) {
        state[8 + node] = coordinates[node].x;
        state[28 + node] = 0.2 * coordinates[node].y;
    }
    const auto conduction = fuelsim::compute_c3d20_thermoelastic(conduction_data, geometry, state);
    bool passed = true;
    for (std::size_t node = 0; node < 8; ++node) {
        const double sign_x = 2.0 * coordinates[node].x - 1.0, sign_y = 2.0 * coordinates[node].y - 1.0,
                     sign_z = 2.0 * coordinates[node].z - 1.0;
        const double expected =
            4.0 * 2.4 * 0.25
            * (sign_x * gradient_x / (1.5 * 1.5) + sign_y * gradient_y / (1.1 * 1.1) + sign_z * gradient_z);
        passed = check(near(conduction[node], expected, 2.0e-13),
                     "finite-strain HEX20 conduction uses midpoint gradients and the current quadratic volume")
                 && passed;
    }

    fuelsim::Hex20LocalValues uniform = state;
    for (std::size_t node = 0; node < 8; ++node)
        uniform[node] = 300.0;
    const fuelsim::CartesianThermoelasticData source_data{fuelsim::IsotropicThermoelasticMaterial(properties),
        10.0,
        1.0,
        fuelsim::StrainFormulation::finite};
    const auto source = fuelsim::compute_c3d20_thermoelastic(source_data, geometry, uniform);
    const double source_sum = std::accumulate(source.begin(), source.begin() + 8, 0.0);
    passed =
        check(near(source_sum, -24.0, 2.0e-13), "finite-strain HEX20 body source uses the current eight-corner volume")
        && passed;

    fuelsim::Hex20LocalValues midside = uniform;
    for (std::size_t node = 0; node < 20; ++node) {
        midside[8 + node] = 0.0;
        midside[28 + node] = 0.0;
    }
    midside[8 + 8] = 0.05;
    const auto midside_source = fuelsim::compute_c3d20_thermoelastic(source_data, geometry, midside);
    const double midside_source_sum = std::accumulate(midside_source.begin(), midside_source.begin() + 8, 0.0);
    passed = check(near(midside_source_sum, -10.0, 2.0e-13),
                 "finite-strain HEX20 body source excludes displacement midside nodes from its volume")
             && passed;
    const fuelsim::CartesianThermoelasticData small_source_data{fuelsim::IsotropicThermoelasticMaterial(properties),
        10.0,
        1.0,
        fuelsim::StrainFormulation::small};
    const auto small_source = fuelsim::compute_c3d20_thermoelastic(small_source_data, geometry, midside);
    passed = check(near(std::accumulate(small_source.begin(), small_source.begin() + 8, 0.0), -10.0, 2.0e-13),
                 "small-strain HEX20 body source retains the reference-volume operator")
             && passed;

    fuelsim::Hex20LocalValues old = state;
    for (std::size_t node = 0; node < 8; ++node) {
        old[node] = 300.0;
        state[node] = 301.0 + static_cast<double>(node);
    }
    for (std::size_t node = 0; node < 20; ++node) {
        old[8 + node] = 0.0;
        old[28 + node] = 0.0;
        old[48 + node] = 0.0;
    }
    const fuelsim::CartesianMaterialHistory history(27);
    const auto transient =
        fuelsim::compute_c3d20_transient(conduction_data, geometry, state, old, history, 1.0, nullptr, true);
    const auto steady =
        fuelsim::compute_c3d20_transient(conduction_data, geometry, state, old, history, 1.0, nullptr, false);
    for (std::size_t node = 0; node < 8; ++node) {
        double expected_capacity = 0.0;
        for (std::size_t other = 0; other < 8; ++other) {
            double mass = 1.0;
            for (std::size_t direction = 0; direction < 3; ++direction) {
                const double left = direction == 0   ? coordinates[node].x
                                    : direction == 1 ? coordinates[node].y
                                                     : coordinates[node].z;
                const double right = direction == 0   ? coordinates[other].x
                                     : direction == 1 ? coordinates[other].y
                                                      : coordinates[other].z;
                mass *= left == right ? 2.0 : 1.0;
            }
            expected_capacity += mass * static_cast<double>(other + 1);
        }
        expected_capacity *= 6.0 * 2.4 / 216.0;
        passed = check(near(transient[node] - steady[node], expected_capacity, 3.0e-12),
                     "finite-strain HEX20 capacity uses the consistent current-volume matrix")
                 && passed;
    }

    fuelsim::Hex20LocalValues invalid = uniform;
    for (std::size_t node = 0; node < 20; ++node) {
        invalid[8 + node] = -2.0 * coordinates[node].x;
        invalid[28 + node] = -2.0 * coordinates[node].y;
    }
    bool midpoint_rejected = false;
    try {
        (void)fuelsim::compute_c3d20_thermoelastic(conduction_data, geometry, invalid);
    } catch (const std::domain_error& error) {
        midpoint_rejected = std::string(error.what()).find("midpoint") != std::string::npos;
    }
    return check(midpoint_rejected, "finite-strain HEX20 thermal integration rejects a singular midpoint geometry")
           && passed;
}

bool test_jacobian_and_transient_history() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::make_c3d20_geometry(coordinates);
    const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(material(true)),
        4.0e5,
        1.0,
        fuelsim::StrainFormulation::small};
    fuelsim::Hex20LocalValues old{}, state{};
    for (std::size_t node = 0; node < 8; ++node)
        old[node] = state[node] = 300.0 + 2.0 * static_cast<double>(node);
    for (std::size_t node = 0; node < 20; ++node) {
        const auto& point = coordinates[node];
        state[8 + node] = 0.02 * point.x + 0.003 * point.y;
        state[28 + node] = 0.003 * point.x - 0.004 * point.y;
        state[48 + node] = 0.002 * point.z;
    }
    const fuelsim::CartesianMaterialHistory history(27);
    const double small_error = directional_jacobian_error(data, geometry, state, old, history);
    const double small_residual_error = residual_path_error(data, geometry, state, old, history);
    const fuelsim::CartesianMaterialHistory update =
        fuelsim::compute_c3d20_transient_update(data, geometry, state, old, history, 0.5);
    bool active = false;
    for (const auto& point : update)
        active = active || point.equivalent_plastic_strain > 0.0 || point.equivalent_creep_strain > 0.0;
    const fuelsim::CartesianThermoelasticData finite_data{fuelsim::IsotropicThermoelasticMaterial(material(true)),
        4.0e5,
        1.0,
        fuelsim::StrainFormulation::finite};
    const double finite_error = directional_jacobian_error(finite_data, geometry, state, old, history);
    const double finite_thermal_error = directional_jacobian_error(finite_data, geometry, state, old, history, 0, 8);
    const double finite_residual_error = residual_path_error(finite_data, geometry, state, old, history);
    fuelsim::Hex20LocalValues invalid = state;
    for (std::size_t node = 0; node < 20; ++node)
        invalid[8 + node] = -2.0 * coordinates[node].x;
    bool invalid_rejected = false;
    try {
        fuelsim::validate_hex20_deformation(geometry.mechanical_points[0], invalid);
    } catch (const std::domain_error&) {
        invalid_rejected = true;
    }
    std::cout << "hex20_directional_jacobian_relative_error=" << small_error << '\n'
              << "hex20_finite_directional_jacobian_relative_error=" << finite_error << '\n'
              << "hex20_finite_thermal_directional_jacobian_relative_error=" << finite_thermal_error << '\n'
              << "hex20_residual_path_relative_error=" << small_residual_error << '\n'
              << "hex20_finite_residual_path_relative_error=" << finite_residual_error << '\n';
    return check(small_error < 2.0e-6,
               "68-DOF narrow automatic-differentiation Jacobian matches a centered directional difference")
           && check(finite_error < 3.0e-6, "finite-strain 68-DOF Jacobian matches a centered directional difference")
           && check(finite_thermal_error < 3.0e-6,
               "finite-strain HEX20 thermal rows match a centered directional difference")
           && check(small_residual_error < 2.0e-14 && finite_residual_error < 2.0e-14,
               "ordinary-double HEX20 residual matches the Jacobian-call residual")
           && check(invalid_rejected, "finite-strain HEX20 rejects a nonpositive deformation Jacobian")
           && check(update.size() == 27 && active,
               "HEX20 transient update commits 27 active inelastic material points");
}

bool test_quadratic_face() {
    const auto cube = unit_cube();
    const fuelsim::Quad8FaceCoordinates coordinates =
        {cube[1], cube[2], cube[6], cube[5], cube[9], cube[14], cube[17], cube[13]};
    const fuelsim::Quad8FaceGeometry geometry = fuelsim::make_quad8_face_geometry(coordinates);
    double thermal_area = 0.0, mechanical_area = 0.0;
    for (const auto& point : geometry.thermal_points)
        thermal_area += point.weighted_measure;
    for (const auto& point : geometry.mechanical_points) {
        const auto& tangent_xi = point.tangent_xi;
        const auto& tangent_eta = point.tangent_eta;
        const double area_x = tangent_xi.y * tangent_eta.z - tangent_xi.z * tangent_eta.y;
        const double area_y = tangent_xi.z * tangent_eta.x - tangent_xi.x * tangent_eta.z;
        const double area_z = tangent_xi.x * tangent_eta.y - tangent_xi.y * tangent_eta.x;
        mechanical_area += std::sqrt(area_x * area_x + area_y * area_y + area_z * area_z) * point.quadrature_weight;
    }
    if (!check(geometry.thermal_points.size() == 4 && geometry.mechanical_points.size() == 9,
            "Quad8 boundary separates 2x2 thermal and 3x3 mechanical integration rules")
        || !check(near(thermal_area, 1.0, 2.0e-14) && near(mechanical_area, 1.0, 2.0e-14),
            "Quad8 thermal and mechanical face rules recover the unit face area"))
        return false;
    fuelsim::Quad8FaceLocalValues state{};
    for (std::size_t node = 0; node < 4; ++node)
        state[node] = 350.0;
    const fuelsim::Quad4FaceBoundaryData pressure = {fuelsim::Quad4FaceBoundaryKind::pressure,
        fuelsim::CartesianTractionComponent::x,
        5.0,
        0.0};
    const auto pressure_residual = fuelsim::compute_quad8_face_boundary(pressure, geometry, state);
    double force_x = 0.0, force_y = 0.0, force_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        force_x += pressure_residual[4 + node];
        force_y += pressure_residual[12 + node];
        force_z += pressure_residual[20 + node];
    }
    const fuelsim::Quad4FaceBoundaryData convection = {fuelsim::Quad4FaceBoundaryKind::convection,
        fuelsim::CartesianTractionComponent::x,
        20.0,
        300.0};
    fuelsim::Quad8FaceLocalJacobian tangent{};
    const auto heat_residual = fuelsim::compute_quad8_face_boundary(convection, geometry, state, &tangent);
    double heat = 0.0, tangent_sum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        heat += heat_residual[row];
        for (std::size_t column = 0; column < 4; ++column)
            tangent_sum += tangent[row * fuelsim::quad8_face_local_dof_count + column];
    }
    return check(near(force_x, 5.0, 2.0e-14) && near(force_y, 0.0, 2.0e-14) && near(force_z, 0.0, 2.0e-14),
               "quadratic face pressure recovers the exact resultant")
           && check(near(heat, 1000.0, 2.0e-14) && near(tangent_sum, 20.0, 2.0e-14),
               "linear face temperature convection recovers exact heat rate and tangent");
}

bool test_warped_geometry() {
    auto coordinates = unit_cube();
    coordinates[8].x = 0.56;
    coordinates[9].y = 0.46;
    coordinates[14].z = 0.54;
    const fuelsim::Hex20Geometry geometry = fuelsim::make_c3d20_geometry(coordinates);
    double thermal_volume = 0.0;
    double mechanical_volume = 0.0;
    for (const auto& point : geometry.thermal_points)
        thermal_volume += point.weighted_measure;
    for (const auto& point : geometry.mechanical_points)
        mechanical_volume += point.weighted_measure;
    return check(thermal_volume > 0.0 && mechanical_volume > 0.0,
               "warped HEX20 geometry has positive thermal and mechanical measures")
           && check(std::isfinite(thermal_volume) && std::isfinite(mechanical_volume),
               "warped HEX20 geometry measures remain finite");
}

bool test_hex20_heat_patch() {
    const auto cube = unit_cube();
    const fuelsim::Quad8FaceCoordinates face =
        {cube[1], cube[2], cube[6], cube[5], cube[9], cube[14], cube[17], cube[13]};
    const auto geometry = fuelsim::make_quad8_face_geometry(face);
    std::vector<fuelsim::Quad8HeatPatchSample> samples;
    std::vector<double> state(84, 0.0), coefficients(84, 0.0);
    double area = 0.0;
    for (std::size_t sample_index = 0; sample_index < 2; ++sample_index) {
        const auto& point = geometry.thermal_points[sample_index];
        fuelsim::Quad8HeatPatchSample sample{{face,
                                                 face,
                                                 point.temperature_shape,
                                                 point.displacement_shape,
                                                 point.derivative_xi,
                                                 point.derivative_eta,
                                                 point.quadrature_weight,
                                                 -1.0},
            {}};
        for (std::size_t entry = 0; entry < sample.local_dofs.size(); ++entry)
            sample.local_dofs[entry] = entry;
        if (sample_index == 1) {
            for (std::size_t node = 0; node < 4; ++node)
                sample.local_dofs[node] = 56 + node;
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t node = 0; node < 8; ++node)
                    sample.local_dofs[8 + component * 16 + node] = 60 + component * 8 + node;
        }
        for (std::size_t node = 0; node < 4; ++node) {
            state[sample.local_dofs[node]] = 350.0 + 20.0 * static_cast<double>(node + sample_index);
            state[sample.local_dofs[4 + node]] = 300.0 + 3.0 * static_cast<double>(node);
            coefficients[sample.local_dofs[node]] += point.weighted_measure * point.temperature_shape[node];
            coefficients[sample.local_dofs[4 + node]] -= point.weighted_measure * point.temperature_shape[node];
        }
        area += point.weighted_measure;
        samples.push_back(sample);
    }
    fuelsim::GapHeatProperties properties{2.0, 0.1};
    std::vector<double> jacobian;
    const auto result = fuelsim::compute_quad8_gap_heat_patch(properties, samples, state, &jacobian);
    double jump_integral = 0.0, analytic_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row)
        jump_integral += coefficients[row] * state[row];
    for (std::size_t row = 0; row < state.size(); ++row)
        analytic_error =
            std::max(analytic_error, std::abs(result[row] - coefficients[row] * 20.0 * jump_integral / area));
    bool passed = check(analytic_error < 1e-10,
        "HEX20 patch averages the temperature difference before transferring heat across shared nodes");
    for (std::size_t entry = 8; entry < 56; ++entry)
        state[entry] = 0.0002 * std::sin(static_cast<double>(entry));
    for (std::size_t entry = 60; entry < 84; ++entry)
        state[entry] = 0.0003 * std::cos(static_cast<double>(entry));
    for (const auto& sample : samples)
        for (std::size_t node = 0; node < 8; ++node)
            state[sample.local_dofs[8 + node]] -= 0.01;
    properties.law = fuelsim::GapHeatConductanceLaw::affine;
    properties.conductance = 10.0;
    properties.clearance_derivative = 2.0;
    properties.temperature_derivative = 0.03;
    properties.reference_temperature = 300.0;
    const auto active = fuelsim::compute_quad8_gap_heat_patch(properties, samples, state, &jacobian);
    const auto passive = fuelsim::compute_quad8_gap_heat_patch(properties, samples, state);
    passed =
        check(active == passive, "HEX20 averaged heat residual is identical with and without its Jacobian") && passed;
    auto plus = state, minus = state;
    std::vector<double> direction(state.size());
    constexpr double step = 1e-6;
    for (std::size_t column = 0; column < state.size(); ++column) {
        direction[column] = 0.1 * std::sin(static_cast<double>(column + 1));
        plus[column] += step * direction[column];
        minus[column] -= step * direction[column];
    }
    const auto forward = fuelsim::compute_quad8_gap_heat_patch(properties, samples, plus);
    const auto backward = fuelsim::compute_quad8_gap_heat_patch(properties, samples, minus);
    double error = 0.0, scale = 0.0, balance = 0.0, column_balance = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0, sum = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column) {
            tangent += jacobian[row * state.size() + column] * direction[column];
            sum += jacobian[column * state.size() + row];
        }
        const double numerical = (forward[row] - backward[row]) / (2.0 * step);
        error = std::max(error, std::abs(tangent - numerical));
        scale = std::max(scale, std::abs(tangent));
        column_balance = std::max(column_balance, std::abs(sum));
        balance += active[row];
    }
    std::cout << "hex20_heat_patch_directional_error=" << error / scale << '\n';
    return check(error / scale < 2e-6, "HEX20 averaged heat geometric and temperature tangent matches differences")
           && check(std::abs(balance) < 1e-10 && column_balance < 1e-10,
               "HEX20 averaged heat residual and each tangent column conserve heat")
           && passed;
}

bool test_hex20_contact_kernels() {
    const auto cube = unit_cube();
    const fuelsim::Quad8FaceCoordinates face =
        {cube[1], cube[2], cube[6], cube[5], cube[9], cube[14], cube[17], cube[13]};
    const fuelsim::Quad8FaceGeometry geometry = fuelsim::make_quad8_face_geometry(face);
    fuelsim::Quad8SurfaceContactLocalValues state{}, committed{};
    for (std::size_t node = 0; node < 4; ++node) {
        state[node] = 400.0;
        state[4 + node] = 300.0;
    }
    for (std::size_t node = 0; node < 8; ++node)
        state[8 + node] = -0.01;
    fuelsim::Quad8ToQuad8HeatGeometry heat_geometry{face,
        face,
        geometry.thermal_points[0].temperature_shape,
        geometry.thermal_points[0].displacement_shape,
        geometry.thermal_points[0].derivative_xi,
        geometry.thermal_points[0].derivative_eta,
        geometry.thermal_points[0].quadrature_weight,
        -1.0};
    fuelsim::Quad8SurfaceContactLocalJacobian heat_jacobian{};
    const fuelsim::GapHeatProperties heat_properties{2.0, 1.0e-6};
    const auto heat_residual =
        fuelsim::compute_quad8_to_quad8_gap_heat(heat_properties, heat_geometry, state, &heat_jacobian);
    double heat_balance = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        heat_balance += heat_residual[node];
    bool passed = check(std::abs(heat_balance) < 1.0e-12 * std::max(1.0, std::abs(heat_residual[0])),
        "HEX20 Q8/Q4 thermal contact residual is exactly conservative at a 2x2 quadrature point");
    const auto heat_value = fuelsim::compute_quad8_to_quad8_gap_heat_value(heat_properties, heat_geometry, state);
    passed = check(heat_value.projected && heat_value.weighted_measure > 0.0 && heat_value.heat_flux > 0.0,
                 "HEX20 thermal contact projects through the quadratic face and produces finite heat flux")
             && passed;
    auto heat_plus = state, heat_minus = state;
    constexpr double thermal_step = 1.0e-7;
    heat_plus[0] += thermal_step;
    heat_minus[0] -= thermal_step;
    const auto heat_plus_residual = fuelsim::compute_quad8_to_quad8_gap_heat(heat_properties, heat_geometry, heat_plus);
    const auto heat_minus_residual =
        fuelsim::compute_quad8_to_quad8_gap_heat(heat_properties, heat_geometry, heat_minus);
    double thermal_jacobian_error = 0.0, thermal_jacobian_scale = 0.0;
    for (std::size_t row = 0; row < heat_plus_residual.size(); ++row) {
        const double numerical = (heat_plus_residual[row] - heat_minus_residual[row]) / (2.0 * thermal_step);
        const double analytic = heat_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count];
        thermal_jacobian_error = std::max(thermal_jacobian_error, std::abs(analytic - numerical));
        thermal_jacobian_scale = std::max({thermal_jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(thermal_jacobian_error / thermal_jacobian_scale < 2.0e-6,
                 "HEX20 thermal contact automatic-differentiation Jacobian matches a centered difference")
             && passed;

    fuelsim::NodeToQuad8ContactGeometry mechanical_geometry{face, face, {}, {}, {}, {}, 0, -1.0};
    for (std::size_t q = 0; q < 9; ++q) {
        mechanical_geometry.secondary_shapes[q] = geometry.mechanical_points[q].displacement_shape;
        mechanical_geometry.secondary_derivatives_xi[q] = geometry.mechanical_points[q].derivative_xi;
        mechanical_geometry.secondary_derivatives_eta[q] = geometry.mechanical_points[q].derivative_eta;
        mechanical_geometry.secondary_quadrature_weights[q] = geometry.mechanical_points[q].quadrature_weight;
    }
    const fuelsim::NormalContactProperties mechanical_properties{1.0e5, 0.0, false};
    fuelsim::Quad8SurfaceContactLocalJacobian mechanical_jacobian{};
    const auto mechanical_residual = fuelsim::compute_node_to_quad8_contact(mechanical_properties,
        mechanical_geometry,
        state,
        committed,
        {},
        &mechanical_jacobian);
    double force_balance_x = 0.0, force_balance_y = 0.0, force_balance_z = 0.0;
    for (std::size_t node = 0; node < 16; ++node) {
        force_balance_x += mechanical_residual[8 + node];
        force_balance_y += mechanical_residual[24 + node];
        force_balance_z += mechanical_residual[40 + node];
    }
    const auto mechanical_value =
        fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, mechanical_geometry, state, committed, {});
    std::array<double, 8> positive_areas{};
    for (std::size_t node = 0; node < positive_areas.size(); ++node) {
        auto node_geometry = mechanical_geometry;
        node_geometry.secondary_local_node = node;
        positive_areas[node] =
            fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, node_geometry, state, committed, {})
                .tributary_area;
    }
    double positive_area_sum = 0.0;
    for (const double area : positive_areas)
        positive_area_sum += area;
    auto consistent_geometry = mechanical_geometry;
    consistent_geometry.nodal_area_rule = fuelsim::Quad8NodalAreaRule::consistent_shape;
    const auto consistent_corner =
        fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, consistent_geometry, state, committed, {});
    consistent_geometry.secondary_local_node = 4;
    const auto consistent_midpoint =
        fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, consistent_geometry, state, committed, {});
    passed = check(mechanical_value.projected && mechanical_value.pressure > 0.0
                       && near(mechanical_value.tributary_area, 3.0 / 76.0, 1.0e-13),
                 "HEX20 Q8 mechanical contact uses the exact positive-lumped corner area")
             && check(near(positive_areas[4], 4.0 / 19.0, 1.0e-13) && near(positive_area_sum, 1.0, 1.0e-13),
                 "HEX20 positive-lumped edge areas are exact and sum to the current face area")
             && check(near(consistent_corner.tributary_area, -1.0 / 12.0, 1.0e-13)
                          && near(consistent_midpoint.tributary_area, 1.0 / 3.0, 1.0e-13),
                 "HEX20 retains the signed consistent-shape area rule for explicit comparisons")
             && check(std::abs(force_balance_x) < 1.0e-12 && std::abs(force_balance_y) < 1.0e-12
                          && std::abs(force_balance_z) < 1.0e-12,
                 "HEX20 Q8 mechanical contact residual is action-reaction conservative")
             && passed;
    fuelsim::NormalContactProperties friction_properties{1.0e5, 0.2, false};
    fuelsim::ContactPointHistory accumulated_slip_history;
    accumulated_slip_history.cartesian_total_tangential_slip = {0.0, 2.0e-4, -3.0e-4};
    const auto open_node = fuelsim::compute_node_to_quad8_contact_value(friction_properties,
        mechanical_geometry,
        committed,
        committed,
        accumulated_slip_history);
    bool open_node_slip_is_retained = open_node.projected && open_node.pressure == 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        open_node_slip_is_retained = open_node_slip_is_retained
                                     && near(open_node.tangential_slip[component],
                                         accumulated_slip_history.cartesian_total_tangential_slip[component],
                                         1.0e-12);
    passed = check(open_node_slip_is_retained,
                 "HEX20 node-to-face accumulated total tangential slip remains constant while contact is open")
             && passed;
    auto sliding_state = state;
    sliding_state[24] = 0.01;
    const auto sliding_value = fuelsim::compute_node_to_quad8_contact_value(friction_properties,
        mechanical_geometry,
        sliding_state,
        committed,
        {});
    passed = check(sliding_value.sliding
                       && sliding_value.tangential_traction
                              <= friction_properties.friction_coefficient * sliding_value.pressure * (1.0 + 1.0e-12),
                 "HEX20 positive-lumped corner contact enters sliding with a bounded tangential traction")
             && passed;
    double jacobian_error = 0.0, jacobian_scale = 0.0;
    constexpr double step = 1.0e-7;
    auto plus = state, minus = state;
    plus[8] += step;
    minus[8] -= step;
    const auto plus_residual =
        fuelsim::compute_node_to_quad8_contact(mechanical_properties, mechanical_geometry, plus, committed, {});
    const auto minus_residual =
        fuelsim::compute_node_to_quad8_contact(mechanical_properties, mechanical_geometry, minus, committed, {});
    for (std::size_t row = 0; row < plus_residual.size(); ++row) {
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double analytic = mechanical_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 8];
        jacobian_error = std::max(jacobian_error, std::abs(analytic - numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(jacobian_error / jacobian_scale < 2.0e-5,
                 "HEX20 node-to-surface mechanical-contact Jacobian matches a centered difference")
             && passed;

    fuelsim::Quad8SurfaceContactLocalResidual surface_residual{}, surface_plus_residual{}, surface_minus_residual{},
        surface_rotate_plus_residual{}, surface_rotate_minus_residual{};
    fuelsim::Quad8SurfaceContactLocalJacobian surface_jacobian{};
    auto rotate_plus = state, rotate_minus = state;
    rotate_plus[17] += step;
    rotate_minus[17] -= step;
    for (const fuelsim::Quad8FaceMechanicalQuadraturePoint& quadrature : geometry.mechanical_points) {
        const fuelsim::Quad8ReferenceProjectionValue reference =
            fuelsim::compute_quad8_reference_projection(face, face, quadrature.displacement_shape, -1.0);
        const fuelsim::Quad8ToQuad8MechanicalGeometry surface_geometry{face,
            face,
            quadrature.displacement_shape,
            quadrature.derivative_xi,
            quadrature.derivative_eta,
            reference.primary_shape,
            reference.primary_derivative_xi,
            reference.primary_derivative_eta,
            quadrature.quadrature_weight,
            -1.0};
        fuelsim::Quad8SurfaceContactLocalJacobian point_jacobian{};
        const auto point_residual = fuelsim::compute_quad8_to_quad8_contact(mechanical_properties,
            surface_geometry,
            state,
            committed,
            {},
            &point_jacobian);
        const auto point_plus =
            fuelsim::compute_quad8_to_quad8_contact(mechanical_properties, surface_geometry, plus, committed, {});
        const auto point_minus =
            fuelsim::compute_quad8_to_quad8_contact(mechanical_properties, surface_geometry, minus, committed, {});
        const auto point_rotate_plus = fuelsim::compute_quad8_to_quad8_contact(mechanical_properties,
            surface_geometry,
            rotate_plus,
            committed,
            {});
        const auto point_rotate_minus = fuelsim::compute_quad8_to_quad8_contact(mechanical_properties,
            surface_geometry,
            rotate_minus,
            committed,
            {});
        for (std::size_t entry = 0; entry < surface_residual.size(); ++entry) {
            surface_residual[entry] += point_residual[entry];
            surface_plus_residual[entry] += point_plus[entry];
            surface_minus_residual[entry] += point_minus[entry];
            surface_rotate_plus_residual[entry] += point_rotate_plus[entry];
            surface_rotate_minus_residual[entry] += point_rotate_minus[entry];
        }
        for (std::size_t entry = 0; entry < surface_jacobian.size(); ++entry)
            surface_jacobian[entry] += point_jacobian[entry];
    }
    force_balance_x = 0.0;
    force_balance_y = 0.0;
    force_balance_z = 0.0;
    for (std::size_t node = 0; node < 16; ++node) {
        force_balance_x += surface_residual[8 + node];
        force_balance_y += surface_residual[24 + node];
        force_balance_z += surface_residual[40 + node];
    }
    const double total_secondary_force =
        std::accumulate(surface_residual.begin() + 8, surface_residual.begin() + 16, 0.0);
    passed = check(near(std::abs(total_secondary_force), 1.0e3, 1.0e-12),
                 "HEX20 surface-to-surface contact integrates constant pressure over the current face")
             && check(near(surface_residual[8], -total_secondary_force / 12.0, 2.0e-12)
                          && near(surface_residual[12], total_secondary_force / 3.0, 2.0e-12),
                 "HEX20 surface-to-surface contact retains the exact signed Q8 consistent nodal forces")
             && check(std::abs(force_balance_x) < 1.0e-12 && std::abs(force_balance_y) < 1.0e-12
                          && std::abs(force_balance_z) < 1.0e-12,
                 "HEX20 surface-to-surface contact is exactly action-reaction conservative")
             && passed;
    auto surface_sliding_state = state;
    for (std::size_t node = 0; node < 8; ++node)
        surface_sliding_state[24 + node] += 0.01;
    const fuelsim::Quad8FaceMechanicalQuadraturePoint& center = geometry.mechanical_points[4];
    const fuelsim::Quad8ReferenceProjectionValue center_reference =
        fuelsim::compute_quad8_reference_projection(face, face, center.displacement_shape, -1.0);
    const fuelsim::Quad8ToQuad8MechanicalGeometry center_geometry{face,
        face,
        center.displacement_shape,
        center.derivative_xi,
        center.derivative_eta,
        center_reference.primary_shape,
        center_reference.primary_derivative_xi,
        center_reference.primary_derivative_eta,
        center.quadrature_weight,
        -1.0};
    auto finite_sliding_geometry = center_geometry;
    finite_sliding_geometry.finite_sliding = true;
    auto finite_sliding_state = state;
    for (std::size_t node = 0; node < 8; ++node)
        finite_sliding_state[24 + node] += 0.25;
    fuelsim::Quad8SurfaceContactLocalJacobian finite_sliding_jacobian{};
    const auto finite_sliding_residual = fuelsim::compute_quad8_to_quad8_contact(mechanical_properties,
        finite_sliding_geometry,
        finite_sliding_state,
        committed,
        {},
        &finite_sliding_jacobian);
    double primary_force = 0.0, primary_force_y_moment = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        primary_force += finite_sliding_residual[16 + node];
        primary_force_y_moment += face[node].y * finite_sliding_residual[16 + node];
    }
    passed = check(std::abs(primary_force_y_moment / primary_force - 0.75) < 1.0e-12,
                 "HEX20 finite sliding moves the primary test-field projection to the current closest point")
             && passed;
    auto finite_plus = finite_sliding_state, finite_minus = finite_sliding_state;
    finite_plus[25] += step;
    finite_minus[25] -= step;
    const auto finite_plus_residual = fuelsim::compute_quad8_to_quad8_contact(mechanical_properties,
        finite_sliding_geometry,
        finite_plus,
        committed,
        {});
    const auto finite_minus_residual = fuelsim::compute_quad8_to_quad8_contact(mechanical_properties,
        finite_sliding_geometry,
        finite_minus,
        committed,
        {});
    double finite_jacobian_error = 0.0, finite_jacobian_scale = 0.0;
    for (std::size_t row = 0; row < finite_plus_residual.size(); ++row) {
        const double numerical = (finite_plus_residual[row] - finite_minus_residual[row]) / (2.0 * step);
        const double analytic = finite_sliding_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 25];
        finite_jacobian_error = std::max(finite_jacobian_error, std::abs(analytic - numerical));
        finite_jacobian_scale = std::max({finite_jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(finite_jacobian_error / finite_jacobian_scale < 2.0e-5,
                 "HEX20 finite-sliding closest-point Jacobian matches a centered directional difference")
             && passed;
    auto outside_state = state;
    for (std::size_t node = 0; node < 8; ++node)
        outside_state[24 + node] += 0.6;
    passed =
        check(!fuelsim::compute_quad8_to_quad8_contact_projection(finite_sliding_geometry, outside_state).projected,
            "HEX20 finite sliding rejects a closest point outside the complete primary face")
        && passed;
    const auto surface_sliding = fuelsim::compute_quad8_to_quad8_contact_value(friction_properties,
        center_geometry,
        surface_sliding_state,
        committed,
        {});
    const auto open_surface = fuelsim::compute_quad8_to_quad8_contact_value(friction_properties,
        center_geometry,
        committed,
        committed,
        accumulated_slip_history);
    bool open_surface_slip_is_retained = open_surface.projected && open_surface.pressure == 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        open_surface_slip_is_retained = open_surface_slip_is_retained
                                        && near(open_surface.tangential_slip[component],
                                            accumulated_slip_history.cartesian_total_tangential_slip[component],
                                            1.0e-12);
    passed = check(surface_sliding.sliding
                       && near(surface_sliding.tangential_traction,
                           friction_properties.friction_coefficient * surface_sliding.pressure,
                           1.0e-12),
                 "HEX20 surface-to-surface Coulomb friction caps the integration-point tangential traction")
             && check(open_surface_slip_is_retained,
                 "HEX20 surface-to-surface accumulated total tangential slip remains constant while contact is open")
             && passed;
    auto biaxial_sliding_state = state;
    for (std::size_t node = 0; node < 8; ++node) {
        biaxial_sliding_state[24 + node] += 0.006;
        biaxial_sliding_state[40 + node] += 0.008;
    }
    const fuelsim::NormalContactProperties controlled_slip_properties{1.0e5, 0.2, false, 4.0e-4};
    const auto biaxial_sliding = fuelsim::compute_quad8_to_quad8_contact_value(controlled_slip_properties,
        center_geometry,
        biaxial_sliding_state,
        committed,
        {});
    const double biaxial_elastic_slip = std::hypot(biaxial_sliding.elastic_tangential_slip[0],
        std::hypot(biaxial_sliding.elastic_tangential_slip[1], biaxial_sliding.elastic_tangential_slip[2]));
    passed = check(biaxial_sliding.sliding && std::abs(biaxial_sliding.tangential_traction_vector[1]) > 0.0
                       && std::abs(biaxial_sliding.tangential_traction_vector[2]) > 0.0,
                 "HEX20 surface-to-surface friction enters sliding with two simultaneous tangential components")
             && check(near(biaxial_sliding.tangential_traction,
                          controlled_slip_properties.friction_coefficient * biaxial_sliding.pressure,
                          1.0e-12)
                          && near(biaxial_elastic_slip, controlled_slip_properties.maximum_elastic_slip, 1.0e-12),
                 "HEX20 biaxial sliding satisfies the Coulomb circle and requested maximum elastic slip")
             && passed;

    constexpr double curved_angle = 0.25;
    const fuelsim::Quad8FaceCoordinates curved_face = {
        fuelsim::CartesianPoint3{std::cos(-curved_angle), std::sin(-curved_angle), 0.0},
        fuelsim::CartesianPoint3{std::cos(curved_angle), std::sin(curved_angle), 0.0},
        fuelsim::CartesianPoint3{std::cos(curved_angle), std::sin(curved_angle), 1.0},
        fuelsim::CartesianPoint3{std::cos(-curved_angle), std::sin(-curved_angle), 1.0},
        fuelsim::CartesianPoint3{1.0, 0.0, 0.0},
        fuelsim::CartesianPoint3{std::cos(curved_angle), std::sin(curved_angle), 0.5},
        fuelsim::CartesianPoint3{1.0, 0.0, 1.0},
        fuelsim::CartesianPoint3{std::cos(-curved_angle), std::sin(-curved_angle), 0.5}};
    const fuelsim::Quad8FaceGeometry curved_geometry = fuelsim::make_quad8_face_geometry(curved_face);
    const fuelsim::Quad8FaceMechanicalQuadraturePoint& curved_center = curved_geometry.mechanical_points[4];
    const fuelsim::Quad8ReferenceProjectionValue curved_reference =
        fuelsim::compute_quad8_reference_projection(curved_face, curved_face, curved_center.displacement_shape, -1.0);
    const fuelsim::Quad8ToQuad8MechanicalGeometry curved_contact_geometry{curved_face,
        curved_face,
        curved_center.displacement_shape,
        curved_center.derivative_xi,
        curved_center.derivative_eta,
        curved_reference.primary_shape,
        curved_reference.primary_derivative_xi,
        curved_reference.primary_derivative_eta,
        curved_center.quadrature_weight,
        -1.0};
    fuelsim::Quad8SurfaceContactLocalValues curved_committed{}, curved_rotated{};
    constexpr double penetration = 0.01, rotation = 0.35;
    constexpr std::array<double, 3> committed_slip = {0.0, 6.0e-4, 8.0e-4};
    const double rotation_cosine = std::cos(rotation), rotation_sine = std::sin(rotation);
    for (std::size_t node = 0; node < 16; ++node) {
        const fuelsim::CartesianPoint3& coordinate = curved_face[node % 8];
        curved_rotated[8 + node] = rotation_cosine * coordinate.x - rotation_sine * coordinate.y - coordinate.x;
        curved_rotated[24 + node] = rotation_sine * coordinate.x + rotation_cosine * coordinate.y - coordinate.y;
        if (node < 8) {
            curved_committed[8 + node] -= penetration;
            curved_committed[24 + node] += committed_slip[1];
            curved_committed[40 + node] += committed_slip[2];
            curved_rotated[8 + node] += -rotation_cosine * penetration - rotation_sine * committed_slip[1];
            curved_rotated[24 + node] += -rotation_sine * penetration + rotation_cosine * committed_slip[1];
            curved_rotated[40 + node] += committed_slip[2];
        }
    }
    fuelsim::ContactPointHistory curved_history{};
    curved_history.cartesian_elastic_tangential_slip = committed_slip;
    const fuelsim::NormalContactProperties curved_friction_properties{1.0e5, 0.5, false};
    fuelsim::Quad8SurfaceContactLocalJacobian curved_contact_jacobian{};
    const auto curved_history_value = fuelsim::compute_quad8_to_quad8_contact_value(curved_friction_properties,
        curved_contact_geometry,
        curved_rotated,
        curved_committed,
        curved_history);
    (void)fuelsim::compute_quad8_to_quad8_contact(curved_friction_properties,
        curved_contact_geometry,
        curved_rotated,
        curved_committed,
        curved_history,
        &curved_contact_jacobian);
    const std::array<double, 3> rotated_history = {-rotation_sine * 6.0e-4, rotation_cosine * 6.0e-4, 8.0e-4};
    double history_error = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        history_error = std::max(history_error,
            std::abs(curved_history_value.elastic_tangential_slip[component] - rotated_history[component]));
    const double history_normal = curved_history_value.elastic_tangential_slip[0] * curved_history_value.normal[0]
                                  + curved_history_value.elastic_tangential_slip[1] * curved_history_value.normal[1]
                                  + curved_history_value.elastic_tangential_slip[2] * curved_history_value.normal[2];
    passed = check(curved_history_value.projected && !curved_history_value.sliding && history_error < 1.0e-13,
                 "HEX20 curved surface rotates a nonzero two-component elastic-slip history with its tangent plane")
             && check(std::abs(history_normal) < 1.0e-13
                          && near(std::hypot(curved_history_value.elastic_tangential_slip[0],
                                      std::hypot(curved_history_value.elastic_tangential_slip[1],
                                          curved_history_value.elastic_tangential_slip[2])),
                              1.0e-3,
                              1.0e-13),
                 "HEX20 curved tangent-plane transport preserves tangency and elastic-slip magnitude")
             && passed;

    auto curved_plus = curved_rotated, curved_minus = curved_rotated;
    curved_plus[17] += step;
    curved_minus[17] -= step;
    const auto curved_plus_residual = fuelsim::compute_quad8_to_quad8_contact(curved_friction_properties,
        curved_contact_geometry,
        curved_plus,
        curved_committed,
        curved_history);
    const auto curved_minus_residual = fuelsim::compute_quad8_to_quad8_contact(curved_friction_properties,
        curved_contact_geometry,
        curved_minus,
        curved_committed,
        curved_history);
    double curved_jacobian_error = 0.0, curved_jacobian_scale = 0.0;
    for (std::size_t row = 0; row < curved_plus_residual.size(); ++row) {
        const double numerical = (curved_plus_residual[row] - curved_minus_residual[row]) / (2.0 * step);
        const double analytic = curved_contact_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 17];
        curved_jacobian_error = std::max(curved_jacobian_error, std::abs(analytic - numerical));
        curved_jacobian_scale = std::max({curved_jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(curved_jacobian_error / curved_jacobian_scale < 2.0e-5,
                 "HEX20 curved history-rotation Jacobian matches a centered directional difference")
             && passed;
    jacobian_error = 0.0;
    jacobian_scale = 0.0;
    for (std::size_t row = 0; row < surface_residual.size(); ++row) {
        const double numerical = (surface_plus_residual[row] - surface_minus_residual[row]) / (2.0 * step);
        const double analytic = surface_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 8];
        jacobian_error = std::max(jacobian_error, std::abs(analytic - numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(analytic), std::abs(numerical)});
        const double rotating_numerical =
                         (surface_rotate_plus_residual[row] - surface_rotate_minus_residual[row]) / (2.0 * step),
                     rotating_analytic = surface_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 17];
        jacobian_error = std::max(jacobian_error, std::abs(rotating_analytic - rotating_numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(rotating_analytic), std::abs(rotating_numerical)});
    }
    return check(jacobian_error / jacobian_scale < 2.0e-5,
               "HEX20 small-sliding contact gap and rotating-normal Jacobian match centered differences")
           && passed;
}

bool test_reduced_volume_kernel() {
    const auto coordinates = unit_cube();
    const auto geometry = fuelsim::make_c3d20_geometry(coordinates, fuelsim::Hex20ElementFormulation::c3d20rt);
    bool passed = check(geometry.mechanical_points.size() == 8, "C3D20RT has eight active material points");
    double volume = 0.0;
    for (const auto& point : geometry.mechanical_points)
        volume += point.weighted_measure;
    passed = check(near(volume, 1.0, 1.0e-14), "C3D20RT integrates the unit cube volume") && passed;
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        for (const unsigned mechanism : {0U, 1U, 2U, 3U}) {
            auto properties = material();
            if ((mechanism & 1U) != 0)
                properties = fuelsim::test::with_norton(std::move(properties), 1.0e-6, 10.0, 3.0, 300.0);
            if ((mechanism & 2U) != 0)
                properties = fuelsim::test::with_plasticity(std::move(properties), 20.0, 10.0, 300.0);
            const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties),
                4.0,
                1.0,
                formulation};
            fuelsim::Hex20LocalValues old{}, state{};
            for (std::size_t node = 0; node < 8; ++node) {
                old[node] = 300.0;
                state[node] = 302.0 + coordinates[node].x + 2.0 * coordinates[node].y;
            }
            for (std::size_t node = 0; node < 20; ++node) {
                state[8 + node] = 0.002 * coordinates[node].x + 0.0003 * coordinates[node].y;
                state[28 + node] = -0.0004 * coordinates[node].y;
                state[48 + node] = 0.0002 * coordinates[node].z;
            }
            const fuelsim::CartesianMaterialHistory history(8);
            passed = check(directional_jacobian_error(data, geometry, state, old, history) < 2.0e-6,
                         "C3D20RT elastic/creep/plastic/coupled small/finite tangent matches centered differences")
                     && passed;
            passed = check(residual_path_error(data, geometry, state, old, history) < 1.0e-13,
                         "C3D20RT independent residual matches the Jacobian residual")
                     && passed;
            const auto updated = fuelsim::compute_c3d20_transient_update(data, geometry, state, old, history, 0.5);
            passed = check(updated.size() == 8, "C3D20RT updates exactly eight material histories") && passed;
            if (mechanism == 0)
                passed = check(fuelsim::compute_c3d20_stress(data, geometry, state).size() == 8,
                             "C3D20RT stress recovery returns only active points")
                         && passed;
        }
    }
    return passed;
}

bool test_primary_projection_shape_derivatives() {
    const std::array<std::array<double, 2>, 8> natural = {
        {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {1, 0}, {0, 1}, {-1, 0}}};
    for (std::size_t sample = 0; sample < 8; ++sample) {
        fuelsim::Quad8ToQuad8HeatGeometry geometry{};
        for (std::size_t node = 0; node < 8; ++node) {
            const double x = natural[node][0], y = natural[node][1];
            geometry.primary_coordinates[node] = {x, y, 0.03 * x * y + 0.01 * static_cast<double>(sample) * x * x};
            geometry.secondary_coordinates[node] = {x + 0.02, y - 0.01, 0.12 + 0.02 * x * y};
        }
        geometry.normal_orientation = 1.0;
        geometry.secondary_displacement_shape =
            fuelsim::make_quad8_face_mechanical_point(geometry.secondary_coordinates,
                -0.4 + 0.1 * static_cast<double>(sample),
                0.23,
                1.0)
                .displacement_shape;
        fuelsim::Quad8SurfaceContactLocalValues state{}, direction{};
        for (std::size_t i = 0; i < state.size(); ++i) {
            state[i] = 0.002 * std::sin(static_cast<double>(i + sample));
            direction[i] = 0.1 * std::cos(0.37 * static_cast<double>(i + 1));
        }
        const auto analytic = fuelsim::compute_quad8_primary_shape_derivatives(geometry, state);
        const auto projected = [&](double step) {
            auto secondary = geometry.secondary_coordinates, primary = geometry.primary_coordinates;
            for (std::size_t node = 0; node < 16; ++node) {
                auto& point = node < 8 ? secondary[node] : primary[node - 8];
                point.x += state[8 + node] + step * direction[8 + node];
                point.y += state[24 + node] + step * direction[24 + node];
                point.z += state[40 + node] + step * direction[40 + node];
            }
            return fuelsim::compute_quad8_reference_projection(secondary,
                primary,
                geometry.secondary_displacement_shape,
                1.0);
        };
        constexpr double h = 1e-6;
        const auto plus = projected(h), minus = projected(-h);
        if (!check(plus.projected && minus.projected,
                "Curved projection derivative samples remain inside the primary face"))
            return false;
        double error = 0.0, scale = 0.0;
        for (std::size_t node = 0; node < 8; ++node) {
            double tangent = 0.0;
            for (std::size_t i = 0; i < state.size(); ++i)
                tangent += analytic[node].derivative(i) * direction[i];
            const double numerical = (plus.primary_shape[node] - minus.primary_shape[node]) / (2 * h);
            error = std::hypot(error, tangent - numerical);
            scale = std::hypot(scale, numerical);
        }
        if (!check(error < 1e-7 * scale,
                "Q8 implicit projection shape derivatives match independent centered projections"))
            return false;
    }
    return true;
}

bool test_disk_transfer() {
    const std::array<std::array<double, 2>, 8> natural = {
        {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {1, 0}, {0, 1}, {-1, 0}}};
    fuelsim::Quad8ToQuad8HeatGeometry geometry{};
    for (std::size_t i = 0; i < 8; ++i) {
        const double x = 0.5 * (natural[i][0] + 1), y = 0.5 * (natural[i][1] + 1);
        geometry.secondary_coordinates[i] = {x, y, 0.01};
        geometry.primary_coordinates[i] = {0.25 + 0.125 * x, 0.125 * y, 0.0};
    }
    geometry.normal_orientation = 1.0;
    auto point = fuelsim::make_quad8_face_mechanical_point(geometry.secondary_coordinates,
        2 * 0.38603095325083064 - 1,
        2 * 0.047248941508978792 - 1,
        1.0);
    geometry.secondary_displacement_shape = point.displacement_shape;
    geometry.secondary_derivative_xi = point.derivative_xi;
    geometry.secondary_derivative_eta = point.derivative_eta;
    auto values = fuelsim::compute_quad8_disk_transfer(geometry, {});
    double fraction = 0.0;
    for (double value : values)
        fraction += value;
    bool passed = check(near(fraction, 0.2580985357387237, 1e-8),
        "Native finite contact distributes the sampling disk onto an extrapolated neighboring face");
    auto neighbor = geometry;
    for (auto& coordinate : neighbor.primary_coordinates)
        coordinate.x += 0.125;
    const auto next = fuelsim::compute_quad8_disk_transfer(neighbor, {});
    double sum = fraction;
    for (double value : next)
        sum += value;
    passed = check(near(sum, 1.0, 1e-12), "Adjacent primary disk fractions preserve total force") && passed;
    passed = check(values[4] < 0.0 || values[5] < 0.0 || values[6] < 0.0 || values[7] < 0.0,
                 "Extrapolated primary midpoint shapes retain their signed coefficients")
             && passed;
    fuelsim::Quad8SurfaceContactLocalValues state{}, direction{};
    for (std::size_t i = 8; i < state.size(); ++i) {
        state[i] = 2e-4 * std::sin(0.73 * static_cast<double>(i));
        direction[i] = 0.1 * std::cos(0.37 * static_cast<double>(i));
    }
    std::array<adlite::Scalar, 8> derivatives{};
    values = fuelsim::compute_quad8_disk_transfer(geometry, state, &derivatives);
    auto plus = state, minus = state;
    constexpr double h = 1e-6;
    for (std::size_t i = 0; i < state.size(); ++i) {
        plus[i] += h * direction[i];
        minus[i] -= h * direction[i];
    }
    const auto positive = fuelsim::compute_quad8_disk_transfer(geometry, plus),
               negative = fuelsim::compute_quad8_disk_transfer(geometry, minus),
               repeated = fuelsim::compute_quad8_disk_transfer(geometry, state);
    double error = 0.0, scale = 0.0;
    for (std::size_t row = 0; row < 8; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            analytic += derivatives[row].derivative(column) * direction[column];
        const double numerical = (positive[row] - negative[row]) / (2 * h);
        error = std::hypot(error, analytic - numerical);
        scale = std::hypot(scale, numerical);
        passed = check(values[row] == repeated[row] && derivatives[row].value() == values[row],
                     "Disk transfer ordinary and derivative calls have identical values")
                 && passed;
    }
    passed =
        check(error < 1e-6 * scale, "Disk area, radius, plane and extrapolation derivatives match centered differences")
        && passed;
    std::vector<fuelsim::Quad8HeatPatchSample> heat_samples;
    const double a = 0.38603095325083064, b = 0.047248941508978792;
    for (std::size_t face = 0; face < 2; ++face) {
        fuelsim::Quad8HeatPatchSample sample{face == 0 ? geometry : neighbor, {}, true, 0};
        sample.geometry.quadrature_weight = 1.0;
        sample.geometry.secondary_temperature_shape = {(1 - a) * (1 - b), a * (1 - b), a * b, (1 - a) * b};
        for (std::size_t i = 0; i < 4; ++i) {
            sample.local_dofs[i] = i;
            sample.local_dofs[4 + i] = 4 + 4 * face + i;
        }
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t i = 0; i < 8; ++i) {
                sample.local_dofs[8 + 16 * c + i] = 12 + 24 * c + i;
                sample.local_dofs[16 + 16 * c + i] = 20 + 24 * c + 8 * face + i;
            }
        heat_samples.push_back(sample);
    }
    std::vector<double> heat_state(84), heat_direction(84), heat_plus(84), heat_minus(84), tangent;
    for (std::size_t i = 0; i < 84; ++i) {
        heat_state[i] = i < 12 ? 300 + static_cast<double>(i) : 2e-4 * std::sin(.73 * static_cast<double>(i));
        heat_direction[i] =
            i < 12 ? std::cos(.37 * static_cast<double>(i)) : .1 * std::cos(.37 * static_cast<double>(i));
        heat_plus[i] = heat_state[i] + h * heat_direction[i];
        heat_minus[i] = heat_state[i] - h * heat_direction[i];
    }
    const fuelsim::GapHeatProperties heat_properties{.001, 1e-5};
    const auto heat = fuelsim::compute_quad8_gap_heat_patch(heat_properties, heat_samples, heat_state, &tangent),
               heat_positive = fuelsim::compute_quad8_gap_heat_patch(heat_properties, heat_samples, heat_plus),
               heat_negative = fuelsim::compute_quad8_gap_heat_patch(heat_properties, heat_samples, heat_minus);
    error = 0.;
    scale = 0.;
    double total_heat = 0.;
    for (std::size_t row = 0; row < 84; ++row) {
        total_heat += heat[row];
        double analytic = 0.;
        for (std::size_t column = 0; column < 84; ++column)
            analytic += tangent[row * 84 + column] * heat_direction[column];
        const double numerical = (heat_positive[row] - heat_negative[row]) / (2 * h);
        error = std::hypot(error, analytic - numerical);
        scale = std::hypot(scale, numerical);
    }
    return check(error < 1e-6 * scale, "Thermal disk transfer retains cross-face geometric and temperature derivatives")
           && check(std::abs(total_heat) < 1e-12, "Thermal disk transfer conserves heat across both primary faces")
           && passed;
}

} // namespace

int run_c3d20t_tests() {
    const bool passed = test_geometry_and_constant_strain() && test_material_value_paths()
                        && test_finite_thermal_operators() && test_jacobian_and_transient_history()
                        && test_quadratic_face() && test_warped_geometry() && test_hex20_heat_patch()
                        && test_hex20_contact_kernels() && test_primary_projection_shape_derivatives()
                        && test_disk_transfer();
    if (passed)
        std::cout << "All HEX20-U2/T1 kernel tests passed\n";
    return passed ? 0 : 1;
}

int run_c3d20rt_tests() {
    return test_reduced_volume_kernel() ? 0 : 1;
}
