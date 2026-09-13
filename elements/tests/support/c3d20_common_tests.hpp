#pragma once
#include "c3d20_types.hpp"
#include "contact_types.hpp"
#include "quad4_face.hpp"
#include "quad8_face.hpp"
#include "support/test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>

namespace fuelsim::test::c3d20 {
inline bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

inline bool near(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance * std::max({1.0, std::abs(actual), std::abs(expected)});
}

inline fuelsim::Hex20Coordinates unit_cube() {
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

inline fuelsim::ThermoelasticProperties material(bool inelastic = false) {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e5, 0.25, 1.2e-5, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
    if (inelastic) {
        result = fuelsim::test::with_norton(std::move(result), 1.0e-6, 10.0, 3.0, 300.0);
        result = fuelsim::test::with_plasticity(std::move(result), 20.0, 10.0, 300.0);
    }
    return result;
}

inline double directional_jacobian_error(const fuelsim::CartesianTestData& data,
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

inline double residual_path_error(const fuelsim::CartesianTestData& data,
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

inline bool test_reference_mass_capacity(fuelsim::Hex20ElementFormulation element_formulation) {
    const auto coordinates = unit_cube();
    const auto geometry = fuelsim::test::make_c3d20_geometry(coordinates, element_formulation);
    const auto registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>();
    functions->name = "c3d20_reference_mass";
    functions->thermal = registry.bind_thermal("constant_thermophysical",
        {{"conductivity", 4.0}, {"density", 2.0}, {"specific_heat", 3.0}});
    functions->thermal.function = [](const fuelsim::ThermoelasticFunctionInput& input,
                                      fuelsim::ThermalPropertyOutput& output) {
        output.conductivity = 4.0;
        output.density = 2.0 + 0.01 * (input.temperature - 300.0) + 5.0 * input.context.time + 0.2 * input.context.x;
        output.specific_heat = 3.0 + 0.02 * (input.temperature - 300.0) + 0.1 * input.context.time;
    };
    functions->elasticity =
        registry.bind_elasticity("constant_isotropic", {{"young_modulus", 2.0e5}, {"poisson_ratio", 0.25}});
    const fuelsim::IsotropicThermoelasticMaterial active_material({functions, 2.0e5});
    bool passed = true;
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        const fuelsim::CartesianTestData data{active_material,
            2.0,
            2.0,
            formulation,
            fuelsim::Hex8ElementFormulation::c3d8t,
            320.0};
        fuelsim::Hex20LocalValues old{}, state{};
        for (std::size_t node = 0; node < 8; ++node) {
            old[node] = 350.0;
            state[node] = 360.0;
        }
        for (std::size_t node = 0; node < 20; ++node) {
            const auto& point = coordinates[node];
            old[8 + node] = 0.1 * point.x;
            old[28 + node] = -0.05 * point.y;
            old[48 + node] = 0.03 * point.z;
            state[8 + node] = point.x;
            state[28 + node] = 0.2 * point.y;
            state[48 + node] = 0.1 * point.z;
        }
        const fuelsim::CartesianMaterialHistory history(geometry.mechanical_points.size());
        for (const bool move_midsides : {false, true}) {
            auto current = state;
            if (move_midsides) {
                current[8 + 8] += 0.04;
                current[28 + 18] -= 0.03;
            }
            fuelsim::Hex20LocalJacobian transient_jacobian{}, steady_jacobian{};
            const auto transient =
                fuelsim::compute_c3d20_transient(data, geometry, current, old, history, 0.5, &transient_jacobian, true);
            const auto steady =
                fuelsim::compute_c3d20_transient(data, geometry, current, old, history, 0.5, &steady_jacobian, false);
            const auto values = fuelsim::compute_c3d20_transient(data, geometry, current, old, history, 0.5);
            // rho0(x)=2.2+0.2*x is fixed at T_initial=320 and time zero.
            // cp(360,time=2)=4.4 and dcp/dT=0.02. Integrate the unit-cube
            // polynomials independently of the element's stored quadrature.
            double stored = 0.0;
            for (std::size_t row = 0; row < 8; ++row) {
                const double mean_x = coordinates[row].x == 0.0 ? 1.0 / 3.0 : 2.0 / 3.0;
                const double initial_nodal_mass = (2.2 + 0.2 * mean_x) / 8.0;
                passed = check(near(transient[row] - steady[row], initial_nodal_mass * 4.4 * 20.0, 3.0e-12),
                             "C3D20 capacity fixes initial temperature, time and reference position across deformation")
                         && passed;
                passed = check(near(transient[row], values[row], 1.0e-13),
                             "C3D20 reference-mass thermal residual agrees in double and AD paths")
                         && passed;
                stored += transient[row] - steady[row];
                for (std::size_t column = 0; column < 8; ++column) {
                    const bool same_x = coordinates[row].x == coordinates[column].x;
                    const double integral_x = same_x ? 1.0 / 3.0 : 1.0 / 6.0;
                    const double moment_x = same_x && coordinates[row].x == 1.0 ? 1.0 / 4.0 : 1.0 / 12.0;
                    const double integral_y = coordinates[row].y == coordinates[column].y ? 1.0 / 3.0 : 1.0 / 6.0;
                    const double integral_z = coordinates[row].z == coordinates[column].z ? 1.0 / 3.0 : 1.0 / 6.0;
                    const double mass = (2.2 * integral_x + 0.2 * moment_x) * integral_y * integral_z;
                    const double capacity = mass * (4.4 + 0.02 * 10.0) / 0.5;
                    passed = check(near(transient_jacobian[68 * row + column] - steady_jacobian[68 * row + column],
                                       capacity,
                                       3.0e-12),
                                 "C3D20 consistent capacity retains off-diagonal mass and the specific-heat derivative")
                             && passed;
                }
                for (std::size_t column = 8; column < fuelsim::hex20_local_dof_count; ++column)
                    passed = check(transient_jacobian[68 * row + column] == steady_jacobian[68 * row + column],
                                 "C3D20 initial-mass capacity has exactly zero displacement derivative")
                             && passed;
            }
            passed = check(near(stored, 2.3 * 4.4 * 20.0, 3.0e-12),
                         "C3D20 total heat storage preserves the initial mass under affine and midside deformation")
                     && passed;
            const double source_volume = formulation == fuelsim::StrainFormulation::finite ? 2.64 : 1.0;
            const double source_residual = std::accumulate(steady.begin(), steady.begin() + 8, 0.0);
            passed = check(near(source_residual, -2.0 * source_volume, 3.0e-12),
                         "C3D20 heat source retains its selected geometric volume independently of fixed mass")
                     && passed;
            passed = check(directional_jacobian_error(data, geometry, current, old, history, 0, 8) < 2.0e-6,
                         "C3D20 fixed-mass thermal tangent agrees with a centered temperature/displacement direction")
                     && passed;
        }
    }
    return passed;
}

} // namespace fuelsim::test::c3d20
