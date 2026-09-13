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

} // namespace fuelsim::test::c3d20
