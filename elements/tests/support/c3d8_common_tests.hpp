#pragma once
#include "c3d8_types.hpp"
#include "c3d8rt.hpp"
#include "c3d8t.hpp"
#include "contact_types.hpp"
#include "quad4_face.hpp"
#include "support/c3d_recovery.hpp"
#include "support/test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace fuelsim::test::c3d8 {
inline bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

inline bool near(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance * std::max({1.0, std::abs(actual), std::abs(expected)});
}

inline fuelsim::Hex8Coordinates unit_cube() {
    return {{{0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.0, 1.0, 1.0}}};
}

inline fuelsim::ThermoelasticProperties properties() {
    return fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, -1.0e8, 0.0, 1.0e-8, 2000.0, 3000.0);
}

inline fuelsim::ThermoelasticProperties inelastic_properties(bool creep, bool plasticity) {
    fuelsim::ThermoelasticProperties result = fuelsim::test::thermoelastic(0.0, 1.0, 200.0, 0.25, 1.0e-4, 300.0);
    if (creep)
        result = fuelsim::test::with_norton(std::move(result), 1.0e-4, 10.0, 3.0, 300.0);
    if (plasticity)
        result = fuelsim::test::with_plasticity(std::move(result), 10.0, 20.0, 300.0);
    return result;
}

inline fuelsim::ThermalPropertyEvaluator temperature_dependent_capacity(const fuelsim::MaterialParameters& named) {
    const double conductivity = named.value("conductivity"), density = named.value("density"),
                 specific_heat = named.value("specific_heat"), slope = named.value("specific_heat_slope");
    return [conductivity, density, specific_heat, slope](const fuelsim::ThermoelasticFunctionInput& input,
               fuelsim::ThermalPropertyOutput& output) {
        output.conductivity = conductivity;
        output.density = density;
        output.specific_heat = specific_heat + slope * input.temperature;
    };
}

inline fuelsim::ThermoelasticProperties capacity_properties() {
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    registry.add_thermal("temperature_dependent_capacity",
        {{"conductivity", "W/(m*K)"},
            {"density", "kg/m^3"},
            {"specific_heat", "J/(kg*K)"},
            {"specific_heat_slope", "J/(kg*K^2)"}},
        &temperature_dependent_capacity);
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>();
    functions->name = "temperature_dependent_capacity_test";
    functions->thermal = registry.bind_thermal("temperature_dependent_capacity",
        {{"conductivity", 4.0}, {"density", 2000.0}, {"specific_heat", 1000.0}, {"specific_heat_slope", 2.0}});
    functions->elasticity =
        registry.bind_elasticity("constant_isotropic", {{"young_modulus", 2.0e11}, {"poisson_ratio", 0.25}});
    return {std::move(functions), 2.0e11};
}

inline fuelsim::ThermalPropertyEvaluator reference_mass_thermal_properties(const fuelsim::MaterialParameters&) {
    return [](const fuelsim::ThermoelasticFunctionInput& input, fuelsim::ThermalPropertyOutput& output) {
        const auto& position = input.context;
        output.conductivity = 4.0;
        output.density = 2.0 + 0.01 * input.temperature + 0.2 * position.time + 0.1 * position.x + 0.2 * position.y
                         + 0.3 * position.z;
        output.specific_heat = 3.0 + 0.02 * input.temperature + 0.4 * position.time + 0.3 * position.x
                               + 0.2 * position.y + 0.1 * position.z;
    };
}

inline bool test_fixed_initial_mass_capacity(bool reduced) {
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    registry.add_thermal("reference_mass_thermal", {}, &reference_mass_thermal_properties);
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>();
    functions->name = "reference_mass_thermal_test";
    functions->thermal = registry.bind_thermal("reference_mass_thermal", {});
    functions->elasticity =
        registry.bind_elasticity("constant_isotropic", {{"young_modulus", 200.0}, {"poisson_ratio", 0.25}});
    const fuelsim::IsotropicThermoelasticMaterial material({std::move(functions), 200.0});
    auto coordinates = unit_cube();
    constexpr double raised_corner = 0.4, initial_temperature = 280.0, time_step = 2.0;
    coordinates[6].z += raised_corner;
    const auto geometry = reduced ? fuelsim::elements::make_c3d8rt_geometry(coordinates)
                                  : fuelsim::elements::make_c3d8t_geometry(coordinates);
    fuelsim::Hex8LocalValues previous{}, current{};
    std::array<double, 8> reference_mass{};
    for (std::size_t node = 0; node < 8; ++node) {
        const auto& position = coordinates[node];
        previous[node] = 300.0 + 0.6 * static_cast<double>(node);
        current[node] = 315.0 + 2.0 * static_cast<double>(node);
        previous[8 + node] = 0.08 * position.x;
        previous[16 + node] = 0.03 * position.y;
        previous[24 + node] = -0.02 * position.z;
        current[8 + node] = 0.35 * position.x + 0.04 * position.y * position.z;
        current[16 + node] = -0.12 * position.y + 0.03 * position.x * position.z;
        current[24 + node] = 0.25 * position.z + 0.02 * position.x * position.y;
        // For z = t * (1 + raised_corner * x * y), det(dx/d(x,y,t)) = 1 + raised_corner * x * y.
        // Integrate N_i analytically for C3D8RT, or evaluate at the paired Gauss point for C3D8T.
        const double x_gauss = 0.5 + (position.x - 0.5) / std::sqrt(3.0);
        const double y_gauss = 0.5 + (position.y - 0.5) / std::sqrt(3.0);
        const double weight = reduced ? 0.125 + raised_corner * (1.0 + position.x) * (1.0 + position.y) / 72.0
                                      : (1.0 + raised_corner * x_gauss * y_gauss) / 8.0;
        const double initial_density =
            2.0 + 0.01 * initial_temperature + 0.1 * position.x + 0.2 * position.y + 0.3 * position.z;
        reference_mass[node] = weight * initial_density;
    }
    const fuelsim::CartesianMaterialHistory history(reduced ? 1 : 8);
    const auto evaluate = reduced ? fuelsim::elements::evaluate_c3d8rt : fuelsim::elements::evaluate_c3d8t;
    double residual_error = 0.0, jacobian_error = 0.0, displacement_derivative = 0.0, stored_power_error = 0.0;
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite})
        for (const double time : {7.0, 11.0})
            for (const bool with_history : {false, true}) {
                fuelsim::elements::C3d8Input input{material,
                    geometry,
                    current,
                    previous,
                    with_history ? &history : nullptr,
                    time_step,
                    time,
                    3.0,
                    formulation,
                    true,
                    initial_temperature};
                const auto active = evaluate(input, {true, true, false, false});
                const auto passive = evaluate(input, {true, false, false, false});
                input.include_thermal_time_term = false;
                const auto baseline = evaluate(input, {true, true, false, false});
                const auto passive_baseline = evaluate(input, {true, false, false, false});
                double expected_power = 0.0, actual_power = 0.0;
                for (std::size_t row = 0; row < 32; ++row) {
                    double expected = 0.0, diagonal = 0.0;
                    if (row < 8) {
                        const auto& position = coordinates[row];
                        const double specific_heat = 3.0 + 0.02 * current[row] + 0.4 * time + 0.3 * position.x
                                                     + 0.2 * position.y + 0.1 * position.z;
                        const double increment = current[row] - previous[row];
                        expected = reference_mass[row] * specific_heat * increment / time_step;
                        diagonal = reference_mass[row] * (specific_heat + 0.02 * increment) / time_step;
                        expected_power += expected;
                        actual_power += active.residual[row] - baseline.residual[row];
                    }
                    residual_error = std::max({residual_error,
                        std::abs(active.residual[row] - baseline.residual[row] - expected),
                        std::abs(passive.residual[row] - passive_baseline.residual[row] - expected)});
                    for (std::size_t column = 0; column < 32; ++column) {
                        const double derivative =
                            active.jacobian[row * 32 + column] - baseline.jacobian[row * 32 + column];
                        jacobian_error =
                            std::max(jacobian_error, std::abs(derivative - (row == column ? diagonal : 0.0)));
                        if (column >= 8)
                            displacement_derivative = std::max(displacement_derivative, std::abs(derivative));
                    }
                }
                stored_power_error = std::max(stored_power_error, std::abs(actual_power - expected_power));
            }
    const std::string model = reduced ? "c3d8rt" : "c3d8t";
    std::cout << model << "_reference_mass_capacity_residual_error=" << residual_error << '\n'
              << model << "_reference_mass_capacity_jacobian_error=" << jacobian_error << '\n'
              << model << "_reference_mass_capacity_displacement_derivative=" << displacement_derivative << '\n'
              << model << "_reference_mass_stored_power_error=" << stored_power_error << '\n';
    return check(residual_error < 1.0e-10 && jacobian_error < 1.0e-10 && stored_power_error < 1.0e-10,
               model + " capacity uses initial density at time zero and reference nodes with current specific heat")
           && check(displacement_derivative == 0.0,
               model + " capacity has no displacement derivative under small or finite strain");
}

inline bool test_history_geometry(bool reduced) {
    const auto coordinates = unit_cube();
    const auto geometry = reduced ? fuelsim::elements::make_c3d8rt_geometry(coordinates)
                                  : fuelsim::elements::make_c3d8t_geometry(coordinates);
    fuelsim::Hex8LocalValues old{}, state{};
    const double angle = 0.35, c = std::cos(angle), s = std::sin(angle);
    for (std::size_t n = 0; n < 8; ++n) {
        const auto& x = coordinates[n];
        const double ox = 1.1 * x.x, oy = 0.9 * x.y, oz = 1.2 * x.z;
        const double nx = c * ox - s * oy, ny = s * ox + c * oy;
        old[8 + n] = ox - x.x;
        old[16 + n] = oy - x.y;
        old[24 + n] = oz - x.z;
        old[n] = 300;
        state[8 + n] = nx - x.x;
        state[16 + n] = ny - x.y;
        state[24 + n] = oz - x.z;
        state[n] = 300 + 2 * nx + 3 * ny - 4 * oz;
    }
    const fuelsim::IsotropicThermoelasticMaterial material(properties());
    const fuelsim::CartesianMaterialHistory history(reduced ? 1 : 8);
    const auto evaluate = reduced ? fuelsim::elements::evaluate_c3d8rt : fuelsim::elements::evaluate_c3d8t;
    const auto history_result = [&](const fuelsim::Hex8LocalValues& current, const fuelsim::Hex8LocalValues& previous) {
        const fuelsim::elements::C3d8Input
            input{material, geometry, current, previous, &history, 1.0, 1.0, 0.0, fuelsim::StrainFormulation::finite};
        return evaluate(input, {false, false, true, false});
    };
    const auto updated = history_result(state, old);
    const auto saved_old = old, saved_state = state;
    const auto diagnose = reduced ? fuelsim::test::recover_c3d8rt : fuelsim::test::recover_c3d8t;
    const auto result = diagnose(geometry, state, old, fuelsim::StrainFormulation::finite);
    bool passed =
        check(updated.history.size() == (reduced ? 1u : 8u) && near(updated.current_volume, 1.1 * 0.9 * 1.2, 2e-14)
                  && near(updated.committed_volume, 1.1 * 0.9 * 1.2, 2e-14),
            "history results report active material points and both physical volumes");
    const std::array<double, 9> rotation = {c, -s, 0, s, c, 0, 0, 0, 1};
    constexpr std::array<std::size_t, 8> gauss_to_node = {0, 1, 3, 2, 4, 5, 7, 6};
    double average_temperature = 0;
    for (std::size_t n = 0; n < 8; ++n)
        average_temperature += state[n] / 8;
    for (std::size_t q = 0; q < result.material_point_count; ++q) {
        const auto& point = result.points[q];
        for (std::size_t j = 0; j < 9; ++j)
            passed = check(near(updated.incremental_rotations[q][j], rotation[j], 2e-14)
                               && near(point.rotation[j], rotation[j], 2e-14),
                         "reconstructed rotation is incremental from the committed configuration")
                     && passed;
        for (std::size_t d = 0; d < 3; ++d) {
            double gradient = 0;
            for (std::size_t n = 0; n < 8; ++n)
                gradient += point.thermal_gradient[n][d] * state[n];
            passed = check(near(gradient, std::array<double, 3>{2, 3, -4}[d], 2e-12),
                         "reconstructed thermal gradient reconstructs an affine current-coordinate temperature")
                     && passed;
        }
        passed = check(near(point.temperature, reduced ? average_temperature : state[gauss_to_node[q]], 2e-12),
                     "reconstructed temperature follows the selected model rule")
                 && passed;
    }
    for (std::size_t q = result.material_point_count; q < 8; ++q) {
        passed = check(result.points[q].temperature == 0 && result.points[q].current_weighted_measure == 0,
                     "inactive reconstructed entries remain unused")
                 && passed;
    }
    auto expanded = state;
    for (std::size_t n = 0; n < 8; ++n) {
        expanded[8 + n] = 1.03 * (coordinates[n].x + state[8 + n]) - coordinates[n].x;
        expanded[16 + n] = 1.03 * (coordinates[n].y + state[16 + n]) - coordinates[n].y;
        expanded[24 + n] = 1.03 * (coordinates[n].z + state[24 + n]) - coordinates[n].z;
    }
    const auto volume = history_result(expanded, old);
    passed = check(near(volume.current_volume, 1.1 * 0.9 * 1.2 * 1.03 * 1.03 * 1.03, 2e-14)
                       && near(volume.committed_volume, 1.1 * 0.9 * 1.2, 2e-14),
                 "history results distinguish current and committed volumes")
             && passed;
    for (int invalid = 0; invalid < 3; ++invalid) {
        fuelsim::Hex8LocalValues bad{}, reference{};
        for (std::size_t n = 0; n < 8; ++n) {
            bad[n] = reference[n] = 300.0;
            bad[8 + n] = -2 * coordinates[n].x;
            if (invalid == 2)
                bad[16 + n] = -2 * coordinates[n].y;
        }
        bool rejected = false;
        try {
            (void)history_result(invalid == 1 ? reference : bad, invalid == 1 ? bad : reference);
        } catch (const std::domain_error&) {
            rejected = true;
        }
        passed =
            check(rejected, "history updates reject invalid current, committed and midpoint configurations") && passed;
    }
    passed = check(state == saved_state && old == saved_old, "history updates do not mutate borrowed states") && passed;
    return passed;
}

} // namespace fuelsim::test::c3d8
