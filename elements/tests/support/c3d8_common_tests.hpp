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
