#include "cax2t_gps.hpp"
#include "cax4t.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {
using namespace fuelsim;
using namespace fuelsim::elements;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double young = 2.0e11, poisson = 0.3, expansion = 1.0e-5;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

double error(double actual, double expected) {
    return std::abs(actual - expected) / std::max({1.0, std::abs(actual), std::abs(expected)});
}

IsotropicThermoelasticMaterial make_material(bool inelastic = false, bool varying = false) {
    const auto registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "cax2t_gps_test";
    functions->thermal = registry.bind_thermal("constant_thermophysical",
        {{"conductivity", 5.0}, {"density", 1000.0}, {"specific_heat", 500.0}});
    if (varying)
        functions->thermal.function = [](const ThermoelasticFunctionInput& input, ThermalPropertyOutput& output) {
            output.conductivity = 3.0 + 120.0 / input.temperature;
            output.density = 1000.0 + 10.0 * input.context.x + 20.0 * input.context.z + input.context.time
                             + 0.2 * (input.temperature - 600.0);
            output.specific_heat = 500.0 + 0.1 * (input.temperature - 600.0);
        };
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", young},
            {"poisson_ratio", poisson},
            {"reference_temperature", 600.0},
            {"young_modulus_temperature_coefficient", varying ? -8.0e7 : 0.0},
            {"poisson_ratio_temperature_coefficient", varying ? 2.0e-5 : 0.0}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", expansion},
            {"reference_temperature", 600.0},
            {"thermal_expansion_temperature_coefficient", varying ? 3.0e-9 : 0.0}}));
    if (inelastic) {
        functions->creep = registry.bind_creep("linear_temperature_norton",
            {{"coefficient", 1.0e-5},
                {"reference_stress", 1.0e8},
                {"stress_exponent", 3.0},
                {"reference_temperature", 600.0},
                {"coefficient_temperature_coefficient", varying ? 1.0e-7 : 0.0},
                {"reference_stress_temperature_coefficient", 0.0},
                {"stress_exponent_temperature_coefficient", varying ? 4.0e-3 : 0.0}});
        functions->plasticity = registry.bind_plasticity("linear_temperature_isotropic_hardening",
            {{"yield_stress", 2.0e8},
                {"hardening_modulus", 1.0e9},
                {"reference_temperature", 600.0},
                {"yield_stress_temperature_coefficient", varying ? 2.0e5 : 0.0},
                {"hardening_temperature_coefficient", varying ? 2.0e7 : 0.0}});
    }
    return IsotropicThermoelasticMaterial({functions, young});
}

bool same_history(const Cax2tGpsMaterialHistory& first, const Cax2tGpsMaterialHistory& second) {
    for (std::size_t q = 0; q < first.size(); ++q) {
        const auto& a = first[q];
        const auto& b = second[q];
        if (a.elastic_strain != b.elastic_strain || a.plastic_strain != b.plastic_strain
            || a.creep_strain != b.creep_strain || a.equivalent_plastic_strain != b.equivalent_plastic_strain
            || a.equivalent_creep_strain != b.equivalent_creep_strain || a.stress.rr != b.stress.rr
            || a.stress.zz != b.stress.zz || a.stress.hoop != b.stress.hoop || a.stress.rz != b.stress.rz)
            return false;
    }
    return true;
}

void check_analytic_mechanics() {
    const auto material = make_material();
    const auto solid = make_cax2t_gps_geometry({0.0, 2.0}, 1.0, 3.0);
    const double strain = 20.0 * expansion;
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues free{620.0, 620.0, 0.0, 2.0 * strain, 0.0, 2.0 * strain};
    const auto expanded = evaluate_cax2t_gps({material, solid, free, initial}, {true, true, true, true});
    for (const auto& stress : expanded.stress)
        require(std::max({std::abs(stress.rr), std::abs(stress.zz), std::abs(stress.hoop), std::abs(stress.rz)}) < 1e-7,
            "Uniform thermal expansion must be stress free, including an element touching the axis");
    for (const double residual : expanded.residual)
        require(std::abs(residual) < 1e-6, "Free thermal expansion must have zero residual");
    for (const auto& history : expanded.history)
        for (const double value : history.elastic_strain)
            require(std::abs(value) < 1e-15, "Steady elastic output must remove the thermal eigenstrain");

    const auto ring = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 3.0);
    const Cax2tGpsLocalValues stretched{600.0, 600.0, 0.0, 0.0, 0.0, 0.003};
    const auto response = evaluate_cax2t_gps({material, ring, stretched, initial}, {true, true, true, true});
    const double lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    const double shear = young / (2.0 * (1.0 + poisson));
    const double axial_stress = (lambda + 2.0 * shear) * 0.001;
    const double radial_stress = lambda * 0.001;
    const double axial_force = axial_stress * 3.0 * pi;
    for (const auto& history : response.history)
        require(history.elastic_strain == std::array<double, 4>{{0.0, 0.001, 0.0, 0.0}},
            "Steady axial loading must return its nonzero elastic history for result output");
    auto finite_input = Cax2tGpsInput{material, ring, stretched, initial};
    finite_input.strain_formulation = StrainFormulation::finite;
    const auto finite_response = evaluate_cax2t_gps(finite_input, {true, true, true, true});
    for (const auto& history : finite_response.history)
        require(std::abs(history.elastic_strain[1] - 0.001 / 1.0005) < 1e-15,
            "Steady finite elastic output must retain the midpoint strain rather than a zero default");
    for (const auto& stress : response.stress)
        require(error(stress.zz, axial_stress) < 1e-14 && error(stress.rr, radial_stress) < 1e-14
                    && error(stress.hoop, radial_stress) < 1e-14 && stress.rz == 0.0,
            "Constant generalized axial strain must produce the analytical constrained stresses without shear");
    require(error(response.residual[4], -axial_force) < 1e-14 && error(response.residual[5], axial_force) < 1e-14,
        "The end sections must receive opposite complete cross-section axial forces");
    require(error(response.residual[2], -2.0 * pi * 1.0 * 3.0 * radial_stress) < 1e-14
                && error(response.residual[3], 2.0 * pi * 2.0 * 3.0 * radial_stress) < 1e-14,
        "Radial force must include the complete circumferential and axial measure");
    auto translated = stretched;
    translated[4] += 0.125;
    translated[5] += 0.125;
    const auto moved = evaluate_cax2t_gps({material, ring, translated, initial}, {true, true, true, true});
    for (std::size_t i = 0; i < stretched.size(); ++i) {
        require(error(response.residual[i], moved.residual[i]) < 1e-13,
            "A shared axial translation must not change the strain or internal force");
        require(response.jacobian[6 * i + 4] == -response.jacobian[6 * i + 5],
            "Rigid axial translation must be an exact null direction of the body tangent");
    }
}

void check_thermal_operators() {
    const auto material = make_material();
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, -1.0, 1.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues heated{620.0, 620.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsMaterialHistory history{};
    const Cax2tGpsInput
        input{material, geometry, heated, initial, &history, 0.5, 0.5, 2.0e7, StrainFormulation::small, true};
    const auto response = evaluate_cax2t_gps(input, {true, true, true, true});
    const double volume = 6.0 * pi;
    require(error(response.stored_heat_rate, volume * 2.0e7) < 1e-14
                && error(response.generated_heat_rate, volume * 2.0e7) < 1e-14,
        "Uniform heating must reproduce the independent cylindrical volume and heat rate");
    require(std::abs(response.residual[0]) < 1e-7 && std::abs(response.residual[1]) < 1e-7,
        "Uniform heating must balance a uniform source at both radial nodes");
    auto stationary_input = input;
    stationary_input.include_thermal_time_term = false;
    const auto stationary = evaluate_cax2t_gps(stationary_input, {true, true, false, false});
    const std::array<double, 2> nodal_volumes = {8.0 * pi / 3.0, 10.0 * pi / 3.0};
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t column = 0; column < 2; ++column) {
            const double capacity = row == column ? nodal_volumes[row] * 5.0e5 / 0.5 : 0.0;
            const double conduction = (row == column ? 1.0 : -1.0) * pi * 2.0 * 5.0 * 3.0;
            require(error(response.jacobian[6 * row + column] - stationary.jacobian[6 * row + column], capacity)
                        < 1e-14,
                "Lumped heat capacity must use cylindrical row-sum weights and an exactly diagonal temperature block");
            require(error(stationary.jacobian[6 * row + column], conduction) < 1e-14,
                "Radial conduction must match the analytical cylindrical stiffness");
        }
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t column = 2; column < 6; ++column)
            require(response.jacobian[6 * row + column] == 0.0,
                "Small-strain reference thermal operators must be independent of displacement");

    const Cax2tGpsLocalValues nonuniform{610.0, 650.0, 0.0, 0.0, 0.0, 0.0};
    const auto varying = evaluate_cax2t_gps({material, geometry, nonuniform, initial}, {true, true, false, true});
    const double stress = -young * expansion * 30.0 / (1.0 - 2.0 * poisson);
    for (const auto& point : varying.stress) {
        require(error(point.zz, stress) < 1e-13,
            "Thermal expansion must use the arithmetic mean of the two radial nodal temperatures");
    }
    const auto variable_material = make_material(false, true);
    const Cax2tGpsInput variable_input{variable_material,
        geometry,
        nonuniform,
        initial,
        &history,
        0.5,
        0.5,
        0.0,
        StrainFormulation::small,
        true};
    const auto nonlinear_thermal = evaluate_cax2t_gps(variable_input, {true, true, true, true});
    auto conduction_input = variable_input;
    conduction_input.include_thermal_time_term = false;
    const auto nonlinear_conduction = evaluate_cax2t_gps(conduction_input, {true, true, false, false});
    std::array<double, 2> expected_heat{}, interpolated_heat{}, conductivity_derivative{};
    double integrated_conductivity = 0.0;
    for (std::size_t q = 0; q < 2; ++q) {
        const double station = (q == 0 ? -1.0 : 1.0) / std::sqrt(3.0);
        const std::array<double, 2> shape = {0.5 * (1.0 - station), 0.5 * (1.0 + station)};
        const double radius = 1.5 + 0.5 * station;
        const double temperature = shape[0] * nonuniform[0] + shape[1] * nonuniform[1];
        const double density = 1000.0 + 10.0 * radius + 0.5 + 0.2 * (temperature - 600.0);
        const double temperature_rate = (temperature - 600.0) / 0.5;
        const double measure = 2.0 * pi * radius;
        integrated_conductivity += measure * (3.0 + 120.0 / nonuniform[q]);
        conductivity_derivative[q] = -measure * 120.0 / (nonuniform[q] * nonuniform[q]);
        for (std::size_t node = 0; node < 2; ++node) {
            const double gradient = node == 0 ? -1.0 : 1.0;
            interpolated_heat[node] +=
                measure
                * ((3.0 + 120.0 / temperature) * gradient * 40.0
                    + shape[node] * density * (500.0 + 0.1 * (temperature - 600.0)) * temperature_rate);
        }
    }
    double expected_storage = 0.0;
    for (std::size_t node = 0; node < 2; ++node) {
        const double change = nonuniform[node] - 600.0;
        const double density = 1000.0 + 10.0 * geometry.radii[node] + 0.5 + 0.2 * change;
        const double heat_capacity = 500.0 + 0.1 * change;
        const double storage = nodal_volumes[node] * density * heat_capacity * change / 0.5;
        expected_storage += storage;
        const double sign = node == 0 ? -1.0 : 1.0;
        expected_heat[node] = sign * integrated_conductivity * 40.0 + storage;
        require(error(nonlinear_thermal.residual[node], expected_heat[node]) < 1e-13,
            "Small-strain conductivity, density, specific heat and temperature rate must use paired nodes");
        require(std::abs(expected_heat[node] - interpolated_heat[node]) > 1e3,
            "Nodal heat operators must be distinguished from the former interpolated consistent operator");
        for (std::size_t column = 0; column < 2; ++column) {
            const double conduction =
                sign * ((column == 0 ? -1.0 : 1.0) * integrated_conductivity + 40.0 * conductivity_derivative[column]);
            const double capacity =
                node == column ? nodal_volumes[node]
                                     * ((0.2 * heat_capacity + 0.1 * density) * change + density * heat_capacity) / 0.5
                               : 0.0;
            require(error(nonlinear_conduction.jacobian[6 * node + column], conduction) < 1e-13,
                "Each nodal conductivity derivative must enter the appropriate temperature column");
            require(
                error(nonlinear_thermal.jacobian[6 * node + column] - nonlinear_conduction.jacobian[6 * node + column],
                    capacity)
                    < 1e-13,
                "Temperature-dependent lumped capacity must remain diagonal and differentiate both rho and cp");
        }
        for (std::size_t column = 2; column < nonuniform.size(); ++column)
            require(nonlinear_thermal.jacobian[6 * node + column] == 0.0,
                "Nonuniform small-strain thermal operators must have exactly zero displacement derivatives");
    }
    require(error(nonlinear_thermal.stored_heat_rate, expected_storage) < 1e-13
                && error(nonlinear_thermal.residual[0] + nonlinear_thermal.residual[1], expected_storage) < 1e-13,
        "Nodal heat storage diagnostics must equal the complete summed heat residual");
}

void check_mean_thermal_expansion(StrainFormulation formulation) {
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues first{610.0, 630.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues state{620.0, 700.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsMaterialHistory initial_history{};
    double maximum_error = 0.0;
    for (const bool temperature_dependent : {false, true}) {
        const auto material = make_material(false, temperature_dependent);
        const auto first_history =
            evaluate_cax2t_gps({material, geometry, first, initial, &initial_history, 0.1, 0.1, 0.0, formulation})
                .history;
        const Cax2tGpsInput input{material, geometry, state, first, &first_history, 0.1, 0.2, 0.0, formulation};
        const auto response = evaluate_cax2t_gps(input, {true, true, true, true});
        const auto steady =
            evaluate_cax2t_gps({material, geometry, state, initial, nullptr, 0.0, 0.2, 0.0, formulation},
                {true, true, true, true});
        auto redistributed = state;
        redistributed[0] += 12.0;
        redistributed[1] -= 12.0;
        const auto same_mean =
            evaluate_cax2t_gps({material, geometry, redistributed, first, &first_history, 0.1, 0.2, 0.0, formulation},
                {true, true, true, true});
        const double mean_change = 60.0;
        const double eigenstrain = (expansion + (temperature_dependent ? 3.0e-9 * mean_change : 0.0)) * mean_change;
        for (std::size_t q = 0; q < response.stress.size(); ++q) {
            for (const bool perturb_mean_preserving : {false, true}) {
                const auto& checked = perturb_mean_preserving ? same_mean : response;
                const double material_change = (perturb_mean_preserving ? redistributed[q] : state[q]) - 600.0;
                const double modulus = young + (temperature_dependent ? -8.0e7 * material_change : 0.0);
                const double ratio = poisson + (temperature_dependent ? 2.0e-5 * material_change : 0.0);
                const double expected = -modulus * eigenstrain / (1.0 - 2.0 * ratio);
                const auto& stress = checked.stress[q];
                require(error(stress.rr, expected) < 1e-13 && error(stress.zz, expected) < 1e-13
                            && error(stress.hoop, expected) < 1e-13 && stress.rz == 0.0,
                    "Mean thermal expansion must coexist with paired-endpoint elastic properties");
                for (std::size_t component = 0; component < 3; ++component)
                    require(std::abs(checked.history[q].elastic_strain[component] + eigenstrain) < 1e-15,
                        "Nonuniform temperatures with the same arithmetic mean must give the same thermal strain");
            }
            require(error(steady.stress[q].rr, response.stress[q].rr) < 1e-13
                        && error(steady.stress[q].zz, response.stress[q].zz) < 1e-13
                        && error(steady.stress[q].hoop, response.stress[q].hoop) < 1e-13,
                "Current and committed mean temperatures must recover the same constrained elastic stress as steady "
                "loading");
            for (std::size_t component = 0; component < 3; ++component)
                require(std::abs(steady.history[q].elastic_strain[component] + eigenstrain) < 1e-15,
                    "Steady elastic history must remove the mean-temperature eigenstrain");
        }
        if (temperature_dependent) {
            require(std::abs(response.stress[0].rr - response.stress[1].rr) > 1.0e6,
                "A mean temperature for expansion must not replace the local temperature for elasticity");
        } else {
            for (std::size_t row = 2; row < state.size(); ++row) {
                require(error(response.residual[row], same_mean.residual[row]) < 1e-14,
                    "A mean-preserving temperature redistribution must leave constant-material mechanical forces "
                    "unchanged");
                require(error(response.jacobian[6 * row], response.jacobian[6 * row + 1]) < 1e-14,
                    "The mean thermal strain must contribute equal derivatives from both radial nodal temperatures");
            }
        }
        const std::array<std::array<double, 2>, 3> directions = {{{1.0, 0.0}, {0.0, 1.0}, {1.0, -1.0}}};
        for (const auto& direction : directions) {
            constexpr double step = 1.0e-3;
            auto plus = state, minus = state;
            for (std::size_t column = 0; column < 2; ++column) {
                plus[column] += step * direction[column];
                minus[column] -= step * direction[column];
            }
            const auto rp =
                evaluate_cax2t_gps({material, geometry, plus, first, &first_history, 0.1, 0.2, 0.0, formulation});
            const auto rm =
                evaluate_cax2t_gps({material, geometry, minus, first, &first_history, 0.1, 0.2, 0.0, formulation});
            for (std::size_t row = 0; row < state.size(); ++row) {
                const double tangent =
                    response.jacobian[6 * row] * direction[0] + response.jacobian[6 * row + 1] * direction[1];
                maximum_error =
                    std::max(maximum_error, error(tangent, (rp.residual[row] - rm.residual[row]) / (2.0 * step)));
            }
        }
    }
    require(maximum_error < 2e-7,
        "Mean expansion and local elastic/thermal property derivatives must match independent centered differences");
    std::cout << "cax2t_gps_" << (formulation == StrainFormulation::finite ? "finite" : "small")
              << "_mean_temperature_jacobian_error=" << maximum_error << '\n';
}

double small_mean_hoop_energy(const Cax2tGpsLocalValues& state) {
    // Integrate the constant assumed strain over the independently known annulus.
    const double thermal = expansion * (0.5 * (state[0] + state[1]) - 600.0);
    const std::array<double, 3> elastic = {state[3] - state[2] - thermal,
        0.5 * (state[5] - state[4]) - thermal,
        (state[2] + state[3]) / 3.0 - thermal};
    const double trace = elastic[0] + elastic[1] + elastic[2];
    const double lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    const double shear = young / (2.0 * (1.0 + poisson));
    return 6.0 * pi
           * (0.5 * lambda * trace * trace
               + shear * (elastic[0] * elastic[0] + elastic[1] * elastic[1] + elastic[2] * elastic[2]));
}

void check_small_mean_hoop_mechanics() {
    const auto material = make_material();
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues state{620.0, 700.0, 0.003, 0.008, 0.002, 0.014};
    const auto response = evaluate_cax2t_gps({material, geometry, state, initial}, {true, true, true, true});
    const auto passive = evaluate_cax2t_gps({material, geometry, state, initial});
    require(response.residual == passive.residual && same_history(response.history, passive.history),
        "Small mean-hoop mechanics must preserve residual-only and tangent-call agreement");
    const double thermal = expansion * 60.0;
    const std::array<double, 3> elastic = {0.005 - thermal, 0.006 - thermal, 0.011 / 3.0 - thermal};
    const double lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    const double shear = young / (2.0 * (1.0 + poisson));
    const double pressure = lambda * (elastic[0] + elastic[1] + elastic[2]);
    const std::array<double, 3> stress = {pressure + 2.0 * shear * elastic[0],
        pressure + 2.0 * shear * elastic[1],
        pressure + 2.0 * shear * elastic[2]};
    const std::array<double, 4> force = {6.0 * pi * (-stress[0] + stress[2] / 3.0),
        6.0 * pi * (stress[0] + stress[2] / 3.0),
        -3.0 * pi * stress[1],
        3.0 * pi * stress[1]};
    double integrated_energy = 0.0;
    for (std::size_t q = 0; q < response.stress.size(); ++q) {
        const auto& actual = response.stress[q];
        require(error(actual.rr, stress[0]) < 1e-13 && error(actual.zz, stress[1]) < 1e-13
                    && error(actual.hoop, stress[2]) < 1e-13 && actual.rz == 0.0,
            "Nonuniform radial motion must use one volume-averaged hoop strain at both material points");
        for (std::size_t component = 0; component < elastic.size(); ++component)
            require(std::abs(response.history[q].elastic_strain[component] - elastic[component]) < 1e-15,
                "Material history must contain the assumed hoop strain, not an averaged output stress");
        integrated_energy += 1.5 * pi * (actual.rr * elastic[0] + actual.zz * elastic[1] + actual.hoop * elastic[2]);
    }
    require(error(integrated_energy, small_mean_hoop_energy(state)) < 1e-13,
        "Assumed-strain elastic energy must not retain the pointwise hoop-strain variance energy");
    for (std::size_t row = 2; row < state.size(); ++row) {
        require(error(response.residual[row], force[row - 2]) < 1e-13,
            "Both radial reactions and axial end forces must follow the assumed-strain virtual work");
        auto plus = state, minus = state;
        constexpr double step = 1e-7;
        plus[row] += step;
        minus[row] -= step;
        const double energy_derivative = (small_mean_hoop_energy(plus) - small_mean_hoop_energy(minus)) / (2.0 * step);
        require(error(response.residual[row], energy_derivative) < 2e-9,
            "Every mechanical reaction must equal an independent energy derivative at fixed nodal temperature");
        for (std::size_t column = 2; column < state.size(); ++column)
            require(error(response.jacobian[6 * row + column], response.jacobian[6 * column + row]) < 1e-14,
                "The elastic small-strain mechanical tangent must remain symmetric");
    }
    double maximum_error = 0.0;
    for (std::size_t column = 0; column < state.size(); ++column) {
        const double step = column < 2 ? 1e-3 : 1e-7;
        auto plus = state, minus = state;
        plus[column] += step;
        minus[column] -= step;
        const auto rp = evaluate_cax2t_gps({material, geometry, plus, initial});
        const auto rm = evaluate_cax2t_gps({material, geometry, minus, initial});
        for (std::size_t row = 0; row < state.size(); ++row)
            maximum_error = std::max(maximum_error,
                error(response.jacobian[6 * row + column], (rp.residual[row] - rm.residual[row]) / (2.0 * step)));
    }
    require(maximum_error < 2e-7,
        "Every small-strain tangent column must differentiate the averaged hoop strain and virtual work");
    std::cout << "cax2t_gps_small_mean_hoop_jacobian_error=" << maximum_error << '\n';
}

void check_finite_mean_hoop_history_and_geometry() {
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues first{620.0, 660.0, 0.1, 0.15, 0.02, 0.10};
    const Cax2tGpsLocalValues state{650.0, 730.0, 0.2, 0.27, 0.04, 0.22};
    const Cax2tGpsMaterialHistory initial_history{};
    const std::array<double, 3> first_stretch = {1.05, 1.04, 1.0 + 0.25 / 3.0};
    const std::array<double, 3> current_stretch = {1.07, 1.09, 1.0 + 0.47 / 3.0};
    double maximum_error = 0.0;
    for (const bool temperature_dependent : {false, true}) {
        const auto material = make_material(false, temperature_dependent);
        const auto first_response = evaluate_cax2t_gps(
            {material, geometry, first, initial, &initial_history, 0.1, 0.1, 0.0, StrainFormulation::finite},
            {true, true, true, true});
        const auto history = first_response.history;
        const auto saved_history = history;
        const Cax2tGpsInput
            input{material, geometry, state, first, &history, 0.1, 0.2, 2.0e6, StrainFormulation::finite, true};
        const auto response = evaluate_cax2t_gps(input, {true, true, true, true});
        const auto passive = evaluate_cax2t_gps(input);
        require(response.residual == passive.residual && same_history(response.history, passive.history),
            "Nonuniform finite hoop averaging must preserve exact residual-only and tangent-call agreement");
        Cax2tGpsLocalResidual expected{};
        double stored = 0.0, generated = 0.0, wrong_mechanical_force = 0.0, wrong_thermal_residual = 0.0;
        std::array<double, 2> nodal_volumes{}, expected_conduction{};
        const double midpoint_span = 0.5 * (first_stretch[0] + current_stretch[0]);
        for (std::size_t q = 0; q < response.stress.size(); ++q) {
            const double station = (q == 0 ? -1.0 : 1.0) / std::sqrt(3.0);
            const std::array<double, 2> shape = {0.5 * (1.0 - station), 0.5 * (1.0 + station)};
            const double reference_radius = 1.5 + 0.5 * station;
            const double material_change = state[q] - 600.0;
            const double modulus = young + (temperature_dependent ? -8.0e7 * material_change : 0.0);
            const double ratio = poisson + (temperature_dependent ? 2.0e-5 * material_change : 0.0);
            const double thermal = (expansion + (temperature_dependent ? 3.0e-9 * 90.0 : 0.0)) * 90.0;
            const double old_thermal = (expansion + (temperature_dependent ? 3.0e-9 * 40.0 : 0.0)) * 40.0;
            std::array<double, 3> elastic{};
            for (std::size_t component = 0; component < elastic.size(); ++component) {
                const double first_increment =
                    2.0 * (first_stretch[component] - 1.0) / (first_stretch[component] + 1.0);
                elastic[component] = first_increment
                                     + 2.0 * (current_stretch[component] - first_stretch[component])
                                           / (current_stretch[component] + first_stretch[component])
                                     - thermal;
                require(std::abs(history[q].elastic_strain[component] - (first_increment - old_thermal)) < 1e-14
                            && std::abs(response.history[q].elastic_strain[component] - elastic[component]) < 1e-14,
                    "Two accepted nonuniform finite increments must accumulate midpoint increments of mean hoop "
                    "stretch");
            }
            const double lambda = modulus * ratio / ((1.0 + ratio) * (1.0 - 2.0 * ratio));
            const double shear = modulus / (2.0 * (1.0 + ratio));
            const double pressure = lambda * (elastic[0] + elastic[1] + elastic[2]);
            const std::array<double, 3> stress = {pressure + 2.0 * shear * elastic[0],
                pressure + 2.0 * shear * elastic[1],
                pressure + 2.0 * shear * elastic[2]};
            require(error(response.stress[q].rr, stress[0]) < 1e-13 && error(response.stress[q].zz, stress[1]) < 1e-13
                        && error(response.stress[q].hoop, stress[2]) < 1e-13 && response.stress[q].rz == 0.0,
                "Finite mean hoop stretch must coexist with endpoint-temperature elastic properties and mean "
                "expansion");
            const double mechanical_measure =
                2.0 * pi * reference_radius * current_stretch[0] * current_stretch[1] * current_stretch[2];
            // Capacity and source use the actual current annulus. Conduction
            // uses its whole-volume ratio and the incremental midpoint gradient.
            const double current_radius = shape[0] * 1.2 + shape[1] * 2.27;
            const double thermal_measure = pi * current_radius * (2.27 - 1.2) * (2.22 - 0.04);
            const double conductivity = temperature_dependent ? 3.0 + 120.0 / state[q] : 5.0;
            generated += thermal_measure * 2.0e6;
            for (std::size_t node = 0; node < 2; ++node) {
                const double radial_gradient = (node == 0 ? -1.0 : 1.0) / (2.27 - 1.2);
                const double force_density = radial_gradient * stress[0] + stress[2] / (1.2 + 2.27);
                expected[2 + node] += mechanical_measure * force_density;
                expected[4 + node] += mechanical_measure * (node == 0 ? -1.0 : 1.0) * stress[1] / (2.22 - 0.04);
                const double change = state[node] - 600.0;
                const double density =
                    temperature_dependent ? 1000.0 + 10.0 * geometry.radii[node] + 20.0 + 0.2 + 0.2 * change : 1000.0;
                const double capacity = density * (temperature_dependent ? 500.0 + 0.1 * change : 500.0);
                const double storage = capacity * (state[node] - first[node]) / 0.1;
                const double conduction = mechanical_measure * conductivity * (node == 0 ? -1.0 : 1.0) * (730.0 - 650.0)
                                          / (midpoint_span * midpoint_span);
                expected_conduction[node] += conduction;
                nodal_volumes[node] += thermal_measure * shape[node];
                stored += thermal_measure * shape[node] * storage;
                expected[node] += conduction + thermal_measure * shape[node] * (storage - 2.0e6);
                if (node == 0) {
                    wrong_mechanical_force += thermal_measure * force_density;
                    wrong_thermal_residual += conduction + mechanical_measure * shape[node] * (storage - 2.0e6);
                }
            }
        }
        for (std::size_t row = 0; row < state.size(); ++row)
            require(error(response.residual[row], expected[row]) < 1e-13,
                "Finite conduction must use midpoint gradients and whole-volume weights, with nodal current-volume "
                "storage and source");
        require(error(response.stored_heat_rate, stored) < 1e-13
                    && error(response.generated_heat_rate, generated) < 1e-13,
            "Finite lumped storage and source diagnostics must retain actual current-volume nodal weights");
        require(error(response.residual[0] + response.residual[1], stored - generated) < 1e-13,
            "Finite nodal heat residuals must conserve stored minus generated heat");
        require(std::abs(expected[0] - wrong_thermal_residual) > 1e6,
            "The thermal geometry regression must distinguish actual point volumes from averaged mechanical weights");
        if (temperature_dependent)
            require(std::abs(expected[2] - wrong_mechanical_force) > 1e6,
                "Point-dependent elasticity must distinguish averaged mechanical measure from thermal point volumes");
        auto stationary_input = input;
        stationary_input.include_thermal_time_term = false;
        const auto stationary = evaluate_cax2t_gps(stationary_input, {true, true, false, false});
        auto conduction_input = stationary_input;
        conduction_input.volumetric_heat_source = 0.0;
        const auto conduction_only = evaluate_cax2t_gps(conduction_input, {true, true, false, false});
        double capacity_geometry = 0.0, conduction_geometry = 0.0;
        for (std::size_t row = 0; row < 2; ++row) {
            require(error(conduction_only.residual[row], expected_conduction[row]) < 1e-13,
                "Finite conduction alone must use paired conductivity, midpoint gradients and whole-volume weights");
            const double change = state[row] - 600.0;
            const double density =
                temperature_dependent ? 1000.0 + 10.0 * geometry.radii[row] + 20.0 + 0.2 + 0.2 * change : 1000.0;
            const double cp = temperature_dependent ? 500.0 + 0.1 * change : 500.0;
            const double slope = temperature_dependent ? 0.2 * cp + 0.1 * density : 0.0;
            for (std::size_t column = 0; column < 2; ++column) {
                const double capacity =
                    row == column ? nodal_volumes[row] * (density * cp + slope * (state[row] - first[row])) / 0.1 : 0.0;
                require(error(response.jacobian[6 * row + column] - stationary.jacobian[6 * row + column], capacity)
                            < 1e-13,
                    "Finite lumped capacity must retain its nodal temperature diagonal");
            }
            for (std::size_t column = 2; column < state.size(); ++column) {
                capacity_geometry = std::max(capacity_geometry,
                    std::abs(response.jacobian[6 * row + column] - stationary.jacobian[6 * row + column]));
                conduction_geometry =
                    std::max(conduction_geometry, std::abs(conduction_only.jacobian[6 * row + column]));
            }
        }
        require(capacity_geometry > 1.0,
            "Current nodal heat-capacity weights must have nonzero displacement derivatives");
        require(conduction_geometry > 1.0
                    && std::abs(conduction_only.residual[0] + conduction_only.residual[1]) < 1e-10,
            "Finite conduction must preserve heat and have nonzero midpoint and volume geometry derivatives");
        for (std::size_t column = 0; column < state.size(); ++column) {
            const double step = column < 2 ? 1e-3 : 1e-6;
            auto plus = state, minus = state;
            plus[column] += step;
            minus[column] -= step;
            const auto rp = evaluate_cax2t_gps(
                {material, geometry, plus, first, &history, 0.1, 0.2, 2.0e6, StrainFormulation::finite, true});
            const auto rm = evaluate_cax2t_gps(
                {material, geometry, minus, first, &history, 0.1, 0.2, 2.0e6, StrainFormulation::finite, true});
            // Isolate the small off-diagonal heat derivative from nodal storage
            // and the temperature-independent source. The exact zero storage
            // entries were checked above; all 36 full tangent entries remain checked.
            Cax2tGpsLocalResidual stationary_plus{}, stationary_minus{};
            if (column < 2) {
                stationary_plus = evaluate_cax2t_gps(
                    {material, geometry, plus, first, &history, 0.1, 0.2, 0.0, StrainFormulation::finite})
                                      .residual;
                stationary_minus = evaluate_cax2t_gps(
                    {material, geometry, minus, first, &history, 0.1, 0.2, 0.0, StrainFormulation::finite})
                                       .residual;
            }
            for (std::size_t row = 0; row < state.size(); ++row) {
                const bool isolated = row < 2 && column < 2 && row != column;
                const double numerical = isolated ? (stationary_plus[row] - stationary_minus[row]) / (2.0 * step)
                                                  : (rp.residual[row] - rm.residual[row]) / (2.0 * step);
                maximum_error = std::max(maximum_error, error(response.jacobian[6 * row + column], numerical));
            }
        }
        require(same_history(history, saved_history),
            "All finite mean-hoop trial evaluations must leave the accepted material histories unchanged");
    }
    require(maximum_error < 2e-7,
        "Every finite tangent column must differentiate mean hoop mechanics and actual thermal geometry consistently");
    std::cout << "cax2t_gps_finite_mean_hoop_jacobian_error=" << maximum_error << '\n';
}

void check_cax4t_projection(StrainFormulation formulation) {
    // The complete radial thermal/mechanical fields expand to the rectangle.
    // Use identical radial material fields without an axial-coordinate dependence.
    auto functions = std::make_shared<MaterialFunctionSet>(make_material(false, true).functions());
    functions->name = "gps_cax4t_complete_projection";
    functions->thermal.function = [](const ThermoelasticFunctionInput& input, ThermalPropertyOutput& output) {
        output.conductivity = 3.0 + 120.0 / input.temperature;
        output.density = 1000.0 + 10.0 * input.context.x + input.context.time + 0.2 * (input.temperature - 600.0);
        output.specific_heat = 500.0 + 0.1 * (input.temperature - 600.0);
    };
    const IsotropicThermoelasticMaterial material({functions, young});
    const auto radial_geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const auto quad_geometry = make_cax4t_geometry({{{1.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {1.0, 2.0}}});
    Cax2tGpsLocalValues old{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const std::array<Cax2tGpsLocalValues, 2> states = {
        {{620.0, 660.0, 0.1, 0.15, 0.02, 0.10}, {650.0, 730.0, 0.2, 0.27, 0.04, 0.22}}};
    // P expands the six radial/end-section values into the twelve rectangular CAX4T values.
    const std::array<std::size_t, 12> radial_dof = {0, 1, 1, 0, 2, 3, 3, 2, 4, 4, 5, 5};
    const std::array<std::size_t, 4> radial_point = {0, 1, 1, 0};
    Cax2tGpsMaterialHistory radial_history{};
    Quad4MaterialHistory quad_history{};
    double maximum_error = 0.0;
    for (std::size_t increment = 0; increment < states.size(); ++increment) {
        const auto& state = states[increment];
        Cax4LocalValues quad_old{}, quad_state{};
        for (std::size_t i = 0; i < radial_dof.size(); ++i) {
            quad_old[i] = old[radial_dof[i]];
            quad_state[i] = state[radial_dof[i]];
        }
        const double time = 0.1 * static_cast<double>(increment + 1);
        const auto radial = evaluate_cax2t_gps(
            {material, radial_geometry, state, old, &radial_history, 0.1, time, 2e6, formulation, true},
            {true, true, true, true});
        const auto quad = evaluate_cax4t(
            {material, quad_geometry, quad_state, quad_old, &quad_history, 0.1, time, 2e6, formulation, true},
            {true, true, true, true});
        Cax2tGpsLocalResidual projected_residual{};
        Cax2tGpsLocalJacobian projected_tangent{};
        for (std::size_t row = 0; row < radial_dof.size(); ++row) {
            projected_residual[radial_dof[row]] += quad.residual[row];
            for (std::size_t column = 0; column < radial_dof.size(); ++column)
                projected_tangent[6 * radial_dof[row] + radial_dof[column]] += quad.jacobian[12 * row + column];
        }
        for (std::size_t q = 0; q < radial_point.size(); ++q) {
            const auto& actual = radial.stress[radial_point[q]];
            const auto& expected = quad.stress[q];
            maximum_error = std::max({maximum_error,
                error(actual.rr, expected.rr),
                error(actual.zz, expected.zz),
                error(actual.hoop, expected.hoop)});
            for (std::size_t component = 0; component < 4; ++component)
                require(std::abs(radial.history[radial_point[q]].elastic_strain[component]
                                 - quad.history[q].elastic_strain[component])
                            < 1e-14,
                    "Reduced CAX4T and radial elements must accumulate identical material strain history");
        }
        for (std::size_t row = 0; row < state.size(); ++row) {
            maximum_error = std::max(maximum_error, error(radial.residual[row], projected_residual[row]));
            for (std::size_t column = 0; column < state.size(); ++column)
                maximum_error = std::max(maximum_error,
                    error(radial.jacobian[6 * row + column], projected_tangent[6 * row + column]));
        }
        require(error(radial.stored_heat_rate, quad.stored_heat_rate) < 2e-11
                    && error(radial.generated_heat_rate, quad.generated_heat_rate) < 2e-11,
            "Both temperature increments must preserve CAX4T heat-storage and source diagnostics");
        radial_history = radial.history;
        quad_history = quad.history;
        old = state;
    }
    require(maximum_error < 2e-11,
        "Complete nonuniform CAX4T thermal/mechanical residual and 6x6 tangent must reduce to the radial element");
    std::cout << "cax2t_gps_" << (formulation == StrainFormulation::finite ? "finite" : "small")
              << "_cax4t_complete_projection_error=" << maximum_error << '\n';
}

void check_paired_material_temperatures(StrainFormulation formulation, bool inelastic) {
    const auto material = make_material(inelastic, true);
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const std::array<std::array<double, 2>, 2> temperatures = {{{620.0, 740.0}, {730.0, 650.0}}};
    Cax2tGpsLocalValues old{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    Cax2tGpsMaterialHistory history{};
    double previous_radial_stretch = 1.0, previous_axial_stretch = 1.0;
    double accumulated_trace = 0.0, maximum_error = 0.0;
    constexpr double dt = 0.25;
    const auto relative = [](double actual, double expected) {
        return std::abs(actual / expected - 1.0);
    };
    for (std::size_t increment = 0; increment < temperatures.size(); ++increment) {
        const double radial_strain = 0.0001 * static_cast<double>(increment + 1);
        const double axial_strain = 0.004 * static_cast<double>(increment + 1);
        const Cax2tGpsLocalValues state{temperatures[increment][0],
            temperatures[increment][1],
            radial_strain,
            2.0 * radial_strain,
            0.0,
            2.0 * axial_strain};
        const double radial_stretch = 1.0 + radial_strain, axial_stretch = 1.0 + axial_strain;
        if (formulation == StrainFormulation::finite)
            accumulated_trace +=
                4.0 * (radial_stretch - previous_radial_stretch) / (radial_stretch + previous_radial_stretch)
                + 2.0 * (axial_stretch - previous_axial_stretch) / (axial_stretch + previous_axial_stretch);
        else
            accumulated_trace = 2.0 * radial_strain + axial_strain;
        const auto saved_history = history;
        const auto saved_state = state, saved_old = old;
        // Separate analytical tests cover lumped heat storage.  This test
        // resolves the small conduction terms alongside material derivatives.
        const Cax2tGpsInput input{material,
            geometry,
            state,
            old,
            &history,
            dt,
            dt * static_cast<double>(increment + 1),
            0.0,
            formulation};
        const auto active = evaluate_cax2t_gps(input, {true, true, true, true});
        const auto passive = evaluate_cax2t_gps(input);
        require(active.residual == passive.residual && same_history(active.history, passive.history),
            "Paired-temperature material paths must return identical residual and trial history");
        const double mean_change = 0.5 * (state[0] + state[1]) - 600.0;
        const double thermal_strain = (expansion + 3.0e-9 * mean_change) * mean_change;
        for (std::size_t q = 0; q < 2; ++q) {
            const auto& point = active.history[q];
            const auto& strain = point.elastic_strain;
            const auto& stress = point.stress;
            const double trace = strain[0] + strain[1] + strain[2];
            const double shear = (stress.zz - stress.rr) / (2.0 * (strain[1] - strain[0]));
            const double bulk = (stress.rr + stress.zz + stress.hoop) / (3.0 * trace);
            const double inferred_young = 9.0 * bulk * shear / (3.0 * bulk + shear);
            const double inferred_poisson = (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear));
            const double change = state[q] - 600.0;
            const double station = (q == 0 ? -1.0 : 1.0) / std::sqrt(3.0);
            const double interpolated_change = 0.5 * ((1.0 - station) * state[0] + (1.0 + station) * state[1]) - 600.0;
            require(relative(inferred_young, young - 8.0e7 * change) < 1e-11,
                "GPS Young modulus must use the paired radial endpoint temperature");
            require(relative(inferred_poisson, poisson + 2.0e-5 * change) < 1e-11,
                "GPS Poisson ratio must use the paired radial endpoint temperature");
            require(relative(inferred_young, young - 8.0e7 * interpolated_change) > 1e-4
                        && relative(inferred_poisson, poisson + 2.0e-5 * interpolated_change) > 1e-4,
                "The elastic probe must independently distinguish both parameter temperatures");
            require(std::abs(trace - accumulated_trace + 3.0 * thermal_strain) < 1e-12,
                "GPS mean thermal expansion must persist across changing endpoint temperature histories");
            if (inelastic) {
                const double equivalent = std::sqrt(0.5
                                                        * ((stress.rr - stress.zz) * (stress.rr - stress.zz)
                                                            + (stress.zz - stress.hoop) * (stress.zz - stress.hoop)
                                                            + (stress.hoop - stress.rr) * (stress.hoop - stress.rr))
                                                    + 3.0 * stress.rz * stress.rz);
                const double yield = 2.0e8 + 2.0e5 * change, hardening = 1.0e9 + 2.0e7 * change;
                const double plastic = point.equivalent_plastic_strain;
                const double creep_rate = (point.equivalent_creep_strain - history[q].equivalent_creep_strain) / dt;
                require(plastic > history[q].equivalent_plastic_strain && creep_rate > 0.0,
                    "Both mechanisms must activate at both GPS points in every temperature increment");
                require(relative(equivalent, yield + hardening * plastic) < 1e-10,
                    "GPS yield stress and hardening must use paired radial temperatures");
                require(relative(equivalent, 2.0e8 + 2.0e5 * interpolated_change + hardening * plastic) > 1e-4
                            && relative(equivalent, yield + (1.0e9 + 2.0e7 * interpolated_change) * plastic) > 1e-4,
                    "The plastic probe must distinguish yield and hardening temperature separately");
                const double coefficient = 1.0e-5 + 1.0e-7 * change, exponent = 3.0 + 4.0e-3 * change;
                require(relative(creep_rate, coefficient * std::pow(equivalent / 1.0e8, exponent)) < 1e-10,
                    "GPS backward Euler creep coefficient and exponent must use paired radial temperatures");
                require(relative(creep_rate,
                            (1.0e-5 + 1.0e-7 * interpolated_change) * std::pow(equivalent / 1.0e8, exponent))
                                > 1e-4
                            && relative(creep_rate,
                                   coefficient * std::pow(equivalent / 1.0e8, 3.0 + 4.0e-3 * interpolated_change))
                                   > 1e-4,
                    "The creep probe must distinguish coefficient and exponent temperature separately");
            }
        }
        constexpr double step = 1e-3;
        for (std::size_t node = 0; node < 2; ++node) {
            auto plus = state, minus = state;
            plus[node] += step;
            minus[node] -= step;
            const auto rp =
                evaluate_cax2t_gps({material, geometry, plus, old, &history, dt, input.time, 0.0, formulation});
            const auto rm =
                evaluate_cax2t_gps({material, geometry, minus, old, &history, dt, input.time, 0.0, formulation});
            for (std::size_t row = 0; row < state.size(); ++row)
                maximum_error = std::max(maximum_error,
                    error(active.jacobian[6 * row + node], (rp.residual[row] - rm.residual[row]) / (2.0 * step)));
        }
        require(maximum_error < 2e-7, "GPS endpoint-temperature tangent columns must agree with centered differences");
        require(same_history(history, saved_history) && state == saved_state && old == saved_old,
            "Endpoint-temperature perturbations must preserve accepted material history and nodal fields");
        history = active.history;
        old = state;
        previous_radial_stretch = radial_stretch;
        previous_axial_stretch = axial_stretch;
    }
    std::cout << "cax2t_gps_" << (formulation == StrainFormulation::finite ? "finite" : "small")
              << (inelastic ? "_coupled" : "_elastic") << "_paired_temperature_jacobian_error=" << maximum_error
              << '\n';
}

void check_nonlinear_transaction(StrainFormulation formulation) {
    const auto material = make_material(true, true);
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.5, 2.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues old{610.0, 630.0, 0.0015, 0.002, 0.0, 0.0005};
    const Cax2tGpsLocalValues state{640.0, 680.0, 0.002, 0.006, 0.001, 0.015};
    const Cax2tGpsMaterialHistory initial_history{};
    const auto history =
        evaluate_cax2t_gps({material, geometry, old, initial, &initial_history, 0.1, 0.1, 0.0, formulation}).history;
    const auto saved_history = history;
    const auto saved_state = state, saved_old = old;
    const Cax2tGpsInput input{material, geometry, state, old, &history, 0.1, 0.2, 2.0e6, formulation, true};
    const auto active = evaluate_cax2t_gps(input, {true, true, true, true});
    const auto passive = evaluate_cax2t_gps(input);
    require(active.residual == passive.residual && same_history(active.history, passive.history),
        "Tangent and residual-only calls must return identical residuals and trial histories");
    const auto tangent_only = evaluate_cax2t_gps(input, {false, true, false, true});
    const auto history_only = evaluate_cax2t_gps(input, {false, false, true, false});
    const auto stress_only = evaluate_cax2t_gps(input, {false, false, false, true});
    require(tangent_only.residual == active.residual && tangent_only.jacobian == active.jacobian
                && same_history(tangent_only.history, {}) && same_history(stress_only.history, {}),
        "Requesting stress or tangent must not force trial history output");
    require(history_only.residual == Cax2tGpsLocalResidual{} && stress_only.residual == Cax2tGpsLocalResidual{}
                && same_history(history_only.history, active.history),
        "Output-only requests must not expose unrequested residuals");
    for (std::size_t q = 0; q < history.size(); ++q) {
        require(active.history[q].equivalent_plastic_strain > history[q].equivalent_plastic_strain
                    && active.history[q].equivalent_creep_strain > history[q].equivalent_creep_strain,
            "Each integration point must produce simultaneous plastic and creep increments");
        const auto& point = active.history[q];
        require(std::abs(point.plastic_strain[0] + point.plastic_strain[1] + point.plastic_strain[2]) < 1e-14
                    && std::abs(point.creep_strain[0] + point.creep_strain[1] + point.creep_strain[2]) < 1e-14,
            "Plastic and creep histories must remain trace free");
        const auto& stress = active.stress[q];
        require(stress_only.stress[q].rr == stress.rr && stress_only.stress[q].zz == stress.zz
                    && stress_only.stress[q].hoop == stress.hoop && stress_only.stress[q].rz == stress.rz,
            "Stress output must be independent of requesting material history");
    }
    require(error(active.residual[0] + active.residual[1], active.stored_heat_rate - active.generated_heat_rate)
                < 1e-14,
        "Summed thermal residual must equal stored minus generated heat");
    require(history_only.stored_heat_rate == active.stored_heat_rate
                && history_only.generated_heat_rate == active.generated_heat_rate,
        "Request flags must preserve thermal conservation diagnostics");
    const Cax2tGpsLocalValues direction{30.0, -20.0, 0.001, -0.0005, 0.0007, -0.0009};
    double maximum_error = 0.0;
    for (const double step : {1.0e-5, 3.0e-6}) {
        auto plus = state, minus = state;
        for (std::size_t column = 0; column < state.size(); ++column) {
            plus[column] += step * direction[column];
            minus[column] -= step * direction[column];
        }
        const auto rp = evaluate_cax2t_gps({material, geometry, plus, old, &history, 0.1, 0.2, 2e6, formulation, true});
        const auto rm =
            evaluate_cax2t_gps({material, geometry, minus, old, &history, 0.1, 0.2, 2e6, formulation, true});
        for (std::size_t row = 0; row < state.size(); ++row) {
            double tangent = 0.0;
            for (std::size_t column = 0; column < state.size(); ++column)
                tangent += active.jacobian[6 * row + column] * direction[column];
            maximum_error =
                std::max(maximum_error, error(tangent, (rp.residual[row] - rm.residual[row]) / (2.0 * step)));
        }
    }
    require(maximum_error < 2e-7, "Coupled inelastic directional tangent must match independent centered differences");
    auto invalid = state;
    invalid[0] = -1.0;
    bool caught = false;
    try {
        (void)evaluate_cax2t_gps({material, geometry, invalid, old, &history, 0.1, 0.2, 2e6, formulation, true});
    } catch (const std::domain_error&) {
        caught = true;
    }
    require(caught, "Invalid trial temperature must raise a domain error");
    const auto retry = evaluate_cax2t_gps(input, {true, true, true, true});
    require(retry.residual == active.residual && retry.jacobian == active.jacobian
                && same_history(retry.history, active.history) && same_history(history, saved_history)
                && state == saved_state && old == saved_old,
        "Discarding a failed trial must preserve inputs, committed history and exact repeatability");
    std::cout << "cax2t_gps_" << (formulation == StrainFormulation::finite ? "finite" : "small")
              << "_coupled_directional_jacobian_error=" << maximum_error << '\n';
}

void check_finite_mechanics_and_thermal_history() {
    const auto material = make_material();
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsMaterialHistory history{};
    const double radial_stretch = 1.1, axial_stretch = 1.04;
    const Cax2tGpsLocalValues stretched{620.0,
        620.0,
        radial_stretch - 1.0,
        2.0 * (radial_stretch - 1.0),
        0.0,
        2.0 * (axial_stretch - 1.0)};
    const auto result = evaluate_cax2t_gps(
        {material, geometry, stretched, initial, &history, 0.5, 0.5, 0.0, StrainFormulation::finite, true},
        {true, true, true, true});
    const double radial_increment = 2.0 * (radial_stretch - 1.0) / (radial_stretch + 1.0);
    const double axial_increment = 2.0 * (axial_stretch - 1.0) / (axial_stretch + 1.0);
    const double eigenstrain = 20.0 * expansion;
    const double lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    const double shear = young / (2.0 * (1.0 + poisson));
    const double pressure = lambda * (2.0 * radial_increment + axial_increment - 3.0 * eigenstrain);
    const double sigma_r = pressure + 2.0 * shear * (radial_increment - eigenstrain);
    const double sigma_z = pressure + 2.0 * shear * (axial_increment - eigenstrain);
    const double area = 3.0 * pi * radial_stretch * radial_stretch;
    const double height = 2.0 * axial_stretch;
    for (const auto& stress : result.stress)
        require(error(stress.rr, sigma_r) < 1e-13 && error(stress.hoop, sigma_r) < 1e-13
                    && error(stress.zz, sigma_z) < 1e-13 && stress.rz == 0.0,
            "Finite homogeneous stretches must reproduce the independent diagonal midpoint material update");
    require(error(result.residual[4], -area * sigma_z) < 1e-13 && error(result.residual[5], area * sigma_z) < 1e-13,
        "Finite axial reactions must use the current cross-section area");
    require(error(result.residual[2], -2.0 * pi * radial_stretch * height * sigma_r) < 1e-13
                && error(result.residual[3], 4.0 * pi * radial_stretch * height * sigma_r) < 1e-13,
        "Finite radial reactions must use the current radius and axial length");
    require(error(result.stored_heat_rate, area * height * 5.0e5 * 20.0 / 0.5) < 1e-13,
        "Finite lumped heat storage must use the complete current volume");

    const double free_increment = expansion * 1000.0;
    const double free_stretch = (2.0 + free_increment) / (2.0 - free_increment);
    const auto solid = make_cax2t_gps_geometry({0.0, 2.0}, 0.0, 2.0);
    const Cax2tGpsLocalValues free{1600.0, 1600.0, 0.0, 2.0 * (free_stretch - 1.0), 0.0, 2.0 * (free_stretch - 1.0)};
    const auto expanded =
        evaluate_cax2t_gps({material, solid, free, initial, &history, 0.1, 0.1, 0.0, StrainFormulation::finite},
            {true, true, true, true});
    for (const auto& stress : expanded.stress)
        require(std::max({std::abs(stress.rr), std::abs(stress.zz), std::abs(stress.hoop)}) < 1e-3,
            "Finite axis-touching free expansion must be stress free under its midpoint integration rule");

    const Cax2tGpsLocalValues first{620.0, 620.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues second{650.0, 650.0, 0.0, 0.0, 0.0, 0.0};
    const auto first_history =
        evaluate_cax2t_gps({material, geometry, first, initial, &history, 0.1, 0.1, 0.0, StrainFormulation::finite})
            .history;
    const auto second_result = evaluate_cax2t_gps(
        {material, geometry, second, first, &first_history, 0.1, 0.2, 0.0, StrainFormulation::finite},
        {true, true, true, true});
    const double expected = -young * expansion * 50.0 / (1.0 - 2.0 * poisson);
    for (const auto& stress : second_result.stress)
        require(error(stress.rr, expected) < 1e-13 && error(stress.zz, expected) < 1e-13
                    && error(stress.hoop, expected) < 1e-13,
            "A second finite increment must retain the difference between old and new thermal eigenstrain");

    for (std::size_t failure = 0; failure < 6; ++failure) {
        auto bad_current = initial, bad_old = initial;
        auto& bad = failure < 3 ? bad_current : bad_old;
        const std::size_t component = failure % 3;
        if (component == 0)
            bad[3] = -2.0;
        if (component == 1)
            bad[5] = -3.0;
        if (component == 2)
            bad[2] = bad[3] = -3.0;
        bool caught = false;
        try {
            (void)evaluate_cax2t_gps(
                {material, geometry, bad_current, bad_old, &history, 0.1, 0.1, 0.0, StrainFormulation::finite});
        } catch (const std::domain_error&) {
            caught = true;
        }
        require(caught, "Finite current and committed radial, axial and hoop stretches must each remain positive");
    }
}

void check_invalid_input() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::array<Cax2tGpsGeometry, 6> invalid_geometry = {{{{1.0, 1.0}, 0.0, 1.0},
        {{2.0, 1.0}, 0.0, 1.0},
        {{-1.0, 1.0}, 0.0, 1.0},
        {{0.0, 1.0}, 1.0, 1.0},
        {{0.0, 1.0}, 0.0, nan},
        {{0.0, nan}, 0.0, 1.0}}};
    for (const auto& geometry : invalid_geometry) {
        bool caught = false;
        try {
            (void)make_cax2t_gps_geometry(geometry.radii, geometry.z_lower, geometry.z_upper);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        require(caught, "Invalid reference radius or height must be rejected");
    }
    const auto material = make_material();
    const auto geometry = make_cax2t_gps_geometry({0.0, 1.0}, 0.0, 1.0);
    const Cax2tGpsLocalValues state{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsMaterialHistory history{};
    for (std::size_t failure = 0; failure < 6; ++failure) {
        Cax2tGpsInput input{material, geometry, state, state};
        if (failure == 0)
            input.strain_formulation = static_cast<StrainFormulation>(99);
        if (failure == 1)
            input.include_thermal_time_term = true;
        if (failure == 2)
            input.committed_history = &history;
        if (failure == 3)
            input.time_step = -0.1;
        if (failure == 4)
            input.volumetric_heat_source = nan;
        if (failure == 5)
            input.time = nan;
        bool caught = false;
        try {
            (void)evaluate_cax2t_gps(input);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        require(caught, "Unsupported formulation and malformed time/source settings must be rejected explicitly");
    }
    const auto annulus = make_cax2t_gps_geometry({1.0, 2.0}, 0.0, 2.0);
    const auto solid = make_cax2t_gps_geometry({0.0, 2.0}, 0.0, 2.0);
    for (const auto formulation : {StrainFormulation::small, StrainFormulation::finite})
        for (const bool with_history : {false, true})
            for (const bool corrupt_committed : {false, true}) {
                // A small-strain steady evaluation does not consume old values.
                if (corrupt_committed && !with_history && formulation == StrainFormulation::small)
                    continue;
                for (std::size_t failure = 0; failure < 5; ++failure) {
                    auto current = state, old = state;
                    auto& invalid = corrupt_committed ? old : current;
                    const auto& test_geometry = failure == 1 || failure == 2 ? solid : annulus;
                    if (failure == 0) {
                        invalid[2] = invalid[3] = -1.1;
                        require(1.0 + 0.5 * (1.0 - 1.0 / std::sqrt(3.0)) + invalid[2] > 0.0,
                            "The crossed-axis regression must still have two positive Gauss-point radii");
                    } else if (failure == 1 || failure == 2)
                        invalid[2] = failure == 1 ? 0.1 : -0.1;
                    else if (failure == 3)
                        invalid[3] = -1.0;
                    else
                        invalid[5] = -2.0;
                    bool caught = false;
                    try {
                        (void)evaluate_cax2t_gps({material,
                            test_geometry,
                            current,
                            old,
                            with_history ? &history : nullptr,
                            with_history ? 0.1 : 0.0,
                            with_history ? 0.1 : 0.0,
                            0.0,
                            formulation});
                    } catch (const std::domain_error&) {
                        caught = true;
                    }
                    require(caught,
                        "Small and finite states must reject crossed-axis endpoints, moving axis nodes and collapsed "
                        "spans");
                }
            }
}
} // namespace

int main() {
    try {
        std::cout << std::scientific << std::setprecision(12);
        check_analytic_mechanics();
        check_thermal_operators();
        check_mean_thermal_expansion(StrainFormulation::small);
        check_mean_thermal_expansion(StrainFormulation::finite);
        check_small_mean_hoop_mechanics();
        check_finite_mean_hoop_history_and_geometry();
        check_cax4t_projection(StrainFormulation::small);
        check_cax4t_projection(StrainFormulation::finite);
        for (const auto formulation : {StrainFormulation::small, StrainFormulation::finite})
            for (const bool inelastic : {false, true})
                check_paired_material_temperatures(formulation, inelastic);
        check_nonlinear_transaction(StrainFormulation::small);
        check_nonlinear_transaction(StrainFormulation::finite);
        check_finite_mechanics_and_thermal_history();
        check_invalid_input();
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
