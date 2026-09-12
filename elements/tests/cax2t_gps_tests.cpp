#include "cax2t_gps.hpp"
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
            output.density = 1000.0 + 10.0 * input.context.x + 20.0 * input.context.z + input.context.time;
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
        functions->creep = registry.bind_creep("norton",
            {{"coefficient", 1.0e-5}, {"reference_stress", 1.0e8}, {"stress_exponent", 3.0}});
        functions->plasticity = registry.bind_plasticity("linear_isotropic_hardening",
            {{"yield_stress", 2.0e8}, {"hardening_modulus", 1.0e9}});
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
    const std::array<std::array<double, 2>, 2> mass_weights = {{{5.0, 3.0}, {3.0, 7.0}}};
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t column = 0; column < 2; ++column) {
            const double capacity = 2.0 * pi * 2.0 / 12.0 * mass_weights[row][column] * 5.0e5 / 0.5;
            const double conduction = (row == column ? 1.0 : -1.0) * pi * 2.0 * 5.0 * 3.0;
            require(error(response.jacobian[6 * row + column] - stationary.jacobian[6 * row + column], capacity)
                        < 1e-14,
                "Consistent capacity must match the analytical cylindrical two-node matrix including its off-diagonal");
            require(error(stationary.jacobian[6 * row + column], conduction) < 1e-14,
                "Radial conduction must match the analytical cylindrical stiffness");
        }
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t column = 2; column < 6; ++column)
            require(response.jacobian[6 * row + column] == 0.0,
                "Small-strain reference thermal operators must be independent of displacement");

    const Cax2tGpsLocalValues nonuniform{610.0, 650.0, 0.0, 0.0, 0.0, 0.0};
    const auto varying = evaluate_cax2t_gps({material, geometry, nonuniform, initial}, {true, true, false, true});
    const double gauss = 1.0 / std::sqrt(3.0);
    const std::array<double, 2> temperatures = {630.0 - 20.0 * gauss, 630.0 + 20.0 * gauss};
    for (std::size_t q = 0; q < temperatures.size(); ++q) {
        const double stress = -young * expansion * (temperatures[q] - 600.0) / (1.0 - 2.0 * poisson);
        require(error(varying.stress[q].zz, stress) < 1e-13,
            "Thermal expansion must use each radial integration-point temperature");
    }
}

void check_nonlinear_transaction(StrainFormulation formulation) {
    const auto material = make_material(true, true);
    const auto geometry = make_cax2t_gps_geometry({1.0, 2.0}, 0.5, 2.0);
    const Cax2tGpsLocalValues initial{600.0, 600.0, 0.0, 0.0, 0.0, 0.0};
    const Cax2tGpsLocalValues old{610.0, 630.0, 0.001, 0.002, 0.0, 0.0005};
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
        "Finite consistent heat storage must use the complete current volume");

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
