#include "cax4t.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {
using namespace fuelsim;
using namespace fuelsim::elements;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

IsotropicThermoelasticMaterial material(bool inelastic) {
    const auto registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "cax4t_library_test";
    functions->thermal = registry.bind_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", 120.0},
            {"conductivity_constant", 3.0},
            {"density", 1000.0},
            {"specific_heat", 500.0}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", 2e11},
            {"poisson_ratio", 0.3},
            {"reference_temperature", 600.0},
            {"young_modulus_temperature_coefficient", 0.0},
            {"poisson_ratio_temperature_coefficient", 0.0}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", 1e-5},
            {"reference_temperature", 600.0},
            {"thermal_expansion_temperature_coefficient", 0.0}}));
    if (inelastic) {
        functions->creep = registry.bind_creep("linear_temperature_norton",
            {{"coefficient", 1e-5},
                {"reference_stress", 1e8},
                {"stress_exponent", 3.0},
                {"reference_temperature", 600.0},
                {"coefficient_temperature_coefficient", 0.0},
                {"reference_stress_temperature_coefficient", 0.0},
                {"stress_exponent_temperature_coefficient", 0.0}});
        functions->plasticity = registry.bind_plasticity("linear_temperature_isotropic_hardening",
            {{"yield_stress", 2e8},
                {"hardening_modulus", 1e9},
                {"reference_temperature", 600.0},
                {"yield_stress_temperature_coefficient", 0.0},
                {"hardening_temperature_coefficient", 0.0}});
    }
    return IsotropicThermoelasticMaterial({functions, 2e11});
}

bool same_history(const Quad4MaterialHistory& a, const Quad4MaterialHistory& b) {
    for (std::size_t q = 0; q < 4; ++q) {
        const auto &x = a[q], &y = b[q];
        if (x.elastic_strain != y.elastic_strain || x.plastic_strain != y.plastic_strain
            || x.creep_strain != y.creep_strain || x.equivalent_plastic_strain != y.equivalent_plastic_strain
            || x.equivalent_creep_strain != y.equivalent_creep_strain || x.stress.rr != y.stress.rr
            || x.stress.zz != y.stress.zz || x.stress.hoop != y.stress.hoop || x.stress.rz != y.stress.rz)
            return false;
    }
    return true;
}

void check_analytic_and_contract(StrainFormulation form) {
    const auto m = material(false);
    const auto g = make_quad4_rz_geometry({{{1.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}}});
    LocalValues old{}, heated{};
    for (std::size_t n = 0; n < 4; ++n) {
        old[n] = 600.0;
        heated[n] = 620.0;
    }
    const Quad4MaterialHistory history{};
    const Cax4tInput input{m, g, heated, old, &history, 0.5, 0.5, 2e7, form, true};
    const auto result = evaluate_cax4t(input, true);
    const double volume = 3.0 * std::acos(-1.0);
    require(std::abs(result.stored_heat_rate / (2e7 * volume) - 1.0) < 1e-14, "Analytical stored heat rate");
    require(std::abs(result.generated_heat_rate / (2e7 * volume) - 1.0) < 1e-14, "Analytical generated heat rate");
    for (std::size_t n = 0; n < 4; ++n)
        require(std::abs(result.residual[n]) < 1e-7, "Uniform heating balances the body source");
    const double expected = -2e11 * 1e-5 * 20.0 / (1.0 - 2.0 * 0.3);
    for (const auto& h : result.history) {
        for (const double stress : {h.stress.rr, h.stress.zz, h.stress.hoop})
            require(std::abs(stress / expected - 1.0) < 1e-13, "Restrained heating analytical stress");
    }
    const auto steady = evaluate_cax4t({m, g, heated, {}, nullptr, 0.0, 0.5, 0.0, form});
    require(steady.stored_heat_rate == 0.0, "Steady evaluation has no storage");
    require(steady.history[0].stress.rr == result.history[0].stress.rr, "Steady and transient elastic stress");
    for (const auto dt : {0.0, -0.1, std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid = input;
        invalid.time_step = dt;
        bool caught = false;
        try {
            (void)evaluate_cax4t(invalid);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        require(caught, "Invalid time increment must be rejected");
    }
}

void check_nonlinear_transaction(StrainFormulation form) {
    const auto m = material(true);
    const auto g = make_quad4_rz_geometry({{{1.0, 0.0}, {2.1, 0.1}, {2.0, 1.2}, {0.9, 1.0}}});
    LocalValues initial{}, old{}, state{};
    for (std::size_t n = 0; n < 4; ++n) {
        initial[n] = 600.0;
        old[n] = 610.0 + 10.0 * static_cast<double>(n);
        state[n] = 630.0 + 15.0 * static_cast<double>(n);
        old[4 + n] = 0.001 * g.coordinates[n].r;
        old[8 + n] = -0.0005 * g.coordinates[n].z;
        state[4 + n] = 0.004 * g.coordinates[n].r;
        state[8 + n] = -0.002 * g.coordinates[n].z;
    }
    state[5] += 0.001;
    state[10] -= 0.001;
    const Quad4MaterialHistory initial_history{};
    const auto history = evaluate_cax4t({m, g, old, initial, &initial_history, 0.1, 0.1, 0.0, form}).history;
    const auto saved_history = history;
    const auto saved_state = state, saved_old = old;
    const Cax4tInput input{m, g, state, old, &history, 0.1, 0.2, 2e6, form, true};
    const auto active = evaluate_cax4t(input, true);
    const auto passive = evaluate_cax4t(input);
    require(active.residual == passive.residual, "Residual-only and tangent calls must match exactly");
    require(same_history(active.history, passive.history), "Trial history must not depend on requesting tangent");
    require(same_history(history, saved_history) && state == saved_state && old == saved_old, "Inputs are immutable");
    bool plastic = false, creep = false;
    for (std::size_t q = 0; q < 4; ++q) {
        plastic = plastic || active.history[q].equivalent_plastic_strain > history[q].equivalent_plastic_strain;
        creep = creep || active.history[q].equivalent_creep_strain > history[q].equivalent_creep_strain;
    }
    require(plastic && creep, "The test must exercise both plasticity and creep increments");
    double heat_sum = 0.0;
    for (std::size_t n = 0; n < 4; ++n)
        heat_sum += active.residual[n];
    require(std::abs(heat_sum - active.stored_heat_rate + active.generated_heat_rate)
                < 1e-12 * std::abs(active.stored_heat_rate),
        "Thermal diagnostics must match summed residual");
    const LocalValues direction = {0.2, -0.3, 0.4, -0.1, 0.3, -0.5, 0.2, 0.4, -0.2, 0.35, -0.45, 0.25};
    // Two resolved step sizes avoid cancellation in the large thermal residual.
    double error = 0.0;
    for (const double step : {1e-5, 3e-6}) {
        auto plus = state, minus = state;
        for (std::size_t j = 0; j < 12; ++j) {
            plus[j] += step * direction[j];
            minus[j] -= step * direction[j];
        }
        const auto rp = evaluate_cax4t({m, g, plus, old, &history, 0.1, 0.2, 2e6, form, true});
        const auto rm = evaluate_cax4t({m, g, minus, old, &history, 0.1, 0.2, 2e6, form, true});
        for (std::size_t i = 0; i < 12; ++i) {
            double derivative = 0.0;
            for (std::size_t j = 0; j < 12; ++j)
                derivative += active.jacobian[12 * i + j] * direction[j];
            const double numerical = (rp.residual[i] - rm.residual[i]) / (2.0 * step);
            error = std::max(error,
                std::abs(derivative - numerical) / std::max({1.0, std::abs(derivative), std::abs(numerical)}));
        }
        require(error < 2e-7, "Coupled inelastic tangent must agree with centered differences");
    }
    auto invalid = state;
    if (form == StrainFormulation::finite)
        for (std::size_t n = 0; n < 4; ++n)
            invalid[4 + n] = -2.0 * g.coordinates[n].r;
    else
        invalid[0] = -1.0;
    bool caught = false;
    try {
        (void)evaluate_cax4t({m, g, invalid, old, &history, 0.1, 0.2, 2e6, form, true}, true);
    } catch (const std::domain_error&) {
        caught = true;
    }
    require(caught, "Invalid trial must raise a domain error");
    const auto retry = evaluate_cax4t(input, true);
    require(retry.residual == active.residual && retry.jacobian == active.jacobian
                && same_history(retry.history, active.history) && same_history(history, saved_history),
        "Discarding a failed trial must preserve exact repeatability");
    std::cout << (form == StrainFormulation::finite ? "finite" : "small")
              << " coupled directional derivative error=" << error << '\n';
}
} // namespace

int main() {
    try {
        for (const auto form : {StrainFormulation::small, StrainFormulation::finite}) {
            check_analytic_and_contract(form);
            check_nonlinear_transaction(form);
        }
        std::cout << "CAX4T independent library contracts passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
