#include "support/cax4_common_tests.hpp"

namespace {
using namespace fuelsim::test::cax4;

bool test_temperature_active_thermoelastic_properties() {
    const fuelsim::ThermoelasticProperties active_properties =
        fuelsim::test::thermoelastic(3824.0, 0.61, 2.0e11, 0.316, 1.0e-5, 600.0, -8.0e7, 2.0e-5, 3.0e-9);
    const fuelsim::IsotropicThermoelasticMaterial material(active_properties);
    constexpr double temperature = 725.0;
    const adlite::Scalar active_temperature = adlite::Scalar::independent(temperature, 0, 1);
    const fuelsim::AxisymmetricStress active = material.stress(1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, active_temperature);
    constexpr double step = 1.0e-3;
    const fuelsim::AxisymmetricStress plus = material.stress(1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, temperature + step);
    const fuelsim::AxisymmetricStress minus = material.stress(1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, temperature - step);
    const std::array<double, 4> analytic = {active.rr.derivative(0),
        active.zz.derivative(0),
        active.hoop.derivative(0),
        active.rz.derivative(0)};
    const std::array<double, 4> finite_difference = {(plus.rr.value() - minus.rr.value()) / (2.0 * step),
        (plus.zz.value() - minus.zz.value()) / (2.0 * step),
        (plus.hoop.value() - minus.hoop.value()) / (2.0 * step),
        (plus.rz.value() - minus.rz.value()) / (2.0 * step)};
    double maximum_error = 0.0;
    for (std::size_t component = 0; component < analytic.size(); ++component)
        maximum_error = std::max(maximum_error, scaled_error(analytic[component], finite_difference[component]));
    std::cout << "active_thermoelastic_temperature_tangent_error=" << maximum_error << '\n';
    return check(maximum_error < 1.0e-8,
        "temperature-dependent thermoelastic AD tangent matches "
        "centered differences");
}

int run_cax4t_contract_tests() {
    std::cout << std::scientific << std::setprecision(12);
    return test_cax_kinematics_and_jacobian(false) && test_temperature_active_thermoelastic_properties() ? 0 : 1;
}

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
        // Exercise nodal capacity with distinct spatial values and a nonzero
        // temperature derivative in the existing nonlinear transaction test.
        functions->thermal.function = [](const ThermoelasticFunctionInput& input, ThermalPropertyOutput& output) {
            output.conductivity = 120.0 / input.temperature + 3.0;
            output.density = 1000.0 + 10.0 * input.context.x + 20.0 * input.context.z + input.context.time;
            output.specific_heat = 500.0 + 0.1 * (input.temperature - 600.0);
        };
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

IsotropicThermoelasticMaterial corner_temperature_material(bool inelastic) {
    const auto registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "cax4t_corner_temperature_test";
    functions->thermal = registry.bind_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", 120.0},
            {"conductivity_constant", 3.0},
            {"density", 1000.0},
            {"specific_heat", 500.0}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", 1e9},
            {"poisson_ratio", 0.25},
            {"reference_temperature", 600.0},
            {"young_modulus_temperature_coefficient", -1e6},
            {"poisson_ratio_temperature_coefficient", 3e-4}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", 1e-5},
            {"reference_temperature", 600.0},
            {"thermal_expansion_temperature_coefficient", 2e-9}}));
    if (inelastic) {
        functions->plasticity = registry.bind_plasticity("linear_temperature_isotropic_hardening",
            {{"yield_stress", 2e6},
                {"hardening_modulus", 2e7},
                {"reference_temperature", 600.0},
                {"yield_stress_temperature_coefficient", 1e4},
                {"hardening_temperature_coefficient", 4e5}});
        functions->creep = registry.bind_creep("linear_temperature_norton",
            {{"coefficient", 1e-4},
                {"reference_stress", 2e6},
                {"stress_exponent", 2.0},
                {"reference_temperature", 600.0},
                {"coefficient_temperature_coefficient", 1e-6},
                {"reference_stress_temperature_coefficient", 0.0},
                {"stress_exponent_temperature_coefficient", 4e-3}});
    }
    return IsotropicThermoelasticMaterial({functions, 1e9});
}

void check_corner_temperature_materials(StrainFormulation form, bool inelastic) {
    const auto m = corner_temperature_material(inelastic);
    const auto g = make_cax4t_geometry({{{1.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}}});
    const std::array<std::array<double, 4>, 2> temperatures = {
        {{610.0, 640.0, 700.0, 750.0}, {750.0, 700.0, 650.0, 620.0}}};
    Cax4LocalValues old{};
    std::fill_n(old.begin(), 4, 600.0);
    Quad4MaterialHistory history{};
    double accumulated_trace = 0.0;
    double tangent_error = 0.0;
    double previous_radial_stretch = 1.0, previous_axial_stretch = 1.0;
    constexpr double dt = 0.25;
    for (std::size_t increment = 0; increment < temperatures.size(); ++increment) {
        Cax4LocalValues state{};
        const double radial_strain = 0.001 * static_cast<double>(increment + 1);
        const double axial_strain = 0.02 * static_cast<double>(increment + 1);
        double mean_temperature = 0.0;
        for (std::size_t node = 0; node < 4; ++node) {
            state[node] = temperatures[increment][node];
            state[4 + node] = radial_strain * g.coordinates[node].r;
            state[8 + node] = axial_strain * g.coordinates[node].z;
            mean_temperature += 0.25 * state[node];
        }
        const double radial_stretch = 1.0 + radial_strain, axial_stretch = 1.0 + axial_strain;
        if (form == StrainFormulation::finite)
            accumulated_trace +=
                4.0 * (radial_stretch - previous_radial_stretch) / (radial_stretch + previous_radial_stretch)
                + 2.0 * (axial_stretch - previous_axial_stretch) / (axial_stretch + previous_axial_stretch);
        else
            accumulated_trace = 2.0 * radial_strain + axial_strain;
        const auto saved_history = history;
        const auto saved_state = state, saved_old = old;
        // The existing transaction test covers nodal heat capacity.  Isolate
        // conduction here so its small off-diagonal temperature derivatives
        // remain resolved when perturbing the material-temperature columns.
        const Cax4Input
            input{m, g, state, old, &history, dt, dt * static_cast<double>(increment + 1), 0.0, form, false};
        const auto active = evaluate_cax4t(input, {true, true, true, false});
        const auto passive = evaluate_cax4t(input);
        require(active.residual == passive.residual && same_history(active.history, passive.history),
            "Corner-temperature material paths must return identical residual and trial history");
        for (std::size_t q = 0; q < 4; ++q) {
            const auto& point = active.history[q];
            const auto& strain = point.elastic_strain;
            const auto& stress = point.stress;
            const double trace = strain[0] + strain[1] + strain[2];
            const double shear = (stress.zz - stress.rr) / (2.0 * (strain[1] - strain[0]));
            const double bulk = (stress.rr + stress.zz + stress.hoop) / (3.0 * trace);
            const double inferred_young = 9.0 * bulk * shear / (3.0 * bulk + shear);
            const double inferred_poisson = (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear));
            const double temperature_change = state[q] - 600.0;
            double interpolated_temperature = 0.0;
            for (std::size_t node = 0; node < 4; ++node)
                interpolated_temperature += g.points[q].shape[node] * state[node];
            const double interpolated_change = interpolated_temperature - 600.0;
            require(relative_difference(inferred_young, 1e9 - 1e6 * temperature_change) < 1e-11,
                "Young modulus must use the paired corner temperature");
            require(relative_difference(inferred_poisson, 0.25 + 3e-4 * temperature_change) < 1e-11,
                "Poisson ratio must use the paired corner temperature");
            require(relative_difference(inferred_young, 1e9 - 1e6 * interpolated_change) > 1e-4
                        && relative_difference(inferred_poisson, 0.25 + 3e-4 * interpolated_change) > 1e-4,
                "Elastic probe must distinguish interpolated temperature for each parameter");
            const double mean_change = mean_temperature - 600.0;
            const double thermal_strain = (1e-5 + 2e-9 * mean_change) * mean_change;
            require(std::abs(trace - accumulated_trace + 3.0 * thermal_strain) < 1e-12,
                "Thermal expansion must retain mean temperature across changing corner histories");
            if (inelastic) {
                const double equivalent = std::sqrt(0.5
                                                        * ((stress.rr - stress.zz) * (stress.rr - stress.zz)
                                                            + (stress.zz - stress.hoop) * (stress.zz - stress.hoop)
                                                            + (stress.hoop - stress.rr) * (stress.hoop - stress.rr))
                                                    + 3.0 * stress.rz * stress.rz);
                const double yield = 2e6 + 1e4 * temperature_change;
                const double hardening = 2e7 + 4e5 * temperature_change;
                const double plastic = point.equivalent_plastic_strain;
                const double creep_rate = (point.equivalent_creep_strain - history[q].equivalent_creep_strain) / dt;
                require(plastic > history[q].equivalent_plastic_strain && creep_rate > 0.0,
                    "Every corner probe point must activate plasticity and creep in every increment");
                require(relative_difference(equivalent, yield + hardening * plastic) < 1e-10,
                    "Plastic yield stress and hardening must use paired corner temperature");
                require(relative_difference(equivalent, 2e6 + 1e4 * interpolated_change + hardening * plastic) > 1e-4
                            && relative_difference(equivalent, yield + (2e7 + 4e5 * interpolated_change) * plastic)
                                   > 1e-4,
                    "Plastic probe must independently distinguish yield and hardening temperature");
                const double coefficient = 1e-4 + 1e-6 * temperature_change;
                const double exponent = 2.0 + 4e-3 * temperature_change;
                require(relative_difference(creep_rate, coefficient * std::pow(equivalent / 2e6, exponent)) < 1e-10,
                    "Backward Euler creep coefficient and exponent must use paired corner temperature");
                require(relative_difference(creep_rate,
                            (1e-4 + 1e-6 * interpolated_change) * std::pow(equivalent / 2e6, exponent))
                                > 1e-4
                            && relative_difference(creep_rate,
                                   coefficient * std::pow(equivalent / 2e6, 2.0 + 4e-3 * interpolated_change))
                                   > 1e-4,
                    "Creep probe must independently distinguish coefficient and exponent temperature");
            }
        }
        // Each temperature column resolves corner material derivatives and the
        // independent arithmetic-mean thermal-expansion chain simultaneously.
        constexpr double step = 1e-3;
        for (std::size_t node = 0; node < 4; ++node) {
            auto plus = state, minus = state;
            plus[node] += step;
            minus[node] -= step;
            const auto rp = evaluate_cax4t({m, g, plus, old, &history, dt, input.time, 0.0, form, false});
            const auto rm = evaluate_cax4t({m, g, minus, old, &history, dt, input.time, 0.0, form, false});
            for (std::size_t row = 0; row < 12; ++row) {
                const double numerical = (rp.residual[row] - rm.residual[row]) / (2.0 * step);
                const double analytic = active.jacobian[12 * row + node];
                const double error =
                    std::abs(analytic - numerical) / std::max({1.0, std::abs(analytic), std::abs(numerical)});
                tangent_error = std::max(tangent_error, error);
                if (error >= 2e-7)
                    std::cerr << "corner temperature derivative row=" << row << " node=" << node
                              << " analytic=" << analytic << " numerical=" << numerical << " error=" << error << '\n';
            }
        }
        require(tangent_error < 2e-7, "Corner-temperature tangent columns must agree with centered differences");
        require(same_history(history, saved_history) && old == saved_old && state == saved_state,
            "Corner-temperature perturbations must not change accepted history or nodal inputs");
        history = active.history;
        old = state;
        previous_radial_stretch = radial_stretch;
        previous_axial_stretch = axial_stretch;
    }
    std::cout << (form == StrainFormulation::finite ? "finite" : "small") << (inelastic ? " coupled" : " elastic")
              << " corner material temperature tangent error=" << tangent_error << '\n';
}

void check_analytic_and_contract(StrainFormulation form) {
    const auto m = material(false);
    const auto g = fuelsim::elements::make_cax4t_geometry({{{1.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}}});
    Cax4LocalValues old{}, heated{};
    for (std::size_t n = 0; n < 4; ++n) {
        old[n] = 600.0;
        heated[n] = 620.0;
    }
    const Quad4MaterialHistory history{};
    const Cax4Input input{m, g, heated, old, &history, 0.5, 0.5, 2e7, form, true};
    const auto result = evaluate_cax4t(input, {true, true, true, false});
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
    const auto g = fuelsim::elements::make_cax4t_geometry({{{1.0, 0.0}, {2.1, 0.1}, {2.0, 1.2}, {0.9, 1.0}}});
    Cax4LocalValues initial{}, old{}, state{};
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
    const Cax4Input input{m, g, state, old, &history, 0.1, 0.2, 2e6, form, true};
    const auto active = evaluate_cax4t(input, {true, true, true, false});
    const auto passive = evaluate_cax4t(input);
    require(active.residual == passive.residual, "Residual-only and tangent calls must match exactly");
    require(same_history(active.history, passive.history), "Trial history must not depend on requesting tangent");
    const auto tangent_only = evaluate_cax4t(input, {false, true, false, true});
    const auto residual_only = evaluate_cax4t(input, {true, false, false, false});
    const auto history_only = evaluate_cax4t(input, {false, false, true, false});
    const auto stress_only = evaluate_cax4t(input, {false, false, false, true});
    require(tangent_only.residual == active.residual && tangent_only.jacobian == active.jacobian
                && residual_only.residual == active.residual,
        "Omitting trial history must preserve the exact residual and tangent");
    require(same_history(tangent_only.history, {}) && same_history(residual_only.history, {})
                && same_history(stress_only.history, {}),
        "Unrequested trial history must be empty");
    require(same_history(history_only.history, active.history) && history_only.residual == Cax4LocalResidual{}
                && stress_only.residual == Cax4LocalResidual{},
        "History and stress requests must not require residual output");
    for (std::size_t q = 0; q < 4; ++q) {
        const auto& expected = active.history[q].stress;
        for (const auto& actual : {tangent_only.stress[q], stress_only.stress[q]})
            require(actual.rr == expected.rr && actual.zz == expected.zz && actual.hoop == expected.hoop
                        && actual.rz == expected.rz,
                "Stress output must not require trial history output");
    }
    require(residual_only.stored_heat_rate == active.stored_heat_rate
                && residual_only.generated_heat_rate == active.generated_heat_rate
                && history_only.stored_heat_rate == active.stored_heat_rate
                && history_only.generated_heat_rate == active.generated_heat_rate,
        "Request flags must preserve thermal conservation diagnostics");
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
    auto without_capacity = input;
    without_capacity.include_thermal_time_term = false;
    const auto stationary = evaluate_cax4t(without_capacity, {true, false, false, false});
    for (std::size_t n = 0; n < 4; ++n) {
        double weight = 0.0;
        for (const auto& point : g.points)
            weight += point.shape[n] * point.weighted_measure;
        const double density = 1000.0 + 10.0 * g.coordinates[n].r + 20.0 * g.coordinates[n].z;
        const double rate = density * (500.0 + 0.1 * (state[n] - 600.0)) * (state[n] - old[n]) / input.time_step;
        require(scaled_error(active.residual[n] - stationary.residual[n], weight * rate) < 1e-12,
            "Capacity must use initial density at time zero, reference nodal mass and current temperature");
    }
    // Resolve the fixed-mass thermal derivative independently of its opposing
    // current-volume source contribution. All thermal and mechanical entries
    // remain active in this mixed direction and both difference step sizes.
    const Cax4LocalValues direction = {2.0, -3.0, 4.0, -1.0, 0.3, -0.5, 0.2, 0.4, -0.2, 0.35, -0.45, 0.25};
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
            const double row_error =
                std::abs(derivative - numerical) / std::max({1.0, std::abs(derivative), std::abs(numerical)});
            error = std::max(error, row_error);
            if (row_error >= 2e-7)
                std::cerr << "coupled capacity derivative step=" << step << " row=" << i << " analytic=" << derivative
                          << " numerical=" << numerical << " error=" << row_error << '\n';
        }
    }
    require(error < 2e-7, "Coupled inelastic tangent must agree with centered differences");
    auto invalid = state;
    if (form == StrainFormulation::finite)
        for (std::size_t n = 0; n < 4; ++n)
            invalid[4 + n] = -2.0 * g.coordinates[n].r;
    else
        invalid[0] = -1.0;
    bool caught = false;
    try {
        (void)evaluate_cax4t({m, g, invalid, old, &history, 0.1, 0.2, 2e6, form, true}, {true, true, true, false});
    } catch (const std::domain_error&) {
        caught = true;
    }
    require(caught, "Invalid trial must raise a domain error");
    const auto retry = evaluate_cax4t(input, {true, true, true, false});
    require(retry.residual == active.residual && retry.jacobian == active.jacobian
                && same_history(retry.history, active.history) && same_history(history, saved_history),
        "Discarding a failed trial must preserve exact repeatability");
    std::cout << (form == StrainFormulation::finite ? "finite" : "small")
              << " coupled directional derivative error=" << error << '\n';
}
} // namespace

int main() {
    if (run_cax4t_contract_tests() != 0)
        return 1;
    try {
        for (const auto form : {StrainFormulation::small, StrainFormulation::finite}) {
            check_analytic_and_contract(form);
            check_nonlinear_transaction(form);
            check_corner_temperature_materials(form, false);
            check_corner_temperature_materials(form, true);
        }
        std::cout << "CAX4T independent library contracts passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
