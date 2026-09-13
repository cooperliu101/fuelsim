#pragma once
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "contact_types.hpp"
#include "line2_rz.hpp"
#include "quad4_face.hpp"
#include "support/test_support.hpp"

#include "material_types.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuelsim::test::cax4 {
constexpr double pi = 3.141592653589793238462643383279502884;

inline bool check(bool condition, const std::string& message) {
    if (!condition)
        std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

inline double scaled_error(double actual, double expected) {
    return std::abs(actual - expected) / (1.0 + std::max(std::abs(actual), std::abs(expected)));
}

inline double relative_difference(double actual, double expected) {
    return std::abs(actual - expected) / std::abs(expected);
}

inline fuelsim::ThermoelasticProperties properties() {
    return fuelsim::test::thermoelastic(3824.0, 0.61, 2.0e11, 0.316, 1.0e-5, 600.0);
}

inline bool test_cax_kinematics_and_jacobian(bool reduced) {
    const std::string name = reduced ? "CAX4RT" : "CAX4T";
    const std::size_t points = reduced ? 1 : 4;
    const fuelsim::Quad4Coordinates coordinates = {{{1.0, 0.0}, {2.1, 0.1}, {2.0, 1.2}, {0.9, 1.0}}};
    const auto geometry = (reduced ? fuelsim::elements::make_cax4rt_geometry(coordinates)
                                   : fuelsim::elements::make_cax4t_geometry(coordinates));
    fuelsim::AxisymmetricTestData data{fuelsim::IsotropicThermoelasticMaterial(properties())};
    data.element_formulation = reduced ? fuelsim::RzElementFormulation::cax4rt : fuelsim::RzElementFormulation::cax4t;
    const fuelsim::Cax4LocalValues direction = {0.2, -0.3, 0.4, -0.1, 0.3, -0.5, 0.2, 0.4, -0.2, 0.35, -0.45, 0.25};
    bool passed = true;
    const auto evaluate = [&](const fuelsim::Cax4LocalValues& state,
                              const fuelsim::Cax4LocalValues& committed,
                              const fuelsim::Quad4MaterialHistory& history,
                              fuelsim::elements::ElementRequest request,
                              bool thermal_time = true) {
        const fuelsim::elements::Cax4Input input{data.material,
            geometry,
            state,
            committed,
            &history,
            0.1,
            data.time,
            data.volumetric_heat_source,
            data.strain_formulation,
            thermal_time,
            data.initial_temperature};
        return reduced ? fuelsim::elements::evaluate_cax4rt(input, request)
                       : fuelsim::elements::evaluate_cax4t(input, request);
    };
    {
        fuelsim::Cax4LocalValues heated{}, committed{};
        for (std::size_t n = 0; n < 4; ++n) {
            heated[n] = 610.0;
            committed[n] = 600.0;
        }
        fuelsim::elements::Cax4Input input{data.material, geometry, heated, committed};
        input.include_thermal_time_term = true;
        for (double step : {0.1, 0.0, std::numeric_limits<double>::quiet_NaN()}) {
            input.time_step = step;
            for (bool jacobian : {false, true}) {
                bool rejected = false;
                try {
                    if (reduced)
                        (void)fuelsim::elements::evaluate_cax4rt(input, {true, jacobian, false, false});
                    else
                        (void)fuelsim::elements::evaluate_cax4t(input, {true, jacobian, false, false});
                } catch (const std::invalid_argument&) {
                    rejected = true;
                }
                passed = check(rejected, name + " heat capacity rejects missing committed history") && passed;
            }
        }
        input.include_thermal_time_term = false;
        input.time_step = 0.0;
        const auto steady =
            reduced ? fuelsim::elements::evaluate_cax4rt(input) : fuelsim::elements::evaluate_cax4t(input);
        passed = check(steady.stored_heat_rate == 0.0, name + " steady evaluation permits missing history") && passed;
    }
    if (!reduced) {
        // Fully restrained heating has the same hydrostatic stress at every
        // point, even on a distorted element with nonuniform corner temperatures.
        for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
            data.strain_formulation = formulation;
            fuelsim::Cax4LocalValues initial{}, old{}, heated{};
            for (std::size_t n = 0; n < 4; ++n)
                initial[n] = 600.0;
            old[0] = 610.0;
            old[1] = 630.0;
            old[2] = 650.0;
            old[3] = 710.0;
            heated[0] = 650.0;
            heated[1] = 690.0;
            heated[2] = 710.0;
            heated[3] = 750.0;
            const auto history = evaluate(old, initial, {}, {false, false, true, false}).history;
            const auto next = evaluate(heated, old, history, {false, false, true, false}).history;
            for (const auto& point : next) {
                const std::array<double, 3> stress = {point.stress.rr, point.stress.zz, point.stress.hoop};
                for (std::size_t component = 0; component < 3; ++component) {
                    passed = check(std::abs(point.elastic_strain[component] + 1e-3) < 1e-14,
                                 "CAX4T expansion uses arithmetic corner temperature across thermal history")
                             && passed;
                    const double expected = -2e11 * 1e-3 / (1.0 - 2.0 * 0.316);
                    passed = check(relative_difference(stress[component], expected) < 1e-12,
                                 "CAX4T restrained nonuniform heating agrees with analytical hydrostatic stress")
                             && passed;
                }
            }
        }
    }
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        data.strain_formulation = formulation;
        data.volumetric_heat_source = formulation == fuelsim::StrainFormulation::finite ? 2e6 : 0.0;
        fuelsim::Cax4LocalValues initial{}, old{}, state{};
        for (std::size_t n = 0; n < 4; ++n) {
            initial[n] = old[n] = state[n] = 600.0;
            old[4 + n] = 0.03 * coordinates[n].r;
            old[8 + n] = -0.02 * coordinates[n].z;
            state[4 + n] = 0.08 * coordinates[n].r;
            state[8 + n] = -0.04 * coordinates[n].z;
        }
        const auto old_history = evaluate(old, initial, {}, {false, false, true, false}).history;
        const auto affine = evaluate(state, old, old_history, {false, false, true, false}).history;
        const bool finite = formulation == fuelsim::StrainFormulation::finite;
        const double radial = finite ? 2.0 * 0.03 / 2.03 + 2.0 * 0.05 / 2.11 : 0.08;
        const double axial = finite ? -2.0 * 0.02 / 1.98 - 2.0 * 0.02 / 1.94 : -0.04;
        double affine_error = 0.0;
        for (std::size_t q = 0; q < points; ++q) {
            const auto& h = affine[q];
            affine_error = std::max({affine_error,
                std::abs(h.elastic_strain[0] - radial),
                std::abs(h.elastic_strain[1] - axial),
                std::abs(h.elastic_strain[2] - radial),
                std::abs(h.elastic_strain[3])});
        }
        passed = check(affine_error < 2e-14, name + " homogeneous stretch follows the analytical incremental strain")
                 && passed;
        if (reduced)
            passed = check(fuelsim::elements::cax4rt_hourglass_energy({data.material,
                               geometry,
                               state,
                               state,
                               nullptr,
                               0.0,
                               data.time,
                               data.volumetric_heat_source,
                               data.strain_formulation,
                               false,
                               data.initial_temperature})
                               < 1e-20,
                         "CAX4RT affine deformation has no artificial hourglass energy")
                     && passed;
        state[5] += 0.025;
        state[6] += 0.015;
        state[10] -= 0.020;
        state[11] += 0.010;
        if (reduced) {
            auto heated = state;
            auto unheated = state;
            for (std::size_t n = 0; n < 4; ++n) {
                unheated[n] = 600.0;
                heated[n] = 601.0;
            }
            const double saved_source = data.volumetric_heat_source;
            data.volumetric_heat_source = data.material.heat_capacity(601.0, {}).value() / 0.1;
            const auto residual = evaluate(heated, unheated, {}, {true, false, false, false}).residual;
            for (std::size_t n = 0; n < 4; ++n)
                passed = check(std::abs(residual[n]) < 1e-10 * data.volumetric_heat_source,
                             "CAX4RT uniform heating balances uniform source on a distorted element")
                         && passed;
            data.volumetric_heat_source = saved_source;
        }
        state[0] += 3.0;
        state[2] -= 2.0;
        fuelsim::Cax4LocalJacobian jacobian{};
        const auto active_result = evaluate(state, old, old_history, {true, true, false, false});
        const auto& active = active_result.residual;
        jacobian = active_result.jacobian;
        const auto passive = evaluate(state, old, old_history, {true, false, false, false}).residual;
        if (!reduced) {
            const auto no_capacity = evaluate(state, old, old_history, {true, false, false, false}, false).residual;
            auto thermal_coordinates = coordinates;
            if (finite)
                for (std::size_t n = 0; n < 4; ++n) {
                    thermal_coordinates[n].r += state[4 + n];
                    thermal_coordinates[n].z += state[8 + n];
                }
            const auto thermal_geometry = (reduced ? fuelsim::elements::make_cax4rt_geometry(thermal_coordinates)
                                                   : fuelsim::elements::make_cax4t_geometry(thermal_coordinates));
            for (std::size_t n = 0; n < 4; ++n) {
                double weight = 0.0;
                for (const auto& p : thermal_geometry.points)
                    weight += p.weighted_measure * p.shape[n];
                const auto& x = coordinates[n];
                const double expected = weight * data.material.heat_capacity(state[n], {data.time, x.r, 0, x.z}).value()
                                        * (state[n] - old[n]) / 0.1;
                passed = check(scaled_error(active[n] - no_capacity[n], expected) < 1e-11,
                             "CAX4T nodal capacity uses row-sum weight and each node's temperature rate")
                         && passed;
            }
        }
        passed =
            check(active == passive, name + " residual and Jacobian evaluations return identical residuals") && passed;
        constexpr double step = 1e-5;
        auto plus = state, minus = state;
        for (std::size_t j = 0; j < 12; ++j) {
            plus[j] += step * direction[j];
            minus[j] -= step * direction[j];
        }
        const auto rp = evaluate(plus, old, old_history, {true, false, false, false}).residual;
        const auto rm = evaluate(minus, old, old_history, {true, false, false, false}).residual;
        double derivative_error = 0.0;
        for (std::size_t i = 0; i < 12; ++i) {
            double ad = 0.0;
            for (std::size_t j = 0; j < 12; ++j)
                ad += jacobian[12 * i + j] * direction[j];
            derivative_error = std::max(derivative_error, scaled_error(ad, (rp[i] - rm[i]) / (2.0 * step)));
        }
        passed = check(derivative_error < 2e-7, name + " coupled volume and hoop Jacobian matches centered differences")
                 && passed;
        std::cout << name + "_" << (finite ? "finite" : "small") << "_directional_jacobian_error=" << derivative_error
                  << '\n';
        // Axial rigid translation is an exact axisymmetric rigid motion, including with committed stress.
        auto translated = state;
        for (std::size_t n = 0; n < 4; ++n)
            translated[8 + n] += 0.5;
        const auto shifted = evaluate(translated, old, old_history, {true, false, false, false}).residual;
        double translation_error = 0.0;
        for (std::size_t i = 0; i < 12; ++i)
            translation_error = std::max(translation_error, scaled_error(shifted[i], passive[i]));
        passed = check(translation_error < 1e-12, name + " internal forces are invariant under axial rigid translation")
                 && passed;
        if (reduced) {
            const double energy = fuelsim::elements::cax4rt_hourglass_energy({data.material,
                geometry,
                state,
                state,
                nullptr,
                0.0,
                data.time,
                data.volumetric_heat_source,
                data.strain_formulation,
                false,
                data.initial_temperature});
            const double shifted_energy = fuelsim::elements::cax4rt_hourglass_energy({data.material,
                geometry,
                translated,
                translated,
                nullptr,
                0.0,
                data.time,
                data.volumetric_heat_source,
                data.strain_formulation,
                false,
                data.initial_temperature});
            passed = check(energy > 0.0 && std::abs(energy - shifted_energy) < 1e-12 * energy,
                         "CAX4RT hourglass energy is positive and invariant under axial translation")
                     && passed;
        }
        if (finite) {
            auto invalid_midpoint = initial;
            for (std::size_t n = 0; n < 4; ++n) {
                invalid_midpoint[4 + n] = 5.0 - 3.0 * coordinates[n].r;
                invalid_midpoint[8 + n] = -1.5 * coordinates[n].z;
            }
            // det(F_new)=1 and all current radii are positive, but det(F_mid)<0.
            bool rejected = false;
            try {
                (void)evaluate(invalid_midpoint, initial, {}, {true, false, false, false}).residual;
            } catch (const std::domain_error& error) {
                rejected = std::string(error.what()).find("midpoint") != std::string::npos;
            }
            passed = check(rejected, name + " rejects an invalid midpoint even when its current volume is positive")
                     && passed;
        }
    }
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        data.strain_formulation = formulation;
        data.volumetric_heat_source = reduced && formulation == fuelsim::StrainFormulation::finite ? 2e6 : 0.0;
        for (int mechanism = 0; mechanism < 3; ++mechanism) {
            auto material = properties();
            if (mechanism != 0)
                material = fuelsim::test::with_norton(material, 1e-4, 1e8, 3.0);
            if (mechanism != 1)
                material = fuelsim::test::with_plasticity(material, 1e8, 2e10);
            data.material = fuelsim::IsotropicThermoelasticMaterial(material);
            fuelsim::Cax4LocalValues initial{}, old{}, state{};
            for (std::size_t n = 0; n < 4; ++n) {
                initial[n] = old[n] = state[n] = 600.0;
                old[4 + n] = 0.001 * coordinates[n].r;
                old[8 + n] = -0.001 * coordinates[n].z;
                state[4 + n] = 0.004 * coordinates[n].r;
                state[8 + n] = -0.005 * coordinates[n].z;
            }
            state[5] += 0.0003;
            state[10] -= 0.0002;
            const auto history = evaluate(old, initial, {}, {false, false, true, false}).history;
            fuelsim::Cax4LocalJacobian jacobian{};
            const auto residual_result = evaluate(state, old, history, {true, true, false, false});
            const auto& residual = residual_result.residual;
            jacobian = residual_result.jacobian;
            passed = check(residual == evaluate(state, old, history, {true, false, false, false}).residual,
                         name + " inelastic material residual agrees exactly between passive and Jacobian paths")
                     && passed;
            {
                const auto full = evaluate(state, old, history, {true, true, true, true});
                const auto tangent = evaluate(state, old, history, {true, true, false, false});
                const auto residual_request = evaluate(state, old, history, {true, false, false, false});
                const auto history_request = evaluate(state, old, history, {false, false, true, false});
                const auto stress_request = evaluate(state, old, history, {false, false, false, true});
                if (full.residual != tangent.residual || full.jacobian != tangent.jacobian
                    || full.residual != residual_request.residual
                    || history_request.residual != decltype(full.residual){}
                    || stress_request.residual != decltype(full.residual){})
                    throw std::runtime_error("Element request changed residual or tangent");
                for (std::size_t q = 0; q < full.history.size(); ++q) {
                    const auto &a = full.history[q], &b = history_request.history[q];
                    if (!fuelsim::test::same_material_state(a, b)
                        || !fuelsim::test::same_material_state(tangent.history[q], {})
                        || !fuelsim::test::same_material_state(residual_request.history[q], {})
                        || !fuelsim::test::same_material_state(stress_request.history[q], {})
                        || stress_request.stress[q].rr != a.stress.rr || stress_request.stress[q].zz != a.stress.zz
                        || stress_request.stress[q].hoop != a.stress.hoop || stress_request.stress[q].rz != a.stress.rz)
                        throw std::runtime_error("Element request changed history or stress");
                }
                if (full.generated_heat_rate != residual_request.generated_heat_rate
                    || full.stored_heat_rate != residual_request.stored_heat_rate
                    || full.generated_heat_rate != history_request.generated_heat_rate
                    || full.stored_heat_rate != history_request.stored_heat_rate)
                    throw std::runtime_error("Element request changed heat diagnostics");
            }
            auto plus = state, minus = state;
            constexpr double step = 1e-5;
            for (std::size_t j = 0; j < 12; ++j) {
                plus[j] += step * direction[j];
                minus[j] -= step * direction[j];
            }
            const auto rp = evaluate(plus, old, history, {true, false, false, false}).residual;
            const auto rm = evaluate(minus, old, history, {true, false, false, false}).residual;
            double error = 0.0;
            for (std::size_t i = 0; i < 12; ++i) {
                double ad = 0.0;
                for (std::size_t j = 0; j < 12; ++j)
                    ad += jacobian[12 * i + j] * direction[j];
                error = std::max(error, scaled_error(ad, (rp[i] - rm[i]) / (2.0 * step)));
            }
            const auto next = evaluate(state, old, history, {false, false, true, false}).history;
            for (std::size_t q = 0; q < points; ++q) {
                const auto& point = next[q];
                if (mechanism != 0)
                    passed = check(point.equivalent_creep_strain > 1e-6, name + " creep branch is active") && passed;
                if (mechanism != 1)
                    passed =
                        check(point.equivalent_plastic_strain > 1e-6, name + " plastic branch is active") && passed;
            }
            passed = check(error < 2e-6, name + " active inelastic Jacobian matches centered differences") && passed;
            std::cout << name + "_inelastic_" << static_cast<int>(formulation) << '_' << mechanism
                      << "_jacobian_error=" << error << '\n';
        }
    }
    return passed;
}

using HeatPointGeometries =
    std::array<fuelsim::Line2RzHeatPointGeometry, fuelsim::line2_interface_quadrature_point_count>;

} // namespace fuelsim::test::cax4
