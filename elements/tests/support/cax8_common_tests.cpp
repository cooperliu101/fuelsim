#include "cax8_types.hpp"
#include "contact_types.hpp"
#include "line3_rz.hpp"
#include "quad4_face.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

int run_cax8_tests(fuelsim::RzElementFormulation selected) {
    try {
        const auto evaluate = [selected](const fuelsim::elements::Cax8Input& input,
                                  fuelsim::elements::ElementRequest request) {
            return selected == fuelsim::RzElementFormulation::cax8rt
                       ? fuelsim::elements::evaluate_cax8rt(input, request)
                       : fuelsim::elements::evaluate_cax8t(input, request);
        };
        const fuelsim::Quad8RzCoordinates coordinates = {
            {{1, 0}, {2.1, .1}, {1.9, 1.2}, {.9, 1}, {1.55, .03}, {2.02, .65}, {1.4, 1.12}, {.93, .5}}};

        auto base = fuelsim::test::thermoelastic(0, 10, 1e6, .25, 1e-5, 600, 0, 0, 0, 1000, 100);
        fuelsim::Quad8RzValues state = {600,
                                   620,
                                   590,
                                   610,
                                   .01,
                                   .03,
                                   .06,
                                   -.02,
                                   .025,
                                   .038,
                                   .022,
                                   -.002,
                                   .02,
                                   -.01,
                                   .07,
                                   .05,
                                   .003,
                                   .034,
                                   .065,
                                   .032},
                               old{};
        for (std::size_t n = 0; n < 4; ++n)
            old[n] = 600;
        double maximum = 0;
        // An affine annulus has exact volume and uniform strain under either rule.
        const fuelsim::Quad8RzCoordinates annulus = {
            {{1, 0}, {2, 0}, {2, 1}, {1, 1}, {1.5, 0}, {2, .5}, {1.5, 1}, {1, .5}}};
        for (auto element : {selected}) {
            const auto geometry = fuelsim::test::make_cax8_geometry(annulus, element);
            fuelsim::AxisymmetricTestData data{fuelsim::IsotropicThermoelasticMaterial(base)};
            data.volumetric_heat_source = 1e4;
            fuelsim::Quad8RzValues affine{};
            for (std::size_t n = 0; n < 4; ++n)
                affine[n] = 600;
            for (std::size_t n = 0; n < 8; ++n) {
                affine[4 + n] = .001 * annulus[n].r;
                affine[12 + n] = -.002 * annulus[n].z;
            }
            const auto result = evaluate({data.material,
                                             geometry,
                                             affine,
                                             {},
                                             nullptr,
                                             0,
                                             data.time,
                                             data.volumetric_heat_source,
                                             data.strain_formulation,
                                             false},
                {true, false, true, false});
            double heat = 0;
            for (std::size_t n = 0; n < 4; ++n)
                heat += result.residual[n];
            if (std::abs(heat + 3 * std::acos(-1.0) * 1e4) > 1e-8)
                throw std::runtime_error("QUAD8 source does not integrate the exact annular volume");
            for (std::size_t q = 0; q < geometry.point_count; ++q) {
                const auto& s = result.history[q].stress;
                if (std::abs(s.rr - 800) > 1e-8 || std::abs(s.zz + 1600) > 1e-8 || std::abs(s.hoop - 800) > 1e-8
                    || std::abs(s.rz) > 1e-8)
                    throw std::runtime_error("QUAD8 affine elastic patch mismatch");
            }
            auto heated = affine;
            for (std::size_t n = 0; n < 4; ++n)
                heated[n] = 630;
            fuelsim::elements::Cax8Input input{data.material, geometry, heated, affine};
            input.include_thermal_time_term = true;
            for (double step : {0.1, 0.0, std::numeric_limits<double>::quiet_NaN()}) {
                input.time_step = step;
                for (bool jacobian : {false, true}) {
                    bool rejected = false;
                    try {
                        if (geometry.point_count == 4)
                            (void)fuelsim::elements::evaluate_cax8rt(input, {true, jacobian, false, false});
                        else
                            (void)fuelsim::elements::evaluate_cax8t(input, {true, jacobian, false, false});
                    } catch (const std::invalid_argument&) {
                        rejected = true;
                    }
                    if (!rejected)
                        throw std::runtime_error("CAX8 heat capacity must reject missing committed history");
                }
            }
            input.include_thermal_time_term = false;
            input.time_step = 0.0;
            const auto steady = geometry.point_count == 4 ? fuelsim::elements::evaluate_cax8rt(input)
                                                          : fuelsim::elements::evaluate_cax8t(input);
            if (steady.stored_heat_rate != 0.0)
                throw std::runtime_error("CAX8 steady evaluation without history must have zero stored heat rate");
            const fuelsim::Quad8MaterialHistory history{};
            const auto transient = evaluate({data.material,
                                                geometry,
                                                heated,
                                                affine,
                                                &history,
                                                .1,
                                                data.time,
                                                data.volumetric_heat_source,
                                                data.strain_formulation,
                                                true},
                {true, false, true, false});
            const double capacity = 3 * std::acos(-1.0) * 1000 * 100 * 30 / .1;
            if (std::abs(transient.stored_heat_rate / capacity - 1) > 1e-13)
                throw std::runtime_error("QUAD8 capacity violates the uniform-heating energy balance");
        }
        for (auto element : {selected}) {
            const auto geometry = fuelsim::test::make_cax8_geometry(coordinates, element);
            if (geometry.point_count != (element == fuelsim::RzElementFormulation::cax8rt ? 4U : 9U))
                throw std::runtime_error("QUAD8 integration count mismatch");
            for (int mechanism = 0; mechanism < 4; ++mechanism) {
                auto properties = base;
                if (mechanism & 1)
                    properties = fuelsim::test::with_norton(properties, .01, 1e4, 3);
                if (mechanism & 2)
                    properties = fuelsim::test::with_plasticity(properties, 1e4, 1e5);
                for (auto form : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
                    fuelsim::AxisymmetricTestData data{fuelsim::IsotropicThermoelasticMaterial(properties)};
                    data.strain_formulation = form;
                    data.volumetric_heat_source = 1e4;
                    const fuelsim::Quad8MaterialHistory history{};
                    const auto active = evaluate({data.material,
                                                     geometry,
                                                     state,
                                                     old,
                                                     &history,
                                                     .1,
                                                     data.time,
                                                     data.volumetric_heat_source,
                                                     data.strain_formulation,
                                                     true},
                        {true, true, true, false});
                    const auto passive = evaluate({data.material,
                                                      geometry,
                                                      state,
                                                      old,
                                                      &history,
                                                      .1,
                                                      data.time,
                                                      data.volumetric_heat_source,
                                                      data.strain_formulation,
                                                      true},
                        {true, false, true, false});
                    if (active.residual != passive.residual)
                        throw std::runtime_error("QUAD8 residual depends on Jacobian request");
                    {
                        const fuelsim::elements::Cax8Input input{data.material,
                            geometry,
                            state,
                            old,
                            &history,
                            .1,
                            data.time,
                            data.volumetric_heat_source,
                            form,
                            true};
                        for (const auto quadrature : {fuelsim::elements::Cax8Quadrature::full,
                                 fuelsim::elements::Cax8Quadrature::reduced,
                                 static_cast<fuelsim::elements::Cax8Quadrature>(-1)}) {
                            const bool matches =
                                quadrature
                                == (geometry.point_count == 4 ? fuelsim::elements::Cax8Quadrature::reduced
                                                              : fuelsim::elements::Cax8Quadrature::full);
                            bool rejected = false;
                            try {
                                fuelsim::elements::evaluate_cax8t(input, {false, false, false, false}, quadrature);
                            } catch (const std::invalid_argument&) {
                                rejected = true;
                            }
                            if (rejected == matches)
                                throw std::runtime_error("CAX8 quadrature validation failed");
                        }
                        const auto full = evaluate(input, {true, true, true, true});
                        const auto tangent = evaluate(input, {true, true, false, false});
                        const auto residual_request = evaluate(input, {true, false, false, false});
                        const auto history_request = evaluate(input, {false, false, true, false});
                        const auto stress_request = evaluate(input, {false, false, false, true});
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
                                || stress_request.stress[q].rr != a.stress.rr
                                || stress_request.stress[q].zz != a.stress.zz
                                || stress_request.stress[q].hoop != a.stress.hoop
                                || stress_request.stress[q].rz != a.stress.rz)
                                throw std::runtime_error("Element request changed history or stress");
                        }
                        if (full.generated_heat_rate != residual_request.generated_heat_rate
                            || full.stored_heat_rate != residual_request.stored_heat_rate
                            || full.generated_heat_rate != history_request.generated_heat_rate
                            || full.stored_heat_rate != history_request.stored_heat_rate)
                            throw std::runtime_error("Element request changed heat diagnostics");
                    }
                    for (std::size_t column = 0; column < 20; ++column) {
                        auto plus = state, minus = state;
                        const double h = column < 4 ? 1e-3 : 1e-7;
                        plus[column] += h;
                        minus[column] -= h;
                        const auto rp = evaluate({data.material,
                                                     geometry,
                                                     plus,
                                                     old,
                                                     &history,
                                                     .1,
                                                     data.time,
                                                     data.volumetric_heat_source,
                                                     data.strain_formulation,
                                                     true},
                            {true, false, true, false});
                        const auto rm = evaluate({data.material,
                                                     geometry,
                                                     minus,
                                                     old,
                                                     &history,
                                                     .1,
                                                     data.time,
                                                     data.volumetric_heat_source,
                                                     data.strain_formulation,
                                                     true},
                            {true, false, true, false});
                        double scale = 0, error = 0;
                        for (std::size_t row = 0; row < 20; ++row) {
                            const auto fd = (rp.residual[row] - rm.residual[row]) / (2 * h),
                                       ad = active.jacobian[20 * row + column];
                            error += std::pow(fd - ad, 2);
                            scale += fd * fd;
                        }
                        const double relative = std::sqrt(error / std::max(scale, 1.0));
                        maximum = std::max(maximum, relative);
                        if (relative > 2e-6)
                            throw std::runtime_error("QUAD8 centered derivative mismatch: " + std::to_string(column)
                                                     + " " + std::to_string(relative));
                    }
                    if (!mechanism) {
                        data.volumetric_heat_source = form == fuelsim::StrainFormulation::small ? 0 : 1e6;
                        std::cout << "probe_form " << (form == fuelsim::StrainFormulation::small ? "small" : "finite")
                                  << '\n';
                        const auto probe = evaluate({data.material,
                                                        geometry,
                                                        state,
                                                        old,
                                                        &history,
                                                        .1,
                                                        data.time,
                                                        data.volumetric_heat_source,
                                                        data.strain_formulation,
                                                        true},
                            {true, false, true, false});
                        std::cout << std::setprecision(17);
                        for (std::size_t q = 0; q < geometry.point_count; ++q) {
                            const auto& s = probe.history[q].stress;
                            std::cout << "probe_point " << q + 1 << ' ' << s.rr << ' ' << s.zz << ' ' << s.hoop << ' '
                                      << s.rz << '\n';
                        }
                        for (std::size_t n = 0; n < 8; ++n)
                            std::cout << "probe_node " << n + 1 << ' ' << probe.residual[4 + n] << ' '
                                      << probe.residual[12 + n] << ' ' << (n < 4 ? probe.residual[n] : 0) << '\n';
                    }
                }
            }
        }
        std::cout << "maximum centered derivative error=" << maximum << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
