#include "fuelsim/core/rz_quad8.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        const fuelsim::Quad8RzCoordinates coordinates = {
            {{1, 0}, {2.1, .1}, {1.9, 1.2}, {.9, 1}, {1.55, .03}, {2.02, .65}, {1.4, 1.12}, {.93, .5}}};

        auto base = fuelsim::test::thermoelastic(0, 10, 1e6, .25, 1e-5, 600, 0, 0, 0, 1000, 100);
        fuelsim::Quad8RzValues state = {600, 620, 590, 610, .01, .03, .06, -.02, .025, .038, .022, -.002, .02, -.01,
                                   .07, .05, .003, .034, .065, .032},
                               old{};
        for (std::size_t n = 0; n < 4; ++n) old[n] = 600;
        double maximum = 0;
        // An affine annulus has exact volume and uniform strain under either rule.
        const fuelsim::Quad8RzCoordinates annulus = {
            {{1, 0}, {2, 0}, {2, 1}, {1, 1}, {1.5, 0}, {2, .5}, {1.5, 1}, {1, .5}}};
        for (auto element : {fuelsim::RzElementFormulation::cax8t, fuelsim::RzElementFormulation::cax8rt}) {
            const auto geometry = fuelsim::make_quad8_rz_geometry(annulus, element);
            fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(base)};
            data.volumetric_heat_source = 1e4;
            fuelsim::Quad8RzValues affine{};
            for (std::size_t n = 0; n < 4; ++n) affine[n] = 600;
            for (std::size_t n = 0; n < 8; ++n) {
                affine[4 + n] = .001 * annulus[n].r;
                affine[12 + n] = -.002 * annulus[n].z;
            }
            const auto result = fuelsim::compute_quad8_rz(data, geometry, affine, {}, nullptr, 0, false, false);
            double heat = 0;
            for (std::size_t n = 0; n < 4; ++n) heat += result.residual[n];
            if (std::abs(heat + 3 * std::acos(-1.0) * 1e4) > 1e-8)
                throw std::runtime_error("QUAD8 source does not integrate the exact annular volume");
            for (std::size_t q = 0; q < geometry.point_count; ++q) {
                const auto& s = result.history[q].stress;
                if (std::abs(s.rr - 800) > 1e-8 || std::abs(s.zz + 1600) > 1e-8 || std::abs(s.hoop - 800) > 1e-8 ||
                    std::abs(s.rz) > 1e-8)
                    throw std::runtime_error("QUAD8 affine elastic patch mismatch");
            }
            auto heated = affine;
            for (std::size_t n = 0; n < 4; ++n) heated[n] = 630;
            const fuelsim::Quad8MaterialHistory history{};
            const auto transient = fuelsim::compute_quad8_rz(data, geometry, heated, affine, &history, .1, false);
            const double capacity = 3 * std::acos(-1.0) * 1000 * 100 * 30 / .1;
            if (std::abs(transient.stored_heat_rate / capacity - 1) > 1e-13)
                throw std::runtime_error("QUAD8 capacity violates the uniform-heating energy balance");
        }
        for (auto element : {fuelsim::RzElementFormulation::cax8t, fuelsim::RzElementFormulation::cax8rt}) {
            const auto geometry = fuelsim::make_quad8_rz_geometry(coordinates, element);
            if (geometry.point_count != (element == fuelsim::RzElementFormulation::cax8rt ? 4U : 9U))
                throw std::runtime_error("QUAD8 integration count mismatch");
            for (int mechanism = 0; mechanism < 4; ++mechanism) {
                auto properties = base;
                if (mechanism & 1) properties = fuelsim::test::with_norton(properties, .01, 1e4, 3);
                if (mechanism & 2) properties = fuelsim::test::with_plasticity(properties, 1e4, 1e5);
                for (auto form : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
                    fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(properties)};
                    data.strain_formulation = form;
                    data.volumetric_heat_source = 1e4;
                    const fuelsim::Quad8MaterialHistory history{};
                    const auto active = fuelsim::compute_quad8_rz(data, geometry, state, old, &history, .1, true);
                    const auto passive = fuelsim::compute_quad8_rz(data, geometry, state, old, &history, .1, false);
                    if (active.residual != passive.residual)
                        throw std::runtime_error("QUAD8 residual depends on Jacobian request");
                    for (std::size_t column = 0; column < 20; ++column) {
                        auto plus = state, minus = state;
                        const double h = column < 4 ? 1e-3 : 1e-7;
                        plus[column] += h;
                        minus[column] -= h;
                        const auto rp = fuelsim::compute_quad8_rz(data, geometry, plus, old, &history, .1, false);
                        const auto rm = fuelsim::compute_quad8_rz(data, geometry, minus, old, &history, .1, false);
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
                            throw std::runtime_error("QUAD8 centered derivative mismatch: " + std::to_string(column) +
                                                     " " + std::to_string(relative));
                    }
                    if (!mechanism) {
                        data.volumetric_heat_source = form == fuelsim::StrainFormulation::small ? 0 : 1e6;
                        std::cout << "probe_form " << (form == fuelsim::StrainFormulation::small ? "small" : "finite")
                                  << '\n';
                        const auto probe = fuelsim::compute_quad8_rz(data, geometry, state, old, &history, .1, false);
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
