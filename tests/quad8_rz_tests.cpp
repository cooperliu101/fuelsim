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
        const auto geometry = fuelsim::make_quad8_rz_geometry(coordinates);
        auto base = fuelsim::test::thermoelastic(0, 10, 1e6, .25, 1e-5, 600, 0, 0, 0, 1000, 100);
        fuelsim::Quad8RzValues state = {600, 620, 590, 610, .01, .03, .06, -.02, .025, .038, .022, -.002, .02, -.01,
                                   .07, .05, .003, .034, .065, .032},
                               old{};
        for (std::size_t n = 0; n < 4; ++n) old[n] = 600;
        double maximum = 0;
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
                    throw std::runtime_error("CAX8T residual depends on Jacobian request");
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
                        throw std::runtime_error("CAX8T centered derivative mismatch: " + std::to_string(column) + " " +
                                                 std::to_string(relative));
                }
                if (!mechanism) {
                    data.volumetric_heat_source = form == fuelsim::StrainFormulation::small ? 0 : 1e6;
                    std::cout << "probe_form " << (form == fuelsim::StrainFormulation::small ? "small" : "finite")
                              << '\n';
                    const auto probe = fuelsim::compute_quad8_rz(data, geometry, state, old, &history, .1, false);
                    std::cout << std::setprecision(17);
                    for (std::size_t q = 0; q < 9; ++q) {
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
        std::cout << "maximum centered derivative error=" << maximum << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
