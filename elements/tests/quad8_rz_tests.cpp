#include "line3_rz_boundary.hpp"
#include "line3_rz_contact.hpp"
#include "rz_quad8.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
void test_quadratic_boundary() {
    using fuelsim::elements::compute_line3_rz_boundary;
    using fuelsim::elements::Line3RzBoundaryKind;
    const std::array<fuelsim::RzPoint, 3> curved = {{{1, 0}, {1.2, 2}, {1.3, .8}}};
    const std::vector<double> state = {500, 650, .02, -.01, .03, .01, .04, -.02};
    for (const bool current : {false, true})
        for (const auto kind : {Line3RzBoundaryKind::pressure,
                 Line3RzBoundaryKind::traction,
                 Line3RzBoundaryKind::heat_flux,
                 Line3RzBoundaryKind::convection})
            for (const auto component : {fuelsim::TractionComponent::radial, fuelsim::TractionComponent::axial}) {
                const fuelsim::elements::Line3RzBoundaryData data{kind, component, 7.0, 300.0, current};
                const auto active = compute_line3_rz_boundary(data, curved, state, true);
                const auto passive = compute_line3_rz_boundary(data, curved, state, false);
                if (active.residual != passive.residual)
                    throw std::runtime_error("Quadratic boundary active/passive residual mismatch");
                for (std::size_t column = 0; column < 8; ++column) {
                    auto plus = state, minus = state;
                    const double step = column < 2 ? 1e-3 : 1e-6;
                    plus[column] += step;
                    minus[column] -= step;
                    const auto p = compute_line3_rz_boundary(data, curved, plus);
                    const auto m = compute_line3_rz_boundary(data, curved, minus);
                    for (std::size_t row = 0; row < 8; ++row) {
                        const double fd = (p.residual[row] - m.residual[row]) / (2 * step);
                        if (std::abs(fd - active.jacobian[8 * row + column]) > 2e-6 * (1 + std::abs(fd)))
                            throw std::runtime_error("Quadratic boundary centered derivative mismatch");
                    }
                }
            }
    const std::array<fuelsim::RzPoint, 3> straight = {{{2, 0}, {2, 3}, {2, 1.5}}};
    const std::vector<double> constant = {400, 400, 0, 0, 0, 0, 0, 0};
    const double area = 12 * std::acos(-1.0);
    for (const auto kind : {Line3RzBoundaryKind::pressure,
             Line3RzBoundaryKind::traction,
             Line3RzBoundaryKind::heat_flux,
             Line3RzBoundaryKind::convection}) {
        const fuelsim::elements::Line3RzBoundaryData data{kind, fuelsim::TractionComponent::radial, 7, 300, false};
        const auto result = compute_line3_rz_boundary(data, straight, constant);
        const bool mechanical = kind == Line3RzBoundaryKind::pressure || kind == Line3RzBoundaryKind::traction;
        const double sum = mechanical ? result.residual[2] + result.residual[3] + result.residual[4]
                                      : result.residual[0] + result.residual[1];
        const double expected = area
                                * (kind == Line3RzBoundaryKind::pressure      ? 7
                                    : kind == Line3RzBoundaryKind::convection ? 700
                                                                              : -7);
        if (std::abs(sum - expected) > 1e-13 * std::abs(expected))
            throw std::runtime_error("Quadratic cylinder boundary resultant mismatch");
    }
}
} // namespace

int main() {
    try {
        test_quadratic_boundary();
        // Nonuniform temperature and genuinely curved geometry distinguish NTS
        // heat transfer from surface integration and from mechanical nodal area.
        for (std::size_t node = 0; node < 2; ++node) {
            fuelsim::rz8::Line3ContactGeometry contact;
            contact.secondary = {{{1, .2}, {1.2, 1.2}, {1.15, .7}}};
            contact.primary = {{{1.8, 2}, {1.4, 0}, {1.7, 1}}};
            contact.primary_first = contact.primary_last = true;
            contact.coordinate = node == 0 ? -1 : 1;
            contact.secondary_node = node;
            contact.nodal_heat = true;
            std::vector<double> v = {500, 650, 300, 350, .01, .02, .03, .01, .02, .03, .03, -.01, .02, -.02, .01, .01};
            const fuelsim::GapHeatProperties heat{.2, 1e-5};
            const fuelsim::NormalContactProperties mechanical{1e8};
            const auto active = fuelsim::rz8::compute_line3_contact(contact, heat, mechanical, v, v, {}, true);
            const auto passive = fuelsim::rz8::compute_line3_contact(contact, heat, mechanical, v, v, {}, false);
            if (!active.thermal.projected || active.residual != passive.residual)
                throw std::runtime_error("Quadratic NTS contact projection or passive residual failed");
            double sum = 0, scale = 0;
            for (std::size_t row = 0; row < 4; ++row) {
                sum += active.residual[row];
                scale += std::abs(active.residual[row]);
            }
            if (std::abs(sum) > 1e-13 * scale)
                throw std::runtime_error("Quadratic NTS heat is not conserved");
            for (std::size_t column = 0; column < 16; ++column) {
                auto plus = v, minus = v;
                const double step = column < 4 ? 1e-3 : 1e-6;
                plus[column] += step;
                minus[column] -= step;
                const auto p = fuelsim::rz8::compute_line3_contact(contact, heat, mechanical, plus, v, {}, false);
                const auto m = fuelsim::rz8::compute_line3_contact(contact, heat, mechanical, minus, v, {}, false);
                for (std::size_t row = 0; row < 16; ++row) {
                    const double fd = (p.residual[row] - m.residual[row]) / (2 * step),
                                 ad = active.jacobian[16 * row + column];
                    if (std::abs(fd - ad) > 2e-6 * (1 + std::abs(fd)))
                        throw std::runtime_error("Quadratic NTS heat Jacobian differs from centered differences");
                }
            }
            auto moved = v;
            moved[6] += .02;
            moved[9] -= .01;
            const auto midnode = fuelsim::rz8::compute_line3_contact(contact, heat, mechanical, moved, v, {}, false);
            if (midnode.residual != passive.residual)
                throw std::runtime_error("NTS corner heat transfer incorrectly depends on secondary midnode");
        }
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
        for (auto element : {fuelsim::RzElementFormulation::cax8t, fuelsim::RzElementFormulation::cax8rt}) {
            const auto geometry = fuelsim::make_quad8_rz_geometry(annulus, element);
            fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(base)};
            data.volumetric_heat_source = 1e4;
            fuelsim::Quad8RzValues affine{};
            for (std::size_t n = 0; n < 4; ++n)
                affine[n] = 600;
            for (std::size_t n = 0; n < 8; ++n) {
                affine[4 + n] = .001 * annulus[n].r;
                affine[12 + n] = -.002 * annulus[n].z;
            }
            const auto result = fuelsim::compute_quad8_rz(data, geometry, affine, {}, nullptr, 0, false, false);
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
                if (mechanism & 1)
                    properties = fuelsim::test::with_norton(properties, .01, 1e4, 3);
                if (mechanism & 2)
                    properties = fuelsim::test::with_plasticity(properties, 1e4, 1e5);
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
                            throw std::runtime_error("QUAD8 centered derivative mismatch: " + std::to_string(column)
                                                     + " " + std::to_string(relative));
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
