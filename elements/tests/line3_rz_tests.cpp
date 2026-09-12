#include "cax8_types.hpp"
#include "contact_types.hpp"
#include "line3_rz.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
void test_nodal_heat_contact() {
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
}

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
        test_nodal_heat_contact();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
