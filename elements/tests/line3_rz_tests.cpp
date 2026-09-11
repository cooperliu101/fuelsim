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
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
