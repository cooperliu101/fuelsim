#include "c3d_common.hpp"
#include "cax_common.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>

namespace {
using namespace fuelsim;
using namespace fuelsim::cartesian_detail;
using namespace fuelsim::c3d8_detail;
using adlite::Scalar;
int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        ++failures;
        std::cerr << message << '\n';
    }
}

bool near(double a, double b, double tolerance = 1e-8) {
    return std::abs(a - b) <= tolerance * std::max(1.0, std::abs(b));
}

void rejects(const std::function<void()>& operation) {
    bool rejected = false;
    try {
        operation();
    } catch (const std::domain_error&) {
        rejected = true;
    }
    check(rejected, "Invalid local geometry must throw a domain error");
}

void central_gradient() {
    const Matrix3 old = {{{1.02, 0.01, -0.02}, {0.03, 0.99, 0.01}, {0.0, 0.02, 1.01}}};
    const Matrix3 current = {{{1.12, 0.07, -0.01}, {0.02, 0.95, 0.04}, {0.03, -0.01, 1.04}}};
    Matrix3 plus = current, minus = current, midpoint{}, direction{};
    ActiveMatrix3 active{}, active_midpoint{};
    constexpr double step = 1e-6;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const double d = 0.1 * std::sin(static_cast<double>(3 * i + j + 1));
            direction[i][j] = d;
            active[i][j] = Scalar::seeded(current[i][j], &d, 1);
            plus[i][j] += step * d;
            minus[i][j] -= step * d;
        }
    const auto gradient = central_increment_gradient(current, old, &midpoint);
    const auto ad = central_increment_gradient(active, old, &active_midpoint);
    const auto forward = central_increment_gradient(plus, old), backward = central_increment_gradient(minus, old);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            double recovered = 0.0, inverse_check = 0.0;
            for (std::size_t k = 0; k < 3; ++k) {
                recovered += gradient[i][k] * (current[k][j] + old[k][j]);
                inverse_check += midpoint[i][k] * (current[k][j] + old[k][j]) / 2.0;
            }
            check(near(recovered, 2.0 * (current[i][j] - old[i][j]), 1e-13),
                "Central gradient satisfies its defining matrix equation");
            check(near(inverse_check, i == j ? 1.0 : 0.0, 1e-13), "Returned midpoint inverse matches midpoint mapping");
            check(gradient[i][j] == ad[i][j].value(), "Central gradient uses identical double and AD value arithmetic");
            check(near(ad[i][j].derivative(0), (forward[i][j] - backward[i][j]) / (2 * step)),
                "Central gradient derivative matches centered differences");
        }
    Matrix3 identity{};
    for (std::size_t i = 0; i < 3; ++i)
        identity[i][i] = 1.0;
    Matrix3 singular = identity;
    singular[0][0] = -1.0;
    rejects([&] { central_increment_gradient(singular, identity); });
    singular[0][0] = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { central_increment_gradient(singular, identity); });
}

void hourglass() {
    std::array<std::array<double, 3>, 8> coordinates = hex8_signs, gradient{}, plus_coordinates{}, minus_coordinates{},
                                         plus_gradient{}, minus_gradient{};
    std::array<std::array<Scalar, 3>, 8> active_coordinates{}, active_gradient{};
    constexpr double step = 1e-6;
    for (std::size_t n = 0; n < 8; ++n)
        for (std::size_t c = 0; c < 3; ++c) {
            gradient[n][c] = hex8_signs[n][c] / 8.0;
            const double dx = 0.03 * std::sin(static_cast<double>(3 * n + c + 1)),
                         dg = 0.02 * std::cos(static_cast<double>(3 * n + c + 1));
            active_coordinates[n][c] = Scalar::seeded(coordinates[n][c], &dx, 1);
            active_gradient[n][c] = Scalar::seeded(gradient[n][c], &dg, 1);
            plus_coordinates[n][c] = coordinates[n][c] + step * dx;
            minus_coordinates[n][c] = coordinates[n][c] - step * dx;
            plus_gradient[n][c] = gradient[n][c] + step * dg;
            minus_gradient[n][c] = gradient[n][c] - step * dg;
        }
    const auto shape = hex8_hourglass_shape(coordinates, gradient),
               forward = hex8_hourglass_shape(plus_coordinates, plus_gradient),
               backward = hex8_hourglass_shape(minus_coordinates, minus_gradient);
    const auto ad = hex8_hourglass_shape(active_coordinates, active_gradient);
    for (std::size_t mode = 0; mode < 4; ++mode) {
        double sum = 0.0;
        for (std::size_t n = 0; n < 8; ++n) {
            sum += shape[n][mode];
            check(shape[n][mode] == ad[n][mode].value(), "Hourglass projection retains double and AD value agreement");
            check(near(ad[n][mode].derivative(0), (forward[n][mode] - backward[n][mode]) / (2 * step)),
                "Hourglass shape derivatives match centered differences");
        }
        check(std::abs(sum) < 1e-14, "Hourglass modes annihilate constant displacement");
        for (std::size_t c = 0; c < 3; ++c) {
            double affine = 0.0;
            for (std::size_t n = 0; n < 8; ++n)
                affine += shape[n][mode] * coordinates[n][c];
            check(std::abs(affine) < 1e-14, "Hourglass modes annihilate affine displacement");
        }
    }
    const double dv = 0.4;
    const auto metric = reduced_hex8_metric(gradient);
    const auto coefficients = reduced_hex8_thermal_hourglass_coefficients(metric, 8.0);
    const auto ordinary = reduced_hex8_thermal_hourglass_coefficients(gradient, 8.0);
    const auto ad_coefficients =
        reduced_hex8_thermal_hourglass_coefficients(active_gradient, Scalar::seeded(8.0, &dv, 1));
    const auto cp = reduced_hex8_thermal_hourglass_coefficients(plus_gradient, 8.0 + step * dv),
               cm = reduced_hex8_thermal_hourglass_coefficients(minus_gradient, 8.0 - step * dv);
    for (std::size_t mode = 0; mode < 4; ++mode) {
        check(coefficients[mode] == ordinary[mode] && coefficients[mode] == ad_coefficients[mode].value(),
            "Cached metric and AD thermal coefficients agree");
        check(near(ad_coefficients[mode].derivative(0), (cp[mode] - cm[mode]) / (2 * step)),
            "Thermal coefficient derivatives match centered differences");
    }
    gradient = {};
    rejects([&] { reduced_hex8_metric(gradient); });
}

void axisymmetric() {
    const std::array<Scalar, 4> sum = {2.2, 0.1, 0.04, 1.9}, difference = {0.2, 0.1, 0.04, -0.1};
    const auto result = evaluate_axisymmetric_midpoint_increment(sum, difference, 2.1, 0.1);
    check(near(result.determinant_sum.value(), 2.2 * 1.9 - 0.1 * 0.04, 1e-14),
        "Axisymmetric midpoint determinant preserves physical measure");
    check(near(result.hoop.value(), 0.2 / 2.1, 1e-14),
        "Axisymmetric hoop increment uses supplied native sum and difference");
    rejects([&] { evaluate_axisymmetric_midpoint_increment({0.0, 0.0, 0.0, 0.0}, difference, 2.0, 0.0); });
    rejects([&] { evaluate_axisymmetric_midpoint_increment(sum, difference, 0.0, 0.0); });
}
} // namespace

int main() {
    central_gradient();
    hourglass();
    axisymmetric();
    return failures == 0 ? 0 : 1;
}
