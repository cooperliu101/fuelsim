#include "c3d_common.hpp"
#include "cax_common.hpp"
#include "quad8_shape.hpp"
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

void quad8_polynomials() {
    using namespace fuelsim::quad8_face_detail;
    const std::array<std::array<double, 2>, 8> nodes = {
        {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {1, 0}, {0, 1}, {-1, 0}}};
    for (std::size_t node = 0; node < 8; ++node) {
        DoubleQuad8ShapeValues values;
        double_quad8_shape(nodes[node][0], nodes[node][1], values);
        for (std::size_t i = 0; i < 8; ++i)
            check(values.shape[i] == (i == node ? 1.0 : 0.0), "QUAD8 nodal interpolation is Kronecker delta");
    }
    for (const auto& point : std::array<std::array<double, 2>, 7>{
             {{0, 0}, {0.23, -0.41}, {-0.67, 0.19}, {-1.0, 0.31}, {0.2, 1.0}, {1.2, -0.7}, {-1.3, 1.4}}}) {
        const double xi = point[0], eta = point[1];
        DoubleQuad8ShapeValues value;
        Quad8ShapeValues active;
        double_quad8_shape(xi, eta, value);
        quad8_shape(Scalar::independent(xi, 0, 2), Scalar::independent(eta, 1, 2), active);
        const std::array<const std::array<double, 8>*, 6> ordinary = {&value.shape,
            &value.derivative_xi,
            &value.derivative_eta,
            &value.second_xi,
            &value.second_xi_eta,
            &value.second_eta};
        const std::array<const std::array<Scalar, 8>*, 6> automatic = {&active.shape,
            &active.derivative_xi,
            &active.derivative_eta,
            &active.second_xi,
            &active.second_xi_eta,
            &active.second_eta};
        for (std::size_t field = 0; field < 6; ++field) {
            double sum = 0.0;
            for (std::size_t i = 0; i < 8; ++i) {
                sum += (*ordinary[field])[i];
                check((*ordinary[field])[i] == (*automatic[field])[i].value(),
                    "QUAD8 double and AD polynomial values are identical");
            }
            check(near(sum, field == 0 ? 1.0 : 0.0, 1e-13), "QUAD8 partition of unity and derivative sums hold");
        }
        constexpr double first_step = 1e-6, second_step = 1e-3;
        DoubleQuad8ShapeValues xp, xm, yp, ym, xxp, xxm, yyp, yym, pp, pm, mp, mm;
        double_quad8_shape(xi + first_step, eta, xp);
        double_quad8_shape(xi - first_step, eta, xm);
        double_quad8_shape(xi, eta + first_step, yp);
        double_quad8_shape(xi, eta - first_step, ym);
        double_quad8_shape(xi + second_step, eta, xxp);
        double_quad8_shape(xi - second_step, eta, xxm);
        double_quad8_shape(xi, eta + second_step, yyp);
        double_quad8_shape(xi, eta - second_step, yym);
        double_quad8_shape(xi + second_step, eta + second_step, pp);
        double_quad8_shape(xi + second_step, eta - second_step, pm);
        double_quad8_shape(xi - second_step, eta + second_step, mp);
        double_quad8_shape(xi - second_step, eta - second_step, mm);
        for (std::size_t i = 0; i < 8; ++i) {
            check(near(value.derivative_xi[i], (xp.shape[i] - xm.shape[i]) / (2 * first_step), 2e-9)
                      && near(value.derivative_eta[i], (yp.shape[i] - ym.shape[i]) / (2 * first_step), 2e-9),
                "QUAD8 first derivatives match independent differences of shape values");
            check(near(value.second_xi[i],
                      (xxp.shape[i] - 2 * value.shape[i] + xxm.shape[i]) / (second_step * second_step),
                      2e-8)
                      && near(value.second_eta[i],
                          (yyp.shape[i] - 2 * value.shape[i] + yym.shape[i]) / (second_step * second_step),
                          2e-8)
                      && near(value.second_xi_eta[i],
                          (pp.shape[i] - pm.shape[i] - mp.shape[i] + mm.shape[i]) / (4 * second_step * second_step),
                          2e-8),
                "QUAD8 second derivatives match independent differences of shape values");
            check(near(active.shape[i].derivative(0), value.derivative_xi[i], 1e-13)
                      && near(active.shape[i].derivative(1), value.derivative_eta[i], 1e-13)
                      && near(active.derivative_xi[i].derivative(0), value.second_xi[i], 1e-13)
                      && near(active.derivative_xi[i].derivative(1), value.second_xi_eta[i], 1e-13)
                      && near(active.derivative_eta[i].derivative(0), value.second_xi_eta[i], 1e-13)
                      && near(active.derivative_eta[i].derivative(1), value.second_eta[i], 1e-13),
                "QUAD8 AD first and second derivatives agree with explicit polynomial derivatives");
        }
    }
}

void axisymmetric() {
    const std::array<double, 4> tensor = {2.0, 5.0, 7.0, 3.0};
    const AxisymmetricRotation quarter_turn = {Scalar::independent(0.0, 0, 1), -1.0, 1.0, 0.0, 1.0};
    check(rotate_axisymmetric_tensor_values(tensor, quarter_turn) == std::array<double, 4>{5.0, 2.0, 7.0, -3.0},
        "Value-only tensor rotation swaps in-plane axes and reverses tensor shear");
    MaterialPointState history;
    history.elastic_strain = tensor;
    history.plastic_strain = {1.0, -2.0, 1.0, 3.0};
    history.creep_strain = {-3.0, 2.0, 1.0, -2.0};
    history.equivalent_plastic_strain = 4.0;
    history.equivalent_creep_strain = 6.0;
    rotate_axisymmetric_strain_history(history, quarter_turn);
    check(history.elastic_strain == std::array<double, 4>{5.0, 2.0, 7.0, -3.0}
              && history.plastic_strain == std::array<double, 4>{-2.0, 1.0, 1.0, -3.0}
              && history.creep_strain == std::array<double, 4>{2.0, -3.0, 1.0, 2.0}
              && history.equivalent_plastic_strain == 4.0 && history.equivalent_creep_strain == 6.0,
        "History rotation preserves hoop components and equivalent scalar histories");

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
    quad8_polynomials();
    return failures == 0 ? 0 : 1;
}
