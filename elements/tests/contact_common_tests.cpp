#include "contact_common.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using adlite::Scalar;
using namespace fuelsim;
using namespace fuelsim::contact_common;
int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << message << '\n';
        ++failures;
    }
}

void domain_error(const std::function<void()>& operation, const char* message) {
    try {
        operation();
        check(false, message);
    } catch (const std::domain_error&) {
    }
}

bool near(double a, double b) {
    return std::abs(a - b) <= 1e-8 * std::max(1.0, std::abs(b));
}

double derivative(const Scalar& value) {
    return value.derivative_size() == 0 ? 0.0 : value.derivative(0);
}

void projection() {
    check(projection_increment_converged(0.0, 0.0, 1.0, 2.0, 32), "Passive projection converges");
    check(!projection_increment_converged(1e-6, 0.0, 1.0, 2.0, 32), "Finite large increment does not converge");
    const double zero = 0.0, large = 1e-5;
    const Scalar z = Scalar::seeded(0.0, &zero, 1);
    check(!projection_increment_converged(Scalar::seeded(0.0, &large, 1), z, z, z, 56),
        "Small primal increment cannot hide an unconverged derivative");
    for (double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        for (std::size_t index = 0; index < 4; ++index) {
            std::array<Scalar, 4> values{z, z, z, z};
            values[index] = Scalar::seeded(0.0, &invalid, 1);
            domain_error([&] { projection_increment_converged(values[0], values[1], values[2], values[3], 56); },
                "Nonfinite increment and coordinate derivatives must be rejected");
            values[index] = Scalar::seeded(invalid, &zero, 1);
            domain_error([&] { projection_increment_converged(values[0], values[1], values[2], values[3], 56); },
                "Nonfinite increment and coordinate values must be rejected");
        }
    }
    for (bool oversized : {false, true}) {
        bool rejected = false;
        try {
            projection_increment_converged(z, z, z, oversized ? z : Scalar(0.0), oversized ? 0 : 56);
        } catch (const std::logic_error&) {
            rejected = true;
        }
        check(rejected, "Mixed passive/active or excessive derivative widths must be rejected");
    }
}

void heat() {
    GapHeatProperties properties{2.0, 0.5};
    for (double gap : {0.25, 0.5, 1.0}) {
        const Scalar h = gap_conductance(properties, Scalar::independent(gap, 0, 1), 400.0, 300.0);
        const double step = 1e-6;
        const double centered_derivative = (gap_conductance(properties, gap + step, 400.0, 300.0).value()
                                               - gap_conductance(properties, gap - step, 400.0, 300.0).value())
                                           / (2 * step);
        check(std::abs(derivative(h) - centered_derivative) < 1e-5,
            "Gas gap derivative including corner matches centered slopes");
    }
    properties.law = GapHeatConductanceLaw::affine;
    properties.conductance = 10.0;
    properties.clearance_derivative = 2.0;
    properties.pressure_derivative = 3.0;
    properties.contact_penalty = 4.0;
    properties.temperature_derivative = 0.1;
    properties.reference_temperature = 350.0;
    const Scalar h = gap_conductance(properties, Scalar::independent(0.0, 0, 1), 400.0, 300.0);
    check(near(h.value(), 10.0) && near(derivative(h), -4.0), "Affine corner uses mean of open and closed slopes");
    check(near(gap_conductance(properties, 0.0, Scalar::independent(400.0, 0, 1), 300.0).derivative(0), 0.05),
        "Average temperature retains both temperature derivatives");
    for (double invalid : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        properties.conductance = invalid;
        domain_error([&] { gap_conductance(properties, 0.0, 400.0, 300.0); }, "Invalid affine conductance is rejected");
        GapHeatProperties gas{invalid, 0.5};
        domain_error([&] { gap_conductance(gas, 0.5, 400.0, 300.0); }, "Invalid gas conductance is rejected");
    }
    domain_error([&] { gap_conductance({2.0, 0.5}, std::numeric_limits<double>::quiet_NaN(), 400.0, 300.0); },
        "Max law cannot hide a nonfinite gap");
}

void friction() {
    const ActivePoint3 normal{0.0, 0.0, 1.0}, zero{};
    NormalContactProperties properties{100.0, 0.5};
    ContactPointHistory history;
    for (double slip : {0.01, 0.2}) {
        const auto result = friction_return(properties,
            history,
            zero,
            zero,
            {Scalar::independent(slip, 0, 1), 0.0, 0.0},
            normal,
            10.0,
            2.0);
        const double expected = slip < 0.05 ? 100.0 * slip : 5.0;
        check(near(result.tangential_traction.value(), expected), "Coulomb return satisfies traction limit");
        check(near(result.tangential_force.value(), 2.0 * expected), "Friction force includes tributary area");
        check(near(result.friction_dissipation.value(), slip < 0.05 ? 0.0 : 1.5),
            "Slip return retains dissipated work");
        const double step = 1e-7;
        const auto plus = friction_return(properties, history, zero, zero, {slip + step, 0.0, 0.0}, normal, 10.0, 2.0);
        const auto minus = friction_return(properties, history, zero, zero, {slip - step, 0.0, 0.0}, normal, 10.0, 2.0);
        check(near(derivative(result.tangential_force),
                  (plus.tangential_force.value() - minus.tangential_force.value()) / (2 * step)),
            "Stick and slip tangents match centered differences");
    }
    const auto threshold = friction_return(properties, history, zero, zero, {0.05, 0.0, 0.0}, normal, 10.0, 2.0);
    history.sliding = true;
    check(!threshold.sliding
              && friction_return(properties, history, zero, zero, {0.05, 0.0, 0.0}, normal, 10.0, 2.0).sliding,
        "Equality at Coulomb limit preserves the committed branch");
    const auto open = friction_return(properties, history, zero, {0.3, 0.0, 0.0}, {0.2, 0.0, 0.0}, normal, 0.0, 2.0);
    check(near(open.tangential_slip[0].value(), 0.3) && open.tangential_force.value() == 0.0,
        "Open contact freezes total slip and contributes no friction");
    for (double invalid :
        {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        properties.penalty = invalid;
        domain_error([&] { friction_return(properties, history, zero, zero, {0.01, 0.0, 0.0}, normal, 10.0, 2.0); },
            "Invalid stick stiffness is rejected");
    }
    properties = {100.0, 0.5, false, std::numeric_limits<double>::denorm_min()};
    domain_error([&] { friction_return(properties, history, zero, zero, {0.01, 0.0, 0.0}, normal, 10.0, 2.0); },
        "Finite slip parameter that overflows stiffness is rejected");
}
} // namespace

int main() {
    projection();
    heat();
    friction();
    return failures == 0 ? 0 : 1;
}
