#include "ring_gps.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using namespace fuelsim;
using namespace fuelsim::elements;
int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) {
        std::cerr << message << '\n';
        ++failures;
    }
}

bool near(double a, double b, double tolerance = 2e-10) {
    return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

bool same_history(const ContactPointHistory& a, const ContactPointHistory& b) {
    return a.elastic_tangential_slip == b.elastic_tangential_slip && a.sliding == b.sliding
           && a.normal_multiplier == b.normal_multiplier
           && a.cartesian_elastic_tangential_slip == b.cartesian_elastic_tangential_slip
           && a.cartesian_total_tangential_slip == b.cartesian_total_tangential_slip
           && a.cartesian_tangent_basis_initialized == b.cartesian_tangent_basis_initialized
           && a.cartesian_contact_normal == b.cartesian_contact_normal
           && a.cartesian_contact_tangent_first == b.cartesian_contact_tangent_first
           && a.total_tangential_slip == b.total_tangential_slip;
}

void rejected(const std::function<void()>& operation, const char* message) {
    bool caught = false;
    try {
        operation();
    } catch (const std::domain_error&) {
        caught = true;
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    check(caught, message);
}

struct Fixture final {
    RingGpsGeometry _geometry{1.1, 1.0, 0.0, 2.0};
    RingGpsLocalValues _state{300.0, 400.0, 0.0, 0.2, 0.0, 0.0, 0.01, 0.01};
    RingGpsLocalValues _committed_state{300.0, 400.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0};
    RingGpsHistory _history{};
    StrainFormulation _strain_formulation = StrainFormulation::small;

    RingGpsInput input() const {
        RingGpsInput value{_geometry, _state, _committed_state, _history};
        value.thermal = true;
        value.mechanical = true;
        value.heat = {2.0, 0.05};
        value.normal = {1000.0, 0.5};
        value.strain_formulation = _strain_formulation;
        return value;
    }

    void commit(const RingGpsResult& result) {
        _committed_state = _state;
        _history = result.history;
    }
};

void conservation_and_requests() {
    Fixture data;
    auto input = data.input();
    const double area = 4.0 * std::acos(-1.0);
    const auto result = evaluate_ring_gps(input, {true, true, true, true});
    check(near(result.heat_rate, 4000.0 * area), "Thermal contact uses the common reference cylindrical area");
    check(near(result.normal_force, 100.0 * area), "Normal contact uses the common reference cylindrical area");
    check(near(result.tangential_force, 10.0 * area), "Tangential contact uses both axial Gauss points");
    check(result.residual[0] == -result.residual[1], "Thermal contact transfers exactly opposite heat rates");
    check(result.residual[2] == -result.residual[3], "Normal contact transfers exactly opposite radial forces");
    check(result.residual[4] == -result.residual[6] && result.residual[5] == -result.residual[7],
        "Tangential contact transfers exactly opposite endpoint forces");
    for (std::size_t column = 0; column < 8; ++column) {
        check(result.jacobian[column] == -result.jacobian[8 + column], "Heat tangent conserves transferred heat");
        check(result.jacobian[16 + column] == -result.jacobian[24 + column], "Normal tangent conserves radial force");
        check(result.jacobian[32 + column] == -result.jacobian[48 + column]
                  && result.jacobian[40 + column] == -result.jacobian[56 + column],
            "Friction tangent conserves axial force");
    }
    const auto residual_only = evaluate_ring_gps(input, {true, false, false, false});
    const auto tangent_only = evaluate_ring_gps(input, {false, true, false, false});
    const auto history_only = evaluate_ring_gps(input, {false, false, true, false});
    check(result.residual == residual_only.residual, "Passive and differentiated residual paths agree exactly");
    check(tangent_only.residual == RingGpsLocalResidual{}, "Unrequested residual stays zero");
    check(residual_only.jacobian == RingGpsLocalJacobian{}, "Unrequested tangent stays zero");
    check(tangent_only.jacobian == result.jacobian, "Tangent requests are independent of residual requests");
    check(history_only.residual == RingGpsLocalResidual{} && history_only.jacobian == RingGpsLocalJacobian{},
        "History-only evaluation does not populate residual or tangent");
    check(near(history_only.history[0].elastic_tangential_slip, 0.01),
        "History is available without a residual request");
    check(residual_only.history[0].elastic_tangential_slip == 0.0, "Unrequested trial history stays empty");
    input.mechanical = false;
    const auto thermal = evaluate_ring_gps(input, {true, true, true, false});
    check(thermal.normal_force == 0.0 && thermal.tangential_force == 0.0 && thermal.heat_rate > 0.0,
        "Thermal-only contact contributes no mechanical force");
    input.mechanical = true;
    input.thermal = false;
    const auto mechanical = evaluate_ring_gps(input, {true, true, true, false});
    check(mechanical.heat_rate == 0.0 && mechanical.normal_force > 0.0, "Mechanical-only contact contributes no heat");
}

void tangent(Fixture& data, bool affine, double maximum_elastic_slip) {
    auto input = data.input();
    input.normal.maximum_elastic_slip = maximum_elastic_slip;
    if (affine) {
        input.heat.law = GapHeatConductanceLaw::affine;
        input.heat.conductance = 30.0;
        input.heat.clearance_derivative = 3.0;
        input.heat.pressure_derivative = 0.2;
        input.heat.temperature_derivative = 0.03;
        input.heat.reference_temperature = 350.0;
        input.heat.contact_penalty = input.normal.penalty;
    }
    const auto result = evaluate_ring_gps(input, {true, true, true, false});
    check(result.residual == evaluate_ring_gps(input, {true, false, false, false}).residual,
        "Small and finite contact residuals match differentiated evaluations exactly");
    const double step = 1e-6;
    for (std::size_t column = 0; column < 8; ++column) {
        const double saved = data._state[column];
        data._state[column] = saved + step;
        const auto plus = evaluate_ring_gps(input, {true, false, false, false});
        data._state[column] = saved - step;
        const auto minus = evaluate_ring_gps(input, {true, false, false, false});
        data._state[column] = saved;
        for (std::size_t row = 0; row < 8; ++row)
            check(near(result.jacobian[row * 8 + column],
                      (plus.residual[row] - minus.residual[row]) / (2.0 * step),
                      2e-6),
                "Compact contact tangent matches every centered finite-difference column");
    }
    const RingGpsLocalValues direction{1.3, -0.7, 0.02, -0.03, 0.1, 0.05, -0.04, 0.06};
    const auto state = data._state;
    for (std::size_t i = 0; i < 8; ++i)
        data._state[i] = state[i] + step * direction[i];
    const auto plus = evaluate_ring_gps(input, {true, false, false, false});
    for (std::size_t i = 0; i < 8; ++i)
        data._state[i] = state[i] - step * direction[i];
    const auto minus = evaluate_ring_gps(input, {true, false, false, false});
    data._state = state;
    for (std::size_t row = 0; row < 8; ++row) {
        double product = 0.0;
        for (std::size_t column = 0; column < 8; ++column)
            product += result.jacobian[row * 8 + column] * direction[column];
        check(near(product, (plus.residual[row] - minus.residual[row]) / (2.0 * step), 2e-6),
            "Compact contact tangent matches a coupled temperature-radial-axial direction");
    }
}

void jacobians() {
    for (const auto formulation : {StrainFormulation::small, StrainFormulation::finite})
        for (const bool affine : {false, true})
            for (const double maximum_elastic_slip : {0.0, 0.02}) {
                Fixture data;
                data._strain_formulation = formulation;
                tangent(data, affine, maximum_elastic_slip);
                data._state[6] = 0.2;
                data._state[7] = 0.3;
                tangent(data, affine, maximum_elastic_slip);
                data._state[6] = -0.2;
                data._state[7] = -0.3;
                tangent(data, affine, maximum_elastic_slip);
                data._state[3] = -0.2;
                tangent(data, affine, maximum_elastic_slip);
            }
}

void history_transactions(StrainFormulation formulation) {
    Fixture data;
    data._strain_formulation = formulation;
    const auto old_state = data._committed_state;
    data._state[6] = 0.2;
    data._state[7] = 0.2;
    const auto sliding = evaluate_ring_gps(data.input(), {true, true, true, false});
    check(sliding.points[0].sliding && near(sliding.points[0].signed_tangential_traction, 50.0),
        "Positive sliding saturates at the Coulomb pressure bound");
    const double radius = formulation == StrainFormulation::small ? 1.0 : 1.2;
    check(near(sliding.points[0].elastic_tangential_slip, 0.05)
              && near(sliding.friction_dissipation, 50.0 * 0.15 * 4.0 * std::acos(-1.0) * radius),
        "Sliding retains only elastic slip and records nonnegative dissipated work");
    check(data._history[0].elastic_tangential_slip == 0.0 && data._history[0].total_tangential_slip == 0.0
              && !data._history[0].sliding && data._committed_state == old_state,
        "Trial evaluation does not mutate committed nodal values or contact histories");
    const auto repeated = evaluate_ring_gps(data.input(), {true, true, true, false});
    check(repeated.residual == sliding.residual && repeated.jacobian == sliding.jacobian,
        "Repeated contact evaluation is independent of abandoned trial states");
    data.commit(sliding);
    const auto committed_history = data._history;
    data._state[6] = 0.19;
    data._state[7] = 0.19;
    const auto unloading = evaluate_ring_gps(data.input());
    check(!unloading.points[0].sliding && near(unloading.points[0].signed_tangential_traction, 40.0),
        "Reverse loading initially unloads the stored tangential spring");
    data._state[6] = -0.15;
    data._state[7] = -0.15;
    const auto reverse = evaluate_ring_gps(data.input());
    check(reverse.points[0].sliding && near(reverse.points[0].signed_tangential_traction, -50.0)
              && reverse.friction_dissipation > 0.0,
        "Sufficient reverse motion produces negative sliding with positive dissipation");
    check(same_history(data._history[0], committed_history[0]) && same_history(data._history[1], committed_history[1]),
        "Reverse trial loading preserves every committed scalar, vector, and tangent-basis history field");
    data._state[3] = -0.1;
    data._state[6] = 1.2;
    data._state[7] = 1.2;
    const auto open = evaluate_ring_gps(data.input());
    check(open.normal_force == 0.0 && open.tangential_force == 0.0 && open.friction_dissipation == 0.0,
        "Open contact releases mechanical traction and dissipates no frictional work");
    check(near(open.history[0].total_tangential_slip, 0.2) && open.history[0].elastic_tangential_slip == 0.0
              && !open.history[0].sliding,
        "Opening freezes accumulated slip and clears elastic tangential history");
    data.commit(open);
    data._state[3] = 0.2;
    data._state[6] = 1.21;
    data._state[7] = 1.21;
    const auto recontact = evaluate_ring_gps(data.input());
    check(near(recontact.history[0].total_tangential_slip, 0.21)
              && near(recontact.history[0].elastic_tangential_slip, 0.01)
              && near(recontact.points[0].signed_tangential_traction, 10.0),
        "Recontact does not replay accumulated free motion or the abandoned sliding spring");
    check(near(recontact.history[0].cartesian_total_tangential_slip[2], recontact.history[0].total_tangential_slip)
              && near(recontact.history[0].cartesian_elastic_tangential_slip[2],
                  recontact.history[0].elastic_tangential_slip)
              && recontact.history[0].cartesian_contact_normal == std::array<double, 3>{1.0, 0.0, 0.0}
              && recontact.history[0].cartesian_contact_tangent_first == std::array<double, 3>{0.0, 0.0, 1.0},
        "Trial history records the fixed radial normal and positive axial tangent consistently");
}

void antisymmetric_slip() {
    Fixture data;
    data._state[6] = 0.02;
    data._state[7] = -0.02;
    const auto result = evaluate_ring_gps(data.input(), {true, true, true, false});
    const double expected = 1000.0 * 4.0 * std::acos(-1.0) * 0.02 / 6.0;
    check(near(result.residual[6], expected) && near(result.residual[7], -expected),
        "Two axial quadrature points retain the antisymmetric endpoint slip mode");
    check(result.points[0].signed_tangential_traction > 0.0 && result.points[1].signed_tangential_traction < 0.0,
        "Opposite local slip directions survive a zero net axial friction force");
    check(near(result.tangential_force, 0.0), "Antisymmetric contact has zero integrated tangential force");
    tangent(data, true, 0.0);
}

void invalid_inputs() {
    Fixture data;
    auto input = data.input();
    input.normal.augmented_lagrangian = true;
    rejected([&] { evaluate_ring_gps(input); }, "Augmented ring contact must be explicitly rejected");
    input.normal.augmented_lagrangian = false;
    for (const double invalid :
        {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        input.normal.penalty = invalid;
        rejected([&] { evaluate_ring_gps(input); }, "Invalid normal penalty must be rejected");
    }
    input.normal.penalty = 1000.0;
    input.heat.minimum_gap = 0.0;
    rejected([&] { evaluate_ring_gps(input); }, "Gas conductance requires a positive minimum gap");
    input.heat.minimum_gap = 0.05;
    data._geometry.z_upper = data._geometry.z_lower;
    rejected([&] { evaluate_ring_gps(input); }, "Empty axial contact intervals must be rejected");
    data._geometry.z_upper = 2.0;
    data._state[2] = -2.0;
    rejected([&] { evaluate_ring_gps(input); }, "A nonpositive displaced contact radius must be rejected");
    data._state[2] = 0.0;
    data._history[0].elastic_tangential_slip = std::numeric_limits<double>::quiet_NaN();
    rejected([&] { evaluate_ring_gps(input); }, "Nonfinite committed history must be rejected");
    data._history[0].elastic_tangential_slip = 0.0;
    data._committed_state[5] = std::numeric_limits<double>::infinity();
    rejected([&] { evaluate_ring_gps(input); }, "Nonfinite committed axial states must be rejected");
}

void finite_geometry() {
    Fixture data;
    data._strain_formulation = StrainFormulation::finite;
    data._state[4] = -0.03;
    data._state[5] = 0.04;
    data._state[6] = 0.01;
    data._state[7] = 0.03;
    const auto input = data.input();
    const auto result = evaluate_ring_gps(input, {true, true, true, false});
    const double area = 2.0 * std::acos(-1.0) * 1.2 * 2.02;
    check(near(result.normal_force, 100.0 * area) && near(result.heat_rate, 4000.0 * area),
        "Finite contact uses the current secondary radius and axial length for both sides");
    check(near(result.points[0].weighted_measure + result.points[1].weighted_measure, area),
        "Finite interface quadrature sums to the current common cylindrical area");
    for (std::size_t column = 0; column < 8; ++column)
        check(result.jacobian[column] == -result.jacobian[8 + column]
                  && result.jacobian[16 + column] == -result.jacobian[24 + column]
                  && result.jacobian[32 + column] == -result.jacobian[48 + column]
                  && result.jacobian[40 + column] == -result.jacobian[56 + column],
            "Finite thermal, radial and axial geometry tangents conserve both sides exactly");
    check(result.residual[0] == -result.residual[1] && result.residual[2] == -result.residual[3]
              && result.residual[4] == -result.residual[6] && result.residual[5] == -result.residual[7],
        "Finite contact conserves heat and all forces with unequal current axial lengths");
    check(result.jacobian[6] > 0.0 && result.jacobian[7] < 0.0 && result.jacobian[3] < 0.0,
        "Thermal residual retains secondary radial and both axial geometric derivatives");
    tangent(data, true, 0.0);
    for (auto* state : {&data._state, &data._committed_state})
        for (const std::size_t upper : {std::size_t(5), std::size_t(7)}) {
            const double saved = (*state)[upper];
            const double saved_lower = (*state)[upper - 1];
            (*state)[upper - 1] = 0.0;
            (*state)[upper] = -2.0;
            rejected([&] { evaluate_ring_gps(input); },
                "Finite contact rejects zero current and committed axial lengths on either side");
            (*state)[upper] = saved;
            (*state)[upper - 1] = saved_lower;
        }
}
} // namespace

int main() {
    conservation_and_requests();
    jacobians();
    history_transactions(StrainFormulation::small);
    history_transactions(StrainFormulation::finite);
    antisymmetric_slip();
    invalid_inputs();
    finite_geometry();
    return failures == 0 ? 0 : 1;
}
