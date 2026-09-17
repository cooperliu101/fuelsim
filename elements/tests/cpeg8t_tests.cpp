#include "cpeg8t.hpp"
#include "line3_plane.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require_close(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance)
        throw std::runtime_error(message);
}

void check_contact() {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    Line3PlaneValues state{410, 430, 300, 320};
    state[16] = .01;
    state[17] = .03;
    state[18] = -.02;
    Line3PlaneContactInput input{{{{-.01, -.0001}, {.01, -.0001}, {0, -.00008}}},
        {{{.012, 0}, {-.012, 0}, {0, .00001}}},
        {},
        {},
        .1,
        .1,
        true,
        true,
        state,
        .2,
        8.0 / 9.0,
        {0, 1e-6, GapHeatConductanceLaw::affine, 1000, 10, .001, .2, 300, 1e9},
        1e9,
        true,
        true};
    const auto active = evaluate_line3_plane_contact(input, true);
    if (!active.projected || !(active.pressure > 0))
        throw std::runtime_error("Curved contact must be active");
    const auto passive = evaluate_line3_plane_contact(input, false);
    Line3PlaneValues direction{}, plus = state, minus = state;
    constexpr double epsilon = 1e-5;
    for (std::size_t i = 0; i < 22; ++i) {
        direction[i] = std::sin(static_cast<double>(i + 1)) * (i < 4 ? 10 : .0001);
        plus[i] += epsilon * direction[i];
        minus[i] -= epsilon * direction[i];
        require_close(active.residual[i], passive.residual[i], 1e-12, "Contact residual agreement");
    }
    auto upper = input, lower = input;
    // Reference members intentionally bind separate perturbed local states.
    const Line3PlaneContactInput upper_input{upper.secondary,
        upper.primary,
        {},
        {},
        .1,
        .1,
        true,
        true,
        plus,
        upper.coordinate,
        upper.weight,
        upper.heat,
        upper.penalty,
        true,
        true};
    const Line3PlaneContactInput lower_input{lower.secondary,
        lower.primary,
        {},
        {},
        .1,
        .1,
        true,
        true,
        minus,
        lower.coordinate,
        lower.weight,
        lower.heat,
        lower.penalty,
        true,
        true};
    const auto rp = evaluate_line3_plane_contact(upper_input, false);
    const auto rm = evaluate_line3_plane_contact(lower_input, false);
    for (std::size_t i = 0; i < 22; ++i) {
        double exact = 0;
        for (std::size_t j = 0; j < 22; ++j)
            exact += active.jacobian[i * 22 + j] * direction[j];
        require_close(exact,
            (rp.residual[i] - rm.residual[i]) / (2 * epsilon),
            2e-6 * std::max(1.0, std::abs(exact)),
            "Curved contact directional tangent");
    }
    for (std::size_t column = 0; column <= 22; ++column) {
        const auto entry = [&](std::size_t row) {
            return column == 22 ? active.residual[row] : active.jacobian[row * 22 + column];
        };
        require_close(entry(0) + entry(1) + entry(2) + entry(3), 0, 1e-9, "Contact heat conservation");
        for (std::size_t c = 0; c < 2; ++c) {
            double force = 0;
            for (std::size_t n = 0; n < 3; ++n)
                force += entry(4 + 3 * c + n) + entry(10 + 3 * c + n);
            require_close(force, 0, 1e-9, "Contact force conservation");
        }
    }
    for (std::size_t n = 0; n < 3; ++n)
        state[7 + n] = .002;
    const auto open = evaluate_line3_plane_contact(input, true);
    require_close(open.pressure, 0, 0, "Open contact pressure");
    for (std::size_t n = 0; n < 3; ++n)
        state[4 + n] = .1;
    if (evaluate_line3_plane_contact(input, false).projected)
        throw std::runtime_error("Projection outside the primary edge must be invalid");
}
} // namespace

int main() {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    try {
        check_contact();
        const Cpeg8Coordinates coordinates{{{-0.01, -0.005},
            {0.01, -0.005},
            {0.01, 0.005},
            {-0.01, 0.005},
            {0.0, -0.005},
            {0.01, 0.0},
            {0.0, 0.005},
            {-0.01, 0.0}}};
        const auto geometry = make_cpeg8t_geometry(coordinates, 0.1);
        const IsotropicThermoelasticMaterial material(
            test::thermoelastic(0, 10, 1e6, .25, 1e-5, 300, 0, 0, 0, 1000, 100));
        Cpeg8Values state{}, previous{};
        for (std::size_t n = 0; n < 4; ++n) {
            state[n] = 400.0;
            previous[n] = 300.0;
        }
        state[20] = 0.0001;
        state[21] = 0.002;
        state[22] = -0.003;
        const Cpeg8History history{};
        Cpeg8Input input{material, geometry, state, previous, &history, 1.0, 1.0};
        input.include_thermal_time_term = true;
        const auto result = evaluate_cpeg8t(input, {true, true, true, true});
        require_close(result.residual[20], -0.16, 1e-13, "Uniform extension reaction");
        require_close(result.residual[21], 4e-5, 1e-15, "First bending reaction");
        require_close(result.residual[22], -2.4e-4, 1e-15, "Second bending reaction");
        require_close(result.stored_heat_rate, 200.0, 1e-10, "Initial mass heat storage");
        for (std::size_t q = 0; q < 9; ++q) {
            const auto& p = geometry.points[q];
            const double ezz = 0.001 + 0.02 * p.y + 0.03 * p.x;
            require_close(result.stress[q].zz, 1.2e6 * ezz - 2000.0, 1e-9, "Analytical axial stress");
            require_close(result.stress[q].xx, 0.4e6 * ezz - 2000.0, 1e-9, "Analytical transverse stress");
        }
        state[20] = 0.01;
        state[21] = 0.2;
        state[22] = -0.3;
        input.strain_formulation = StrainFormulation::finite;
        const auto finite_probe = evaluate_cpeg8t(input, {true, true, true, true});
        // Independently extracted native one-increment Abaqus 2025 values.
        const std::array<double, 9> native_axial_stress{78056.7560071207,
            103823.51684521978,
            129024.20091959235,
            86709.74668762767,
            112285.71428571422,
            137301.84656868933,
            95298.42888339,
            120685.71457289465,
            145519.3149743173};
        for (std::size_t q = 0; q < 9; ++q)
            require_close(finite_probe.stress[q].zz, native_axial_stress[q], 1e-7, "Finite axial increment probe");
        state[20] = 0.0001;
        state[21] = 0.002;
        state[22] = -0.003;
        // Curvature, temperature gradients and non-affine displacements exercise
        // the complete point-to-node chain, including shared control columns.
        state[1] += 13.0;
        state[2] -= 7.0;
        state[3] += 3.0;
        state[8] = 0.0002;
        state[17] = -0.0001;
        input.volumetric_heat_source = 35.0;
        Cpeg8Values direction{};
        for (std::size_t i = 0; i < 23; ++i)
            direction[i] = std::sin(static_cast<double>(i + 1)) * (i < 4 ? 10.0 : 0.001);
        const auto base = test::thermoelastic(0, 10, 1e6, .25, 1e-5, 300, 0, 0, 0, 1000, 100);
        const std::array<ThermoelasticProperties, 4> materials{base,
            test::with_plasticity(base, 500, 10000),
            test::with_norton(base, 1e-5, 1000, 3),
            test::with_norton(test::with_plasticity(base, 500, 10000), 1e-5, 1000, 3)};
        for (const auto& properties : materials)
            for (auto formulation : {StrainFormulation::small, StrainFormulation::finite}) {
                const IsotropicThermoelasticMaterial tested_material(properties);
                const Cpeg8Input
                    tested{tested_material, geometry, state, previous, &history, 1.0, 1.0, 35.0, formulation, true};
                const auto active = evaluate_cpeg8t(tested, {true, true, true, true});
                const auto passive = evaluate_cpeg8t(tested, {true, false, true, true});
                Cpeg8Values plus = state, minus = state;
                constexpr double epsilon = 1e-5;
                for (std::size_t i = 0; i < 23; ++i) {
                    plus[i] += epsilon * direction[i];
                    minus[i] -= epsilon * direction[i];
                    require_close(active.residual[i],
                        passive.residual[i],
                        1e-12 * std::max(1.0, std::abs(passive.residual[i])),
                        "Residual path agreement");
                }
                const Cpeg8Input
                    upper{tested_material, geometry, plus, previous, &history, 1.0, 1.0, 35.0, formulation, true};
                const Cpeg8Input
                    lower{tested_material, geometry, minus, previous, &history, 1.0, 1.0, 35.0, formulation, true};
                const auto upper_result = evaluate_cpeg8t(upper), lower_result = evaluate_cpeg8t(lower);
                for (std::size_t i = 0; i < 23; ++i) {
                    double exact = 0.0;
                    for (std::size_t j = 0; j < 23; ++j)
                        exact += active.jacobian[23 * i + j] * direction[j];
                    const double difference = (upper_result.residual[i] - lower_result.residual[i]) / (2.0 * epsilon);
                    require_close(exact, difference, 1e-6 * std::max(1.0, std::abs(exact)), "Directional tangent");
                }
                bool plastic_active = false, creep_active = false;
                for (const auto& point : active.history) {
                    plastic_active = plastic_active || point.equivalent_plastic_strain > 0;
                    creep_active = creep_active || point.equivalent_creep_strain > 0;
                    require_close(point.plastic_strain[0] + point.plastic_strain[1] + point.plastic_strain[2],
                        0,
                        1e-14,
                        "Plastic incompressibility");
                    require_close(point.creep_strain[0] + point.creep_strain[1] + point.creep_strain[2],
                        0,
                        1e-14,
                        "Creep incompressibility");
                }
                if (properties.functions->has_plasticity() != plastic_active
                    || properties.functions->has_creep() != creep_active)
                    throw std::runtime_error("Requested inelastic mechanism was not exercised");
                if (formulation == StrainFormulation::finite) {
                    const double cosine = std::cos(.3), sine = std::sin(.3);
                    auto rotated = state;
                    for (std::size_t n = 0; n < 8; ++n) {
                        const double x = coordinates[n][0] + state[4 + n], y = coordinates[n][1] + state[12 + n];
                        rotated[4 + n] = cosine * x - sine * y - coordinates[n][0];
                        rotated[12 + n] = sine * x + cosine * y - coordinates[n][1];
                    }
                    rotated[21] = cosine * state[21] - sine * state[22];
                    rotated[22] = cosine * state[22] + sine * state[21];
                    const Cpeg8Input rotation_input{tested_material,
                        geometry,
                        rotated,
                        state,
                        &active.history,
                        0,
                        1,
                        0,
                        StrainFormulation::finite,
                        false};
                    const auto objective = evaluate_cpeg8t(rotation_input, {true, false, true, true});
                    require_close(objective.elastic_energy_change, 0, 1e-12, "Rigid rotation preserves elastic energy");
                    require_close(objective.plastic_dissipation_increment,
                        0,
                        1e-12,
                        "Rigid rotation has no plastic work");
                    require_close(objective.creep_dissipation_increment, 0, 1e-12, "Rigid rotation has no creep work");
                    for (std::size_t q = 0; q < 9; ++q) {
                        const auto& s = active.stress[q];
                        require_close(objective.stress[q].xx,
                            cosine * cosine * s.xx + sine * sine * s.yy - 2 * cosine * sine * s.xy,
                            1e-7,
                            "Rigid rotation transforms stress objectively");
                        require_close(objective.history[q].equivalent_plastic_strain,
                            active.history[q].equivalent_plastic_strain,
                            1e-14,
                            "Plastic scalar is rotation invariant");
                        require_close(objective.history[q].equivalent_creep_strain,
                            active.history[q].equivalent_creep_strain,
                            1e-14,
                            "Creep scalar is rotation invariant");
                    }
                }
            }
        for (auto kind : {Cpeg8BoundaryKind::pressure,
                 Cpeg8BoundaryKind::traction_x,
                 Cpeg8BoundaryKind::traction_y,
                 Cpeg8BoundaryKind::heat_flux,
                 Cpeg8BoundaryKind::convection})
            for (bool current : {false, true}) {
                Cpeg8Values plus = state, minus = state;
                constexpr double epsilon = 1e-5;
                for (std::size_t i = 0; i < 23; ++i) {
                    plus[i] += epsilon * direction[i];
                    minus[i] -= epsilon * direction[i];
                }
                const auto active = evaluate_cpeg8t_boundary({geometry, state, 2, kind, 1200, 350, current}, true);
                const auto upper = evaluate_cpeg8t_boundary({geometry, plus, 2, kind, 1200, 350, current}, false);
                const auto lower = evaluate_cpeg8t_boundary({geometry, minus, 2, kind, 1200, 350, current}, false);
                for (std::size_t i = 0; i < 23; ++i) {
                    double exact = 0;
                    for (std::size_t j = 0; j < 23; ++j)
                        exact += active.jacobian[23 * i + j] * direction[j];
                    require_close(exact,
                        (upper.residual[i] - lower.residual[i]) / (2 * epsilon),
                        1e-6 * std::max(1.0, std::abs(exact)),
                        "Boundary directional tangent");
                }
            }
        auto invalid = state;
        invalid[20] = -1;
        bool rejected = false;
        try {
            (void)evaluate_cpeg8t(
                {material, geometry, invalid, previous, &history, 1, 1, 0, StrainFormulation::finite});
        } catch (const std::domain_error&) {
            rejected = true;
        }
        if (!rejected)
            throw std::runtime_error("Negative thickness must reject the trial state");
        std::cout << "CPEG8T analytical, inelastic, boundary and curved contact tangent checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
