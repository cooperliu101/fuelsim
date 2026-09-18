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

void check_averaged_contact(bool curved, bool crossing, bool shared = false) {
    using namespace fuelsim::elements;
    PlaneAveragedContactGeometry geometry;
    geometry.penalty = 1e9;
    geometry.coordinates = {{{-.01, -.0001}}, {{.01, -.0001}}, {{0, -.0001}}, {{.015, 0}}, {{-.015, 0}}, {{0, 0}}};
    geometry.secondary = {{{0, 1, 2}, 0, .1}};
    geometry.primary = {{{3, 4, 5}}};
    if (curved) {
        geometry.coordinates[2][1] += .00002;
        geometry.coordinates[5][1] += .00001;
    }
    if (crossing) {
        geometry.coordinates = {{{-.005, -.0001}},
            {{.005, -.0001}},
            {{0, -.0001}},
            {{0, 0}},
            {{-.02, 0}},
            {{-.01, 0}},
            {{.02, 0}},
            {{.01, 0}}};
        geometry.primary = {{{3, 4, 5}}, {{6, 3, 7}}};
    }
    if (shared) {
        geometry.coordinates = {{{-.01, -.0002}},
            {{0, -.00025}},
            {{.01, -.00015}},
            {{-.005, -.00022}},
            {{.005, -.00021}},
            {{.015, 0}},
            {{-.015, 0}},
            {{0, .00001}}};
        geometry.secondary = {{{0, 1, 3}, 1, .1}, {{1, 2, 4}, 0, .1}};
        geometry.primary = {{{5, 6, 7}}};
    }
    const std::size_t size = 2 * geometry.coordinates.size();
    std::vector<double> state(size), direction(size), plus(size), minus(size), total(size);
    for (std::size_t i = 0; i < size; ++i) {
        direction[i] = .0001 * std::sin(static_cast<double>(i + 1));
        plus[i] = 1e-5 * direction[i];
        minus[i] = -plus[i];
    }
    for (std::size_t node = 0; node < (shared ? 1U : 3U); ++node) {
        if (!shared)
            geometry.secondary[0].local_node = node;
        const auto result = evaluate_plane_averaged_contact(geometry, state, true);
        const auto passive = evaluate_plane_averaged_contact(geometry, state, false);
        const auto upper = evaluate_plane_averaged_contact(geometry, plus, false);
        const auto lower = evaluate_plane_averaged_contact(geometry, minus, false);
        auto rotated = geometry;
        const double cosine = std::cos(.37), sine = std::sin(.37);
        for (auto& point : rotated.coordinates) {
            const auto original = point;
            point = {cosine * original[0] - sine * original[1] + .017,
                sine * original[0] + cosine * original[1] - .004};
        }
        const auto rigid = evaluate_plane_averaged_contact(rotated, state, false);
        require_close(rigid.point.gap, result.point.gap, 1e-14, "Averaged contact rigid-motion gap invariance");
        for (std::size_t n = 0; n < size / 2; ++n) {
            require_close(rigid.residual[2 * n],
                cosine * result.residual[2 * n] - sine * result.residual[2 * n + 1],
                1e-8,
                "Averaged contact objective force x");
            require_close(rigid.residual[2 * n + 1],
                sine * result.residual[2 * n] + cosine * result.residual[2 * n + 1],
                1e-8,
                "Averaged contact objective force y");
        }
        for (std::size_t i = 0; i < size; ++i) {
            total[i] += result.residual[i];
            require_close(result.residual[i], passive.residual[i], 1e-12, "Averaged contact residual agreement");
            double derivative = 0.0;
            for (std::size_t j = 0; j < size; ++j)
                derivative += result.jacobian[i * size + j] * direction[j];
            require_close(derivative,
                (upper.residual[i] - lower.residual[i]) / 2e-5,
                3e-6 * std::max(1.0, std::abs(derivative)),
                "Averaged contact full geometric tangent");
        }
        double moment = 0.0;
        for (std::size_t n = 0; n < size / 2; ++n)
            moment += geometry.coordinates[n][0] * result.residual[2 * n + 1]
                      - geometry.coordinates[n][1] * result.residual[2 * n];
        require_close(moment, 0.0, 1e-12, "Averaged contact moment conservation");
        for (std::size_t column = 0; column <= size; ++column)
            for (std::size_t c = 0; c < 2; ++c) {
                double sum = 0.0;
                for (std::size_t n = 0; n < size / 2; ++n)
                    sum += column == size ? result.residual[2 * n + c] : result.jacobian[(2 * n + c) * size + column];
                require_close(sum, 0.0, column == size ? 1e-10 : 1e-7, "Averaged contact force conservation");
            }
    }
    if (crossing) {
        // Independent native fully prescribed sliding reference at t=1 s.
        const std::array<std::size_t, 5> nodes{4, 5, 3, 7, 6};
        const std::array<double, 5> reference{-4.110544158648,
            20.608843301259,
            67.003401714779,
            20.608843301259,
            -4.110544158648};
        for (std::size_t n = 0; n < nodes.size(); ++n)
            require_close(total[2 * nodes[n] + 1], reference[n], 1e-9, "Averaged contact native primary transfer");
    }
    for (const auto& edge : geometry.secondary)
        for (auto node : edge.nodes)
            state[2 * node + 1] = .002;
    const auto open = evaluate_plane_averaged_contact(geometry, state, true);
    for (auto value : open.residual)
        require_close(value, 0.0, 0.0, "Open averaged contact has zero force");
    for (auto value : open.jacobian)
        require_close(value, 0.0, 0.0, "Open averaged contact has zero tangent");
}

void check_thermal_averaged_contact(double slide = 0.0) {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    PlaneAveragedContactGeometry geometry;
    // Both sides are curved; the primary tangent jumps at its shared endpoint.
    geometry.coordinates = {{{-.006, -.00012}},
        {{0, -.00008}},
        {{.006, -.0001}},
        {{-.003, -.00007}},
        {{.003, -.00006}},
        {{-.02, 0}},
        {{0, .00002}},
        {{.02, 0}},
        {{-.01, .00006}},
        {{.01, -.00001}}};
    geometry.secondary = {{{0, 1, 3}, 1, .1}, {{1, 2, 4}, 0, .1}};
    geometry.primary = {{{6, 5, 8}}, {{7, 6, 9}}};
    geometry.temperature_nodes = {0, 1, 2, 5, 6, 7};
    geometry.penalty = 1e9;
    geometry.heat = {0.02, 1e-5, GapHeatConductanceLaw::affine, 1000, 1e4, .001, .2, 300, 1e9};
    const auto offset = 2 * geometry.coordinates.size();
    std::vector<double> state(offset + geometry.temperature_nodes.size()), direction(state.size());
    for (std::size_t node = 0; node < 5; ++node)
        state[2 * node] = slide;
    for (std::size_t i = 0; i < direction.size(); ++i)
        direction[i] = (i < offset ? .0001 : 10.0) * std::sin(static_cast<double>(i + 1));
    for (std::size_t i = offset; i < state.size(); ++i)
        state[i] = i < offset + 3 ? 400.0 + 7.0 * static_cast<double>(i - offset)
                                  : 300.0 + 11.0 * static_cast<double>(i - offset - 3);
    for (auto law : {GapHeatConductanceLaw::affine, GapHeatConductanceLaw::gas_gap})
        for (bool open : {false, true}) {
            geometry.heat.law = law;
            for (std::size_t node = 0; node < 5; ++node)
                state[2 * node + 1] = open ? .0004 : 0.0;
            const auto exact = evaluate_plane_averaged_contact(geometry, state, true);
            const auto passive = evaluate_plane_averaged_contact(geometry, state, false);
            auto plus = state, minus = state;
            constexpr double epsilon = 1e-5;
            for (std::size_t i = 0; i < state.size(); ++i) {
                plus[i] += epsilon * direction[i];
                minus[i] -= epsilon * direction[i];
            }
            const auto upper = evaluate_plane_averaged_contact(geometry, plus, false);
            const auto lower = evaluate_plane_averaged_contact(geometry, minus, false);
            for (std::size_t i = 0; i < state.size(); ++i) {
                require_close(exact.residual[i], passive.residual[i], 0.0, "Thermal averaged residual agreement");
                double derivative = 0.0;
                for (std::size_t j = 0; j < state.size(); ++j)
                    derivative += exact.jacobian[i * state.size() + j] * direction[j];
                require_close(derivative,
                    (upper.residual[i] - lower.residual[i]) / (2 * epsilon),
                    2e-6 * std::max(1.0, std::abs(derivative)),
                    "Thermal averaged geometric and material tangent");
                if (i < offset)
                    require_close(exact.residual[i], 0.0, 0.0, "Thermal contact must not apply mechanical force");
            }
            for (std::size_t column = 0; column <= state.size(); ++column) {
                double balance = 0.0;
                for (std::size_t i = offset; i < state.size(); ++i)
                    balance += column == state.size() ? exact.residual[i] : exact.jacobian[i * state.size() + column];
                require_close(balance, 0.0, 1e-8, "Thermal averaged heat and tangent conservation");
            }
        }
    for (std::size_t node = 0; node < 5; ++node)
        state[2 * node] = 1.0;
    bool rejected = false;
    try {
        (void)evaluate_plane_averaged_contact(geometry, state, false);
    } catch (const std::domain_error&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("A thermal neighborhood with no primary projection must reject the state");
}

} // namespace

int main() {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    try {
        check_averaged_contact(false, false);
        check_averaged_contact(true, false);
        check_averaged_contact(false, true);
        check_averaged_contact(true, false, true);
        check_thermal_averaged_contact();
        check_thermal_averaged_contact(.0179);
        check_thermal_averaged_contact(-.0179);
        const Cpeg8Coordinates coordinates{{{-0.01, -0.005},
            {0.01, -0.005},
            {0.01, 0.005},
            {-0.01, 0.005},
            {0.0, -0.005},
            {0.01, 0.0},
            {0.0, 0.005},
            {-0.01, 0.0}}};
        const auto geometry = make_cpeg8t_geometry(coordinates, 0.1);
        Cpeg8Values film_state{};
        film_state[2] = 400;
        film_state[3] = 350;
        const auto film =
            evaluate_cpeg8t_boundary({geometry, film_state, 2, Cpeg8BoundaryKind::convection, 1000, 300, false}, true);
        require_close(film.residual[2], 100, 1e-12, "Native nodal film heat at first endpoint");
        require_close(film.residual[3], 50, 1e-12, "Native nodal film heat at second endpoint");
        require_close(film.jacobian[23 * 2 + 2], 1, 1e-14, "Film diagonal temperature derivative");
        require_close(film.jacobian[23 * 2 + 3], 0, 1e-14, "Film has no off-diagonal temperature derivative");
        film_state[20] = .01;
        film_state[21] = .02;
        film_state[22] = -.03;
        const auto stretched_film =
            evaluate_cpeg8t_boundary({geometry, film_state, 2, Cpeg8BoundaryKind::convection, 1000, 300, true}, true);
        require_close(stretched_film.residual[2], 100, 1e-12, "Film uses initial section thickness");
        require_close(stretched_film.residual[3], 50, 1e-12, "Film uses initial section thickness at both endpoints");
        require_close(stretched_film.jacobian[23 * 2 + 20], 0, 1e-14, "Film is independent of section extension");
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
