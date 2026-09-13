#include "c3d20rt.hpp"
#include "c3d20t.hpp"
#include "cax2t_gps.hpp"
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "cax8rt.hpp"
#include "cax8t.hpp"
#include "support/c3d20_common_tests.hpp"
#include "support/c3d8_common_tests.hpp"
#include <limits>

namespace {
bool check(bool condition, const char* message) {
    if (!condition)
        std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

bool check_load(const double* loaded,
    const double* unloaded,
    const double* tangent,
    const double* unloaded_tangent,
    std::size_t temperatures,
    std::size_t nodes,
    std::size_t components,
    double mass,
    const std::array<double, 3>& acceleration) {
    const std::size_t size = temperatures + nodes * components;
    bool passed = true;
    for (std::size_t component = 0; component < components; ++component) {
        double force = 0.0;
        for (std::size_t node = 0; node < nodes; ++node) {
            const auto index = temperatures + component * nodes + node;
            force += loaded[index] - unloaded[index];
        }
        const auto direction = components == 2 && component == 1 ? 2 : component;
        passed = check(std::abs(force + mass * acceleration[direction]) < 1e-8 * mass,
                     "body force equals initial mass times global acceleration")
                 && passed;
    }
    for (std::size_t row = 0; row < temperatures; ++row)
        passed = check(loaded[row] == unloaded[row], "gravity does not add heat") && passed;
    for (std::size_t entry = 0; entry < size * size; ++entry)
        passed = check(tangent[entry] == unloaded_tangent[entry],
                     "current density and volume derivatives cancel in the body force")
                 && passed;
    return passed;
}
} // namespace

int main() {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    const IsotropicThermoelasticMaterial material(test::thermoelastic(0, 4, 1e5, .25, 0, 300, 0, 0, 0, 2000, 1000));
    bool passed = true;
    const auto ratio = adlite::Scalar::independent(1.5, 0, 1);
    const auto density = material.current_density(300, ratio);
    const auto mass_density = density * ratio;
    passed = check(std::abs(density.value() - 2000.0 / 1.5) < 1e-12
                       && std::abs(density.derivative(0) + 2000.0 / 2.25) < 1e-12
                       && std::abs(mass_density.value() - 2000) < 1e-12 && std::abs(mass_density.derivative(0)) < 1e-12,
                 "rho_initial/J conserves mass with its automatic derivative")
             && passed;
    for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity()}) {
        bool rejected = false;
        try {
            (void)material.current_density(300, adlite::Scalar(invalid));
        } catch (const std::domain_error&) {
            rejected = true;
        }
        passed = check(rejected, "invalid density volume ratio is rejected") && passed;
    }
    const double pi = std::acos(-1.0);
    for (const auto formulation : {StrainFormulation::small, StrainFormulation::finite}) {
        const std::array<double, 3> axis_acceleration{2, 0, -9.81};
        const Quad4Coordinates corners{{{1, 0}, {2, 0}, {2, 1}, {1, 1}}};
        const auto geometry = make_cax4t_geometry(corners);
        Cax4LocalValues state{}, old{};
        for (std::size_t node = 0; node < 4; ++node) {
            old[node] = 300;
            state[node] = 400 + 10 * static_cast<double>(node);
            state[4 + node] = .2 * corners[node].r;
            state[8 + node] = .3 * corners[node].z;
        }
        for (bool reduced : {false, true}) {
            Cax4Input input{material, geometry, state, old};
            input.strain_formulation = formulation;
            input.initial_temperature = 300;
            const auto unloaded = reduced ? evaluate_cax4rt(input) : evaluate_cax4t(input);
            input.body_acceleration = axis_acceleration;
            const auto loaded = reduced ? evaluate_cax4rt(input) : evaluate_cax4t(input);
            passed = check_load(loaded.residual.data(),
                         unloaded.residual.data(),
                         loaded.jacobian.data(),
                         unloaded.jacobian.data(),
                         4,
                         4,
                         2,
                         6000 * pi,
                         axis_acceleration)
                     && passed;
        }
        const Quad8RzCoordinates quadratic{{{1, 0}, {2, 0}, {2, 1}, {1, 1}, {1.5, 0}, {2, .5}, {1.5, 1}, {1, .5}}};
        Quad8RzValues state8{}, old8{};
        for (std::size_t node = 0; node < 4; ++node) {
            old8[node] = 300;
            state8[node] = 400 + 10 * static_cast<double>(node);
        }
        for (std::size_t node = 0; node < 8; ++node) {
            state8[4 + node] = .2 * quadratic[node].r;
            state8[12 + node] = .3 * quadratic[node].z;
        }
        for (bool reduced : {false, true}) {
            const auto geometry8 = reduced ? make_cax8rt_geometry(quadratic) : make_cax8t_geometry(quadratic);
            Cax8Input input{material, geometry8, state8, old8};
            input.initial_temperature = 300;
            input.strain_formulation = formulation;
            const auto unloaded = reduced ? evaluate_cax8rt(input) : evaluate_cax8t(input);
            input.body_acceleration = axis_acceleration;
            const auto loaded = reduced ? evaluate_cax8rt(input) : evaluate_cax8t(input);
            passed = check_load(loaded.residual.data(),
                         unloaded.residual.data(),
                         loaded.jacobian.data(),
                         unloaded.jacobian.data(),
                         4,
                         8,
                         2,
                         6000 * pi,
                         axis_acceleration)
                     && passed;
        }
        const auto bar = make_cax2t_gps_geometry({1, 2}, 0, 1);
        const Cax2tGpsLocalValues bar_state{400, 420, .2, .4, 0, .3}, bar_old{300, 300, 0, 0, 0, 0};
        Cax2tGpsInput bar_input{material, bar, bar_state, bar_old};
        bar_input.initial_temperature = 300;
        bar_input.strain_formulation = formulation;
        const auto bar_unloaded = evaluate_cax2t_gps(bar_input);
        bar_input.body_acceleration = axis_acceleration;
        const auto bar_loaded = evaluate_cax2t_gps(bar_input);
        passed = check_load(bar_loaded.residual.data(),
                     bar_unloaded.residual.data(),
                     bar_loaded.jacobian.data(),
                     bar_unloaded.jacobian.data(),
                     2,
                     2,
                     2,
                     6000 * pi,
                     axis_acceleration)
                 && passed;

        const std::array<double, 3> acceleration{2, -3, -9.81};
        const auto cube = test::c3d8::unit_cube();
        Hex8LocalValues cartesian{}, cartesian_old{};
        for (std::size_t node = 0; node < 8; ++node) {
            cartesian[node] = 400 + 10 * static_cast<double>(node);
            cartesian_old[node] = 300;
            cartesian[8 + node] = .2 * cube[node].x;
            cartesian[16 + node] = .3 * cube[node].y;
            cartesian[24 + node] = .4 * cube[node].z;
        }
        for (bool reduced : {false, true}) {
            const auto volume = reduced ? make_c3d8rt_geometry(cube) : make_c3d8t_geometry(cube);
            C3d8Input input{material, volume, cartesian, cartesian_old};
            input.initial_temperature = 300;
            input.strain_formulation = formulation;
            const auto unloaded = reduced ? evaluate_c3d8rt(input) : evaluate_c3d8t(input);
            input.body_acceleration = acceleration;
            const auto loaded = reduced ? evaluate_c3d8rt(input) : evaluate_c3d8t(input);
            passed = check_load(loaded.residual.data(),
                         unloaded.residual.data(),
                         loaded.jacobian.data(),
                         unloaded.jacobian.data(),
                         8,
                         8,
                         3,
                         2000,
                         acceleration)
                     && passed;
        }
        const auto cube20 = test::c3d20::unit_cube();
        Hex20LocalValues state20{}, old20{};
        for (std::size_t node = 0; node < 8; ++node) {
            state20[node] = 400 + 10 * static_cast<double>(node);
            old20[node] = 300;
        }
        for (std::size_t node = 0; node < 20; ++node) {
            state20[8 + node] = .2 * cube20[node].x;
            state20[28 + node] = .3 * cube20[node].y;
            state20[48 + node] = .4 * cube20[node].z;
        }
        for (bool reduced : {false, true}) {
            const auto volume = reduced ? make_c3d20rt_geometry(cube20) : make_c3d20t_geometry(cube20);
            C3d20Input input{material, volume, state20, old20};
            input.initial_temperature = 300;
            input.strain_formulation = formulation;
            const auto unloaded = reduced ? evaluate_c3d20rt(input) : evaluate_c3d20t(input);
            input.body_acceleration = acceleration;
            const auto loaded = reduced ? evaluate_c3d20rt(input) : evaluate_c3d20t(input);
            passed = check_load(loaded.residual.data(),
                         unloaded.residual.data(),
                         loaded.jacobian.data(),
                         unloaded.jacobian.data(),
                         8,
                         20,
                         3,
                         2000,
                         acceleration)
                     && passed;
        }
    }
    if (passed)
        std::cout << "[PASS] Initial mass and gravity for all nine element formulations\n";
    return passed ? 0 : 1;
}
