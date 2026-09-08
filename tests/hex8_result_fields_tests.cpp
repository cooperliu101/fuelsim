#include "io/hex8_result_fields.hpp"
#include "support/material_factory.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void close(double actual, double expected, const char* field) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > 2e-12 * std::max(1.0, std::abs(expected)))
        throw std::runtime_error(std::string("HEX8 derived output disagrees with analytic ") + field);
}

void affine_case(double angle, bool reduced, bool finite) {
    const fuelsim::Hex8Coordinates nodes = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    const auto geometry = fuelsim::make_hex8_geometry(nodes);
    const fuelsim::IsotropicThermoelasticMaterial material(fuelsim::test::thermoelastic(0, 4, 1e9, 0.3, 0, 300));
    const double c = std::cos(angle), s = std::sin(angle);
    constexpr double a = 1.1, b = 0.9, d = 1.05;
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        const auto& p = nodes[node];
        state[node] = 300 + 2 * p.x + 3 * p.y + 5 * p.z;
        state[8 + node] = (c * a - 1) * p.x - s * b * p.y;
        state[16 + node] = s * a * p.x + (c * b - 1) * p.y;
        state[24 + node] = (d - 1) * p.z;
    }
    const auto values = fuelsim::io_detail::hex8_derived_results(geometry,
        state,
        material,
        finite ? fuelsim::StrainFormulation::finite : fuelsim::StrainFormulation::small,
        reduced,
        1.0);
    const std::array<double, 3> gradient =
        finite ? std::array<double, 3>{c * 2 / a - s * 3 / b, s * 2 / a + c * 3 / b, 5 / d}
               : std::array<double, 3>{2, 3, 5};
    const std::array<double, 6> logarithmic = {c * c * std::log(a) + s * s * std::log(b),
        s * s * std::log(a) + c * c * std::log(b),
        std::log(d),
        c * s * (std::log(a) - std::log(b)),
        0,
        0};
    const std::array<double, 6> infinitesimal = {c * a - 1, c * b - 1, d - 1, 0.5 * s * (a - b), 0, 0};
    constexpr std::array<std::size_t, 8> material_node = {0, 1, 3, 2, 4, 5, 7, 6};
    for (std::size_t q = 0; q < 8; ++q) {
        const auto& value = values[q];
        close(value[3], c * a * value[0] - s * b * value[1], "current x");
        close(value[4], s * a * value[0] + c * b * value[1], "current y");
        close(value[5], d * value[2], "current z");
        close(value[6], reduced ? 305 : state[material_node[q]], "material temperature");
        close(value[7], (finite ? a * b * d : 1.0) * (reduced ? 1.0 : 0.125), "integration measure");
        close(value[23], value[7], "local current measure");
        for (std::size_t component = 0; component < 3; ++component)
            close(value[8 + component], -4 * gradient[component], "heat flux");
        for (std::size_t component = 0; component < 6; ++component) {
            close(value[11 + component], logarithmic[component], "logarithmic strain");
            close(value[17 + component], infinitesimal[component], "infinitesimal strain");
        }
    }
}

void singular_small_strain() {
    const fuelsim::Hex8Coordinates nodes = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    const auto geometry = fuelsim::make_hex8_geometry(nodes);
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = 300;
        state[8 + node] = -nodes[node].x;
    }
    const fuelsim::IsotropicThermoelasticMaterial material(fuelsim::test::thermoelastic(0, 4, 1e9, 0.3, 0, 300));
    const auto values =
        fuelsim::io_detail::hex8_derived_results(geometry, state, material, fuelsim::StrainFormulation::small, true, 0);
    close(values[0][17], -1, "linearized strain at singular total deformation");
    if (!std::isnan(values[0][11]))
        throw std::runtime_error("Undefined derived logarithm must be marked missing");
    bool rejected = false;
    try {
        fuelsim::io_detail::hex8_derived_results(geometry,
            state,
            material,
            fuelsim::StrainFormulation::finite,
            true,
            0);
    } catch (const std::domain_error&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("Finite-strain output must reject nonpositive geometry");
}
} // namespace

int main() {
    try {
        for (bool reduced : {false, true})
            for (bool finite : {false, true})
                for (double angle : {0.0, 0.37, 3.14159265358979323846})
                    affine_case(angle, reduced, finite);
        singular_small_strain();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
