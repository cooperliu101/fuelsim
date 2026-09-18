#include "cax2t_gps.hpp"
#include "cpeg8t.hpp"
#include "support/c3d20_common_tests.hpp"
#include "support/c3d8_common_tests.hpp"
#include <limits>

namespace {
using namespace fuelsim;
using namespace fuelsim::elements;

void require(bool value, const char* message) {
    if (!value)
        throw std::runtime_error(message);
}

void near(double actual, double expected) {
    require(std::abs(actual - expected) < 2e-12, "Material point rate sampling differs from independent expectation");
}

double expected(double temperature, double x, double y, double z, double history) {
    return 1e-4 * temperature + 1e-3 * (2.0 + x + 2 * y + 3 * z) + history + 0.01 * std::sqrt(27.0);
}

IsotropicThermoelasticMaterial material() {
    auto base = test::thermoelastic(0, 1, 1e9, 0.3, 0, 300);
    auto functions = std::make_shared<MaterialFunctionSet>(*base.functions);
    functions->creep.builtin.kind = CreepBuiltinParameters::Kind::custom;
    functions->creep.function = [](const CreepRateInput& input) {
        return 1e-4 * input.temperature
               + 1e-3 * (input.context.time + input.context.x + 2 * input.context.y + 3 * input.context.z)
               + input.equivalent_creep_strain + 0.01 * input.equivalent_stress;
    };
    return IsotropicThermoelasticMaterial({functions, 300});
}

void axisymmetric() {
    const auto m = material();
    MaterialPointState point;
    point.stress = {7, 7, 7, 3}; // Pure shear plus hydrostatic stress: q=sqrt(3)*3.
    point.equivalent_creep_strain = 0.2;
    const auto bar = make_cax2t_gps_geometry({1, 2}, 0, 1);
    const auto br = cax2t_gps_creep_rates(m, bar, {300, 400, 0, 0, 0, 0}, {point, point}, 2);
    require(br.size() == 2, "BAR2 active point count");
    for (std::size_t q = 0; q < 2; ++q)
        near(br[q],
            expected(300 + 100 * static_cast<double>(q), 1.5 + (q == 0 ? -1 : 1) / std::sqrt(12.0), 0, 0.5, 0.2));
    const Quad4Coordinates corners{{{1, 0}, {2, 0}, {2, 1}, {1, 1}}};
    const auto geometry = make_cax4t_geometry(corners);
    Cax4LocalValues state{300, 400, 450, 350};
    Quad4MaterialHistory history{point, point, point, point};
    for (std::size_t q = 0; q < history.size(); ++q)
        history[q].equivalent_creep_strain += 0.01 * static_cast<double>(q);
    const auto full = cax4t_creep_rates(m, geometry, state, history, 2);
    require(full.size() == 4, "CAX4T active point count");
    for (std::size_t q = 0; q < 4; ++q)
        near(full[q],
            expected(state[q],
                geometry.points[q].radius,
                0,
                geometry.points[q].axial_coordinate,
                0.2 + 0.01 * static_cast<double>(q)));
    // Inactive slots deliberately contain invalid values and must never be read.
    for (std::size_t q = 1; q < 4; ++q)
        history[q].equivalent_creep_strain = std::numeric_limits<double>::quiet_NaN();
    for (auto formulation : {StrainFormulation::small, StrainFormulation::finite}) {
        auto current = corners;
        for (std::size_t n = 0; n < 4; ++n) {
            state[4 + n] = 0.2 * (corners[n].r - 1) * corners[n].z;
            if (formulation == StrainFormulation::finite)
                current[n].r += state[4 + n];
        }
        const auto displaced = make_cax4t_geometry(current);
        double temperature = 0, volume = 0;
        for (const auto& p : displaced.points) {
            volume += p.weighted_measure;
            for (std::size_t n = 0; n < 4; ++n)
                temperature += p.weighted_measure * p.shape[n] * state[n];
        }
        const auto rates = cax4rt_creep_rates(m, geometry, state, history, 2, formulation);
        require(rates.size() == 1, "CAX4RT skips inactive slots");
        near(rates[0], expected(temperature / volume, 14.0 / 9, 0, 0.5, 0.2));
    }
    const Quad8RzCoordinates nodes{{{1, 0}, {2, 0}, {2, 1}, {1, 1}, {1.5, 0}, {2, 0.5}, {1.5, 1}, {1, 0.5}}};
    Quad8RzValues quadratic{300, 400, 450, 350};
    for (bool reduced : {false, true}) {
        const auto g = reduced ? make_cax8rt_geometry(nodes) : make_cax8t_geometry(nodes);
        Quad8MaterialHistory h;
        h.fill(point);
        for (std::size_t q = 0; q < h.size(); ++q)
            h[q].equivalent_creep_strain += 0.01 * static_cast<double>(q);
        if (reduced)
            for (std::size_t q = 4; q < 9; ++q)
                h[q].equivalent_creep_strain = std::numeric_limits<double>::quiet_NaN();
        const auto rates =
            reduced ? cax8rt_creep_rates(m, g, quadratic, h, 2) : cax8t_creep_rates(m, g, quadratic, h, 2);
        require(rates.size() == (reduced ? 4 : 9), "CAX8 active material point count");
        for (std::size_t q = 0; q < rates.size(); ++q) {
            const auto& p = g.points[q];
            near(rates[q],
                expected(300 + 100 * (p.radius - 1) + 50 * p.axial_coordinate,
                    p.radius,
                    0,
                    p.axial_coordinate,
                    0.2 + 0.01 * static_cast<double>(q)));
        }
    }
}

void cartesian() {
    const auto m = material();
    CartesianMaterialPointState point;
    point.stress = {7, 7, 7, 3, 0, 0};
    point.equivalent_creep_strain = 0.2;
    const auto nodes = test::c3d8::unit_cube();
    const auto g = make_c3d8rt_geometry(nodes);
    Hex8LocalValues state{};
    for (std::size_t n = 0; n < 8; ++n)
        state[n] = 300 + 100 * nodes[n].x + 50 * nodes[n].y + 20 * nodes[n].z;
    for (auto formulation : {StrainFormulation::small, StrainFormulation::finite}) {
        auto current = nodes;
        for (std::size_t n = 0; n < 8; ++n) {
            state[8 + n] = 0.2 * nodes[n].x * nodes[n].y;
            if (formulation == StrainFormulation::finite)
                current[n].x += state[8 + n];
        }
        const auto displaced = make_c3d8t_geometry(current);
        double temperature = 0, volume = 0;
        for (const auto& p : displaced.points) {
            volume += p.weighted_measure;
            for (std::size_t n = 0; n < 8; ++n)
                temperature += p.weighted_measure * p.shape[n] * state[n];
        }
        const auto rates = c3d8rt_creep_rates(m, g, state, {point}, 2, formulation);
        require(rates.size() == 1, "C3D8RT active material point count");
        near(rates[0], expected(temperature / volume, 0.5, 0.5, 0.5, 0.2));
    }
    const auto quadratic_nodes = test::c3d20::unit_cube();
    Hex20LocalValues quadratic{};
    for (std::size_t n = 0; n < 8; ++n)
        quadratic[n] = 300 + 100 * nodes[n].x + 50 * nodes[n].y + 20 * nodes[n].z;
    for (bool reduced : {false, true}) {
        const auto geometry = reduced ? make_c3d20rt_geometry(quadratic_nodes) : make_c3d20t_geometry(quadratic_nodes);
        CartesianMaterialHistory h(reduced ? 8 : 27, point);
        for (std::size_t q = 0; q < h.size(); ++q)
            h[q].equivalent_creep_strain += 0.01 * static_cast<double>(q);
        const auto rates = reduced ? c3d20rt_creep_rates(m, geometry, quadratic, h, 2)
                                   : c3d20t_creep_rates(m, geometry, quadratic, h, 2);
        require(rates.size() == h.size(), "C3D20 active material point count");
        for (std::size_t q = 0; q < rates.size(); ++q) {
            const auto p = geometry.mechanical_points[q].position;
            near(rates[q],
                expected(300 + 100 * p.x + 50 * p.y + 20 * p.z, p.x, p.y, p.z, 0.2 + 0.01 * static_cast<double>(q)));
        }
    }
    const Cpeg8Coordinates plane{{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0}, {1, 0.5}, {0.5, 1}, {0, 0.5}}};
    const auto pg = make_cpeg8t_geometry(plane, 1);
    const Cpeg8Values ps{300, 400, 450, 350};
    const auto pr = cpeg8t_creep_rates(m, pg, ps, CartesianMaterialHistory(9, point), 2);
    require(pr.size() == 9, "CPEG8T active material point count");
    for (std::size_t q = 0; q < 9; ++q) {
        const auto& p = pg.points[q];
        near(pr[q], expected(300 + 100 * p.x + 50 * p.y, p.x, p.y, 0, 0.2));
    }
}
} // namespace

int main() {
    try {
        axisymmetric();
        cartesian();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
