#include "support/c3d20_common_tests.hpp"

namespace {
using namespace fuelsim::test::c3d20;

bool test_reduced_volume_kernel() {
    const auto coordinates = unit_cube();
    const auto geometry = fuelsim::test::make_c3d20_geometry(coordinates, fuelsim::Hex20ElementFormulation::c3d20rt);
    bool passed = check(geometry.mechanical_points.size() == 8, "C3D20RT has eight active material points");
    double volume = 0.0;
    for (const auto& point : geometry.mechanical_points)
        volume += point.weighted_measure;
    passed = check(near(volume, 1.0, 1.0e-14), "C3D20RT integrates the unit cube volume") && passed;
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        for (const unsigned mechanism : {0U, 1U, 2U, 3U}) {
            auto properties = material();
            if ((mechanism & 1U) != 0)
                properties = fuelsim::test::with_norton(std::move(properties), 1.0e-6, 10.0, 3.0, 300.0);
            if ((mechanism & 2U) != 0)
                properties = fuelsim::test::with_plasticity(std::move(properties), 20.0, 10.0, 300.0);
            const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(properties),
                4.0,
                1.0,
                formulation,
                fuelsim::Hex8ElementFormulation::c3d8t,
                300.0};
            fuelsim::Hex20LocalValues old{}, state{};
            for (std::size_t node = 0; node < 8; ++node) {
                old[node] = 300.0;
                state[node] = 302.0 + coordinates[node].x + 2.0 * coordinates[node].y;
            }
            for (std::size_t node = 0; node < 20; ++node) {
                state[8 + node] = 0.002 * coordinates[node].x + 0.0003 * coordinates[node].y;
                state[28 + node] = -0.0004 * coordinates[node].y;
                state[48 + node] = 0.0002 * coordinates[node].z;
            }
            const fuelsim::CartesianMaterialHistory history(8);
            passed = check(directional_jacobian_error(data, geometry, state, old, history) < 2.0e-6,
                         "C3D20RT elastic/creep/plastic/coupled small/finite tangent matches centered differences")
                     && passed;
            passed = check(residual_path_error(data, geometry, state, old, history) < 1.0e-13,
                         "C3D20RT independent residual matches the Jacobian residual")
                     && passed;
            const auto updated = fuelsim::compute_c3d20_transient_update(data, geometry, state, old, history, 0.5);
            passed = check(updated.size() == 8, "C3D20RT updates exactly eight material histories") && passed;
            if (mechanism == 0)
                passed = check(fuelsim::compute_c3d20_stress(data, geometry, state).size() == 8,
                             "C3D20RT stress recovery returns only active points")
                         && passed;
        }
    }
    return passed;
}

int run_c3d20rt_tests() {
    return test_reduced_volume_kernel() && test_reference_mass_capacity(fuelsim::Hex20ElementFormulation::c3d20rt) ? 0
                                                                                                                   : 1;
}
} // namespace

int main() {
    return run_c3d20rt_tests();
}
