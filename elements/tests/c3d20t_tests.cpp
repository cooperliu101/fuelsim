#include "c3d20t.hpp"
#include "support/c3d20_common_tests.hpp"

namespace {
using namespace fuelsim::test::c3d20;

double material_state_error(const fuelsim::CartesianMaterialPointState& first,
    const fuelsim::CartesianMaterialPointState& second) {
    double error = 0.0, scale = 1.0;
    const std::array<const std::array<double, 6>*, 3> first_histories = {&first.elastic_strain,
        &first.plastic_strain,
        &first.creep_strain};
    const std::array<const std::array<double, 6>*, 3> second_histories = {&second.elastic_strain,
        &second.plastic_strain,
        &second.creep_strain};
    for (std::size_t history = 0; history < first_histories.size(); ++history)
        for (std::size_t component = 0; component < 6; ++component) {
            const double first_value = (*first_histories[history])[component];
            const double second_value = (*second_histories[history])[component];
            error = std::max(error, std::abs(first_value - second_value));
            scale = std::max({scale, std::abs(first_value), std::abs(second_value)});
        }
    const std::array<double, 8> first_scalars = {first.stress.xx,
        first.stress.yy,
        first.stress.zz,
        first.stress.xy,
        first.stress.yz,
        first.stress.xz,
        first.equivalent_plastic_strain,
        first.equivalent_creep_strain};
    const std::array<double, 8> second_scalars = {second.stress.xx,
        second.stress.yy,
        second.stress.zz,
        second.stress.xy,
        second.stress.yz,
        second.stress.xz,
        second.equivalent_plastic_strain,
        second.equivalent_creep_strain};
    for (std::size_t component = 0; component < first_scalars.size(); ++component) {
        error = std::max(error, std::abs(first_scalars[component] - second_scalars[component]));
        scale = std::max({scale, std::abs(first_scalars[component]), std::abs(second_scalars[component])});
    }
    return error / scale;
}

bool test_material_value_paths() {
    fuelsim::ThermoelasticProperties elastic =
        fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e5, 0.25, 1.2e-5, 300.0, -10.0, 1.0e-5, 2.0e-8);
    std::array<fuelsim::ThermoelasticProperties, 4> properties = {elastic,
        fuelsim::test::with_plasticity(elastic, 20.0, 10.0, 300.0, -0.01, 0.02),
        fuelsim::test::with_norton(elastic, 1.0e-6, 10.0, 3.0, 300.0, 1.0e-8, 0.01, 1.0e-3),
        fuelsim::test::with_plasticity(
            fuelsim::test::with_norton(elastic, 1.0e-6, 10.0, 3.0, 300.0, 1.0e-8, 0.01, 1.0e-3),
            20.0,
            10.0,
            300.0,
            -0.01,
            0.02)};
    const fuelsim::SymmetricTensor3Values strain{0.02, -0.004, 0.002, 0.003, -0.001, 0.002};
    const fuelsim::SymmetricTensor3 active_strain{strain.xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz};
    fuelsim::CartesianRotation rotation;
    constexpr double angle = 0.31;
    rotation.xx = std::cos(angle);
    rotation.xy = -std::sin(angle);
    rotation.yx = std::sin(angle);
    rotation.yy = std::cos(angle);
    double maximum_error = 0.0;
    for (const fuelsim::ThermoelasticProperties& property : properties) {
        const fuelsim::IsotropicThermoelasticMaterial model(property);
        const fuelsim::CartesianMaterialPointState committed{};
        const fuelsim::CartesianMaterialPointState active =
            model.response(active_strain, adlite::Scalar(315.0), 0.5, committed).trial_state;
        const fuelsim::CartesianMaterialPointState values = model.response_values(strain, 315.0, 0.5, committed);
        maximum_error = std::max(maximum_error, material_state_error(active, values));
        const fuelsim::CartesianMaterialPointState active_incremental =
            model.incremental_response(active_strain, rotation, adlite::Scalar(325.0), 315.0, 0.5, active).trial_state;
        const fuelsim::CartesianMaterialPointState values_incremental =
            model.incremental_response_values(strain, rotation, 325.0, 315.0, 0.5, active);
        maximum_error = std::max(maximum_error, material_state_error(active_incremental, values_incremental));
    }
    std::cout << "hex20_material_value_path_relative_error=" << maximum_error << '\n';
    return check(maximum_error < 2.0e-14,
        "ordinary-double Cartesian material values match passive automatic differentiation for all built-in branches");
}

bool test_geometry_and_constant_strain() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::test::make_c3d20_geometry(coordinates);
    double thermal_volume = 0.0, mechanical_volume = 0.0;
    for (const fuelsim::Hex20ThermalQuadraturePoint& point : geometry.thermal_points) {
        thermal_volume += point.weighted_measure;
        double temperature_sum = 0.0;
        std::array<double, 3> temperature_gradient{};
        for (std::size_t node = 0; node < 8; ++node) {
            temperature_sum += point.temperature_shape[node];
            for (std::size_t direction = 0; direction < 3; ++direction)
                temperature_gradient[direction] += point.temperature_gradient[node][direction];
        }
        if (!check(near(temperature_sum, 1.0, 2.0e-14), "HEX20 thermal shape functions form a partition of unity")
            || !check(std::max({std::abs(temperature_gradient[0]),
                          std::abs(temperature_gradient[1]),
                          std::abs(temperature_gradient[2])})
                          < 2.0e-14,
                "HEX20 thermal shape gradients sum to zero"))
            return false;
    }
    for (const fuelsim::Hex20MechanicalQuadraturePoint& point : geometry.mechanical_points) {
        mechanical_volume += point.weighted_measure;
        double temperature_sum = 0.0, displacement_sum = 0.0;
        std::array<double, 3> displacement_gradient{};
        for (std::size_t node = 0; node < 8; ++node)
            temperature_sum += point.temperature_shape[node];
        for (std::size_t node = 0; node < 20; ++node) {
            displacement_sum += point.displacement_shape[node];
            for (std::size_t direction = 0; direction < 3; ++direction)
                displacement_gradient[direction] += point.displacement_gradient[node][direction];
        }
        if (!check(near(temperature_sum, 1.0, 2.0e-14) && near(displacement_sum, 1.0, 2.0e-14),
                "HEX20 mechanical points carry separate U2 and T1 partitions of unity")
            || !check(std::max({std::abs(displacement_gradient[0]),
                          std::abs(displacement_gradient[1]),
                          std::abs(displacement_gradient[2])})
                          < 2.0e-14,
                "HEX20 mechanical displacement shape gradients sum to zero"))
            return false;
    }
    if (!check(near(thermal_volume, 1.0, 2.0e-14), "HEX20 thermal 8-point integration recovers unit volume")
        || !check(near(mechanical_volume, 1.0, 2.0e-14), "HEX20 mechanical 27-point integration recovers unit volume"))
        return false;
    fuelsim::Hex20LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node)
        state[node] = 300.0;
    for (std::size_t node = 0; node < 20; ++node) {
        const auto& point = coordinates[node];
        state[8 + node] = 0.01 * point.x + 0.004 * point.y + 0.008 * point.z;
        state[28 + node] = 0.004 * point.x - 0.02 * point.y - 0.006 * point.z;
        state[48 + node] = 0.008 * point.x - 0.006 * point.y + 0.03 * point.z;
    }
    const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(material()), 0.0, 0.0};
    const auto stresses = fuelsim::compute_c3d20_stress(data, geometry, state);
    const double lambda = 2.0e5 * 0.25 / (1.25 * 0.5), shear = 2.0e5 / 2.5, trace = 0.02;
    for (const auto& stress : stresses)
        if (!check(near(stress.xx, lambda * trace + 2.0 * shear * 0.01, 3.0e-13)
                       && near(stress.yy, lambda * trace - 2.0 * shear * 0.02, 3.0e-13)
                       && near(stress.zz, lambda * trace + 2.0 * shear * 0.03, 3.0e-13)
                       && near(stress.xy, 2.0 * shear * 0.004, 3.0e-13)
                       && near(stress.yz, -2.0 * shear * 0.006, 3.0e-13)
                       && near(stress.xz, 2.0 * shear * 0.008, 3.0e-13),
                "quadratic displacement interpolation reproduces all constant-strain stresses"))
            return false;
    return true;
}

bool test_finite_thermal_operators() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::test::make_c3d20_geometry(coordinates);
    const fuelsim::ThermoelasticProperties properties =
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e5, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2.0, 3.0);
    const fuelsim::CartesianTestData conduction_data{fuelsim::IsotropicThermoelasticMaterial(properties),
        0.0,
        1.0,
        fuelsim::StrainFormulation::finite,
        fuelsim::Hex8ElementFormulation::c3d8t,
        300.0};
    fuelsim::Hex20LocalValues state{};
    constexpr double gradient_x = 7.0, gradient_y = -5.0, gradient_z = 3.0;
    for (std::size_t node = 0; node < 8; ++node)
        state[node] = 300.0 + gradient_x * coordinates[node].x + gradient_y * coordinates[node].y
                      + gradient_z * coordinates[node].z;
    for (std::size_t node = 0; node < 20; ++node) {
        state[8 + node] = coordinates[node].x;
        state[28 + node] = 0.2 * coordinates[node].y;
    }
    const auto conduction = fuelsim::compute_c3d20_thermoelastic(conduction_data, geometry, state);
    bool passed = true;
    for (std::size_t node = 0; node < 8; ++node) {
        const double sign_x = 2.0 * coordinates[node].x - 1.0, sign_y = 2.0 * coordinates[node].y - 1.0,
                     sign_z = 2.0 * coordinates[node].z - 1.0;
        const double expected =
            4.0 * 2.4 * 0.25
            * (sign_x * gradient_x / (1.5 * 1.5) + sign_y * gradient_y / (1.1 * 1.1) + sign_z * gradient_z);
        passed = check(near(conduction[node], expected, 2.0e-13),
                     "finite-strain HEX20 conduction uses midpoint gradients and the current quadratic volume")
                 && passed;
    }

    fuelsim::Hex20LocalValues uniform = state;
    for (std::size_t node = 0; node < 8; ++node)
        uniform[node] = 300.0;
    const fuelsim::CartesianTestData source_data{fuelsim::IsotropicThermoelasticMaterial(properties),
        10.0,
        1.0,
        fuelsim::StrainFormulation::finite};
    const auto source = fuelsim::compute_c3d20_thermoelastic(source_data, geometry, uniform);
    const double source_sum = std::accumulate(source.begin(), source.begin() + 8, 0.0);
    passed =
        check(near(source_sum, -24.0, 2.0e-13), "finite-strain HEX20 body source uses the current eight-corner volume")
        && passed;

    fuelsim::Hex20LocalValues midside = uniform;
    for (std::size_t node = 0; node < 20; ++node) {
        midside[8 + node] = 0.0;
        midside[28 + node] = 0.0;
    }
    midside[8 + 8] = 0.05;
    const auto midside_source = fuelsim::compute_c3d20_thermoelastic(source_data, geometry, midside);
    const double midside_source_sum = std::accumulate(midside_source.begin(), midside_source.begin() + 8, 0.0);
    passed = check(near(midside_source_sum, -10.0, 2.0e-13),
                 "finite-strain HEX20 body source excludes displacement midside nodes from its volume")
             && passed;
    const fuelsim::CartesianTestData small_source_data{fuelsim::IsotropicThermoelasticMaterial(properties),
        10.0,
        1.0,
        fuelsim::StrainFormulation::small};
    const auto small_source = fuelsim::compute_c3d20_thermoelastic(small_source_data, geometry, midside);
    passed = check(near(std::accumulate(small_source.begin(), small_source.begin() + 8, 0.0), -10.0, 2.0e-13),
                 "small-strain HEX20 body source retains the reference-volume operator")
             && passed;

    fuelsim::Hex20LocalValues old = state;
    for (std::size_t node = 0; node < 8; ++node) {
        old[node] = 300.0;
        state[node] = 301.0 + static_cast<double>(node);
    }
    for (std::size_t node = 0; node < 20; ++node) {
        old[8 + node] = 0.0;
        old[28 + node] = 0.0;
        old[48 + node] = 0.0;
    }
    const fuelsim::CartesianMaterialHistory history(27);
    const auto transient =
        fuelsim::compute_c3d20_transient(conduction_data, geometry, state, old, history, 1.0, nullptr, true);
    const auto steady =
        fuelsim::compute_c3d20_transient(conduction_data, geometry, state, old, history, 1.0, nullptr, false);
    for (std::size_t node = 0; node < 8; ++node) {
        double expected_capacity = 0.0;
        for (std::size_t other = 0; other < 8; ++other) {
            double mass = 1.0;
            for (std::size_t direction = 0; direction < 3; ++direction) {
                const double left = direction == 0   ? coordinates[node].x
                                    : direction == 1 ? coordinates[node].y
                                                     : coordinates[node].z;
                const double right = direction == 0   ? coordinates[other].x
                                     : direction == 1 ? coordinates[other].y
                                                      : coordinates[other].z;
                mass *= left == right ? 2.0 : 1.0;
            }
            expected_capacity += mass * static_cast<double>(other + 1);
        }
        expected_capacity *= 6.0 / 216.0;
        passed = check(near(transient[node] - steady[node], expected_capacity, 3.0e-12),
                     "finite-strain HEX20 capacity uses the consistent initial-mass matrix")
                 && passed;
    }

    fuelsim::Hex20LocalValues invalid = uniform;
    for (std::size_t node = 0; node < 20; ++node) {
        invalid[8 + node] = -2.0 * coordinates[node].x;
        invalid[28 + node] = -2.0 * coordinates[node].y;
    }
    bool midpoint_rejected = false;
    try {
        (void)fuelsim::compute_c3d20_thermoelastic(conduction_data, geometry, invalid);
    } catch (const std::domain_error& error) {
        midpoint_rejected = std::string(error.what()).find("midpoint") != std::string::npos;
    }
    return check(midpoint_rejected, "finite-strain HEX20 thermal integration rejects a singular midpoint geometry")
           && passed;
}

bool test_jacobian_and_transient_history() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::test::make_c3d20_geometry(coordinates);
    const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(material(true)),
        4.0e5,
        1.0,
        fuelsim::StrainFormulation::small,
        fuelsim::Hex8ElementFormulation::c3d8t,
        300.0};
    fuelsim::Hex20LocalValues old{}, state{};
    for (std::size_t node = 0; node < 8; ++node)
        old[node] = state[node] = 300.0 + 2.0 * static_cast<double>(node);
    for (std::size_t node = 0; node < 20; ++node) {
        const auto& point = coordinates[node];
        state[8 + node] = 0.02 * point.x + 0.003 * point.y;
        state[28 + node] = 0.003 * point.x - 0.004 * point.y;
        state[48 + node] = 0.002 * point.z;
    }
    const fuelsim::CartesianMaterialHistory history(27);
    const double small_error = directional_jacobian_error(data, geometry, state, old, history);
    const double small_residual_error = residual_path_error(data, geometry, state, old, history);
    const fuelsim::CartesianMaterialHistory update =
        fuelsim::compute_c3d20_transient_update(data, geometry, state, old, history, 0.5);
    bool active = false;
    for (const auto& point : update)
        active = active || point.equivalent_plastic_strain > 0.0 || point.equivalent_creep_strain > 0.0;
    const fuelsim::CartesianTestData finite_data{fuelsim::IsotropicThermoelasticMaterial(material(true)),
        4.0e5,
        1.0,
        fuelsim::StrainFormulation::finite,
        fuelsim::Hex8ElementFormulation::c3d8t,
        300.0};
    const double finite_error = directional_jacobian_error(finite_data, geometry, state, old, history);
    const double finite_thermal_error = directional_jacobian_error(finite_data, geometry, state, old, history, 0, 8);
    const double finite_residual_error = residual_path_error(finite_data, geometry, state, old, history);
    fuelsim::Hex20LocalValues invalid = state;
    for (std::size_t node = 0; node < 20; ++node)
        invalid[8 + node] = -2.0 * coordinates[node].x;
    bool invalid_rejected = false;
    try {
        fuelsim::elements::validate_c3d20t_deformation(geometry.mechanical_points[0], invalid);
    } catch (const std::domain_error&) {
        invalid_rejected = true;
    }
    std::cout << "hex20_directional_jacobian_relative_error=" << small_error << '\n'
              << "hex20_finite_directional_jacobian_relative_error=" << finite_error << '\n'
              << "hex20_finite_thermal_directional_jacobian_relative_error=" << finite_thermal_error << '\n'
              << "hex20_residual_path_relative_error=" << small_residual_error << '\n'
              << "hex20_finite_residual_path_relative_error=" << finite_residual_error << '\n';
    return check(small_error < 2.0e-6,
               "68-DOF narrow automatic-differentiation Jacobian matches a centered directional difference")
           && check(finite_error < 3.0e-6, "finite-strain 68-DOF Jacobian matches a centered directional difference")
           && check(finite_thermal_error < 3.0e-6,
               "finite-strain HEX20 thermal rows match a centered directional difference")
           && check(small_residual_error < 2.0e-14 && finite_residual_error < 2.0e-14,
               "ordinary-double HEX20 residual matches the Jacobian-call residual")
           && check(invalid_rejected, "finite-strain HEX20 rejects a nonpositive deformation Jacobian")
           && check(update.size() == 27 && active,
               "HEX20 transient update commits 27 active inelastic material points");
}

bool test_warped_geometry() {
    auto coordinates = unit_cube();
    coordinates[8].x = 0.56;
    coordinates[9].y = 0.46;
    coordinates[14].z = 0.54;
    const fuelsim::Hex20Geometry geometry = fuelsim::test::make_c3d20_geometry(coordinates);
    double thermal_volume = 0.0;
    double mechanical_volume = 0.0;
    for (const auto& point : geometry.thermal_points)
        thermal_volume += point.weighted_measure;
    for (const auto& point : geometry.mechanical_points)
        mechanical_volume += point.weighted_measure;
    return check(thermal_volume > 0.0 && mechanical_volume > 0.0,
               "warped HEX20 geometry has positive thermal and mechanical measures")
           && check(std::isfinite(thermal_volume) && std::isfinite(mechanical_volume),
               "warped HEX20 geometry measures remain finite");
}

int run_c3d20t_tests() {
    const bool passed = test_geometry_and_constant_strain() && test_material_value_paths()
                        && test_finite_thermal_operators() && test_jacobian_and_transient_history()
                        && test_warped_geometry()
                        && test_reference_mass_capacity(fuelsim::Hex20ElementFormulation::c3d20t)

        ;
    if (passed)
        std::cout << "All HEX20-U2/T1 kernel tests passed\n";
    return passed ? 0 : 1;
}
} // namespace

int main() {
    using namespace fuelsim::elements;
    const fuelsim::IsotropicThermoelasticMaterial test_material(material());
    const fuelsim::Hex20LocalValues state{};
    for (const auto selected : {C3d20Quadrature::full, C3d20Quadrature::reduced}) {
        const auto geometry = make_c3d20t_geometry(unit_cube(), selected);
        const C3d20Input input{test_material, geometry, state, state};
        for (const auto requested :
            {C3d20Quadrature::full, C3d20Quadrature::reduced, static_cast<C3d20Quadrature>(-1)}) {
            bool rejected = false;
            try {
                evaluate_c3d20t(input, {false, false, false, false}, requested);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            if (!check(rejected == (requested != selected), "C3D20 quadrature validation"))
                return 1;
        }
    }
    return run_c3d20t_tests();
}
