#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/contact.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool near(double actual, double expected, double tolerance) {
    return std::abs(actual - expected) <= tolerance * std::max({1.0, std::abs(actual), std::abs(expected)});
}

fuelsim::Hex8Coordinates unit_cube() {
    return {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}}};
}

fuelsim::ThermoelasticProperties properties() {
    return fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, -1.0e8, 0.0, 1.0e-8, 2000.0, 3000.0);
}

fuelsim::ThermoelasticProperties inelastic_properties(bool creep, bool plasticity) {
    fuelsim::ThermoelasticProperties result = fuelsim::test::thermoelastic(0.0, 1.0, 200.0, 0.25, 0.0, 300.0);
    if (creep) result = fuelsim::test::with_norton(std::move(result), 1.0e-4, 10.0, 3.0, 300.0);
    if (plasticity) result = fuelsim::test::with_plasticity(std::move(result), 10.0, 20.0, 300.0);
    return result;
}

double equivalent_stress(const fuelsim::SymmetricTensor3& stress) {
    const double mean = (stress.xx.value() + stress.yy.value() + stress.zz.value()) / 3.0;
    const double xx = stress.xx.value() - mean, yy = stress.yy.value() - mean, zz = stress.zz.value() - mean;
    return std::sqrt(1.5 * (xx * xx + yy * yy + zz * zz +
                               2.0 * (stress.xy.value() * stress.xy.value() + stress.yz.value() * stress.yz.value() +
                                         stress.xz.value() * stress.xz.value())));
}

bool test_geometry_and_constant_strain() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
    double volume = 0.0;
    for (const fuelsim::Hex8QuadraturePoint& point : geometry.points) {
        volume += point.weighted_measure;
        double shape_sum = 0.0;
        std::array<double, 3> gradient_sum{};
        for (std::size_t node = 0; node < 8; ++node) {
            shape_sum += point.shape[node];
            for (std::size_t component = 0; component < 3; ++component)
                gradient_sum[component] += point.gradient[node][component];
        }
        if (!check(near(shape_sum, 1.0, 1.0e-14), "HEX8 shape functions form a partition of unity") ||
            !check(near(gradient_sum[0], 0.0, 1.0e-14) && near(gradient_sum[1], 0.0, 1.0e-14) &&
                       near(gradient_sum[2], 0.0, 1.0e-14),
                "HEX8 shape gradients sum to zero"))
            return false;
    }
    if (!check(near(volume, 1.0, 1.0e-14), "HEX8 eight-point integration recovers unit volume")) return false;
    fuelsim::Hex8LocalValues state{};
    const double exx = 0.01;
    const double eyy = -0.02;
    const double ezz = 0.03;
    const double exy = 0.004;
    const double eyz = -0.006;
    const double exz = 0.008;
    for (std::size_t node = 0; node < 8; ++node) {
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[node] = 300.0;
        state[8 + node] = exx * point.x + exy * point.y + exz * point.z;
        state[16 + node] = exy * point.x + eyy * point.y + eyz * point.z;
        state[24 + node] = exz * point.x + eyz * point.y + ezz * point.z;
    }
    const fuelsim::Hex8ThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 0.0, 0.0};
    const auto stresses = fuelsim::compute_hex8_stress(data, geometry, state);
    const double lambda = 2.0e11 * 0.25 / (1.25 * 0.5);
    const double shear = 2.0e11 / 2.5;
    const double trace = exx + eyy + ezz;
    for (const fuelsim::SymmetricTensor3Values& stress : stresses)
        if (!check(near(stress.xx, lambda * trace + 2.0 * shear * exx, 2.0e-13) &&
                       near(stress.yy, lambda * trace + 2.0 * shear * eyy, 2.0e-13) &&
                       near(stress.zz, lambda * trace + 2.0 * shear * ezz, 2.0e-13) &&
                       near(stress.xy, 2.0 * shear * exy, 2.0e-13) && near(stress.yz, 2.0 * shear * eyz, 2.0e-13) &&
                       near(stress.xz, 2.0 * shear * exz, 2.0e-13),
                "HEX8 reproduces all six constant-strain stress components"))
            return false;
    return true;
}

bool test_free_thermal_expansion_and_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
    const fuelsim::Hex8ThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 7.0e5, 0.0};
    fuelsim::Hex8LocalValues state{};
    const double temperature = 650.0;
    const double active_alpha = 1.2e-5 + 1.0e-8 * (temperature - 300.0);
    const double strain = active_alpha * (temperature - 300.0);
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = temperature;
        state[8 + node] = strain * coordinates[node].x;
        state[16 + node] = strain * coordinates[node].y;
        state[24 + node] = strain * coordinates[node].z;
    }
    for (const fuelsim::SymmetricTensor3Values& stress : fuelsim::compute_hex8_stress(data, geometry, state))
        if (!check(std::max({std::abs(stress.xx), std::abs(stress.yy), std::abs(stress.zz), std::abs(stress.xy),
                       std::abs(stress.yz), std::abs(stress.xz)}) < 1.0e-4,
                "uniform three-dimensional thermal expansion is stress free"))
            return false;
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        state[dof] += dof < 8 ? 2.0 * static_cast<double>(dof) : 1.0e-5 * static_cast<double>(dof + 1);
    std::array<double, 32> direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.37 * static_cast<double>(dof + 1));
    fuelsim::Hex8LocalJacobian jacobian{};
    (void)fuelsim::compute_hex8_thermoelastic(data, geometry, state, nullptr, 0.0, &jacobian);
    const double epsilon = 1.0e-7;
    fuelsim::Hex8LocalValues plus = state;
    fuelsim::Hex8LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += epsilon * direction[dof];
        minus[dof] -= epsilon * direction[dof];
    }
    const fuelsim::Hex8LocalResidual plus_residual = fuelsim::compute_hex8_thermoelastic(data, geometry, plus);
    const fuelsim::Hex8LocalResidual minus_residual = fuelsim::compute_hex8_thermoelastic(data, geometry, minus);
    double maximum_error = 0.0;
    double scale = 0.0;
    for (std::size_t row = 0; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column) analytic += jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * epsilon);
        maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    return check(maximum_error / scale < 3.0e-7,
        "full 32-DOF HEX8 automatic-differentiation Jacobian matches a centered directional difference");
}

bool test_transient_capacity_and_faces() {
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
    const fuelsim::Hex8ThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 1.2e7, 0.0};
    fuelsim::Hex8LocalValues old_state{};
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        old_state[node] = 300.0;
        state[node] = 302.0;
    }
    const fuelsim::Hex8LocalResidual residual =
        fuelsim::compute_hex8_thermoelastic(data, geometry, state, &old_state, 1.0);
    for (std::size_t node = 0; node < 8; ++node)
        if (!check(std::abs(residual[node]) < 1.0e-7,
                "Backward Euler consistent heat capacity balances uniform volumetric heating; residual=" +
                    std::to_string(residual[node])))
            return false;
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Quad4FaceCoordinates face_coordinates = {
        {coordinates[1], coordinates[2], coordinates[6], coordinates[5]}};
    const fuelsim::Quad4FaceGeometry face = fuelsim::make_quad4_face_geometry(face_coordinates);
    fuelsim::Quad4FaceLocalValues face_state{};
    const fuelsim::Quad4FaceBoundaryData pressure = {
        fuelsim::Quad4FaceBoundaryKind::pressure, fuelsim::CartesianTractionComponent::x, 5.0, 0.0};
    const fuelsim::Quad4FaceLocalResidual pressure_residual =
        fuelsim::compute_quad4_face_boundary(pressure, face, face_state);
    double force_x = 0.0;
    double force_y = 0.0;
    double force_z = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        force_x += pressure_residual[4 + node];
        force_y += pressure_residual[8 + node];
        force_z += pressure_residual[12 + node];
    }
    if (!check(near(force_x, 5.0, 1.0e-14) && near(force_y, 0.0, 1.0e-14) && near(force_z, 0.0, 1.0e-14),
            "reference pressure uses the outward three-dimensional face area vector and exact total force"))
        return false;
    face_state.fill(0.0);
    for (std::size_t node = 0; node < 4; ++node) {
        face_state[4 + node] = 0.2 * face_coordinates[node].y;
        face_state[8 + node] = 0.1 * face_coordinates[node].y;
        face_state[12 + node] = -0.05 * face_coordinates[node].z;
    }
    const fuelsim::Quad4FaceBoundaryData follower = {
        fuelsim::Quad4FaceBoundaryKind::pressure, fuelsim::CartesianTractionComponent::x, 5.0, 0.0, true};
    fuelsim::Quad4FaceLocalJacobian follower_jacobian{};
    const auto follower_residual = fuelsim::compute_quad4_face_boundary(follower, face, face_state, &follower_jacobian);
    force_x = 0.0;
    force_y = 0.0;
    force_z = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        force_x += follower_residual[4 + node];
        force_y += follower_residual[8 + node];
        force_z += follower_residual[12 + node];
    }
    if (!check(near(force_x, 5.225, 1.0e-14) && near(force_y, -0.95, 1.0e-14) && near(force_z, 0.0, 1.0e-14),
            "three-dimensional follower pressure uses the complete current area vector and exact resultant"))
        return false;
    std::array<double, 16> follower_direction{};
    for (std::size_t dof = 4; dof < 16; ++dof) follower_direction[dof] = std::cos(0.41 * static_cast<double>(dof + 1));
    constexpr double follower_step = 1.0e-7;
    auto follower_plus = face_state, follower_minus = face_state;
    for (std::size_t dof = 0; dof < 16; ++dof) {
        follower_plus[dof] += follower_step * follower_direction[dof];
        follower_minus[dof] -= follower_step * follower_direction[dof];
    }
    const auto follower_plus_residual = fuelsim::compute_quad4_face_boundary(follower, face, follower_plus);
    const auto follower_minus_residual = fuelsim::compute_quad4_face_boundary(follower, face, follower_minus);
    double follower_error = 0.0;
    for (std::size_t row = 4; row < 16; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 16; ++column)
            analytic += follower_jacobian[row * 16 + column] * follower_direction[column];
        const double numerical = (follower_plus_residual[row] - follower_minus_residual[row]) / (2.0 * follower_step);
        follower_error = std::max(follower_error, std::abs(analytic - numerical));
    }
    if (!check(follower_error < 2.0e-9,
            "three-dimensional follower-pressure geometric Jacobian matches centered difference"))
        return false;
    for (std::size_t node = 0; node < 4; ++node) face_state[node] = 350.0;
    const fuelsim::Quad4FaceBoundaryData convection = {
        fuelsim::Quad4FaceBoundaryKind::convection, fuelsim::CartesianTractionComponent::x, 20.0, 300.0};
    fuelsim::Quad4FaceLocalJacobian convection_jacobian{};
    const fuelsim::Quad4FaceLocalResidual convection_residual =
        fuelsim::compute_quad4_face_boundary(convection, face, face_state, &convection_jacobian);
    double heat = 0.0;
    double tangent_sum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        heat += convection_residual[row];
        for (std::size_t column = 0; column < 4; ++column) tangent_sum += convection_jacobian[row * 16 + column];
    }
    return check(near(heat, 1000.0, 1.0e-14) && near(tangent_sum, 20.0, 1.0e-14),
        "three-dimensional convection has the exact face heat rate and consistent temperature tangent");
}

bool test_cartesian_inelastic_material() {
    const fuelsim::SymmetricTensor3 strain{0.20, -0.04, -0.03, 0.02, -0.015, 0.01};
    const fuelsim::CartesianMaterialPointState committed{};
    bool passed = true;
    for (const std::array<bool, 2> branch : {std::array<bool, 2>{false, true}, {true, false}, {true, true}}) {
        const fuelsim::IsotropicThermoelasticMaterial material(inelastic_properties(branch[0], branch[1]));
        const fuelsim::CartesianInelasticStressResponse response = material.response(strain, 300.0, 1.0, committed);
        const fuelsim::CartesianMaterialPointState state = response.trial_state;
        const double plastic_trace = state.plastic_strain[0] + state.plastic_strain[1] + state.plastic_strain[2];
        const double creep_trace = state.creep_strain[0] + state.creep_strain[1] + state.creep_strain[2];
        passed = check((branch[1] ? state.equivalent_plastic_strain > 0.0 : state.equivalent_plastic_strain == 0.0) &&
                           (branch[0] ? state.equivalent_creep_strain > 0.0 : state.equivalent_creep_strain == 0.0) &&
                           std::abs(plastic_trace) < 2.0e-15 && std::abs(creep_trace) < 2.0e-15 &&
                           std::isfinite(equivalent_stress(response.stress)),
                     "Cartesian plastic, creep, and coupled updates activate the requested traceless branches") &&
                 passed;
        const adlite::Scalar active_xx = adlite::Scalar::independent(strain.xx.value(), 0, 1);
        const fuelsim::CartesianInelasticStressResponse active = material.response(
            {active_xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz}, 300.0, 1.0, committed);
        constexpr double step = 1.0e-7;
        const double plus =
            material
                .response({strain.xx.value() + step, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz}, 300.0, 1.0,
                    committed)
                .stress.xx.value();
        const double minus =
            material
                .response({strain.xx.value() - step, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz}, 300.0, 1.0,
                    committed)
                .stress.xx.value();
        const double numerical = (plus - minus) / (2.0 * step);
        passed =
            check(std::abs(active.stress.xx.derivative(0) - numerical) / std::max(1.0, std::abs(numerical)) < 2.0e-7,
                "Cartesian inelastic material automatic-differentiation tangent matches centered difference") &&
            passed;
    }
    return passed;
}

bool test_finite_strain_kinematics_and_coupled_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
    fuelsim::Hex8LocalValues state{};
    const double stretch_x = 1.12, stretch_y = 0.94, stretch_z = 1.03;
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = 300.0;
        state[8 + node] = (stretch_x - 1.0) * coordinates[node].x;
        state[16 + node] = (stretch_y - 1.0) * coordinates[node].y;
        state[24 + node] = (stretch_z - 1.0) * coordinates[node].z;
    }
    fuelsim::Hex8LocalAdValues passive{};
    for (std::size_t dof = 0; dof < state.size(); ++dof) passive[dof] = state[dof];
    const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
        geometry.points[0], passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::finite);
    const auto taylor = [](double stretch) {
        const double cinv_minus_one = 1.0 / (stretch * stretch) - 1.0;
        return -0.5 * cinv_minus_one + 0.25 * cinv_minus_one * cinv_minus_one;
    };
    bool passed = check(near(kinematics.strain_increment.xx.value(), taylor(stretch_x), 2.0e-14) &&
                            near(kinematics.strain_increment.yy.value(), taylor(stretch_y), 2.0e-14) &&
                            near(kinematics.strain_increment.zz.value(), taylor(stretch_z), 2.0e-14) &&
                            near(kinematics.current_weighted_measure.value(),
                                geometry.points[0].weighted_measure * stretch_x * stretch_y * stretch_z, 2.0e-14),
        "finite-strain HEX8 recovers the MOOSE Taylor diagonal increment and current volume measure");
    for (std::size_t node = 0; node < 8; ++node) {
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[8 + node] = 0.20 * point.x + 0.08 * point.y - 0.03 * point.z;
        state[16 + node] = -0.02 * point.x - 0.04 * point.y + 0.06 * point.z;
        state[24 + node] = 0.04 * point.x - 0.05 * point.y - 0.03 * point.z;
    }
    const fuelsim::Hex8ThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(inelastic_properties(true, true)),
        0.0, 1.0, fuelsim::StrainFormulation::finite};
    const fuelsim::Hex8LocalValues committed_state = [] {
        fuelsim::Hex8LocalValues value{};
        for (std::size_t node = 0; node < 8; ++node) value[node] = 300.0;
        return value;
    }();
    const fuelsim::Hex8MaterialHistory committed_material{};
    fuelsim::Hex8LocalJacobian jacobian{};
    (void)fuelsim::compute_hex8_transient(data, geometry, state, committed_state, committed_material, 1.0, &jacobian);
    std::array<double, 32> direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = dof < 8 ? 0.0 : std::sin(0.29 * static_cast<double>(dof + 1));
    constexpr double step = 2.0e-7;
    fuelsim::Hex8LocalValues plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const auto plus_residual =
        fuelsim::compute_hex8_transient(data, geometry, plus, committed_state, committed_material, 1.0);
    const auto minus_residual =
        fuelsim::compute_hex8_transient(data, geometry, minus, committed_state, committed_material, 1.0);
    double maximum_error = 0.0, scale = 0.0;
    for (std::size_t row = 8; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column) analytic += jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(maximum_error / scale < 2.0e-6,
                 "finite-strain coupled HEX8 automatic-differentiation Jacobian matches centered difference") &&
             passed;
    fuelsim::Hex8LocalValues inverted = state;
    for (std::size_t node = 0; node < 8; ++node) inverted[8 + node] = -2.0 * coordinates[node].x;
    try {
        (void)fuelsim::compute_hex8_transient(data, geometry, inverted, committed_state, committed_material, 1.0);
        passed = check(false, "finite-strain HEX8 rejects a nonpositive deformation Jacobian") && passed;
    } catch (const std::domain_error&) {}
    fuelsim::CartesianMaterialPointState history;
    history.elastic_strain = {0.01, -0.02, 0.03, 0.004, -0.005, 0.006};
    history.plastic_strain = {0.02, -0.01, -0.01, 0.007, 0.008, -0.009};
    history.creep_strain = {-0.03, 0.01, 0.02, -0.011, 0.012, 0.013};
    history.equivalent_plastic_strain = 0.04;
    history.equivalent_creep_strain = 0.05;
    const fuelsim::CartesianRotation quarter_turn{0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    const auto rotated =
        fuelsim::IsotropicThermoelasticMaterial(inelastic_properties(false, false))
            .incremental_response({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, quarter_turn, 300.0, 300.0, 1.0, history)
            .trial_state;
    const std::array<std::array<double, 6>, 3> expected = {{{-0.02, 0.01, 0.03, -0.004, 0.006, 0.005},
        {-0.01, 0.02, -0.01, -0.007, -0.009, -0.008}, {0.01, -0.03, 0.02, 0.011, 0.013, -0.012}}};
    const std::array<const std::array<double, 6>*, 3> actual = {
        &rotated.elastic_strain, &rotated.plastic_strain, &rotated.creep_strain};
    bool objective = true;
    for (std::size_t tensor_index = 0; tensor_index < actual.size(); ++tensor_index)
        for (std::size_t component = 0; component < 6; ++component)
            objective =
                objective && near((*actual[tensor_index])[component], expected[tensor_index][component], 1.0e-14);
    passed =
        check(objective && rotated.equivalent_plastic_strain == history.equivalent_plastic_strain &&
                  rotated.equivalent_creep_strain == history.equivalent_creep_strain,
            "finite-strain Cartesian elastic, plastic, and creep tensors rotate objectively while scalars do not") &&
        passed;
    return passed;
}

bool test_cartesian_surface_contact_kernels() {
    const fuelsim::Quad4FaceCoordinates secondary = {{{0.99, 0.0, 0.0}, {0.99, 1.0, 0.0}, {0.99, 1.0, 1.0},
                                            {0.99, 0.0, 1.0}}},
                                        primary = {
                                            {{1.0, -0.5, -0.5}, {1.0, 1.5, -0.5}, {1.0, 1.5, 1.5}, {1.0, -0.5, 1.5}}};
    const fuelsim::Quad4FaceGeometry secondary_face = fuelsim::make_quad4_face_geometry(secondary);
    const fuelsim::Quad4FaceQuadraturePoint& point = secondary_face.points[0];
    const fuelsim::Quad4ToQuad4HeatGeometry heat_geometry{
        secondary, primary, point.shape, point.derivative_xi, point.derivative_eta, 1.0};
    fuelsim::Quad4SurfaceContactLocalValues state{}, committed{};
    for (std::size_t node = 0; node < 4; ++node) {
        state[node] = 400.0;
        state[4 + node] = 300.0;
        committed[node] = 400.0;
        committed[4 + node] = 300.0;
        state[8 + node] = 0.02;
        state[16 + node] = 0.001;
    }
    fuelsim::Quad4SurfaceContactLocalJacobian heat_jacobian{};
    const fuelsim::Quad4SurfaceContactLocalResidual heat =
        fuelsim::compute_quad4_to_quad4_gap_heat({0.2, 0.001}, heat_geometry, state, &heat_jacobian);
    double secondary_heat = 0.0, primary_heat = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        secondary_heat += heat[node];
        primary_heat += heat[4 + node];
    }
    bool passed = check(near(secondary_heat, 5000.0, 1.0e-12) && near(primary_heat, -5000.0, 1.0e-12),
        "three-dimensional gap heat transfer is exactly conservative on the secondary quadrature point");
    std::array<std::array<double, 4>, 4> shapes{}, derivatives_xi{}, derivatives_eta{};
    for (std::size_t q = 0; q < 4; ++q) {
        shapes[q] = secondary_face.points[q].shape;
        derivatives_xi[q] = secondary_face.points[q].derivative_xi;
        derivatives_eta[q] = secondary_face.points[q].derivative_eta;
    }
    const fuelsim::NodeToQuad4ContactGeometry contact_geometry{
        secondary, primary, shapes, derivatives_xi, derivatives_eta, 0, 1.0};
    fuelsim::Quad4SurfaceContactLocalJacobian contact_jacobian{};
    const fuelsim::NormalContactProperties stick_properties{1000.0, 1.0, false};
    const fuelsim::ContactPointHistory history{};
    const fuelsim::Quad4SurfaceContactLocalResidual contact = fuelsim::compute_node_to_quad4_contact(
        stick_properties, contact_geometry, state, committed, history, &contact_jacobian);
    const fuelsim::CartesianContactPointValue stick =
        fuelsim::compute_node_to_quad4_contact_value(stick_properties, contact_geometry, state, committed, history);
    std::array<double, 3> resultant{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 8; ++node) resultant[component] += contact[8 * (component + 1) + node];
    passed =
        check(stick.projected && near(stick.gap, -0.01, 1.0e-12) && near(stick.pressure, 10.0, 1.0e-12) &&
                  near(stick.tributary_area, 0.25, 1.0e-12) && near(stick.contact_force, 2.5, 1.0e-12) &&
                  near(stick.tangential_traction, 1.0, 1.0e-12) && !stick.sliding && near(resultant[0], 0.0, 1.0e-12) &&
                  near(resultant[1], 0.0, 1.0e-12) && near(resultant[2], 0.0, 1.0e-12),
            "three-dimensional node-to-face contact recovers pressure, tributary area, stick traction, and "
            "equal-and-opposite force") &&
        passed;
    const fuelsim::CartesianContactPointValue sliding = fuelsim::compute_node_to_quad4_contact_value(
        {1000.0, 0.05, false}, contact_geometry, state, committed, history);
    passed = check(sliding.sliding && near(sliding.tangential_traction, 0.5, 1.0e-12) &&
                       near(sliding.elastic_tangential_slip[1], 0.0005, 1.0e-12),
                 "three-dimensional Coulomb contact caps sliding traction and stores the vector elastic slip") &&
             passed;
    std::array<double, 32> direction{};
    for (std::size_t dof = 8; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.37 * static_cast<double>(dof + 1));
    constexpr double step = 1.0e-7;
    fuelsim::Quad4SurfaceContactLocalValues plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Quad4SurfaceContactLocalResidual plus_contact = fuelsim::compute_node_to_quad4_contact(
                                                        stick_properties, contact_geometry, plus, committed, history),
                                                    minus_contact = fuelsim::compute_node_to_quad4_contact(
                                                        stick_properties, contact_geometry, minus, committed, history);
    double maximum_error = 0.0, scale = 0.0;
    std::size_t maximum_row = 0;
    double maximum_analytic = 0.0, maximum_numerical = 0.0;
    for (std::size_t row = 8; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column)
            analytic += contact_jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_contact[row] - minus_contact[row]) / (2.0 * step);
        if (std::abs(analytic - numerical) > maximum_error) {
            maximum_error = std::abs(analytic - numerical);
            maximum_row = row;
            maximum_analytic = analytic;
            maximum_numerical = numerical;
        }
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    std::cout << "cartesian_contact_jacobian_relative_error=" << maximum_error / scale << '\n'
              << "cartesian_contact_jacobian_maximum_row=" << maximum_row << '\n'
              << "cartesian_contact_jacobian_maximum_analytic=" << maximum_analytic << '\n'
              << "cartesian_contact_jacobian_maximum_numerical=" << maximum_numerical << '\n';
    passed = check(maximum_error / scale < 2.0e-6,
                 "three-dimensional sticking contact automatic-differentiation Jacobian matches centered difference") &&
             passed;
    return passed;
}
} // namespace

int main() {
    bool passed = true;
    passed = test_geometry_and_constant_strain() && passed;
    passed = test_free_thermal_expansion_and_jacobian() && passed;
    passed = test_transient_capacity_and_faces() && passed;
    passed = test_cartesian_inelastic_material() && passed;
    passed = test_finite_strain_kinematics_and_coupled_jacobian() && passed;
    passed = test_cartesian_surface_contact_kernels() && passed;
    if (!passed) return 1;
    std::cout << "HEX8 thermo-mechanics tests passed\n";
    return 0;
}
