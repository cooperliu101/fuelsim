#include "fuelsim/core/cartesian3d_hex20.hpp"
#include "fuelsim/core/contact.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
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

fuelsim::Hex20Coordinates unit_cube() {
    return {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {0.5, 0.0, 0.0}, {1.0, 0.5, 0.0}, {0.5, 1.0, 0.0}, {0.0, 0.5, 0.0},
        {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {1.0, 1.0, 0.5}, {0.0, 1.0, 0.5}, {0.5, 0.0, 1.0}, {1.0, 0.5, 1.0},
        {0.5, 1.0, 1.0}, {0.0, 0.5, 1.0}}};
}

fuelsim::ThermoelasticProperties material(bool inelastic = false) {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(3000.0, 4.0, 2.0e5, 0.25, 1.2e-5, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
    if (inelastic) {
        result = fuelsim::test::with_norton(std::move(result), 1.0e-6, 10.0, 3.0, 300.0);
        result = fuelsim::test::with_plasticity(std::move(result), 20.0, 10.0, 300.0);
    }
    return result;
}

double directional_jacobian_error(const fuelsim::CartesianThermoelasticData& data,
    const fuelsim::Hex20Geometry& geometry, const fuelsim::Hex20LocalValues& state,
    const fuelsim::Hex20LocalValues& old, const fuelsim::CartesianMaterialHistory& history) {
    fuelsim::Hex20LocalJacobian jacobian{};
    (void)fuelsim::compute_hex20_transient(data, geometry, state, old, history, 0.5, &jacobian);
    fuelsim::Hex20LocalValues direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.29 * static_cast<double>(dof + 1));
    constexpr double step = 2.0e-7;
    auto plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const auto plus_residual = fuelsim::compute_hex20_transient(data, geometry, plus, old, history, 0.5);
    const auto minus_residual = fuelsim::compute_hex20_transient(data, geometry, minus, old, history, 0.5);
    double error = 0.0, scale = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            analytic += jacobian[row * state.size() + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        error = std::max(error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    return error / scale;
}

bool test_geometry_and_constant_strain() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::make_hex20_geometry(coordinates);
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
        if (!check(near(temperature_sum, 1.0, 2.0e-14), "HEX20 thermal shape functions form a partition of unity") ||
            !check(std::max({std::abs(temperature_gradient[0]), std::abs(temperature_gradient[1]),
                       std::abs(temperature_gradient[2])}) < 2.0e-14,
                "HEX20 thermal shape gradients sum to zero"))
            return false;
    }
    for (const fuelsim::Hex20MechanicalQuadraturePoint& point : geometry.mechanical_points) {
        mechanical_volume += point.weighted_measure;
        double temperature_sum = 0.0, displacement_sum = 0.0;
        std::array<double, 3> displacement_gradient{};
        for (std::size_t node = 0; node < 8; ++node) temperature_sum += point.temperature_shape[node];
        for (std::size_t node = 0; node < 20; ++node) {
            displacement_sum += point.displacement_shape[node];
            for (std::size_t direction = 0; direction < 3; ++direction)
                displacement_gradient[direction] += point.displacement_gradient[node][direction];
        }
        if (!check(near(temperature_sum, 1.0, 2.0e-14) && near(displacement_sum, 1.0, 2.0e-14),
                "HEX20 mechanical points carry separate U2 and T1 partitions of unity") ||
            !check(std::max({std::abs(displacement_gradient[0]), std::abs(displacement_gradient[1]),
                       std::abs(displacement_gradient[2])}) < 2.0e-14,
                "HEX20 mechanical displacement shape gradients sum to zero"))
            return false;
    }
    if (!check(near(thermal_volume, 1.0, 2.0e-14), "HEX20 thermal 8-point integration recovers unit volume") ||
        !check(near(mechanical_volume, 1.0, 2.0e-14), "HEX20 mechanical 27-point integration recovers unit volume"))
        return false;
    fuelsim::Hex20LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) state[node] = 300.0;
    for (std::size_t node = 0; node < 20; ++node) {
        const auto& point = coordinates[node];
        state[8 + node] = 0.01 * point.x + 0.004 * point.y + 0.008 * point.z;
        state[28 + node] = 0.004 * point.x - 0.02 * point.y - 0.006 * point.z;
        state[48 + node] = 0.008 * point.x - 0.006 * point.y + 0.03 * point.z;
    }
    const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(material()), 0.0, 0.0};
    const auto stresses = fuelsim::compute_hex20_stress(data, geometry, state);
    const double lambda = 2.0e5 * 0.25 / (1.25 * 0.5), shear = 2.0e5 / 2.5, trace = 0.02;
    for (const auto& stress : stresses)
        if (!check(near(stress.xx, lambda * trace + 2.0 * shear * 0.01, 3.0e-13) &&
                       near(stress.yy, lambda * trace - 2.0 * shear * 0.02, 3.0e-13) &&
                       near(stress.zz, lambda * trace + 2.0 * shear * 0.03, 3.0e-13) &&
                       near(stress.xy, 2.0 * shear * 0.004, 3.0e-13) &&
                       near(stress.yz, -2.0 * shear * 0.006, 3.0e-13) && near(stress.xz, 2.0 * shear * 0.008, 3.0e-13),
                "quadratic displacement interpolation reproduces all constant-strain stresses"))
            return false;
    return true;
}

bool test_jacobian_and_transient_history() {
    const auto coordinates = unit_cube();
    const fuelsim::Hex20Geometry geometry = fuelsim::make_hex20_geometry(coordinates);
    const fuelsim::CartesianThermoelasticData data{
        fuelsim::IsotropicThermoelasticMaterial(material(true)), 4.0e5, 1.0, fuelsim::StrainFormulation::small};
    fuelsim::Hex20LocalValues old{}, state{};
    for (std::size_t node = 0; node < 8; ++node) old[node] = state[node] = 300.0 + 2.0 * static_cast<double>(node);
    for (std::size_t node = 0; node < 20; ++node) {
        const auto& point = coordinates[node];
        state[8 + node] = 0.02 * point.x + 0.003 * point.y;
        state[28 + node] = 0.003 * point.x - 0.004 * point.y;
        state[48 + node] = 0.002 * point.z;
    }
    const fuelsim::CartesianMaterialHistory history(27);
    const double small_error = directional_jacobian_error(data, geometry, state, old, history);
    const fuelsim::CartesianMaterialHistory update =
        fuelsim::compute_hex20_transient_update(data, geometry, state, old, history, 0.5);
    bool active = false;
    for (const auto& point : update)
        active = active || point.equivalent_plastic_strain > 0.0 || point.equivalent_creep_strain > 0.0;
    const fuelsim::CartesianThermoelasticData finite_data{
        fuelsim::IsotropicThermoelasticMaterial(material(true)), 4.0e5, 1.0, fuelsim::StrainFormulation::finite};
    const double finite_error = directional_jacobian_error(finite_data, geometry, state, old, history);
    fuelsim::Hex20LocalValues invalid = state;
    for (std::size_t node = 0; node < 20; ++node) invalid[8 + node] = -2.0 * coordinates[node].x;
    bool invalid_rejected = false;
    try {
        fuelsim::validate_hex20_deformation(geometry.mechanical_points[0], invalid);
    } catch (const std::domain_error&) { invalid_rejected = true; }
    std::cout << "hex20_directional_jacobian_relative_error=" << small_error << '\n'
              << "hex20_finite_directional_jacobian_relative_error=" << finite_error << '\n';
    return check(small_error < 2.0e-6,
               "68-DOF narrow automatic-differentiation Jacobian matches a centered directional difference") &&
           check(finite_error < 3.0e-6, "finite-strain 68-DOF Jacobian matches a centered directional difference") &&
           check(invalid_rejected, "finite-strain HEX20 rejects a nonpositive deformation Jacobian") &&
           check(update.size() == 27 && active, "HEX20 transient update commits 27 active inelastic material points");
}

bool test_quadratic_face() {
    const auto cube = unit_cube();
    const fuelsim::Quad8FaceCoordinates coordinates = {
        cube[1], cube[2], cube[6], cube[5], cube[9], cube[14], cube[17], cube[13]};
    const fuelsim::Quad8FaceGeometry geometry = fuelsim::make_quad8_face_geometry(coordinates);
    double thermal_area = 0.0, mechanical_area = 0.0;
    for (const auto& point : geometry.thermal_points) thermal_area += point.weighted_measure;
    for (const auto& point : geometry.mechanical_points) {
        const auto& tangent_xi = point.tangent_xi;
        const auto& tangent_eta = point.tangent_eta;
        const double area_x = tangent_xi.y * tangent_eta.z - tangent_xi.z * tangent_eta.y;
        const double area_y = tangent_xi.z * tangent_eta.x - tangent_xi.x * tangent_eta.z;
        const double area_z = tangent_xi.x * tangent_eta.y - tangent_xi.y * tangent_eta.x;
        mechanical_area += std::sqrt(area_x * area_x + area_y * area_y + area_z * area_z) * point.quadrature_weight;
    }
    if (!check(geometry.thermal_points.size() == 4 && geometry.mechanical_points.size() == 9,
            "Quad8 boundary separates 2x2 thermal and 3x3 mechanical integration rules") ||
        !check(near(thermal_area, 1.0, 2.0e-14) && near(mechanical_area, 1.0, 2.0e-14),
            "Quad8 thermal and mechanical face rules recover the unit face area"))
        return false;
    fuelsim::Quad8FaceLocalValues state{};
    for (std::size_t node = 0; node < 4; ++node) state[node] = 350.0;
    const fuelsim::Quad4FaceBoundaryData pressure = {
        fuelsim::Quad4FaceBoundaryKind::pressure, fuelsim::CartesianTractionComponent::x, 5.0, 0.0};
    const auto pressure_residual = fuelsim::compute_quad8_face_boundary(pressure, geometry, state);
    double force_x = 0.0, force_y = 0.0, force_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        force_x += pressure_residual[4 + node];
        force_y += pressure_residual[12 + node];
        force_z += pressure_residual[20 + node];
    }
    const fuelsim::Quad4FaceBoundaryData convection = {
        fuelsim::Quad4FaceBoundaryKind::convection, fuelsim::CartesianTractionComponent::x, 20.0, 300.0};
    fuelsim::Quad8FaceLocalJacobian tangent{};
    const auto heat_residual = fuelsim::compute_quad8_face_boundary(convection, geometry, state, &tangent);
    double heat = 0.0, tangent_sum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        heat += heat_residual[row];
        for (std::size_t column = 0; column < 4; ++column)
            tangent_sum += tangent[row * fuelsim::quad8_face_local_dof_count + column];
    }
    return check(near(force_x, 5.0, 2.0e-14) && near(force_y, 0.0, 2.0e-14) && near(force_z, 0.0, 2.0e-14),
               "quadratic face pressure recovers the exact resultant") &&
           check(near(heat, 1000.0, 2.0e-14) && near(tangent_sum, 20.0, 2.0e-14),
               "linear face temperature convection recovers exact heat rate and tangent");
}

bool test_warped_geometry() {
    auto coordinates = unit_cube();
    coordinates[8].x = 0.56;
    coordinates[9].y = 0.46;
    coordinates[14].z = 0.54;
    const fuelsim::Hex20Geometry geometry = fuelsim::make_hex20_geometry(coordinates);
    double thermal_volume = 0.0;
    double mechanical_volume = 0.0;
    for (const auto& point : geometry.thermal_points) thermal_volume += point.weighted_measure;
    for (const auto& point : geometry.mechanical_points) mechanical_volume += point.weighted_measure;
    return check(thermal_volume > 0.0 && mechanical_volume > 0.0,
               "warped HEX20 geometry has positive thermal and mechanical measures") &&
           check(std::isfinite(thermal_volume) && std::isfinite(mechanical_volume),
               "warped HEX20 geometry measures remain finite");
}

bool test_hex20_contact_kernels() {
    const auto cube = unit_cube();
    const fuelsim::Quad8FaceCoordinates face = {
        cube[1], cube[2], cube[6], cube[5], cube[9], cube[14], cube[17], cube[13]};
    const fuelsim::Quad8FaceGeometry geometry = fuelsim::make_quad8_face_geometry(face);
    fuelsim::Quad8SurfaceContactLocalValues state{}, committed{};
    for (std::size_t node = 0; node < 4; ++node) {
        state[node] = 400.0;
        state[4 + node] = 300.0;
    }
    for (std::size_t node = 0; node < 8; ++node) state[8 + node] = -0.01;
    fuelsim::Quad8ToQuad8HeatGeometry heat_geometry{face, face, geometry.thermal_points[0].temperature_shape,
        geometry.thermal_points[0].displacement_shape, geometry.thermal_points[0].derivative_xi,
        geometry.thermal_points[0].derivative_eta, geometry.thermal_points[0].quadrature_weight, -1.0};
    fuelsim::Quad8SurfaceContactLocalJacobian heat_jacobian{};
    const fuelsim::GapHeatProperties heat_properties{2.0, 1.0e-6};
    const auto heat_residual =
        fuelsim::compute_quad8_to_quad8_gap_heat(heat_properties, heat_geometry, state, &heat_jacobian);
    double heat_balance = 0.0;
    for (std::size_t node = 0; node < 8; ++node) heat_balance += heat_residual[node];
    bool passed = check(std::abs(heat_balance) < 1.0e-12 * std::max(1.0, std::abs(heat_residual[0])),
        "HEX20 Q8/Q4 thermal contact residual is exactly conservative at a 2x2 quadrature point");
    const auto heat_value = fuelsim::compute_quad8_to_quad8_gap_heat_value(heat_properties, heat_geometry, state);
    passed = check(heat_value.projected && heat_value.weighted_measure > 0.0 && heat_value.heat_flux > 0.0,
                 "HEX20 thermal contact projects through the quadratic face and produces finite heat flux") &&
             passed;
    auto heat_plus = state, heat_minus = state;
    constexpr double thermal_step = 1.0e-7;
    heat_plus[0] += thermal_step;
    heat_minus[0] -= thermal_step;
    const auto heat_plus_residual = fuelsim::compute_quad8_to_quad8_gap_heat(heat_properties, heat_geometry, heat_plus);
    const auto heat_minus_residual =
        fuelsim::compute_quad8_to_quad8_gap_heat(heat_properties, heat_geometry, heat_minus);
    double thermal_jacobian_error = 0.0, thermal_jacobian_scale = 0.0;
    for (std::size_t row = 0; row < heat_plus_residual.size(); ++row) {
        const double numerical = (heat_plus_residual[row] - heat_minus_residual[row]) / (2.0 * thermal_step);
        const double analytic = heat_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count];
        thermal_jacobian_error = std::max(thermal_jacobian_error, std::abs(analytic - numerical));
        thermal_jacobian_scale = std::max({thermal_jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(thermal_jacobian_error / thermal_jacobian_scale < 2.0e-6,
                 "HEX20 thermal contact automatic-differentiation Jacobian matches a centered difference") &&
             passed;

    fuelsim::NodeToQuad8ContactGeometry mechanical_geometry{face, face, {}, {}, {}, {}, 0, -1.0};
    for (std::size_t q = 0; q < 9; ++q) {
        mechanical_geometry.secondary_shapes[q] = geometry.mechanical_points[q].displacement_shape;
        mechanical_geometry.secondary_derivatives_xi[q] = geometry.mechanical_points[q].derivative_xi;
        mechanical_geometry.secondary_derivatives_eta[q] = geometry.mechanical_points[q].derivative_eta;
        mechanical_geometry.secondary_quadrature_weights[q] = geometry.mechanical_points[q].quadrature_weight;
    }
    const fuelsim::NormalContactProperties mechanical_properties{1.0e5, 0.0, false};
    fuelsim::Quad8SurfaceContactLocalJacobian mechanical_jacobian{};
    const auto mechanical_residual = fuelsim::compute_node_to_quad8_contact(
        mechanical_properties, mechanical_geometry, state, committed, {}, &mechanical_jacobian);
    double force_balance_x = 0.0, force_balance_y = 0.0, force_balance_z = 0.0;
    for (std::size_t node = 0; node < 16; ++node) {
        force_balance_x += mechanical_residual[8 + node];
        force_balance_y += mechanical_residual[24 + node];
        force_balance_z += mechanical_residual[40 + node];
    }
    const auto mechanical_value =
        fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, mechanical_geometry, state, committed, {});
    std::array<double, 8> positive_areas{};
    for (std::size_t node = 0; node < positive_areas.size(); ++node) {
        auto node_geometry = mechanical_geometry;
        node_geometry.secondary_local_node = node;
        positive_areas[node] =
            fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, node_geometry, state, committed, {})
                .tributary_area;
    }
    double positive_area_sum = 0.0;
    for (const double area : positive_areas) positive_area_sum += area;
    auto consistent_geometry = mechanical_geometry;
    consistent_geometry.nodal_area_rule = fuelsim::Quad8NodalAreaRule::consistent_shape;
    const auto consistent_corner =
        fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, consistent_geometry, state, committed, {});
    consistent_geometry.secondary_local_node = 4;
    const auto consistent_midpoint =
        fuelsim::compute_node_to_quad8_contact_value(mechanical_properties, consistent_geometry, state, committed, {});
    passed = check(mechanical_value.projected && mechanical_value.pressure > 0.0 &&
                       near(mechanical_value.tributary_area, 3.0 / 76.0, 1.0e-13),
                 "HEX20 Q8 mechanical contact uses the exact positive-lumped corner area") &&
             check(near(positive_areas[4], 4.0 / 19.0, 1.0e-13) && near(positive_area_sum, 1.0, 1.0e-13),
                 "HEX20 positive-lumped edge areas are exact and sum to the current face area") &&
             check(near(consistent_corner.tributary_area, -1.0 / 12.0, 1.0e-13) &&
                       near(consistent_midpoint.tributary_area, 1.0 / 3.0, 1.0e-13),
                 "HEX20 retains the signed consistent-shape area rule for explicit comparisons") &&
             check(std::abs(force_balance_x) < 1.0e-12 && std::abs(force_balance_y) < 1.0e-12 &&
                       std::abs(force_balance_z) < 1.0e-12,
                 "HEX20 Q8 mechanical contact residual is action-reaction conservative") &&
             passed;
    fuelsim::NormalContactProperties friction_properties{1.0e5, 0.2, false};
    auto sliding_state = state;
    sliding_state[24] = 0.01;
    const auto sliding_value = fuelsim::compute_node_to_quad8_contact_value(
        friction_properties, mechanical_geometry, sliding_state, committed, {});
    passed = check(sliding_value.sliding &&
                       sliding_value.tangential_traction <=
                           friction_properties.friction_coefficient * sliding_value.pressure * (1.0 + 1.0e-12),
                 "HEX20 positive-lumped corner contact enters sliding with a bounded tangential traction") &&
             passed;
    double jacobian_error = 0.0, jacobian_scale = 0.0;
    constexpr double step = 1.0e-7;
    auto plus = state, minus = state;
    plus[8] += step;
    minus[8] -= step;
    const auto plus_residual =
        fuelsim::compute_node_to_quad8_contact(mechanical_properties, mechanical_geometry, plus, committed, {});
    const auto minus_residual =
        fuelsim::compute_node_to_quad8_contact(mechanical_properties, mechanical_geometry, minus, committed, {});
    for (std::size_t row = 0; row < plus_residual.size(); ++row) {
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double analytic = mechanical_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 8];
        jacobian_error = std::max(jacobian_error, std::abs(analytic - numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(jacobian_error / jacobian_scale < 2.0e-5,
                 "HEX20 node-to-surface mechanical-contact Jacobian matches a centered difference") &&
             passed;

    fuelsim::Quad8SurfaceContactLocalResidual surface_residual{}, surface_plus_residual{}, surface_minus_residual{},
        surface_rotate_plus_residual{}, surface_rotate_minus_residual{};
    fuelsim::Quad8SurfaceContactLocalJacobian surface_jacobian{};
    auto rotate_plus = state, rotate_minus = state;
    rotate_plus[17] += step;
    rotate_minus[17] -= step;
    for (const fuelsim::Quad8FaceMechanicalQuadraturePoint& quadrature : geometry.mechanical_points) {
        const fuelsim::Quad8ReferenceProjectionValue reference =
            fuelsim::compute_quad8_reference_projection(face, face, quadrature.displacement_shape, -1.0);
        const fuelsim::Quad8ToQuad8MechanicalGeometry surface_geometry{face, face, quadrature.displacement_shape,
            quadrature.derivative_xi, quadrature.derivative_eta, reference.primary_shape,
            reference.primary_derivative_xi, reference.primary_derivative_eta, quadrature.quadrature_weight, -1.0};
        fuelsim::Quad8SurfaceContactLocalJacobian point_jacobian{};
        const auto point_residual = fuelsim::compute_quad8_to_quad8_contact(
            mechanical_properties, surface_geometry, state, committed, {}, &point_jacobian);
        const auto point_plus =
            fuelsim::compute_quad8_to_quad8_contact(mechanical_properties, surface_geometry, plus, committed, {});
        const auto point_minus =
            fuelsim::compute_quad8_to_quad8_contact(mechanical_properties, surface_geometry, minus, committed, {});
        const auto point_rotate_plus = fuelsim::compute_quad8_to_quad8_contact(
            mechanical_properties, surface_geometry, rotate_plus, committed, {});
        const auto point_rotate_minus = fuelsim::compute_quad8_to_quad8_contact(
            mechanical_properties, surface_geometry, rotate_minus, committed, {});
        for (std::size_t entry = 0; entry < surface_residual.size(); ++entry) {
            surface_residual[entry] += point_residual[entry];
            surface_plus_residual[entry] += point_plus[entry];
            surface_minus_residual[entry] += point_minus[entry];
            surface_rotate_plus_residual[entry] += point_rotate_plus[entry];
            surface_rotate_minus_residual[entry] += point_rotate_minus[entry];
        }
        for (std::size_t entry = 0; entry < surface_jacobian.size(); ++entry)
            surface_jacobian[entry] += point_jacobian[entry];
    }
    force_balance_x = 0.0;
    force_balance_y = 0.0;
    force_balance_z = 0.0;
    for (std::size_t node = 0; node < 16; ++node) {
        force_balance_x += surface_residual[8 + node];
        force_balance_y += surface_residual[24 + node];
        force_balance_z += surface_residual[40 + node];
    }
    const double total_secondary_force =
        std::accumulate(surface_residual.begin() + 8, surface_residual.begin() + 16, 0.0);
    passed = check(near(std::abs(total_secondary_force), 1.0e3, 1.0e-12),
                 "HEX20 surface-to-surface contact integrates constant pressure over the current face") &&
             check(near(surface_residual[8], -total_secondary_force / 12.0, 2.0e-12) &&
                       near(surface_residual[12], total_secondary_force / 3.0, 2.0e-12),
                 "HEX20 surface-to-surface contact retains the exact signed Q8 consistent nodal forces") &&
             check(std::abs(force_balance_x) < 1.0e-12 && std::abs(force_balance_y) < 1.0e-12 &&
                       std::abs(force_balance_z) < 1.0e-12,
                 "HEX20 surface-to-surface contact is exactly action-reaction conservative") &&
             passed;
    auto surface_sliding_state = state;
    for (std::size_t node = 0; node < 8; ++node) surface_sliding_state[24 + node] += 0.01;
    const fuelsim::Quad8FaceMechanicalQuadraturePoint& center = geometry.mechanical_points[4];
    const fuelsim::Quad8ReferenceProjectionValue center_reference =
        fuelsim::compute_quad8_reference_projection(face, face, center.displacement_shape, -1.0);
    const fuelsim::Quad8ToQuad8MechanicalGeometry center_geometry{face, face, center.displacement_shape,
        center.derivative_xi, center.derivative_eta, center_reference.primary_shape,
        center_reference.primary_derivative_xi, center_reference.primary_derivative_eta, center.quadrature_weight,
        -1.0};
    auto finite_sliding_geometry = center_geometry;
    finite_sliding_geometry.finite_sliding = true;
    auto finite_sliding_state = state;
    for (std::size_t node = 0; node < 8; ++node) finite_sliding_state[24 + node] += 0.25;
    fuelsim::Quad8SurfaceContactLocalJacobian finite_sliding_jacobian{};
    const auto finite_sliding_residual = fuelsim::compute_quad8_to_quad8_contact(
        mechanical_properties, finite_sliding_geometry, finite_sliding_state, committed, {}, &finite_sliding_jacobian);
    double primary_force = 0.0, primary_force_y_moment = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        primary_force += finite_sliding_residual[16 + node];
        primary_force_y_moment += face[node].y * finite_sliding_residual[16 + node];
    }
    passed = check(std::abs(primary_force_y_moment / primary_force - 0.75) < 1.0e-12,
                 "HEX20 finite sliding moves the primary test-field projection to the current closest point") &&
             passed;
    auto finite_plus = finite_sliding_state, finite_minus = finite_sliding_state;
    finite_plus[25] += step;
    finite_minus[25] -= step;
    const auto finite_plus_residual = fuelsim::compute_quad8_to_quad8_contact(
        mechanical_properties, finite_sliding_geometry, finite_plus, committed, {});
    const auto finite_minus_residual = fuelsim::compute_quad8_to_quad8_contact(
        mechanical_properties, finite_sliding_geometry, finite_minus, committed, {});
    double finite_jacobian_error = 0.0, finite_jacobian_scale = 0.0;
    for (std::size_t row = 0; row < finite_plus_residual.size(); ++row) {
        const double numerical = (finite_plus_residual[row] - finite_minus_residual[row]) / (2.0 * step);
        const double analytic = finite_sliding_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 25];
        finite_jacobian_error = std::max(finite_jacobian_error, std::abs(analytic - numerical));
        finite_jacobian_scale = std::max({finite_jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(finite_jacobian_error / finite_jacobian_scale < 2.0e-5,
                 "HEX20 finite-sliding closest-point Jacobian matches a centered directional difference") &&
             passed;
    auto outside_state = state;
    for (std::size_t node = 0; node < 8; ++node) outside_state[24 + node] += 0.6;
    passed =
        check(!fuelsim::compute_quad8_to_quad8_contact_projection(finite_sliding_geometry, outside_state).projected,
            "HEX20 finite sliding rejects a closest point outside the complete primary face") &&
        passed;
    const auto surface_sliding = fuelsim::compute_quad8_to_quad8_contact_value(
        friction_properties, center_geometry, surface_sliding_state, committed, {});
    passed = check(surface_sliding.sliding &&
                       near(surface_sliding.tangential_traction,
                           friction_properties.friction_coefficient * surface_sliding.pressure, 1.0e-12),
                 "HEX20 surface-to-surface Coulomb friction caps the integration-point tangential traction") &&
             passed;
    auto biaxial_sliding_state = state;
    for (std::size_t node = 0; node < 8; ++node) {
        biaxial_sliding_state[24 + node] += 0.006;
        biaxial_sliding_state[40 + node] += 0.008;
    }
    const fuelsim::NormalContactProperties controlled_slip_properties{1.0e5, 0.2, false, 4.0e-4};
    const auto biaxial_sliding = fuelsim::compute_quad8_to_quad8_contact_value(
        controlled_slip_properties, center_geometry, biaxial_sliding_state, committed, {});
    const double biaxial_elastic_slip = std::hypot(biaxial_sliding.elastic_tangential_slip[0],
        std::hypot(biaxial_sliding.elastic_tangential_slip[1], biaxial_sliding.elastic_tangential_slip[2]));
    passed = check(biaxial_sliding.sliding && std::abs(biaxial_sliding.tangential_traction_vector[1]) > 0.0 &&
                       std::abs(biaxial_sliding.tangential_traction_vector[2]) > 0.0,
                 "HEX20 surface-to-surface friction enters sliding with two simultaneous tangential components") &&
             check(near(biaxial_sliding.tangential_traction,
                       controlled_slip_properties.friction_coefficient * biaxial_sliding.pressure, 1.0e-12) &&
                       near(biaxial_elastic_slip, controlled_slip_properties.maximum_elastic_slip, 1.0e-12),
                 "HEX20 biaxial sliding satisfies the Coulomb circle and requested maximum elastic slip") &&
             passed;

    constexpr double curved_angle = 0.25;
    const fuelsim::Quad8FaceCoordinates curved_face = {
        fuelsim::CartesianPoint3{std::cos(-curved_angle), std::sin(-curved_angle), 0.0},
        fuelsim::CartesianPoint3{std::cos(curved_angle), std::sin(curved_angle), 0.0},
        fuelsim::CartesianPoint3{std::cos(curved_angle), std::sin(curved_angle), 1.0},
        fuelsim::CartesianPoint3{std::cos(-curved_angle), std::sin(-curved_angle), 1.0},
        fuelsim::CartesianPoint3{1.0, 0.0, 0.0},
        fuelsim::CartesianPoint3{std::cos(curved_angle), std::sin(curved_angle), 0.5},
        fuelsim::CartesianPoint3{1.0, 0.0, 1.0},
        fuelsim::CartesianPoint3{std::cos(-curved_angle), std::sin(-curved_angle), 0.5}};
    const fuelsim::Quad8FaceGeometry curved_geometry = fuelsim::make_quad8_face_geometry(curved_face);
    const fuelsim::Quad8FaceMechanicalQuadraturePoint& curved_center = curved_geometry.mechanical_points[4];
    const fuelsim::Quad8ReferenceProjectionValue curved_reference =
        fuelsim::compute_quad8_reference_projection(curved_face, curved_face, curved_center.displacement_shape, -1.0);
    const fuelsim::Quad8ToQuad8MechanicalGeometry curved_contact_geometry{curved_face, curved_face,
        curved_center.displacement_shape, curved_center.derivative_xi, curved_center.derivative_eta,
        curved_reference.primary_shape, curved_reference.primary_derivative_xi, curved_reference.primary_derivative_eta,
        curved_center.quadrature_weight, -1.0};
    fuelsim::Quad8SurfaceContactLocalValues curved_committed{}, curved_rotated{};
    constexpr double penetration = 0.01, rotation = 0.35;
    constexpr std::array<double, 3> committed_slip = {0.0, 6.0e-4, 8.0e-4};
    const double rotation_cosine = std::cos(rotation), rotation_sine = std::sin(rotation);
    for (std::size_t node = 0; node < 16; ++node) {
        const fuelsim::CartesianPoint3& coordinate = curved_face[node % 8];
        curved_rotated[8 + node] = rotation_cosine * coordinate.x - rotation_sine * coordinate.y - coordinate.x;
        curved_rotated[24 + node] = rotation_sine * coordinate.x + rotation_cosine * coordinate.y - coordinate.y;
        if (node < 8) {
            curved_committed[8 + node] -= penetration;
            curved_committed[24 + node] += committed_slip[1];
            curved_committed[40 + node] += committed_slip[2];
            curved_rotated[8 + node] += -rotation_cosine * penetration - rotation_sine * committed_slip[1];
            curved_rotated[24 + node] += -rotation_sine * penetration + rotation_cosine * committed_slip[1];
            curved_rotated[40 + node] += committed_slip[2];
        }
    }
    fuelsim::ContactPointHistory curved_history{};
    curved_history.cartesian_elastic_tangential_slip = committed_slip;
    const fuelsim::NormalContactProperties curved_friction_properties{1.0e5, 0.5, false};
    fuelsim::Quad8SurfaceContactLocalJacobian curved_contact_jacobian{};
    const auto curved_history_value = fuelsim::compute_quad8_to_quad8_contact_value(
        curved_friction_properties, curved_contact_geometry, curved_rotated, curved_committed, curved_history);
    (void)fuelsim::compute_quad8_to_quad8_contact(curved_friction_properties, curved_contact_geometry, curved_rotated,
        curved_committed, curved_history, &curved_contact_jacobian);
    const std::array<double, 3> rotated_history = {-rotation_sine * 6.0e-4, rotation_cosine * 6.0e-4, 8.0e-4};
    double history_error = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        history_error = std::max(history_error,
            std::abs(curved_history_value.elastic_tangential_slip[component] - rotated_history[component]));
    const double history_normal = curved_history_value.elastic_tangential_slip[0] * curved_history_value.normal[0] +
                                  curved_history_value.elastic_tangential_slip[1] * curved_history_value.normal[1] +
                                  curved_history_value.elastic_tangential_slip[2] * curved_history_value.normal[2];
    passed =
        check(curved_history_value.projected && !curved_history_value.sliding && history_error < 1.0e-13,
            "HEX20 curved surface rotates a nonzero two-component elastic-slip history with its tangent plane") &&
        check(std::abs(history_normal) < 1.0e-13 && near(std::hypot(curved_history_value.elastic_tangential_slip[0],
                                                             std::hypot(curved_history_value.elastic_tangential_slip[1],
                                                                 curved_history_value.elastic_tangential_slip[2])),
                                                        1.0e-3, 1.0e-13),
            "HEX20 curved tangent-plane transport preserves tangency and elastic-slip magnitude") &&
        passed;

    auto curved_plus = curved_rotated, curved_minus = curved_rotated;
    curved_plus[17] += step;
    curved_minus[17] -= step;
    const auto curved_plus_residual = fuelsim::compute_quad8_to_quad8_contact(
        curved_friction_properties, curved_contact_geometry, curved_plus, curved_committed, curved_history);
    const auto curved_minus_residual = fuelsim::compute_quad8_to_quad8_contact(
        curved_friction_properties, curved_contact_geometry, curved_minus, curved_committed, curved_history);
    double curved_jacobian_error = 0.0, curved_jacobian_scale = 0.0;
    for (std::size_t row = 0; row < curved_plus_residual.size(); ++row) {
        const double numerical = (curved_plus_residual[row] - curved_minus_residual[row]) / (2.0 * step);
        const double analytic = curved_contact_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 17];
        curved_jacobian_error = std::max(curved_jacobian_error, std::abs(analytic - numerical));
        curved_jacobian_scale = std::max({curved_jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    passed = check(curved_jacobian_error / curved_jacobian_scale < 2.0e-5,
                 "HEX20 curved history-rotation Jacobian matches a centered directional difference") &&
             passed;
    jacobian_error = 0.0;
    jacobian_scale = 0.0;
    for (std::size_t row = 0; row < surface_residual.size(); ++row) {
        const double numerical = (surface_plus_residual[row] - surface_minus_residual[row]) / (2.0 * step);
        const double analytic = surface_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 8];
        jacobian_error = std::max(jacobian_error, std::abs(analytic - numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(analytic), std::abs(numerical)});
        const double rotating_numerical =
                         (surface_rotate_plus_residual[row] - surface_rotate_minus_residual[row]) / (2.0 * step),
                     rotating_analytic = surface_jacobian[row * fuelsim::quad8_surface_contact_local_dof_count + 17];
        jacobian_error = std::max(jacobian_error, std::abs(rotating_analytic - rotating_numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(rotating_analytic), std::abs(rotating_numerical)});
    }
    return check(jacobian_error / jacobian_scale < 2.0e-5,
               "HEX20 small-sliding contact gap and rotating-normal Jacobian match centered differences") &&
           passed;
}
} // namespace

int main() {
    const bool passed = test_geometry_and_constant_strain() && test_jacobian_and_transient_history() &&
                        test_quadratic_face() && test_warped_geometry() && test_hex20_contact_kernels();
    if (passed) std::cout << "All HEX20-U2/T1 kernel tests passed\n";
    return passed ? 0 : 1;
}
