#include "fuelsim/core/cartesian3d_hex20.hpp"
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
} // namespace

int main() {
    const bool passed =
        test_geometry_and_constant_strain() && test_jacobian_and_transient_history() && test_quadratic_face();
    if (passed) std::cout << "All HEX20-U2/T1 kernel tests passed\n";
    return passed ? 0 : 1;
}
