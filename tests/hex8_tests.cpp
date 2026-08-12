#include "fuelsim/hex8_thermoelastic.hpp"

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

bool near(double actual, double expected, double tolerance) { return std::abs(actual - expected) <= tolerance * std::max({1.0, std::abs(actual), std::abs(expected)}); }

fuelsim::Hex8Coordinates unit_cube() { return {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}}}; }

fuelsim::ThermoelasticProperties properties() { return {3000.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, -1.0e8, 0.0, 1.0e-8}; }

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
            for (std::size_t component = 0; component < 3; ++component) gradient_sum[component] += point.gradient[node][component];
        }
        if (!check(near(shape_sum, 1.0, 1.0e-14), "HEX8 shape functions form a partition of unity") ||
            !check(near(gradient_sum[0], 0.0, 1.0e-14) && near(gradient_sum[1], 0.0, 1.0e-14) && near(gradient_sum[2], 0.0, 1.0e-14), "HEX8 shape gradients sum to zero"))
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
    fuelsim::Hex8ThermoelasticKernel kernel(fuelsim::IsotropicThermoelasticMaterial(properties()), 0.0);
    const auto stresses = kernel.stress_values(geometry, state);
    const double lambda = 2.0e11 * 0.25 / (1.25 * 0.5);
    const double shear = 2.0e11 / 2.5;
    const double trace = exx + eyy + ezz;
    for (const fuelsim::SymmetricTensor3Values& stress : stresses) {
        if (!check(near(stress.xx, lambda * trace + 2.0 * shear * exx, 2.0e-13) && near(stress.yy, lambda * trace + 2.0 * shear * eyy, 2.0e-13) && near(stress.zz, lambda * trace + 2.0 * shear * ezz, 2.0e-13) && near(stress.xy, 2.0 * shear * exy, 2.0e-13) &&
                       near(stress.yz, 2.0 * shear * eyz, 2.0e-13) && near(stress.xz, 2.0 * shear * exz, 2.0e-13),
                   "HEX8 reproduces all six constant-strain stress components"))
            return false;
    }
    return true;
}

bool test_free_thermal_expansion_and_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
    fuelsim::Hex8ThermoelasticKernel kernel(fuelsim::IsotropicThermoelasticMaterial(properties()), 7.0e5);
    fuelsim::Hex8LocalValues state{};
    const double temperature = 650.0;
    const double active_alpha = properties().thermal_expansion + properties().thermal_expansion_temperature_coefficient * (temperature - properties().reference_temperature);
    const double strain = active_alpha * (temperature - properties().reference_temperature);
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = temperature;
        state[8 + node] = strain * coordinates[node].x;
        state[16 + node] = strain * coordinates[node].y;
        state[24 + node] = strain * coordinates[node].z;
    }
    for (const fuelsim::SymmetricTensor3Values& stress : kernel.stress_values(geometry, state)) {
        if (!check(std::max({std::abs(stress.xx), std::abs(stress.yy), std::abs(stress.zz), std::abs(stress.xy), std::abs(stress.yz), std::abs(stress.xz)}) < 1.0e-4, "uniform three-dimensional thermal expansion is stress free")) return false;
    }

    for (std::size_t dof = 0; dof < state.size(); ++dof) state[dof] += dof < 8 ? 2.0 * static_cast<double>(dof) : 1.0e-5 * static_cast<double>(dof + 1);
    std::array<double, 32> direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof) direction[dof] = std::sin(0.37 * static_cast<double>(dof + 1));
    const fuelsim::Hex8LocalSystem system = kernel.linearize(geometry, state);
    const double epsilon = 1.0e-7;
    fuelsim::Hex8LocalValues plus = state;
    fuelsim::Hex8LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += epsilon * direction[dof];
        minus[dof] -= epsilon * direction[dof];
    }
    const fuelsim::Hex8LocalResidual plus_residual = kernel.residual(geometry, plus);
    const fuelsim::Hex8LocalResidual minus_residual = kernel.residual(geometry, minus);
    double maximum_error = 0.0;
    double scale = 0.0;
    for (std::size_t row = 0; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column) analytic += system.jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * epsilon);
        maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    return check(maximum_error / scale < 3.0e-7, "full 32-DOF HEX8 automatic-differentiation Jacobian matches a centered directional difference");
}

bool test_transient_capacity_and_faces() {
    const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(unit_cube());
    fuelsim::Hex8ThermoelasticKernel kernel(fuelsim::IsotropicThermoelasticMaterial(properties()), 6.0e6, 1.2e7);
    fuelsim::Hex8LocalValues old_state{};
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        old_state[node] = 300.0;
        state[node] = 302.0;
    }
    const fuelsim::Hex8LocalResidual residual = kernel.residual(geometry, state, old_state, 1.0);
    for (std::size_t node = 0; node < 8; ++node) {
        if (!check(std::abs(residual[node]) < 1.0e-7, "Backward Euler consistent heat capacity balances uniform volumetric heating; residual=" + std::to_string(residual[node]))) return false;
    }

    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Quad4FaceCoordinates face_coordinates = {{coordinates[1], coordinates[2], coordinates[6], coordinates[5]}};
    const fuelsim::Quad4FaceGeometry face = fuelsim::make_quad4_face_geometry(face_coordinates);
    fuelsim::Quad4FaceLocalValues face_state{};
    fuelsim::Quad4FacePressureKernel pressure(5.0);
    const fuelsim::Quad4FaceLocalResidual pressure_residual = pressure.residual(face, face_state);
    double force_x = 0.0;
    double force_y = 0.0;
    double force_z = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        force_x += pressure_residual[4 + node];
        force_y += pressure_residual[8 + node];
        force_z += pressure_residual[12 + node];
    }
    if (!check(near(force_x, 5.0, 1.0e-14) && near(force_y, 0.0, 1.0e-14) && near(force_z, 0.0, 1.0e-14), "reference pressure uses the outward three-dimensional face area vector and exact total force")) return false;

    for (std::size_t node = 0; node < 4; ++node) face_state[node] = 350.0;
    fuelsim::Quad4FaceConvectionKernel convection(20.0, 300.0);
    const fuelsim::Quad4FaceLocalSystem convection_system = convection.linearize(face, face_state);
    double heat = 0.0;
    double tangent_sum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        heat += convection_system.residual[row];
        for (std::size_t column = 0; column < 4; ++column) tangent_sum += convection_system.jacobian[row * 16 + column];
    }
    return check(near(heat, 1000.0, 1.0e-14) && near(tangent_sum, 20.0, 1.0e-14), "three-dimensional convection has the exact face heat rate and consistent temperature tangent");
}

} // namespace

int main() {
    bool passed = true;
    passed = test_geometry_and_constant_strain() && passed;
    passed = test_free_thermal_expansion_and_jacobian() && passed;
    passed = test_transient_capacity_and_faces() && passed;
    if (!passed) return 1;
    std::cout << "HEX8 thermo-mechanics tests passed\n";
    return 0;
}
