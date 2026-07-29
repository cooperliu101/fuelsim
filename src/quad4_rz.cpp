#include "fuelsim/quad4_rz.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

std::array<double, 4> shape_functions(double xi, double eta) {
    return {
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    };
}

std::array<double, 4> shape_derivative_xi(double eta) {
    return {
        -0.25 * (1.0 - eta),
        0.25 * (1.0 - eta),
        0.25 * (1.0 + eta),
        -0.25 * (1.0 + eta),
    };
}

std::array<double, 4> shape_derivative_eta(double xi) {
    return {
        -0.25 * (1.0 - xi),
        -0.25 * (1.0 + xi),
        0.25 * (1.0 + xi),
        0.25 * (1.0 - xi),
    };
}

adlite::Scalar interpolate(const std::array<double, 4>& coefficients,
                           const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 4; ++node)
        result += coefficients[node] * state[offset + node];
    return result;
}

AxisymmetricStress
stress_at_point(const RzQuadraturePoint& point, const LocalAdValues& state,
                const IsotropicThermoelasticMaterial& material) {
    const adlite::Scalar temperature = interpolate(point.shape, state, 0);
    const adlite::Scalar radial_displacement =
        interpolate(point.shape, state, 4);
    const adlite::Scalar strain_rr = interpolate(point.gradient_r, state, 4);
    const adlite::Scalar strain_zz = interpolate(point.gradient_z, state, 8);
    const adlite::Scalar strain_hoop = radial_displacement / point.radius;
    const adlite::Scalar strain_rz =
        0.5 * (interpolate(point.gradient_z, state, 4) +
               interpolate(point.gradient_r, state, 8));

    return material.stress(strain_rr, strain_zz, strain_hoop, strain_rz,
                           temperature);
}

} // namespace

Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates) {
    constexpr double gp = 0.577350269189625764509148780501957456;
    const std::array<std::array<double, 2>, 4> locations = {{
        {{-gp, -gp}},
        {{gp, -gp}},
        {{gp, gp}},
        {{-gp, gp}},
    }};

    Quad4RzGeometry geometry{};
    for (std::size_t q = 0; q < locations.size(); ++q) {
        const double xi = locations[q][0];
        const double eta = locations[q][1];
        const std::array<double, 4> shape = shape_functions(xi, eta);
        const std::array<double, 4> derivative_xi = shape_derivative_xi(eta);
        const std::array<double, 4> derivative_eta = shape_derivative_eta(xi);

        double dr_dxi = 0.0;
        double dr_deta = 0.0;
        double dz_dxi = 0.0;
        double dz_deta = 0.0;
        double radius = 0.0;
        for (std::size_t node = 0; node < 4; ++node) {
            dr_dxi += derivative_xi[node] * coordinates[node].r;
            dr_deta += derivative_eta[node] * coordinates[node].r;
            dz_dxi += derivative_xi[node] * coordinates[node].z;
            dz_deta += derivative_eta[node] * coordinates[node].z;
            radius += shape[node] * coordinates[node].r;
        }

        const double determinant = dr_dxi * dz_deta - dr_deta * dz_dxi;
        if (!(determinant > 0.0))
            throw std::invalid_argument(
                "Quad4RzGeometry requires positive Jacobian determinant");
        if (!(radius > 0.0))
            throw std::invalid_argument(
                "Quad4RzGeometry requires positive quadrature radius");

        RzQuadraturePoint& point = geometry.points[q];
        point.shape = shape;
        point.radius = radius;
        point.weighted_measure = 2.0 * pi * radius * determinant;

        for (std::size_t node = 0; node < 4; ++node) {
            point.gradient_r[node] = (dz_deta * derivative_xi[node] -
                                      dz_dxi * derivative_eta[node]) /
                                     determinant;
            point.gradient_z[node] = (-dr_deta * derivative_xi[node] +
                                      dr_dxi * derivative_eta[node]) /
                                     determinant;
        }
    }

    return geometry;
}

Quad4RzThermoelasticKernel::Quad4RzThermoelasticKernel(
    IsotropicThermoelasticMaterial material, double volumetric_heat_source)
    : _material(material), _volumetric_heat_source(volumetric_heat_source) {}

LocalResidual
Quad4RzThermoelasticKernel::residual(const Quad4RzGeometry& geometry,
                                     const LocalValues& state) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];

    LocalAdValues passive_residual{};
    residual_ad(geometry, passive_state, passive_residual);

    LocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = passive_residual[row].value();
    return result;
}

LocalSystem
Quad4RzThermoelasticKernel::linearize(const Quad4RzGeometry& geometry,
                                      const LocalValues& state) const {
    LocalAdValues active_state{};
    adlite::seed_identity(state.data(), state.size(), active_state.data());

    LocalAdValues active_residual{};
    residual_ad(geometry, active_state, active_residual);

    LocalSystem result{};
    adlite::extract_jacobian(active_residual.data(), active_residual.size(),
                             active_state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

std::array<AxisymmetricStressValues, 4>
Quad4RzThermoelasticKernel::stress_values(const Quad4RzGeometry& geometry,
                                          const LocalValues& state) const {
    LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        passive_state[dof] = state[dof];

    std::array<AxisymmetricStressValues, 4> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const AxisymmetricStress stress =
            stress_at_point(geometry.points[q], passive_state, _material);
        result[q] = {stress.rr.value(), stress.zz.value(), stress.hoop.value(),
                     stress.rz.value()};
    }
    return result;
}

void Quad4RzThermoelasticKernel::residual_ad(const Quad4RzGeometry& geometry,
                                             const LocalAdValues& state,
                                             LocalAdValues& residual) const {
    std::fill(residual.begin(), residual.end(), adlite::Scalar(0.0));

    for (const RzQuadraturePoint& point : geometry.points) {
        const adlite::Scalar temperature = interpolate(point.shape, state, 0);
        const adlite::Scalar gradient_temperature_r =
            interpolate(point.gradient_r, state, 0);
        const adlite::Scalar gradient_temperature_z =
            interpolate(point.gradient_z, state, 0);
        const adlite::Scalar conductivity = _material.conductivity(temperature);
        const AxisymmetricStress stress =
            stress_at_point(point, state, _material);

        for (std::size_t node = 0; node < 4; ++node) {
            residual[node] +=
                point.weighted_measure *
                (conductivity *
                     (point.gradient_r[node] * gradient_temperature_r +
                      point.gradient_z[node] * gradient_temperature_z) -
                 _volumetric_heat_source * point.shape[node]);

            residual[4 + node] +=
                point.weighted_measure *
                (stress.rr * point.gradient_r[node] +
                 stress.hoop * point.shape[node] / point.radius +
                 stress.rz * point.gradient_z[node]);

            residual[8 + node] +=
                point.weighted_measure * (stress.zz * point.gradient_z[node] +
                                          stress.rz * point.gradient_r[node]);
        }
    }
}

} // namespace fuelsim
