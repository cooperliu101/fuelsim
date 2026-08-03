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
                const IsotropicThermoelasticMaterial& material,
                StrainFormulation strain_formulation) {
    const adlite::Scalar temperature = interpolate(point.shape, state, 0);
    const AxisymmetricKinematics kinematics =
        evaluate_axisymmetric_kinematics(point, state, strain_formulation);
    return material.stress(
        kinematics.strain_rr, kinematics.strain_zz,
        kinematics.strain_hoop, kinematics.strain_rz, temperature);
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

AxisymmetricKinematics evaluate_axisymmetric_kinematics(
    const RzQuadraturePoint& point, const LocalAdValues& state,
    StrainFormulation strain_formulation) {
    const adlite::Scalar radial_displacement =
        interpolate(point.shape, state, 4);
    const adlite::Scalar displacement_gradient_rr =
        interpolate(point.gradient_r, state, 4);
    const adlite::Scalar displacement_gradient_rz =
        interpolate(point.gradient_z, state, 4);
    const adlite::Scalar displacement_gradient_zr =
        interpolate(point.gradient_r, state, 8);
    const adlite::Scalar displacement_gradient_zz =
        interpolate(point.gradient_z, state, 8);

    AxisymmetricKinematics result{};
    if (strain_formulation == StrainFormulation::small) {
        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            result.gradient_r[node] = point.gradient_r[node];
            result.gradient_z[node] = point.gradient_z[node];
        }
        result.radius = point.radius;
        result.weighted_measure = point.weighted_measure;
        result.strain_rr = displacement_gradient_rr;
        result.strain_zz = displacement_gradient_zz;
        result.strain_hoop = radial_displacement / point.radius;
        result.strain_rz =
            0.5 * (displacement_gradient_rz + displacement_gradient_zr);
        return result;
    }

    const adlite::Scalar deformation_rr = 1.0 + displacement_gradient_rr;
    const adlite::Scalar deformation_rz = displacement_gradient_rz;
    const adlite::Scalar deformation_zr = displacement_gradient_zr;
    const adlite::Scalar deformation_zz = 1.0 + displacement_gradient_zz;
    const adlite::Scalar deformation_hoop =
        1.0 + radial_displacement / point.radius;
    const adlite::Scalar determinant_rz =
        deformation_rr * deformation_zz -
        deformation_rz * deformation_zr;
    const adlite::Scalar current_radius = point.radius + radial_displacement;
    if (!std::isfinite(determinant_rz.value()) ||
        !(determinant_rz.value() > 0.0) ||
        !std::isfinite(deformation_hoop.value()) ||
        !(deformation_hoop.value() > 0.0) ||
        !std::isfinite(current_radius.value()) ||
        !(current_radius.value() > 0.0))
        throw std::domain_error(
            "Finite-strain Quad4 RZ deformation must preserve positive "
            "Jacobian and radius");

    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        result.gradient_r[node] =
            (deformation_zz * point.gradient_r[node] -
             deformation_zr * point.gradient_z[node]) /
            determinant_rz;
        result.gradient_z[node] =
            (-deformation_rz * point.gradient_r[node] +
             deformation_rr * point.gradient_z[node]) /
            determinant_rz;
    }
    result.radius = current_radius;
    result.weighted_measure =
        point.weighted_measure * determinant_rz * deformation_hoop;

    const adlite::Scalar left_cauchy_rr =
        deformation_rr * deformation_rr + deformation_rz * deformation_rz;
    const adlite::Scalar left_cauchy_zz =
        deformation_zr * deformation_zr + deformation_zz * deformation_zz;
    const adlite::Scalar left_cauchy_rz =
        deformation_rr * deformation_zr + deformation_rz * deformation_zz;
    const adlite::Scalar mean =
        0.5 * (left_cauchy_rr + left_cauchy_zz);
    const adlite::Scalar radius = adlite::hypot(
        0.5 * (left_cauchy_rr - left_cauchy_zz), left_cauchy_rz);
    const adlite::Scalar eigenvalue_high = mean + radius;
    const adlite::Scalar eigenvalue_low = mean - radius;
    if (!std::isfinite(eigenvalue_low.value()) ||
        !(eigenvalue_low.value() > 0.0) ||
        !std::isfinite(eigenvalue_high.value()) ||
        !(eigenvalue_high.value() > 0.0))
        throw std::domain_error(
            "Finite-strain Quad4 RZ left Cauchy-Green tensor must be "
            "positive definite");

    const adlite::Scalar logarithmic_mean =
        0.25 * (adlite::log(eigenvalue_high) +
                adlite::log(eigenvalue_low));
    adlite::Scalar logarithmic_slope;
    if (radius.value() <= 1.0e-10 * mean.value())
        logarithmic_slope = 0.5 / mean;
    else
        logarithmic_slope =
            adlite::atanh(radius / mean) / (2.0 * radius);

    result.strain_rr =
        logarithmic_mean +
        logarithmic_slope * (left_cauchy_rr - mean);
    result.strain_zz =
        logarithmic_mean +
        logarithmic_slope * (left_cauchy_zz - mean);
    result.strain_rz = logarithmic_slope * left_cauchy_rz;
    result.strain_hoop = adlite::log(deformation_hoop);
    return result;
}

Quad4RzThermoelasticKernel::Quad4RzThermoelasticKernel(
    IsotropicThermoelasticMaterial material, double volumetric_heat_source,
    StrainFormulation strain_formulation)
    : _material(material), _volumetric_heat_source(0.0),
      _strain_formulation(strain_formulation) {
    set_volumetric_heat_source(volumetric_heat_source);
}

double Quad4RzThermoelasticKernel::volumetric_heat_source() const noexcept {
    return _volumetric_heat_source;
}

void Quad4RzThermoelasticKernel::set_volumetric_heat_source(
    double volumetric_heat_source) {
    if (!std::isfinite(volumetric_heat_source) ||
        !(volumetric_heat_source >= 0.0))
        throw std::invalid_argument(
            "Quad4RzThermoelasticKernel volumetric heat source must be finite "
            "and nonnegative");
    _volumetric_heat_source = volumetric_heat_source;
}

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
            stress_at_point(geometry.points[q], passive_state, _material,
                            _strain_formulation);
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
        const AxisymmetricKinematics kinematics =
            evaluate_axisymmetric_kinematics(point, state,
                                             _strain_formulation);
        const AxisymmetricStress stress = _material.stress(
            kinematics.strain_rr, kinematics.strain_zz,
            kinematics.strain_hoop, kinematics.strain_rz, temperature);

        for (std::size_t node = 0; node < 4; ++node) {
            residual[node] +=
                point.weighted_measure *
                (conductivity *
                     (point.gradient_r[node] * gradient_temperature_r +
                      point.gradient_z[node] * gradient_temperature_z) -
                 _volumetric_heat_source * point.shape[node]);

            residual[4 + node] +=
                kinematics.weighted_measure *
                (stress.rr * kinematics.gradient_r[node] +
                 stress.hoop * point.shape[node] / kinematics.radius +
                 stress.rz * kinematics.gradient_z[node]);

            residual[8 + node] +=
                kinematics.weighted_measure *
                (stress.zz * kinematics.gradient_z[node] +
                 stress.rz * kinematics.gradient_r[node]);
        }
    }
}

} // namespace fuelsim
