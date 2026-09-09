#include "fuelsim/elements/rz_geometry.hpp"
#include "fuelsim/elements/detail/rz_point.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss = 0.577350269189625764509148780501957456;

std::array<double, quad4_node_count> shape_functions(double xi, double eta) {
    return {
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    };
}

std::array<double, quad4_node_count> shape_derivative_xi(double eta) {
    return {
        -0.25 * (1.0 - eta),
        0.25 * (1.0 - eta),
        0.25 * (1.0 + eta),
        -0.25 * (1.0 + eta),
    };
}

std::array<double, quad4_node_count> shape_derivative_eta(double xi) {
    return {
        -0.25 * (1.0 - xi),
        -0.25 * (1.0 + xi),
        0.25 * (1.0 + xi),
        0.25 * (1.0 - xi),
    };
}
} // namespace

Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates) {
    const std::array<std::array<double, 2>, 4> locations = {{
        {{-gauss, -gauss}},
        {{gauss, -gauss}},
        {{gauss, gauss}},
        {{-gauss, gauss}},
    }};
    Quad4RzGeometry geometry{};
    geometry.coordinates = coordinates;
    for (std::size_t q = 0; q < locations.size(); ++q) {
        const double xi = locations[q][0], eta = locations[q][1];
        const std::array<double, quad4_node_count> shape = shape_functions(xi, eta);
        const std::array<double, quad4_node_count> derivative_xi = shape_derivative_xi(eta);
        const std::array<double, quad4_node_count> derivative_eta = shape_derivative_eta(xi);
        double dr_dxi = 0.0, dr_deta = 0.0, dz_dxi = 0.0, dz_deta = 0.0, radius = 0.0, axial_coordinate = 0.0;
        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            dr_dxi += derivative_xi[node] * coordinates[node].r;
            dr_deta += derivative_eta[node] * coordinates[node].r;
            dz_dxi += derivative_xi[node] * coordinates[node].z;
            dz_deta += derivative_eta[node] * coordinates[node].z;
            radius += shape[node] * coordinates[node].r;
            axial_coordinate += shape[node] * coordinates[node].z;
        }
        const double determinant = dr_dxi * dz_deta - dr_deta * dz_dxi;
        if (!(determinant > 0.0))
            throw std::invalid_argument("Quad4RzGeometry requires positive Jacobian determinant");
        if (!(radius > 0.0))
            throw std::invalid_argument("Quad4RzGeometry requires positive quadrature radius");
        RzQuadraturePoint& point = geometry.points[q];
        point.shape = shape;
        point.radius = radius;
        point.axial_coordinate = axial_coordinate;
        point.weighted_measure = 2.0 * pi * radius * determinant;
        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            point.gradient_r[node] = (dz_deta * derivative_xi[node] - dz_dxi * derivative_eta[node]) / determinant;
            point.gradient_z[node] = (-dr_deta * derivative_xi[node] + dr_dxi * derivative_eta[node]) / determinant;
        }
    }
    return geometry;
}

// Kinematics core driven by the quadrature-point radial displacement and the four in-plane
// displacement-gradient components. These point quantities may carry any ADlite seeding
// (passive, or independent variables), which lets callers choose the derivative width.
AxisymmetricKinematics evaluate_axisymmetric_kinematics_from_point(const RzQuadraturePoint& point,
    const adlite::Scalar& radial_displacement,
    const adlite::Scalar& displacement_gradient_rr,
    const adlite::Scalar& displacement_gradient_rz,
    const adlite::Scalar& displacement_gradient_zr,
    const adlite::Scalar& displacement_gradient_zz,
    const LocalValues& committed_state,
    StrainFormulation strain_formulation,
    bool cax4t) {
    AxisymmetricKinematics result{};
    if (strain_formulation == StrainFormulation::small) {
        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            result.gradient_r[node] = point.gradient_r[node];
            result.gradient_z[node] = point.gradient_z[node];
        }
        result.radius = point.radius;
        result.weighted_measure = point.weighted_measure;
        result.midpoint_weighted_measure = point.weighted_measure;
        result.strain_rr = displacement_gradient_rr;
        result.strain_zz = displacement_gradient_zz;
        result.strain_hoop = radial_displacement / point.radius;
        result.strain_rz = 0.5 * (displacement_gradient_rz + displacement_gradient_zr);
        return result;
    }
    const adlite::Scalar deformation_rr = 1.0 + displacement_gradient_rr, deformation_rz = displacement_gradient_rz,
                         deformation_zr = displacement_gradient_zr, deformation_zz = 1.0 + displacement_gradient_zz,
                         deformation_hoop = 1.0 + radial_displacement / point.radius,
                         determinant_rz = deformation_rr * deformation_zz - deformation_rz * deformation_zr,
                         current_radius = point.radius + radial_displacement;
    if (!std::isfinite(determinant_rz.value()) || !(determinant_rz.value() > 0.0))
        throw std::domain_error("Finite-strain Quad4 RZ deformation must preserve a positive in-plane Jacobian");
    if (!std::isfinite(deformation_hoop.value()) || !(deformation_hoop.value() > 0.0)
        || !std::isfinite(current_radius.value()) || !(current_radius.value() > 0.0))
        throw std::domain_error("Finite-strain RZ deformation requires positive hoop stretch and radius");
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        result.gradient_r[node] =
            (deformation_zz * point.gradient_r[node] - deformation_zr * point.gradient_z[node]) / determinant_rz;
        result.gradient_z[node] =
            (-deformation_rz * point.gradient_r[node] + deformation_rr * point.gradient_z[node]) / determinant_rz;
    }
    result.radius = current_radius;
    result.weighted_measure = point.weighted_measure * determinant_rz * deformation_hoop;
    const double old_radial_displacement = quad4_rz_detail::interpolate(point.shape, committed_state, 4),
                 old_deformation_rr = 1.0 + quad4_rz_detail::interpolate(point.gradient_r, committed_state, 4),
                 old_deformation_rz = quad4_rz_detail::interpolate(point.gradient_z, committed_state, 4),
                 old_deformation_zr = quad4_rz_detail::interpolate(point.gradient_r, committed_state, 8),
                 old_deformation_zz = 1.0 + quad4_rz_detail::interpolate(point.gradient_z, committed_state, 8),
                 old_deformation_hoop = 1.0 + old_radial_displacement / point.radius,
                 old_determinant_rz = old_deformation_rr * old_deformation_zz - old_deformation_rz * old_deformation_zr;
    if (!std::isfinite(old_determinant_rz) || !(old_determinant_rz > 0.0) || !std::isfinite(old_deformation_hoop)
        || !(old_deformation_hoop > 0.0))
        throw std::domain_error("Committed finite-strain RZ state requires positive Jacobian and hoop stretch");
    const double old_inverse_rr = old_deformation_zz / old_determinant_rz,
                 old_inverse_rz = -old_deformation_rz / old_determinant_rz,
                 old_inverse_zr = -old_deformation_zr / old_determinant_rz,
                 old_inverse_zz = old_deformation_rr / old_determinant_rz;
    const adlite::Scalar incremental_rr = deformation_rr * old_inverse_rr + deformation_rz * old_inverse_zr,
                         incremental_rz = deformation_rr * old_inverse_rz + deformation_rz * old_inverse_zz,
                         incremental_zr = deformation_zr * old_inverse_rr + deformation_zz * old_inverse_zr,
                         incremental_zz = deformation_zr * old_inverse_rz + deformation_zz * old_inverse_zz,
                         incremental_hoop = deformation_hoop / old_deformation_hoop,
                         incremental_determinant = incremental_rr * incremental_zz - incremental_rz * incremental_zr;
    if (!std::isfinite(incremental_determinant.value()) || !(incremental_determinant.value() > 0.0)
        || !std::isfinite(incremental_hoop.value()) || !(incremental_hoop.value() > 0.0))
        throw std::domain_error("Incremental finite-strain RZ state requires positive Jacobian and hoop stretch");
    if (cax4t) {
        const adlite::Scalar sum_rr = deformation_rr + old_deformation_rr, sum_rz = deformation_rz + old_deformation_rz,
                             sum_zr = deformation_zr + old_deformation_zr, sum_zz = deformation_zz + old_deformation_zz,
                             sum_hoop = deformation_hoop + old_deformation_hoop,
                             determinant_sum = sum_rr * sum_zz - sum_rz * sum_zr;
        if (!(determinant_sum.value() > 0.0) || !(sum_hoop.value() > 0.0))
            throw std::domain_error("CAX4T midpoint configuration must preserve positive volume");
        result.midpoint_weighted_measure = point.weighted_measure * determinant_sum * sum_hoop / 8.0;
        const adlite::Scalar
            hrr =
                2.0 * ((deformation_rr - old_deformation_rr) * sum_zz - (deformation_rz - old_deformation_rz) * sum_zr)
                / determinant_sum,
            hrz = 2.0
                  * (-(deformation_rr - old_deformation_rr) * sum_rz + (deformation_rz - old_deformation_rz) * sum_rr)
                  / determinant_sum,
            hzr = 2.0
                  * ((deformation_zr - old_deformation_zr) * sum_zz - (deformation_zz - old_deformation_zz) * sum_zr)
                  / determinant_sum,
            hzz = 2.0
                  * (-(deformation_zr - old_deformation_zr) * sum_rz + (deformation_zz - old_deformation_zz) * sum_rr)
                  / determinant_sum,
            shear = 0.5 * (hrz + hzr), spin = 0.25 * (hrz - hzr), denom = 1.0 + spin * spin,
            cosine = (1.0 - spin * spin) / denom, sine = 2.0 * spin / denom;
        result.rotation = {cosine, sine, -sine, cosine, adlite::Scalar(1.0)};
        result.strain_rr = cosine * cosine * hrr - 2.0 * cosine * sine * shear + sine * sine * hzz;
        result.strain_zz = sine * sine * hrr + 2.0 * cosine * sine * shear + cosine * cosine * hzz;
        result.strain_rz = cosine * sine * (hrr - hzz) + (cosine * cosine - sine * sine) * shear;
        result.strain_hoop = 2.0 * (deformation_hoop - old_deformation_hoop) / sum_hoop;
        return result;
    }
    const adlite::Scalar inverse_rr = incremental_zz / incremental_determinant,
                         inverse_rz = -incremental_rz / incremental_determinant,
                         inverse_zr = -incremental_zr / incremental_determinant,
                         inverse_zz = incremental_rr / incremental_determinant, inverse_hoop = 1.0 / incremental_hoop,
                         cinv_rr = inverse_rr * inverse_rr + inverse_rz * inverse_rz - 1.0,
                         cinv_zz = inverse_zr * inverse_zr + inverse_zz * inverse_zz - 1.0,
                         cinv_rz = inverse_rr * inverse_zr + inverse_rz * inverse_zz,
                         cinv_hoop = inverse_hoop * inverse_hoop - 1.0;
    result.strain_rr = -0.5 * cinv_rr + 0.25 * (cinv_rr * cinv_rr + cinv_rz * cinv_rz);
    result.strain_zz = -0.5 * cinv_zz + 0.25 * (cinv_rz * cinv_rz + cinv_zz * cinv_zz);
    result.strain_rz = -0.5 * cinv_rz + 0.25 * cinv_rz * (cinv_rr + cinv_zz);
    result.strain_hoop = -0.5 * cinv_hoop + 0.25 * cinv_hoop * cinv_hoop;
    const adlite::Scalar axial_rotation = inverse_rz - inverse_zr, q = 0.25 * axial_rotation * axial_rotation,
                         trace_minus_one = inverse_rr + inverse_zz + inverse_hoop - 1.0,
                         p = 0.25 * trace_minus_one * trace_minus_one, sum = p + q;
    if (!std::isfinite(sum.value()) || !(sum.value() > 0.0))
        throw std::domain_error("MOOSE Taylor finite-strain rotation has invalid p+q");
    const adlite::Scalar p2 = p * p, p3 = p2 * p, sum2 = sum * sum, sum3 = sum2 * sum,
                         c1_squared = p + 3.0 * p2 * (1.0 - sum) / sum2 - 2.0 * p3 * (1.0 - sum) / sum3;
    if (!std::isfinite(c1_squared.value()) || !(c1_squared.value() > 0.0))
        throw std::domain_error("MOOSE Taylor finite-strain rotation has nonpositive C1 squared");
    const adlite::Scalar c1 = adlite::sqrt(c1_squared);
    adlite::Scalar c2;
    if (q.value() > 0.01) {
        c2 = (1.0 - c1) / (4.0 * q);
    } else {
        const adlite::Scalar q2 = q * q, q3 = q2 * q, p4 = p3 * p;
        c2 = 0.125 + q * 0.03125 * (p2 - 12.0 * (p - 1.0)) / p2 + q2 * (p - 2.0) * (p2 - 10.0 * p + 32.0) / p3
             + q3 * (1104.0 - 992.0 * p + 376.0 * p2 - 72.0 * p3 + 5.0 * p4) / (512.0 * p4);
    }
    const adlite::Scalar c3_test = (p * q * (3.0 - q) + p3 + q * q) / sum3;
    if (!std::isfinite(c3_test.value()) || !(c3_test.value() > 0.0))
        throw std::domain_error("MOOSE Taylor finite-strain rotation has nonpositive C3 test");
    const adlite::Scalar c3 = 0.5 * adlite::sqrt(c3_test);
    result.rotation.rr = c1;
    result.rotation.rz = -c3 * axial_rotation;
    result.rotation.zr = c3 * axial_rotation;
    result.rotation.zz = c1;
    result.rotation.hoop = c1 + c2 * axial_rotation * axial_rotation;
    return result;
}

// Evaluates the constitutive relation with AD seeded only on the four strain components and
// the temperature (width 5), returning the stress values, the consistent material tangent
// d(stress)/d(strain), and the thermal coupling d(stress)/dT.
AxisymmetricStressTangent evaluate_axisymmetric_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 4>& fed_strain,
    double temperature,
    double time_step,
    const MaterialPointState* committed_material,
    MaterialFunctionContext context) {
    const std::array<double, 5> seeds = {fed_strain[0], fed_strain[1], fed_strain[2], fed_strain[3], temperature};
    std::array<adlite::Scalar, 5> active{};
    adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    const AxisymmetricStress stress =
        committed_material == nullptr ? material.stress(active[0], active[1], active[2], active[3], active[4], context)
                                      : material
                                            .response(active[0],
                                                active[1],
                                                active[2],
                                                active[3],
                                                active[4],
                                                time_step,
                                                *committed_material,
                                                context)
                                            .stress;
    const std::array<const adlite::Scalar*, 4> components = {&stress.rr, &stress.zz, &stress.hoop, &stress.rz};
    AxisymmetricStressTangent result{};
    result.stress = {stress.rr.value(), stress.zz.value(), stress.hoop.value(), stress.rz.value()};
    std::array<double, 5> derivatives{};
    for (std::size_t row = 0; row < 4; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 4; ++column)
            result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[4];
    }
    return result;
}

AxisymmetricKinematics evaluate_axisymmetric_kinematics(const RzQuadraturePoint& point,
    const LocalAdValues& state,
    StrainFormulation strain_formulation) {
    const LocalValues undeformed{};
    return evaluate_axisymmetric_incremental_kinematics(point, state, undeformed, strain_formulation);
}

AxisymmetricKinematics evaluate_axisymmetric_incremental_kinematics(const RzQuadraturePoint& point,
    const LocalAdValues& current_state,
    const LocalValues& committed_state,
    StrainFormulation strain_formulation) {
    return evaluate_axisymmetric_kinematics_from_point(point,
        quad4_rz_detail::interpolate(point.shape, current_state, 4),
        quad4_rz_detail::interpolate(point.gradient_r, current_state, 4),
        quad4_rz_detail::interpolate(point.gradient_z, current_state, 4),
        quad4_rz_detail::interpolate(point.gradient_r, current_state, 8),
        quad4_rz_detail::interpolate(point.gradient_z, current_state, 8),
        committed_state,
        strain_formulation);
}

} // namespace fuelsim
