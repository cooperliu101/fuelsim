#include "ad_local_system.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/quad4_rz.hpp"
#include <adlite/adlite.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
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
        if (!(determinant > 0.0)) throw std::invalid_argument("Quad4RzGeometry requires positive Jacobian determinant");
        if (!(radius > 0.0)) throw std::invalid_argument("Quad4RzGeometry requires positive quadrature radius");
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
namespace quad4_rz_detail {
inline adlite::Scalar interpolate(
    const std::array<double, quad4_node_count>& coefficients, const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node) result += coefficients[node] * state[offset + node];
    return result;
}
inline double interpolate(
    const std::array<double, quad4_node_count>& coefficients, const LocalValues& state, std::size_t offset) {
    double result = 0.0;
    for (std::size_t node = 0; node < quad4_node_count; ++node) result += coefficients[node] * state[offset + node];
    return result;
}
inline LocalAdValues ad_state(const LocalValues& state, bool active = false) {
    LocalAdValues result{};
    if (active)
        ad_local_system::make_active(state.data(), state.size(), result.data());
    else
        ad_local_system::make_passive(state.data(), state.size(), result.data());
    return result;
}
inline LocalResidual values(
    const LocalAdValues& state, const LocalAdValues& residual, LocalJacobian* jacobian = nullptr) {
    LocalResidual result{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), result.data(), jacobian->data());
    return result;
}
inline void add_mechanical_point_residual(const RzQuadraturePoint& point, const AxisymmetricKinematics& kinematics,
    const AxisymmetricStress& stress, LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        residual[4 + node] += kinematics.weighted_measure * (stress.rr * kinematics.gradient_r[node] +
                                                                stress.hoop * point.shape[node] / kinematics.radius +
                                                                stress.rz * kinematics.gradient_z[node]);
        residual[8 + node] += kinematics.weighted_measure *
                              (stress.zz * kinematics.gradient_z[node] + stress.rz * kinematics.gradient_r[node]);
    }
}
inline void add_point_residual(const RzQuadraturePoint& point, const adlite::Scalar& gradient_temperature_r,
    const adlite::Scalar& gradient_temperature_z, const AxisymmetricKinematics& kinematics,
    const adlite::Scalar& heat_capacity, const adlite::Scalar& temperature_rate, const adlite::Scalar& conductivity,
    double volumetric_heat_source, const AxisymmetricStress& stress, LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        residual[node] +=
            point.weighted_measure * (heat_capacity * point.shape[node] * temperature_rate +
                                         conductivity * (point.gradient_r[node] * gradient_temperature_r +
                                                            point.gradient_z[node] * gradient_temperature_z) -
                                         volumetric_heat_source * point.shape[node]);
    }
    add_mechanical_point_residual(point, kinematics, stress, residual);
}
} // namespace quad4_rz_detail
AxisymmetricKinematics evaluate_axisymmetric_kinematics(
    const RzQuadraturePoint& point, const LocalAdValues& state, StrainFormulation strain_formulation) {
    const LocalValues undeformed{};
    return evaluate_axisymmetric_incremental_kinematics(point, state, undeformed, strain_formulation);
}
AxisymmetricKinematics evaluate_axisymmetric_incremental_kinematics(const RzQuadraturePoint& point,
    const LocalAdValues& current_state, const LocalValues& committed_state, StrainFormulation strain_formulation) {
    const adlite::Scalar radial_displacement = quad4_rz_detail::interpolate(point.shape, current_state, 4),
                         displacement_gradient_rr = quad4_rz_detail::interpolate(point.gradient_r, current_state, 4),
                         displacement_gradient_rz = quad4_rz_detail::interpolate(point.gradient_z, current_state, 4),
                         displacement_gradient_zr = quad4_rz_detail::interpolate(point.gradient_r, current_state, 8),
                         displacement_gradient_zz = quad4_rz_detail::interpolate(point.gradient_z, current_state, 8);
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
    if (!std::isfinite(deformation_hoop.value()) || !(deformation_hoop.value() > 0.0) ||
        !std::isfinite(current_radius.value()) || !(current_radius.value() > 0.0))
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
    if (!std::isfinite(old_determinant_rz) || !(old_determinant_rz > 0.0) || !std::isfinite(old_deformation_hoop) ||
        !(old_deformation_hoop > 0.0))
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
    if (!std::isfinite(incremental_determinant.value()) || !(incremental_determinant.value() > 0.0) ||
        !std::isfinite(incremental_hoop.value()) || !(incremental_hoop.value() > 0.0))
        throw std::domain_error("Incremental finite-strain RZ state requires positive Jacobian and hoop stretch");
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
        c2 = 0.125 + q * 0.03125 * (p2 - 12.0 * (p - 1.0)) / p2 + q2 * (p - 2.0) * (p2 - 10.0 * p + 32.0) / p3 +
             q3 * (1104.0 - 992.0 * p + 376.0 * p2 - 72.0 * p3 + 5.0 * p4) / (512.0 * p4);
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
namespace {
MaterialFunctionContext rz_material_context(double time, const RzQuadraturePoint& point) {
    return {time, point.radius, 0.0, point.axial_coordinate};
}
struct ThermoelasticPointResponse final {
    adlite::Scalar temperature, gradient_temperature_r, gradient_temperature_z;
    AxisymmetricKinematics kinematics;
    AxisymmetricStress stress;
};
ThermoelasticPointResponse point_response(const RzQuadraturePoint& point, const LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material, StrainFormulation strain_formulation, double time) {
    const adlite::Scalar temperature = quad4_rz_detail::interpolate(point.shape, state, 0);
    const AxisymmetricKinematics kinematics = evaluate_axisymmetric_kinematics(point, state, strain_formulation);
    AxisymmetricStress stress = material.stress(kinematics.strain_rr, kinematics.strain_zz, kinematics.strain_hoop,
        kinematics.strain_rz, temperature, rz_material_context(time, point));
    if (strain_formulation == StrainFormulation::finite)
        stress = rotate_axisymmetric_tensor(stress, kinematics.rotation);
    return {
        temperature,
        quad4_rz_detail::interpolate(point.gradient_r, state, 0),
        quad4_rz_detail::interpolate(point.gradient_z, state, 0),
        kinematics,
        stress,
    };
}
} // namespace
LocalResidual compute_quad4_rz_thermoelastic(
    const Quad4RzData& data, const Quad4RzGeometry& geometry, const LocalValues& state, LocalJacobian* jacobian) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state, jacobian != nullptr);
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    for (const RzQuadraturePoint& point : geometry.points) {
        const ThermoelasticPointResponse response =
            point_response(point, ad_state, data.material, data.strain_formulation, data.time);
        const adlite::Scalar conductivity =
            data.material.conductivity(response.temperature, rz_material_context(data.time, point));
        quad4_rz_detail::add_point_residual(point, response.gradient_temperature_r, response.gradient_temperature_z,
            response.kinematics, 0.0, 0.0, conductivity, data.volumetric_heat_source, response.stress, ad_residual);
    }
    return quad4_rz_detail::values(ad_state, ad_residual, jacobian);
}
std::array<AxisymmetricStressValues, 4> compute_quad4_rz_thermoelastic_stress(
    const Quad4RzData& data, const Quad4RzGeometry& geometry, const LocalValues& state) {
    const LocalAdValues passive_state = quad4_rz_detail::ad_state(state);
    std::array<AxisymmetricStressValues, 4> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const AxisymmetricStress& stress =
            point_response(geometry.points[q], passive_state, data.material, data.strain_formulation, data.time).stress;
        result[q] = {stress.rr.value(), stress.zz.value(), stress.hoop.value(), stress.rz.value()};
    }
    return result;
}
namespace {
struct TransientPointResponse final {
    adlite::Scalar temperature, gradient_temperature_r, gradient_temperature_z;
    AxisymmetricKinematics kinematics;
    double old_temperature;
    InelasticStressResponse response;
};
TransientPointResponse transient_point_response(const RzQuadraturePoint& point, const LocalAdValues& state,
    const LocalValues& committed_state, const IsotropicThermoelasticMaterial& material,
    const MaterialPointState& committed_material, double time_step, StrainFormulation strain_formulation, double time) {
    const adlite::Scalar temperature = quad4_rz_detail::interpolate(point.shape, state, 0);
    const AxisymmetricKinematics kinematics =
        evaluate_axisymmetric_incremental_kinematics(point, state, committed_state, strain_formulation);
    const double old_temperature = quad4_rz_detail::interpolate(point.shape, committed_state, 0);
    const MaterialFunctionContext context = rz_material_context(time, point);
    InelasticStressResponse response =
        strain_formulation == StrainFormulation::finite
            ? material.incremental_response(kinematics.strain_rr, kinematics.strain_zz, kinematics.strain_hoop,
                  kinematics.strain_rz, kinematics.rotation, temperature, old_temperature, time_step,
                  committed_material, context)
            : material.response(kinematics.strain_rr, kinematics.strain_zz, kinematics.strain_hoop,
                  kinematics.strain_rz, temperature, time_step, committed_material, context);
    return {temperature, quad4_rz_detail::interpolate(point.gradient_r, state, 0),
        quad4_rz_detail::interpolate(point.gradient_z, state, 0), kinematics, old_temperature, std::move(response)};
}
void validate_time_step(double time_step) {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("Transient Quad4 time_step must be finite and positive");
}
void validate_committed_state(const LocalValues& committed_state) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        const double temperature = committed_state[node];
        if (!std::isfinite(temperature) || !(temperature > 0.0))
            throw std::invalid_argument("Transient Quad4 committed temperatures must be finite and positive");
    }
}
} // namespace
LocalResidual compute_quad4_rz_transient(const Quad4RzData& data, const Quad4RzGeometry& geometry,
    const LocalValues& current_state, const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step, LocalJacobian* jacobian) {
    validate_time_step(time_step);
    validate_committed_state(committed_state);
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(current_state, jacobian != nullptr);
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const RzQuadraturePoint& point = geometry.points[q];
        const TransientPointResponse evaluation = transient_point_response(point, ad_state, committed_state,
            data.material, committed_material[q], time_step, data.strain_formulation, data.time);
        const adlite::Scalar temperature_rate = (evaluation.temperature - evaluation.old_temperature) / time_step;
        const MaterialFunctionContext context = rz_material_context(data.time, point);
        const adlite::Scalar conductivity = data.material.conductivity(evaluation.temperature, context);
        const adlite::Scalar heat_capacity = data.material.heat_capacity(evaluation.temperature, context);
        quad4_rz_detail::add_point_residual(point, evaluation.gradient_temperature_r, evaluation.gradient_temperature_z,
            evaluation.kinematics, heat_capacity, temperature_rate, conductivity, data.volumetric_heat_source,
            evaluation.response.stress, ad_residual);
    }
    return quad4_rz_detail::values(ad_state, ad_residual, jacobian);
}
Quad4MaterialUpdate compute_quad4_rz_transient_update(const Quad4RzData& data, const Quad4RzGeometry& geometry,
    const LocalValues& converged_state, const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material, double time_step) {
    validate_time_step(time_step);
    const LocalAdValues passive_state = quad4_rz_detail::ad_state(converged_state);
    Quad4MaterialUpdate result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const TransientPointResponse evaluation = transient_point_response(geometry.points[q], passive_state,
            committed_state, data.material, committed_material[q], time_step, data.strain_formulation, data.time);
        if (!std::isfinite(evaluation.temperature.value()) || !(evaluation.temperature.value() > 0.0))
            throw std::domain_error("Transient Quad4 trial temperature must be finite and positive");
        result.history[q] = IsotropicThermoelasticMaterial::state_values(evaluation.response.trial_state);
        result.stress[q] = {
            evaluation.response.stress.rr.value(),
            evaluation.response.stress.zz.value(),
            evaluation.response.stress.hoop.value(),
            evaluation.response.stress.rz.value(),
        };
    }
    return result;
}
namespace {
bool finite_point(const RzPoint& point) { return std::isfinite(point.r) && std::isfinite(point.z); }
void validate_line(const Line2InterfaceSideCoordinates& coordinates, const char* name) {
    for (const RzPoint& point : coordinates) {
        if (!finite_point(point)) throw std::invalid_argument(std::string(name) + " requires finite coordinates");
        if (!(point.r >= 0.0)) throw std::invalid_argument(std::string(name) + " requires nonnegative radii");
    }
    const double dr = coordinates[1].r - coordinates[0].r, dz = coordinates[1].z - coordinates[0].z;
    if (!(std::hypot(dr, dz) > 0.0)) throw std::invalid_argument(std::string(name) + " requires a nonzero line length");
}
void validate_edge(
    const std::array<RzPoint, 2>& coordinates, const std::array<std::size_t, 2>& local_nodes, const char* name) {
    if (local_nodes[0] >= 4 || local_nodes[1] >= 4 || local_nodes[0] == local_nodes[1])
        throw std::invalid_argument(std::string(name) + " requires two distinct Quad4 local nodes");
    validate_line(coordinates, name);
}
void displaced_edge_coordinates(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes, const LocalAdValues& state, bool use_displaced_geometry,
    std::array<adlite::Scalar, 2>& radius, std::array<adlite::Scalar, 2>& axial) {
    for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
        const std::size_t local = local_nodes[edge_node];
        radius[edge_node] = coordinates[edge_node].r;
        axial[edge_node] = coordinates[edge_node].z;
        if (use_displaced_geometry) {
            radius[edge_node] += state[4 + local];
            axial[edge_node] += state[8 + local];
        }
    }
}
} // namespace
Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(
    const std::array<RzPoint, 2>& coordinates, const std::array<std::size_t, 2>& local_nodes) {
    validate_edge(coordinates, local_nodes, "Boundary edge");
    return {coordinates, local_nodes};
}
namespace {
void compute_line2_rz_mechanical_residual_ad(const Line2RzBoundaryData& data, const Line2RzBoundaryGeometry& geometry,
    const LocalAdValues& state, LocalAdValues& residual) {
    residual.fill(adlite::Scalar(0.0));
    std::array<adlite::Scalar, 2> radius{};
    std::array<adlite::Scalar, 2> axial{};
    displaced_edge_coordinates(
        geometry.coordinates, geometry.local_nodes, state, data.use_displaced_geometry, radius, axial);
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const double xi : locations) {
        const std::array<double, 2> shape = {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
        const adlite::Scalar current_radius = shape[0] * radius[0] + shape[1] * radius[1];
        const adlite::Scalar dr_dxi = 0.5 * (radius[1] - radius[0]), dz_dxi = 0.5 * (axial[1] - axial[0]);
        const adlite::Scalar measure = data.kind == Line2RzBoundaryKind::pressure
                                           ? 2.0 * pi * current_radius
                                           : 2.0 * pi * current_radius * adlite::hypot(dr_dxi, dz_dxi);
        if (!std::isfinite(measure.value()) || !(measure.value() > 0.0))
            throw std::domain_error("Mechanical boundary current measure must be finite and positive");
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node) {
            const std::size_t local = geometry.local_nodes[edge_node];
            if (data.kind == Line2RzBoundaryKind::pressure) {
                residual[4 + local] += measure * data.load * shape[edge_node] * dz_dxi;
                residual[8 + local] -= measure * data.load * shape[edge_node] * dr_dxi;
            } else {
                const std::size_t offset = data.component == TractionComponent::radial ? 4 : 8;
                residual[offset + local] -= measure * data.load * shape[edge_node];
            }
        }
    }
}
void compute_line2_rz_convection_residual_ad(const Line2RzBoundaryData& data, const Line2RzBoundaryGeometry& geometry,
    const LocalAdValues& state, LocalAdValues& residual) {
    residual.fill(adlite::Scalar(0.0));
    const double dr = geometry.coordinates[1].r - geometry.coordinates[0].r,
                 dz = geometry.coordinates[1].z - geometry.coordinates[0].z, line_jacobian = 0.5 * std::hypot(dr, dz);
    const std::array<double, 2> locations = {-gauss, gauss};
    for (const double xi : locations) {
        const std::array<double, 2> shape = {0.5 * (1.0 - xi), 0.5 * (1.0 + xi)};
        const double radius = shape[0] * geometry.coordinates[0].r + shape[1] * geometry.coordinates[1].r;
        adlite::Scalar temperature = 0.0;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node)
            temperature += shape[edge_node] * state[geometry.local_nodes[edge_node]];
        const adlite::Scalar heat_flux = data.load * (temperature - data.ambient);
        const double measure = 2.0 * pi * radius * line_jacobian;
        for (std::size_t edge_node = 0; edge_node < 2; ++edge_node)
            residual[geometry.local_nodes[edge_node]] += measure * shape[edge_node] * heat_flux;
    }
}
} // namespace
LocalResidual compute_line2_rz_boundary(const Line2RzBoundaryData& data, const Line2RzBoundaryGeometry& geometry,
    const LocalValues& state, LocalJacobian* jacobian) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state, jacobian != nullptr);
    LocalAdValues ad_residual{};
    if (data.kind == Line2RzBoundaryKind::convection)
        compute_line2_rz_convection_residual_ad(data, geometry, ad_state, ad_residual);
    else
        compute_line2_rz_mechanical_residual_ad(data, geometry, ad_state, ad_residual);
    return quad4_rz_detail::values(ad_state, ad_residual, jacobian);
}
namespace {
struct HeatAdQuadratureValue final {
    bool projected = false;
    adlite::Scalar gap = 0.0, heat_flux = 0.0, weighted_measure = 0.0, primary_shape_0 = 0.0, primary_shape_1 = 0.0;
};
struct ContactAdValue final {
    bool projected = false;
    adlite::Scalar gap = 0.0, pressure = 0.0, tributary_area = 0.0, tributary_length = 0.0, contact_force = 0.0,
                   primary_shape_0 = 0.0, primary_shape_1 = 0.0, normal_r = 0.0, normal_z = 0.0, tangent_r = 0.0,
                   tangent_z = 0.0, tangential_traction = 0.0, tangential_force = 0.0, elastic_tangential_slip = 0.0;
    bool sliding = false;
};
double reference_projection_fraction(
    const RzPoint& secondary, const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double tangent_r = primary_coordinates[1].r - primary_coordinates[0].r,
                 tangent_z = primary_coordinates[1].z - primary_coordinates[0].z,
                 length_squared = tangent_r * tangent_r + tangent_z * tangent_z;
    return ((secondary.r - primary_coordinates[0].r) * tangent_r +
               (secondary.z - primary_coordinates[0].z) * tangent_z) /
           length_squared;
}
double reference_normal_orientation(const RzPoint& secondary, const Line2InterfaceSideCoordinates& primary_coordinates,
    double primary_fraction, double zero_gap_orientation_hint) {
    const double tangent_r = primary_coordinates[1].r - primary_coordinates[0].r,
                 tangent_z = primary_coordinates[1].z - primary_coordinates[0].z,
                 length = std::hypot(tangent_r, tangent_z),
                 delta_r = primary_coordinates[0].r + primary_fraction * tangent_r - secondary.r,
                 delta_z = primary_coordinates[0].z + primary_fraction * tangent_z - secondary.z,
                 raw_gap = delta_r * tangent_z / length - delta_z * tangent_r / length;
    if (raw_gap == 0.0) {
        constexpr double endpoint_tolerance = 1.0e-12;
        const double projection = reference_projection_fraction(secondary, primary_coordinates);
        if (projection >= -endpoint_tolerance && projection <= 1.0 + endpoint_tolerance) {
            if (zero_gap_orientation_hint == 0.0)
                throw std::invalid_argument("Contact requires a positive gap or nonzero material-side hint");
            return zero_gap_orientation_hint < 0.0 ? 1.0 : -1.0;
        }
        const double normal_r = tangent_z / length, normal_z = -tangent_r / length,
                     center_r = 0.5 * (primary_coordinates[0].r + primary_coordinates[1].r),
                     center_z = 0.5 * (primary_coordinates[0].z + primary_coordinates[1].z),
                     side = (secondary.r - center_r) * normal_r + (secondary.z - center_z) * normal_z;
        return side >= 0.0 ? 1.0 : -1.0;
    }
    return raw_gap > 0.0 ? 1.0 : -1.0;
}
adlite::Scalar interpolate(
    const std::array<double, line2_interface_side_node_count>& shape, const LocalAdValues& state, std::size_t offset) {
    adlite::Scalar value = 0.0;
    for (std::size_t node = 0; node < line2_interface_side_node_count; ++node)
        value += shape[node] * state[offset + node];
    return value;
}
bool projection_is_inside(double fraction, bool includes_second_endpoint) {
    if (fraction < 0.0) return false;
    if (includes_second_endpoint) return fraction <= 1.0;
    return fraction < 1.0;
}
HeatAdQuadratureValue evaluate_heat_quadrature(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, const Line2RzHeatQuadraturePoint& point,
    const LocalAdValues& state, const GapHeatProperties& properties, bool includes_second_endpoint) {
    const adlite::Scalar secondary_radius_0 = secondary_coordinates[0].r + state[4],
                         secondary_radius_1 = secondary_coordinates[1].r + state[5],
                         secondary_axial_0 = secondary_coordinates[0].z + state[8],
                         secondary_axial_1 = secondary_coordinates[1].z + state[9];
    const adlite::Scalar secondary_radius =
        point.secondary_shape[0] * secondary_radius_0 + point.secondary_shape[1] * secondary_radius_1;
    const adlite::Scalar secondary_axial =
        point.secondary_shape[0] * secondary_axial_0 + point.secondary_shape[1] * secondary_axial_1;
    const adlite::Scalar primary_radius_0 = primary_coordinates[0].r + state[6],
                         primary_radius_1 = primary_coordinates[1].r + state[7],
                         primary_axial_0 = primary_coordinates[0].z + state[10],
                         primary_axial_1 = primary_coordinates[1].z + state[11],
                         tangent_r = primary_radius_1 - primary_radius_0, tangent_z = primary_axial_1 - primary_axial_0,
                         tangent_length = adlite::hypot(tangent_r, tangent_z);
    adlite::Scalar primary_fraction =
        ((secondary_radius - primary_radius_0) * tangent_r + (secondary_axial - primary_axial_0) * tangent_z) /
        (tangent_length * tangent_length);
    if (!projection_is_inside(primary_fraction.value(), includes_second_endpoint)) return {};
    const adlite::Scalar primary_shape_0 = 1.0 - primary_fraction, primary_shape_1 = primary_fraction,
                         primary_radius = primary_shape_0 * primary_radius_0 + primary_shape_1 * primary_radius_1,
                         primary_axial = primary_shape_0 * primary_axial_0 + primary_shape_1 * primary_axial_1,
                         normal_r = point.normal_orientation * tangent_z / tangent_length,
                         normal_z = -point.normal_orientation * tangent_r / tangent_length;
    const adlite::Scalar gap =
        (primary_radius - secondary_radius) * normal_r + (primary_axial - secondary_axial) * normal_z;
    const adlite::Scalar secondary_temperature = interpolate(point.secondary_shape, state, 0),
                         primary_temperature = primary_shape_0 * state[2] + primary_shape_1 * state[3],
                         thermal_gap = adlite::max(gap, adlite::Scalar(properties.minimum_gap)),
                         conductance = properties.gap_conductivity / thermal_gap,
                         heat_flux = conductance * (secondary_temperature - primary_temperature),
                         dr_dxi = 0.5 * (secondary_radius_1 - secondary_radius_0),
                         dz_dxi = 0.5 * (secondary_axial_1 - secondary_axial_0),
                         surface_jacobian = adlite::sqrt(dr_dxi * dr_dxi + dz_dxi * dz_dxi),
                         weighted_measure = 2.0 * pi * secondary_radius * surface_jacobian * point.integration_weight;
    return {true, gap, heat_flux, weighted_measure, primary_shape_0, primary_shape_1};
}
bool clamp_owned_chain_endpoint(
    adlite::Scalar& fraction, double reference_fraction, bool primary_segment_is_first, bool includes_second_endpoint) {
    constexpr double endpoint_tolerance = 1.0e-12;
    const double value = fraction.value();
    if (primary_segment_is_first && value < 0.0 && std::abs(reference_fraction) <= endpoint_tolerance) {
        fraction = 0.0;
        return true;
    }
    if (includes_second_endpoint && value > 1.0 && std::abs(reference_fraction - 1.0) <= endpoint_tolerance) {
        fraction = 1.0;
        return true;
    }
    return false;
}
void evaluate_friction(ContactAdValue& value, const NodeToLineRzContactGeometry& geometry, const LocalAdValues& state,
    const LocalValues& committed_state, const ContactPointHistory& history, const NormalContactProperties& properties) {
    if (properties.friction_coefficient == 0.0 || !(value.pressure.value() > 0.0)) return;
    const std::size_t secondary = geometry.secondary_local_node;
    const adlite::Scalar secondary_increment_r = state[4 + secondary] - committed_state[4 + secondary],
                         secondary_increment_z = state[8 + secondary] - committed_state[8 + secondary];
    const adlite::Scalar primary_increment_r = value.primary_shape_0 * (state[6] - committed_state[6]) +
                                               value.primary_shape_1 * (state[7] - committed_state[7]);
    const adlite::Scalar primary_increment_z = value.primary_shape_0 * (state[10] - committed_state[10]) +
                                               value.primary_shape_1 * (state[11] - committed_state[11]);
    const adlite::Scalar tangential_increment = (secondary_increment_r - primary_increment_r) * value.tangent_r +
                                                (secondary_increment_z - primary_increment_z) * value.tangent_z;
    const adlite::Scalar trial_slip = history.elastic_tangential_slip + tangential_increment,
                         trial_traction = properties.penalty * trial_slip,
                         sliding_limit = properties.friction_coefficient * value.pressure;
    const double trial_magnitude = std::abs(trial_traction.value());
    if (trial_magnitude < sliding_limit.value() || (trial_magnitude == sliding_limit.value() && !history.sliding)) {
        value.tangential_traction = trial_traction;
        value.elastic_tangential_slip = trial_slip;
    } else {
        const double direction = trial_traction.value() < 0.0 ? -1.0 : 1.0;
        value.tangential_traction = direction * sliding_limit;
        value.elastic_tangential_slip = value.tangential_traction / properties.penalty;
        value.sliding = true;
    }
    value.tangential_force = value.tangential_traction * value.tributary_area;
}
ContactAdValue evaluate_contact(const NodeToLineRzContactGeometry& geometry, const LocalAdValues& state,
    const LocalValues& committed_state, const ContactPointHistory& history, const NormalContactProperties& properties) {
    const std::size_t secondary = geometry.secondary_local_node, other = secondary == 0 ? 1 : 0;
    const adlite::Scalar secondary_radius = geometry.secondary_edge_coordinates[secondary].r + state[4 + secondary],
                         secondary_z = geometry.secondary_edge_coordinates[secondary].z + state[8 + secondary],
                         primary_radius_0 = geometry.primary_segment_coordinates[0].r + state[6],
                         primary_radius_1 = geometry.primary_segment_coordinates[1].r + state[7],
                         primary_z_0 = geometry.primary_segment_coordinates[0].z + state[10],
                         primary_z_1 = geometry.primary_segment_coordinates[1].z + state[11],
                         tangent_r = primary_radius_1 - primary_radius_0, tangent_z = primary_z_1 - primary_z_0,
                         tangent_length = adlite::hypot(tangent_r, tangent_z);
    adlite::Scalar primary_fraction =
        ((secondary_radius - primary_radius_0) * tangent_r + (secondary_z - primary_z_0) * tangent_z) /
        (tangent_length * tangent_length);
    bool projected = projection_is_inside(primary_fraction.value(), geometry.primary_segment_includes_second_endpoint);
    if (!projected)
        projected = clamp_owned_chain_endpoint(primary_fraction, geometry.reference_primary_fraction,
            geometry.primary_segment_is_first, geometry.primary_segment_includes_second_endpoint);
    if (!projected) return {};
    const adlite::Scalar primary_shape_0 = 1.0 - primary_fraction, primary_shape_1 = primary_fraction,
                         primary_radius = primary_shape_0 * primary_radius_0 + primary_shape_1 * primary_radius_1,
                         primary_z = primary_shape_0 * primary_z_0 + primary_shape_1 * primary_z_1,
                         normal_r = geometry.normal_orientation * tangent_z / tangent_length,
                         normal_z = -geometry.normal_orientation * tangent_r / tangent_length,
                         gap = (primary_radius - secondary_radius) * normal_r + (primary_z - secondary_z) * normal_z;
    const adlite::Scalar multiplier =
        properties.augmented_lagrangian ? adlite::Scalar(history.normal_multiplier) : adlite::Scalar(0.0);
    const adlite::Scalar pressure = adlite::max(multiplier - properties.penalty * gap, adlite::Scalar(0.0)),
                         other_radius = geometry.secondary_edge_coordinates[other].r + state[4 + other],
                         other_z = geometry.secondary_edge_coordinates[other].z + state[8 + other],
                         dr = other_radius - secondary_radius, dz = other_z - secondary_z,
                         edge_length = adlite::sqrt(dr * dr + dz * dz), tributary_length = 0.5 * edge_length,
                         tributary_area = 2.0 * pi * 0.5 * edge_length * (2.0 * secondary_radius + other_radius) / 3.0,
                         contact_force = pressure * tributary_area;
    ContactAdValue result;
    result.projected = true;
    result.gap = gap;
    result.pressure = pressure;
    result.tributary_area = tributary_area;
    result.tributary_length = tributary_length;
    result.contact_force = contact_force;
    result.primary_shape_0 = primary_shape_0;
    result.primary_shape_1 = primary_shape_1;
    result.normal_r = normal_r;
    result.normal_z = normal_z;
    result.tangent_r = tangent_r / tangent_length;
    result.tangent_z = tangent_z / tangent_length;
    evaluate_friction(result, geometry, state, committed_state, history, properties);
    return result;
}
} // namespace
std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> make_line2_rz_heat_quadrature(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, double secondary_coordinate_lower,
    double secondary_coordinate_upper, double zero_gap_orientation_hint) {
    if (!std::isfinite(secondary_coordinate_lower) || !std::isfinite(secondary_coordinate_upper) ||
        secondary_coordinate_lower < -1.0 || secondary_coordinate_upper > 1.0 ||
        !(secondary_coordinate_upper > secondary_coordinate_lower))
        throw std::invalid_argument("Heat quadrature requires a nonempty secondary interval in [-1,1]");
    const std::array<double, line2_interface_quadrature_point_count> locations = {-gauss, gauss};
    std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> result{};
    for (std::size_t q = 0; q < result.size(); ++q) {
        const double secondary_xi = 0.5 * ((1.0 - locations[q]) * secondary_coordinate_lower +
                                              (1.0 + locations[q]) * secondary_coordinate_upper);
        const std::array<double, 2> secondary_shape = {
            0.5 * (1.0 - secondary_xi),
            0.5 * (1.0 + secondary_xi),
        };
        result[q] = make_line2_rz_heat_point_geometry(secondary_coordinates, primary_coordinates, secondary_shape,
            0.5 * (secondary_coordinate_upper - secondary_coordinate_lower), true, zero_gap_orientation_hint)
                        .point;
        const double primary_fraction = result[q].primary_shape[1];
        if (primary_fraction < 0.0 || primary_fraction > 1.0)
            throw std::invalid_argument("Heat-contact Gauss point lies outside the primary segment");
    }
    return result;
}
Line2RzHeatPointGeometry make_line2_rz_heat_point_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const std::array<double, line2_interface_side_node_count>& secondary_shape, double integration_weight,
    bool primary_segment_includes_second_endpoint, double zero_gap_orientation_hint) {
    validate_line(secondary_coordinates, "Line2RzHeatPointGeometry secondary");
    validate_line(primary_coordinates, "Line2RzHeatPointGeometry primary");
    const RzPoint secondary_point = {
        secondary_shape[0] * secondary_coordinates[0].r + secondary_shape[1] * secondary_coordinates[1].r,
        secondary_shape[0] * secondary_coordinates[0].z + secondary_shape[1] * secondary_coordinates[1].z,
    };
    const double primary_fraction = reference_projection_fraction(secondary_point, primary_coordinates),
                 closest_fraction = std::max(0.0, std::min(1.0, primary_fraction));
    return {
        secondary_coordinates,
        primary_coordinates,
        {secondary_shape, {1.0 - primary_fraction, primary_fraction}, integration_weight,
            reference_normal_orientation(
                secondary_point, primary_coordinates, closest_fraction, zero_gap_orientation_hint)},
        primary_segment_includes_second_endpoint,
    };
}
LocalResidual compute_line2_rz_gap_heat(const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry,
    const LocalValues& state, LocalJacobian* jacobian) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state, jacobian != nullptr);
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    const HeatAdQuadratureValue value =
        evaluate_heat_quadrature(geometry.secondary_coordinates, geometry.primary_coordinates, geometry.point, ad_state,
            properties, geometry.primary_segment_includes_second_endpoint);
    if (value.projected) {
        for (std::size_t node = 0; node < line2_interface_side_node_count; ++node) {
            ad_residual[node] += value.weighted_measure * geometry.point.secondary_shape[node] * value.heat_flux;
            const adlite::Scalar primary_shape = node == 0 ? value.primary_shape_0 : value.primary_shape_1;
            ad_residual[2 + node] -= value.weighted_measure * primary_shape * value.heat_flux;
        }
    }
    return quad4_rz_detail::values(ad_state, ad_residual, jacobian);
}
HeatQuadratureValue compute_line2_rz_gap_heat_value(
    const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry, const LocalValues& state) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state);
    const HeatAdQuadratureValue value =
        evaluate_heat_quadrature(geometry.secondary_coordinates, geometry.primary_coordinates, geometry.point, ad_state,
            properties, geometry.primary_segment_includes_second_endpoint);
    return {value.projected, value.gap.value(), value.heat_flux.value(), value.weighted_measure.value()};
}
NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates, std::size_t secondary_local_node,
    bool primary_segment_is_first, bool primary_segment_includes_upper_endpoint, double zero_gap_orientation_hint) {
    validate_line(secondary_edge_coordinates, "NodeToLineRzContactGeometry secondary");
    validate_line(primary_segment_coordinates, "NodeToLineRzContactGeometry primary");
    if (secondary_local_node >= line2_interface_side_node_count)
        throw std::invalid_argument("NodeToLineRzContactGeometry secondary node is out of range");
    const RzPoint secondary_point = secondary_edge_coordinates[secondary_local_node];
    const double primary_fraction = reference_projection_fraction(secondary_point, primary_segment_coordinates),
                 closest_fraction = std::max(0.0, std::min(1.0, primary_fraction));
    const double orientation = reference_normal_orientation(
        secondary_point, primary_segment_coordinates, closest_fraction, zero_gap_orientation_hint);
    return {
        secondary_edge_coordinates,
        primary_segment_coordinates,
        secondary_local_node,
        primary_segment_is_first,
        primary_segment_includes_upper_endpoint,
        orientation,
        primary_fraction,
    };
}
LocalResidual compute_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history, LocalJacobian* jacobian) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state, jacobian != nullptr);
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    const ContactAdValue value = evaluate_contact(geometry, ad_state, committed_state, history, properties);
    if (value.projected) {
        const std::size_t secondary = geometry.secondary_local_node;
        ad_residual[4 + secondary] += value.contact_force * value.normal_r;
        ad_residual[6] -= value.primary_shape_0 * value.contact_force * value.normal_r;
        ad_residual[7] -= value.primary_shape_1 * value.contact_force * value.normal_r;
        ad_residual[8 + secondary] += value.contact_force * value.normal_z;
        ad_residual[10] -= value.primary_shape_0 * value.contact_force * value.normal_z;
        ad_residual[11] -= value.primary_shape_1 * value.contact_force * value.normal_z;
        if (properties.friction_coefficient != 0.0) {
            ad_residual[4 + secondary] += value.tangential_force * value.tangent_r;
            ad_residual[6] -= value.primary_shape_0 * value.tangential_force * value.tangent_r;
            ad_residual[7] -= value.primary_shape_1 * value.tangential_force * value.tangent_r;
            ad_residual[8 + secondary] += value.tangential_force * value.tangent_z;
            ad_residual[10] -= value.primary_shape_0 * value.tangential_force * value.tangent_z;
            ad_residual[11] -= value.primary_shape_1 * value.tangential_force * value.tangent_z;
        }
    }
    return quad4_rz_detail::values(ad_state, ad_residual, jacobian);
}
ContactPointValue compute_node_to_line_rz_contact_value(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state);
    const ContactAdValue result = evaluate_contact(geometry, ad_state, committed_state, history, properties);
    return {
        result.projected,
        result.gap.value(),
        result.pressure.value(),
        result.tributary_area.value(),
        result.tributary_length.value(),
        result.contact_force.value(),
        result.tangential_traction.value(),
        result.tangential_force.value(),
        result.elastic_tangential_slip.value(),
        result.sliding,
    };
}
} // namespace fuelsim
