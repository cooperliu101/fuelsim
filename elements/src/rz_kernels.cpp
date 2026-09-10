#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "contact.hpp"
#include "detail/ad_local_system.hpp"
#include "detail/rz_point.hpp"
#include "rz_quad4.hpp"
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
} // namespace

namespace quad4_rz_detail {
inline LocalAdValues ad_state(const LocalValues& state, bool active = false) {
    LocalAdValues result{};
    if (active)
        ad_local_system::make_active(state.data(), state.size(), result.data());
    else
        ad_local_system::make_passive(state.data(), state.size(), result.data());
    return result;
}

inline LocalResidual
values(const LocalAdValues& state, const LocalAdValues& residual, LocalJacobian* jacobian = nullptr) {
    LocalResidual result{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), result.data(), jacobian->data());
    return result;
}

inline void add_mechanical_point_residual(const RzQuadraturePoint& point,
    const AxisymmetricKinematics& kinematics,
    const AxisymmetricStress& stress,
    LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        residual[4 + node] +=
            kinematics.weighted_measure
            * (stress.rr * kinematics.gradient_r[node] + stress.hoop * point.shape[node] / kinematics.radius
                + stress.rz * kinematics.gradient_z[node]);
        residual[8 + node] += kinematics.weighted_measure
                              * (stress.zz * kinematics.gradient_z[node] + stress.rz * kinematics.gradient_r[node]);
    }
}

// Null heat-capacity inputs keep steady evaluations free of transient AD work.
inline void add_point_residual(const RzQuadraturePoint& point,
    const adlite::Scalar& gradient_temperature_r,
    const adlite::Scalar& gradient_temperature_z,
    const AxisymmetricKinematics& kinematics,
    const adlite::Scalar* heat_capacity,
    const adlite::Scalar* temperature_rate,
    const adlite::Scalar& conductivity,
    double volumetric_heat_source,
    const AxisymmetricStress& stress,
    LocalAdValues& residual) {
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        const adlite::Scalar thermal = heat_capacity == nullptr
                                           ? conductivity
                                                 * (point.gradient_r[node] * gradient_temperature_r
                                                     + point.gradient_z[node] * gradient_temperature_z)
                                           : *heat_capacity * point.shape[node] * *temperature_rate
                                                 + conductivity
                                                       * (point.gradient_r[node] * gradient_temperature_r
                                                           + point.gradient_z[node] * gradient_temperature_z);
        residual[node] += point.weighted_measure * (thermal - volumetric_heat_source * point.shape[node]);
    }
    add_mechanical_point_residual(point, kinematics, stress, residual);
}
} // namespace quad4_rz_detail

namespace {
MaterialFunctionContext rz_material_context(double time, const RzQuadraturePoint& point) {
    return {time, point.radius, 0.0, point.axial_coordinate};
}

struct ThermoelasticPointResponse final {
    adlite::Scalar temperature, gradient_temperature_r, gradient_temperature_z;
    AxisymmetricKinematics kinematics;
    AxisymmetricStress stress;
};

ThermoelasticPointResponse point_response(const RzQuadraturePoint& point,
    const LocalAdValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time) {
    const adlite::Scalar temperature = quad4_rz_detail::interpolate(point.shape, state, 0);
    const AxisymmetricKinematics kinematics = evaluate_axisymmetric_kinematics(point, state, strain_formulation);
    AxisymmetricStress stress = material.stress(kinematics.strain_rr,
        kinematics.strain_zz,
        kinematics.strain_hoop,
        kinematics.strain_rz,
        temperature,
        rz_material_context(time, point));
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

// Assembles one quadrature point's residual and exact 12-by-12 Jacobian with narrow AD. The
// kinematics chain is seeded on the four in-plane displacement-gradient components, the
// quadrature-point radial displacement, and the quadrature-point temperature (width 6), while
// the constitutive evaluation uses width 5 and is reattached with adlite::compose, so stress
// rotation and current-configuration geometry derivatives remain exact. The final chain from
// the point seeds to the 12 local DOFs is linear and applied in closed form:
// d(gradient_rr)/dur[b] = dN_b/dr, d(gradient_rz)/dur[b] = dN_b/dz, d(gradient_zr)/duz[b] =
// dN_b/dr, d(gradient_zz)/duz[b] = dN_b/dz, d(radial_displacement)/dur[b] = N_b, and
// d(temperature)/dT[b] = N_b.
void add_quad4_rz_point_system(const RzQuadraturePoint& point,
    const LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation strain_formulation,
    double time,
    double volumetric_heat_source,
    const LocalValues* committed_state,
    const MaterialPointState* committed_material,
    double time_step,
    LocalAdValues& residual,
    LocalJacobian& jacobian,
    bool include_thermal_time_term) {
    constexpr std::size_t point_width = 6, radial_index = 4, temperature_index = 5;
    std::array<adlite::Scalar, point_width> active{};
    active[0] = adlite::Scalar::independent(quad4_rz_detail::interpolate(point.gradient_r, state, 4), 0, point_width);
    active[1] = adlite::Scalar::independent(quad4_rz_detail::interpolate(point.gradient_z, state, 4), 1, point_width);
    active[2] = adlite::Scalar::independent(quad4_rz_detail::interpolate(point.gradient_r, state, 8), 2, point_width);
    active[3] = adlite::Scalar::independent(quad4_rz_detail::interpolate(point.gradient_z, state, 8), 3, point_width);
    active[radial_index] =
        adlite::Scalar::independent(quad4_rz_detail::interpolate(point.shape, state, 4), radial_index, point_width);
    const double temperature_value = quad4_rz_detail::interpolate(point.shape, state, 0);
    active[temperature_index] = adlite::Scalar::independent(temperature_value, temperature_index, point_width);
    const LocalValues undeformed{};
    const LocalValues& old_state = committed_state == nullptr ? undeformed : *committed_state;
    const AxisymmetricKinematics kinematics = evaluate_axisymmetric_kinematics_from_point(point,
        active[radial_index],
        active[0],
        active[1],
        active[2],
        active[3],
        old_state,
        strain_formulation);
    const MaterialFunctionContext context = rz_material_context(time, point);
    const std::array<const adlite::Scalar*, 4> strain_components = {&kinematics.strain_rr,
        &kinematics.strain_zz,
        &kinematics.strain_hoop,
        &kinematics.strain_rz};
    std::array<double, 4> fed_strain{};
    if (committed_material != nullptr && strain_formulation == StrainFormulation::finite) {
        const double old_temperature = quad4_rz_detail::interpolate(point.shape, old_state, 0);
        if (!std::isfinite(old_temperature) || !(old_temperature > 0.0))
            throw std::domain_error("Incremental material committed temperature must be finite and positive");
        MaterialFunctionContext old_context = context;
        old_context.time -= time_step;
        const AxisymmetricStrain old_imposed = material.eigenstrain_rz(adlite::Scalar(old_temperature), old_context);
        const std::array<double, 4> imposed = {old_imposed.rr.value(),
            old_imposed.zz.value(),
            old_imposed.hoop.value(),
            old_imposed.rz.value()};
        for (std::size_t component = 0; component < 4; ++component)
            fed_strain[component] = committed_material->elastic_strain[component]
                                    + strain_components[component]->value() + imposed[component]
                                    + committed_material->plastic_strain[component]
                                    + committed_material->creep_strain[component];
    } else {
        for (std::size_t component = 0; component < 4; ++component)
            fed_strain[component] = strain_components[component]->value();
    }
    const AxisymmetricStressTangent tangent = evaluate_axisymmetric_stress_tangent(material,
        fed_strain,
        temperature_value,
        time_step,
        committed_material,
        context);
    const std::array<adlite::Scalar, 5> compose_inputs = {*strain_components[0],
        *strain_components[1],
        *strain_components[2],
        *strain_components[3],
        active[temperature_index]};
    const std::array<double, 4> stress_values = {tangent.stress.rr,
        tangent.stress.zz,
        tangent.stress.hoop,
        tangent.stress.rz};
    std::array<double, 5> partials{};
    std::array<adlite::Scalar, 4> composed{};
    for (std::size_t component = 0; component < 4; ++component) {
        for (std::size_t column = 0; column < 4; ++column)
            partials[column] = tangent.tangent[component][column];
        partials[4] = tangent.thermal[component];
        composed[component] =
            adlite::compose(stress_values[component], compose_inputs.data(), partials.data(), compose_inputs.size());
    }
    AxisymmetricStress stress{composed[0], composed[1], composed[2], composed[3]};
    if (strain_formulation == StrainFormulation::finite)
        stress = rotate_axisymmetric_tensor(stress, kinematics.rotation);
    const adlite::Scalar conductivity = material.conductivity(active[temperature_index], context);
    adlite::Scalar temperature_rate = 0.0, heat_capacity = 0.0;
    if (committed_state != nullptr && include_thermal_time_term) {
        const double old_temperature = quad4_rz_detail::interpolate(point.shape, old_state, 0);
        temperature_rate = (active[temperature_index] - old_temperature) / time_step;
        heat_capacity = material.heat_capacity(active[temperature_index], context);
    }
    const double gradient_temperature_r = quad4_rz_detail::interpolate(point.gradient_r, state, 0),
                 gradient_temperature_z = quad4_rz_detail::interpolate(point.gradient_z, state, 0);
    LocalAdValues point_residual{};
    point_residual.fill(adlite::Scalar(0.0));
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        point_residual[node] +=
            point.weighted_measure
            * (conductivity
                    * (point.gradient_r[node] * gradient_temperature_r
                        + point.gradient_z[node] * gradient_temperature_z)
                + point.shape[node] * heat_capacity * temperature_rate - point.shape[node] * volumetric_heat_source);
        point_residual[4 + node] +=
            kinematics.weighted_measure
            * (stress.rr * kinematics.gradient_r[node] + stress.hoop * point.shape[node] / kinematics.radius
                + stress.rz * kinematics.gradient_z[node]);
        point_residual[8 + node] +=
            kinematics.weighted_measure
            * (stress.zz * kinematics.gradient_z[node] + stress.rz * kinematics.gradient_r[node]);
    }
    std::array<double, point_width> derivatives{};
    const double conductivity_value = conductivity.value();
    for (std::size_t node = 0; node < quad4_node_count; ++node) {
        point_residual[node].copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t other = 0; other < quad4_node_count; ++other)
            jacobian[node * local_dof_count + other] += point.weighted_measure * conductivity_value
                                                            * (point.gradient_r[node] * point.gradient_r[other]
                                                                + point.gradient_z[node] * point.gradient_z[other])
                                                        + point.shape[other] * derivatives[temperature_index];
        for (std::size_t equation = 0; equation < 2; ++equation) {
            const std::size_t row = 4 * (equation + 1) + node;
            point_residual[row].copy_derivatives(derivatives.data(), derivatives.size());
            for (std::size_t other = 0; other < quad4_node_count; ++other) {
                jacobian[row * local_dof_count + other] += derivatives[temperature_index] * point.shape[other];
                jacobian[row * local_dof_count + 4 + other] += derivatives[0] * point.gradient_r[other]
                                                               + derivatives[1] * point.gradient_z[other]
                                                               + derivatives[radial_index] * point.shape[other];
                jacobian[row * local_dof_count + 8 + other] +=
                    derivatives[2] * point.gradient_r[other] + derivatives[3] * point.gradient_z[other];
            }
        }
    }
    for (std::size_t row = 0; row < residual.size(); ++row)
        residual[row] += point_residual[row];
}
} // namespace

LocalResidual compute_quad4_rz_thermoelastic(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian) {
    if (data.element_formulation == RzElementFormulation::cax4rt) {
        const auto result = rz::compute_cax4rt(data, geometry, state, {}, nullptr, 0.0, jacobian != nullptr, false);
        if (jacobian)
            *jacobian = result.jacobian;
        return result.residual;
    }
    if (data.element_formulation == RzElementFormulation::cax4t) {
        const auto result = elements::evaluate_cax4t({data.material,
                                                         geometry,
                                                         state,
                                                         {},
                                                         nullptr,
                                                         0.0,
                                                         data.time,
                                                         data.volumetric_heat_source,
                                                         data.strain_formulation,
                                                         false},
            jacobian != nullptr);
        if (jacobian)
            *jacobian = result.jacobian;
        return result.residual;
    }
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    if (jacobian != nullptr) {
        jacobian->fill(0.0);
        for (const RzQuadraturePoint& point : geometry.points)
            add_quad4_rz_point_system(point,
                state,
                data.material,
                data.strain_formulation,
                data.time,
                data.volumetric_heat_source,
                nullptr,
                nullptr,
                0.0,
                ad_residual,
                *jacobian,
                true);
        LocalResidual result{};
        ad_local_system::extract_residual(ad_residual.data(), ad_residual.size(), result.data());
        return result;
    }
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state);
    for (const RzQuadraturePoint& point : geometry.points) {
        const ThermoelasticPointResponse response =
            point_response(point, ad_state, data.material, data.strain_formulation, data.time);
        const adlite::Scalar conductivity =
            data.material.conductivity(response.temperature, rz_material_context(data.time, point));
        quad4_rz_detail::add_point_residual(point,
            response.gradient_temperature_r,
            response.gradient_temperature_z,
            response.kinematics,
            nullptr,
            nullptr,
            conductivity,
            data.volumetric_heat_source,
            response.stress,
            ad_residual);
    }
    return quad4_rz_detail::values(ad_state, ad_residual);
}

std::array<AxisymmetricStressValues, 4> compute_quad4_rz_thermoelastic_stress(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& state) {
    if (data.element_formulation == RzElementFormulation::cax4rt) {
        const auto result = rz::compute_cax4rt(data, geometry, state, {}, nullptr, 0.0, false, false);
        std::array<AxisymmetricStressValues, 4> stress{};
        stress[0] = result.history[0].stress;
        return stress;
    }
    if (data.element_formulation == RzElementFormulation::cax4t) {
        const auto result = elements::evaluate_cax4t({data.material,
                                                         geometry,
                                                         state,
                                                         {},
                                                         nullptr,
                                                         0.0,
                                                         data.time,
                                                         data.volumetric_heat_source,
                                                         data.strain_formulation,
                                                         false},
            false);
        std::array<AxisymmetricStressValues, 4> stress;
        for (std::size_t q = 0; q < 4; ++q)
            stress[q] = result.history[q].stress;
        return stress;
    }
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

TransientPointResponse transient_point_response(const RzQuadraturePoint& point,
    const LocalAdValues& state,
    const LocalValues& committed_state,
    const IsotropicThermoelasticMaterial& material,
    const MaterialPointState& committed_material,
    double time_step,
    StrainFormulation strain_formulation,
    double time) {
    const adlite::Scalar temperature = quad4_rz_detail::interpolate(point.shape, state, 0);
    const AxisymmetricKinematics kinematics =
        evaluate_axisymmetric_incremental_kinematics(point, state, committed_state, strain_formulation);
    const double old_temperature = quad4_rz_detail::interpolate(point.shape, committed_state, 0);
    const MaterialFunctionContext context = rz_material_context(time, point);
    InelasticStressResponse response = strain_formulation == StrainFormulation::finite
                                           ? material.incremental_response(kinematics.strain_rr,
                                                 kinematics.strain_zz,
                                                 kinematics.strain_hoop,
                                                 kinematics.strain_rz,
                                                 kinematics.rotation,
                                                 temperature,
                                                 old_temperature,
                                                 time_step,
                                                 committed_material,
                                                 context)
                                           : material.response(kinematics.strain_rr,
                                                 kinematics.strain_zz,
                                                 kinematics.strain_hoop,
                                                 kinematics.strain_rz,
                                                 temperature,
                                                 time_step,
                                                 committed_material,
                                                 context);
    return {temperature,
        quad4_rz_detail::interpolate(point.gradient_r, state, 0),
        quad4_rz_detail::interpolate(point.gradient_z, state, 0),
        kinematics,
        old_temperature,
        std::move(response)};
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

LocalResidual compute_quad4_rz_transient(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& current_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material,
    double time_step,
    LocalJacobian* jacobian,
    bool include_thermal_time_term) {
    validate_time_step(time_step);
    validate_committed_state(committed_state);
    if (data.element_formulation == RzElementFormulation::cax4rt) {
        const auto result = rz::compute_cax4rt(data,
            geometry,
            current_state,
            committed_state,
            &committed_material,
            time_step,
            jacobian != nullptr,
            include_thermal_time_term);
        if (jacobian)
            *jacobian = result.jacobian;
        return result.residual;
    }
    if (data.element_formulation == RzElementFormulation::cax4t) {
        const auto result = elements::evaluate_cax4t({data.material,
                                                         geometry,
                                                         current_state,
                                                         committed_state,
                                                         &committed_material,
                                                         time_step,
                                                         data.time,
                                                         data.volumetric_heat_source,
                                                         data.strain_formulation,
                                                         include_thermal_time_term},
            jacobian != nullptr);
        if (jacobian)
            *jacobian = result.jacobian;
        return result.residual;
    }
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    if (jacobian != nullptr) {
        jacobian->fill(0.0);
        for (std::size_t q = 0; q < geometry.points.size(); ++q)
            add_quad4_rz_point_system(geometry.points[q],
                current_state,
                data.material,
                data.strain_formulation,
                data.time,
                data.volumetric_heat_source,
                &committed_state,
                &committed_material[q],
                time_step,
                ad_residual,
                *jacobian,
                include_thermal_time_term);
        LocalResidual result{};
        ad_local_system::extract_residual(ad_residual.data(), ad_residual.size(), result.data());
        return result;
    }
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(current_state);
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const RzQuadraturePoint& point = geometry.points[q];
        const TransientPointResponse evaluation = transient_point_response(point,
            ad_state,
            committed_state,
            data.material,
            committed_material[q],
            time_step,
            data.strain_formulation,
            data.time);
        const MaterialFunctionContext context = rz_material_context(data.time, point);
        const adlite::Scalar conductivity = data.material.conductivity(evaluation.temperature, context);
        adlite::Scalar temperature_rate = 0.0, heat_capacity = 0.0;
        if (include_thermal_time_term) {
            temperature_rate = (evaluation.temperature - evaluation.old_temperature) / time_step;
            heat_capacity = data.material.heat_capacity(evaluation.temperature, context);
        }
        quad4_rz_detail::add_point_residual(point,
            evaluation.gradient_temperature_r,
            evaluation.gradient_temperature_z,
            evaluation.kinematics,
            &heat_capacity,
            &temperature_rate,
            conductivity,
            data.volumetric_heat_source,
            evaluation.response.stress,
            ad_residual);
    }
    return quad4_rz_detail::values(ad_state, ad_residual);
}

Quad4MaterialHistory compute_quad4_rz_transient_update(const Quad4RzData& data,
    const Quad4RzGeometry& geometry,
    const LocalValues& converged_state,
    const LocalValues& committed_state,
    const Quad4MaterialHistory& committed_material,
    double time_step) {
    validate_time_step(time_step);
    if (data.element_formulation == RzElementFormulation::cax4rt)
        return rz::compute_cax4rt(data,
            geometry,
            converged_state,
            committed_state,
            &committed_material,
            time_step,
            false,
            false)
            .history;
    if (data.element_formulation == RzElementFormulation::cax4t)
        return elements::evaluate_cax4t({data.material,
                                            geometry,
                                            converged_state,
                                            committed_state,
                                            &committed_material,
                                            time_step,
                                            data.time,
                                            data.volumetric_heat_source,
                                            data.strain_formulation,
                                            false},
            false)
            .history;
    const LocalAdValues passive_state = quad4_rz_detail::ad_state(converged_state);
    Quad4MaterialHistory result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const TransientPointResponse evaluation = transient_point_response(geometry.points[q],
            passive_state,
            committed_state,
            data.material,
            committed_material[q],
            time_step,
            data.strain_formulation,
            data.time);
        if (!std::isfinite(evaluation.temperature.value()) || !(evaluation.temperature.value() > 0.0))
            throw std::domain_error("Transient Quad4 trial temperature must be finite and positive");
        result[q] = IsotropicThermoelasticMaterial::state_values(evaluation.response.trial_state);
        result[q].stress = {
            evaluation.response.stress.rr.value(),
            evaluation.response.stress.zz.value(),
            evaluation.response.stress.hoop.value(),
            evaluation.response.stress.rz.value(),
        };
    }
    return result;
}

namespace {
bool finite_point(const RzPoint& point) {
    return std::isfinite(point.r) && std::isfinite(point.z);
}

void validate_line(const Line2InterfaceSideCoordinates& coordinates, const char* name) {
    for (const RzPoint& point : coordinates) {
        if (!finite_point(point))
            throw std::invalid_argument(std::string(name) + " requires finite coordinates");
        if (!(point.r >= 0.0))
            throw std::invalid_argument(std::string(name) + " requires nonnegative radii");
    }
    const double dr = coordinates[1].r - coordinates[0].r, dz = coordinates[1].z - coordinates[0].z;
    if (!(std::hypot(dr, dz) > 0.0))
        throw std::invalid_argument(std::string(name) + " requires a nonzero line length");
}

void validate_edge(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes,
    const char* name) {
    if (local_nodes[0] >= 4 || local_nodes[1] >= 4 || local_nodes[0] == local_nodes[1])
        throw std::invalid_argument(std::string(name) + " requires two distinct Quad4 local nodes");
    validate_line(coordinates, name);
}

void displaced_edge_coordinates(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes,
    const LocalAdValues& state,
    bool use_displaced_geometry,
    std::array<adlite::Scalar, 2>& radius,
    std::array<adlite::Scalar, 2>& axial) {
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

Line2RzBoundaryGeometry make_line2_rz_boundary_geometry(const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes) {
    validate_edge(coordinates, local_nodes, "Boundary edge");
    return {coordinates, local_nodes};
}

namespace {
void compute_line2_rz_mechanical_residual_ad(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const LocalAdValues& state,
    LocalAdValues& residual) {
    residual.fill(adlite::Scalar(0.0));
    std::array<adlite::Scalar, 2> radius{};
    std::array<adlite::Scalar, 2> axial{};
    displaced_edge_coordinates(geometry.coordinates,
        geometry.local_nodes,
        state,
        data.use_displaced_geometry,
        radius,
        axial);
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

void compute_line2_rz_convection_residual_ad(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const LocalAdValues& state,
    LocalAdValues& residual) {
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

LocalResidual compute_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian) {
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
                   tangent_z = 0.0, tangential_traction = 0.0, tangential_force = 0.0, elastic_tangential_slip = 0.0,
                   total_tangential_slip = 0.0;
    bool sliding = false;
};

double reference_projection_fraction(const RzPoint& secondary,
    const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double tangent_r = primary_coordinates[1].r - primary_coordinates[0].r,
                 tangent_z = primary_coordinates[1].z - primary_coordinates[0].z,
                 length_squared = tangent_r * tangent_r + tangent_z * tangent_z;
    return ((secondary.r - primary_coordinates[0].r) * tangent_r + (secondary.z - primary_coordinates[0].z) * tangent_z)
           / length_squared;
}

double reference_normal_orientation(const RzPoint& secondary,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    double primary_fraction,
    double zero_gap_orientation_hint) {
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

adlite::Scalar interpolate(const std::array<double, line2_interface_side_node_count>& shape,
    const LocalAdValues& state,
    std::size_t offset) {
    adlite::Scalar value = 0.0;
    for (std::size_t node = 0; node < line2_interface_side_node_count; ++node)
        value += shape[node] * state[offset + node];
    return value;
}

bool projection_is_inside(double fraction, bool includes_second_endpoint, bool first = false, double extension = 0) {
    if (fraction < (first ? -extension : 0.0))
        return false;
    if (includes_second_endpoint)
        return fraction <= 1.0 + extension;
    return fraction < 1.0;
}

ContactProjectionValue project_to_current_line(double secondary_r,
    double secondary_z,
    double primary_r_0,
    double primary_z_0,
    double primary_r_1,
    double primary_z_1,
    double normal_orientation,
    bool primary_segment_is_first,
    bool includes_second_endpoint,
    double extension) {
    const double tangent_r = primary_r_1 - primary_r_0, tangent_z = primary_z_1 - primary_z_0,
                 length_squared = tangent_r * tangent_r + tangent_z * tangent_z;
    double fraction =
        ((secondary_r - primary_r_0) * tangent_r + (secondary_z - primary_z_0) * tangent_z) / length_squared;
    const bool projected =
        projection_is_inside(fraction, includes_second_endpoint, primary_segment_is_first, extension);
    if (!projected)
        return {};
    const double length = std::sqrt(length_squared), primary_r = primary_r_0 + fraction * tangent_r,
                 primary_z = primary_z_0 + fraction * tangent_z, normal_r = normal_orientation * tangent_z / length,
                 normal_z = -normal_orientation * tangent_r / length;
    return {true, (primary_r - secondary_r) * normal_r + (primary_z - secondary_z) * normal_z};
}

HeatAdQuadratureValue evaluate_heat_quadrature(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const Line2RzHeatQuadraturePoint& point,
    const LocalAdValues& state,
    const GapHeatProperties& properties,
    bool includes_second_endpoint,
    bool first) {
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
        ((secondary_radius - primary_radius_0) * tangent_r + (secondary_axial - primary_axial_0) * tangent_z)
        / (tangent_length * tangent_length);
    if (!projection_is_inside(primary_fraction.value(), includes_second_endpoint, first, point.nodal ? .1 : 0))
        return {};
    const adlite::Scalar primary_shape_0 = 1.0 - primary_fraction, primary_shape_1 = primary_fraction,
                         primary_radius = primary_shape_0 * primary_radius_0 + primary_shape_1 * primary_radius_1,
                         primary_axial = primary_shape_0 * primary_axial_0 + primary_shape_1 * primary_axial_1,
                         normal_r = point.normal_orientation * tangent_z / tangent_length,
                         normal_z = -point.normal_orientation * tangent_r / tangent_length;
    const adlite::Scalar gap =
        (primary_radius - secondary_radius) * normal_r + (primary_axial - secondary_axial) * normal_z;
    const adlite::Scalar secondary_temperature = interpolate(point.secondary_shape, state, 0),
                         primary_temperature = primary_shape_0 * state[2] + primary_shape_1 * state[3];
    adlite::Scalar conductance;
    if (properties.law == GapHeatConductanceLaw::gas_gap) {
        const adlite::Scalar thermal_gap = adlite::max(gap, adlite::Scalar(properties.minimum_gap));
        conductance = properties.gap_conductivity / thermal_gap;
    } else {
        const adlite::Scalar pressure = adlite::max(-properties.contact_penalty * gap, adlite::Scalar(0.0));
        const adlite::Scalar average_temperature = 0.5 * (secondary_temperature + primary_temperature);
        conductance = properties.conductance + properties.clearance_derivative * gap
                      + properties.pressure_derivative * pressure
                      + properties.temperature_derivative * (average_temperature - properties.reference_temperature);
        if (!std::isfinite(conductance.value()) || conductance.value() < 0.0)
            throw std::domain_error("Axisymmetric affine gap conductance must be finite and nonnegative");
    }
    const adlite::Scalar heat_flux = conductance * (secondary_temperature - primary_temperature),
                         dr_dxi = 0.5 * (secondary_radius_1 - secondary_radius_0),
                         dz_dxi = 0.5 * (secondary_axial_1 - secondary_axial_0),
                         surface_jacobian = adlite::sqrt(dr_dxi * dr_dxi + dz_dxi * dz_dxi),
                         area_radius = point.nodal ? ((1.0 + point.secondary_shape[0]) * secondary_radius_0
                                                         + (1.0 + point.secondary_shape[1]) * secondary_radius_1)
                                                         / 3.0
                                                   : secondary_radius,
                         weighted_measure = 2.0 * pi * area_radius * surface_jacobian * point.integration_weight;
    return {true, gap, heat_flux, weighted_measure, primary_shape_0, primary_shape_1};
}

void evaluate_friction(ContactAdValue& value,
    const NodeToLineRzContactGeometry& geometry,
    const LocalAdValues& state,
    const LocalValues& committed_state,
    const ContactPointHistory& history,
    const NormalContactProperties& properties) {
    value.total_tangential_slip = history.total_tangential_slip;
    if (properties.friction_coefficient == 0.0 || !(value.pressure.value() > 0.0))
        return;
    const std::size_t secondary = geometry.secondary_local_node;
    const adlite::Scalar secondary_increment_r = state[4 + secondary] - committed_state[4 + secondary],
                         secondary_increment_z = state[8 + secondary] - committed_state[8 + secondary];
    const adlite::Scalar primary_increment_r = value.primary_shape_0 * (state[6] - committed_state[6])
                                               + value.primary_shape_1 * (state[7] - committed_state[7]);
    const adlite::Scalar primary_increment_z = value.primary_shape_0 * (state[10] - committed_state[10])
                                               + value.primary_shape_1 * (state[11] - committed_state[11]);
    const adlite::Scalar tangential_increment = (secondary_increment_r - primary_increment_r) * value.tangent_r
                                                + (secondary_increment_z - primary_increment_z) * value.tangent_z;
    value.total_tangential_slip += tangential_increment;
    const adlite::Scalar trial_slip = history.elastic_tangential_slip + tangential_increment,
                         sliding_limit = properties.friction_coefficient * value.pressure,
                         stick_stiffness = properties.maximum_elastic_slip > 0.0
                                               ? sliding_limit / properties.maximum_elastic_slip
                                               : adlite::Scalar(properties.penalty),
                         trial_traction = properties.maximum_elastic_slip > 0.0 ? stick_stiffness * trial_slip
                                                                                : properties.penalty * trial_slip;
    const double trial_magnitude = std::abs(trial_traction.value());
    if (trial_magnitude < sliding_limit.value() || (trial_magnitude == sliding_limit.value() && !history.sliding)) {
        value.tangential_traction = trial_traction;
        value.elastic_tangential_slip = trial_slip;
    } else {
        const double direction = trial_traction.value() < 0.0 ? -1.0 : 1.0;
        value.tangential_traction = direction * sliding_limit;
        value.elastic_tangential_slip = properties.maximum_elastic_slip > 0.0
                                            ? value.tangential_traction / stick_stiffness
                                            : value.tangential_traction / properties.penalty;
        value.sliding = true;
    }
    value.tangential_force = value.tangential_traction * value.tributary_area;
}

ContactAdValue evaluate_contact(const NodeToLineRzContactGeometry& geometry,
    const LocalAdValues& state,
    const LocalValues& committed_state,
    const ContactPointHistory& history,
    const NormalContactProperties& properties) {
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
        ((secondary_radius - primary_radius_0) * tangent_r + (secondary_z - primary_z_0) * tangent_z)
        / (tangent_length * tangent_length);
    const bool projected = projection_is_inside(primary_fraction.value(),
        geometry.primary_segment_includes_second_endpoint,
        geometry.primary_segment_is_first,
        .1);
    if (!projected)
        return {};
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
    const Line2InterfaceSideCoordinates& primary_coordinates,
    double secondary_coordinate_lower,
    double secondary_coordinate_upper,
    double zero_gap_orientation_hint) {
    if (!std::isfinite(secondary_coordinate_lower) || !std::isfinite(secondary_coordinate_upper)
        || secondary_coordinate_lower < -1.0 || secondary_coordinate_upper > 1.0
        || !(secondary_coordinate_upper > secondary_coordinate_lower))
        throw std::invalid_argument("Heat quadrature requires a nonempty secondary interval in [-1,1]");
    const std::array<double, line2_interface_quadrature_point_count> locations = {-gauss, gauss};
    std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> result{};
    for (std::size_t q = 0; q < result.size(); ++q) {
        const double secondary_xi =
            0.5
            * ((1.0 - locations[q]) * secondary_coordinate_lower + (1.0 + locations[q]) * secondary_coordinate_upper);
        const std::array<double, 2> secondary_shape = {
            0.5 * (1.0 - secondary_xi),
            0.5 * (1.0 + secondary_xi),
        };
        result[q] = make_line2_rz_heat_point_geometry(secondary_coordinates,
            primary_coordinates,
            secondary_shape,
            0.5 * (secondary_coordinate_upper - secondary_coordinate_lower),
            true,
            zero_gap_orientation_hint)
                        .point;
        const double primary_fraction = result[q].primary_shape[1];
        if (primary_fraction < 0.0 || primary_fraction > 1.0)
            throw std::invalid_argument("Heat-contact Gauss point lies outside the primary segment");
    }
    return result;
}

Line2RzHeatPointGeometry make_line2_rz_heat_point_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const std::array<double, line2_interface_side_node_count>& secondary_shape,
    double integration_weight,
    bool primary_segment_includes_second_endpoint,
    double zero_gap_orientation_hint) {
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
        {secondary_shape,
            {1.0 - primary_fraction, primary_fraction},
            integration_weight,
            reference_normal_orientation(secondary_point,
                primary_coordinates,
                closest_fraction,
                zero_gap_orientation_hint)},
        primary_segment_includes_second_endpoint,
    };
}

LocalResidual compute_line2_rz_gap_heat(const GapHeatProperties& properties,
    const Line2RzHeatPointGeometry& geometry,
    const LocalValues& state,
    LocalJacobian* jacobian) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state, jacobian != nullptr);
    LocalAdValues ad_residual{};
    ad_residual.fill(adlite::Scalar(0.0));
    const HeatAdQuadratureValue value = evaluate_heat_quadrature(geometry.secondary_coordinates,
        geometry.primary_coordinates,
        geometry.point,
        ad_state,
        properties,
        geometry.primary_segment_includes_second_endpoint,
        geometry.primary_segment_is_first);
    if (value.projected) {
        for (std::size_t node = 0; node < line2_interface_side_node_count; ++node) {
            ad_residual[node] += value.weighted_measure * geometry.point.secondary_shape[node] * value.heat_flux;
            const adlite::Scalar primary_shape = node == 0 ? value.primary_shape_0 : value.primary_shape_1;
            ad_residual[2 + node] -= value.weighted_measure * primary_shape * value.heat_flux;
        }
    }
    return quad4_rz_detail::values(ad_state, ad_residual, jacobian);
}

HeatQuadratureValue compute_line2_rz_gap_heat_value(const GapHeatProperties& properties,
    const Line2RzHeatPointGeometry& geometry,
    const LocalValues& state) {
    const LocalAdValues ad_state = quad4_rz_detail::ad_state(state);
    const HeatAdQuadratureValue value = evaluate_heat_quadrature(geometry.secondary_coordinates,
        geometry.primary_coordinates,
        geometry.point,
        ad_state,
        properties,
        geometry.primary_segment_includes_second_endpoint,
        geometry.primary_segment_is_first);
    return {value.projected, value.gap.value(), value.heat_flux.value(), value.weighted_measure.value()};
}

ContactProjectionValue compute_line2_rz_heat_projection(const Line2RzHeatPointGeometry& geometry,
    const LocalValues& state) {
    const double secondary_r_0 = geometry.secondary_coordinates[0].r + state[4],
                 secondary_r_1 = geometry.secondary_coordinates[1].r + state[5],
                 secondary_z_0 = geometry.secondary_coordinates[0].z + state[8],
                 secondary_z_1 = geometry.secondary_coordinates[1].z + state[9],
                 secondary_r = geometry.point.secondary_shape[0] * secondary_r_0
                               + geometry.point.secondary_shape[1] * secondary_r_1,
                 secondary_z = geometry.point.secondary_shape[0] * secondary_z_0
                               + geometry.point.secondary_shape[1] * secondary_z_1;
    return project_to_current_line(secondary_r,
        secondary_z,
        geometry.primary_coordinates[0].r + state[6],
        geometry.primary_coordinates[0].z + state[10],
        geometry.primary_coordinates[1].r + state[7],
        geometry.primary_coordinates[1].z + state[11],
        geometry.point.normal_orientation,
        geometry.primary_segment_is_first,
        geometry.primary_segment_includes_second_endpoint,
        geometry.point.nodal ? .1 : 0);
}

NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates,
    std::size_t secondary_local_node,
    bool primary_segment_is_first,
    bool primary_segment_includes_upper_endpoint,
    double zero_gap_orientation_hint) {
    validate_line(secondary_edge_coordinates, "NodeToLineRzContactGeometry secondary");
    validate_line(primary_segment_coordinates, "NodeToLineRzContactGeometry primary");
    if (secondary_local_node >= line2_interface_side_node_count)
        throw std::invalid_argument("NodeToLineRzContactGeometry secondary node is out of range");
    const RzPoint secondary_point = secondary_edge_coordinates[secondary_local_node];
    const double primary_fraction = reference_projection_fraction(secondary_point, primary_segment_coordinates),
                 closest_fraction = std::max(0.0, std::min(1.0, primary_fraction));
    const double orientation = reference_normal_orientation(secondary_point,
        primary_segment_coordinates,
        closest_fraction,
        zero_gap_orientation_hint);
    return {
        secondary_edge_coordinates,
        primary_segment_coordinates,
        secondary_local_node,
        primary_segment_is_first,
        primary_segment_includes_upper_endpoint,
        orientation,
    };
}

LocalResidual compute_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry,
    const LocalValues& state,
    const LocalValues& committed_state,
    const ContactPointHistory& history,
    LocalJacobian* jacobian) {
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
    const NodeToLineRzContactGeometry& geometry,
    const LocalValues& state,
    const LocalValues& committed_state,
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
        result.total_tangential_slip.value(),
    };
}

ContactProjectionValue compute_node_to_line_rz_contact_projection(const NodeToLineRzContactGeometry& geometry,
    const LocalValues& state) {
    const std::size_t secondary = geometry.secondary_local_node;
    return project_to_current_line(geometry.secondary_edge_coordinates[secondary].r + state[4 + secondary],
        geometry.secondary_edge_coordinates[secondary].z + state[8 + secondary],
        geometry.primary_segment_coordinates[0].r + state[6],
        geometry.primary_segment_coordinates[0].z + state[10],
        geometry.primary_segment_coordinates[1].r + state[7],
        geometry.primary_segment_coordinates[1].z + state[11],
        geometry.normal_orientation,
        geometry.primary_segment_is_first,
        geometry.primary_segment_includes_second_endpoint,
        .1);
}
} // namespace fuelsim
