#include "cax4t.hpp"
#include "detail/rz_point.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

struct AxisymmetricKinematics final {
    std::array<adlite::Scalar, quad4_node_count> gradient_r, gradient_z;
    adlite::Scalar radius, weighted_measure, strain_rr, strain_zz, strain_hoop, strain_rz;
    AxisymmetricRotation rotation;
    adlite::Scalar midpoint_weighted_measure = 0.0;
};

AxisymmetricKinematics evaluate_axisymmetric_kinematics_from_point(const RzQuadraturePoint& point,
    const adlite::Scalar& radial_displacement,
    const adlite::Scalar& displacement_gradient_rr,
    const adlite::Scalar& displacement_gradient_rz,
    const adlite::Scalar& displacement_gradient_zr,
    const adlite::Scalar& displacement_gradient_zz,
    const Cax4LocalValues& committed_state,
    StrainFormulation strain_formulation);

} // namespace
} // namespace fuelsim

namespace fuelsim::elements {
namespace {
// CAX4T uses point-local width-six kinematics and width-five constitutive AD.
// Cross-point volume/pressure dependencies are assembled through explicit nodal chains.
Cax4LocalValues cax4t_nodal_chain(const adlite::Scalar& value, const RzQuadraturePoint& point) {
    std::array<double, 6> d{};
    value.copy_derivatives(d.data(), d.size());
    Cax4LocalValues result{};
    for (std::size_t n = 0; n < 4; ++n) {
        result[n] = d[5] * point.shape[n];
        result[4 + n] = d[0] * point.gradient_r[n] + d[1] * point.gradient_z[n] + d[4] * point.shape[n];
        result[8 + n] = d[2] * point.gradient_r[n] + d[3] * point.gradient_z[n];
    }
    return result;
}

struct Cax4tPointSystem final {
    AxisymmetricKinematics kinematics;
    adlite::Scalar temperature;
    AxisymmetricStress stress;
    std::array<Cax4LocalValues, 4> stress_derivatives{};
    MaterialPointState history;
};

} // namespace

Cax4Result evaluate_cax4t(const Cax4Input& input, ElementRequest request) {
    const bool jacobian = request.jacobian;
    const auto& data = input;
    const auto& geometry = input.geometry;
    const auto& state = input.state;
    const auto& old_state = input.committed_state;
    const auto* old_history = input.committed_history;
    const double time_step = input.time_step;
    const bool thermal_time = input.include_thermal_time_term;
    if (!std::isfinite(input.time) || !std::isfinite(input.volumetric_heat_source))
        throw std::invalid_argument("CAX4T time and heat source must be finite");
    if (old_history && (!std::isfinite(time_step) || !(time_step > 0.0)))
        throw std::invalid_argument("CAX4T transient time_step must be finite and positive");
    if (!old_history && thermal_time)
        throw std::invalid_argument("CAX4T heat capacity requires committed material history");
    for (std::size_t n = 0; n < 12; ++n) {
        if (!std::isfinite(state[n]) || (n < 4 && !(state[n] > 0.0)))
            throw std::domain_error("CAX4T trial values must be finite and temperatures positive");
        if (!std::isfinite(old_state[n]) || (old_history && n < 4 && !(old_state[n] > 0.0)))
            throw std::invalid_argument("CAX4T committed values must be finite and temperatures positive");
    }
    std::array<Cax4tPointSystem, 4> systems;
    double numerator = 0.0, midpoint_volume = 0.0, current_volume = 0.0, reference_volume = 0.0;
    Cax4LocalValues numerator_derivatives{}, midpoint_derivatives{}, volume_derivatives{}, hoop_shape{};
    double hoop_new = 0.0, hoop_old = 0.0;
    for (std::size_t q = 0; q < 4; ++q) {
        const auto& point = geometry.points[q];
        auto& s = systems[q];
        const std::array<double, 6> seeds = {quad4_rz_detail::interpolate(point.gradient_r, state, 4),
            quad4_rz_detail::interpolate(point.gradient_z, state, 4),
            quad4_rz_detail::interpolate(point.gradient_r, state, 8),
            quad4_rz_detail::interpolate(point.gradient_z, state, 8),
            quad4_rz_detail::interpolate(point.shape, state, 4),
            quad4_rz_detail::interpolate(point.shape, state, 0)};
        std::array<adlite::Scalar, 6> active;
        for (std::size_t i = 0; i < 6; ++i)
            active[i] = jacobian ? adlite::Scalar::independent(seeds[i], i, 6) : adlite::Scalar(seeds[i]);
        s.temperature = active[5];
        s.kinematics = evaluate_axisymmetric_kinematics_from_point(point,
            active[4],
            active[0],
            active[1],
            active[2],
            active[3],
            old_state,
            data.strain_formulation);
        const auto& k = s.kinematics;
        const adlite::Scalar weighted_trace = k.midpoint_weighted_measure * (k.strain_rr + k.strain_zz + k.strain_hoop);
        numerator += weighted_trace.value();
        midpoint_volume += k.midpoint_weighted_measure.value();
        current_volume += k.weighted_measure.value();
        reference_volume += point.weighted_measure;
        hoop_new += point.weighted_measure * (1.0 + seeds[4] / point.radius);
        hoop_old +=
            point.weighted_measure * (1.0 + quad4_rz_detail::interpolate(point.shape, old_state, 4) / point.radius);
        for (std::size_t n = 0; n < 4; ++n)
            hoop_shape[4 + n] += point.weighted_measure * point.shape[n] / point.radius;
        if (jacobian) {
            const auto a = cax4t_nodal_chain(weighted_trace, point);
            const auto b = cax4t_nodal_chain(k.midpoint_weighted_measure, point);
            const auto c = cax4t_nodal_chain(k.weighted_measure, point);
            for (std::size_t j = 0; j < 12; ++j) {
                numerator_derivatives[j] += a[j];
                midpoint_derivatives[j] += b[j];
                volume_derivatives[j] += c[j];
            }
        }
    }
    if (!std::isfinite(midpoint_volume) || !std::isfinite(current_volume) || !std::isfinite(reference_volume)
        || !(midpoint_volume > 0.0) || !(current_volume > 0.0) || !(reference_volume > 0.0))
        throw std::domain_error("CAX4T requires positive reference, midpoint and current volumes");
    const double average_trace = numerator / midpoint_volume;
    Cax4LocalValues average_derivatives{};
    for (std::size_t j = 0; j < 12; ++j)
        average_derivatives[j] = (numerator_derivatives[j] - average_trace * midpoint_derivatives[j]) / midpoint_volume;
    hoop_new /= reference_volume;
    hoop_old /= reference_volume;
    for (double& d : hoop_shape)
        d /= reference_volume;
    const bool finite_strain = data.strain_formulation == StrainFormulation::finite;
    // Abaqus averages F_hoop over reference volume before forming its central increment.
    // The independently averaged full trace then supplies only the in-plane correction.
    const double average_hoop = finite_strain ? 2.0 * (hoop_new - hoop_old) / (hoop_new + hoop_old) : hoop_new - 1.0;
    Cax4LocalValues hoop_derivatives{};
    for (std::size_t j = 0; j < 12; ++j)
        hoop_derivatives[j] =
            hoop_shape[j] * (finite_strain ? 4.0 * hoop_old / ((hoop_new + hoop_old) * (hoop_new + hoop_old)) : 1.0);
    double pressure = 0.0, hoop_stress = 0.0;
    Cax4LocalValues pressure_derivatives{}, hoop_stress_derivatives{};
    Cax4Result result;
    for (std::size_t q = 0; q < 4; ++q) {
        const auto& point = geometry.points[q];
        auto& s = systems[q];
        auto& k = s.kinematics;
        const adlite::Scalar correction = (average_trace - average_hoop - k.strain_rr - k.strain_zz) / 2.0;
        k.strain_rr += correction;
        k.strain_zz += correction;
        k.strain_hoop = average_hoop;
        std::array<adlite::Scalar, 5> inputs = {k.strain_rr, k.strain_zz, k.strain_hoop, k.strain_rz, s.temperature};
        const MaterialFunctionContext context = {data.time, point.radius, 0.0, point.axial_coordinate};
        // CAX4T uses the arithmetic corner temperature for expansion, while
        // constitutive properties retain the point temperature. CAX4RT differs.
        const double center_temperature = (state[0] + state[1] + state[2] + state[3]) / 4.0;
        const auto point_eigen = data.material.eigenstrain_rz(s.temperature, context);
        const auto center_eigen = data.material.eigenstrain_rz(
            jacobian ? adlite::Scalar::independent(center_temperature, 0, 1) : adlite::Scalar(center_temperature),
            context);
        const std::array<adlite::Scalar, 4> pe = {point_eigen.rr, point_eigen.zz, point_eigen.hoop, point_eigen.rz};
        const std::array<adlite::Scalar, 4> ce = {center_eigen.rr, center_eigen.zz, center_eigen.hoop, center_eigen.rz};
        std::array<double, 4> center_derivative{};
        for (std::size_t i = 0; i < 4; ++i) {
            inputs[i] += pe[i] - ce[i].value();
            if (jacobian)
                ce[i].copy_derivatives(&center_derivative[i], 1);
        }
        if (old_history && finite_strain) {
            auto old_context = context;
            old_context.time -= time_step;
            const auto op =
                data.material.eigenstrain_rz(adlite::Scalar(quad4_rz_detail::interpolate(point.shape, old_state, 0)),
                    old_context);
            const auto oc = data.material.eigenstrain_rz(
                adlite::Scalar((old_state[0] + old_state[1] + old_state[2] + old_state[3]) / 4.0),
                old_context);
            const std::array<double, 4> difference = {op.rr.value() - oc.rr.value(),
                op.zz.value() - oc.zz.value(),
                op.hoop.value() - oc.hoop.value(),
                op.rz.value() - oc.rz.value()};
            for (std::size_t i = 0; i < 4; ++i)
                inputs[i] -= difference[i];
        }
        std::array<double, 4> fed = {inputs[0].value(), inputs[1].value(), inputs[2].value(), inputs[3].value()};
        if (old_history && data.strain_formulation == StrainFormulation::finite) {
            auto old_context = context;
            old_context.time -= time_step;
            const auto eigen =
                data.material.eigenstrain_rz(adlite::Scalar(quad4_rz_detail::interpolate(point.shape, old_state, 0)),
                    old_context);
            const std::array<double, 4> imposed = {eigen.rr.value(),
                eigen.zz.value(),
                eigen.hoop.value(),
                eigen.rz.value()};
            for (std::size_t i = 0; i < 4; ++i)
                fed[i] = (*old_history)[q].elastic_strain[i] + inputs[i].value() + imposed[i]
                         + (*old_history)[q].plastic_strain[i] + (*old_history)[q].creep_strain[i];
        }
        AxisymmetricStressTangent tangent;
        if (jacobian) {
            tangent = evaluate_axisymmetric_stress_tangent(data.material,
                fed,
                s.temperature.value(),
                time_step,
                old_history ? &(*old_history)[q] : nullptr,
                context);
        } else {
            const auto sigma =
                old_history ? data.material
                                  .response(fed[0],
                                      fed[1],
                                      fed[2],
                                      fed[3],
                                      s.temperature.value(),
                                      time_step,
                                      (*old_history)[q],
                                      context)
                                  .stress
                            : data.material.stress(fed[0], fed[1], fed[2], fed[3], s.temperature.value(), context);
            tangent.stress = {sigma.rr.value(), sigma.zz.value(), sigma.hoop.value(), sigma.rz.value()};
        }
        const std::array<double, 4> values = {tangent.stress.rr,
            tangent.stress.zz,
            tangent.stress.hoop,
            tangent.stress.rz};
        std::array<adlite::Scalar, 4> sigma;
        std::array<double, 4> trace_response{}, hoop_response{}, center_response{};
        for (std::size_t i = 0; i < 4; ++i) {
            std::array<double, 5> partials{};
            for (std::size_t j = 0; j < 4; ++j)
                partials[j] = tangent.tangent[i][j];
            partials[4] = tangent.thermal[i];
            sigma[i] =
                jacobian ? adlite::compose(values[i], inputs.data(), partials.data(), 5) : adlite::Scalar(values[i]);
            trace_response[i] = (partials[0] + partials[1]) / 2.0;
            hoop_response[i] = partials[2] - trace_response[i];
            for (std::size_t j = 0; j < 4; ++j)
                center_response[i] -= partials[j] * center_derivative[j];
        }
        s.stress = {sigma[0], sigma[1], sigma[2], sigma[3]};
        AxisymmetricStress trace_stress = {trace_response[0], trace_response[1], trace_response[2], trace_response[3]};
        AxisymmetricStress hoop_tangent = {hoop_response[0], hoop_response[1], hoop_response[2], hoop_response[3]};
        AxisymmetricStress center_tangent = {center_response[0],
            center_response[1],
            center_response[2],
            center_response[3]};
        if (data.strain_formulation == StrainFormulation::finite) {
            s.stress = rotate_axisymmetric_tensor(s.stress, k.rotation);
            trace_stress = rotate_axisymmetric_tensor(trace_stress, k.rotation);
            hoop_tangent = rotate_axisymmetric_tensor(hoop_tangent, k.rotation);
            center_tangent = rotate_axisymmetric_tensor(center_tangent, k.rotation);
        }
        const std::array<adlite::Scalar, 4> final_sigma = {s.stress.rr, s.stress.zz, s.stress.hoop, s.stress.rz};
        trace_response = {trace_stress.rr.value(),
            trace_stress.zz.value(),
            trace_stress.hoop.value(),
            trace_stress.rz.value()};
        hoop_response = {hoop_tangent.rr.value(),
            hoop_tangent.zz.value(),
            hoop_tangent.hoop.value(),
            hoop_tangent.rz.value()};
        center_response = {center_tangent.rr.value(),
            center_tangent.zz.value(),
            center_tangent.hoop.value(),
            center_tangent.rz.value()};
        if (jacobian)
            for (std::size_t i = 0; i < 4; ++i) {
                s.stress_derivatives[i] = cax4t_nodal_chain(final_sigma[i], point);
                for (std::size_t j = 0; j < 12; ++j)
                    s.stress_derivatives[i][j] +=
                        trace_response[i] * average_derivatives[j] + hoop_response[i] * hoop_derivatives[j];
                for (std::size_t j = 0; j < 4; ++j)
                    s.stress_derivatives[i][j] += center_response[i] / 4.0;
            }
        const double weight = point.weighted_measure / reference_volume;
        pressure += weight * (s.stress.rr.value() + s.stress.zz.value()) / 2.0;
        hoop_stress += weight * s.stress.hoop.value();
        for (std::size_t j = 0; j < 12; ++j) {
            pressure_derivatives[j] += weight * (s.stress_derivatives[0][j] + s.stress_derivatives[1][j]) / 2.0;
            hoop_stress_derivatives[j] += weight * s.stress_derivatives[2][j];
        }
        if (old_history) {
            const AxisymmetricRotation rotation = {k.rotation.rr.value(),
                k.rotation.rz.value(),
                k.rotation.zr.value(),
                k.rotation.zz.value(),
                k.rotation.hoop.value()};
            const auto response = data.strain_formulation == StrainFormulation::finite
                                      ? data.material.incremental_response(inputs[0].value(),
                                            inputs[1].value(),
                                            inputs[2].value(),
                                            inputs[3].value(),
                                            rotation,
                                            s.temperature.value(),
                                            quad4_rz_detail::interpolate(point.shape, old_state, 0),
                                            time_step,
                                            (*old_history)[q],
                                            context)
                                      : data.material.response(fed[0],
                                            fed[1],
                                            fed[2],
                                            fed[3],
                                            s.temperature.value(),
                                            time_step,
                                            (*old_history)[q],
                                            context);
            s.history = IsotropicThermoelasticMaterial::state_values(response.trial_state);
        }
        s.history.stress = {s.stress.rr.value(), s.stress.zz.value(), s.stress.hoop.value(), s.stress.rz.value()};
        result.history[q] = s.history;
    }
    for (std::size_t q = 0; q < 4; ++q) {
        const auto& point = geometry.points[q];
        const auto& s = systems[q];
        const auto& k = s.kinematics;
        const double w = k.weighted_measure.value();
        const double point_pressure = (s.stress.rr.value() + s.stress.zz.value()) / 2.0;
        const bool finite = data.strain_formulation == StrainFormulation::finite;
        const double scale = finite ? point.weighted_measure * current_volume / (reference_volume * w) : 1.0;
        const std::array<double, 4> deviator = {s.stress.rr.value() - point_pressure,
            s.stress.zz.value() - point_pressure,
            0.0,
            s.stress.rz.value()};
        const std::array<double, 4> stress = {scale * deviator[0] + pressure,
            scale * deviator[1] + pressure,
            pressure,
            scale * deviator[3]};
        const auto wd = jacobian ? cax4t_nodal_chain(k.weighted_measure, point) : Cax4LocalValues{};
        std::array<Cax4LocalValues, 4> sd{};
        if (jacobian)
            for (std::size_t j = 0; j < 12; ++j) {
                const double pd = (s.stress_derivatives[0][j] + s.stress_derivatives[1][j]) / 2.0;
                const double ds = finite ? scale * (volume_derivatives[j] / current_volume - wd[j] / w) : 0.0;
                for (std::size_t i = 0; i < 4; ++i)
                    sd[i][j] = i == 2 ? pressure_derivatives[j]
                                      : ds * deviator[i] + scale * (s.stress_derivatives[i][j] - (i < 2 ? pd : 0.0))
                                            + (i < 2 ? pressure_derivatives[j] : 0.0);
            }
        const MaterialFunctionContext context = {data.time, point.radius, 0.0, point.axial_coordinate};
        // CAX4T uses the temperature corner associated with this Gauss station
        // for conductivity; mechanical properties retain point temperature.
        const auto conductivity = data.material.conductivity(jacobian ? adlite::Scalar::independent(state[q], 0, 1)
                                                                      : adlite::Scalar(state[q]),
            context);
        double conductivity_derivative = 0.0;
        if (jacobian)
            conductivity.copy_derivatives(&conductivity_derivative, 1);
        std::array<adlite::Scalar, 4> thermal_gr{}, thermal_gz{};
        adlite::Scalar fa = 1.0, fb = 0.0, fc = 0.0, fd = 1.0;
        if (finite) {
            std::array<adlite::Scalar, 4> gradient{};
            for (std::size_t component = 0; component < 4; ++component) {
                const auto& coefficients = component % 2 == 0 ? point.gradient_r : point.gradient_z;
                const std::size_t offset = component < 2 ? 4 : 8;
                const double current = quad4_rz_detail::interpolate(coefficients, state, offset);
                gradient[component] =
                    jacobian ? adlite::Scalar::independent(current, component, 6) : adlite::Scalar(current);
                gradient[component] =
                    .5 * (gradient[component] + quad4_rz_detail::interpolate(coefficients, old_state, offset));
            }
            fa += gradient[0];
            fb = gradient[1];
            fc = gradient[2];
            fd += gradient[3];
        }
        const auto thermal_determinant = fa * fd - fb * fc;
        if (!(thermal_determinant.value() > 0))
            throw std::domain_error("CAX4T thermal midpoint must remain positive");
        adlite::Scalar tr = 0.0, tz = 0.0;
        for (std::size_t n = 0; n < 4; ++n) {
            thermal_gr[n] = (fd * point.gradient_r[n] - fc * point.gradient_z[n]) / thermal_determinant;
            thermal_gz[n] = (fa * point.gradient_z[n] - fb * point.gradient_r[n]) / thermal_determinant;
            tr += thermal_gr[n] * state[n];
            tz += thermal_gz[n] * state[n];
        }
        const double conduction_measure = point.weighted_measure * (finite ? current_volume / reference_volume : 1.0);
        const auto thermal_measure = finite ? k.weighted_measure : adlite::Scalar(point.weighted_measure);
        result.generated_heat_rate += thermal_measure.value() * data.volumetric_heat_source;
        for (std::size_t n = 0; n < 4; ++n) {
            const adlite::Scalar hoop = point.shape[n] / k.radius;
            const std::array<adlite::Scalar, 4> radial = {k.gradient_r[n], adlite::Scalar(0.0), hoop, k.gradient_z[n]};
            const std::array<adlite::Scalar, 4> axial = {adlite::Scalar(0.0),
                k.gradient_z[n],
                adlite::Scalar(0.0),
                k.gradient_r[n]};
            for (std::size_t equation = 0; equation < 2; ++equation) {
                const auto& b = equation == 0 ? radial : axial;
                const auto row = 4 * (equation + 1) + n;
                double contraction = 0.0;
                for (std::size_t i = 0; i < 4; ++i)
                    contraction += b[i].value() * stress[i];
                result.residual[row] += w * contraction;
                if (jacobian) {
                    for (std::size_t j = 0; j < 12; ++j)
                        result.jacobian[row * 12 + j] += wd[j] * contraction;
                    for (std::size_t i = 0; i < 4; ++i) {
                        const auto bd = cax4t_nodal_chain(b[i], point);
                        for (std::size_t j = 0; j < 12; ++j)
                            result.jacobian[row * 12 + j] += w * (bd[j] * stress[i] + b[i].value() * sd[i][j]);
                    }
                }
            }
            const auto conduction = thermal_gr[n] * tr + thermal_gz[n] * tz;
            auto thermal = conduction_measure * conductivity.value() * conduction
                           - thermal_measure * data.volumetric_heat_source * point.shape[n];
            // CAX4T capacity is row-sum lumped: integrate N_i over the
            // selected volume, then evaluate rho*cp and temperature rate at i.
            if (thermal_time && old_history) {
                const auto temperature =
                    jacobian ? adlite::Scalar::independent(state[n], 0, 1) : adlite::Scalar(state[n]);
                const auto& coordinate = geometry.coordinates[n];
                const auto rate = data.material.heat_capacity(temperature, {data.time, coordinate.r, 0.0, coordinate.z})
                                  * (temperature - old_state[n]) / time_step;
                thermal += thermal_measure * point.shape[n] * rate.value();
                result.stored_heat_rate += thermal_measure.value() * point.shape[n] * rate.value();
                if (jacobian) {
                    double derivative = 0.0;
                    rate.copy_derivatives(&derivative, 1);
                    result.jacobian[n * 12 + n] += thermal_measure.value() * point.shape[n] * derivative;
                }
            }
            result.residual[n] += thermal.value();
            if (jacobian) {
                const auto td = cax4t_nodal_chain(thermal, point);
                for (std::size_t j = 0; j < 12; ++j)
                    result.jacobian[n * 12 + j] += td[j];
                result.jacobian[n * 12 + q] += conduction_measure * conductivity_derivative * conduction.value();
                if (finite)
                    for (std::size_t j = 0; j < 12; ++j)
                        result.jacobian[n * 12 + j] += point.weighted_measure / reference_volume * volume_derivatives[j]
                                                       * conductivity.value() * conduction.value();
                for (std::size_t j = 0; j < 4; ++j)
                    result.jacobian[n * 12 + j] += conduction_measure * conductivity.value()
                                                   * (thermal_gr[n].value() * thermal_gr[j].value()
                                                       + thermal_gz[n].value() * thermal_gz[j].value());
            }
        }
    }
    // The virtual-work hoop term is p_plane*N/r plus
    // (average_sigma_hoop-p_plane)*average_reference(N/R)/average_F_hoop.
    // Its nonlocal part cannot be assembled with the unmodified pointwise hoop B entry.
    const double hoop_factor = current_volume / (finite_strain ? hoop_new : 1.0);
    for (std::size_t n = 0; n < 4; ++n) {
        const auto row = 4 + n;
        result.residual[row] += hoop_factor * (hoop_stress - pressure) * hoop_shape[row];
        if (jacobian)
            for (std::size_t j = 0; j < 12; ++j) {
                const double factor_derivative =
                    finite_strain
                        ? volume_derivatives[j] / hoop_new - current_volume * hoop_shape[j] / (hoop_new * hoop_new)
                        : 0.0;
                result.jacobian[row * 12 + j] +=
                    hoop_shape[row]
                    * (factor_derivative * (hoop_stress - pressure)
                        + hoop_factor * (hoop_stress_derivatives[j] - pressure_derivatives[j]));
            }
    }
    if (request.stress)
        for (std::size_t q = 0; q < result.history.size(); ++q)
            result.stress[q] = result.history[q].stress;
    if (!request.residual && !request.jacobian)
        result.residual.fill(0.0);
    if (!request.history)
        result.history = {};
    return result;
}
} // namespace fuelsim::elements

namespace fuelsim {
namespace {
AxisymmetricKinematics evaluate_axisymmetric_kinematics_from_point(const RzQuadraturePoint& point,
    const adlite::Scalar& radial_displacement,
    const adlite::Scalar& displacement_gradient_rr,
    const adlite::Scalar& displacement_gradient_rz,
    const adlite::Scalar& displacement_gradient_zr,
    const adlite::Scalar& displacement_gradient_zz,
    const Cax4LocalValues& committed_state,
    StrainFormulation strain_formulation) {
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
    {
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
}

} // namespace

} // namespace fuelsim
