#include "cax4rt.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
using Gradient = std::array<std::array<double, 2>, 4>;
using GradientDerivative = std::array<std::array<Cax4LocalValues, 2>, 4>;
constexpr std::array<double, 4> mode = {1.0, -1.0, 1.0, -1.0};
constexpr std::array<std::array<double, 2>, 4> signs = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};

struct ReducedGeometry final {
    double volume = 0.0;
    Cax4LocalValues volume_derivative{};
    Gradient gradient{};
    GradientDerivative gradient_derivative{};
    std::array<double, 4> hoop{}, measures{}, gamma{};
    std::array<Cax4LocalValues, 4> hoop_derivative{}, measure_derivative{}, gamma_derivative{};
    double thermal_coefficient = 0.0;
    Cax4LocalValues thermal_derivative{};
};

ReducedGeometry reduce_geometry(const Quad4RzGeometry& reference, const Cax4LocalValues& state, double chain_scale) {
    ReducedGeometry result;
    for (const auto& p : reference.points) {
        double frr = 1.0, frz = 0.0, fzr = 0.0, fzz = 1.0, radius = p.radius;
        for (std::size_t n = 0; n < 4; ++n) {
            frr += state[4 + n] * p.gradient_r[n];
            frz += state[4 + n] * p.gradient_z[n];
            fzr += state[8 + n] * p.gradient_r[n];
            fzz += state[8 + n] * p.gradient_z[n];
            radius += state[4 + n] * p.shape[n];
        }
        const double det = frr * fzz - frz * fzr;
        if (!std::isfinite(det) || !(det > 0.0) || !std::isfinite(radius) || !(radius > 0.0))
            throw std::domain_error("CAX4RT reference, committed, midpoint and current geometry must remain positive");
        const double w = p.weighted_measure * det * radius / p.radius;
        Gradient b{};
        for (std::size_t n = 0; n < 4; ++n)
            b[n] = {(fzz * p.gradient_r[n] - fzr * p.gradient_z[n]) / det,
                (-frz * p.gradient_r[n] + frr * p.gradient_z[n]) / det};
        result.volume += w;
        for (std::size_t n = 0; n < 4; ++n) {
            result.measures[n] += w * p.shape[n];
            result.hoop[n] += w * p.shape[n] / radius;
            for (std::size_t d = 0; d < 2; ++d)
                result.gradient[n][d] += w * b[n][d];
        }
        if (chain_scale == 0.0)
            continue;
        for (std::size_t j = 0; j < 8; ++j) {
            const auto a = j % 4, c = j / 4, column = 4 + j;
            const double dr = c == 0 ? p.shape[a] : 0.0;
            const double dw = chain_scale * w * (b[a][c] + dr / radius);
            result.volume_derivative[column] += dw;
            for (std::size_t n = 0; n < 4; ++n) {
                result.measure_derivative[n][column] += dw * p.shape[n];
                result.hoop_derivative[n][column] +=
                    dw * p.shape[n] / radius - chain_scale * w * p.shape[n] * dr / (radius * radius);
                for (std::size_t d = 0; d < 2; ++d)
                    result.gradient_derivative[n][d][column] += dw * b[n][d] - chain_scale * w * b[a][d] * b[n][c];
            }
        }
    }
    for (std::size_t n = 0; n < 4; ++n) {
        result.hoop[n] /= result.volume;
        for (std::size_t d = 0; d < 2; ++d)
            result.gradient[n][d] /= result.volume;
        for (std::size_t j = 0; j < 12; ++j) {
            result.hoop_derivative[n][j] =
                (result.hoop_derivative[n][j] - result.hoop[n] * result.volume_derivative[j]) / result.volume;
            for (std::size_t d = 0; d < 2; ++d)
                result.gradient_derivative[n][d][j] =
                    (result.gradient_derivative[n][d][j] - result.gradient[n][d] * result.volume_derivative[j])
                    / result.volume;
        }
    }
    std::array<double, 2> projected{}, center{};
    double jrr = 0.0, jrz = 0.0, jzr = 0.0, jzz = 0.0;
    for (std::size_t n = 0; n < 4; ++n) {
        const double r = reference.coordinates[n].r + state[4 + n], z = reference.coordinates[n].z + state[8 + n];
        projected[0] += mode[n] * r;
        projected[1] += mode[n] * z;
        center[0] += r / 4.0;
        center[1] += z / 4.0;
        jrr += r * signs[n][0] / 4.0;
        jrz += r * signs[n][1] / 4.0;
        jzr += z * signs[n][0] / 4.0;
        jzz += z * signs[n][1] / 4.0;
    }
    const double det_center = jrr * jzz - jrz * jzr;
    if (!(det_center > 0.0) || !(center[0] > 0.0))
        throw std::domain_error("CAX4RT center geometry must be positive");
    for (std::size_t n = 0; n < 4; ++n) {
        result.gamma[n] = mode[n] - result.gradient[n][0] * projected[0] - result.gradient[n][1] * projected[1];
        for (std::size_t j = 4; j < 12; ++j) {
            const auto a = (j - 4) % 4, c = (j - 4) / 4;
            result.gamma_derivative[n][j] = -chain_scale * result.gradient[n][c] * mode[a];
            for (std::size_t d = 0; d < 2; ++d)
                result.gamma_derivative[n][j] -= result.gradient_derivative[n][d][j] * projected[d];
        }
    }
    // Abaqus/Standard CAX4RT thermal stabilization uses the first two local-node
    // mean gradients and the per-radian volume. Native operator identification
    // includes radial translation, aspect ratio, skew and a 1000-fold scale change.
    double norm = 0.0;
    for (std::size_t n = 0; n < 2; ++n)
        for (double g : result.gradient[n])
            norm += g * g;
    const double denominator = 12.0 * std::acos(-1.0);
    result.thermal_coefficient = result.volume * norm / denominator;
    for (std::size_t j = 4; j < 12; ++j) {
        double dnorm = 0.0;
        for (std::size_t n = 0; n < 2; ++n)
            for (std::size_t d = 0; d < 2; ++d)
                dnorm += 2.0 * result.gradient[n][d] * result.gradient_derivative[n][d][j];
        result.thermal_derivative[j] = (result.volume_derivative[j] * norm + result.volume * dnorm) / denominator;
    }
    return result;
}

Cax4LocalValues chain(const adlite::Scalar& x, const std::array<Cax4LocalValues, 6>& seeds) {
    std::array<double, 6> dx{};
    x.copy_derivatives(dx.data(), dx.size());
    Cax4LocalValues result{};
    for (std::size_t j = 0; j < 12; ++j)
        for (std::size_t k = 0; k < 6; ++k)
            result[j] += dx[k] * seeds[k][j];
    return result;
}

struct HourglassState final {
    double coefficient = 0.0;
    std::array<double, 2> amplitude{}, transported{};
    std::array<std::array<double, 2>, 2> deformation = {{{1, 0}, {0, 1}}};
};

HourglassState hourglass_state(const elements::Cax4Input& data,
    const Quad4RzGeometry& geometry,
    const ReducedGeometry& reference,
    const Cax4LocalValues& state) {
    HourglassState result;
    double norm = 0.0, normalization = 0.0, radius = 0.0, axial = 0.0;
    for (const auto& point : geometry.points) {
        radius += point.weighted_measure * point.radius / reference.volume;
        axial += point.weighted_measure * point.axial_coordinate / reference.volume;
    }
    for (std::size_t n = 0; n < 4; ++n) {
        normalization += reference.gamma[n] * mode[n];
        for (double b : reference.gradient[n])
            norm += b * b;
        for (std::size_t c = 0; c < 2; ++c) {
            result.amplitude[c] += reference.gamma[n] * state[4 + 4 * c + n];
            if (data.strain_formulation == StrainFormulation::finite)
                for (std::size_t d = 0; d < 2; ++d)
                    result.deformation[c][d] += reference.gradient[n][d] * state[4 + 4 * c + n];
        }
    }
    const auto initial = data.material.active_properties(data.initial_temperature, {0.0, radius, 0.0, axial});
    result.coefficient =
        .005 * initial.shear_modulus.value() * reference.volume * norm / (normalization * normalization);
    if (!std::isfinite(result.coefficient) || !(result.coefficient > 0.0))
        throw std::domain_error("CAX4RT requires a positive finite initial hourglass stiffness");
    for (std::size_t d = 0; d < 2; ++d)
        for (std::size_t c = 0; c < 2; ++c)
            result.transported[d] += result.deformation[c][d] * result.amplitude[c];
    return result;
}
} // namespace

elements::Cax4Result compute_cax4rt(const elements::Cax4Input& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    const Quad4MaterialHistory* history,
    double time_step,
    bool jacobian,
    bool thermal_time) {
    const bool finite = data.strain_formulation == StrainFormulation::finite;
    const auto reference = reduce_geometry(geometry, {}, 0.0);
    const auto current = finite ? reduce_geometry(geometry, state, jacobian ? 1.0 : 0.0) : reference;
    const auto old = finite ? reduce_geometry(geometry, committed, 0.0) : reference;
    Cax4LocalValues midpoint_state{};
    for (std::size_t j = 0; j < 12; ++j)
        midpoint_state[j] = (state[j] + committed[j]) / 2.0;
    const auto midpoint = finite ? reduce_geometry(geometry, midpoint_state, jacobian ? 0.5 : 0.0) : reference;
    std::array<double, 6> values{}, old_values{};
    std::array<Cax4LocalValues, 6> seeds{};
    for (std::size_t n = 0; n < 4; ++n) {
        for (std::size_t c = 0; c < 2; ++c)
            for (std::size_t d = 0; d < 2; ++d) {
                values[2 * c + d] += reference.gradient[n][d] * state[4 + 4 * c + n];
                old_values[2 * c + d] += reference.gradient[n][d] * committed[4 + 4 * c + n];
                seeds[2 * c + d][4 + 4 * c + n] = reference.gradient[n][d];
            }
        values[4] += reference.hoop[n] * state[4 + n];
        old_values[4] += reference.hoop[n] * committed[4 + n];
        seeds[4][4 + n] = reference.hoop[n];
        values[5] += current.measures[n] * state[n] / current.volume;
        old_values[5] += old.measures[n] * committed[n] / old.volume;
        seeds[5][n] = current.measures[n] / current.volume;
    }
    for (std::size_t j = 4; j < 12; ++j) {
        for (std::size_t n = 0; n < 4; ++n)
            seeds[5][j] += current.measure_derivative[n][j] * state[n] / current.volume;
        seeds[5][j] -= values[5] * current.volume_derivative[j] / current.volume;
    }
    double trace = 0.0;
    Cax4LocalValues trace_derivative{};
    for (std::size_t n = 0; n < 4; ++n) {
        const double ur = state[4 + n] - (finite ? committed[4 + n] : 0.0),
                     uz = state[8 + n] - (finite ? committed[8 + n] : 0.0);
        trace += ur * (midpoint.gradient[n][0] + midpoint.hoop[n]) + uz * midpoint.gradient[n][1];
        trace_derivative[4 + n] += midpoint.gradient[n][0] + midpoint.hoop[n];
        trace_derivative[8 + n] += midpoint.gradient[n][1];
        for (std::size_t j = 4; j < 12; ++j)
            trace_derivative[j] += ur * (midpoint.gradient_derivative[n][0][j] + midpoint.hoop_derivative[n][j])
                                   + uz * midpoint.gradient_derivative[n][1][j];
    }
    std::array<adlite::Scalar, 6> active;
    for (std::size_t i = 0; i < 6; ++i)
        active[i] = jacobian ? adlite::Scalar::independent(values[i], i, 6) : adlite::Scalar(values[i]);
    adlite::Scalar rr = active[0], zz = active[3], rz = (active[1] + active[2]) / 2.0, hoop = active[4];
    AxisymmetricRotation rotation = {1.0, 0.0, 0.0, 1.0, 1.0};
    const adlite::Scalar frr = 1.0 + active[0], frz = active[1], fzr = active[2], fzz = 1.0 + active[3],
                         fh = 1.0 + active[4];
    const adlite::Scalar fd = frr * fzz - frz * fzr;
    if (finite) {
        const double old_det = (1.0 + old_values[0]) * (1.0 + old_values[3]) - old_values[1] * old_values[2];
        if (!(fd.value() > 0.0) || !(fh.value() > 0.0) || !(old_det > 0.0) || !(1.0 + old_values[4] > 0.0))
            throw std::domain_error("CAX4RT averaged deformation must be positive");
        const adlite::Scalar a = 2.0 + active[0] + old_values[0], b = active[1] + old_values[1],
                             c = active[2] + old_values[2], d = 2.0 + active[3] + old_values[3], det = a * d - b * c;
        if (!(det.value() > 0.0))
            throw std::domain_error("CAX4RT averaged midpoint must be positive");
        const adlite::Scalar hrr = 2.0 * ((active[0] - old_values[0]) * d - (active[1] - old_values[1]) * c) / det,
                             hrz = 2.0 * (-(active[0] - old_values[0]) * b + (active[1] - old_values[1]) * a) / det,
                             hzr = 2.0 * ((active[2] - old_values[2]) * d - (active[3] - old_values[3]) * c) / det,
                             hzz = 2.0 * (-(active[2] - old_values[2]) * b + (active[3] - old_values[3]) * a) / det,
                             spin = (hrz - hzr) / 4.0, den = 1.0 + spin * spin, cs = (1.0 - spin * spin) / den,
                             sn = 2.0 * spin / den, shear = (hrz + hzr) / 2.0;
        rotation = {cs, sn, -sn, cs, 1.0};
        rr = cs * cs * hrr - 2.0 * cs * sn * shear + sn * sn * hzz;
        zz = sn * sn * hrr + 2.0 * cs * sn * shear + cs * cs * hzz;
        rz = cs * sn * (hrr - hzz) + (cs * cs - sn * sn) * shear;
        hoop = 2.0 * (active[4] - old_values[4]) / (2.0 + active[4] + old_values[4]);
    }
    const adlite::Scalar correction = (trace - hoop - rr - zz) / 2.0;
    rr += correction;
    zz += correction;
    const std::array<adlite::Scalar, 5> inputs = {rr, zz, hoop, rz, active[5]};
    std::array<double, 5> fed = {rr.value(), zz.value(), hoop.value(), rz.value(), active[5].value()};
    double radius = 0.0, axial = 0.0;
    for (const auto& p : geometry.points) {
        radius += p.weighted_measure * p.radius / reference.volume;
        axial += p.weighted_measure * p.axial_coordinate / reference.volume;
    }
    const MaterialFunctionContext context = {data.time, radius, 0.0, axial};
    if (finite && history) {
        auto old_context = context;
        old_context.time -= time_step;
        const auto eigen = data.material.eigenstrain_rz(old_values[5], old_context);
        const std::array<double, 4> imposed = {eigen.rr.value(),
            eigen.zz.value(),
            eigen.hoop.value(),
            eigen.rz.value()};
        for (std::size_t i = 0; i < 4; ++i)
            fed[i] += (*history)[0].elastic_strain[i] + (*history)[0].plastic_strain[i] + (*history)[0].creep_strain[i]
                      + imposed[i];
    }
    std::array<adlite::Scalar, 5> material_inputs;
    for (std::size_t i = 0; i < 5; ++i)
        material_inputs[i] = jacobian ? adlite::Scalar::independent(fed[i], i, 5) : adlite::Scalar(fed[i]);
    const auto raw = history ? data.material
                                   .response(material_inputs[0],
                                       material_inputs[1],
                                       material_inputs[2],
                                       material_inputs[3],
                                       material_inputs[4],
                                       time_step,
                                       (*history)[0],
                                       context)
                                   .stress
                             : data.material.stress(material_inputs[0],
                                   material_inputs[1],
                                   material_inputs[2],
                                   material_inputs[3],
                                   material_inputs[4],
                                   context);
    const std::array<adlite::Scalar, 4> raw_components = {raw.rr, raw.zz, raw.hoop, raw.rz};
    std::array<adlite::Scalar, 4> composed;
    std::array<double, 4> trace_response{};
    for (std::size_t i = 0; i < 4; ++i) {
        std::array<double, 5> partials{};
        raw_components[i].copy_derivatives(partials.data(), partials.size());
        composed[i] = jacobian ? adlite::compose(raw_components[i].value(), inputs.data(), partials.data(), 5)
                               : adlite::Scalar(raw_components[i].value());
        trace_response[i] = (partials[0] + partials[1]) / 2.0;
    }
    AxisymmetricStress sigma = {composed[0], composed[1], composed[2], composed[3]},
                       ds = {trace_response[0], trace_response[1], trace_response[2], trace_response[3]};
    if (finite) {
        sigma = rotate_axisymmetric_tensor(sigma, rotation);
        ds = rotate_axisymmetric_tensor(ds, rotation);
    }
    const std::array<adlite::Scalar, 4> stress = {sigma.rr, sigma.zz, sigma.hoop, sigma.rz};
    trace_response = {ds.rr.value(), ds.zz.value(), ds.hoop.value(), ds.rz.value()};
    std::array<Cax4LocalValues, 4> stress_derivative{};
    for (std::size_t i = 0; i < 4; ++i) {
        stress_derivative[i] = chain(stress[i], seeds);
        for (std::size_t j = 0; j < 12; ++j)
            stress_derivative[i][j] += trace_response[i] * trace_derivative[j];
    }
    elements::Cax4Result result;
    if (history) {
        const AxisymmetricRotation rot = {rotation.rr.value(),
            rotation.rz.value(),
            rotation.zr.value(),
            rotation.zz.value(),
            1.0};
        const auto response =
            finite ? data.material.incremental_response(rr.value(),
                         zz.value(),
                         hoop.value(),
                         rz.value(),
                         rot,
                         values[5],
                         old_values[5],
                         time_step,
                         (*history)[0],
                         context)
                   : data.material.response(fed[0], fed[1], fed[2], fed[3], fed[4], time_step, (*history)[0], context);
        result.history[0] = IsotropicThermoelasticMaterial::state_values(response.trial_state);
    }
    result.history[0].stress = {sigma.rr.value(), sigma.zz.value(), sigma.hoop.value(), sigma.rz.value()};
    const double p = (sigma.rr.value() + sigma.zz.value()) / 2.0;
    Cax4LocalValues dp{};
    for (std::size_t j = 0; j < 12; ++j)
        dp[j] = (stress_derivative[0][j] + stress_derivative[1][j]) / 2.0;
    for (std::size_t n = 0; n < 4; ++n) {
        const adlite::Scalar br = finite ? (reference.gradient[n][0] * fzz - reference.gradient[n][1] * fzr) / fd
                                         : adlite::Scalar(reference.gradient[n][0]);
        const adlite::Scalar bz = finite ? (-reference.gradient[n][0] * frz + reference.gradient[n][1] * frr) / fd
                                         : adlite::Scalar(reference.gradient[n][1]);
        const adlite::Scalar bh = finite ? reference.hoop[n] / fh : adlite::Scalar(reference.hoop[n]);
        const auto dr = chain(br, seeds), dz = chain(bz, seeds), dh = chain(bh, seeds);
        for (std::size_t c = 0; c < 2; ++c) {
            const auto row = 4 + 4 * c + n;
            const double dev = stress[c].value() - p;
            const double bmain = c == 0 ? br.value() : bz.value(), bcross = c == 0 ? bz.value() : br.value();
            const double value = dev * bmain + sigma.rz.value() * bcross + p * current.gradient[n][c]
                                 + (c == 0 ? p * current.hoop[n] + (sigma.hoop.value() - p) * bh.value() : 0.0);
            result.residual[row] = current.volume * value;
            if (!jacobian)
                continue;
            for (std::size_t j = 0; j < 12; ++j) {
                double dv = (stress_derivative[c][j] - dp[j]) * bmain + dev * (c == 0 ? dr[j] : dz[j])
                            + stress_derivative[3][j] * bcross + sigma.rz.value() * (c == 0 ? dz[j] : dr[j])
                            + dp[j] * current.gradient[n][c] + p * current.gradient_derivative[n][c][j];
                if (c == 0)
                    dv += dp[j] * current.hoop[n] + p * current.hoop_derivative[n][j]
                          + (stress_derivative[2][j] - dp[j]) * bh.value() + (sigma.hoop.value() - p) * dh[j];
                result.jacobian[row * 12 + j] = current.volume_derivative[j] * value + current.volume * dv;
            }
        }
    }
    // Total-stiffness hourglass energy: 0.5*C*|F^T*a|^2, including both derivatives of F.
    const auto hourglass = hourglass_state(data, geometry, reference, state);
    const double coefficient = hourglass.coefficient;
    const auto &amplitude = hourglass.amplitude, &transported = hourglass.transported;
    const auto& F = hourglass.deformation;
    for (std::size_t n = 0; n < 4; ++n)
        for (std::size_t c = 0; c < 2; ++c) {
            const auto row = 4 + 4 * c + n;
            double spatial = 0.0, gnb = 0.0;
            for (std::size_t d = 0; d < 2; ++d) {
                spatial += F[c][d] * transported[d];
                gnb += reference.gradient[n][d] * transported[d];
            }
            result.residual[row] += coefficient * (reference.gamma[n] * spatial + (finite ? amplitude[c] * gnb : 0.0));
            if (!jacobian)
                continue;
            for (std::size_t j = 0; j < 8; ++j) {
                const auto a = j % 4, q = j / 4;
                double spatial_derivative = 0.0, gn_derivative = 0.0;
                for (std::size_t d = 0; d < 2; ++d) {
                    const double db =
                        F[q][d] * reference.gamma[a] + (finite ? reference.gradient[a][d] * amplitude[q] : 0.0);
                    spatial_derivative +=
                        F[c][d] * db + (finite && c == q ? reference.gradient[a][d] * transported[d] : 0.0);
                    gn_derivative += reference.gradient[n][d] * db;
                }
                result.jacobian[row * 12 + 4 + j] +=
                    coefficient
                    * (reference.gamma[n] * spatial_derivative
                        + (finite ? (c == q ? reference.gamma[a] * gnb : 0.0) + amplitude[c] * gn_derivative : 0.0));
            }
        }
    const auto k =
        data.material.conductivity(jacobian ? adlite::Scalar::independent(values[5], 0, 1) : adlite::Scalar(values[5]),
            context);
    const double dk = k.is_active() ? k.derivative(0) : 0.0;
    std::array<double, 2> gradT{};
    double modalT = 0.0;
    for (std::size_t n = 0; n < 4; ++n) {
        modalT += current.gamma[n] * state[n];
        for (std::size_t d = 0; d < 2; ++d)
            gradT[d] += current.gradient[n][d] * state[n];
    }
    for (std::size_t n = 0; n < 4; ++n) {
        double uniform = 0.0;
        for (std::size_t d = 0; d < 2; ++d)
            uniform += current.gradient[n][d] * gradT[d];
        const double conduction = current.volume * uniform + current.thermal_coefficient * current.gamma[n] * modalT;
        result.residual[n] += k.value() * conduction - data.volumetric_heat_source * current.measures[n];
        if (jacobian)
            for (std::size_t j = 0; j < 12; ++j) {
                double du = 0.0, dm = j < 4 ? current.gamma[j] : 0.0;
                for (std::size_t m = 0; m < 4; ++m)
                    dm += current.gamma_derivative[m][j] * state[m];
                for (std::size_t d = 0; d < 2; ++d) {
                    double dg = j < 4 ? current.gradient[j][d] : 0.0;
                    for (std::size_t m = 0; m < 4; ++m)
                        dg += current.gradient_derivative[m][d][j] * state[m];
                    du += current.gradient_derivative[n][d][j] * gradT[d] + current.gradient[n][d] * dg;
                }
                const double dc =
                    current.volume_derivative[j] * uniform + current.volume * du
                    + current.thermal_derivative[j] * current.gamma[n] * modalT
                    + current.thermal_coefficient * (current.gamma_derivative[n][j] * modalT + current.gamma[n] * dm);
                result.jacobian[n * 12 + j] += dk * seeds[5][j] * conduction + k.value() * dc
                                               - data.volumetric_heat_source * current.measure_derivative[n][j];
            }
        if (history && thermal_time) {
            const auto cap = data.material.heat_capacity(jacobian ? adlite::Scalar::independent(state[n], 0, 1)
                                                                  : adlite::Scalar(state[n]),
                {data.time, geometry.coordinates[n].r, 0.0, geometry.coordinates[n].z});
            const double rate = (state[n] - committed[n]) / time_step;
            result.residual[n] += current.measures[n] * cap.value() * rate;
            if (jacobian) {
                for (std::size_t j = 4; j < 12; ++j)
                    result.jacobian[n * 12 + j] += current.measure_derivative[n][j] * cap.value() * rate;
                result.jacobian[n * 12 + n] +=
                    current.measures[n]
                    * (cap.value() / time_step + (cap.is_active() ? cap.derivative(0) : 0.0) * rate);
            }
        }
    }
    return result;
}

std::array<double, 2> cax4rt_thermal_rates(const elements::Cax4Input& data,
    const Quad4RzGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed,
    double time_step,
    bool thermal_time) {
    const auto current = reduce_geometry(geometry,
        data.strain_formulation == StrainFormulation::finite ? state : Cax4LocalValues{},
        0.0);
    std::array<double, 2> rates = {0.0, data.volumetric_heat_source * current.volume};
    if (thermal_time)
        for (std::size_t n = 0; n < 4; ++n)
            rates[0] +=
                current.measures[n]
                * data.material
                      .heat_capacity(state[n], {data.time, geometry.coordinates[n].r, 0.0, geometry.coordinates[n].z})
                      .value()
                * (state[n] - committed[n]) / time_step;
    return rates;
}

double cax4rt_hourglass_energy(const Cax4Input& data) {
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    const auto hourglass = hourglass_state(data, geometry, reduce_geometry(geometry, {}, 0.0), state);
    return 0.5 * hourglass.coefficient
           * (hourglass.transported[0] * hourglass.transported[0]
               + hourglass.transported[1] * hourglass.transported[1]);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Cax4Result evaluate_cax4rt(const Cax4Input& input, ElementRequest request) {
    const bool jacobian = request.jacobian;
    if (input.committed_history && (!std::isfinite(input.time_step) || !(input.time_step > 0)))
        throw std::invalid_argument("CAX4RT history update requires a positive finite time step");
    const auto& data = input;
    auto result = compute_cax4rt(data,
        input.geometry,
        input.state,
        input.committed_state,
        input.committed_history,
        input.time_step,
        jacobian,
        input.include_thermal_time_term);
    const auto rates = cax4rt_thermal_rates(data,
        input.geometry,
        input.state,
        input.committed_state,
        input.time_step,
        input.include_thermal_time_term);
    result.stored_heat_rate = rates[0];
    result.generated_heat_rate = rates[1];
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
