#include "cax8_assembly.hpp"
#include "cax8_types.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
struct PointKinematics final {
    std::array<adlite::Scalar, 6> active;
    std::array<Quad8RzValues, 6> chain{};
    std::array<double, 6> old{};
    std::array<adlite::Scalar, 4> strain;
    AxisymmetricRotation rotation;
    adlite::Scalar frr, frz, fzr, fzz, det, radius, measure;
};

PointKinematics kinematics(const Quad8RzPoint& p,
    const Quad8RzValues& state,
    const Quad8RzValues& committed,
    bool finite,
    bool jacobian) {
    PointKinematics k;
    std::array<double, 6> value{};
    for (std::size_t n = 0; n < 8; ++n) {
        for (std::size_t c = 0; c < 2; ++c) {
            const auto index = 4 + 8 * c + n;
            value[2 * c] += state[index] * p.gradient_r[n];
            value[2 * c + 1] += state[index] * p.gradient_z[n];
            k.old[2 * c] += committed[index] * p.gradient_r[n];
            k.old[2 * c + 1] += committed[index] * p.gradient_z[n];
            k.chain[2 * c][index] = p.gradient_r[n];
            k.chain[2 * c + 1][index] = p.gradient_z[n];
        }
        value[4] += state[4 + n] * p.shape[n];
        k.old[4] += committed[4 + n] * p.shape[n];
        k.chain[4][4 + n] = p.shape[n];
    }
    for (std::size_t n = 0; n < 4; ++n) {
        value[5] += state[n] * p.temperature_shape[n];
        k.old[5] += committed[n] * p.temperature_shape[n];
        k.chain[5][n] = p.temperature_shape[n];
    }
    for (std::size_t i = 0; i < 6; ++i)
        k.active[i] = jacobian ? adlite::Scalar::independent(value[i], i, 6) : adlite::Scalar(value[i]);
    const auto& v = k.active;
    k.frr = 1 + v[0];
    k.frz = v[1];
    k.fzr = v[2];
    k.fzz = 1 + v[3];
    k.det = k.frr * k.fzz - k.frz * k.fzr;
    k.radius = p.radius + v[4];
    k.measure = finite ? p.weighted_measure * k.det * k.radius / p.radius : adlite::Scalar(p.weighted_measure);
    k.strain = {v[0], v[3], v[4] / p.radius, (v[1] + v[2]) / 2};
    if (!finite)
        return k;
    const double old_det = (1 + k.old[0]) * (1 + k.old[3]) - k.old[1] * k.old[2];
    if (!(k.det.value() > 0) || !(old_det > 0) || !(k.radius.value() > 0) || !(p.radius + k.old[4] > 0))
        throw std::domain_error("CAX8T committed and current deformation and radius must remain positive");
    const adlite::Scalar a = 2 + v[0] + k.old[0], b = v[1] + k.old[1], c = v[2] + k.old[2], d = 2 + v[3] + k.old[3],
                         det = a * d - b * c;
    if (!(det.value() > 0))
        throw std::domain_error("CAX8T midpoint deformation must remain positive");
    const adlite::Scalar hrr = 2 * ((v[0] - k.old[0]) * d - (v[1] - k.old[1]) * c) / det,
                         hrz = 2 * (-(v[0] - k.old[0]) * b + (v[1] - k.old[1]) * a) / det,
                         hzr = 2 * ((v[2] - k.old[2]) * d - (v[3] - k.old[3]) * c) / det,
                         hzz = 2 * (-(v[2] - k.old[2]) * b + (v[3] - k.old[3]) * a) / det, spin = (hrz - hzr) / 4,
                         den = 1 + spin * spin, cs = (1 - spin * spin) / den, sn = 2 * spin / den,
                         shear = (hrz + hzr) / 2;
    k.rotation = {cs, sn, -sn, cs, 1.0};
    k.strain = {cs * cs * hrr - 2 * cs * sn * shear + sn * sn * hzz,
        sn * sn * hrr + 2 * cs * sn * shear + cs * cs * hzz,
        2 * (v[4] - k.old[4]) / (2 * p.radius + v[4] + k.old[4]),
        cs * sn * (hrr - hzz) + (cs * cs - sn * sn) * shear};
    return k;
}

void add_row(elements::Cax8Result& result,
    std::size_t row,
    const adlite::Scalar& value,
    const PointKinematics& k,
    bool jacobian) {
    result.residual[row] += value.value();
    if (!jacobian)
        return;
    std::array<double, 6> d{};
    value.copy_derivatives(d.data(), 6);
    for (std::size_t j = 0; j < 20; ++j)
        for (std::size_t a = 0; a < 6; ++a)
            result.jacobian[20 * row + j] += d[a] * k.chain[a][j];
}
} // namespace

elements::Cax8Result compute_quad8_rz(const elements::Cax8Input& data,
    const Quad8RzGeometry& geometry,
    const Quad8RzValues& state,
    const Quad8RzValues& committed,
    const Quad8MaterialHistory* history,
    double dt,
    bool jacobian,
    bool thermal_time) {
    if (history && (!(dt > 0) || !std::isfinite(dt)))
        throw std::invalid_argument("CAX8T material update requires positive finite time step");
    const bool finite = data.strain_formulation == StrainFormulation::finite;
    elements::Cax8Result result;
    for (std::size_t q = 0; q < geometry.point_count; ++q) {
        const auto& p = geometry.points[q];
        const auto k = kinematics(p, state, committed, finite, jacobian);
        const auto& t = k.active[5];
        const MaterialFunctionContext context = {data.time, p.radius, 0, p.axial_coordinate};
        std::array<double, 5> fed = {k.strain[0].value(),
            k.strain[1].value(),
            k.strain[2].value(),
            k.strain[3].value(),
            t.value()};
        if (finite && history) {
            auto old_context = context;
            old_context.time -= dt;
            const auto eigen = data.material.eigenstrain_rz(k.old[5], old_context);
            const std::array<double, 4> imposed = {eigen.rr.value(),
                eigen.zz.value(),
                eigen.hoop.value(),
                eigen.rz.value()};
            for (std::size_t c = 0; c < 4; ++c)
                fed[c] += (*history)[q].elastic_strain[c] + (*history)[q].plastic_strain[c]
                          + (*history)[q].creep_strain[c] + imposed[c];
        }
        std::array<adlite::Scalar, 5> material;
        for (std::size_t i = 0; i < 5; ++i)
            material[i] = jacobian ? adlite::Scalar::independent(fed[i], i, 5) : adlite::Scalar(fed[i]);
        const auto raw =
            history ? data.material
                          .response(material[0],
                              material[1],
                              material[2],
                              material[3],
                              material[4],
                              dt,
                              (*history)[q],
                              context)
                          .stress
                    : data.material.stress(material[0], material[1], material[2], material[3], material[4], context);
        const std::array<adlite::Scalar, 4> components = {raw.rr, raw.zz, raw.hoop, raw.rz};
        const std::array<adlite::Scalar, 5> inputs = {k.strain[0], k.strain[1], k.strain[2], k.strain[3], t};
        std::array<adlite::Scalar, 4> composed;
        for (std::size_t c = 0; c < 4; ++c) {
            std::array<double, 5> partials{};
            if (jacobian)
                components[c].copy_derivatives(partials.data(), 5);
            composed[c] = jacobian ? adlite::compose(components[c].value(), inputs.data(), partials.data(), 5)
                                   : adlite::Scalar(components[c].value());
        }
        AxisymmetricStress stress = {composed[0], composed[1], composed[2], composed[3]};
        if (finite)
            stress = rotate_axisymmetric_tensor(stress, k.rotation);
        if (history) {
            const AxisymmetricRotation rotation = {k.rotation.rr.value(),
                k.rotation.rz.value(),
                k.rotation.zr.value(),
                k.rotation.zz.value(),
                1.0};
            const auto response =
                finite ? data.material.incremental_response(k.strain[0].value(),
                             k.strain[1].value(),
                             k.strain[2].value(),
                             k.strain[3].value(),
                             rotation,
                             t.value(),
                             k.old[5],
                             dt,
                             (*history)[q],
                             context)
                       : data.material.response(fed[0], fed[1], fed[2], fed[3], fed[4], dt, (*history)[q], context);
            result.history[q] = IsotropicThermoelasticMaterial::state_values(response.trial_state);
        }
        result.history[q].stress = {stress.rr.value(), stress.zz.value(), stress.hoop.value(), stress.rz.value()};
        for (std::size_t n = 0; n < 8; ++n) {
            const adlite::Scalar br = finite ? (p.gradient_r[n] * k.fzz - p.gradient_z[n] * k.fzr) / k.det
                                             : adlite::Scalar(p.gradient_r[n]),
                                 bz = finite ? (-p.gradient_r[n] * k.frz + p.gradient_z[n] * k.frr) / k.det
                                             : adlite::Scalar(p.gradient_z[n]),
                                 bh = finite ? p.shape[n] / k.radius : adlite::Scalar(p.shape[n] / p.radius);
            add_row(result, 4 + n, k.measure * (br * stress.rr + bz * stress.rz + bh * stress.hoop), k, jacobian);
            add_row(result, 12 + n, k.measure * (bz * stress.zz + br * stress.rz), k, jacobian);
        }
        // Four temperature shape functions share the full quadratic geometric map.
        const adlite::Scalar a = 1 + (k.active[0] + k.old[0]) / 2, b = (k.active[1] + k.old[1]) / 2,
                             c = (k.active[2] + k.old[2]) / 2, d = 1 + (k.active[3] + k.old[3]) / 2,
                             det = a * d - b * c;
        std::array<adlite::Scalar, 4> gr, gz;
        adlite::Scalar tr = 0, tz = 0;
        for (std::size_t n = 0; n < 4; ++n) {
            gr[n] = finite ? (p.temperature_gradient_r[n] * d - p.temperature_gradient_z[n] * c) / det
                           : adlite::Scalar(p.temperature_gradient_r[n]);
            gz[n] = finite ? (-p.temperature_gradient_r[n] * b + p.temperature_gradient_z[n] * a) / det
                           : adlite::Scalar(p.temperature_gradient_z[n]);
            tr += gr[n] * state[n];
            tz += gz[n] * state[n];
        }
        const auto conductivity = data.material.conductivity(t, context);
        const adlite::Scalar capacity =
            history && thermal_time ? data.material.heat_capacity(t, context) * (t - k.old[5]) / dt : adlite::Scalar(0);
        double source_measure = p.weighted_measure;
        Quad8RzValues source_derivative{};
        if (finite) {
            double fa = 1, fb = 0, fc = 0, fd = 1, radius = p.source_radius;
            for (std::size_t n = 0; n < 4; ++n) {
                fa += state[4 + n] * p.source_gradient_r[n];
                fb += state[4 + n] * p.source_gradient_z[n];
                fc += state[12 + n] * p.source_gradient_r[n];
                fd += state[12 + n] * p.source_gradient_z[n];
                radius += state[4 + n] * p.temperature_shape[n];
            }
            const double determinant = fa * fd - fb * fc;
            if (!(determinant > 0) || !(radius > 0))
                throw std::domain_error("CAX8T current linear source geometry must be positive");
            source_measure = p.source_measure * determinant * radius / p.source_radius;
            for (std::size_t n = 0; n < 4; ++n) {
                source_derivative[4 + n] = source_measure
                                           * ((fd * p.source_gradient_r[n] - fc * p.source_gradient_z[n]) / determinant
                                               + p.temperature_shape[n] / radius);
                source_derivative[12 + n] =
                    source_measure * ((fa * p.source_gradient_z[n] - fb * p.source_gradient_r[n]) / determinant);
            }
        }
        for (std::size_t n = 0; n < 4; ++n) {
            const auto row = k.measure * (conductivity * (gr[n] * tr + gz[n] * tz) + p.temperature_shape[n] * capacity);
            add_row(result, n, row, k, jacobian);
            const double source = p.temperature_shape[n] * data.volumetric_heat_source;
            result.residual[n] -= source * source_measure;
            if (jacobian)
                for (std::size_t j = 4; j < 20; ++j)
                    result.jacobian[20 * n + j] -= source * source_derivative[j];
            if (jacobian)
                for (std::size_t j = 0; j < 4; ++j)
                    result.jacobian[20 * n + j] += k.measure.value() * conductivity.value()
                                                   * (gr[n].value() * gr[j].value() + gz[n].value() * gz[j].value());
        }
        result.generated_heat_rate += source_measure * data.volumetric_heat_source;
        result.stored_heat_rate += k.measure.value() * capacity.value();
    }
    return result;
}
} // namespace fuelsim
