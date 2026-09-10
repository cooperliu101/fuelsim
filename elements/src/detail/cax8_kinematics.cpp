#include "cax8_kinematics.hpp"
#include <stdexcept>

namespace fuelsim::cax8_detail {

PointKinematics evaluate_kinematics(const Quad8RzPoint& p,
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

MechanicalGradient mechanical_gradient(const Quad8RzPoint& p, const PointKinematics& k, std::size_t n, bool finite) {
    const adlite::Scalar br = finite ? (p.gradient_r[n] * k.fzz - p.gradient_z[n] * k.fzr) / k.det
                                     : adlite::Scalar(p.gradient_r[n]),
                         bz = finite ? (-p.gradient_r[n] * k.frz + p.gradient_z[n] * k.frr) / k.det
                                     : adlite::Scalar(p.gradient_z[n]),
                         bh = finite ? p.shape[n] / k.radius : adlite::Scalar(p.shape[n] / p.radius);
    return {br, bz, bh};
}

ThermalGeometry thermal_geometry(const Quad8RzPoint& p, const PointKinematics& k, bool finite) {
    const adlite::Scalar a = 1 + (k.active[0] + k.old[0]) / 2, b = (k.active[1] + k.old[1]) / 2,
                         c = (k.active[2] + k.old[2]) / 2, d = 1 + (k.active[3] + k.old[3]) / 2, det = a * d - b * c;
    std::array<adlite::Scalar, 4> gr, gz;
    for (std::size_t n = 0; n < 4; ++n) {
        gr[n] = finite ? (p.temperature_gradient_r[n] * d - p.temperature_gradient_z[n] * c) / det
                       : adlite::Scalar(p.temperature_gradient_r[n]);
        gz[n] = finite ? (-p.temperature_gradient_r[n] * b + p.temperature_gradient_z[n] * a) / det
                       : adlite::Scalar(p.temperature_gradient_z[n]);
    }
    return {gr, gz};
}

SourceGeometry source_geometry(const Quad8RzPoint& p, const Quad8RzValues& state, bool finite) {
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
    return {source_measure, source_derivative};
}
} // namespace fuelsim::cax8_detail
