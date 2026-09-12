#include "cax8t.hpp"
#include "cax_common.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
Quad8RzPoint evaluate_quad8_rz_point(const Quad8RzCoordinates& coordinates, double x, double y, double weight) {
    constexpr std::array<std::array<double, 2>, 4> signs = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
    Quad8RzPoint p;
    std::array<double, 8> dx{}, dy{};
    std::array<double, 4> tx{}, ty{};
    for (std::size_t n = 0; n < 4; ++n) {
        const double a = signs[n][0], b = signs[n][1];
        p.shape[n] = .25 * (1 + a * x) * (1 + b * y) * (a * x + b * y - 1);
        dx[n] = .25 * a * (1 + b * y) * (2 * a * x + b * y);
        dy[n] = .25 * b * (1 + a * x) * (a * x + 2 * b * y);
        p.temperature_shape[n] = .25 * (1 + a * x) * (1 + b * y);
        tx[n] = .25 * a * (1 + b * y);
        ty[n] = .25 * b * (1 + a * x);
    }
    p.shape[4] = .5 * (1 - x * x) * (1 - y);
    dx[4] = -x * (1 - y);
    dy[4] = -.5 * (1 - x * x);
    p.shape[5] = .5 * (1 + x) * (1 - y * y);
    dx[5] = .5 * (1 - y * y);
    dy[5] = -(1 + x) * y;
    p.shape[6] = .5 * (1 - x * x) * (1 + y);
    dx[6] = -x * (1 + y);
    dy[6] = .5 * (1 - x * x);
    p.shape[7] = .5 * (1 - x) * (1 - y * y);
    dx[7] = -.5 * (1 - y * y);
    dy[7] = -(1 - x) * y;
    double a = 0, b = 0, c = 0, d = 0;
    for (std::size_t n = 0; n < 8; ++n) {
        a += dx[n] * coordinates[n].r;
        b += dy[n] * coordinates[n].r;
        c += dx[n] * coordinates[n].z;
        d += dy[n] * coordinates[n].z;
        p.radius += p.shape[n] * coordinates[n].r;
        p.axial_coordinate += p.shape[n] * coordinates[n].z;
    }
    const double det = a * d - b * c;
    if (!std::isfinite(det) || !(det > 0) || !std::isfinite(p.radius) || !(p.radius > 0))
        throw std::domain_error("CAX8 requires positive reference Jacobian and integration-point radius");
    p.weighted_measure = 2 * std::acos(-1.0) * p.radius * det * weight;
    for (std::size_t n = 0; n < 8; ++n) {
        p.gradient_r[n] = (dx[n] * d - dy[n] * c) / det;
        p.gradient_z[n] = (dy[n] * a - dx[n] * b) / det;
    }
    for (std::size_t n = 0; n < 4; ++n) {
        p.temperature_gradient_r[n] = (tx[n] * d - ty[n] * c) / det;
        p.temperature_gradient_z[n] = (ty[n] * a - tx[n] * b) / det;
    }
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    for (std::size_t n = 0; n < 4; ++n) {
        a += tx[n] * coordinates[n].r;
        b += ty[n] * coordinates[n].r;
        c += tx[n] * coordinates[n].z;
        d += ty[n] * coordinates[n].z;
        p.source_radius += p.temperature_shape[n] * coordinates[n].r;
    }
    const double source_det = a * d - b * c;
    if (!(source_det > 0) || !(p.source_radius > 0))
        throw std::domain_error("CAX8 linear source geometry must be positive");
    p.source_measure = 2 * std::acos(-1.0) * p.source_radius * source_det * weight;
    for (std::size_t n = 0; n < 4; ++n) {
        p.source_gradient_r[n] = (tx[n] * d - ty[n] * c) / source_det;
        p.source_gradient_z[n] = (ty[n] * a - tx[n] * b) / source_det;
    }
    return p;
}

struct PointKinematics final {
    std::array<adlite::Scalar, 6> active;
    std::array<Quad8RzValues, 6> chain{};
    std::array<double, 6> old{};
    std::array<adlite::Scalar, 4> strain;
    AxisymmetricRotation rotation;
    adlite::Scalar frr, frz, fzr, fzz, det, radius, measure;
};

struct MechanicalGradient final {
    adlite::Scalar radial, axial, hoop;
};

struct ThermalGeometry final {
    std::array<adlite::Scalar, 4> radial, axial;
};

struct SourceGeometry final {
    double measure = 0.0;
    Quad8RzValues derivative{};
};

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
        throw std::domain_error("CAX8 committed and current deformation and radius must remain positive");
    const auto midpoint = evaluate_axisymmetric_midpoint_increment(
        {2 + v[0] + k.old[0], v[1] + k.old[1], v[2] + k.old[2], 2 + v[3] + k.old[3]},
        {v[0] - k.old[0], v[1] - k.old[1], v[2] - k.old[2], v[3] - k.old[3]},
        2 * p.radius + v[4] + k.old[4],
        v[4] - k.old[4]);
    k.rotation = midpoint.in_plane.rotation;
    k.strain = {midpoint.in_plane.rr, midpoint.in_plane.zz, midpoint.hoop, midpoint.in_plane.rz};
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
            throw std::domain_error("CAX8 current linear source geometry must be positive");
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

Cax8Result evaluate_cax8t(const Cax8Input& data, ElementRequest request, Cax8Quadrature quadrature) {
    if (quadrature != Cax8Quadrature::full && quadrature != Cax8Quadrature::reduced)
        throw std::invalid_argument("Invalid CAX8 quadrature");
    const std::size_t expected = quadrature == Cax8Quadrature::full ? 9 : 4;
    if (data.geometry.point_count != expected)
        throw std::invalid_argument("CAX8 geometry does not match the selected quadrature");
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    const auto& committed = data.committed_state;
    const auto* history = data.committed_history;
    const double dt = data.time_step;
    const bool jacobian = request.jacobian;
    const bool thermal_time = data.include_thermal_time_term;
    validate_cax_time_input(history != nullptr, dt, thermal_time);
    const bool finite = data.strain_formulation == StrainFormulation::finite;
    elements::Cax8Result result;
    for (std::size_t q = 0; q < geometry.point_count; ++q) {
        const auto& p = geometry.points[q];
        const auto k = evaluate_kinematics(p, state, committed, finite, jacobian);
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
        const auto tangent = evaluate_axisymmetric_material_response(data.material,
            {fed[0], fed[1], fed[2], fed[3]},
            fed[4],
            dt,
            history ? &(*history)[q] : nullptr,
            context,
            jacobian,
            request.history);
        const std::array<adlite::Scalar, 5> inputs = {k.strain[0], k.strain[1], k.strain[2], k.strain[3], t};
        AxisymmetricStress stress = compose_axisymmetric_stress(tangent, inputs, jacobian);
        if (finite)
            stress = rotate_axisymmetric_tensor(stress, k.rotation);
        if (history && request.history) {
            result.history[q] = tangent.history;
            if (finite)
                rotate_axisymmetric_strain_history(result.history[q], k.rotation);
        }
        result.history[q].stress = {stress.rr.value(), stress.zz.value(), stress.hoop.value(), stress.rz.value()};
        for (std::size_t n = 0; n < 8; ++n) {
            const auto gradient = mechanical_gradient(p, k, n, finite);
            const auto& br = gradient.radial;
            const auto& bz = gradient.axial;
            const auto& bh = gradient.hoop;
            add_row(result, 4 + n, k.measure * (br * stress.rr + bz * stress.rz + bh * stress.hoop), k, jacobian);
            add_row(result, 12 + n, k.measure * (bz * stress.zz + br * stress.rz), k, jacobian);
        }
        // Four temperature shape functions share the full quadratic geometric map.
        const auto thermal_map = thermal_geometry(p, k, finite);
        const auto& gr = thermal_map.radial;
        const auto& gz = thermal_map.axial;
        adlite::Scalar tr = 0, tz = 0;
        for (std::size_t n = 0; n < 4; ++n) {
            tr += gr[n] * state[n];
            tz += gz[n] * state[n];
        }
        const auto conductivity = data.material.conductivity(t, context);
        const adlite::Scalar capacity =
            history && thermal_time ? data.material.heat_capacity(t, context) * (t - k.old[5]) / dt : adlite::Scalar(0);
        const auto source_map = source_geometry(p, state, finite);
        const double source_measure = source_map.measure;
        const auto& source_derivative = source_map.derivative;
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
    if (request.stress)
        for (std::size_t q = 0; q < result.history.size(); ++q)
            result.stress[q] = result.history[q].stress;
    if (!request.residual && !request.jacobian)
        result.residual.fill(0.0);
    if (!request.history)
        result.history = {};
    return result;
}

Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates, Cax8Quadrature quadrature) {
    if (quadrature != Cax8Quadrature::full && quadrature != Cax8Quadrature::reduced)
        throw std::invalid_argument("Invalid CAX8 quadrature");
    Quad8RzGeometry result{coordinates, {}};
    if (quadrature == Cax8Quadrature::reduced) {
        result.point_count = 4;
        const double g = 1 / std::sqrt(3.0);
        const std::array<double, 4> xi = {-g, g, g, -g}, eta = {-g, -g, g, g};
        for (std::size_t q = 0; q < 4; ++q)
            result.points[q] = evaluate_quad8_rz_point(coordinates, xi[q], eta[q], 1);
        return result;
    }
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> q = {-g, 0, g}, w = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            result.points[3 * j + i] = evaluate_quad8_rz_point(coordinates, q[i], q[j], w[i] * w[j]);
    return result;
}

Cax8Result evaluate_cax8t(const Cax8Input& input, ElementRequest request) {
    return evaluate_cax8t(input, request, Cax8Quadrature::full);
}

Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates) {
    return make_cax8t_geometry(coordinates, Cax8Quadrature::full);
}
} // namespace fuelsim::elements
