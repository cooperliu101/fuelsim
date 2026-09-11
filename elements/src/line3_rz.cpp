#include "line3_rz.hpp"
#include "contact_common.hpp"
#include "contact_types.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::elements {
Line3RzBoundaryResult compute_line3_rz_boundary(const Line3RzBoundaryData& data,
    const std::array<RzPoint, 3>& coordinates,
    const std::vector<double>& state,
    bool jacobian) {
    if (state.size() != 8)
        throw std::invalid_argument("CAX8T quadratic boundary requires eight local degrees of freedom");
    std::array<adlite::Scalar, 8> v;
    for (std::size_t i = 0; i < 8; ++i)
        v[i] = jacobian ? adlite::Scalar::independent(state[i], i, 8) : adlite::Scalar(state[i]);
    std::array<adlite::Scalar, 8> rows{};
    const bool current = data.use_displaced_geometry;
    const double load = data.load;
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> points = {-g, 0, g}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    for (std::size_t q = 0; q < 3; ++q) {
        const double x = points[q];
        const std::array<double, 3> shape = {x * (x - 1) / 2, x * (x + 1) / 2, 1 - x * x},
                                    derivative = {x - .5, x + .5, -2 * x};
        const std::array<double, 2> thermal = {(1 - x) / 2, (1 + x) / 2};
        adlite::Scalar radius = 0, tr = 0, tz = 0;
        for (std::size_t n = 0; n < 3; ++n) {
            const auto& point = coordinates[n];
            const adlite::Scalar r = current ? point.r + v[2 + n] : adlite::Scalar(point.r),
                                 z = current ? point.z + v[5 + n] : adlite::Scalar(point.z);
            radius += shape[n] * r;
            tr += derivative[n] * r;
            tz += derivative[n] * z;
        }
        const auto length = adlite::hypot(tr, tz), factor = 2 * std::acos(-1.0) * weights[q] * radius,
                   measure = factor * length;
        if (!(length.value() > 0) || !(radius.value() >= 0))
            throw std::domain_error("CAX8T boundary geometry is invalid");
        if (data.kind == Line3RzBoundaryKind::pressure)
            for (std::size_t n = 0; n < 3; ++n) {
                rows[2 + n] += shape[n] * load * factor * tz;
                rows[5 + n] -= shape[n] * load * factor * tr;
            }
        else if (data.kind == Line3RzBoundaryKind::traction) {
            const auto offset = data.component == TractionComponent::radial ? 2U : 5U;
            for (std::size_t n = 0; n < 3; ++n)
                rows[offset + n] -= shape[n] * load * measure;
        } else {
            adlite::Scalar flux = -load;
            if (data.kind == Line3RzBoundaryKind::convection) {
                flux = load * (thermal[0] * v[0] + thermal[1] * v[1] - data.ambient);
            }
            for (std::size_t n = 0; n < 2; ++n)
                rows[n] += thermal[n] * flux * measure;
        }
    }
    Line3RzBoundaryResult result;
    for (std::size_t i = 0; i < 8; ++i) {
        result.residual[i] = rows[i].value();
        if (jacobian)
            rows[i].copy_derivatives(result.jacobian.data() + 8 * i, 8);
    }
    return result;
}
} // namespace fuelsim::elements

namespace fuelsim::rz8 {
namespace {
struct ActivePoint final {
    adlite::Scalar r, z;
};

struct CurvePoint final {
    ActivePoint position, tangent, second;
    std::array<adlite::Scalar, 3> shape;
};

CurvePoint curve(const std::array<ActivePoint, 3>& points, const adlite::Scalar& x) {
    CurvePoint result;
    result.shape = {x * (x - 1) / 2, x * (x + 1) / 2, 1 - x * x};
    const std::array<adlite::Scalar, 3> derivative = {x - .5, x + .5, -2 * x};
    constexpr std::array<double, 3> second = {1, 1, -2};
    for (std::size_t n = 0; n < 3; ++n) {
        result.position.r += result.shape[n] * points[n].r;
        result.position.z += result.shape[n] * points[n].z;
        result.tangent.r += derivative[n] * points[n].r;
        result.tangent.z += derivative[n] * points[n].z;
        result.second.r += second[n] * points[n].r;
        result.second.z += second[n] * points[n].z;
    }
    return result;
}

struct Projection final {
    bool projected = false, clamped = false;
    double x = 0;
};

Projection closest(const std::array<ActivePoint, 3>& points,
    const ActivePoint& secondary,
    bool first,
    bool last,
    bool mechanical) {
    const double ar = (points[0].r.value() + points[1].r.value()) / 2 - points[2].r.value(),
                 az = (points[0].z.value() + points[1].z.value()) / 2 - points[2].z.value(),
                 br = (points[1].r.value() - points[0].r.value()) / 2,
                 bz = (points[1].z.value() - points[0].z.value()) / 2, cr = points[2].r.value() - secondary.r.value(),
                 cz = points[2].z.value() - secondary.z.value();
    const std::array<double, 4> f = {br * cr + bz * cz,
        br * br + bz * bz + 2 * (ar * cr + az * cz),
        3 * (ar * br + az * bz),
        2 * (ar * ar + az * az)};
    const auto evaluate = [&](double x) {
        return ((f[3] * x + f[2]) * x + f[1]) * x + f[0];
    };
    const double lower = first && mechanical ? -1.2 : -1.0, upper = last && mechanical ? 1.2 : 1.0;
    std::vector<double> cuts = {lower, upper};
    const double discriminant = 4 * f[2] * f[2] - 12 * f[3] * f[1];
    if (f[3] > 0 && discriminant >= 0)
        for (double sign : {-1., 1.}) {
            const double x = (-2 * f[2] + sign * std::sqrt(discriminant)) / (6 * f[3]);
            if (x > lower && x < upper)
                cuts.push_back(x);
        }
    std::sort(cuts.begin(), cuts.end());
    std::vector<double> roots;
    const double tolerance = 1e-12 * (std::abs(f[0]) + std::abs(f[1]) + std::abs(f[2]) + std::abs(f[3]));
    for (double x : cuts)
        if (std::abs(evaluate(x)) <= tolerance)
            roots.push_back(x);
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        double left = cuts[i - 1], right = cuts[i], fl = evaluate(left), fr = evaluate(right);
        if (fl * fr >= 0)
            continue;
        for (int iteration = 0; iteration < 60; ++iteration) {
            const double middle = (left + right) / 2, fm = evaluate(middle);
            if ((fm > 0) == (fl > 0)) {
                left = middle;
                fl = fm;
            } else
                right = middle;
        }
        roots.push_back((left + right) / 2);
    }
    Projection result;
    double distance = std::numeric_limits<double>::infinity();
    const auto consider = [&](double x, bool clamped) {
        if (x >= 1 - 1e-12 && !last)
            return;
        const double r = (ar * x + br) * x + cr, z = (az * x + bz) * x + cz, d = r * r + z * z;
        if (d < distance) {
            distance = d;
            result = {true, clamped, x};
        }
    };
    for (double x : roots)
        if ((3 * f[3] * x + 2 * f[2]) * x + f[1] > 0)
            consider(x, false);
    return result;
}
} // namespace

std::pair<bool, double> project_line3(const std::array<RzPoint, 3>& primary, RzPoint point) {
    std::array<ActivePoint, 3> coordinates;
    for (std::size_t n = 0; n < 3; ++n)
        coordinates[n] = {primary[n].r, primary[n].z};
    const auto result = closest(coordinates, {point.r, point.z}, true, true, false);
    return {result.projected, result.x};
}

Line3ContactResult compute_line3_contact(const Line3ContactGeometry& geometry,
    const GapHeatProperties& heat,
    const NormalContactProperties& mechanical,
    const std::vector<double>& state,
    const std::vector<double>& old,
    const ContactPointHistory& history,
    bool jacobian) {
    if (state.size() != 16 || old.size() != 16)
        throw std::invalid_argument("Quadratic RZ contact requires sixteen local degrees of freedom");
    std::array<adlite::Scalar, 16> v;
    for (std::size_t i = 0; i < 16; ++i)
        v[i] = jacobian ? adlite::Scalar::independent(state[i], i, 16) : adlite::Scalar(state[i]);
    std::array<ActivePoint, 3> secondary, primary;
    for (std::size_t n = 0; n < 3; ++n) {
        secondary[n] = {geometry.secondary[n].r + v[4 + n], geometry.secondary[n].z + v[7 + n]};
        primary[n] = {geometry.primary[n].r + v[10 + n], geometry.primary[n].z + v[13 + n]};
    }
    const double xs = geometry.mechanical ? (geometry.secondary_node == 0      ? -1
                                                : geometry.secondary_node == 1 ? 1
                                                                               : 0)
                                          : geometry.coordinate;
    const auto sp = curve(secondary, xs);
    const auto projection = closest(primary,
        sp.position,
        geometry.primary_first,
        geometry.primary_last,
        geometry.mechanical || geometry.nodal_heat);
    Line3ContactResult result;
    if (!projection.projected)
        return result;
    const auto trial = curve(primary, projection.x);
    const auto dr = trial.position.r - sp.position.r, dz = trial.position.z - sp.position.z;
    const auto equation = dr * trial.tangent.r + dz * trial.tangent.z;
    const double derivative = (trial.tangent.r * trial.tangent.r + trial.tangent.z * trial.tangent.z
                               + dr * trial.second.r + dz * trial.second.z)
                                  .value();
    if (!(derivative > 0) || !std::isfinite(derivative))
        throw std::domain_error("Quadratic RZ projection is singular");
    const adlite::Scalar x =
        projection.clamped ? adlite::Scalar(projection.x) : projection.x - (equation - equation.value()) / derivative;
    const auto pp = curve(primary, x);
    const auto length = adlite::hypot(pp.tangent.r, pp.tangent.z);
    if (!(length.value() > 0) || !(sp.position.r.value() >= 0))
        throw std::domain_error("Quadratic RZ contact geometry is invalid");
    const adlite::Scalar nr = -pp.tangent.z / length, nz = pp.tangent.r / length, tr = pp.tangent.r / length,
                         tz = pp.tangent.z / length,
                         gap = (pp.position.r - sp.position.r) * nr + (pp.position.z - sp.position.z) * nz;
    result.primary_coordinate = x.value();
    std::array<adlite::Scalar, 16> rows{};
    if (!geometry.mechanical) {
        const adlite::Scalar ts = (1 - xs) / 2 * v[0] + (1 + xs) / 2 * v[1],
                             tp = (1 - x) / 2 * v[2] + (1 + x) / 2 * v[3];
        const adlite::Scalar h = contact_common::gap_conductance(heat, gap, ts, tp);
        // Native NTS heat transfer uses only temperature corner nodes. Its area
        // integrates linear shape functions on the current corner chord, even
        // when mechanical projection uses the full quadratic primary geometry.
        const auto measure =
            geometry.nodal_heat
                ? std::acos(-1.0) * adlite::hypot(secondary[1].r - secondary[0].r, secondary[1].z - secondary[0].z)
                      * ((2 * secondary[geometry.secondary_node].r + secondary[1 - geometry.secondary_node].r) / 3)
                : 2 * std::acos(-1.0) * sp.position.r * adlite::hypot(sp.tangent.r, sp.tangent.z) * geometry.weight;
        const auto flux = h * (ts - tp), rate = flux * measure;
        rows[0] += (1 - xs) / 2 * rate;
        rows[1] += (1 + xs) / 2 * rate;
        rows[2] -= (1 - x) / 2 * rate;
        rows[3] -= (1 + x) / 2 * rate;
        result.thermal = {true, gap.value(), flux.value(), measure.value()};
    } else {
        const adlite::Scalar trial_pressure =
            (mechanical.augmented_lagrangian ? history.normal_multiplier : 0) - mechanical.penalty * gap;
        const adlite::Scalar pressure = trial_pressure.value() > 0 ? trial_pressure : adlite::Scalar(0);
        adlite::Scalar area = 0, tributary_length = 0;
        const double g = std::sqrt(3.0 / 5.0);
        const std::array<double, 3> qs = {-g, 0, g}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9};
        for (std::size_t q = 0; q < 3; ++q) {
            const auto p = curve(secondary, qs[q]);
            const auto line = p.shape[geometry.secondary_node] * adlite::hypot(p.tangent.r, p.tangent.z) * weights[q];
            area += 2 * std::acos(-1.0) * p.position.r * line;
            tributary_length += line;
        }
        if (!(area.value() > 0))
            throw std::domain_error("Quadratic RZ contact node requires positive tributary area");
        adlite::Scalar traction = 0, elastic = 0, total = history.total_tangential_slip;
        bool sliding = false;
        if (mechanical.friction_coefficient > 0 && pressure.value() > 0) {
            adlite::Scalar relative_r = v[4 + geometry.secondary_node] - old[4 + geometry.secondary_node],
                           relative_z = v[7 + geometry.secondary_node] - old[7 + geometry.secondary_node];
            for (std::size_t n = 0; n < 3; ++n) {
                relative_r -= pp.shape[n] * (v[10 + n] - old[10 + n]);
                relative_z -= pp.shape[n] * (v[13 + n] - old[13 + n]);
            }
            const auto increment = relative_r * tr + relative_z * tz;
            elastic = history.elastic_tangential_slip + increment;
            total += increment;
            double penalty = mechanical.penalty;
            adlite::Scalar trial_traction = penalty * elastic, limit = mechanical.friction_coefficient * pressure;
            if (mechanical.maximum_elastic_slip > 0) {
                const auto stiffness = limit / mechanical.maximum_elastic_slip;
                trial_traction = stiffness * elastic;
                const double magnitude = std::abs(trial_traction.value());
                if (magnitude > limit.value() || (magnitude == limit.value() && history.sliding)) {
                    sliding = true;
                    traction = trial_traction.value() > 0 ? limit : -limit;
                    elastic = traction / stiffness;
                } else
                    traction = trial_traction;
            } else if (std::abs(trial_traction.value()) > limit.value()) {
                sliding = true;
                traction = trial_traction.value() > 0 ? limit : -limit;
                elastic = traction / penalty;
            } else
                traction = trial_traction;
        }
        const auto fr = area * (pressure * nr + traction * tr), fz = area * (pressure * nz + traction * tz);
        rows[4 + geometry.secondary_node] += fr;
        rows[7 + geometry.secondary_node] += fz;
        for (std::size_t n = 0; n < 3; ++n) {
            rows[10 + n] -= pp.shape[n] * fr;
            rows[13 + n] -= pp.shape[n] * fz;
        }
        result.mechanical = {true,
            gap.value(),
            pressure.value(),
            area.value(),
            tributary_length.value(),
            (pressure * area).value(),
            traction.value(),
            (traction * area).value(),
            elastic.value(),
            sliding,
            total.value()};
    }
    for (std::size_t i = 0; i < 16; ++i) {
        result.residual[i] = rows[i].value();
        if (jacobian)
            rows[i].copy_derivatives(result.jacobian.data() + 16 * i, 16);
    }
    return result;
}
} // namespace fuelsim::rz8
