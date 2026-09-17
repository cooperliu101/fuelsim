#include "line3_plane.hpp"
#include "contact_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace fuelsim::elements {
namespace {
using Point = std::array<adlite::Scalar, 2>;

struct Curve final {
    Point position{}, tangent{}, second{};
    std::array<adlite::Scalar, 3> shape;
};

Curve curve(const std::array<Point, 3>& nodes, const adlite::Scalar& xi) {
    Curve result;
    result.shape = {0.5 * xi * (xi - 1), 0.5 * xi * (xi + 1), 1 - xi * xi};
    const std::array<adlite::Scalar, 3> derivative{xi - 0.5, xi + 0.5, -2 * xi};
    constexpr std::array<double, 3> second{1, 1, -2};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t c = 0; c < 2; ++c) {
            result.position[c] += result.shape[i] * nodes[i][c];
            result.tangent[c] += derivative[i] * nodes[i][c];
            result.second[c] += second[i] * nodes[i][c];
        }
    return result;
}

std::pair<bool, double> closest(const std::array<Point, 3>& nodes, const Point& point) {
    std::array<double, 4> polynomial{};
    for (std::size_t i = 0; i < 2; ++i) {
        const double a = 0.5 * (nodes[0][i].value() + nodes[1][i].value()) - nodes[2][i].value();
        const double b = 0.5 * (nodes[1][i].value() - nodes[0][i].value());
        const double c = nodes[2][i].value() - point[i].value();
        polynomial[0] += b * c;
        polynomial[1] += b * b + 2 * a * c;
        polynomial[2] += 3 * a * b;
        polynomial[3] += 2 * a * a;
    }
    const auto value = [&](double x) {
        return ((polynomial[3] * x + polynomial[2]) * x + polynomial[1]) * x + polynomial[0];
    };
    std::vector<double> cuts{-1.0, 1.0}, roots;
    const double discriminant = 4 * polynomial[2] * polynomial[2] - 12 * polynomial[3] * polynomial[1];
    if (polynomial[3] > 0 && discriminant >= 0)
        for (double sign : {-1.0, 1.0}) {
            const double x = (-2 * polynomial[2] + sign * std::sqrt(discriminant)) / (6 * polynomial[3]);
            if (x > -1 && x < 1)
                cuts.push_back(x);
        }
    std::sort(cuts.begin(), cuts.end());
    double scale = 0.0;
    for (auto coefficient : polynomial)
        scale += std::abs(coefficient);
    if (!(scale > 0) || !std::isfinite(scale))
        throw std::domain_error("Plane contact curve projection is degenerate");
    for (auto x : cuts)
        if (std::abs(value(x)) <= 1e-13 * scale)
            roots.push_back(x);
    for (std::size_t i = 1; i < cuts.size(); ++i) {
        double left = cuts[i - 1], right = cuts[i], fl = value(left);
        if ((fl > 0) == (value(right) > 0))
            continue;
        for (int iteration = 0; iteration < 60; ++iteration) {
            const double mid = 0.5 * (left + right), fm = value(mid);
            if ((fm > 0) == (fl > 0)) {
                left = mid;
                fl = fm;
            } else
                right = mid;
        }
        roots.push_back(0.5 * (left + right));
    }
    double best = std::numeric_limits<double>::infinity(), coordinate = 0.0;
    for (auto root : roots) {
        if ((3 * polynomial[3] * root + 2 * polynomial[2]) * root + polynomial[1] <= 0)
            continue;
        const auto projected = curve(nodes, root);
        const double dx = projected.position[0].value() - point[0].value();
        const double dy = projected.position[1].value() - point[1].value();
        const double distance = dx * dx + dy * dy;
        if (distance < best) {
            best = distance;
            coordinate = root;
        }
    }
    return {std::isfinite(best), coordinate};
}
} // namespace

Line3PlaneContactResult evaluate_line3_plane_contact(const Line3PlaneContactInput& input, bool jacobian) {
    std::array<adlite::Scalar, 22> v, rows{};
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (!std::isfinite(input.state[i]))
            throw std::domain_error("Plane contact state must be finite");
        v[i] = jacobian ? adlite::Scalar::independent(input.state[i], i, v.size()) : adlite::Scalar(input.state[i]);
    }
    std::array<Point, 3> secondary, primary;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t c = 0; c < 2; ++c) {
            secondary[i][c] = input.secondary[i][c] + v[4 + 3 * c + i];
            primary[i][c] = input.primary[i][c] + v[10 + 3 * c + i];
        }
    const auto sp = curve(secondary, input.coordinate);
    const auto projection = closest(primary, sp.position);
    Line3PlaneContactResult result;
    if (!projection.first)
        return result;
    const auto trial = curve(primary, projection.second);
    adlite::Scalar equation = 0.0;
    double derivative = 0.0;
    for (std::size_t c = 0; c < 2; ++c) {
        const auto distance = trial.position[c] - sp.position[c];
        equation += distance * trial.tangent[c];
        derivative += trial.tangent[c].value() * trial.tangent[c].value() + distance.value() * trial.second[c].value();
    }
    if (!(derivative > 0) || !std::isfinite(derivative))
        throw std::domain_error("Plane contact projection derivative must be positive");
    const auto coordinate = projection.second - (equation - equation.value()) / derivative;
    const auto pp = curve(primary, coordinate);
    if ((pp.tangent[0] * sp.tangent[0] + pp.tangent[1] * sp.tangent[1]).value() >= 0.0)
        return result;
    const auto length = adlite::hypot(pp.tangent[0], pp.tangent[1]);
    const auto secondary_length = adlite::hypot(sp.tangent[0], sp.tangent[1]);
    if (!(length.value() > 0) || !(secondary_length.value() > 0))
        throw std::domain_error("Plane contact edge has zero length");
    const Point normal{pp.tangent[1] / length, -pp.tangent[0] / length};
    const auto dx = sp.position[0] - pp.position[0], dy = sp.position[1] - pp.position[1];
    const auto gap = dx * normal[0] + dy * normal[1];
    adlite::Scalar thickness = input.secondary_thickness;
    if (input.secondary_finite)
        thickness += v[16] + v[17] * (sp.position[1] - input.secondary_reference[1])
                     - v[18] * (sp.position[0] - input.secondary_reference[0]);
    adlite::Scalar primary_thickness = input.primary_thickness;
    if (input.primary_finite)
        primary_thickness += v[19] + v[20] * (pp.position[1] - input.primary_reference[1])
                             - v[21] * (pp.position[0] - input.primary_reference[0]);
    if (!(thickness.value() > 0) || !(primary_thickness.value() > 0))
        throw std::domain_error("Plane contact thickness must remain positive");
    const auto area = input.weight * secondary_length * thickness;
    const adlite::Scalar pressure = input.mechanical && gap.value() < 0 ? -input.penalty * gap : adlite::Scalar(0);
    if (input.mechanical)
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t c = 0; c < 2; ++c) {
                const auto force = pressure * area * normal[c];
                rows[4 + 3 * c + i] -= sp.shape[i] * force;
                rows[10 + 3 * c + i] += pp.shape[i] * force;
            }
    if (input.thermal) {
        const std::array<adlite::Scalar, 2> primary_shape{0.5 * (1 - coordinate), 0.5 * (1 + coordinate)};
        const std::array<double, 2> secondary_shape{0.5 * (1 - input.coordinate), 0.5 * (1 + input.coordinate)};
        const auto ts = secondary_shape[0] * v[0] + secondary_shape[1] * v[1];
        const auto tp = primary_shape[0] * v[2] + primary_shape[1] * v[3];
        const auto rate = contact_common::gap_conductance(input.heat, gap, ts, tp) * (ts - tp) * area;
        for (std::size_t i = 0; i < 2; ++i) {
            rows[i] += secondary_shape[i] * rate;
            rows[2 + i] -= primary_shape[i] * rate;
        }
        result.heat_rate = rate.value();
    }
    result.projected = true;
    result.gap = gap.value();
    result.pressure = pressure.value();
    result.area = area.value();
    result.force = pressure.value() * area.value();
    result.distance_squared = dx.value() * dx.value() + dy.value() * dy.value();
    result.primary_coordinate = coordinate.value();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        result.residual[i] = rows[i].value();
        if (jacobian)
            rows[i].copy_derivatives(result.jacobian.data() + 22 * i, 22);
    }
    return result;
}
} // namespace fuelsim::elements
