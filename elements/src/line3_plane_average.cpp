#include "line3_plane.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
using Active = adlite::Scalar;
using Point = std::array<Active, 2>;

struct Sample final {
    double coordinate, weight;
};

std::array<Sample, 11> sampling_rule() {
    std::array<Sample, 11> result{};
    std::size_t index = 0;
    // Equal-weight degree-three rules in the consistent nodal-area intervals.
    // These formulas reproduce the independent native gap and transfer probes;
    // no measured contact matrix coefficients enter the implementation.
    for (std::size_t patch = 0; patch < 3; ++patch) {
        const double left = patch == 0 ? -1.0 : (patch == 1 ? -2.0 / 3.0 : 2.0 / 3.0);
        const double right = patch == 0 ? -2.0 / 3.0 : (patch == 1 ? 2.0 / 3.0 : 1.0);
        const std::size_t count = patch == 1 ? 5 : 3;
        const double n = static_cast<double>(count);
        for (std::size_t q = 0; q < count; ++q)
            result[index++] = {0.5 * (left + right)
                                   + (right - left) * (static_cast<double>(q) - 0.5 * (n - 1.0))
                                         / std::sqrt(n * n - 1.0),
                (right - left) / n};
    }
    return result;
}

double test_shape(std::size_t node, double coordinate) {
    if (node == 0)
        return std::max(0.0, -coordinate);
    if (node == 1)
        return std::max(0.0, coordinate);
    if (node == 2)
        return 1.0 - std::abs(coordinate);
    throw std::invalid_argument("Plane contact local constraint node must be less than three");
}

struct EdgePoint final {
    Point position{}, tangent{};
    std::array<Active, 3> shape;
};

EdgePoint interpolate(const std::vector<Point>& nodes, const std::array<std::size_t, 3>& edge, const Active& x) {
    EdgePoint result;
    result.shape = {0.5 * x * (x - 1.0), 0.5 * x * (x + 1.0), 1.0 - x * x};
    const std::array<Active, 3> derivative{x - 0.5, x + 0.5, -2.0 * x};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t c = 0; c < 2; ++c) {
            result.position[c] += result.shape[i] * nodes.at(edge[i])[c];
            result.tangent[c] += derivative[i] * nodes.at(edge[i])[c];
        }
    return result;
}

Active project(const std::vector<Point>& nodes,
    const std::array<std::size_t, 3>& edge,
    const Point& point,
    const Point& tangent) {
    Active a = 0.0, b = 0.0, c = 0.0;
    for (std::size_t d = 0; d < 2; ++d) {
        a += (0.5 * (nodes[edge[0]][d] + nodes[edge[1]][d]) - nodes[edge[2]][d]) * tangent[d];
        b += 0.5 * (nodes[edge[1]][d] - nodes[edge[0]][d]) * tangent[d];
        c += (nodes[edge[2]][d] - point[d]) * tangent[d];
    }
    const double av = a.value(), bv = b.value(), cv = c.value();
    const double scale = std::abs(av) + std::abs(bv) + std::abs(cv);
    double root = 0.0;
    if (!(scale > 0.0) || !std::isfinite(scale))
        throw std::domain_error("Averaged plane contact projection is degenerate");
    if (std::abs(av) < 1e-14 * scale) {
        if (!(bv < 0.0))
            throw std::domain_error("Averaged plane contact surfaces have inconsistent orientation");
        root = -cv / bv;
    } else {
        const double discriminant = bv * bv - 4.0 * av * cv;
        if (!(discriminant > 0.0))
            throw std::domain_error("Averaged plane contact normal does not intersect the primary curve");
        // Select the root whose tangent opposes the secondary tangent.
        root = bv < 0.0 ? -2.0 * cv / (bv - std::sqrt(discriminant)) : (-bv - std::sqrt(discriminant)) / (2.0 * av);
    }
    const double derivative = 2.0 * av * root + bv;
    if (!(derivative < -1e-12 * scale))
        throw std::domain_error("Averaged plane contact projection is singular");
    const Active equation = (a * root + b) * root + c;
    return root - (equation - equation.value()) / derivative;
}
} // namespace

PlaneSurfaceContactResult evaluate_plane_averaged_contact(const PlaneAveragedContactGeometry& geometry,
    const std::vector<double>& state,
    bool jacobian) {
    const std::size_t count = geometry.coordinates.size(), size = 2 * count;
    if (state.size() != size || geometry.secondary.empty() || geometry.primary.empty() || !(geometry.penalty > 0.0)
        || !std::isfinite(geometry.penalty))
        throw std::invalid_argument("Invalid averaged plane contact input");
    std::vector<Point> nodes(count);
    for (std::size_t n = 0; n < count; ++n)
        for (std::size_t c = 0; c < 2; ++c) {
            const auto i = 2 * n + c;
            if (!std::isfinite(state[i]))
                throw std::domain_error("Averaged plane contact state must be finite");
            nodes[n][c] =
                geometry.coordinates[n][c] + (jacobian ? Active::independent(state[i], i, size) : Active(state[i]));
        }
    const auto samples = sampling_rule();
    Active area = 0.0, gap_integral = 0.0;
    std::vector<Point> secondary_vectors(count), primary_vectors(count);
    for (const auto& edge : geometry.secondary) {
        if (!(edge.thickness > 0.0))
            throw std::invalid_argument("Averaged plane contact requires positive initial thickness");
        // The first equal-weight sample is the physical smoothing half-width.
        const Active radius = (1.0 - 1.0 / std::sqrt(2.0)) / 12.0
                              * adlite::hypot(nodes[edge.nodes[1]][0] - nodes[edge.nodes[0]][0],
                                  nodes[edge.nodes[1]][1] - nodes[edge.nodes[0]][1]);
        if (!(radius.value() > 0.0))
            throw std::domain_error("Averaged plane contact secondary edge is degenerate");
        for (const auto& sample : samples) {
            const double test = test_shape(edge.local_node, sample.coordinate);
            if (test == 0.0)
                continue;
            const auto s = interpolate(nodes, edge.nodes, sample.coordinate);
            const Active length = adlite::hypot(s.tangent[0], s.tangent[1]);
            if (!(length.value() > 0.0))
                throw std::domain_error("Averaged plane contact has an undefined secondary normal");
            const Point tangent{s.tangent[0] / length, s.tangent[1] / length};
            const Point normal{-tangent[1], tangent[0]};
            const Active weight = sample.weight * test * edge.thickness * length;
            Active fraction_sum = 0.0;
            for (const auto& primary : geometry.primary) {
                Active left = 0.0, right = 0.0;
                for (std::size_t c = 0; c < 2; ++c) {
                    left += (nodes[primary[1]][c] - s.position[c]) * tangent[c];
                    right += (nodes[primary[0]][c] - s.position[c]) * tangent[c];
                }
                if (!(right.value() > left.value()))
                    throw std::domain_error("Averaged plane contact primary edge is reversed or folded");
                if (left.value() >= radius.value() || right.value() <= -radius.value())
                    continue;
                const Active a = left.value() > -radius.value() ? left : -radius;
                const Active b = right.value() < radius.value() ? right : radius;
                const Active fraction = (b - a) / (2.0 * radius);
                fraction_sum += fraction;
                const auto p = interpolate(nodes, primary, project(nodes, primary, s.position, tangent));
                Active gap = 0.0;
                for (std::size_t c = 0; c < 2; ++c)
                    gap += (s.position[c] - p.position[c]) * normal[c];
                const Active measure = weight * fraction;
                area += measure;
                gap_integral += measure * gap;
                // Average shape times normal, retaining the normal field inside
                // the integral. Factoring out one mean normal loses curvature.
                for (std::size_t n = 0; n < 3; ++n)
                    for (std::size_t c = 0; c < 2; ++c) {
                        secondary_vectors[edge.nodes[n]][c] += measure * s.shape[n] * normal[c];
                        primary_vectors[primary[n]][c] += measure * p.shape[n] * normal[c];
                    }
            }
            if (fraction_sum.value() > 1.0 + 1e-10)
                throw std::domain_error("Averaged plane contact primary intervals overlap");
        }
    }
    if (!(area.value() > 0.0))
        throw std::domain_error("Averaged plane contact constraint has no primary projection");
    const Active gap = gap_integral / area;
    const Active pressure = gap.value() < 0.0 ? -geometry.penalty * gap : Active(0.0);
    PlaneSurfaceContactResult result;
    result.residual.resize(size);
    if (jacobian)
        result.jacobian.resize(size * size);
    for (std::size_t n = 0; n < count; ++n)
        for (std::size_t c = 0; c < 2; ++c) {
            const Active row = pressure * (primary_vectors[n][c] - secondary_vectors[n][c]);
            const auto i = 2 * n + c;
            result.residual[i] = row.value();
            if (jacobian)
                row.copy_derivatives(result.jacobian.data() + i * size, size);
        }
    result.point.projected = true;
    result.point.gap = gap.value();
    result.point.pressure = pressure.value();
    result.point.area = area.value();
    result.point.force = pressure.value() * area.value();
    return result;
}
} // namespace fuelsim::elements
