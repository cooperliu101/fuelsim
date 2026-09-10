#include "cax8_geometry.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
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
        throw std::domain_error("CAX8T requires positive reference Jacobian and integration-point radius");
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
        throw std::domain_error("CAX8T linear source geometry must be positive");
    p.source_measure = 2 * std::acos(-1.0) * p.source_radius * source_det * weight;
    for (std::size_t n = 0; n < 4; ++n) {
        p.source_gradient_r[n] = (tx[n] * d - ty[n] * c) / source_det;
        p.source_gradient_z[n] = (ty[n] * a - tx[n] * b) / source_det;
    }
    return p;
}

Quad8RzGeometry make_quad8_rz_geometry(const Quad8RzCoordinates& coordinates, std::size_t order) {
    Quad8RzGeometry result{coordinates, {}};
    if (order == 2) {
        result.point_count = 4;
        const double g = 1 / std::sqrt(3.0);
        const std::array<double, 4> xi = {-g, g, g, -g}, eta = {-g, -g, g, g};
        for (std::size_t q = 0; q < 4; ++q)
            result.points[q] = evaluate_quad8_rz_point(coordinates, xi[q], eta[q], 1);
        return result;
    }
    if (order != 3)
        throw std::invalid_argument("QUAD8 geometry requires cax8t or cax8rt");
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> q = {-g, 0, g}, w = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            result.points[3 * j + i] = evaluate_quad8_rz_point(coordinates, q[i], q[j], w[i] * w[j]);
    return result;
}

} // namespace fuelsim
