#pragma once
#include "c3d8_types.hpp"
#include <cmath>
#include <stdexcept>

// Independent, ordinary-double reconstruction from nodal results for verification.
// This header is test support and must never be included by production code.
namespace fuelsim::test {
using RecoveryMatrix = std::array<std::array<double, 3>, 3>;

inline double recovery_determinant(const RecoveryMatrix& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

inline RecoveryMatrix recovery_inverse(const RecoveryMatrix& a) {
    const double d = recovery_determinant(a);
    if (!std::isfinite(d) || d == 0)
        throw std::domain_error("Singular reconstruction matrix");
    return {{{(a[1][1] * a[2][2] - a[1][2] * a[2][1]) / d,
                 (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / d,
                 (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / d},
        {(a[1][2] * a[2][0] - a[1][0] * a[2][2]) / d,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / d,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / d},
        {(a[1][0] * a[2][1] - a[1][1] * a[2][0]) / d,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / d,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / d}}};
}

inline RecoveryMatrix recovery_product(const RecoveryMatrix& a, const RecoveryMatrix& b) {
    RecoveryMatrix c{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                c[i][j] += a[i][k] * b[k][j];
    return c;
}

struct RecoveredC3d8Point final {
    SymmetricTensor3Values strain_increment{};
    std::array<double, 9> rotation{};
    std::array<std::array<double, 3>, 8> thermal_gradient{};
    double temperature = 0, current_weighted_measure = 0;
};

struct RecoveredC3d8 final {
    double current_volume = 0, committed_volume = 0;
    std::size_t material_point_count = 0;
    std::array<RecoveredC3d8Point, 8> points{};
};

inline RecoveredC3d8Point recover_c3d8_point(const Hex8QuadraturePoint& point,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    StrainFormulation formulation) {
    RecoveryMatrix f{}, old{}, gradient{}, rotation{};
    for (std::size_t i = 0; i < 3; ++i) {
        f[i][i] = old[i][i] = rotation[i][i] = 1;
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t n = 0; n < 8; ++n) {
                gradient[i][j] += state[8 * (i + 1) + n] * point.gradient[n][j];
                old[i][j] += committed[8 * (i + 1) + n] * point.gradient[n][j];
            }
    }
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            f[i][j] += gradient[i][j];
    const bool finite = formulation == StrainFormulation::finite;
    RecoveryMatrix h = gradient, inverse{};
    double determinant = 1;
    if (finite) {
        determinant = recovery_determinant(f);
        const double old_determinant = recovery_determinant(old);
        if (!(determinant > 0) || !std::isfinite(determinant) || !(old_determinant > 0)
            || !std::isfinite(old_determinant))
            throw std::domain_error("Invalid reconstruction configuration");
        inverse = recovery_inverse(f);
        RecoveryMatrix sum{}, difference{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                sum[i][j] = f[i][j] + old[i][j];
                difference[i][j] = 2 * (f[i][j] - old[i][j]);
            }
        if (!(recovery_determinant(sum) > 0))
            throw std::domain_error("Invalid reconstruction midpoint");
        h = recovery_product(difference, recovery_inverse(sum));
        RecoveryMatrix numerator = rotation, denominator = rotation;
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const double spin = 0.25 * (h[i][j] - h[j][i]);
                numerator[i][j] += spin;
                denominator[i][j] -= spin;
            }
        rotation = recovery_product(numerator, recovery_inverse(denominator));
    } else
        inverse = rotation;
    RecoveryMatrix strain{}, transpose{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            strain[i][j] = 0.5 * (h[i][j] + h[j][i]);
            transpose[i][j] = rotation[j][i];
        }
    strain = recovery_product(transpose, recovery_product(strain, rotation));
    RecoveredC3d8Point result;
    result.strain_increment = {strain[0][0], strain[1][1], strain[2][2], strain[0][1], strain[1][2], strain[0][2]};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            result.rotation[3 * i + j] = rotation[i][j];
    for (std::size_t n = 0; n < 8; ++n)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                result.thermal_gradient[n][j] += point.gradient[n][k] * inverse[k][j];
    result.current_weighted_measure = point.weighted_measure * determinant;
    return result;
}

inline RecoveredC3d8 recover_c3d8(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    StrainFormulation formulation,
    bool reduced) {
    RecoveredC3d8 result;
    result.material_point_count = reduced ? 1 : 8;
    std::array<RecoveredC3d8Point, 8> points;
    for (std::size_t q = 0; q < 8; ++q) {
        points[q] = recover_c3d8_point(geometry.points[q], state, committed, formulation);
        result.current_volume += points[q].current_weighted_measure;
        result.committed_volume +=
            recover_c3d8_point(geometry.points[q], committed, committed, formulation).current_weighted_measure;
    }
    if (formulation == StrainFormulation::small)
        result.current_volume = result.committed_volume = geometry.reference_volume;
    constexpr std::array<std::size_t, 8> order = {0, 1, 3, 2, 4, 5, 7, 6};
    for (std::size_t q = 0; q < result.material_point_count; ++q) {
        result.points[q] =
            reduced ? recover_c3d8_point(geometry.reduced_point, state, committed, formulation) : points[q];
        auto& p = result.points[q];
        if (!reduced)
            p.temperature = state[order[q]];
        else {
            p.thermal_gradient = {};
            for (std::size_t g = 0; g < 8; ++g)
                for (std::size_t n = 0; n < 8; ++n) {
                    p.temperature += points[g].current_weighted_measure * geometry.points[g].shape[n] * state[n]
                                     / result.current_volume;
                    for (std::size_t d = 0; d < 3; ++d)
                        p.thermal_gradient[n][d] += points[g].current_weighted_measure
                                                    * points[g].thermal_gradient[n][d] / result.current_volume;
                }
        }
    }
    return result;
}

inline RecoveredC3d8 recover_c3d8t(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    StrainFormulation formulation) {
    return recover_c3d8(geometry, state, committed, formulation, false);
}

inline RecoveredC3d8 recover_c3d8rt(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    StrainFormulation formulation) {
    return recover_c3d8(geometry, state, committed, formulation, true);
}
} // namespace fuelsim::test
