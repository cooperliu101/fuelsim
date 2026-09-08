#include "io/hex8_result_fields.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::io_detail {
namespace {
using Matrix3 = std::array<std::array<double, 3>, 3>;

SymmetricTensor3Values logarithmic_strain(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    Matrix3 deformation = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                deformation[component][direction] +=
                    state[8 * (component + 1) + node] * point.gradient[node][direction];
    Matrix3 left{};
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t column = 0; column < 3; ++column)
            for (std::size_t inner = 0; inner < 3; ++inner)
                left[row][column] += deformation[row][inner] * deformation[column][inner];
    Matrix3 vectors = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (std::size_t sweep = 0; sweep < 32; ++sweep) {
        std::size_t p = 0, q = 1;
        double largest = std::abs(left[p][q]);
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t column = row + 1; column < 3; ++column)
                if (std::abs(left[row][column]) > largest) {
                    p = row;
                    q = column;
                    largest = std::abs(left[row][column]);
                }
        const double scale = std::max({1.0, std::abs(left[0][0]), std::abs(left[1][1]), std::abs(left[2][2])});
        if (largest <= 1.0e-15 * scale)
            break;
        const double tau = (left[q][q] - left[p][p]) / (2.0 * left[p][q]);
        const double tangent = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent), sine = tangent * cosine;
        const double app = left[p][p], aqq = left[q][q], apq = left[p][q];
        left[p][p] = app - tangent * apq;
        left[q][q] = aqq + tangent * apq;
        left[p][q] = left[q][p] = 0.0;
        for (std::size_t row = 0; row < 3; ++row) {
            if (row == p || row == q)
                continue;
            const double arp = left[row][p], arq = left[row][q];
            left[row][p] = left[p][row] = cosine * arp - sine * arq;
            left[row][q] = left[q][row] = sine * arp + cosine * arq;
        }
        for (std::size_t row = 0; row < 3; ++row) {
            const double vrp = vectors[row][p], vrq = vectors[row][q];
            vectors[row][p] = cosine * vrp - sine * vrq;
            vectors[row][q] = sine * vrp + cosine * vrq;
        }
    }
    Matrix3 result{};
    for (std::size_t mode = 0; mode < 3; ++mode) {
        if (!(left[mode][mode] > 0.0))
            throw std::domain_error("HEX8 logarithmic strain output requires a positive left stretch tensor");
        const double value = 0.5 * std::log(left[mode][mode]);
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t column = 0; column < 3; ++column)
                result[row][column] += value * vectors[row][mode] * vectors[column][mode];
    }
    return {result[0][0], result[1][1], result[2][2], result[0][1], result[1][2], result[0][2]};
}

Matrix3 deformation_gradient(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    Matrix3 result = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t node = 0; node < 8; ++node)
                result[i][j] += state[8 * (i + 1) + node] * point.gradient[node][j];
    return result;
}

double determinant(const Matrix3& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

Matrix3 inverse(const Matrix3& a, double det) {
    Matrix3 result{};
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t column = 0; column < 3; ++column) {
            const auto i = (column + 1) % 3, j = (column + 2) % 3;
            const auto k = (row + 1) % 3, l = (row + 2) % 3;
            result[row][column] = (a[i][k] * a[j][l] - a[i][l] * a[j][k]) / det;
        }
    return result;
}

struct PointGeometry final {
    std::array<std::array<double, 3>, 8> gradient{};
    double measure;
};

PointGeometry current_geometry(const Hex8QuadraturePoint& point, const Hex8LocalValues& state, bool finite) {
    PointGeometry result{point.gradient, point.weighted_measure};
    if (!finite)
        return result;
    const auto f = deformation_gradient(point, state);
    const double det = determinant(f);
    if (!(det > 0.0) || !std::isfinite(det))
        throw std::domain_error("HEX8 output requires a positive finite current Jacobian");
    const auto inv = inverse(f, det);
    result.measure *= det;
    result.gradient = {};
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t current = 0; current < 3; ++current)
            for (std::size_t reference = 0; reference < 3; ++reference)
                result.gradient[node][current] += point.gradient[node][reference] * inv[reference][current];
    return result;
}

std::array<double, 6> components(const SymmetricTensor3Values& tensor) {
    return {tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz};
}
} // namespace

std::array<std::array<double, 24>, 8> hex8_derived_results(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation formulation,
    bool reduced,
    double time) {
    const bool finite = formulation == StrainFormulation::finite;
    std::array<PointGeometry, 8> current{};
    double volume = 0.0, average_trace = 0.0;
    std::array<std::array<double, 3>, 8> average_gradient{};
    double average_temperature = 0.0;
    for (std::size_t q = 0; q < 8; ++q) {
        const auto& point = geometry.points[q];
        current[q] = current_geometry(point, state, finite);
        volume += current[q].measure;
        for (std::size_t node = 0; node < 8; ++node) {
            average_temperature += current[q].measure * point.shape[node] * state[node];
            for (std::size_t component = 0; component < 3; ++component) {
                average_gradient[node][component] += current[q].measure * current[q].gradient[node][component];
                average_trace +=
                    point.weighted_measure * point.gradient[node][component] * state[8 * (component + 1) + node];
            }
        }
    }
    if (!(volume > 0.0) || !std::isfinite(volume))
        throw std::domain_error("HEX8 output requires a positive finite integration volume");
    average_temperature /= volume;
    average_trace /= geometry.reference_volume;
    for (auto& gradient : average_gradient)
        for (double& component : gradient)
            component /= volume;
    constexpr std::array<std::size_t, 8> material_node = {0, 1, 3, 2, 4, 5, 7, 6};
    std::array<std::array<double, 24>, 8> result{};
    for (std::size_t q = 0; q < 8; ++q) {
        const auto& point = reduced ? geometry.reduced_point : geometry.points[q];
        auto& values = result[q];
        values[0] = values[3] = point.position.x;
        values[1] = values[4] = point.position.y;
        values[2] = values[5] = point.position.z;
        for (std::size_t node = 0; node < 8; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                values[3 + component] += point.shape[node] * state[8 * (component + 1) + node];
        const double temperature = reduced ? average_temperature : state[material_node[q]];
        values[6] = temperature;
        values[7] = reduced ? volume : point.weighted_measure / geometry.reference_volume * volume;
        const auto& gradient = reduced ? average_gradient : current[q].gradient;
        const double conductivity =
            material
                .conductivity(adlite::Scalar(temperature), {time, point.position.x, point.position.y, point.position.z})
                .value();
        for (std::size_t node = 0; node < 8; ++node)
            for (std::size_t component = 0; component < 3; ++component)
                values[8 + component] -= conductivity * gradient[node][component] * state[node];
        std::array<double, 6> logarithmic{};
        try {
            logarithmic = components(logarithmic_strain(point, state));
        } catch (const std::domain_error&) {
            if (finite)
                throw;
            // The linearized small-strain problem does not constrain the
            // total deformation gradient. An undefined derived logarithm
            // must not make an otherwise valid small-strain output fail.
            logarithmic.fill(std::numeric_limits<double>::quiet_NaN());
        }
        if (finite && !reduced) {
            const double correction =
                (std::log(volume / geometry.reference_volume) - logarithmic[0] - logarithmic[1] - logarithmic[2]) / 3.0;
            for (std::size_t component = 0; component < 3; ++component)
                logarithmic[component] += correction;
        }
        const auto f = deformation_gradient(point, state);
        std::array<double, 6> infinitesimal = {f[0][0] - 1,
            f[1][1] - 1,
            f[2][2] - 1,
            0.5 * (f[0][1] + f[1][0]),
            0.5 * (f[1][2] + f[2][1]),
            0.5 * (f[0][2] + f[2][0])};
        if (!reduced) {
            const double correction = (average_trace - infinitesimal[0] - infinitesimal[1] - infinitesimal[2]) / 3.0;
            for (std::size_t component = 0; component < 3; ++component)
                infinitesimal[component] += correction;
        }
        for (std::size_t component = 0; component < 6; ++component) {
            values[11 + component] = logarithmic[component];
            values[17 + component] = infinitesimal[component];
        }
        values[23] = current_geometry(point, state, finite).measure;
    }
    return result;
}
} // namespace fuelsim::io_detail
