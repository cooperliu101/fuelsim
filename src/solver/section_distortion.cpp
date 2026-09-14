#include "solver/section_distortion.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace fuelsim {
namespace {
double dot(const std::vector<double>& a, const std::vector<double>& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

std::vector<double> section_action(const CrossSection& section, const std::vector<double>& field, bool stiffness) {
    const auto n = section.nodes().size();
    std::vector<double> result(2 * n, 0.0);
    for (const auto& point : section.points()) {
        double ux = 0.0, uy = 0.0, exx = 0.0, eyy = 0.0, gamma = 0.0;
        for (std::size_t i = 0; i < 8; ++i) {
            const double x = field[point.nodes[i]], y = field[n + point.nodes[i]];
            ux += point.shape[i] * x;
            uy += point.shape[i] * y;
            exx += point.gradient[i][0] * x;
            eyy += point.gradient[i][1] * y;
            gamma += point.gradient[i][1] * x + point.gradient[i][0] * y;
        }
        const auto& region = section.regions()[point.region];
        const auto properties =
            region.material.active_properties(region.temperature, material_context(0.0, point.position));
        const double mu = properties.shear_modulus.value(), lambda = properties.lame_lambda.value();
        const double sx = (lambda + 2.0 * mu) * exx + lambda * eyy;
        const double sy = lambda * exx + (lambda + 2.0 * mu) * eyy;
        const double tau = mu * gamma;
        for (std::size_t i = 0; i < 8; ++i) {
            result[point.nodes[i]] +=
                point.weight
                * (stiffness ? point.gradient[i][0] * sx + point.gradient[i][1] * tau : point.shape[i] * ux);
            result[n + point.nodes[i]] +=
                point.weight
                * (stiffness ? point.gradient[i][1] * sy + point.gradient[i][0] * tau : point.shape[i] * uy);
        }
    }
    return result;
}

// Twice-reorthogonalized modified Gram-Schmidt in the section displacement inner product.
bool append_section_direction(const CrossSection& section,
    std::vector<double> candidate,
    std::vector<std::vector<double>>& basis,
    std::vector<std::vector<double>>& weighted) {
    auto image = section_action(section, candidate, false);
    const double initial = dot(candidate, image);
    for (int pass = 0; pass < 2; ++pass)
        for (std::size_t j = 0; j < basis.size(); ++j) {
            const double projection = dot(weighted[j], candidate);
            for (std::size_t i = 0; i < candidate.size(); ++i) {
                candidate[i] -= projection * basis[j][i];
                image[i] -= projection * weighted[j][i];
            }
        }
    // Recompute the image instead of trusting subtraction near a dependent direction.
    image = section_action(section, candidate, false);
    const double norm = dot(candidate, image);
    if (norm <= 1.0e-20 * initial)
        return false;
    for (std::size_t i = 0; i < candidate.size(); ++i) {
        candidate[i] /= std::sqrt(norm);
        image[i] /= std::sqrt(norm);
    }
    basis.push_back(std::move(candidate));
    weighted.push_back(std::move(image));
    return true;
}

struct ProjectedSpectrum final {
    std::vector<double> values, vectors;
};

// Specialized small offline section pencil: S v = mu D v, with D positive definite.
// Cholesky triangular solves whiten D; cyclic symmetric Jacobi rotations solve the
// resulting real symmetric problem. No explicit matrix inverse is constructed.
ProjectedSpectrum solve_section_pencil(std::vector<double> shear, std::vector<double> distortion, std::size_t n) {
    double scale = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        scale = std::max(scale, distortion[n * i + i]);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double value = distortion[n * i + j];
            for (std::size_t k = 0; k < j; ++k)
                value -= distortion[n * i + k] * distortion[n * j + k];
            if (i == j) {
                if (!std::isfinite(value) || value <= 1.0e-13 * scale)
                    throw std::runtime_error("Projected section distortion matrix is singular or ill-conditioned");
                distortion[n * i + j] = std::sqrt(value);
            } else
                distortion[n * i + j] = value / distortion[n * j + j];
        }
    for (std::size_t column = 0; column < n; ++column)
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < i; ++k)
                shear[n * i + column] -= distortion[n * i + k] * shear[n * k + column];
            shear[n * i + column] /= distortion[n * i + i];
        }
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < j; ++k)
                shear[n * row + j] -= distortion[n * j + k] * shear[n * row + k];
            shear[n * row + j] /= distortion[n * j + j];
        }
    ProjectedSpectrum result;
    result.vectors.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        result.vectors[n * i + i] = 1.0;
        for (std::size_t j = 0; j < i; ++j) {
            const double average = 0.5 * (shear[n * i + j] + shear[n * j + i]);
            shear[n * i + j] = average;
            shear[n * j + i] = average;
        }
    }
    bool converged = false;
    for (int sweep = 0; sweep < 80; ++sweep) {
        double diagonal = 0.0, off = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            diagonal = std::hypot(diagonal, shear[n * i + i]);
            for (std::size_t j = 0; j < i; ++j)
                off = std::hypot(off, shear[n * i + j]);
        }
        if (off <= 1.0e-13 * diagonal) {
            converged = true;
            break;
        }
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = p + 1; q < n; ++q) {
                const double apq = shear[n * p + q];
                if (std::abs(apq) < 1.0e-16 * diagonal)
                    continue;
                const double delta = (shear[n * q + q] - shear[n * p + p]) / (2.0 * apq);
                const double t = std::copysign(1.0, delta) / (std::abs(delta) + std::hypot(1.0, delta));
                const double c = 1.0 / std::hypot(1.0, t), s = t * c;
                shear[n * p + p] -= t * apq;
                shear[n * q + q] += t * apq;
                shear[n * p + q] = shear[n * q + p] = 0.0;
                for (std::size_t k = 0; k < n; ++k) {
                    if (k != p && k != q) {
                        const double kp = shear[n * k + p], kq = shear[n * k + q];
                        shear[n * k + p] = shear[n * p + k] = c * kp - s * kq;
                        shear[n * k + q] = shear[n * q + k] = s * kp + c * kq;
                    }
                    const double kp = result.vectors[n * k + p], kq = result.vectors[n * k + q];
                    result.vectors[n * k + p] = c * kp - s * kq;
                    result.vectors[n * k + q] = s * kp + c * kq;
                }
            }
    }
    if (!converged)
        throw std::runtime_error("Projected section eigensolver did not converge");
    result.values.resize(n);
    for (std::size_t j = 0; j < n; ++j) {
        result.values[j] = shear[n * j + j];
        for (std::size_t offset = 0; offset < n; ++offset) {
            const std::size_t i = n - 1 - offset;
            for (std::size_t k = i + 1; k < n; ++k)
                result.vectors[n * i + j] -= distortion[n * k + i] * result.vectors[n * k + j];
            result.vectors[n * i + j] /= distortion[n * i + i];
        }
    }
    return result;
}
} // namespace

SectionDistortionSpectrum
build_section_distortion_modes(const CrossSection& section, std::size_t count, std::size_t shear_free_count) {
    const auto nodes = section.nodes().size(), full = 2 * nodes;
    if (full > 512 || count + shear_free_count == 0 || count + shear_free_count > full - 3)
        throw std::invalid_argument("Section distortion requires 1..(2*n-3) modes and at most 512 transverse DOFs");
    std::vector<std::vector<double>> basis, weighted;
    for (std::size_t global = 0; global < 3; ++global) {
        std::vector<double> field(full);
        for (std::size_t i = 0; i < nodes; ++i) {
            if (global == 0)
                field[i] = 1.0;
            else if (global == 1)
                field[nodes + i] = 1.0;
            else {
                field[i] = -(section.nodes()[i].y - section.centroid().y);
                field[nodes + i] = section.nodes()[i].x - section.centroid().x;
            }
        }
        if (!append_section_direction(section, std::move(field), basis, weighted))
            throw std::runtime_error("Classical section rigid space is linearly dependent");
    }
    for (std::size_t i = 0; i < full; ++i) {
        std::vector<double> field(full);
        field[i] = 1.0;
        (void)append_section_direction(section, std::move(field), basis, weighted);
    }
    if (basis.size() != full)
        throw std::runtime_error("Section complement failed the rank check");
    const auto n = full - 3;
    std::vector<double> d(n * n), s(n * n);
    SectionWarpingSolver warping(section);
    for (std::size_t column = 0; column < n; ++column) {
        const auto di = section_action(section, basis[3 + column], true);
        const auto si = warping.solve(basis[3 + column]).condensed_force;
        for (std::size_t row = 0; row < n; ++row) {
            d[n * row + column] = dot(basis[3 + row], di);
            s[n * row + column] = dot(basis[3 + row], si);
        }
    }
    const auto eigen = solve_section_pencil(s, d, n);
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&eigen](std::size_t a, std::size_t b) {
        return eigen.values[a] > eigen.values[b];
    });
    const double maximum = eigen.values[order.front()];
    SectionDistortionSpectrum result;
    for (double value : eigen.values) {
        if (value < -1.0e-9 * maximum)
            throw std::runtime_error("Condensed section shear operator is not positive semidefinite");
        if (value <= 1.0e-10 * maximum)
            ++result.shear_kernel_dimension;
    }
    if (count > n - result.shear_kernel_dimension)
        throw std::invalid_argument("Requested distortion count exceeds the finite shear spectrum");
    for (std::size_t selected = 0; selected < count; ++selected) {
        const auto j = order[selected];
        SectionDistortionMode mode;
        mode.eigenvalue = 1.0 / eigen.values[j];
        mode.distortion_energy = mode.eigenvalue;
        mode.transverse.resize(full);
        for (std::size_t k = 0; k < n; ++k)
            for (std::size_t i = 0; i < full; ++i)
                mode.transverse[i] += basis[3 + k][i] * eigen.vectors[n * k + j];
        mode.warping = warping.solve(mode.transverse);
        const double norm = std::sqrt(mode.warping.stiffness);
        const auto pivot = std::max_element(mode.transverse.begin(), mode.transverse.end(), [](double a, double b) {
            return std::abs(a) < std::abs(b);
        });
        const double factor = std::copysign(1.0 / norm, *pivot);
        for (double& value : mode.transverse)
            value *= factor;
        mode.warping = warping.solve(mode.transverse);
        mode.shear_norm = std::sqrt(mode.warping.stiffness);
        mode.section_norm = std::sqrt(dot(mode.transverse, section_action(section, mode.transverse, false)));
        for (std::size_t global = 0; global < 3; ++global)
            mode.classic_projection =
                std::max(mode.classic_projection, std::abs(dot(weighted[global], mode.transverse)) / mode.section_norm);
        const auto di = section_action(section, mode.transverse, true);
        double residual = 0.0, denominator = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            const double kd = dot(basis[3 + k], di);
            const double ks = mode.eigenvalue * dot(basis[3 + k], mode.warping.condensed_force);
            residual = std::hypot(residual, kd - ks);
            denominator = std::hypot(denominator, std::abs(kd) + std::abs(ks));
        }
        mode.eigen_relative_residual = residual / denominator;
        for (const auto& previous : result.modes)
            mode.orthogonality_error =
                std::max(mode.orthogonality_error, std::abs(dot(previous.transverse, mode.warping.condensed_force)));
        if (mode.eigen_relative_residual > 1.0e-7 || mode.classic_projection > 1.0e-9
            || mode.orthogonality_error > 1.0e-7)
            throw std::runtime_error("Section eigenmode failed residual or orthogonality verification");
        result.modes.push_back(std::move(mode));
    }
    if (shear_free_count > result.shear_kernel_dimension)
        throw std::invalid_argument("Requested shear-free count exceeds the section shear kernel");
    if (shear_free_count != 0) {
        // KS has non-rigid null vectors: gradients which the scalar warping space
        // can cancel exactly. Discarding them removes, for example, lateral
        // dilation needed to relax restrained Poisson contraction. Their pencil
        // eigenvalue is infinite, not a spurious small finite eigenvalue.
        // Within this already identified kernel ONLY, order by KD/W. The pencil
        // eigenvectors are D-orthonormal, while the complement basis is W-orthonormal.
        const auto k = result.shear_kernel_dimension;
        std::vector<double> mass(k * k), identity(k * k);
        for (std::size_t a = 0; a < k; ++a) {
            identity[k * a + a] = 1.0;
            for (std::size_t b = 0; b < k; ++b)
                for (std::size_t i = 0; i < n; ++i)
                    mass[k * a + b] +=
                        eigen.vectors[n * i + order[n - k + a]] * eigen.vectors[n * i + order[n - k + b]];
        }
        const auto kernel = solve_section_pencil(mass, identity, k);
        std::vector<std::size_t> kernel_order(k);
        std::iota(kernel_order.begin(), kernel_order.end(), 0);
        std::stable_sort(kernel_order.begin(), kernel_order.end(), [&kernel](std::size_t a, std::size_t b) {
            return kernel.values[a] > kernel.values[b];
        });
        std::vector<std::vector<double>> kernel_fields, kernel_weights;
        for (std::size_t selected = 0; selected < k; ++selected) {
            const auto j = kernel_order[selected];
            SectionDistortionMode mode;
            mode.eigenvalue = std::numeric_limits<double>::infinity();
            mode.distortion_energy = 1.0 / kernel.values[j];
            mode.transverse.resize(full);
            for (std::size_t i = 0; i < n; ++i) {
                double coefficient = 0.0;
                for (std::size_t a = 0; a < k; ++a)
                    coefficient += eigen.vectors[n * i + order[n - k + a]] * kernel.vectors[k * a + j];
                for (std::size_t node = 0; node < full; ++node)
                    mode.transverse[node] += basis[3 + i][node] * coefficient;
            }
            const double norm = std::sqrt(dot(mode.transverse, section_action(section, mode.transverse, false)));
            const auto pivot = std::max_element(mode.transverse.begin(), mode.transverse.end(), [](double a, double b) {
                return std::abs(a) < std::abs(b);
            });
            const double factor = std::copysign(1.0 / norm, *pivot);
            for (double& value : mode.transverse)
                value *= factor;
            kernel_fields.push_back(std::move(mode.transverse));
            kernel_weights.push_back(section_action(section, kernel_fields.back(), false));
        }
        std::vector<std::vector<double>> selected_fields, selected_weights;
        for (const auto& field : kernel_fields)
            if (selected_fields.size() < shear_free_count)
                (void)append_section_direction(section, field, selected_fields, selected_weights);
        if (selected_fields.size() != shear_free_count)
            throw std::runtime_error("Shear-free selection failed its rank check");
        for (auto& field : selected_fields) {
            SectionDistortionMode mode;
            mode.eigenvalue = std::numeric_limits<double>::infinity();
            mode.transverse = std::move(field);
            const auto weight = section_action(section, mode.transverse, false);
            mode.distortion_energy = dot(mode.transverse, section_action(section, mode.transverse, true));
            mode.section_norm = std::sqrt(dot(mode.transverse, weight));
            mode.warping = warping.solve(mode.transverse);
            mode.shear_norm = std::sqrt(mode.warping.stiffness);
            for (std::size_t global = 0; global < 3; ++global)
                mode.classic_projection =
                    std::max(mode.classic_projection, std::abs(dot(weighted[global], mode.transverse)));
            for (const auto& previous : result.shear_free_modes)
                mode.orthogonality_error =
                    std::max(mode.orthogonality_error, std::abs(dot(previous.transverse, weight)));
            // A scale-independent original-operator check, not just a check on
            // eigenvalues of the numerically transformed pencil.
            const auto di = section_action(section, mode.transverse, true);
            mode.eigen_relative_residual = std::sqrt(dot(mode.warping.condensed_force, mode.warping.condensed_force))
                                           / (maximum * std::sqrt(dot(di, di)));
            if (mode.classic_projection > 1.0e-9 || mode.orthogonality_error > 1.0e-7
                || mode.eigen_relative_residual > 1.0e-7)
                throw std::runtime_error("Shear-free section mode failed the original-operator check");
            result.shear_free_modes.push_back(std::move(mode));
        }
    }
    return result;
}
} // namespace fuelsim
