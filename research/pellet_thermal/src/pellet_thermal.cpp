#include "pellet_thermal.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
ExactCondensedPellet::ExactCondensedPellet(const std::vector<double>& stiffness,
    const std::vector<double>& source,
    const std::vector<std::size_t>& boundary) {
    const auto n = source.size(), b = boundary.size();
    if (b == 0 || b > n || stiffness.size() != n * n)
        throw std::invalid_argument("Invalid pellet condensation dimensions");
    for (double v : stiffness)
        if (!std::isfinite(v))
            throw std::invalid_argument("Nonfinite pellet stiffness");
    for (double v : source)
        if (!std::isfinite(v))
            throw std::invalid_argument("Nonfinite pellet source");
    std::vector<bool> selected(n, false);
    for (auto i : boundary) {
        if (i >= n || selected[i])
            throw std::invalid_argument("Invalid or repeated pellet surface node");
        selected[i] = true;
    }
    std::vector<std::size_t> interior;
    for (std::size_t i = 0; i < n; ++i)
        if (!selected[i])
            interior.push_back(i);
    const auto m = interior.size();
    std::vector<double> lower(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double v = stiffness[interior[i] * n + interior[j]];
            for (std::size_t k = 0; k < j; ++k)
                v -= lower[i * m + k] * lower[j * m + k];
            if (i == j) {
                if (!(v > 0.0) || !std::isfinite(v))
                    throw std::invalid_argument("Pellet interior conductivity matrix is not positive definite");
                lower[i * m + j] = std::sqrt(v);
            } else
                lower[i * m + j] = v / lower[j * m + j];
        }
    _stiffness.resize(b * b);
    _source.resize(b);
    // Solve K_II X = [K_IG, f_I], without forming an inverse.
    for (std::size_t col = 0; col <= b; ++col) {
        std::vector<double> x(m);
        for (std::size_t i = 0; i < m; ++i) {
            x[i] = col == b ? source[interior[i]] : stiffness[interior[i] * n + boundary[col]];
            for (std::size_t j = 0; j < i; ++j)
                x[i] -= lower[i * m + j] * x[j];
            x[i] /= lower[i * m + i];
        }
        for (std::size_t i = m; i-- > 0;) {
            for (std::size_t j = i + 1; j < m; ++j)
                x[i] -= lower[j * m + i] * x[j];
            x[i] /= lower[i * m + i];
        }
        for (std::size_t i = 0; i < b; ++i) {
            double v = col == b ? source[boundary[i]] : stiffness[boundary[i] * n + boundary[col]];
            for (std::size_t j = 0; j < m; ++j)
                v -= stiffness[boundary[i] * n + interior[j]] * x[j];
            if (col == b)
                _source[i] = v;
            else
                _stiffness[i * b + col] = v;
        }
    }
}

void ExactCondensedPellet::evaluate(const std::vector<double>& temperature,
    double heat_source,
    std::vector<double>& residual) const {
    const auto n = _source.size();
    if (temperature.size() != n || n == 0 || !std::isfinite(heat_source))
        throw std::invalid_argument("Invalid pellet response input");
    for (double v : temperature)
        if (!std::isfinite(v))
            throw std::domain_error("Nonfinite pellet temperature");
    residual.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        residual[i] = -heat_source * _source[i];
        for (std::size_t j = 0; j < n; ++j)
            residual[i] += _stiffness[i * n + j] * temperature[j];
    }
}

void ExactCondensedPellet::evaluate_with_jacobian(const std::vector<double>& temperature,
    double heat_source,
    std::vector<double>& residual,
    std::vector<double>& jacobian) const {
    evaluate(temperature, heat_source, residual);
    jacobian = _stiffness;
}
} // namespace fuelsim::elements
