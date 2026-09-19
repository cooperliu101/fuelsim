#include "pellet_mlp.hpp"
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace fuelsim::elements {
std::uint64_t pellet_mesh_signature(const std::vector<CartesianPoint3>& coordinates,
    const std::vector<std::vector<std::size_t>>& connectivity,
    double conductivity) {
    std::uint64_t hash = 14695981039346656037ULL;
    const auto integer = [&](std::uint64_t value) {
        for (int i = 0; i < 8; ++i) {
            hash ^= value & 255U;
            hash *= 1099511628211ULL;
            value >>= 8;
        }
    };
    const auto real = [&](double value) {
        std::uint64_t bits;
        static_assert(sizeof(bits) == sizeof(value), "Pellet signature requires binary64");
        std::memcpy(&bits, &value, sizeof(bits));
        integer(bits);
    };
    real(conductivity);
    integer(connectivity.size());
    for (const auto& element : connectivity) {
        integer(element.size());
        for (auto node : element) {
            integer(node);
            const auto& p = coordinates.at(node);
            real(p.x);
            real(p.y);
            real(p.z);
        }
    }
    return hash;
}

SurrogatePellet::SurrogatePellet(const std::string& path) {
    std::ifstream input(path);
    std::string magic;
    int version = 0;
    std::size_t n = 0, layers = 0;
    if (!(input >> magic >> version >> n >> layers >> _volume >> _signature) || magic != "FUELSIM_PELLET_MLP"
        || version != 1 || n == 0 || n > 4096 || layers == 0 || layers > 16 || !(_volume > 0)
        || !std::isfinite(_volume))
        throw std::invalid_argument("Invalid pellet model header: " + path);
    _nodes.resize(n);
    _coordinates.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        auto& p = _coordinates[i];
        if (!(input >> _nodes[i] >> p.x >> p.y >> p.z) || !std::isfinite(p.x) || !std::isfinite(p.y)
            || !std::isfinite(p.z) || (i > 0 && _nodes[i] <= _nodes[i - 1]))
            throw std::invalid_argument("Invalid pellet model node ordering");
    }
    const auto read = [&](std::vector<double>& values, std::size_t count, bool positive) {
        values.resize(count);
        for (auto& v : values)
            if (!(input >> v) || !std::isfinite(v) || (positive && !(v > 0)))
                throw std::invalid_argument("Invalid pellet model coefficient or normalization");
    };
    read(_input_mean, n + 1, false);
    read(_input_scale, n + 1, true);
    read(_output_mean, n, false);
    read(_output_scale, n, true);
    std::size_t width = n + 1;
    for (std::size_t l = 0; l < layers; ++l) {
        Layer layer{};
        std::string activation;
        if (!(input >> layer.input >> layer.output >> activation) || layer.input != width || layer.output == 0
            || layer.output > 4096 || (activation != "tanh" && activation != "linear")
            || (l + 1 == layers && (layer.output != n || activation != "linear")))
            throw std::invalid_argument("Invalid pellet network layer");
        layer.tanh = activation == "tanh";
        read(layer.weights, layer.input * layer.output, false);
        read(layer.biases, layer.output, false);
        width = layer.output;
        _layers.push_back(std::move(layer));
    }
    std::string extra;
    if (input >> extra)
        throw std::invalid_argument("Unexpected trailing pellet model data");
}

void SurrogatePellet::validate_mesh(const std::vector<std::size_t>& nodes,
    const std::vector<CartesianPoint3>& coordinates,
    std::uint64_t signature) const {
    if (nodes != _nodes || signature != _signature)
        throw std::invalid_argument("Pellet model mesh, material or node order mismatch");
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const auto& p = coordinates.at(nodes[i]);
        const auto& q = _coordinates[i];
        if (p.x != q.x || p.y != q.y || p.z != q.z)
            throw std::invalid_argument("Pellet model surface coordinates mismatch");
    }
}

void SurrogatePellet::compute(const std::vector<double>& temperature,
    double heat_source,
    std::vector<double>& residual,
    std::vector<double>* jacobian) const {
    const auto n = _nodes.size();
    if (temperature.size() != n || !std::isfinite(heat_source))
        throw std::invalid_argument("Invalid pellet network input");
    std::vector<double> values = temperature;
    values.push_back(heat_source);
    std::vector<double> derivatives;
    if (jacobian)
        derivatives.assign((n + 1) * n, 0);
    for (std::size_t i = 0; i <= n; ++i) {
        if (!std::isfinite(values[i]))
            throw std::domain_error("Nonfinite pellet network temperature");
        values[i] = (values[i] - _input_mean[i]) / _input_scale[i];
        if (jacobian && i < n)
            derivatives[i * n + i] = 1 / _input_scale[i];
    }
    for (const auto& layer : _layers) {
        std::vector<double> next = layer.biases, chain;
        if (jacobian)
            chain.assign(layer.output * n, 0);
        for (std::size_t i = 0; i < layer.output; ++i) {
            for (std::size_t j = 0; j < layer.input; ++j) {
                const double w = layer.weights[i * layer.input + j];
                next[i] += w * values[j];
                if (jacobian)
                    for (std::size_t k = 0; k < n; ++k)
                        chain[i * n + k] += w * derivatives[j * n + k];
            }
            if (!std::isfinite(next[i]))
                throw std::domain_error("Nonfinite pellet network activation");
            if (layer.tanh) {
                next[i] = std::tanh(next[i]);
                if (jacobian)
                    for (std::size_t k = 0; k < n; ++k)
                        chain[i * n + k] *= 1 - next[i] * next[i];
            }
        }
        values = std::move(next);
        derivatives = std::move(chain);
    }
    residual.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        residual[i] = values[i] * _output_scale[i] + _output_mean[i];
        if (!std::isfinite(residual[i]))
            throw std::domain_error("Nonfinite pellet network residual");
        if (jacobian)
            for (std::size_t j = 0; j < n; ++j)
                derivatives[i * n + j] *= _output_scale[i];
    }
    if (jacobian) {
        for (double value : derivatives)
            if (!std::isfinite(value))
                throw std::domain_error("Nonfinite pellet network tangent");
        *jacobian = std::move(derivatives);
    }
}

void SurrogatePellet::evaluate(const std::vector<double>& temperature,
    double heat_source,
    std::vector<double>& residual) const {
    compute(temperature, heat_source, residual, nullptr);
}

void SurrogatePellet::evaluate_with_jacobian(const std::vector<double>& temperature,
    double heat_source,
    std::vector<double>& residual,
    std::vector<double>& jacobian) const {
    compute(temperature, heat_source, residual, &jacobian);
}
} // namespace fuelsim::elements
