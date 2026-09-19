#pragma once
#include "element_types.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim::elements {
std::uint64_t pellet_mesh_signature(const std::vector<CartesianPoint3>& coordinates,
    const std::vector<std::vector<std::size_t>>& connectivity,
    double conductivity);

class SurrogatePellet final {
  public:
    explicit SurrogatePellet(const std::string& path);
    void validate_mesh(const std::vector<std::size_t>& nodes,
        const std::vector<CartesianPoint3>& coordinates,
        std::uint64_t signature) const;
    void evaluate(const std::vector<double>& temperature, double heat_source, std::vector<double>& residual) const;
    void evaluate_with_jacobian(const std::vector<double>& temperature,
        double heat_source,
        std::vector<double>& residual,
        std::vector<double>& jacobian) const;

    double volume() const noexcept { return _volume; }

  private:
    struct Layer final {
        std::size_t input, output;
        bool tanh;
        std::vector<double> weights, biases;
    };

    void compute(const std::vector<double>& temperature,
        double heat_source,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const;
    std::vector<Layer> _layers;
    std::vector<std::size_t> _nodes;
    std::vector<CartesianPoint3> _coordinates;
    std::vector<double> _input_mean, _input_scale, _output_mean, _output_scale;
    std::uint64_t _signature = 0;
    double _volume = 0;
};
} // namespace fuelsim::elements
