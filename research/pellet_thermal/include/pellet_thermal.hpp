#pragma once
#include <cstddef>
#include <vector>

namespace fuelsim::elements {
// Matrices use row-major storage; source is the nodal load for unit W/m^3.
class ExactCondensedPellet final {
  public:
    ExactCondensedPellet() = default;
    ExactCondensedPellet(const std::vector<double>& stiffness,
        const std::vector<double>& source,
        const std::vector<std::size_t>& boundary);
    void evaluate(const std::vector<double>& temperature, double heat_source, std::vector<double>& residual) const;
    void evaluate_with_jacobian(const std::vector<double>& temperature,
        double heat_source,
        std::vector<double>& residual,
        std::vector<double>& jacobian) const;

    const std::vector<double>& stiffness() const noexcept { return _stiffness; }

    const std::vector<double>& source() const noexcept { return _source; }

  private:
    std::vector<double> _stiffness, _source;
};
} // namespace fuelsim::elements
