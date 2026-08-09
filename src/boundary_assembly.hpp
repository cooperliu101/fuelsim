#ifndef FUELSIM_BOUNDARY_ASSEMBLY_HPP
#define FUELSIM_BOUNDARY_ASSEMBLY_HPP

#include "fuelsim/boundary.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "spatial_layout.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim {

class SpatialAssembly;

class BoundaryAssembly final {
  public:
    BoundaryAssembly() = default;

    void set_load_factor(double load_factor, const SpatialLayout& layout);
    double load_factor() const noexcept;
    void set_time(double time, const SpatialLayout& layout);
    double region_heat_source(std::size_t region_index,
                              const SpatialLayout& layout) const;

  private:
    friend class SpatialAssembly;

    double function_value(const std::string& name,
                          const SpatialLayout& layout) const;
    double load_multiplier(bool scale_with_load,
                           const std::string& function,
                           const SpatialLayout& layout) const;
    void refresh_controlled_values(const SpatialLayout& layout);

    struct PressureLoad final {
        double pressure;
        bool scale_with_load;
        std::string function;
    };
    struct TractionLoad final {
        double traction;
        bool scale_with_load;
        std::string function;
    };
    struct ControlledDirichlet final {
        std::size_t dof;
        double value;
        bool scale_with_load;
        std::string function;
    };
    struct ConvectionLoad final {
        double heat_transfer_coefficient;
        double ambient_temperature;
        std::string coefficient_function;
        std::string ambient_temperature_function;
    };
    struct PressureContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzPressureGeometry geometry;
    };
    struct TractionContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzTractionGeometry geometry;
    };
    struct ConvectionContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzConvectionGeometry geometry;
    };

    std::vector<Line2RzConvectionKernel> _convection_kernels;
    std::vector<ConvectionLoad> _convection_loads;
    std::vector<ConvectionContribution> _convection_contributions;
    std::vector<Line2RzPressureKernel> _pressure_kernels;
    std::vector<PressureContribution> _pressure_contributions;
    std::vector<Line2RzTractionKernel> _traction_kernels;
    std::vector<TractionContribution> _traction_contributions;
    std::vector<DirichletCondition> _dirichlet_conditions;
    std::vector<ControlledDirichlet> _controlled_dirichlet_conditions;
    std::vector<PressureLoad> _pressure_loads;
    std::vector<TractionLoad> _traction_loads;
    double _load_factor = 1.0;
    double _time = 0.0;
};

} // namespace fuelsim

#endif
