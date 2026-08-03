#ifndef FUELSIM_BOUNDARY_HPP
#define FUELSIM_BOUNDARY_HPP

#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"

#include <array>
#include <cstddef>

namespace fuelsim {

struct Line2RzConvectionGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};

struct ConvectionProperties final {
    double heat_transfer_coefficient;
    double ambient_temperature;
};

struct Line2RzPressureGeometry final {
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};

struct PressureProperties final {
    double pressure;
    bool use_displaced_geometry;
};

Line2RzConvectionGeometry make_line2_rz_convection_geometry(
    const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes);

Line2RzPressureGeometry make_line2_rz_pressure_geometry(
    const std::array<RzPoint, 2>& coordinates,
    const std::array<std::size_t, 2>& local_nodes);

class Line2RzPressureKernel final {
  public:
    explicit Line2RzPressureKernel(PressureProperties properties);

    const PressureProperties& properties() const noexcept;
    void set_properties(PressureProperties properties);
    LocalResidual residual(const Line2RzPressureGeometry& geometry,
                           const LocalValues& state) const;
    LocalSystem linearize(const Line2RzPressureGeometry& geometry,
                          const LocalValues& state) const;

  private:
    void residual_ad(const Line2RzPressureGeometry& geometry,
                     const LocalAdValues& state, LocalAdValues& residual) const;

    PressureProperties _properties;
};

class Line2RzConvectionKernel final {
  public:
    explicit Line2RzConvectionKernel(ConvectionProperties properties);

    const ConvectionProperties& properties() const noexcept;
    void set_properties(ConvectionProperties properties);
    LocalResidual residual(const Line2RzConvectionGeometry& geometry,
                           const LocalValues& state) const;
    LocalSystem linearize(const Line2RzConvectionGeometry& geometry,
                          const LocalValues& state) const;

  private:
    void residual_ad(const Line2RzConvectionGeometry& geometry,
                     const LocalAdValues& state, LocalAdValues& residual) const;

    ConvectionProperties _properties;
};

} // namespace fuelsim

#endif
