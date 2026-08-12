#ifndef FUELSIM_HEX8_THERMOELASTIC_HPP
#define FUELSIM_HEX8_THERMOELASTIC_HPP

#include "fuelsim/hex8.hpp"
#include "fuelsim/material.hpp"

#include <array>

namespace fuelsim {

class Hex8ThermoelasticKernel final {
  public:
    Hex8ThermoelasticKernel(IsotropicThermoelasticMaterial material, double volumetric_heat_source);
    Hex8ThermoelasticKernel(IsotropicThermoelasticMaterial material, double constant_heat_capacity, double volumetric_heat_source);

    double volumetric_heat_source() const noexcept;
    double heat_capacity(double temperature, double x, double y, double z) const;
    void set_volumetric_heat_source(double volumetric_heat_source) noexcept;
    void set_time(double time) noexcept;

    Hex8LocalResidual residual(const Hex8Geometry& geometry, const Hex8LocalValues& state) const;
    Hex8LocalSystem linearize(const Hex8Geometry& geometry, const Hex8LocalValues& state) const;
    Hex8LocalResidual residual(const Hex8Geometry& geometry, const Hex8LocalValues& current_state, const Hex8LocalValues& committed_state, double time_step) const;
    Hex8LocalSystem linearize(const Hex8Geometry& geometry, const Hex8LocalValues& current_state, const Hex8LocalValues& committed_state, double time_step) const;
    std::array<SymmetricTensor3Values, 8> stress_values(const Hex8Geometry& geometry, const Hex8LocalValues& state) const;

  private:
    void residual_ad(const Hex8Geometry& geometry, const Hex8LocalAdValues& state, const Hex8LocalValues* committed_state, double time_step, Hex8LocalAdValues& residual) const;

    IsotropicThermoelasticMaterial _material;
    double _constant_heat_capacity;
    double _volumetric_heat_source;
    double _time;
};

enum class CartesianTractionComponent { x, y, z };

class Quad4FacePressureKernel final {
  public:
    explicit Quad4FacePressureKernel(double pressure);
    void set_pressure(double pressure) noexcept;
    Quad4FaceLocalResidual residual(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;
    Quad4FaceLocalSystem linearize(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;

  private:
    double _pressure;
};

class Quad4FaceTractionKernel final {
  public:
    Quad4FaceTractionKernel(CartesianTractionComponent component, double traction);
    void set_traction(double traction) noexcept;
    Quad4FaceLocalResidual residual(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;
    Quad4FaceLocalSystem linearize(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;

  private:
    CartesianTractionComponent _component;
    double _traction;
};

class Quad4FaceConvectionKernel final {
  public:
    Quad4FaceConvectionKernel(double heat_transfer_coefficient, double ambient_temperature);
    void set_properties(double heat_transfer_coefficient, double ambient_temperature) noexcept;
    Quad4FaceLocalResidual residual(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;
    Quad4FaceLocalSystem linearize(const Quad4FaceGeometry& geometry, const Quad4FaceLocalValues& state) const;

  private:
    double _heat_transfer_coefficient;
    double _ambient_temperature;
};

} // namespace fuelsim

#endif
