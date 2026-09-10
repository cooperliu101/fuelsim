#pragma once
#include "c3d20_types.hpp"

namespace fuelsim {
Hex20LocalResidual compute_hex20_thermoelastic(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues* committed_state,
    double time_step,
    Hex20LocalJacobian* jacobian,
    bool include_thermal_time_term);
Hex20LocalResidual compute_hex20_transient(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex20LocalJacobian* jacobian,
    bool include_thermal_time_term);
CartesianMaterialHistory compute_hex20_transient_update(const elements::C3d20Input& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step);
std::vector<SymmetricTensor3Values>
compute_hex20_stress(const elements::C3d20Input& data, const Hex20Geometry& geometry, const Hex20LocalValues& state);
} // namespace fuelsim
