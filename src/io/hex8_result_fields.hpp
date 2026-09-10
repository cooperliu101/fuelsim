#pragma once
#include "c3d8_types.hpp"
#include <array>

namespace fuelsim::io_detail {
inline constexpr std::array<const char*, 24> hex8_derived_field_names = {"reference_x",
    "reference_y",
    "reference_z",
    "current_x",
    "current_y",
    "current_z",
    "material_temperature",
    "integration_measure",
    "heat_flux_x",
    "heat_flux_y",
    "heat_flux_z",
    "logarithmic_strain_xx",
    "logarithmic_strain_yy",
    "logarithmic_strain_zz",
    "logarithmic_strain_xy",
    "logarithmic_strain_yz",
    "logarithmic_strain_xz",
    "infinitesimal_strain_xx",
    "infinitesimal_strain_yy",
    "infinitesimal_strain_zz",
    "infinitesimal_strain_xy",
    "infinitesimal_strain_yz",
    "infinitesimal_strain_xz",
    "current_measure"};

// Postprocessing only: no constitutive update, history mutation, or seeded derivatives.
std::array<std::array<double, 24>, 8> hex8_derived_results(const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const IsotropicThermoelasticMaterial& material,
    StrainFormulation formulation,
    bool reduced,
    double time);
} // namespace fuelsim::io_detail
