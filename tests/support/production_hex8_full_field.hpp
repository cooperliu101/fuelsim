#pragma once
#include <cstddef>
#include <string>

namespace fuelsim::test {
struct ProductionHex8FullFieldOptions final {
    std::string case_name;
    std::string reference_prefix;
    std::size_t expected_steps = 0;
    double time_step = 0.0;
    double bulk_relative_tolerance = 1.0e-2;
    double bulk_pointwise_relative_tolerance = 0.0;
    double reaction_heat_flux_pointwise_relative_tolerance = 0.0;
    double reaction_heat_flux_pointwise_absolute_tolerance = 0.0;
    double displacement_pointwise_relative_tolerance = 0.0;
    double displacement_pointwise_absolute_tolerance = 0.0;
    double reaction_pointwise_relative_tolerance = 0.0;
    double reaction_pointwise_absolute_tolerance = 0.0;
    double stress_pointwise_relative_tolerance = 0.0;
    double stress_pointwise_absolute_tolerance = 0.0;
    double logarithmic_strain_pointwise_relative_tolerance = 0.0;
    double logarithmic_strain_pointwise_absolute_tolerance = 0.0;
    double elastic_strain_pointwise_relative_tolerance = 0.0;
    double elastic_strain_pointwise_absolute_tolerance = 0.0;
    double inelastic_pointwise_relative_tolerance = 0.0;
    double reaction_zero_absolute_tolerance = 1.0;
    double energy_relative_tolerance = 1.0e-2;
    double energy_pointwise_relative_tolerance = 0.0;
    double hourglass_energy_pointwise_absolute_tolerance = 0.0;
    double external_work_pointwise_absolute_tolerance = 0.0;
    double coordinate_tolerance = 1.0e-7;
    bool reduced_integration = false;
    bool gate_reaction_heat_flux = true;
};

bool compare_production_hex8_full_field(const std::string& output_path, const ProductionHex8FullFieldOptions& options);
} // namespace fuelsim::test
