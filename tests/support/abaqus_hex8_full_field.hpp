#ifndef FUELSIM_TEST_ABAQUS_HEX8_FULL_FIELD_HPP
#define FUELSIM_TEST_ABAQUS_HEX8_FULL_FIELD_HPP

#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include <string>
#include <vector>

namespace fuelsim::test {
struct AbaqusHex8StepSnapshot final {
    double time = 0.0;
    std::vector<double> state;
    std::vector<std::vector<CartesianMaterialPointState>> material_by_source_element;
    std::vector<CartesianContactNodeSummary> contact;
    std::vector<ContactPointHistory> contact_history;
    TransientConservationSummary conservation;
};

class AbaqusHex8SnapshotObserver final : public TransientStepObserver {
  public:
    void accepted_step(const TransientProblem& problem, const TransientAcceptedStep& step) override;

    const std::vector<AbaqusHex8StepSnapshot>& snapshots() const noexcept { return _snapshots; }

  private:
    std::vector<AbaqusHex8StepSnapshot> _snapshots;
};

struct AbaqusHex8FullFieldOptions final {
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
    double contact_relative_tolerance = 1.25e-2;
    double contact_pointwise_relative_tolerance = 0.0;
    double contact_slip_pointwise_relative_tolerance = 0.0;
    double contact_slip_pointwise_absolute_tolerance = 0.0;
    double contact_replayed_heat_rate_relative_tolerance = 0.0;
    double contact_replayed_heat_rate_pointwise_relative_tolerance = 0.0;
    double contact_replayed_heat_rate_pointwise_absolute_tolerance = 0.0;
    double contact_total_heat_rate_relative_tolerance = 0.0;
    double energy_relative_tolerance = 1.0e-2;
    double energy_pointwise_relative_tolerance = 0.0;
    double external_work_pointwise_absolute_tolerance = 0.0;
    double coordinate_tolerance = 1.0e-7;
    double minimum_contact_state_match_fraction = 1.0;
    double tangent_basis_tolerance = 1.0e-9;
    bool use_contact_summary_total_slip = false;
    bool gate_contact_slip = true;
    bool gate_contact_pressure = true;
    bool gate_contact_state = true;
};

bool compare_abaqus_hex8_full_field(const TransientProblem& solved_problem, const SpatialDefinition& definition,
    const UnstructuredHex8Mesh& mesh, const std::vector<AbaqusHex8StepSnapshot>& snapshots,
    const AbaqusHex8FullFieldOptions& options);
} // namespace fuelsim::test

#endif
