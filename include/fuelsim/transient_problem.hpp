#pragma once
#include "fuelsim/nonlinear_problem.hpp"
#include <array>
#include <cstddef>
#include <memory>
#include <vector>
namespace fuelsim {
struct TransientTimeErrorEstimate;
struct TransientTimeOptions;
struct SpatialDefinition;
class SpatialProblemStorage;
class BackendAccess;
struct TransientStepInput final {
    double end_time, load_factor;
};
struct RegionStateSummary final {
    double maximum_temperature, maximum_equivalent_plastic_strain, maximum_equivalent_creep_strain;
};
struct TransientConservationSummary final {
    double generated_heat_rate = 0.0, stored_heat_rate = 0.0, convection_heat_rate = 0.0,
           interface_heat_imbalance = 0.0, dirichlet_heat_input_rate = 0.0, global_thermal_balance = 0.0,
           relative_thermal_balance = 0.0, unconstrained_thermal_residual_l2 = 0.0,
           internal_mechanical_work_increment = 0.0, pressure_traction_work_increment = 0.0,
           dirichlet_reaction_work_increment = 0.0, contact_work_increment = 0.0, mechanical_work_balance = 0.0,
           relative_mechanical_work_balance = 0.0, unconstrained_mechanical_residual_l2 = 0.0,
           elastic_energy_change = 0.0, plastic_dissipation_increment = 0.0, creep_dissipation_increment = 0.0;
};
struct TransientConservationField final {
    const char* name;
    double TransientConservationSummary::* member;
};
inline constexpr std::array<TransientConservationField, 18> transient_conservation_fields = {{
    {"generated_heat_rate", &TransientConservationSummary::generated_heat_rate},
    {"stored_heat_rate", &TransientConservationSummary::stored_heat_rate},
    {"convection_heat_rate", &TransientConservationSummary::convection_heat_rate},
    {"interface_heat_imbalance", &TransientConservationSummary::interface_heat_imbalance},
    {"dirichlet_heat_input_rate", &TransientConservationSummary::dirichlet_heat_input_rate},
    {"global_thermal_balance", &TransientConservationSummary::global_thermal_balance},
    {"relative_thermal_balance", &TransientConservationSummary::relative_thermal_balance},
    {"unconstrained_thermal_residual_l2", &TransientConservationSummary::unconstrained_thermal_residual_l2},
    {"internal_mechanical_work_increment", &TransientConservationSummary::internal_mechanical_work_increment},
    {"pressure_traction_work_increment", &TransientConservationSummary::pressure_traction_work_increment},
    {"dirichlet_reaction_work_increment", &TransientConservationSummary::dirichlet_reaction_work_increment},
    {"contact_work_increment", &TransientConservationSummary::contact_work_increment},
    {"mechanical_work_balance", &TransientConservationSummary::mechanical_work_balance},
    {"relative_mechanical_work_balance", &TransientConservationSummary::relative_mechanical_work_balance},
    {"unconstrained_mechanical_residual_l2", &TransientConservationSummary::unconstrained_mechanical_residual_l2},
    {"elastic_energy_change", &TransientConservationSummary::elastic_energy_change},
    {"plastic_dissipation_increment", &TransientConservationSummary::plastic_dissipation_increment},
    {"creep_dissipation_increment", &TransientConservationSummary::creep_dissipation_increment},
}};
class TransientProblem final : public NonlinearProblem {
  public:
    TransientProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh);
    TransientProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);
    ~TransientProblem() override;
    bool is_cartesian_3d() const noexcept;
    const SpatialDefinition& definition() const noexcept;
    const std::vector<double>& committed_solution() const noexcept;
    double committed_time() const noexcept;
    double committed_load_factor() const noexcept;
    bool time_step_active() const noexcept;
    std::vector<double> time_events() const;
    RegionStateSummary summarize_region(std::size_t region) const;
    ProblemStateSnapshot capture_state() const;
    void restore_state(const ProblemStateSnapshot& snapshot);
    TransientTimeErrorEstimate step_doubling_error(const ProblemStateSnapshot& full_step,
        const ProblemStateSnapshot& two_half_steps, const TransientTimeOptions& options) const;
    void combine_last_half_step_conservation(const TransientConservationSummary& first_half);
    void begin_time_step(const TransientStepInput& input);
    void commit_time_step(const std::vector<double>& converged_solution);
    void rollback_time_step() noexcept;
    bool uses_augmented_contact() const noexcept override;
    AugmentedContactUpdate update_augmented_contact_multipliers(
        const std::vector<double>& state, std::size_t completed_updates) override;
    const TransientConservationSummary& last_conservation_summary() const noexcept;
    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    const std::vector<FieldDescriptor>& field_layout() const noexcept override;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept override;
    void validate_state(const std::vector<double>& state) const override;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;
    void compute_contribution(std::size_t index, const std::vector<double>& state, std::vector<double>& residual,
        std::vector<double>* jacobian) const override;

  private:
    friend class BackendAccess;
    void apply_spatial_controls(double time, double load_factor);
    std::vector<double> accumulate_contribution_conservation(
        const std::vector<double>& solution, TransientConservationSummary& summary) const;
    void clear_active_time_step() noexcept;
    void require_active_time_step() const;
    std::unique_ptr<SpatialProblemStorage> _impl;
};
} // namespace fuelsim
