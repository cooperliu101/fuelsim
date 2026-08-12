#ifndef FUELSIM_TRANSIENT_PROBLEM_HPP
#define FUELSIM_TRANSIENT_PROBLEM_HPP

#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace fuelsim {

namespace rz {
class ProblemAccess;
class TransientConservationCalculator;
} // namespace rz
namespace cartesian3d {
class ProblemAccess;
}
struct TransientTimeErrorEstimate;
struct TransientTimeOptions;

struct TransientRegionDefinition final {
    std::string region;
    TransientInelasticProperties material;
};

struct TransientProblemDefinition final {
    SpatialDefinition spatial;
    std::vector<TransientRegionDefinition> regions;
};

struct TransientStepInput final {
    double end_time;
    double load_factor;
};

struct RegionInelasticSummary final {
    double maximum_equivalent_plastic_strain;
    double maximum_equivalent_creep_strain;
};

struct TransientConservationSummary final {
    double generated_heat_rate = 0.0;
    double stored_heat_rate = 0.0;
    double convection_heat_rate = 0.0;
    double interface_heat_imbalance = 0.0;
    double dirichlet_heat_input_rate = 0.0;
    double global_thermal_balance = 0.0;
    double relative_thermal_balance = 0.0;
    double unconstrained_thermal_residual_l2 = 0.0;

    double internal_mechanical_work_increment = 0.0;
    double pressure_traction_work_increment = 0.0;
    double dirichlet_reaction_work_increment = 0.0;
    double contact_work_increment = 0.0;
    double mechanical_work_balance = 0.0;
    double relative_mechanical_work_balance = 0.0;
    double unconstrained_mechanical_residual_l2 = 0.0;
    double elastic_energy_change = 0.0;
    double plastic_dissipation_increment = 0.0;
    double creep_dissipation_increment = 0.0;
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

// Opaque, immutable transaction snapshot used by the common time integrator.
// Its concrete material and contact layout remains owned by TransientProblem.
class TransientStateSnapshot final {
  public:
    TransientStateSnapshot();
    ~TransientStateSnapshot();
    TransientStateSnapshot(const TransientStateSnapshot& other);
    TransientStateSnapshot& operator=(const TransientStateSnapshot& other);
    TransientStateSnapshot(TransientStateSnapshot&& other) noexcept;
    TransientStateSnapshot& operator=(TransientStateSnapshot&& other) noexcept;

    bool empty() const noexcept;

  private:
    std::shared_ptr<const void> snapshot_owner() const noexcept;

    struct Storage;
    explicit TransientStateSnapshot(std::shared_ptr<const Storage> storage);

    std::shared_ptr<const Storage> _storage;
    friend class TransientProblem;
};

class TransientProblem final : public NonlinearProblem {
  public:
    TransientProblem(TransientProblemDefinition definition, const UnstructuredQuad4Mesh& source_mesh);
    TransientProblem(TransientProblemDefinition definition, const UnstructuredHex8Mesh& source_mesh);
    ~TransientProblem() override;

    const std::vector<double>& committed_solution() const noexcept;
    double committed_time() const noexcept;
    double committed_load_factor() const noexcept;
    bool time_step_active() const noexcept;
    std::vector<double> time_events() const;

    TransientStateSnapshot capture_state() const;
    void restore_state(const TransientStateSnapshot& snapshot);
    TransientTimeErrorEstimate step_doubling_error(const TransientStateSnapshot& full_step,
                                                   const TransientStateSnapshot& two_half_steps,
                                                   const TransientTimeOptions& options) const;
    void combine_last_half_step_conservation(const TransientConservationSummary& first_half);

    void begin_time_step(const TransientStepInput& input);
    void commit_time_step(const std::vector<double>& converged_solution);
    void rollback_time_step() noexcept;
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                std::size_t completed_updates);

    const TransientConservationSummary& last_conservation_summary() const noexcept;

    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    const std::vector<FieldDescriptor>& field_layout() const noexcept override;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept override;
    void validate_state(const std::vector<double>& state) const override;
    std::vector<std::size_t> required_state_dofs(std::size_t contribution_begin,
                                                 std::size_t contribution_end) const override;
    void validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                              const GlobalStateView& state) const override;
    std::size_t contribution_dof_count(std::size_t contribution_index) const override;
    void fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const override;
    void compute_contribution_residual(std::size_t contribution_index, const std::vector<double>& state,
                                       std::vector<double>& residual) const override;
    void compute_contribution_system(std::size_t contribution_index, const std::vector<double>& state,
                                     std::vector<double>& residual, std::vector<double>& jacobian) const override;

  private:
    friend class rz::ProblemAccess;
    friend class cartesian3d::ProblemAccess;
    friend class rz::TransientConservationCalculator;

    class Implementation;
    void apply_spatial_controls(double time, double load_factor);
    void clear_active_time_step() noexcept;
    void refresh_region_heat_sources();
    void require_active_time_step() const;

    std::unique_ptr<Implementation> _implementation;
};

} // namespace fuelsim

#endif
