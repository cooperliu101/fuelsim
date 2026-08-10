#ifndef FUELSIM_TRANSIENT_PROBLEM_HPP
#define FUELSIM_TRANSIENT_PROBLEM_HPP

#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz_transient.hpp"
#include "fuelsim/spatial_definition.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace fuelsim {

class SpatialAssembly;
class TransientConservationCalculator;

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

struct TransientCommittedState final {
    std::vector<double> solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> material_stresses;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    TransientConservationSummary conservation;
    double time = 0.0;
    double load_factor = 0.0;
};

class TransientProblem final : public NonlinearProblem {
  public:
    TransientProblem(TransientProblemDefinition definition, const UnstructuredQuad4Mesh& source_mesh);
    ~TransientProblem() override;

    const TransientProblemDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;
    std::size_t region_count() const noexcept;
    std::size_t region_index(const std::string& name) const;
    std::size_t region_node_offset(std::size_t region_index) const;
    const RegionDefinition& region(std::size_t region_index) const;
    const RegionMesh& region_mesh(std::size_t region_index) const;
    const Quad4RzTransientKernel& region_kernel(std::size_t region_index) const;
    const Quad4RzGeometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;

    const std::vector<double>& committed_solution() const noexcept;
    double committed_time() const noexcept;
    double committed_load_factor() const noexcept;
    bool time_step_active() const noexcept;
    std::vector<double> time_events() const;

    TransientCommittedState committed_state() const;
    void restore_committed_state(TransientCommittedState state);

    void begin_time_step(const TransientStepInput& input);
    void commit_time_step(const std::vector<double>& converged_solution);
    void rollback_time_step() noexcept;
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                std::size_t completed_updates);

    const Quad4MaterialHistory& material_history(std::size_t region_index, std::size_t element_index) const;
    const std::array<AxisymmetricStressValues, 4>& material_stress(std::size_t region_index,
                                                                   std::size_t element_index) const;
    RegionInelasticSummary summarize_region_history(std::size_t region_index) const;
    const TransientConservationSummary& last_conservation_summary() const noexcept;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;
    std::vector<ContactNodeSummary> summarize_contact_nodes(std::size_t contact_index,
                                                            const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;

    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept override;
    void validate_state(const std::vector<double>& state) const override;
    std::vector<std::size_t> required_state_dofs(std::size_t contribution_begin,
                                                 std::size_t contribution_end) const override;
    void validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                              const GlobalStateView& state) const override;
    LocalDofs contribution_dofs(std::size_t contribution_index) const override;
    LocalResidual contribution_residual(std::size_t contribution_index, const LocalValues& state) const override;
    LocalSystem linearize_contribution(std::size_t contribution_index, const LocalValues& state) const override;

  private:
    friend class TransientConservationCalculator;

    void apply_spatial_controls(double time, double load_factor);
    void clear_active_time_step() noexcept;
    void refresh_region_heat_sources();
    void require_active_time_step() const;

    TransientProblemDefinition _definition;
    std::unique_ptr<SpatialAssembly> _spatial;
    std::vector<Quad4RzTransientKernel> _region_kernels;
    std::vector<std::vector<Quad4MaterialHistory>> _material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> _material_stresses;
    TransientConservationSummary _last_conservation_summary;
    std::vector<double> _committed_solution;
    double _committed_time;
    double _committed_load_factor;
    double _active_time_step;
    double _active_end_time;
    double _active_load_factor;
    std::vector<std::vector<ContactPointHistory>> _active_contact_histories;
    bool _time_step_active;
};

} // namespace fuelsim

#endif
