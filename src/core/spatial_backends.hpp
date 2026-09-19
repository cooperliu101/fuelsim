#pragma once
#include "core/problem_backend_access.hpp"
#include "cpeg8t.hpp"
#include "spatial_backend.hpp"

namespace fuelsim {
class ThermalBackend final : public SpatialBackend {
  public:
    ThermalBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredQuad4Mesh& mesh,
        bool transient);
    ThermalBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredQuad8Mesh& mesh,
        bool transient);
    ThermalBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredHex8Mesh& mesh,
        bool transient);
    ThermalBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredHex20Mesh& mesh,
        bool transient);

    const thermal::SpatialAssembly& spatial() const noexcept { return _spatial; }

    const spatial_detail::SpatialLayout& layout() const noexcept override { return _spatial; }

    std::size_t contribution_count() const noexcept override { return _spatial.contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const override {
        return _spatial.contribution_type(index);
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    void validate_state(const std::vector<double>& state) const override { _spatial.validate_state(state); }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override {
        return _spatial.required_state_dofs(first, last);
    }

    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override {
        _spatial.validate_local_state(first, last, state);
    }

    void set_time(double value) override;
    void set_load_factor(double value) override;
    void set_heat_source_interval(double begin, double end) override;
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) override;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const override;
    void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) override;
    void publish_material_history() noexcept override;
    void export_histories(TransientCommittedState& state) const override;
    void restore_histories(TransientCommittedState& state) override;
    RegionStateSummary summarize_region(std::size_t region) const override;
    std::vector<double> committed_creep_rates() const override;

    bool contribution_metadata_is_fixed() const noexcept override { return true; }

  private:
    thermal::SpatialAssembly _spatial;
};

class RadialBackend final : public SpatialBackend {
  public:
    RadialBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredBar2Mesh& mesh,
        bool transient);

    const radial::SpatialAssembly& spatial() const noexcept { return _spatial; }

    const spatial_detail::SpatialLayout& layout() const noexcept override { return _spatial; }

    std::size_t contribution_count() const noexcept override { return _spatial.contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const override {
        return _spatial.contribution_type(index);
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    void validate_state(const std::vector<double>& state) const override { _spatial.validate_state(state); }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override {
        return _spatial.required_state_dofs(first, last);
    }

    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override {
        _spatial.validate_local_state(first, last, state);
    }

    void set_time(double value) override;
    void set_load_factor(double value) override;
    void set_heat_source_interval(double begin, double end) override;
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) override;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const override;
    void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) override;
    void publish_material_history() noexcept override;
    void export_histories(TransientCommittedState& state) const override;
    void restore_histories(TransientCommittedState& state) override;
    RegionStateSummary summarize_region(std::size_t region) const override;
    std::vector<double> committed_creep_rates() const override;

    bool contribution_metadata_is_fixed() const noexcept override { return true; }

    bool combined_contact_contribution() const noexcept override { return true; }

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept override {
        return _spatial.committed_contact_histories();
    }

    void commit_contact_state(const std::vector<double>& state) override;

    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories) override {
        _spatial.restore_contact_state(state, std::move(histories));
    }

    const std::vector<std::vector<Cax2tGpsMaterialHistory>>& histories() const noexcept { return _histories; }

    void stage_contact_history(const std::vector<double>& state) override {
        _staged_contact_solution = state;
        _staged_friction_dissipation = _spatial.prepare_contact_state(state, _staged_contact_histories);
    }

    void publish_contact_history() noexcept override {
        _spatial.publish_contact_state(_staged_contact_solution, _staged_contact_histories);
    }

    void complete_contact_diagnostics(TransientConservationSummary& conservation) const override {
        conservation.friction_dissipation_increment = _staged_friction_dissipation;
    }

  private:
    radial::SpatialAssembly _spatial;
    std::vector<std::vector<Cax2tGpsMaterialHistory>> _histories, _staged;
    double _time = 0.0;
};

class Rz4Backend final : public SpatialBackend {
  public:
    Rz4Backend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredQuad4Mesh& mesh,
        bool transient);

    const rz::SpatialAssembly& spatial() const noexcept { return _spatial; }

    const spatial_detail::SpatialLayout& layout() const noexcept override { return _spatial; }

    std::size_t contribution_count() const noexcept override { return _spatial.contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const override {
        return _spatial.contribution_type(index);
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    void validate_state(const std::vector<double>& state) const override { _spatial.validate_state(state); }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override {
        return _spatial.required_state_dofs(first, last);
    }

    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override {
        _spatial.validate_local_state(first, last, state);
    }

    void set_time(double value) override;
    void set_load_factor(double value) override;
    void set_heat_source_interval(double begin, double end) override;
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) override;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const override;
    void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) override;
    void publish_material_history() noexcept override;
    void export_histories(TransientCommittedState& state) const override;
    void restore_histories(TransientCommittedState& state) override;
    RegionStateSummary summarize_region(std::size_t region) const override;
    std::vector<double> committed_creep_rates() const override;

    std::size_t sparsity_contribution_count() const noexcept override { return _spatial.sparsity_contribution_count(); }

    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept override {
        return _spatial.committed_contact_histories();
    }

    void commit_contact_state(const std::vector<double>& state) override;

    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories) override {
        _spatial.restore_contact_state(state, std::move(histories));
    }

    bool uses_augmented_contact() const noexcept override { return _spatial.uses_augmented_contact(); }

    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
        std::size_t completed_updates) override {
        return _spatial.update_augmented_contact_multipliers(state, completed_updates);
    }

    const std::vector<AxisymmetricRegionData>& kernel_data() const noexcept { return _kernel_data; }

    const std::vector<std::vector<Quad4MaterialHistory>>& histories() const noexcept { return _histories; }

    void stage_contact_history(const std::vector<double>& state) override {
        _staged_contact_solution = state;
        _spatial.prepare_contact_state(state, _staged_contact_histories);
    }

    void publish_contact_history() noexcept override {
        _spatial.publish_contact_state(_staged_contact_solution, _staged_contact_histories);
    }

  private:
    rz::SpatialAssembly _spatial;
    std::vector<std::vector<Quad4MaterialHistory>> _histories, _staged;
    std::vector<AxisymmetricRegionData> _kernel_data;
};

class Rz8Backend final : public SpatialBackend {
  public:
    Rz8Backend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredQuad8Mesh& mesh,
        bool transient);

    const rz8::SpatialAssembly& spatial() const noexcept { return _spatial; }

    const spatial_detail::SpatialLayout& layout() const noexcept override { return _spatial; }

    std::size_t contribution_count() const noexcept override { return _spatial.contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const override {
        return _spatial.contribution_type(index);
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    void validate_state(const std::vector<double>& state) const override { _spatial.validate_state(state); }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override {
        return _spatial.required_state_dofs(first, last);
    }

    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override {
        _spatial.validate_local_state(first, last, state);
    }

    void set_time(double value) override;
    void set_load_factor(double value) override;
    void set_heat_source_interval(double begin, double end) override;
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) override;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const override;
    void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) override;
    void publish_material_history() noexcept override;
    void export_histories(TransientCommittedState& state) const override;
    void restore_histories(TransientCommittedState& state) override;
    RegionStateSummary summarize_region(std::size_t region) const override;
    std::vector<double> committed_creep_rates() const override;

    std::size_t sparsity_contribution_count() const noexcept override { return _spatial.sparsity_contribution_count(); }

    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept override {
        return _spatial.committed_contact_histories();
    }

    void commit_contact_state(const std::vector<double>& state) override;

    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories) override {
        _spatial.restore_contact_state(state, std::move(histories));
    }

    bool uses_augmented_contact() const noexcept override { return _spatial.uses_augmented_contact(); }

    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
        std::size_t completed_updates) override {
        return _spatial.update_augmented_contact_multipliers(state, completed_updates);
    }

    const std::vector<AxisymmetricRegionData>& kernel_data() const noexcept { return _kernel_data; }

    const std::vector<std::vector<Quad8MaterialHistory>>& histories() const noexcept { return _histories; }

    void stage_contact_history(const std::vector<double>& state) override {
        _staged_contact_solution = state;
        _spatial.prepare_contact_state(state, _staged_contact_histories);
    }

    void publish_contact_history() noexcept override {
        _spatial.publish_contact_state(_staged_contact_solution, _staged_contact_histories);
    }

  private:
    rz8::SpatialAssembly _spatial;
    std::vector<std::vector<Quad8MaterialHistory>> _histories, _staged;
    std::vector<AxisymmetricRegionData> _kernel_data;
};

class PlaneBackend final : public SpatialBackend {
  public:
    PlaneBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredPlaneQuad8Mesh& mesh,
        bool transient);

    const plane::SpatialAssembly& spatial() const noexcept { return _spatial; }

    const spatial_detail::SpatialLayout& layout() const noexcept override { return _spatial; }

    std::size_t contribution_count() const noexcept override { return _spatial.contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const override {
        return _spatial.contribution_type(index);
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    void validate_state(const std::vector<double>& state) const override { _spatial.validate_state(state); }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override {
        return _spatial.required_state_dofs(first, last);
    }

    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override {
        _spatial.validate_local_state(first, last, state);
    }

    void set_time(double value) override;
    void set_load_factor(double value) override;
    void set_heat_source_interval(double begin, double end) override;
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) override;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const override;
    void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) override;
    void publish_material_history() noexcept override;
    void export_histories(TransientCommittedState& state) const override;
    void restore_histories(TransientCommittedState& state) override;
    RegionStateSummary summarize_region(std::size_t region) const override;
    std::vector<double> committed_creep_rates() const override;

    bool combined_contact_contribution() const noexcept override { return true; }

    const std::vector<std::vector<CartesianMaterialHistory>>& histories() const noexcept { return _histories; }

  private:
    plane::SpatialAssembly _spatial;
    std::vector<std::vector<CartesianMaterialHistory>> _histories, _staged;
    double _time = 0.0;
    elements::Cpeg8Result evaluate_plane(std::size_t index,
        const std::vector<double>& local,
        bool transient,
        bool jacobian,
        bool history) const;
};

class CartesianBackend final : public SpatialBackend {
  public:
    CartesianBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredHex8Mesh& mesh,
        bool transient);
    CartesianBackend(SpatialTimeState& state,
        SpatialDefinition definition,
        const UnstructuredHex20Mesh& mesh,
        bool transient);

    const cartesian::SpatialAssembly& spatial() const noexcept { return _spatial; }

    const spatial_detail::SpatialLayout& layout() const noexcept override { return _spatial; }

    std::size_t contribution_count() const noexcept override { return _spatial.contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const override {
        return _spatial.contribution_type(index);
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    void validate_state(const std::vector<double>& state) const override { _spatial.validate_state(state); }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override {
        return _spatial.required_state_dofs(first, last);
    }

    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override {
        _spatial.validate_local_state(first, last, state);
    }

    void set_time(double value) override;
    void set_load_factor(double value) override;
    void set_heat_source_interval(double begin, double end) override;
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) override;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const override;
    void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) override;
    void publish_material_history() noexcept override;
    void export_histories(TransientCommittedState& state) const override;
    void restore_histories(TransientCommittedState& state) override;
    RegionStateSummary summarize_region(std::size_t region) const override;
    std::vector<double> committed_creep_rates() const override;

    std::size_t sparsity_contribution_count() const noexcept override { return _spatial.sparsity_contribution_count(); }

    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;

    bool jacobian_sparsity_is_state_dependent() const noexcept override {
        return _spatial.jacobian_sparsity_is_state_dependent();
    }

    std::pair<std::size_t, std::size_t> contribution_partition(std::size_t partition,
        std::size_t count) const override {
        return _spatial.contribution_partition(partition, count);
    }

    void contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const override {
        _spatial.contribution_jacobian_pattern(index, pattern);
    }

    void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const override {
        _spatial.sparsity_contribution_jacobian_pattern(index, pattern);
    }

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept override {
        return _spatial.committed_contact_histories();
    }

    void commit_contact_state(const std::vector<double>& state) override;

    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories) override {
        _spatial.restore_contact_state(state, std::move(histories));
    }

    const std::vector<std::vector<CartesianMaterialHistory>>& histories() const noexcept { return _histories; }

    void stage_contact_history(const std::vector<double>& state) override {
        _staged_contact_solution = state;
        _staged_friction_dissipation = _spatial.prepare_contact_state(state, _staged_contact_histories);
    }

    void publish_contact_history() noexcept override {
        _spatial.publish_contact_state(_staged_contact_solution, _staged_contact_histories);
    }

    void complete_contact_diagnostics(TransientConservationSummary& conservation) const override {
        conservation.friction_dissipation_increment = _staged_friction_dissipation;
    }

  private:
    cartesian::SpatialAssembly _spatial;
    std::vector<std::vector<CartesianMaterialHistory>> _histories, _staged;
};
} // namespace fuelsim
