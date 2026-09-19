#pragma once
#include "contact_types.hpp"
#include "core/transient_problem.hpp"
#include "spatial_layout.hpp"
#include <exception>
#include <functional>

namespace fuelsim {
struct TransientCommittedState;

struct SpatialTimeState final {
    TransientConservationSummary last_conservation_summary;
    std::vector<double> committed_solution, previous_committed_solution, committed_raw_residual,
        committed_external_load_residual;
    double committed_time = 0.0, committed_load_factor = 0.0, active_time_step = 0.0, active_end_time = 0.0,
           active_load_factor = 0.0, previous_committed_time = 0.0;
    bool time_step_active = false, include_thermal_time_term = true, track_previous_committed_solution = false;
};

class SpatialBackend {
  public:
    explicit SpatialBackend(SpatialTimeState& state) : _state(state) {}

    virtual ~SpatialBackend() = default;
    SpatialBackend(const SpatialBackend&) = delete;
    SpatialBackend& operator=(const SpatialBackend&) = delete;
    SpatialBackend(SpatialBackend&&) = delete;
    SpatialBackend& operator=(SpatialBackend&&) = delete;
    virtual const spatial_detail::SpatialLayout& layout() const noexcept = 0;
    virtual std::size_t contribution_count() const noexcept = 0;

    virtual std::size_t sparsity_contribution_count() const noexcept { return contribution_count(); }

    virtual bool jacobian_sparsity_is_state_dependent() const noexcept { return false; }

    virtual bool contribution_metadata_is_fixed() const noexcept { return layout().definition().contacts.empty(); }

    virtual std::pair<std::size_t, std::size_t> contribution_partition(std::size_t partition,
        std::size_t partition_count) const;
    virtual void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const = 0;
    virtual void contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const;

    virtual void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
        contribution_dofs(index, dofs);
    }

    virtual void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const;
    virtual SpatialContributionType contribution_type(std::size_t index) const = 0;

    virtual bool combined_contact_contribution() const noexcept { return false; }

    virtual void validate_state(const std::vector<double>& state) const = 0;
    virtual std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const = 0;
    virtual void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const = 0;
    virtual void set_load_factor(double value) = 0;
    virtual void set_time(double value) = 0;
    virtual void set_region_strain_formulation(std::size_t region, StrainFormulation formulation) = 0;
    virtual void set_heat_source_interval(double begin, double end) = 0;
    virtual void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool transient) const = 0;
    virtual const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept;
    virtual void commit_contact_state(const std::vector<double>& state);
    virtual void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories);

    virtual bool uses_augmented_contact() const noexcept { return false; }

    virtual AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
        std::size_t completed_updates);
    virtual void prepare_time_step(const NonlinearProblem& problem,
        const std::vector<double>& converged_solution,
        std::size_t first_contribution,
        std::size_t last_contribution,
        const std::function<void(std::vector<double>&)>& sum_partitions,
        TransientConservationSummary& conservation,
        std::vector<double>& raw_residual,
        std::vector<double>& external_load_residual,
        std::exception_ptr partition_failure) = 0;
    virtual void publish_material_history() noexcept = 0;
    void prepare_contacts(const std::vector<double>& state,
        std::exception_ptr& failure,
        const std::function<void(std::vector<double>&)>& sum_partitions);

    virtual void publish_contact_history() noexcept {}

    virtual void complete_contact_diagnostics(TransientConservationSummary&) const {}

    void capture_active_contact_state();
    void restore_active_contact_state();

    void clear_active_contact_state() noexcept { _active_contact_histories.clear(); }

    virtual void export_histories(TransientCommittedState& state) const = 0;
    virtual void restore_histories(TransientCommittedState& state) = 0;
    virtual RegionStateSummary summarize_region(std::size_t region) const = 0;
    virtual std::vector<double> committed_creep_rates() const = 0;

  protected:
    virtual void stage_contact_history(const std::vector<double>&) {}

    // The owning storage is immovable and outlives this borrowed time state.
    const SpatialTimeState& _state;
    std::vector<std::vector<ContactPointHistory>> _staged_contact_histories;
    std::vector<double> _staged_contact_solution;
    double _staged_friction_dissipation = 0.0;

  private:
    std::vector<std::vector<ContactPointHistory>> _active_contact_histories;
};

} // namespace fuelsim
