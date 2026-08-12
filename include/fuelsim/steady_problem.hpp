#ifndef FUELSIM_STEADY_PROBLEM_HPP
#define FUELSIM_STEADY_PROBLEM_HPP

#include "fuelsim/dof_map.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz_thermoelastic.hpp"
#include "fuelsim/spatial_definition.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {

namespace rz {
class SpatialAssembly;
}

class SteadyStateSnapshot final {
  public:
    SteadyStateSnapshot();
    ~SteadyStateSnapshot();
    SteadyStateSnapshot(const SteadyStateSnapshot& other);
    SteadyStateSnapshot& operator=(const SteadyStateSnapshot& other);
    SteadyStateSnapshot(SteadyStateSnapshot&& other) noexcept;
    SteadyStateSnapshot& operator=(SteadyStateSnapshot&& other) noexcept;

    bool empty() const noexcept;

  private:
    std::shared_ptr<const void> snapshot_owner() const noexcept;

    struct Storage;
    explicit SteadyStateSnapshot(std::shared_ptr<const Storage> storage);

    std::shared_ptr<const Storage> _storage;
    friend class SteadyProblem;
};

class SteadyProblem final : public NonlinearProblem {
  public:
    SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh);
    ~SteadyProblem() override;

    const SpatialDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;

    std::size_t region_count() const noexcept;
    std::size_t region_index(const std::string& name) const;
    const RegionDefinition& region(std::size_t region_index) const;
    const RegionMesh& region_mesh(std::size_t region_index) const;
    const Quad4RzThermoelasticKernel& region_kernel(std::size_t region_index) const;
    std::size_t region_node_offset(std::size_t region_index) const;
    std::size_t region_element_count(std::size_t region_index) const;
    std::size_t region_element_offset(std::size_t region_index) const;
    std::size_t volume_contribution_count() const noexcept;
    SpatialContributionType contribution_type(std::size_t contribution_index) const;
    const Quad4RzGeometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;

    std::size_t contact_count() const noexcept;
    const ContactDefinition& contact(std::size_t contact_index) const;
    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept;
    void commit_contact_state(const std::vector<double>& state);
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                std::size_t completed_updates);
    void restore_contact_state(const std::vector<double>& state,
                               std::vector<std::vector<ContactPointHistory>> histories);

    void set_load_factor(double load_factor);
    double load_factor() const noexcept;
    void set_time(double time);

    std::vector<double> initial_state() const;
    std::vector<ContactNodeSummary> summarize_contact_nodes(std::size_t contact_index,
                                                            const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;
    SteadyStateSnapshot capture_internal_state() const;
    void restore_internal_state(const SteadyStateSnapshot& snapshot, const std::vector<double>& state);
    void commit_internal_state(const std::vector<double>& state);
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

    LocalDofs contribution_dofs(std::size_t contribution_index) const;
    LocalValues contribution_state(std::size_t contribution_index, const std::vector<double>& global_state) const;
    LocalValues contribution_state(std::size_t contribution_index, const GlobalStateView& global_state) const;
    LocalResidual contribution_residual(std::size_t contribution_index, const LocalValues& state) const;
    LocalSystem linearize_contribution(std::size_t contribution_index, const LocalValues& state) const;

  private:
    void refresh_region_heat_sources();

    std::unique_ptr<rz::SpatialAssembly> _spatial;
    std::vector<Quad4RzThermoelasticKernel> _region_kernels;
};

} // namespace fuelsim

#endif
