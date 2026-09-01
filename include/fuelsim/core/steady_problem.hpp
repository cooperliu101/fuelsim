#pragma once
#include "fuelsim/core/nonlinear_problem.hpp"
#include <cstddef>
#include <memory>
#include <vector>

namespace fuelsim {
struct AugmentedContactUpdate;
struct SpatialDefinition;
class UnstructuredQuad4Mesh;
class UnstructuredHex8Mesh;
class UnstructuredHex20Mesh;
class SpatialProblemStorage;
class BackendAccess;

class SteadyProblem final : public NonlinearProblem {
  public:
    SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh);
    SteadyProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);
    SteadyProblem(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh);
    ~SteadyProblem() override;
    bool uses_augmented_contact() const noexcept override;
    AugmentedContactUpdate update_augmented_contact_multipliers(
        const std::vector<double>& state, std::size_t completed_updates) override;
    void set_load_factor(double load_factor);
    double load_factor() const noexcept;
    void set_time(double time);
    std::vector<double> initial_state() const;
    ProblemStateSnapshot capture_internal_state() const;
    void restore_internal_state(const ProblemStateSnapshot& snapshot, const std::vector<double>& state);
    void commit_internal_state(const std::vector<double>& state);
    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    std::size_t sparsity_contribution_count() const noexcept override;
    bool jacobian_sparsity_is_state_dependent() const noexcept override;
    std::pair<std::size_t, std::size_t> contribution_partition(
        std::size_t partition, std::size_t partition_count) const override;
    const std::vector<FieldDescriptor>& field_layout() const noexcept override;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept override;
    void validate_state(const std::vector<double>& state) const override;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const override;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const override;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;
    void contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const override;
    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const override;
    void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const override;
    void compute_contribution(std::size_t index, const std::vector<double>& state, std::vector<double>& residual,
        std::vector<double>* jacobian) const override;

  private:
    friend class BackendAccess;
    std::unique_ptr<SpatialProblemStorage> _impl;
};
} // namespace fuelsim
