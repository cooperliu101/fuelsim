#ifndef FUELSIM_STEADY_PROBLEM_HPP
#define FUELSIM_STEADY_PROBLEM_HPP

#include "fuelsim/nonlinear_problem.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace fuelsim {

struct AugmentedContactUpdate;
struct SpatialDefinition;
class UnstructuredQuad4Mesh;
class UnstructuredHex8Mesh;

namespace rz {
class ProblemAccess;
}
namespace cartesian3d {
class ProblemAccess;
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
    SteadyProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);
    ~SteadyProblem() override;

    bool is_cartesian_3d() const noexcept;
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state, std::size_t completed_updates);

    void set_load_factor(double load_factor);
    double load_factor() const noexcept;
    void set_time(double time);

    std::vector<double> initial_state() const;
    SteadyStateSnapshot capture_internal_state() const;
    void restore_internal_state(const SteadyStateSnapshot& snapshot, const std::vector<double>& state);
    void commit_internal_state(const std::vector<double>& state);
    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    const std::vector<FieldDescriptor>& field_layout() const noexcept override;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept override;
    void validate_state(const std::vector<double>& state) const override;
    std::vector<std::size_t> required_state_dofs(std::size_t contribution_begin, std::size_t contribution_end) const override;
    void validate_local_state(std::size_t contribution_begin, std::size_t contribution_end, const GlobalStateView& state) const override;
    std::size_t contribution_dof_count(std::size_t contribution_index) const override;
    void fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const override;
    void compute_contribution_residual(std::size_t contribution_index, const std::vector<double>& state, std::vector<double>& residual) const override;
    void compute_contribution_system(std::size_t contribution_index, const std::vector<double>& state, std::vector<double>& residual, std::vector<double>& jacobian) const override;

  private:
    friend class rz::ProblemAccess;
    friend class cartesian3d::ProblemAccess;

    class Implementation;
    void refresh_region_heat_sources();

    std::unique_ptr<Implementation> _implementation;
};

} // namespace fuelsim

#endif
