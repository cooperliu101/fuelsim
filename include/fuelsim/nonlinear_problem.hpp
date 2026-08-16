#pragma once
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
struct AugmentedContactUpdate;
class SteadyProblem;
class TransientProblem;

class ProblemStateSnapshot final {
  public:
    ProblemStateSnapshot() = default;

    bool empty() const noexcept { return _state == nullptr; }

  private:
    ProblemStateSnapshot(const std::shared_ptr<const void>& owner, const std::shared_ptr<const void>& state)
        : _owner(owner), _state(state) {}

    std::shared_ptr<const void> _owner, _state;
    friend class SteadyProblem;
    friend class TransientProblem;
};
enum class FieldCategory {
    thermal,
    mechanical,
};

struct FieldDescriptor final {
    std::string name;
    std::size_t begin, end;
    FieldCategory category;
};

struct ContributionWorkspace final {
    void reserve(std::size_t maximum_dof_count);
    void resize(std::size_t dof_count, bool include_jacobian);
    std::vector<std::size_t> dofs;
    std::vector<double> state, residual;
    std::vector<double> jacobian;
};

struct DirichletCondition final {
    std::size_t dof;
    double value;
};

class NonlinearProblem {
  public:
    NonlinearProblem() = default;
    virtual ~NonlinearProblem() = default;
    NonlinearProblem(const NonlinearProblem&) = delete;
    NonlinearProblem& operator=(const NonlinearProblem&) = delete;
    NonlinearProblem(NonlinearProblem&&) = delete;
    NonlinearProblem& operator=(NonlinearProblem&&) = delete;

    // The identity keeps field ranges, contribution mappings, constraints, and shadow DOFs immutable while load,
    // time, trial state, and residual values may change.
    std::shared_ptr<const void> discretization_identity() const noexcept { return _discretization_identity; }

    virtual std::size_t dof_count() const noexcept = 0;
    virtual std::size_t contribution_count() const noexcept = 0;
    virtual std::size_t sparsity_contribution_count() const noexcept;
    virtual std::pair<std::size_t, std::size_t> contribution_partition(
        std::size_t partition, std::size_t partition_count) const;
    virtual const std::vector<FieldDescriptor>& field_layout() const noexcept = 0;
    virtual void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const = 0;
    virtual void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    virtual void compute_contribution(std::size_t index, const std::vector<double>& state,
        std::vector<double>& residual, std::vector<double>* jacobian) const = 0;
    virtual const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept = 0;
    virtual bool uses_augmented_contact() const noexcept;
    virtual AugmentedContactUpdate update_augmented_contact_multipliers(
        const std::vector<double>& state, std::size_t completed_updates);
    virtual void validate_state(const std::vector<double>& state) const;
    virtual std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    virtual void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void validate_discretization() const;
    std::size_t field_index(std::size_t dof) const;
    void evaluate_contribution(std::size_t index, const std::vector<double>& global_state,
        ContributionWorkspace& workspace, bool linearize) const;

  private:
    std::shared_ptr<const void> _discretization_identity = std::make_shared<unsigned char>(0);
};
} // namespace fuelsim
