#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
namespace fuelsim {
struct AugmentedContactUpdate;
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
class GlobalStateView final {
  public:
    explicit GlobalStateView(const std::vector<double>& dense_values);
    GlobalStateView(
        std::size_t global_size, const std::vector<std::uint32_t>& global_dofs, const std::vector<double>& values);
    std::size_t global_size() const noexcept;
    std::size_t local_size() const noexcept;
    bool contains(std::size_t global_dof) const;
    double value(std::size_t global_dof) const;

  private:
    std::size_t _global_size;
    const std::vector<double>* _dense_values;
    const std::vector<std::uint32_t>* _global_dofs;
    const std::vector<double>* _sparse_values;
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
    virtual const std::vector<FieldDescriptor>& field_layout() const noexcept = 0;
    virtual void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const = 0;
    virtual void compute_contribution_residual(
        std::size_t index, const std::vector<double>& state, std::vector<double>& residual) const = 0;
    virtual void compute_contribution_system(std::size_t index, const std::vector<double>& state,
        std::vector<double>& residual, std::vector<double>& jacobian) const = 0;
    virtual const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept = 0;
    virtual bool uses_augmented_contact() const noexcept;
    virtual AugmentedContactUpdate update_augmented_contact_multipliers(
        const std::vector<double>& state, std::size_t completed_updates);
    virtual void validate_state(const std::vector<double>& state) const;
    virtual std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    virtual void validate_local_state(std::size_t first, std::size_t last, const GlobalStateView& state) const;
    void validate_discretization() const;
    std::size_t field_index(std::size_t dof) const;
    void evaluate_contribution_residual(
        std::size_t index, const GlobalStateView& global_state, ContributionWorkspace& workspace) const;
    void evaluate_contribution_system(
        std::size_t index, const GlobalStateView& global_state, ContributionWorkspace& workspace) const;
    void assemble_residual(const std::vector<double>& state, std::vector<double>& residual) const;

  private:
    void gather_contribution_state(std::size_t index, const GlobalStateView& global_state,
        ContributionWorkspace& workspace, bool include_jacobian) const;
    std::shared_ptr<const void> _discretization_identity = std::make_shared<unsigned char>(0);
};
struct FieldNorms final {
    std::vector<double> l2, maximum_absolute;
};
struct DirectionalJacobianCheck final {
    FieldNorms residual, analytic_directional_derivative;
    FieldNorms finite_difference_directional_derivative, difference;
};
std::vector<double> constrained_residual(const NonlinearProblem& problem, const std::vector<double>& state);
FieldNorms field_norms(const NonlinearProblem& problem, const std::vector<double>& values);
DirectionalJacobianCheck check_directional_jacobian(const NonlinearProblem& problem, const std::vector<double>& state,
    const std::vector<double>& direction, double step);
} // namespace fuelsim
