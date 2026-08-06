#ifndef FUELSIM_NONLINEAR_PROBLEM_HPP
#define FUELSIM_NONLINEAR_PROBLEM_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include "fuelsim/local_system.hpp"

namespace fuelsim {

struct DirichletCondition final {
    std::size_t dof;
    double value;
};

class GlobalStateView final {
  public:
    explicit GlobalStateView(const std::vector<double>& dense_values);
    GlobalStateView(std::size_t global_size,
                    const std::vector<std::uint32_t>& global_dofs,
                    const std::vector<double>& values);

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
    virtual ~NonlinearProblem() = default;

    virtual std::size_t dof_count() const noexcept = 0;
    virtual std::size_t contribution_count() const noexcept = 0;

    virtual LocalDofs
    contribution_dofs(std::size_t contribution_index) const = 0;
    virtual LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const LocalValues& state) const = 0;
    virtual LocalSystem
    linearize_contribution(std::size_t contribution_index,
                           const LocalValues& state) const = 0;

    virtual const std::vector<DirichletCondition>&
    dirichlet_conditions() const noexcept = 0;

    virtual void validate_state(const std::vector<double>& state) const;
    virtual std::vector<std::size_t> required_state_dofs(
        std::size_t contribution_begin,
        std::size_t contribution_end) const;
    virtual void validate_local_state(
        std::size_t contribution_begin, std::size_t contribution_end,
        const GlobalStateView& state) const;

    LocalValues
    contribution_state(std::size_t contribution_index,
                       const std::vector<double>& global_state) const;
    LocalValues contribution_state(std::size_t contribution_index,
                                   const GlobalStateView& global_state) const;

    void assemble_residual(const std::vector<double>& state,
                           std::vector<double>& residual) const;
    void assemble_state_independent_residual(
        std::vector<double>& residual) const;

  protected:
    virtual void
    add_state_independent_residual(std::vector<double>& residual) const = 0;
};

} // namespace fuelsim

#endif
