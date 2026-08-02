#ifndef FUELSIM_NONLINEAR_PROBLEM_HPP
#define FUELSIM_NONLINEAR_PROBLEM_HPP

#include <cstddef>
#include <vector>

#include "fuelsim/local_system.hpp"

namespace fuelsim {

struct DirichletCondition final {
    std::size_t dof;
    double value;
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

    LocalValues
    contribution_state(std::size_t contribution_index,
                       const std::vector<double>& global_state) const;

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
