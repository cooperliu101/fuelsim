#ifndef FUELSIM_DIAGNOSTICS_HPP
#define FUELSIM_DIAGNOSTICS_HPP
#include "fuelsim/nonlinear_problem.hpp"
#include <vector>
namespace fuelsim {
struct FieldNorms final {
    std::vector<double> l2;
    std::vector<double> maximum_absolute;
};
struct DirectionalJacobianCheck final {
    FieldNorms residual;
    FieldNorms analytic_directional_derivative;
    FieldNorms finite_difference_directional_derivative;
    FieldNorms difference;
};
std::vector<double> constrained_residual(const NonlinearProblem& problem, const std::vector<double>& state);
FieldNorms field_norms(const NonlinearProblem& problem, const std::vector<double>& values);
DirectionalJacobianCheck check_directional_jacobian(const NonlinearProblem& problem, const std::vector<double>& state,
    const std::vector<double>& direction, double step);
} // namespace fuelsim
#endif
