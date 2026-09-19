#pragma once
#include "solver/petsc_solver.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace fuelsim {
class TransientProblem;
}

namespace fuelsim::solver_detail {
void commit_time_step(TransientProblem& problem, const SolveResult& result);

inline void add_quadratic_time_correction(std::vector<double>& prediction,
    const std::vector<double>& current,
    const std::vector<double>& previous,
    const std::vector<double>& older,
    double current_time,
    double previous_time,
    double older_time,
    double target_time) {
    const double recent_interval = current_time - previous_time;
    const double older_interval = previous_time - older_time;
    const double forward_interval = target_time - current_time;
    if (prediction.size() != current.size() || previous.size() != current.size() || older.size() != current.size()
        || !(recent_interval > 0.0) || !(older_interval > 0.0) || !(forward_interval > 0.0))
        return;
    // Do not extrapolate curvature across abrupt changes in step size.
    const double previous_ratio = recent_interval / older_interval;
    const double forward_ratio = forward_interval / recent_interval;
    if (previous_ratio < 0.5 || previous_ratio > 2.0 || forward_ratio < 0.5 || forward_ratio > 2.0)
        return;
    const double factor = forward_interval * (forward_interval + recent_interval) / (recent_interval + older_interval);
    if (!std::isfinite(factor))
        return;
    for (std::size_t dof = 0; dof < current.size(); ++dof)
        prediction[dof] +=
            factor * ((current[dof] - previous[dof]) / recent_interval - (previous[dof] - older[dof]) / older_interval);
}

using Clock = std::chrono::steady_clock;

inline double seconds_since(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

inline void accumulate_timing(SolveTiming& total, const SolveTiming& step) {
    total.setup_seconds += step.setup_seconds;
    total.nonlinear_solve_seconds += step.nonlinear_solve_seconds;
    total.residual_callback_seconds += step.residual_callback_seconds;
    total.jacobian_callback_seconds += step.jacobian_callback_seconds;
    total.minimum_residual_assembly_seconds += step.minimum_residual_assembly_seconds;
    total.maximum_residual_assembly_seconds += step.maximum_residual_assembly_seconds;
    total.minimum_jacobian_assembly_seconds += step.minimum_jacobian_assembly_seconds;
    total.maximum_jacobian_assembly_seconds += step.maximum_jacobian_assembly_seconds;
    total.local_residual_assembly_seconds += step.local_residual_assembly_seconds;
    total.local_jacobian_assembly_seconds += step.local_jacobian_assembly_seconds;
    total.total_seconds += step.total_seconds;
    total.initial_resident_bytes = std::max(total.initial_resident_bytes, step.initial_resident_bytes);
    total.setup_resident_bytes = std::max(total.setup_resident_bytes, step.setup_resident_bytes);
    total.solve_resident_bytes = std::max(total.solve_resident_bytes, step.solve_resident_bytes);
    total.final_resident_bytes = std::max(total.final_resident_bytes, step.final_resident_bytes);
    total.minimum_peak_resident_bytes = std::max(total.minimum_peak_resident_bytes, step.minimum_peak_resident_bytes);
    total.maximum_peak_resident_bytes = std::max(total.maximum_peak_resident_bytes, step.maximum_peak_resident_bytes);
    total.total_peak_resident_bytes = std::max(total.total_peak_resident_bytes, step.total_peak_resident_bytes);
    total.residual_evaluations += step.residual_evaluations;
    total.jacobian_evaluations += step.jacobian_evaluations;
    total.workspace_setups += step.workspace_setups;
    total.solve_calls += step.solve_calls;
}
} // namespace fuelsim::solver_detail
