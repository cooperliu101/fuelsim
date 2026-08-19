#pragma once
#include "fuelsim/solver/petsc_solver.hpp"
#include <algorithm>
#include <chrono>

namespace fuelsim::solver_detail {
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
