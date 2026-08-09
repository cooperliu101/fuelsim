#ifndef FUELSIM_TRANSIENT_TIME_CONTROL_HPP
#define FUELSIM_TRANSIENT_TIME_CONTROL_HPP

#include "fuelsim/problem_solver.hpp"

namespace fuelsim::time_control {

TransientConservationSummary combine_half_step_conservation(
    const TransientConservationSummary& first,
    const TransientConservationSummary& second);

TransientTimeErrorEstimate step_doubling_error(
    const TransientCommittedState& full_step,
    const TransientCommittedState& two_half_steps,
    const TransientTimeOptions& options);

double step_factor(const TransientTimeOptions& options, double error);

} // namespace fuelsim::time_control

#endif
