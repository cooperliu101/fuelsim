#pragma once
#include "fuelsim/core/steady_problem.hpp"
#include "fuelsim/core/transient_problem.hpp"
#include "fuelsim/solver/petsc_solver.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim {
struct SteadyLoadOptions final {
    std::size_t load_steps = 1;
    double cutback_factor = 0.5;
    std::size_t maximum_cutbacks_per_step = 12;
    double minimum_load_increment = 1.0e-6;
};

struct SteadyRejectedLoadStep final {
    double attempted_load_factor, load_increment;
    std::size_t cutback_index;
    SolveFailureCategory failure_category;
    std::string failure_message;
};

struct SteadyResult final {
    SolveResult solve;
    std::vector<SteadyRejectedLoadStep> rejected_steps;
    std::size_t completed_steps = 0, total_cutbacks = 0;
    bool completed = false;
    int total_nonlinear_iterations = 0, total_linear_iterations = 0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
};

SteadyResult solve_steady(
    SteadyProblem& problem, const SteadyLoadOptions& load_options, const SolverOptions& options = SolverOptions{});

struct TransientTimeOptions final {
    double end_time, initial_time_step, minimum_time_step, maximum_time_step, growth_factor, cutback_factor;
    std::size_t maximum_cutbacks_per_step;
    double load_ramp_time;
    std::size_t target_nonlinear_iterations = 0, iteration_window = 0;
    double time_error_relative_tolerance = 0.0, temperature_time_absolute_tolerance = 1.0e-3,
           displacement_time_absolute_tolerance = 1.0e-10, time_error_safety_factor = 0.9,
           strain_history_time_absolute_tolerance = 1.0e-10, stress_history_time_absolute_tolerance = 1.0;
    bool include_thermal_time_term = true;
};

struct TransientFieldTimeError final {
    std::string name;
    double value = 0.0;
};

struct TransientTimeErrorEstimate final {
    std::vector<TransientFieldTimeError> nodal_fields;
    double elastic_strain = 0.0, plastic_strain = 0.0, creep_strain = 0.0, equivalent_plastic_strain = 0.0,
           equivalent_creep_strain = 0.0, stress = 0.0, contact_friction = 0.0, contact_normal_multiplier = 0.0,
           maximum = 0.0;
};
enum class TransientTerminationReason {
    not_started,
    completed,
    maximum_cutbacks,
    minimum_time_step,
};

struct TransientRejectedStep final {
    double attempted_end_time, time_step;
    std::size_t cutback_index;
    int nonlinear_iterations, linear_iterations, convergence_reason;
    double residual_norm;
    SolveFailureCategory failure_category;
    std::string failure_message;
    double time_error_estimate = 0.0;
    TransientTimeErrorEstimate time_error_components;
};

struct TransientAcceptedStep final {
    double time, time_step, next_time_step, load_factor;
    std::size_t cutbacks;
    int nonlinear_iterations, linear_iterations;
    double time_error_estimate = 0.0;
    TransientTimeErrorEstimate time_error_components;
    TransientConservationSummary conservation;
};

struct TransientResult final {
    SolveResult last_attempt;
    std::vector<double> committed_state;
    std::vector<TransientAcceptedStep> accepted_steps;
    std::vector<TransientRejectedStep> rejected_steps;
    bool completed = false;
    std::size_t total_cutbacks = 0, time_error_rejections = 0;
    int total_nonlinear_iterations = 0, total_linear_iterations = 0;
    double committed_time = 0.0, next_time_step = 0.0, total_seconds = 0.0;
    SolveTiming aggregate_timing;
    TransientTerminationReason termination_reason = TransientTerminationReason::not_started;
};

const char* transient_termination_reason_name(TransientTerminationReason reason) noexcept;

class TransientStepObserver {
  public:
    virtual ~TransientStepObserver() = default;
    virtual void accepted_step(const TransientProblem& problem, const TransientAcceptedStep& step) = 0;
};

TransientResult solve_transient(TransientProblem& problem, const TransientTimeOptions& time_options,
    const SolverOptions& solver_options = SolverOptions{}, TransientStepObserver* observer = nullptr);
} // namespace fuelsim
