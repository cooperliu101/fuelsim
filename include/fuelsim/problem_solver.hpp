#ifndef FUELSIM_PROBLEM_SOLVER_HPP
#define FUELSIM_PROBLEM_SOLVER_HPP

#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

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
    double attempted_load_factor;
    double load_increment;
    std::size_t cutback_index;
    SolveFailureCategory failure_category;
    std::string failure_message;
};

struct SteadyResult final {
    SolveResult solve;
    std::vector<SteadyRejectedLoadStep> rejected_steps;
    std::size_t completed_steps = 0;
    std::size_t total_cutbacks = 0;
    bool completed = false;
    int total_nonlinear_iterations = 0;
    int total_linear_iterations = 0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
};

SteadyResult solve_steady(SteadyProblem& problem, const SteadyLoadOptions& load_options, const SolverOptions& options = SolverOptions{});

struct TransientTimeOptions final {
    double end_time;
    double initial_time_step;
    double minimum_time_step;
    double maximum_time_step;
    double growth_factor;
    double cutback_factor;
    std::size_t maximum_cutbacks_per_step;
    double load_ramp_time;
    std::size_t target_nonlinear_iterations = 0;
    std::size_t iteration_window = 0;
    double time_error_relative_tolerance = 0.0;
    double temperature_time_absolute_tolerance = 1.0e-3;
    double displacement_time_absolute_tolerance = 1.0e-10;
    double time_error_safety_factor = 0.9;
    double strain_history_time_absolute_tolerance = 1.0e-10;
    double stress_history_time_absolute_tolerance = 1.0;
};

struct TransientFieldTimeError final {
    std::string name;
    double value = 0.0;
};

struct TransientTimeErrorEstimate final {
    std::vector<TransientFieldTimeError> nodal_fields;
    double elastic_strain = 0.0;
    double plastic_strain = 0.0;
    double creep_strain = 0.0;
    double equivalent_plastic_strain = 0.0;
    double equivalent_creep_strain = 0.0;
    double stress = 0.0;
    double contact_friction = 0.0;
    double contact_normal_multiplier = 0.0;
    double maximum = 0.0;
};

enum class TransientTerminationReason {
    not_started,
    completed,
    maximum_cutbacks,
    minimum_time_step,
};

struct TransientRejectedStep final {
    double attempted_end_time;
    double time_step;
    std::size_t cutback_index;
    int nonlinear_iterations;
    int linear_iterations;
    int convergence_reason;
    double residual_norm;
    SolveFailureCategory failure_category;
    std::string failure_message;
    double time_error_estimate = 0.0;
    TransientTimeErrorEstimate time_error_components;
};

struct TransientAcceptedStep final {
    double time;
    double time_step;
    double next_time_step;
    double load_factor;
    std::size_t cutbacks;
    int nonlinear_iterations;
    int linear_iterations;
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
    std::size_t total_cutbacks = 0;
    std::size_t time_error_rejections = 0;
    int total_nonlinear_iterations = 0;
    int total_linear_iterations = 0;
    double committed_time = 0.0;
    double next_time_step = 0.0;
    double total_seconds = 0.0;
    SolveTiming aggregate_timing;
    TransientTerminationReason termination_reason = TransientTerminationReason::not_started;
};

const char* transient_termination_reason_name(TransientTerminationReason reason) noexcept;

class TransientStepObserver {
  public:
    virtual ~TransientStepObserver() = default;
    virtual void accepted_step(const TransientProblem& problem, const TransientAcceptedStep& step) = 0;
};

TransientResult solve_transient(TransientProblem& problem, const TransientTimeOptions& time_options, const SolverOptions& solver_options = SolverOptions{}, TransientStepObserver* observer = nullptr);

} // namespace fuelsim

#endif
