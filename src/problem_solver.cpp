#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

using SteadyClock = std::chrono::steady_clock;

double seconds_since(const SteadyClock::time_point& start) {
    return std::chrono::duration<double>(SteadyClock::now() - start).count();
}

void accumulate_timing(SolveTiming& total, const SolveTiming& step) {
    total.setup_seconds += step.setup_seconds;
    total.nonlinear_solve_seconds += step.nonlinear_solve_seconds;
    total.residual_callback_seconds += step.residual_callback_seconds;
    total.jacobian_callback_seconds += step.jacobian_callback_seconds;
    total.total_seconds += step.total_seconds;
    total.residual_evaluations += step.residual_evaluations;
    total.jacobian_evaluations += step.jacobian_evaluations;
    total.workspace_setups += step.workspace_setups;
    total.solve_calls += step.solve_calls;
}

class TimeStepTransaction final {
  public:
    TimeStepTransaction(TransientProblem& problem,
                        const TransientStepInput& input)
        : _problem(problem), _committed(false) {
        _problem.begin_time_step(input);
    }

    ~TimeStepTransaction() {
        if (!_committed)
            _problem.rollback_time_step();
    }

    TimeStepTransaction(const TimeStepTransaction&) = delete;
    TimeStepTransaction& operator=(const TimeStepTransaction&) = delete;

    void commit(const std::vector<double>& solution) {
        _problem.commit_time_step(solution);
        _committed = true;
    }

  private:
    TransientProblem& _problem;
    bool _committed;
};

bool reaches_end(double time, double end_time) {
    const double scale = std::max({1.0, std::abs(time), std::abs(end_time)});
    return end_time - time <=
           16.0 * std::numeric_limits<double>::epsilon() * scale;
}

void validate_time_options(const TransientProblem& problem,
                           const TransientTimeOptions& options) {
    if (problem.time_step_active())
        throw std::logic_error(
            "solve_transient cannot start with an active time step");
    const double time_scale = std::max(
        {1.0, std::abs(problem.committed_time()), std::abs(options.end_time)});
    const double time_tolerance =
        16.0 * std::numeric_limits<double>::epsilon() * time_scale;
    if (!std::isfinite(options.end_time) ||
        options.end_time < problem.committed_time() - time_tolerance)
        throw std::invalid_argument(
            "solve_transient end time must not precede committed time");
    if (!std::isfinite(options.initial_time_step) ||
        !std::isfinite(options.minimum_time_step) ||
        !std::isfinite(options.maximum_time_step) ||
        !(options.minimum_time_step > 0.0) ||
        options.initial_time_step < options.minimum_time_step ||
        options.maximum_time_step < options.initial_time_step)
        throw std::invalid_argument("solve_transient requires 0 < minimum <= "
                                    "initial <= maximum time step");
    if (!std::isfinite(options.growth_factor) || options.growth_factor < 1.0)
        throw std::invalid_argument(
            "solve_transient growth factor must be at least one");
    if (!std::isfinite(options.cutback_factor) ||
        !(options.cutback_factor > 0.0 && options.cutback_factor < 1.0))
        throw std::invalid_argument(
            "solve_transient cutback factor must lie between zero and one");
    if (!std::isfinite(options.load_ramp_time) || options.load_ramp_time < 0.0)
        throw std::invalid_argument(
            "solve_transient ramp time must be finite and nonnegative");
    if (options.target_nonlinear_iterations == 0 &&
        options.iteration_window != 0)
        throw std::invalid_argument(
            "solve_transient iteration window requires a target");
    if (options.target_nonlinear_iterations > 0 &&
        options.iteration_window >= options.target_nonlinear_iterations)
        throw std::invalid_argument(
            "solve_transient iteration window must be smaller than target");
    if (!std::isfinite(options.time_error_relative_tolerance) ||
        options.time_error_relative_tolerance < 0.0 ||
        !std::isfinite(options.temperature_time_absolute_tolerance) ||
        !(options.temperature_time_absolute_tolerance > 0.0) ||
        !std::isfinite(options.displacement_time_absolute_tolerance) ||
        !(options.displacement_time_absolute_tolerance > 0.0) ||
        !std::isfinite(
            options.strain_history_time_absolute_tolerance) ||
        !(options.strain_history_time_absolute_tolerance > 0.0) ||
        !std::isfinite(
            options.stress_history_time_absolute_tolerance) ||
        !(options.stress_history_time_absolute_tolerance > 0.0) ||
        !std::isfinite(options.time_error_safety_factor) ||
        !(options.time_error_safety_factor > 0.0 &&
          options.time_error_safety_factor < 1.0))
        throw std::invalid_argument(
            "solve_transient time-error tolerances must be finite and "
            "nonnegative/positive, and safety factor must lie in (0, 1)");
}

double load_factor_at_time(const TransientTimeOptions& options, double time) {
    if (options.load_ramp_time == 0.0)
        return 1.0;
    return std::min(time / options.load_ramp_time, 1.0);
}

std::vector<double>
initial_guess_with_dirichlet_values(const NonlinearProblem& problem,
                                    const std::vector<double>& state) {
    if (state.size() != problem.dof_count())
        throw std::invalid_argument(
            "Dirichlet initial-guess state size does not match problem");
    std::vector<double> result = state;
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = condition.value;
    return result;
}

double accepted_next_time_step(const TransientTimeOptions& options,
                               double actual_time_step,
                               double controller_time_step,
                               bool event_truncated, std::size_t cutbacks,
                               int nonlinear_iterations) {
    const double base = event_truncated && cutbacks == 0 ? controller_time_step
                                                         : actual_time_step;
    if (options.target_nonlinear_iterations == 0) {
        if (cutbacks > 0)
            return std::clamp(base, options.minimum_time_step,
                              options.maximum_time_step);
        return std::min(options.maximum_time_step,
                        base * options.growth_factor);
    }
    const std::size_t iterations =
        nonlinear_iterations < 0
            ? 0
            : static_cast<std::size_t>(nonlinear_iterations);
    const std::size_t lower =
        options.target_nonlinear_iterations - options.iteration_window;
    const std::size_t upper =
        options.target_nonlinear_iterations >
                std::numeric_limits<std::size_t>::max() -
                    options.iteration_window
            ? std::numeric_limits<std::size_t>::max()
            : options.target_nonlinear_iterations + options.iteration_window;
    if (iterations < lower && cutbacks == 0)
        return std::min(options.maximum_time_step,
                        base * options.growth_factor);
    if (iterations > upper)
        return std::max(options.minimum_time_step,
                        base * options.cutback_factor);
    return std::clamp(base, options.minimum_time_step,
                      options.maximum_time_step);
}

void combine_attempt(SolveResult& aggregate, const SolveResult& addition) {
    const int nonlinear_iterations =
        aggregate.nonlinear_iterations + addition.nonlinear_iterations;
    const int linear_iterations =
        aggregate.linear_iterations + addition.linear_iterations;
    const std::size_t nonlinear_attempts =
        aggregate.nonlinear_attempts + addition.nonlinear_attempts;
    SolveTiming timing = aggregate.timing;
    accumulate_timing(timing, addition.timing);
    const bool used_backtracking = aggregate.used_backtracking_fallback ||
                                   addition.used_backtracking_fallback;
    const SolveFailureCategory basic_failure =
        aggregate.basic_failure_category != SolveFailureCategory::none
            ? aggregate.basic_failure_category
            : addition.basic_failure_category;
    const std::string basic_message =
        !aggregate.basic_failure_message.empty()
            ? aggregate.basic_failure_message
            : addition.basic_failure_message;
    aggregate = addition;
    aggregate.nonlinear_iterations = nonlinear_iterations;
    aggregate.linear_iterations = linear_iterations;
    aggregate.nonlinear_attempts = nonlinear_attempts;
    aggregate.timing = timing;
    aggregate.used_backtracking_fallback = used_backtracking;
    aggregate.basic_failure_category = basic_failure;
    aggregate.basic_failure_message = basic_message;
}

TransientConservationSummary combine_half_step_conservation(
    const TransientConservationSummary& first,
    const TransientConservationSummary& second) {
    TransientConservationSummary result;
    const auto average = [](double left, double right) {
        return 0.5 * (left + right);
    };
    result.generated_heat_rate =
        average(first.generated_heat_rate, second.generated_heat_rate);
    result.stored_heat_rate =
        average(first.stored_heat_rate, second.stored_heat_rate);
    result.convection_heat_rate =
        average(first.convection_heat_rate, second.convection_heat_rate);
    result.interface_heat_imbalance = average(
        first.interface_heat_imbalance, second.interface_heat_imbalance);
    result.dirichlet_heat_input_rate = average(
        first.dirichlet_heat_input_rate, second.dirichlet_heat_input_rate);
    result.global_thermal_balance =
        result.stored_heat_rate + result.convection_heat_rate +
        result.interface_heat_imbalance - result.generated_heat_rate -
        result.dirichlet_heat_input_rate;
    const double thermal_scale =
        std::abs(result.generated_heat_rate) +
        std::abs(result.stored_heat_rate) +
        std::abs(result.convection_heat_rate) +
        std::abs(result.interface_heat_imbalance) +
        std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0
            ? std::abs(result.global_thermal_balance) / thermal_scale
            : 0.0;
    result.unconstrained_thermal_residual_l2 = std::max(
        first.unconstrained_thermal_residual_l2,
        second.unconstrained_thermal_residual_l2);

    result.internal_mechanical_work_increment =
        first.internal_mechanical_work_increment +
        second.internal_mechanical_work_increment;
    result.pressure_traction_work_increment =
        first.pressure_traction_work_increment +
        second.pressure_traction_work_increment;
    result.dirichlet_reaction_work_increment =
        first.dirichlet_reaction_work_increment +
        second.dirichlet_reaction_work_increment;
    result.contact_work_increment = first.contact_work_increment +
                                    second.contact_work_increment;
    result.mechanical_work_balance =
        result.internal_mechanical_work_increment +
        result.contact_work_increment -
        result.pressure_traction_work_increment -
        result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) +
        std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment) +
        std::abs(result.contact_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0
            ? std::abs(result.mechanical_work_balance) / mechanical_scale
            : 0.0;
    result.unconstrained_mechanical_residual_l2 = std::max(
        first.unconstrained_mechanical_residual_l2,
        second.unconstrained_mechanical_residual_l2);
    result.elastic_energy_change = first.elastic_energy_change +
                                   second.elastic_energy_change;
    result.plastic_dissipation_increment =
        first.plastic_dissipation_increment +
        second.plastic_dissipation_increment;
    result.creep_dissipation_increment =
        first.creep_dissipation_increment +
        second.creep_dissipation_increment;
    return result;
}

struct ErrorAccumulator final {
    double difference_squared = 0.0;
    double solution_squared = 0.0;
    std::size_t count = 0;
};

void accumulate_error(ErrorAccumulator& accumulator, double full_step,
                      double two_half_steps) {
    const double difference = two_half_steps - full_step;
    accumulator.difference_squared += difference * difference;
    accumulator.solution_squared += two_half_steps * two_half_steps;
    ++accumulator.count;
}

double normalized_error(const ErrorAccumulator& accumulator,
                        double absolute_tolerance,
                        double relative_tolerance) {
    if (accumulator.count == 0)
        return 0.0;
    const double denominator =
        absolute_tolerance *
            std::sqrt(static_cast<double>(accumulator.count)) +
        relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}

TransientTimeErrorEstimate
step_doubling_error(const TransientCommittedState& full_step,
                    const TransientCommittedState& two_half_steps,
                    const TransientTimeOptions& options) {
    if (full_step.solution.size() != two_half_steps.solution.size() ||
        full_step.solution.size() % 3 != 0)
        throw std::logic_error(
            "step-doubling states must share the [T, ur, uz] layout");
    if (full_step.material_histories.size() !=
            two_half_steps.material_histories.size() ||
        full_step.material_stresses.size() !=
            two_half_steps.material_stresses.size())
        throw std::logic_error(
            "step-doubling material-state region layouts differ");

    const std::size_t node_count = full_step.solution.size() / 3;
    std::array<ErrorAccumulator, 3> nodal{};
    for (std::size_t field = 0; field < 3; ++field) {
        const std::size_t begin = field * node_count;
        const std::size_t end = begin + node_count;
        for (std::size_t dof = begin; dof < end; ++dof)
            accumulate_error(nodal[field], full_step.solution[dof],
                             two_half_steps.solution[dof]);
    }

    ErrorAccumulator elastic;
    ErrorAccumulator plastic;
    ErrorAccumulator creep;
    ErrorAccumulator equivalent_plastic;
    ErrorAccumulator equivalent_creep;
    ErrorAccumulator stress;
    ErrorAccumulator contact_friction;
    bool contact_state_mismatch = false;
    for (std::size_t region = 0;
         region < full_step.material_histories.size(); ++region) {
        const auto& full_history = full_step.material_histories[region];
        const auto& half_history = two_half_steps.material_histories[region];
        const auto& full_stress = full_step.material_stresses[region];
        const auto& half_stress = two_half_steps.material_stresses[region];
        if (full_history.size() != half_history.size() ||
            full_stress.size() != half_stress.size() ||
            full_history.size() != full_stress.size())
            throw std::logic_error(
                "step-doubling material-state element layouts differ");
        for (std::size_t element = 0; element < full_history.size();
             ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const MaterialPointState& full_point =
                    full_history[element][q];
                const MaterialPointState& half_point =
                    half_history[element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_error(elastic,
                                     full_point.elastic_strain[component],
                                     half_point.elastic_strain[component]);
                    accumulate_error(plastic,
                                     full_point.plastic_strain[component],
                                     half_point.plastic_strain[component]);
                    accumulate_error(creep,
                                     full_point.creep_strain[component],
                                     half_point.creep_strain[component]);
                }
                accumulate_error(equivalent_plastic,
                                 full_point.equivalent_plastic_strain,
                                 half_point.equivalent_plastic_strain);
                accumulate_error(equivalent_creep,
                                 full_point.equivalent_creep_strain,
                                 half_point.equivalent_creep_strain);

                const AxisymmetricStressValues& full_value =
                    full_stress[element][q];
                const AxisymmetricStressValues& half_value =
                    half_stress[element][q];
                accumulate_error(stress, full_value.rr, half_value.rr);
                accumulate_error(stress, full_value.zz, half_value.zz);
                accumulate_error(stress, full_value.hoop, half_value.hoop);
                accumulate_error(stress, full_value.rz, half_value.rz);
            }
        }
    }
    if (full_step.contact_histories.size() !=
        two_half_steps.contact_histories.size())
        throw std::logic_error(
            "step-doubling contact-history layouts differ");
    for (std::size_t contact = 0;
         contact < full_step.contact_histories.size(); ++contact) {
        if (full_step.contact_histories[contact].size() !=
            two_half_steps.contact_histories[contact].size())
            throw std::logic_error(
                "step-doubling contact-node history layouts differ");
        for (std::size_t node = 0;
             node < full_step.contact_histories[contact].size(); ++node) {
            const ContactPointHistory& full =
                full_step.contact_histories[contact][node];
            const ContactPointHistory& half =
                two_half_steps.contact_histories[contact][node];
            accumulate_error(contact_friction,
                             full.elastic_tangential_slip,
                             half.elastic_tangential_slip);
            contact_state_mismatch =
                contact_state_mismatch || full.sliding != half.sliding;
        }
    }

    TransientTimeErrorEstimate result;
    result.temperature = normalized_error(
        nodal[0], options.temperature_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.radial_displacement = normalized_error(
        nodal[1], options.displacement_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.axial_displacement = normalized_error(
        nodal[2], options.displacement_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.elastic_strain = normalized_error(
        elastic, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.plastic_strain = normalized_error(
        plastic, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.creep_strain = normalized_error(
        creep, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.equivalent_plastic_strain = normalized_error(
        equivalent_plastic, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.equivalent_creep_strain = normalized_error(
        equivalent_creep, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.stress = normalized_error(
        stress, options.stress_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.contact_friction =
        contact_state_mismatch
            ? std::numeric_limits<double>::infinity()
            : normalized_error(
                  contact_friction,
                  options.displacement_time_absolute_tolerance,
                  options.time_error_relative_tolerance);
    result.maximum = std::max(
        {result.temperature, result.radial_displacement,
         result.axial_displacement, result.elastic_strain,
         result.plastic_strain, result.creep_strain,
         result.equivalent_plastic_strain,
         result.equivalent_creep_strain, result.stress,
         result.contact_friction});
    return result;
}

double time_error_step_factor(const TransientTimeOptions& options,
                              double error) {
    if (!(error > 0.0))
        return options.growth_factor;
    return std::clamp(options.time_error_safety_factor / std::sqrt(error),
                      0.1, options.growth_factor);
}

} // namespace

SteadyResult solve_steady(SteadyProblem& problem,
                          const SteadyLoadOptions& load_options,
                          const SolverOptions& options) {
    if (load_options.load_steps == 0)
        throw std::invalid_argument("solve_steady load_steps must be positive");
    if (!std::isfinite(load_options.cutback_factor) ||
        !(load_options.cutback_factor > 0.0 &&
          load_options.cutback_factor < 1.0))
        throw std::invalid_argument(
            "solve_steady cutback factor must lie between zero and one");
    if (!std::isfinite(load_options.minimum_load_increment) ||
        !(load_options.minimum_load_increment > 0.0) ||
        load_options.minimum_load_increment > 1.0)
        throw std::invalid_argument(
            "solve_steady minimum load increment must lie in (0, 1]");
    const SteadyClock::time_point start = SteadyClock::now();
    SteadyResult result;
    PetscSolver solver;
    std::vector<double> state = problem.initial_state();
    double accepted_load_factor = 0.0;
    for (std::size_t step = 1; step <= load_options.load_steps; ++step) {
        const double target_load_factor =
            static_cast<double>(step) /
            static_cast<double>(load_options.load_steps);
        while (accepted_load_factor < target_load_factor) {
            std::size_t cutbacks = 0;
            double attempted_load_factor = target_load_factor;
            double load_increment =
                attempted_load_factor - accepted_load_factor;
            SolveResult attempt;
            for (;;) {
                attempt = SolveResult{};
                try {
                    problem.set_load_factor(attempted_load_factor);
                    attempt = solver.solve(
                        problem,
                        initial_guess_with_dirichlet_values(problem, state),
                        options);
                    if (attempt.converged)
                        problem.commit_contact_state(attempt.state);
                } catch (const std::domain_error& error) {
                    attempt.converged = false;
                    attempt.failure_category =
                        SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                } catch (const std::overflow_error& error) {
                    attempt.converged = false;
                    attempt.failure_category =
                        SolveFailureCategory::physical_domain;
                    attempt.failure_message = error.what();
                }
                accumulate_timing(result.aggregate_timing, attempt.timing);
                result.total_nonlinear_iterations +=
                    attempt.nonlinear_iterations;
                result.total_linear_iterations += attempt.linear_iterations;
                if (attempt.converged)
                    break;

                result.rejected_steps.push_back(
                    {attempted_load_factor, load_increment, cutbacks,
                     attempt.failure_category, attempt.failure_message});
                result.solve = attempt;
                if (cutbacks >=
                    load_options.maximum_cutbacks_per_step) {
                    result.total_seconds = seconds_since(start);
                    return result;
                }
                const double tolerance =
                    16.0 * std::numeric_limits<double>::epsilon();
                if (load_increment <=
                    load_options.minimum_load_increment + tolerance) {
                    result.total_seconds = seconds_since(start);
                    return result;
                }
                load_increment = std::max(
                    load_increment * load_options.cutback_factor,
                    load_options.minimum_load_increment);
                attempted_load_factor = accepted_load_factor + load_increment;
                ++cutbacks;
                ++result.total_cutbacks;
            }
            accepted_load_factor = attempted_load_factor;
            state = attempt.state;
            result.solve = std::move(attempt);
        }
        result.completed_steps = step;
    }
    result.completed = true;
    result.total_seconds = seconds_since(start);
    return result;
}

TransientResult solve_transient(TransientProblem& problem,
                                const TransientTimeOptions& options,
                                const SolverOptions& solver_options,
                                TransientStepObserver* observer) {
    validate_time_options(problem, options);
    const SteadyClock::time_point start = SteadyClock::now();
    TransientResult result;
    PetscSolver solver;
    const std::vector<double> events = problem.time_events();
    double next_time_step = options.initial_time_step;
    while (!reaches_end(problem.committed_time(), options.end_time)) {
        const double controller_time_step = std::min(
            next_time_step, options.end_time - problem.committed_time());
        double time_step = controller_time_step;
        bool event_truncated = false;
        for (const double event : events) {
            if (reaches_end(problem.committed_time(), event))
                continue;
            if (event >= options.end_time ||
                reaches_end(event, options.end_time))
                break;
            const double event_step = event - problem.committed_time();
            if (event_step < time_step &&
                !reaches_end(event, problem.committed_time() + time_step)) {
                time_step = event_step;
                event_truncated = true;
            }
            break;
        }
        std::size_t cutbacks = 0;
        for (;;) {
            const double end_time = problem.committed_time() + time_step;
            SolveResult attempt;
            double time_error_estimate = 0.0;
            TransientTimeErrorEstimate time_error_components;
            TransientConservationSummary first_half_conservation;
            int controller_nonlinear_iterations = 0;
            const bool error_control =
                options.time_error_relative_tolerance > 0.0;
            TransientCommittedState base_state;
            if (error_control)
                base_state = problem.committed_state();
            try {
                if (!error_control) {
                    TimeStepTransaction transaction(
                        problem,
                        {end_time, load_factor_at_time(options, end_time)});
                    attempt = solver.solve(
                        problem,
                        initial_guess_with_dirichlet_values(
                            problem, problem.committed_solution()),
                        solver_options);
                    controller_nonlinear_iterations =
                        attempt.nonlinear_iterations;
                    if (attempt.converged)
                        transaction.commit(attempt.state);
                } else {
                    SolveResult full_step;
                    TransientCommittedState full_step_state;
                    {
                        TimeStepTransaction transaction(
                            problem,
                            {end_time,
                             load_factor_at_time(options, end_time)});
                        full_step = solver.solve(
                            problem,
                            initial_guess_with_dirichlet_values(
                                problem, problem.committed_solution()),
                            solver_options);
                        if (full_step.converged) {
                            transaction.commit(full_step.state);
                            full_step_state = problem.committed_state();
                        }
                    }
                    attempt = full_step;
                    controller_nonlinear_iterations =
                        full_step.nonlinear_iterations;
                    if (full_step.converged) {
                        problem.restore_committed_state(base_state);
                        const double half_time =
                            base_state.time + 0.5 * time_step;
                        SolveResult first_half;
                        {
                            TimeStepTransaction transaction(
                                problem,
                                {half_time,
                                 load_factor_at_time(options, half_time)});
                            first_half = solver.solve(
                                problem,
                                initial_guess_with_dirichlet_values(
                                    problem, problem.committed_solution()),
                                solver_options);
                            if (first_half.converged)
                                transaction.commit(first_half.state);
                        }
                        if (first_half.converged)
                            first_half_conservation =
                                problem.last_conservation_summary();
                        controller_nonlinear_iterations = std::max(
                            controller_nonlinear_iterations,
                            first_half.nonlinear_iterations);
                        combine_attempt(attempt, first_half);
                        if (!first_half.converged) {
                            problem.restore_committed_state(
                                std::move(base_state));
                        } else {
                            SolveResult second_half;
                            {
                                TimeStepTransaction transaction(
                                    problem,
                                    {end_time,
                                     load_factor_at_time(options, end_time)});
                                second_half = solver.solve(
                                    problem,
                                    initial_guess_with_dirichlet_values(
                                        problem,
                                        problem.committed_solution()),
                                    solver_options);
                                if (second_half.converged)
                                    transaction.commit(second_half.state);
                            }
                            controller_nonlinear_iterations = std::max(
                                controller_nonlinear_iterations,
                                second_half.nonlinear_iterations);
                            combine_attempt(attempt, second_half);
                            if (!second_half.converged) {
                                problem.restore_committed_state(
                                    std::move(base_state));
                            } else {
                                time_error_components = step_doubling_error(
                                    full_step_state,
                                    problem.committed_state(), options);
                                time_error_estimate =
                                    time_error_components.maximum;
                                if (!(time_error_estimate <= 1.0)) {
                                    problem.restore_committed_state(
                                        std::move(base_state));
                                    attempt.converged = false;
                                    attempt.failure_category =
                                        SolveFailureCategory::
                                            time_discretization;
                                    attempt.failure_message =
                                        "Backward-Euler step-doubling error "
                                        "exceeded one";
                                    ++result.time_error_rejections;
                                } else {
                                    TransientCommittedState accepted_state =
                                        problem.committed_state();
                                    accepted_state.conservation =
                                        combine_half_step_conservation(
                                            first_half_conservation,
                                            accepted_state.conservation);
                                    problem.restore_committed_state(
                                        std::move(accepted_state));
                                }
                            }
                        }
                    }
                }
            } catch (const std::domain_error& error) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                attempt.converged = false;
                attempt.failure_category =
                    SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (const std::overflow_error& error) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                attempt.converged = false;
                attempt.failure_category =
                    SolveFailureCategory::physical_domain;
                attempt.failure_message = error.what();
            } catch (...) {
                if (error_control && !base_state.solution.empty())
                    problem.restore_committed_state(std::move(base_state));
                throw;
            }
            accumulate_timing(result.aggregate_timing, attempt.timing);
            result.total_nonlinear_iterations += attempt.nonlinear_iterations;
            result.total_linear_iterations += attempt.linear_iterations;
            result.last_attempt = std::move(attempt);
            if (result.last_attempt.converged) {
                std::vector<RegionInelasticSummary> histories;
                histories.reserve(problem.region_count());
                for (std::size_t region = 0; region < problem.region_count();
                     ++region)
                    histories.push_back(
                        problem.summarize_region_history(region));
                next_time_step = accepted_next_time_step(
                    options, time_step, controller_time_step, event_truncated,
                    cutbacks, controller_nonlinear_iterations);
                if (error_control) {
                    const double error_limited_step =
                        time_step *
                        time_error_step_factor(options, time_error_estimate);
                    next_time_step =
                        std::clamp(std::min(next_time_step, error_limited_step),
                                   options.minimum_time_step,
                                   options.maximum_time_step);
                }
                result.accepted_steps.push_back(
                    {problem.committed_time(), time_step, next_time_step,
                     problem.committed_load_factor(), cutbacks,
                     result.last_attempt.nonlinear_iterations,
                     result.last_attempt.linear_iterations,
                     std::move(histories), time_error_estimate,
                     time_error_components,
                     problem.last_conservation_summary()});
                if (observer != nullptr)
                    observer->accepted_step(problem,
                                            result.accepted_steps.back());
                break;
            }
            result.rejected_steps.push_back(
                {problem.committed_time() + time_step, time_step, cutbacks,
                 result.last_attempt.nonlinear_iterations,
                 result.last_attempt.linear_iterations,
                 result.last_attempt.convergence_reason,
                 result.last_attempt.residual_norm,
                 result.last_attempt.failure_category,
                 result.last_attempt.failure_message, time_error_estimate,
                 time_error_components});
            if (cutbacks >= options.maximum_cutbacks_per_step) {
                result.termination_reason =
                    TransientTerminationReason::maximum_cutbacks;
                break;
            }
            double reduction_factor = options.cutback_factor;
            if (result.last_attempt.failure_category ==
                SolveFailureCategory::time_discretization)
                reduction_factor = std::min(
                    0.9,
                    std::max(options.cutback_factor,
                             time_error_step_factor(options,
                                                    time_error_estimate)));
            const double reduced = time_step * reduction_factor;
            const double minimum_tolerance =
                16.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, options.minimum_time_step);
            if (time_step <= options.minimum_time_step + minimum_tolerance) {
                result.termination_reason =
                    TransientTerminationReason::minimum_time_step;
                break;
            }
            time_step = std::max(reduced, options.minimum_time_step);
            ++cutbacks;
            ++result.total_cutbacks;
        }
        if (!result.last_attempt.converged)
            break;
    }
    result.completed = reaches_end(problem.committed_time(), options.end_time);
    if (result.completed)
        result.termination_reason = TransientTerminationReason::completed;
    result.committed_state = problem.committed_solution();
    result.committed_time = problem.committed_time();
    result.next_time_step = next_time_step;
    result.total_seconds = seconds_since(start);
    return result;
}

const char*
transient_termination_reason_name(TransientTerminationReason reason) noexcept {
    switch (reason) {
    case TransientTerminationReason::not_started:
        return "not_started";
    case TransientTerminationReason::completed:
        return "completed";
    case TransientTerminationReason::maximum_cutbacks:
        return "maximum_cutbacks";
    case TransientTerminationReason::minimum_time_step:
        return "minimum_time_step";
    }
    return "unknown";
}

} // namespace fuelsim
