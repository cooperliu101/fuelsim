#include "transient_output_observer.hpp"

#include "fuelsim/checkpoint_io.hpp"
#include "fuelsim/petsc_solver.hpp"

#include <utility>

namespace fuelsim::app {

TransientOutputObserver::TransientOutputObserver(
    ExodusTransientResultsWriter* results,
    EngineeringHistoryWriter* history,
    std::string checkpoint_file,
    std::size_t exodus_interval,
    std::size_t history_interval,
    std::size_t progress_interval,
    std::size_t checkpoint_interval,
    const PetscSession& session,
    CaseOutput& progress_output)
    : _results(results), _history(history),
      _checkpoint_file(std::move(checkpoint_file)),
      _exodus_interval(exodus_interval), _history_interval(history_interval),
      _progress_interval(progress_interval),
      _checkpoint_interval(checkpoint_interval), _accepted_steps(0),
      _exodus_at_latest(true), _history_at_latest(true),
      _last_time_step(0.0), _last_next_time_step(0.0),
      _last_nonlinear_iterations(0), _checkpoint_at_latest(false),
      _session(session), _progress_output(progress_output) {}

void TransientOutputObserver::accepted_step(
    const TransientProblem& problem, const TransientAcceptedStep& step) {
    ++_accepted_steps;
    _last_time_step = step.time_step;
    _last_next_time_step = step.next_time_step;
    _last_nonlinear_iterations = step.nonlinear_iterations;
    _exodus_at_latest = false;
    _history_at_latest = false;
    _checkpoint_at_latest = false;
    _session.collective_root_action([&]() {
        if (_results != nullptr && _accepted_steps % _exodus_interval == 0) {
            _results->append(problem);
            _exodus_at_latest = true;
        }
        if (_history != nullptr && _accepted_steps % _history_interval == 0) {
            _history->append(problem, step.time_step, step.next_time_step,
                             step.nonlinear_iterations);
            _history_at_latest = true;
        }
        if (_accepted_steps % _progress_interval == 0) {
            _progress_output.value("progress.accepted_steps", _accepted_steps);
            _progress_output.value("progress.time", step.time);
            _progress_output.value("progress.time_step", step.time_step);
            _progress_output.value("progress.next_time_step",
                                   step.next_time_step);
            _progress_output.value("progress.nonlinear_iterations",
                                   step.nonlinear_iterations);
            _progress_output.value("progress.linear_iterations",
                                   step.linear_iterations);
            _progress_output.value("progress.cutbacks", step.cutbacks);
            _progress_output.value("progress.time_error_estimate",
                                   step.time_error_estimate);
            write_time_error_components("progress.time_error.",
                                        step.time_error_components,
                                        _progress_output);
            write_conservation_summary("progress.conservation.",
                                       step.conservation, _progress_output);
        }
        if (!_checkpoint_file.empty() &&
            _accepted_steps % _checkpoint_interval == 0) {
            TransientCheckpointIo::write(_checkpoint_file, problem,
                                         step.next_time_step);
            _checkpoint_at_latest = true;
        }
    });
}

void TransientOutputObserver::finalize(const TransientProblem& problem,
                                       double next_time_step) {
    _session.collective_root_action([&]() {
        if (_results != nullptr && !_exodus_at_latest)
            _results->append(problem);
        if (_history != nullptr && !_history_at_latest)
            _history->append(problem, _last_time_step, _last_next_time_step,
                             _last_nonlinear_iterations);
        if (!_checkpoint_file.empty() && !_checkpoint_at_latest)
            TransientCheckpointIo::write(_checkpoint_file, problem,
                                         next_time_step);
    });
    _exodus_at_latest = true;
    _history_at_latest = true;
    _checkpoint_at_latest = true;
}

} // namespace fuelsim::app
