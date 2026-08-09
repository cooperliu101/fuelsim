#ifndef FUELSIM_TRANSIENT_OUTPUT_OBSERVER_HPP
#define FUELSIM_TRANSIENT_OUTPUT_OBSERVER_HPP

#include "case_output.hpp"

#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"

#include <cstddef>
#include <string>

namespace fuelsim {

class PetscSession;

namespace app {

class TransientOutputObserver final : public TransientStepObserver {
  public:
    TransientOutputObserver(ExodusTransientResultsWriter* results,
                            EngineeringHistoryWriter* history,
                            std::string checkpoint_file,
                            std::size_t exodus_interval,
                            std::size_t history_interval,
                            std::size_t progress_interval,
                            std::size_t checkpoint_interval,
                            const PetscSession& session,
                            CaseOutput& progress_output);

    void accepted_step(const TransientProblem& problem,
                       const TransientAcceptedStep& step) override;
    void finalize(const TransientProblem& problem, double next_time_step);

  private:
    ExodusTransientResultsWriter* _results;
    EngineeringHistoryWriter* _history;
    std::string _checkpoint_file;
    std::size_t _exodus_interval;
    std::size_t _history_interval;
    std::size_t _progress_interval;
    std::size_t _checkpoint_interval;
    std::size_t _accepted_steps;
    bool _exodus_at_latest;
    bool _history_at_latest;
    double _last_time_step;
    double _last_next_time_step;
    int _last_nonlinear_iterations;
    bool _checkpoint_at_latest;
    const PetscSession& _session;
    CaseOutput& _progress_output;
};

} // namespace app
} // namespace fuelsim

#endif
