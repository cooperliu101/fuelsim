#ifndef FUELSIM_TRANSIENT_CONSERVATION_HPP
#define FUELSIM_TRANSIENT_CONSERVATION_HPP

#include "fuelsim/transient_problem.hpp"

namespace fuelsim {

class TransientConservationCalculator final {
  public:
    static TransientConservationSummary summarize(
        const TransientProblem& problem,
        const std::vector<double>& converged_solution,
        const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
        const std::vector<
            std::vector<std::array<AxisymmetricStressValues, 4>>>&
            staged_stresses);
};

} // namespace fuelsim

#endif
