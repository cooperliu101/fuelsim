#pragma once
#include "assembly.hpp"
#include "cartesian3d_assembly.hpp"
#include "fuelsim/quad4_rz.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"
#include <array>
#include <vector>

namespace fuelsim {
struct TransientCommittedState final {
    std::vector<double> solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<Hex8MaterialHistory>> cartesian_material_histories;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    TransientConservationSummary conservation;
    double time = 0.0, load_factor = 0.0;
};

namespace rz {
struct SteadyBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<Quad4RzData>& kernel_data;
};

struct TransientBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<Quad4RzData>& kernel_data;
    const std::vector<std::vector<Quad4MaterialHistory>>& histories;
    const std::vector<double>& committed_solution;
    double active_time_step;
    bool time_step_active;
};
} // namespace rz

class BackendAccess final {
  public:
    static rz::SteadyBackendView steady(const SteadyProblem& problem) noexcept;
    static rz::TransientBackendView transient(const TransientProblem& problem) noexcept;
    static const cartesian::SpatialAssembly& cartesian_spatial(const SteadyProblem& problem) noexcept;
    static const cartesian::SpatialAssembly& cartesian_spatial(const TransientProblem& problem) noexcept;
    static const std::vector<std::vector<Hex8MaterialHistory>>& cartesian_material_histories(
        const TransientProblem& problem) noexcept;
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace fuelsim
