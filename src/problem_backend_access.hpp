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
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> material_stresses;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    TransientConservationSummary conservation;
    double time = 0.0, load_factor = 0.0;
};
namespace rz {
struct SteadyBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<Quad4RzThermoelasticKernel>& kernels;
};
struct TransientBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<Quad4RzTransientKernel>& kernels;
    const std::vector<std::vector<Quad4MaterialHistory>>& histories;
    const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>& stresses;
    const std::vector<double>& committed_solution;
    double active_time_step;
    bool time_step_active;
};
class BackendAccess final {
  public:
    static SteadyBackendView steady(const SteadyProblem& problem) noexcept;
    static TransientBackendView transient(const TransientProblem& problem) noexcept;
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace rz
namespace cartesian {
class BackendAccess final {
  public:
    static const SpatialAssembly& spatial(const SteadyProblem& problem) noexcept;
    static const SpatialAssembly& spatial(const TransientProblem& problem) noexcept;
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace cartesian
} // namespace fuelsim
