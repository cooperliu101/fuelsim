#ifndef FUELSIM_PROBLEM_BACKEND_ACCESS_HPP
#define FUELSIM_PROBLEM_BACKEND_ACCESS_HPP
#include "assembly.hpp"
#include "cartesian3d_assembly.hpp"
#include "fuelsim/quad4_rz_thermoelastic.hpp"
#include "fuelsim/quad4_rz_transient.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"
#include <array>
#include <vector>
namespace fuelsim {
namespace rz {
struct SteadyBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<Quad4RzThermoelasticKernel>& kernels;
};
struct TransientBackendView final {
    const TransientProblemDefinition& definition;
    const SpatialAssembly& spatial;
    const std::vector<Quad4RzTransientKernel>& kernels;
    const std::vector<std::vector<Quad4MaterialHistory>>& histories;
    const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>& stresses;
    const std::vector<double>& committed_solution;
    double active_time_step;
    bool time_step_active;
};
struct TransientCommittedState final {
    std::vector<double> solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> material_stresses;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    TransientConservationSummary conservation;
    double time = 0.0;
    double load_factor = 0.0;
};
class BackendAccess final {
  public:
    static SteadyBackendView steady(const SteadyProblem& problem) noexcept;
    static TransientBackendView transient(const TransientProblem& problem) noexcept;
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace rz
namespace cartesian3d {
struct SteadyBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<Hex8ThermoelasticKernel>& kernels;
};
struct TransientBackendView final {
    const TransientProblemDefinition& definition;
    const SpatialAssembly& spatial;
    const std::vector<Hex8ThermoelasticKernel>& kernels;
    const std::vector<double>& committed_solution;
    double active_time_step;
    bool time_step_active;
};
struct TransientCommittedState final {
    std::vector<double> solution;
    TransientConservationSummary conservation;
    double time = 0.0;
    double load_factor = 0.0;
};
class BackendAccess final {
  public:
    static SteadyBackendView steady(const SteadyProblem& problem) noexcept;
    static TransientBackendView transient(const TransientProblem& problem) noexcept;
    static std::array<SymmetricTensor3Values, 8> stress(
        const TransientProblem& problem, const std::vector<double>& state, std::size_t region, std::size_t element);
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace cartesian3d
} // namespace fuelsim
#endif
