#ifndef FUELSIM_TEST_CARTESIAN3D_PROBLEM_ACCESS_HPP
#define FUELSIM_TEST_CARTESIAN3D_PROBLEM_ACCESS_HPP
#include "problem_backend_access.hpp"
namespace fuelsim::cartesian {
class ProblemAccess final {
  public:
    static SteadyBackendView view(const SteadyProblem& problem) { return BackendAccess::steady(problem); }
    static TransientBackendView view(const TransientProblem& problem) { return BackendAccess::transient(problem); }
    static const DofMap& dof_map(const SteadyProblem& problem) noexcept { return view(problem).spatial.dof_map(); }
    static std::size_t region_count(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.region_count();
    }
    static const RegionDefinition& region(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region(index);
    }
    static const Hex8RegionMesh& region_mesh(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region_mesh(index);
    }
    static std::size_t region_node_offset(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region_node_offset(index);
    }
    static const Hex8Geometry& region_element_geometry(
        const SteadyProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }
    static const TransientProblemDefinition& definition(const TransientProblem& problem) noexcept {
        return view(problem).definition;
    }
    static const DofMap& dof_map(const TransientProblem& problem) noexcept { return view(problem).spatial.dof_map(); }
    static std::size_t region_count(const TransientProblem& problem) noexcept {
        return view(problem).spatial.region_count();
    }
    static const RegionDefinition& region(const TransientProblem& problem, std::size_t index) {
        return view(problem).spatial.region(index);
    }
    static const Hex8RegionMesh& region_mesh(const TransientProblem& problem, std::size_t index) {
        return view(problem).spatial.region_mesh(index);
    }
    static std::size_t region_node_offset(const TransientProblem& problem, std::size_t index) {
        return view(problem).spatial.region_node_offset(index);
    }
    static const Hex8Geometry& region_element_geometry(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }
    static std::array<SymmetricTensor3Values, 8> stress(
        const SteadyProblem& problem, const std::vector<double>& state, std::size_t region, std::size_t element) {
        return view(problem).spatial.stress(region, element, state);
    }
    static std::array<SymmetricTensor3Values, 8> stress(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return BackendAccess::stress(problem, problem.committed_solution(), region, element);
    }
    static TransientCommittedState committed_state(const TransientProblem& problem) {
        return BackendAccess::committed_state(problem);
    }
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
        BackendAccess::restore_committed_state(problem, std::move(state));
    }
};
} // namespace fuelsim::cartesian
#endif
