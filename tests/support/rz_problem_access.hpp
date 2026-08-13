#ifndef FUELSIM_TEST_RZ_PROBLEM_ACCESS_HPP
#define FUELSIM_TEST_RZ_PROBLEM_ACCESS_HPP
#include "problem_backend_access.hpp"
#include <algorithm>
#include <stdexcept>
namespace fuelsim::rz {
inline std::size_t named_region(const SpatialDefinition& definition, const std::string& name) {
    const auto found = std::find_if(definition.regions.begin(), definition.regions.end(),
        [&name](const RegionDefinition& region) { return region.name == name; });
    if (found == definition.regions.end()) throw std::invalid_argument("Unknown region: " + name);
    return static_cast<std::size_t>(found - definition.regions.begin());
}
class ProblemAccess final {
  public:
    static SteadyBackendView view(const SteadyProblem& problem) { return BackendAccess::steady(problem); }
    static TransientBackendView view(const TransientProblem& problem) { return BackendAccess::transient(problem); }
    static const SpatialDefinition& definition(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.definition();
    }
    static const DofMap& dof_map(const SteadyProblem& problem) noexcept { return view(problem).spatial.dof_map(); }
    static std::size_t region_count(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.region_count();
    }
    static std::size_t region_index(const SteadyProblem& problem, const std::string& name) {
        return named_region(view(problem).spatial.definition(), name);
    }
    static const RegionDefinition& region(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region(index);
    }
    static const RegionMesh& region_mesh(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region_mesh(index);
    }
    static const Quad4RzThermoelasticKernel& region_kernel(const SteadyProblem& problem, std::size_t index) {
        return view(problem).kernels.at(index);
    }
    static std::size_t region_node_offset(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region_node_offset(index);
    }
    static std::size_t region_element_count(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region_element_count(index);
    }
    static std::size_t region_element_offset(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.region_element_offset(index);
    }
    static std::size_t volume_contribution_count(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.volume_contribution_count();
    }
    static SpatialContributionType contribution_type(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.contribution_type(index);
    }
    static const Quad4RzGeometry& region_element_geometry(
        const SteadyProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }
    static std::size_t contact_count(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.contact_count();
    }
    static const ContactDefinition& contact(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.contact(index);
    }
    static const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories(
        const SteadyProblem& problem) noexcept {
        return view(problem).spatial.committed_contact_histories();
    }
    static std::vector<ContactNodeSummary> summarize_contact_nodes(
        const SteadyProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_contact_nodes(contact, state);
    }
    static std::vector<std::size_t> contact_secondary_source_nodes(const SteadyProblem& problem, std::size_t contact) {
        return view(problem).spatial.contact_secondary_source_nodes(contact);
    }
    static InterfaceSummary summarize_interface(
        const SteadyProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_interface(contact, state);
    }
    static LocalDofs contribution_dofs(const SteadyProblem& problem, std::size_t contribution) {
        return view(problem).spatial.contribution_dofs(contribution);
    }
    static LocalValues contribution_state(
        const SteadyProblem& problem, std::size_t contribution, const std::vector<double>& state) {
        return contribution_state(problem, contribution, GlobalStateView(state));
    }
    static LocalValues contribution_state(
        const SteadyProblem& problem, std::size_t contribution, const GlobalStateView& state) {
        if (state.global_size() != problem.dof_count())
            throw std::invalid_argument("SteadyProblem contribution state has the wrong global size");
        LocalValues result{};
        const LocalDofs dofs = contribution_dofs(problem, contribution);
        for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = state.value(dofs[local]);
        return result;
    }
    static LocalResidual contribution_residual(
        const SteadyProblem& problem, std::size_t contribution, const LocalValues& state) {
        const SteadyBackendView backend = view(problem);
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.contribution_residual(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return backend.kernels[location.first].residual(
            backend.spatial.region_element_geometry(location.first, location.second), state);
    }
    static LocalSystem linearize_contribution(
        const SteadyProblem& problem, std::size_t contribution, const LocalValues& state) {
        const SteadyBackendView backend = view(problem);
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.linearize_contribution(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return backend.kernels[location.first].linearize(
            backend.spatial.region_element_geometry(location.first, location.second), state);
    }
    static const SpatialDefinition& definition(const TransientProblem& problem) noexcept {
        return problem.definition();
    }
    static const DofMap& dof_map(const TransientProblem& problem) noexcept { return view(problem).spatial.dof_map(); }
    static std::size_t region_count(const TransientProblem& problem) noexcept {
        return view(problem).spatial.region_count();
    }
    static std::size_t region_index(const TransientProblem& problem, const std::string& name) {
        return named_region(view(problem).spatial.definition(), name);
    }
    static std::size_t region_node_offset(const TransientProblem& problem, std::size_t index) {
        return view(problem).spatial.region_node_offset(index);
    }
    static const RegionDefinition& region(const TransientProblem& problem, std::size_t index) {
        return view(problem).spatial.region(index);
    }
    static const RegionMesh& region_mesh(const TransientProblem& problem, std::size_t index) {
        return view(problem).spatial.region_mesh(index);
    }
    static const Quad4RzTransientKernel& region_kernel(const TransientProblem& problem, std::size_t index) {
        return view(problem).kernels.at(index);
    }
    static const Quad4RzGeometry& region_element_geometry(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }
    static TransientCommittedState committed_state(const TransientProblem& problem) {
        return BackendAccess::committed_state(problem);
    }
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
        BackendAccess::restore_committed_state(problem, std::move(state));
    }
    static const Quad4MaterialHistory& material_history(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).histories.at(region).at(element);
    }
    static const std::array<AxisymmetricStressValues, 4>& material_stress(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).stresses.at(region).at(element);
    }
    static RegionStateSummary summarize_region_history(const TransientProblem& problem, std::size_t region) {
        return problem.summarize_region(region);
    }
    static InterfaceSummary summarize_interface(
        const TransientProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_interface(contact, state);
    }
    static std::vector<ContactNodeSummary> summarize_contact_nodes(
        const TransientProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_contact_nodes(contact, state);
    }
    static std::vector<std::size_t> contact_secondary_source_nodes(
        const TransientProblem& problem, std::size_t contact) {
        return view(problem).spatial.contact_secondary_source_nodes(contact);
    }
    static LocalDofs contribution_dofs(const TransientProblem& problem, std::size_t contribution) {
        return view(problem).spatial.contribution_dofs(contribution);
    }
    static LocalValues contribution_state(
        const TransientProblem& problem, std::size_t contribution, const std::vector<double>& state) {
        return contribution_state(problem, contribution, GlobalStateView(state));
    }
    static LocalValues contribution_state(
        const TransientProblem& problem, std::size_t contribution, const GlobalStateView& state) {
        if (state.global_size() != problem.dof_count())
            throw std::invalid_argument("TransientProblem contribution state has the wrong global size");
        LocalValues result{};
        const LocalDofs dofs = contribution_dofs(problem, contribution);
        for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = state.value(dofs[local]);
        return result;
    }
    static LocalResidual contribution_residual(
        const TransientProblem& problem, std::size_t contribution, const LocalValues& state) {
        const TransientBackendView backend = view(problem);
        if (!backend.time_step_active)
            throw std::logic_error("TransientProblem residual evaluation requires an active time step");
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.contribution_residual(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return backend.kernels[location.first].residual(
            backend.spatial.region_element_geometry(location.first, location.second), state,
            contribution_state(problem, contribution, backend.committed_solution),
            backend.histories[location.first][location.second], backend.active_time_step);
    }
    static LocalSystem linearize_contribution(
        const TransientProblem& problem, std::size_t contribution, const LocalValues& state) {
        const TransientBackendView backend = view(problem);
        if (!backend.time_step_active)
            throw std::logic_error("TransientProblem residual evaluation requires an active time step");
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.linearize_contribution(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return backend.kernels[location.first].linearize(
            backend.spatial.region_element_geometry(location.first, location.second), state,
            contribution_state(problem, contribution, backend.committed_solution),
            backend.histories[location.first][location.second], backend.active_time_step);
    }
};
} // namespace fuelsim::rz
#endif
