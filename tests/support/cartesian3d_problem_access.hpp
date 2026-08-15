#ifndef FUELSIM_TEST_CARTESIAN3D_PROBLEM_ACCESS_HPP
#define FUELSIM_TEST_CARTESIAN3D_PROBLEM_ACCESS_HPP
#include "problem_backend_access.hpp"

namespace fuelsim::cartesian {
class ProblemAccess final {
  public:
    static const SpatialAssembly& view(const SteadyProblem& problem) {
        return fuelsim::BackendAccess::cartesian_spatial(problem);
    }

    static const SpatialAssembly& view(const TransientProblem& problem) {
        return fuelsim::BackendAccess::cartesian_spatial(problem);
    }

    static const spatial_detail::SpatialLayout& dof_map(const SteadyProblem& problem) noexcept { return view(problem); }

    static std::size_t region_count(const SteadyProblem& problem) noexcept { return view(problem).region_count(); }

    static const RegionDefinition& region(const SteadyProblem& problem, std::size_t index) {
        return view(problem).region(index);
    }

    static const Hex8RegionMesh& region_mesh(const SteadyProblem& problem, std::size_t index) {
        return view(problem).region_mesh(index);
    }

    static std::size_t region_node_offset(const SteadyProblem& problem, std::size_t index) {
        return view(problem).region_node_offset(index);
    }

    static std::vector<CartesianContactNodeSummary> summarize_contact_nodes(
        const SteadyProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).summarize_contact_nodes(contact, state);
    }

    static std::vector<std::size_t> contact_secondary_source_nodes(const SteadyProblem& problem, std::size_t contact) {
        return view(problem).contact_secondary_source_nodes(contact);
    }

    static InterfaceSummary summarize_interface(
        const SteadyProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).summarize_interface(contact, state);
    }

    static const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories(
        const SteadyProblem& problem) noexcept {
        return view(problem).committed_contact_histories();
    }

    static const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories(
        const TransientProblem& problem) noexcept {
        return view(problem).committed_contact_histories();
    }

    static const Hex8Geometry& region_element_geometry(
        const SteadyProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).region_element_geometry(region, element);
    }

    static const SpatialDefinition& definition(const TransientProblem& problem) noexcept {
        return problem.definition();
    }

    static const spatial_detail::SpatialLayout& dof_map(const TransientProblem& problem) noexcept {
        return view(problem);
    }

    static std::size_t region_count(const TransientProblem& problem) noexcept { return view(problem).region_count(); }

    static const RegionDefinition& region(const TransientProblem& problem, std::size_t index) {
        return view(problem).region(index);
    }

    static const Hex8RegionMesh& region_mesh(const TransientProblem& problem, std::size_t index) {
        return view(problem).region_mesh(index);
    }

    static std::size_t region_node_offset(const TransientProblem& problem, std::size_t index) {
        return view(problem).region_node_offset(index);
    }

    static std::vector<CartesianContactNodeSummary> summarize_contact_nodes(
        const TransientProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).summarize_contact_nodes(contact, state);
    }

    static std::vector<std::size_t> contact_secondary_source_nodes(
        const TransientProblem& problem, std::size_t contact) {
        return view(problem).contact_secondary_source_nodes(contact);
    }

    static InterfaceSummary summarize_interface(
        const TransientProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).summarize_interface(contact, state);
    }

    static const Hex8Geometry& region_element_geometry(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).region_element_geometry(region, element);
    }

    static std::array<SymmetricTensor3Values, 8> stress(
        const SteadyProblem& problem, const std::vector<double>& state, std::size_t region, std::size_t element) {
        return view(problem).stress(region, element, state);
    }

    static std::array<SymmetricTensor3Values, 8> stress(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        std::array<SymmetricTensor3Values, 8> result{};
        const Hex8MaterialHistory& history = material_history(problem, region, element);
        for (std::size_t q = 0; q < result.size(); ++q) result[q] = history[q].stress;
        return result;
    }

    static const Hex8MaterialHistory& material_history(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return fuelsim::BackendAccess::cartesian_material_histories(problem).at(region).at(element);
    }

    static TransientCommittedState committed_state(const TransientProblem& problem) {
        return fuelsim::BackendAccess::committed_state(problem);
    }

    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
        fuelsim::BackendAccess::restore_committed_state(problem, std::move(state));
    }
};
} // namespace fuelsim::cartesian
#endif
