#include "contact_types.hpp"
#include "core/cax_evaluation.hpp"
#include "core/element_region_data.hpp"
#ifndef FUELSIM_TEST_RZ_PROBLEM_ACCESS_HPP
#define FUELSIM_TEST_RZ_PROBLEM_ACCESS_HPP
#include "core/problem_backend_access.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <stdexcept>

namespace fuelsim::rz {
inline std::size_t named_region(const SpatialDefinition& definition, const std::string& name) {
    const auto found = std::find_if(definition.regions.begin(),
        definition.regions.end(),
        [&name](const RegionDefinition& region) { return region.name == name; });
    if (found == definition.regions.end())
        throw std::invalid_argument("Unknown region: " + name);
    return static_cast<std::size_t>(found - definition.regions.begin());
}

class ProblemAccess final {
  public:
    static SteadyBackendView view(const SteadyProblem& problem) { return fuelsim::BackendAccess::steady(problem); }

    static TransientBackendView view(const TransientProblem& problem) {
        return fuelsim::BackendAccess::transient(problem);
    }

    static const SpatialDefinition& definition(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.definition();
    }

    static const spatial_detail::SpatialLayout& dof_map(const SteadyProblem& problem) noexcept {
        return view(problem).spatial;
    }

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

    static const AxisymmetricRegionData& region_kernel_data(const SteadyProblem& problem, std::size_t index) {
        return view(problem).kernel_data.at(index);
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

    static const Quad4RzGeometry&
    region_element_geometry(const SteadyProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }

    static std::size_t contact_count(const SteadyProblem& problem) noexcept {
        return view(problem).spatial.definition().contacts.size();
    }

    static const ContactDefinition& contact(const SteadyProblem& problem, std::size_t index) {
        return view(problem).spatial.definition().contacts.at(index);
    }

    static const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories(
        const SteadyProblem& problem) noexcept {
        return view(problem).spatial.committed_contact_histories();
    }

    static std::vector<ContactNodeSummary>
    summarize_contact_nodes(const SteadyProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_contact_nodes(contact, state);
    }

    static std::vector<std::size_t> contact_secondary_source_nodes(const SteadyProblem& problem, std::size_t contact) {
        return view(problem).spatial.contact_secondary_source_nodes(contact);
    }

    static InterfaceSummary
    summarize_interface(const SteadyProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_interface(contact, state);
    }

    static Cax4LocalDofs contribution_dofs(const SteadyProblem& problem, std::size_t contribution) {
        return view(problem).spatial.contribution_dofs(contribution);
    }

    static Cax4LocalValues
    contribution_state(const SteadyProblem& problem, std::size_t contribution, const std::vector<double>& state) {
        if (state.size() != problem.dof_count())
            throw std::invalid_argument("SteadyProblem contribution state has the wrong global size");
        Cax4LocalValues result{};
        const Cax4LocalDofs dofs = contribution_dofs(problem, contribution);
        for (std::size_t local = 0; local < dofs.size(); ++local)
            result[local] = state.at(dofs[local]);
        return result;
    }

    static Cax4LocalResidual
    contribution_residual(const SteadyProblem& problem, std::size_t contribution, const Cax4LocalValues& state) {
        const SteadyBackendView backend = view(problem);
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.compute_contribution(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return evaluate_cax4(backend.kernel_data[location.first],
            backend.spatial.region_element_geometry(location.first, location.second),
            state,
            {},
            nullptr,
            0.0,
            false,
            {true, false, false, false})
            .residual;
    }

    static LocalLinearization
    linearize_contribution(const SteadyProblem& problem, std::size_t contribution, const Cax4LocalValues& state) {
        const SteadyBackendView backend = view(problem);
        LocalLinearization result{};
        if (contribution >= backend.spatial.volume_contribution_count())
            result.residual = backend.spatial.compute_contribution(contribution, state, &result.jacobian);
        else {
            const auto location = backend.spatial.element_location(contribution);
            const auto volume = evaluate_cax4(backend.kernel_data[location.first],
                backend.spatial.region_element_geometry(location.first, location.second),
                state,
                {},
                nullptr,
                0.0,
                false,
                {true, true, false, false});
            result.residual = volume.residual;
            result.jacobian = volume.jacobian;
        }
        return result;
    }

    static const SpatialDefinition& definition(const TransientProblem& problem) noexcept {
        return problem.definition();
    }

    static const spatial_detail::SpatialLayout& dof_map(const TransientProblem& problem) noexcept {
        return view(problem).spatial;
    }

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

    static const AxisymmetricRegionData& region_kernel_data(const TransientProblem& problem, std::size_t index) {
        return view(problem).kernel_data.at(index);
    }

    static const Quad4RzGeometry&
    region_element_geometry(const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }

    static TransientCommittedState committed_state(const TransientProblem& problem) {
        return fuelsim::BackendAccess::committed_state(problem);
    }

    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
        fuelsim::BackendAccess::restore_committed_state(problem, std::move(state));
    }

    static const Quad4MaterialHistory&
    material_history(const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).histories.at(region).at(element);
    }

    static std::array<AxisymmetricStressValues, 4>
    material_stress(const TransientProblem& problem, std::size_t region, std::size_t element) {
        std::array<AxisymmetricStressValues, 4> result{};
        const Quad4MaterialHistory& history = view(problem).histories.at(region).at(element);
        for (std::size_t q = 0; q < result.size(); ++q)
            result[q] = history[q].stress;
        return result;
    }

    static RegionStateSummary summarize_region_history(const TransientProblem& problem, std::size_t region) {
        return problem.summarize_region(region);
    }

    static InterfaceSummary
    summarize_interface(const TransientProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_interface(contact, state);
    }

    static std::vector<ContactNodeSummary>
    summarize_contact_nodes(const TransientProblem& problem, std::size_t contact, const std::vector<double>& state) {
        return view(problem).spatial.summarize_contact_nodes(contact, state);
    }

    static std::vector<std::size_t> contact_secondary_source_nodes(const TransientProblem& problem,
        std::size_t contact) {
        return view(problem).spatial.contact_secondary_source_nodes(contact);
    }

    static Cax4LocalDofs contribution_dofs(const TransientProblem& problem, std::size_t contribution) {
        return view(problem).spatial.contribution_dofs(contribution);
    }

    static Cax4LocalValues
    contribution_state(const TransientProblem& problem, std::size_t contribution, const std::vector<double>& state) {
        if (state.size() != problem.dof_count())
            throw std::invalid_argument("TransientProblem contribution state has the wrong global size");
        Cax4LocalValues result{};
        const Cax4LocalDofs dofs = contribution_dofs(problem, contribution);
        for (std::size_t local = 0; local < dofs.size(); ++local)
            result[local] = state.at(dofs[local]);
        return result;
    }

    static Cax4LocalResidual
    contribution_residual(const TransientProblem& problem, std::size_t contribution, const Cax4LocalValues& state) {
        const TransientBackendView backend = view(problem);
        if (!backend.time_step_active)
            throw std::logic_error("TransientProblem residual evaluation requires an active time step");
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.compute_contribution(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return evaluate_cax4(backend.kernel_data[location.first],
            backend.spatial.region_element_geometry(location.first, location.second),
            state,
            contribution_state(problem, contribution, backend.committed_solution),
            &backend.histories[location.first][location.second],
            backend.active_time_step,
            true,
            {true, false, false, false})
            .residual;
    }

    static LocalLinearization
    linearize_contribution(const TransientProblem& problem, std::size_t contribution, const Cax4LocalValues& state) {
        const TransientBackendView backend = view(problem);
        if (!backend.time_step_active)
            throw std::logic_error("TransientProblem residual evaluation requires an active time step");
        LocalLinearization result{};
        if (contribution >= backend.spatial.volume_contribution_count())
            result.residual = backend.spatial.compute_contribution(contribution, state, &result.jacobian);
        else {
            const auto location = backend.spatial.element_location(contribution);
            const auto volume = evaluate_cax4(backend.kernel_data[location.first],
                backend.spatial.region_element_geometry(location.first, location.second),
                state,
                contribution_state(problem, contribution, backend.committed_solution),
                &backend.histories[location.first][location.second],
                backend.active_time_step,
                true,
                {true, true, false, false});
            result.residual = volume.residual;
            result.jacobian = volume.jacobian;
        }
        return result;
    }
};
} // namespace fuelsim::rz
#endif
