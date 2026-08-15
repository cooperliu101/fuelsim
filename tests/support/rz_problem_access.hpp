#ifndef FUELSIM_TEST_RZ_PROBLEM_ACCESS_HPP
#define FUELSIM_TEST_RZ_PROBLEM_ACCESS_HPP
#include "problem_backend_access.hpp"
#include <algorithm>
#include <stdexcept>

namespace fuelsim::rz {
struct LocalLinearization final {
    LocalResidual residual;
    LocalJacobian jacobian;
};

inline LocalLinearization linearize_quad4_rz_thermoelastic(
    const Quad4RzData& data, const Quad4RzGeometry& geometry, const LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_quad4_rz_thermoelastic(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_quad4_rz_transient(const Quad4RzData& data, const Quad4RzGeometry& geometry,
    const LocalValues& state, const LocalValues& committed_state, const Quad4MaterialHistory& history,
    double time_step) {
    LocalLinearization result{};
    result.residual =
        compute_quad4_rz_transient(data, geometry, state, committed_state, history, time_step, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_line2_rz_boundary(
    const Line2RzBoundaryData& data, const Line2RzBoundaryGeometry& geometry, const LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_line2_rz_boundary(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_line2_rz_gap_heat(
    const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry, const LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_line2_rz_gap_heat(properties, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history) {
    LocalLinearization result{};
    result.residual =
        compute_node_to_line_rz_contact(properties, geometry, state, committed_state, history, &result.jacobian);
    return result;
}

inline std::size_t named_region(const SpatialDefinition& definition, const std::string& name) {
    const auto found = std::find_if(definition.regions.begin(), definition.regions.end(),
        [&name](const RegionDefinition& region) { return region.name == name; });
    if (found == definition.regions.end()) throw std::invalid_argument("Unknown region: " + name);
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

    static const Quad4RzData& region_kernel_data(const SteadyProblem& problem, std::size_t index) {
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

    static const Quad4RzGeometry& region_element_geometry(
        const SteadyProblem& problem, std::size_t region, std::size_t element) {
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
        if (state.size() != problem.dof_count())
            throw std::invalid_argument("SteadyProblem contribution state has the wrong global size");
        LocalValues result{};
        const LocalDofs dofs = contribution_dofs(problem, contribution);
        for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = state.at(dofs[local]);
        return result;
    }

    static LocalResidual contribution_residual(
        const SteadyProblem& problem, std::size_t contribution, const LocalValues& state) {
        const SteadyBackendView backend = view(problem);
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.compute_contribution(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return compute_quad4_rz_thermoelastic(backend.kernel_data[location.first],
            backend.spatial.region_element_geometry(location.first, location.second), state);
    }

    static LocalLinearization linearize_contribution(
        const SteadyProblem& problem, std::size_t contribution, const LocalValues& state) {
        const SteadyBackendView backend = view(problem);
        LocalLinearization result{};
        if (contribution >= backend.spatial.volume_contribution_count())
            result.residual = backend.spatial.compute_contribution(contribution, state, &result.jacobian);
        else {
            const auto location = backend.spatial.element_location(contribution);
            result.residual = compute_quad4_rz_thermoelastic(backend.kernel_data[location.first],
                backend.spatial.region_element_geometry(location.first, location.second), state, &result.jacobian);
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

    static const Quad4RzData& region_kernel_data(const TransientProblem& problem, std::size_t index) {
        return view(problem).kernel_data.at(index);
    }

    static const Quad4RzGeometry& region_element_geometry(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).spatial.region_element_geometry(region, element);
    }

    static TransientCommittedState committed_state(const TransientProblem& problem) {
        return fuelsim::BackendAccess::committed_state(problem);
    }

    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
        fuelsim::BackendAccess::restore_committed_state(problem, std::move(state));
    }

    static const Quad4MaterialHistory& material_history(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        return view(problem).histories.at(region).at(element);
    }

    static std::array<AxisymmetricStressValues, 4> material_stress(
        const TransientProblem& problem, std::size_t region, std::size_t element) {
        std::array<AxisymmetricStressValues, 4> result{};
        const Quad4MaterialHistory& history = view(problem).histories.at(region).at(element);
        for (std::size_t q = 0; q < result.size(); ++q) result[q] = history[q].stress;
        return result;
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
        if (state.size() != problem.dof_count())
            throw std::invalid_argument("TransientProblem contribution state has the wrong global size");
        LocalValues result{};
        const LocalDofs dofs = contribution_dofs(problem, contribution);
        for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = state.at(dofs[local]);
        return result;
    }

    static LocalResidual contribution_residual(
        const TransientProblem& problem, std::size_t contribution, const LocalValues& state) {
        const TransientBackendView backend = view(problem);
        if (!backend.time_step_active)
            throw std::logic_error("TransientProblem residual evaluation requires an active time step");
        if (contribution >= backend.spatial.volume_contribution_count())
            return backend.spatial.compute_contribution(contribution, state);
        const auto location = backend.spatial.element_location(contribution);
        return compute_quad4_rz_transient(backend.kernel_data[location.first],
            backend.spatial.region_element_geometry(location.first, location.second), state,
            contribution_state(problem, contribution, backend.committed_solution),
            backend.histories[location.first][location.second], backend.active_time_step);
    }

    static LocalLinearization linearize_contribution(
        const TransientProblem& problem, std::size_t contribution, const LocalValues& state) {
        const TransientBackendView backend = view(problem);
        if (!backend.time_step_active)
            throw std::logic_error("TransientProblem residual evaluation requires an active time step");
        LocalLinearization result{};
        if (contribution >= backend.spatial.volume_contribution_count())
            result.residual = backend.spatial.compute_contribution(contribution, state, &result.jacobian);
        else {
            const auto location = backend.spatial.element_location(contribution);
            result.residual = compute_quad4_rz_transient(backend.kernel_data[location.first],
                backend.spatial.region_element_geometry(location.first, location.second), state,
                contribution_state(problem, contribution, backend.committed_solution),
                backend.histories[location.first][location.second], backend.active_time_step, &result.jacobian);
        }
        return result;
    }
};
} // namespace fuelsim::rz
#endif
