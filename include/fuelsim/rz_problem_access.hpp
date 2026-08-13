#ifndef FUELSIM_RZ_PROBLEM_ACCESS_HPP
#define FUELSIM_RZ_PROBLEM_ACCESS_HPP

#include "fuelsim/dof_map.hpp"
#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz_thermoelastic.hpp"
#include "fuelsim/quad4_rz_transient.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim {
namespace rz {

// Concrete RZ checkpoint and verification state. The common time integrator
// uses TransientStateSnapshot and does not depend on this layout.
struct TransientCommittedState final {
    std::vector<double> solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> material_stresses;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    TransientConservationSummary conservation;
    double time = 0.0;
    double load_factor = 0.0;
};

// Explicit access to the concrete RZ backend for I/O, verification, and
// backend-specific diagnostics. Solver and time-control code must use the
// geometry-independent SteadyProblem and TransientProblem ports instead.
class ProblemAccess final {
  public:
    static const SpatialDefinition& definition(const SteadyProblem& problem) noexcept;
    static const DofMap& dof_map(const SteadyProblem& problem) noexcept;
    static std::size_t region_count(const SteadyProblem& problem) noexcept;
    static std::size_t region_index(const SteadyProblem& problem, const std::string& name);
    static const RegionDefinition& region(const SteadyProblem& problem, std::size_t region_index);
    static const RegionMesh& region_mesh(const SteadyProblem& problem, std::size_t region_index);
    static const Quad4RzThermoelasticKernel& region_kernel(const SteadyProblem& problem, std::size_t region_index);
    static std::size_t region_node_offset(const SteadyProblem& problem, std::size_t region_index);
    static std::size_t region_element_count(const SteadyProblem& problem, std::size_t region_index);
    static std::size_t region_element_offset(const SteadyProblem& problem, std::size_t region_index);
    static std::size_t volume_contribution_count(const SteadyProblem& problem) noexcept;
    static SpatialContributionType contribution_type(const SteadyProblem& problem, std::size_t contribution_index);
    static const Quad4RzGeometry& region_element_geometry(const SteadyProblem& problem, std::size_t region_index,
                                                          std::size_t element_index);
    static std::size_t contact_count(const SteadyProblem& problem) noexcept;
    static const ContactDefinition& contact(const SteadyProblem& problem, std::size_t contact_index);
    static const std::vector<std::vector<ContactPointHistory>>&
    committed_contact_histories(const SteadyProblem& problem) noexcept;
    static void commit_contact_state(SteadyProblem& problem, const std::vector<double>& state);
    static void restore_contact_state(SteadyProblem& problem, const std::vector<double>& state,
                                      std::vector<std::vector<ContactPointHistory>> histories);
    static std::vector<ContactNodeSummary>
    summarize_contact_nodes(const SteadyProblem& problem, std::size_t contact_index, const std::vector<double>& state);
    static std::vector<std::size_t> contact_secondary_source_nodes(const SteadyProblem& problem,
                                                                   std::size_t contact_index);
    static InterfaceSummary summarize_interface(const SteadyProblem& problem, std::size_t contact_index,
                                                const std::vector<double>& state);
    static LocalDofs contribution_dofs(const SteadyProblem& problem, std::size_t contribution_index);
    static LocalValues contribution_state(const SteadyProblem& problem, std::size_t contribution_index,
                                          const std::vector<double>& global_state);
    static LocalValues contribution_state(const SteadyProblem& problem, std::size_t contribution_index,
                                          const GlobalStateView& global_state);
    static LocalResidual contribution_residual(const SteadyProblem& problem, std::size_t contribution_index,
                                               const LocalValues& state);
    static LocalSystem linearize_contribution(const SteadyProblem& problem, std::size_t contribution_index,
                                              const LocalValues& state);

    static const TransientProblemDefinition& definition(const TransientProblem& problem) noexcept;
    static const DofMap& dof_map(const TransientProblem& problem) noexcept;
    static std::size_t region_count(const TransientProblem& problem) noexcept;
    static std::size_t region_index(const TransientProblem& problem, const std::string& name);
    static std::size_t region_node_offset(const TransientProblem& problem, std::size_t region_index);
    static const RegionDefinition& region(const TransientProblem& problem, std::size_t region_index);
    static const RegionMesh& region_mesh(const TransientProblem& problem, std::size_t region_index);
    static const Quad4RzTransientKernel& region_kernel(const TransientProblem& problem, std::size_t region_index);
    static const Quad4RzGeometry& region_element_geometry(const TransientProblem& problem, std::size_t region_index,
                                                          std::size_t element_index);
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
    static const Quad4MaterialHistory& material_history(const TransientProblem& problem, std::size_t region_index,
                                                        std::size_t element_index);
    static const std::array<AxisymmetricStressValues, 4>&
    material_stress(const TransientProblem& problem, std::size_t region_index, std::size_t element_index);
    static RegionInelasticSummary summarize_region_history(const TransientProblem& problem, std::size_t region_index);
    static InterfaceSummary summarize_interface(const TransientProblem& problem, std::size_t contact_index,
                                                const std::vector<double>& state);
    static std::vector<ContactNodeSummary> summarize_contact_nodes(const TransientProblem& problem,
                                                                   std::size_t contact_index,
                                                                   const std::vector<double>& state);
    static std::vector<std::size_t> contact_secondary_source_nodes(const TransientProblem& problem,
                                                                   std::size_t contact_index);
    static LocalDofs contribution_dofs(const TransientProblem& problem, std::size_t contribution_index);
    static LocalValues contribution_state(const TransientProblem& problem, std::size_t contribution_index,
                                          const std::vector<double>& global_state);
    static LocalValues contribution_state(const TransientProblem& problem, std::size_t contribution_index,
                                          const GlobalStateView& global_state);
    static LocalResidual contribution_residual(const TransientProblem& problem, std::size_t contribution_index,
                                               const LocalValues& state);
    static LocalSystem linearize_contribution(const TransientProblem& problem, std::size_t contribution_index,
                                              const LocalValues& state);
};

} // namespace rz
} // namespace fuelsim

#endif
