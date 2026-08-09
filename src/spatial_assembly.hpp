#ifndef FUELSIM_SPATIAL_ASSEMBLY_HPP
#define FUELSIM_SPATIAL_ASSEMBLY_HPP

#include "fuelsim/boundary.hpp"
#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "boundary_assembly.hpp"
#include "contact_assembly.hpp"
#include "spatial_layout.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {

// Internal shared spatial assembly. It deliberately does not implement the
// nonlinear-solver port: SteadyProblem and TransientProblem retain ownership
// of their volume physics and expose the two production problem types.
class SpatialAssembly final {
  public:
    SpatialAssembly(SpatialDefinition definition,
                    const UnstructuredQuad4Mesh& source_mesh);

    const SpatialDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;
    std::size_t region_count() const noexcept;
    std::size_t region_index(const std::string& name) const;
    const RegionDefinition& region(std::size_t region_index) const;
    const RegionMesh& region_mesh(std::size_t region_index) const;
    std::size_t region_node_offset(std::size_t region_index) const;
    std::size_t region_element_count(std::size_t region_index) const;
    std::size_t region_element_offset(std::size_t region_index) const;
    std::size_t volume_contribution_count() const noexcept;
    SpatialContributionType
    contribution_type(std::size_t contribution_index) const;
    std::pair<std::size_t, std::size_t>
    element_location(std::size_t contribution_index) const;
    const Quad4RzGeometry&
    region_element_geometry(std::size_t region_index,
                            std::size_t element_index) const;
    double region_heat_source(std::size_t region_index) const;

    std::size_t contact_count() const noexcept;
    const ContactDefinition& contact(std::size_t contact_index) const;
    const std::vector<std::vector<ContactPointHistory>>&
    committed_contact_histories() const noexcept;
    void commit_contact_state(const std::vector<double>& state);
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate
    update_augmented_contact_multipliers(const std::vector<double>& state,
                                         std::size_t completed_updates);
    void restore_contact_state(
        const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories);

    void set_load_factor(double load_factor);
    double load_factor() const noexcept;
    void set_time(double time);

    std::vector<double> initial_state() const;
    std::vector<ContactNodeSummary>
    summarize_contact_nodes(std::size_t contact_index,
                            const std::vector<double>& state) const;
    std::vector<std::size_t>
    contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary
    summarize_interface(std::size_t contact_index,
                        const std::vector<double>& state) const;

    std::size_t dof_count() const noexcept;
    std::size_t contribution_count() const noexcept;
    const std::vector<DirichletCondition>&
    dirichlet_conditions() const noexcept;
    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t>
    required_state_dofs(std::size_t contribution_begin,
                        std::size_t contribution_end) const;
    void validate_local_state(std::size_t contribution_begin,
                              std::size_t contribution_end,
                              const GlobalStateView& state) const;
    LocalDofs contribution_dofs(std::size_t contribution_index) const;
    LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const LocalValues& state) const;
    LocalSystem linearize_contribution(std::size_t contribution_index,
                                       const LocalValues& state) const;

  private:
    friend class ContactAssembly;

    struct ContributionRanges final {
        std::size_t thermal_begin;
        std::size_t mechanical_begin;
        std::size_t pressure_begin;
        std::size_t traction_begin;
        std::size_t convection_begin;
        std::size_t end;
    };
    struct ContributionLocation final {
        SpatialContributionType type;
        std::size_t local_index;
    };
    using ResolvedBoundary = SpatialLayout::ResolvedBoundary;
    using ThermalContribution = ContactAssembly::ThermalContribution;
    using MechanicalContribution = ContactAssembly::MechanicalContribution;
    using PressureContribution = BoundaryAssembly::PressureContribution;
    using TractionContribution = BoundaryAssembly::TractionContribution;
    using ConvectionContribution = BoundaryAssembly::ConvectionContribution;
    using ControlledDirichlet = BoundaryAssembly::ControlledDirichlet;
    using PressureLoad = BoundaryAssembly::PressureLoad;
    using TractionLoad = BoundaryAssembly::TractionLoad;
    using ConvectionLoad = BoundaryAssembly::ConvectionLoad;

    SpatialAssembly(SpatialDefinition definition,
                    const UnstructuredQuad4Mesh& source_mesh,
                    std::vector<std::int64_t> block_ids,
                    std::vector<RegionMesh> meshes);
    static std::vector<std::int64_t>
    resolve_block_ids(const SpatialDefinition& definition,
                      const UnstructuredQuad4Mesh& source_mesh);
    static std::vector<RegionMesh>
    build_meshes(const SpatialDefinition& definition,
                 const UnstructuredQuad4Mesh& source_mesh);

    ResolvedBoundary resolve_boundary(const UnstructuredQuad4Mesh& source_mesh,
                                      const std::string& name) const;
    std::size_t global_node(std::size_t region_index,
                            std::size_t local_node) const;
    void build_contacts(const UnstructuredQuad4Mesh& source_mesh);
    void build_boundary_conditions(const UnstructuredQuad4Mesh& source_mesh);
    void update_mechanical_candidates(const std::vector<double>& state) const;
    void update_mechanical_candidates(std::size_t contribution_begin,
                                      std::size_t contribution_end,
                                      const GlobalStateView& state) const;
    ContributionRanges contribution_ranges() const noexcept;
    ContributionLocation
    locate_contribution(std::size_t contribution_index) const;
    LocalValues contribution_state(
        std::size_t contribution_index,
        const std::vector<double>& global_state) const;
    LocalValues contribution_state(std::size_t contribution_index,
                                   const GlobalStateView& global_state) const;
    SpatialLayout _layout;
    BoundaryAssembly _boundary;
    ContactAssembly _contact;
};

} // namespace fuelsim

#endif
