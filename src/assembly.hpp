#ifndef FUELSIM_ASSEMBLY_HPP
#define FUELSIM_ASSEMBLY_HPP

#include "fuelsim/boundary.hpp"
#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz.hpp"
#include "fuelsim/spatial_definition.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {

class SpatialLayout final {
  public:
    struct ResolvedBoundary final {
        std::size_t region;
        RegionBoundary boundary;
    };

    SpatialLayout(SpatialDefinition definition, std::vector<std::int64_t> block_ids, std::vector<RegionMesh> meshes);

    const SpatialDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;
    std::size_t region_count() const noexcept;
    const RegionDefinition& region(std::size_t index) const;
    const RegionMesh& region_mesh(std::size_t index) const;
    std::size_t region_node_offset(std::size_t index) const;
    std::size_t region_element_offset(std::size_t index) const;
    std::size_t volume_contribution_count() const noexcept;
    const Quad4RzGeometry& element_geometry(std::size_t region, std::size_t element) const;
    ResolvedBoundary resolve_boundary(const UnstructuredQuad4Mesh& source_mesh, const std::string& name) const;
    std::size_t global_node(std::size_t region, std::size_t local_node) const;
    std::pair<std::size_t, std::array<std::size_t, 2>> edge_parent(std::size_t region,
                                                                   const Line2BoundaryElement& edge) const;

  private:
    friend class SpatialAssembly;

    void build_volume_geometries();

    SpatialDefinition _definition;
    std::vector<std::int64_t> _block_ids;
    std::vector<RegionMesh> _meshes;
    std::vector<std::size_t> _node_offsets;
    std::vector<std::size_t> _element_offsets;
    DofMap _dof_map;
    std::vector<std::vector<Quad4RzGeometry>> _region_geometries;
};

class BoundaryAssembly final {
  public:
    BoundaryAssembly() = default;

    void build(const UnstructuredQuad4Mesh& source_mesh, const SpatialLayout& layout);
    void set_load_factor(double load_factor, const SpatialLayout& layout);
    double load_factor() const noexcept;
    void set_time(double time, const SpatialLayout& layout);
    double region_heat_source(std::size_t region_index, const SpatialLayout& layout) const;
    std::size_t pressure_contribution_count() const noexcept;
    std::size_t traction_contribution_count() const noexcept;
    std::size_t convection_contribution_count() const noexcept;
    LocalDofs contribution_dofs(SpatialContributionType type, std::size_t index, const SpatialLayout& layout) const;
    LocalResidual contribution_residual(SpatialContributionType type, std::size_t index,
                                        const LocalValues& state) const;
    LocalSystem linearize_contribution(SpatialContributionType type, std::size_t index, const LocalValues& state) const;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept;

  private:
    double function_value(const std::string& name, const SpatialLayout& layout) const;
    double load_multiplier(bool scale_with_load, const std::string& function, const SpatialLayout& layout) const;
    void refresh_controlled_values(const SpatialLayout& layout);

    struct PressureLoad final {
        double pressure;
        bool scale_with_load;
        std::string function;
    };
    struct TractionLoad final {
        double traction;
        bool scale_with_load;
        std::string function;
    };
    struct ControlledDirichlet final {
        std::size_t dof;
        double value;
        bool scale_with_load;
        std::string function;
    };
    struct ConvectionLoad final {
        double heat_transfer_coefficient;
        double ambient_temperature;
        std::string coefficient_function;
        std::string ambient_temperature_function;
    };
    struct PressureContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzPressureGeometry geometry;
    };
    struct TractionContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzTractionGeometry geometry;
    };
    struct ConvectionContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzConvectionGeometry geometry;
    };

    std::vector<Line2RzConvectionKernel> _convection_kernels;
    std::vector<ConvectionLoad> _convection_loads;
    std::vector<ConvectionContribution> _convection_contributions;
    std::vector<Line2RzPressureKernel> _pressure_kernels;
    std::vector<PressureContribution> _pressure_contributions;
    std::vector<Line2RzTractionKernel> _traction_kernels;
    std::vector<TractionContribution> _traction_contributions;
    std::vector<DirichletCondition> _dirichlet_conditions;
    std::vector<ControlledDirichlet> _controlled_dirichlet_conditions;
    std::vector<PressureLoad> _pressure_loads;
    std::vector<TractionLoad> _traction_loads;
    double _load_factor = 1.0;
    double _time = 0.0;
};

class SpatialAssembly;

class ContactAssembly final {
  public:
    ContactAssembly() = default;

    std::size_t contact_count(const SpatialLayout& layout) const noexcept;
    const ContactDefinition& contact(std::size_t contact_index, const SpatialLayout& layout) const;
    const std::vector<std::vector<ContactPointHistory>>& committed_histories() const noexcept;
    bool uses_augmented_contact(const SpatialLayout& layout) const noexcept;
    AugmentedContactUpdate update_augmented_multipliers(SpatialAssembly& assembly, const std::vector<double>& state,
                                                        std::size_t completed_updates);
    void commit_state(SpatialAssembly& assembly, const std::vector<double>& state);
    void restore_state(const SpatialAssembly& assembly, const std::vector<double>& state,
                       std::vector<std::vector<ContactPointHistory>> histories);

  private:
    friend class SpatialAssembly;

    struct ThermalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        Line2RzHeatGeometry geometry;
    };
    struct MechanicalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        NodeToLineRzContactGeometry geometry;
        std::size_t secondary;
        std::size_t primary;
        mutable bool active = false;
    };

    std::vector<Line2RzGapHeatKernel> _thermal_kernels;
    std::vector<NodeToLineRzContactKernel> _mechanical_kernels;
    std::vector<ThermalContribution> _thermal_contributions;
    std::vector<MechanicalContribution> _mechanical_contributions;
    mutable std::vector<std::vector<bool>> _projected_mechanical_nodes;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<SpatialLayout::ResolvedBoundary> _primary_boundaries;
    std::vector<SpatialLayout::ResolvedBoundary> _secondary_boundaries;
};

// Internal shared spatial assembly. It deliberately does not implement the
// nonlinear-solver port: SteadyProblem and TransientProblem retain ownership
// of their volume physics and expose the two production problem types.
class SpatialAssembly final {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh);

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
    SpatialContributionType contribution_type(std::size_t contribution_index) const;
    std::pair<std::size_t, std::size_t> element_location(std::size_t contribution_index) const;
    const Quad4RzGeometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    double region_heat_source(std::size_t region_index) const;

    std::size_t contact_count() const noexcept;
    const ContactDefinition& contact(std::size_t contact_index) const;
    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept;
    void commit_contact_state(const std::vector<double>& state);
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                std::size_t completed_updates);
    void restore_contact_state(const std::vector<double>& state,
                               std::vector<std::vector<ContactPointHistory>> histories);

    void set_load_factor(double load_factor);
    double load_factor() const noexcept;
    void set_time(double time);

    std::vector<double> initial_state() const;
    std::vector<ContactNodeSummary> summarize_contact_nodes(std::size_t contact_index,
                                                            const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;

    std::size_t dof_count() const noexcept;
    std::size_t contribution_count() const noexcept;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept;
    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t contribution_begin, std::size_t contribution_end) const;
    void validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                              const GlobalStateView& state) const;
    LocalDofs contribution_dofs(std::size_t contribution_index) const;
    LocalResidual contribution_residual(std::size_t contribution_index, const LocalValues& state) const;
    LocalSystem linearize_contribution(std::size_t contribution_index, const LocalValues& state) const;

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

    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh,
                    std::vector<std::int64_t> block_ids, std::vector<RegionMesh> meshes);
    static std::vector<std::int64_t> resolve_block_ids(const SpatialDefinition& definition,
                                                       const UnstructuredQuad4Mesh& source_mesh);
    static std::vector<RegionMesh> build_meshes(const SpatialDefinition& definition,
                                                const UnstructuredQuad4Mesh& source_mesh);

    void build_contacts(const UnstructuredQuad4Mesh& source_mesh);
    void update_mechanical_candidates(const std::vector<double>& state) const;
    void update_mechanical_candidates(std::size_t contribution_begin, std::size_t contribution_end,
                                      const GlobalStateView& state) const;
    std::vector<std::vector<bool>> touched_mechanical_nodes(std::size_t contribution_begin,
                                                            std::size_t contribution_end) const;
    ContributionRanges contribution_ranges() const noexcept;
    ContributionLocation locate_contribution(std::size_t contribution_index) const;
    LocalValues contribution_state(std::size_t contribution_index, const std::vector<double>& global_state) const;
    LocalValues contribution_state(std::size_t contribution_index, const GlobalStateView& global_state) const;
    SpatialLayout _layout;
    BoundaryAssembly _boundary;
    ContactAssembly _contact;
};

} // namespace fuelsim

#endif
