#pragma once
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
namespace fuelsim::rz {
struct ResolvedBoundary final {
    std::size_t region;
    RegionBoundary boundary;
};
class SpatialAssembly final {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh);
    const SpatialDefinition& definition() const noexcept { return _definition; }
    const DofMap& dof_map() const noexcept { return _dof_map; }
    std::size_t region_count() const noexcept { return _definition.regions.size(); }
    std::size_t region_index(const std::string& name) const;
    const RegionDefinition& region(std::size_t index) const { return _definition.regions.at(index); }
    const RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }
    std::size_t region_node_offset(std::size_t region_index) const;
    std::size_t region_element_count(std::size_t region_index) const;
    std::size_t region_element_offset(std::size_t region_index) const;
    std::size_t volume_contribution_count() const noexcept { return _element_offsets.back(); }
    SpatialContributionType contribution_type(std::size_t index) const;
    std::pair<std::size_t, std::size_t> element_location(std::size_t index) const;
    const Quad4RzGeometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    double region_heat_source(std::size_t region_index) const;
    std::size_t contact_count() const noexcept { return _definition.contacts.size(); }
    const ContactDefinition& contact(std::size_t index) const { return _definition.contacts.at(index); }
    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept {
        return _contact_histories;
    }
    void commit_contact_state(const std::vector<double>& state);
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(
        const std::vector<double>& state, std::size_t completed_updates);
    void restore_contact_state(
        const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories);
    void set_load_factor(double load_factor);
    double load_factor() const noexcept { return _load_factor; }
    void set_time(double time);
    std::vector<double> initial_state() const;
    std::vector<ContactNodeSummary> summarize_contact_nodes(
        std::size_t contact_index, const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;
    std::size_t dof_count() const noexcept { return _dof_map.dof_count(); }
    std::size_t contribution_count() const noexcept { return contribution_ranges().end; }
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept { return _dirichlet_conditions; }
    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const GlobalStateView& state) const;
    LocalDofs contribution_dofs(std::size_t index) const;
    LocalResidual contribution_residual(std::size_t index, const LocalValues& state) const;
    LocalSystem linearize_contribution(std::size_t index, const LocalValues& state) const;
    std::size_t global_node(std::size_t region, std::size_t local_node) const;
    std::pair<std::size_t, std::array<std::size_t, 2>> edge_parent(
        std::size_t region, const Line2BoundaryElement& edge) const;

  private:
    struct ContributionRanges final {
        std::size_t thermal_begin, mechanical_begin, pressure_begin, traction_begin, convection_begin, end;
    };
    struct ContributionLocation final {
        SpatialContributionType type;
        std::size_t local_index;
    };
    struct ControlledScalar final {
        double value;
        bool scale_with_load;
        std::string function;
    };
    struct ControlledDirichlet final {
        std::size_t dof;
        ControlledScalar control;
    };
    struct ConvectionLoad final {
        double heat_transfer_coefficient, ambient_temperature;
        std::string coefficient_function, ambient_temperature_function;
    };
    struct PressureContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzBoundaryGeometry geometry;
    };
    struct TractionContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzBoundaryGeometry geometry;
    };
    struct ConvectionContribution final {
        std::size_t load;
        std::array<std::size_t, 4> nodes;
        Line2RzBoundaryGeometry geometry;
    };
    struct ThermalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        Line2RzHeatPointGeometry geometry;
        std::size_t integration_point, primary;
        mutable bool active = false;
    };
    struct MechanicalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        NodeToLineRzContactGeometry geometry;
        std::size_t secondary, primary;
        mutable bool active = false;
    };
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh,
        std::vector<std::int64_t> block_ids, std::vector<RegionMesh> meshes);
    static std::vector<std::int64_t> resolve_block_ids(
        const SpatialDefinition& definition, const UnstructuredQuad4Mesh& source_mesh);
    static std::vector<RegionMesh> build_meshes(
        const SpatialDefinition& definition, const UnstructuredQuad4Mesh& source_mesh);
    void build_volume_geometries();
    ResolvedBoundary resolve_boundary(const UnstructuredQuad4Mesh& source_mesh, const std::string& name) const;
    void build_boundaries(const UnstructuredQuad4Mesh& source_mesh);
    double function_value(const std::string& name) const;
    double controlled_value(const ControlledScalar& control) const;
    void refresh_controlled_values();
    LocalDofs boundary_contribution_dofs(SpatialContributionType type, std::size_t index) const;
    LocalResidual boundary_contribution_residual(
        SpatialContributionType type, std::size_t index, const LocalValues& state) const;
    LocalSystem linearize_boundary_contribution(
        SpatialContributionType type, std::size_t index, const LocalValues& state) const;
    LocalDofs local_dofs(const std::array<std::size_t, 4>& nodes) const;
    void build_contacts(const UnstructuredQuad4Mesh& source_mesh);
    void update_thermal_candidates(const std::vector<double>& state) const;
    void update_thermal_candidates(std::size_t first, std::size_t last, const GlobalStateView& state) const;
    void update_mechanical_candidates(const std::vector<double>& state) const;
    void update_mechanical_candidates(std::size_t first, std::size_t last, const GlobalStateView& state) const;
    std::vector<std::vector<bool>> touched_thermal_points(std::size_t first, std::size_t last) const;
    std::vector<std::vector<bool>> touched_mechanical_nodes(std::size_t first, std::size_t last) const;
    ContributionRanges contribution_ranges() const noexcept;
    ContributionLocation locate_contribution(std::size_t index) const;
    LocalValues contribution_state(std::size_t index, const std::vector<double>& global_state) const;
    LocalValues contribution_state(std::size_t index, const GlobalStateView& global_state) const;
    SpatialDefinition _definition;
    std::vector<std::int64_t> _block_ids;
    std::vector<RegionMesh> _meshes;
    std::vector<std::size_t> _node_offsets, _element_offsets;
    DofMap _dof_map;
    std::vector<std::vector<Quad4RzGeometry>> _region_geometries;
    std::vector<Line2RzBoundaryKernel> _convection_kernels;
    std::vector<ConvectionLoad> _convection_loads;
    std::vector<ConvectionContribution> _convection_contributions;
    std::vector<Line2RzBoundaryKernel> _pressure_kernels;
    std::vector<PressureContribution> _pressure_contributions;
    std::vector<Line2RzBoundaryKernel> _traction_kernels;
    std::vector<TractionContribution> _traction_contributions;
    std::vector<DirichletCondition> _dirichlet_conditions;
    std::vector<ControlledDirichlet> _controlled_dirichlet_conditions;
    std::vector<ControlledScalar> _pressure_loads, _traction_loads;
    std::vector<Line2RzGapHeatKernel> _thermal_kernels;
    std::vector<NodeToLineRzContactKernel> _mechanical_kernels;
    std::vector<ThermalContribution> _thermal_contributions;
    std::vector<MechanicalContribution> _mechanical_contributions;
    mutable std::vector<std::vector<bool>> _projected_thermal_points;
    mutable std::vector<std::vector<bool>> _projected_mechanical_nodes;
    std::vector<std::size_t> _thermal_point_counts;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<ResolvedBoundary> _primary_boundaries, _secondary_boundaries;
    double _load_factor = 1.0, _time = 0.0;
};
} // namespace fuelsim::rz
