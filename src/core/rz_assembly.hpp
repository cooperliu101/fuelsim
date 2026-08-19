#pragma once
#include "fuelsim/core/contact.hpp"
#include "fuelsim/core/mesh.hpp"
#include "fuelsim/core/nonlinear_problem.hpp"
#include "fuelsim/core/rz_quad4.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include "spatial_layout.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim::rz {
struct ResolvedBoundary final {
    std::size_t region;
    RegionBoundary boundary;
};

class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh);

    const RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }

    SpatialContributionType contribution_type(std::size_t index) const;
    const Quad4RzGeometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;

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
    void set_time(double time);
    std::vector<ContactNodeSummary> summarize_contact_nodes(
        std::size_t contact_index, const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;

    std::size_t contribution_count() const noexcept { return contribution_ranges().end; }

    std::size_t sparsity_contribution_count() const noexcept;

    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    LocalDofs contribution_dofs(std::size_t index) const;
    LocalDofs sparsity_contribution_dofs(std::size_t index) const;
    LocalResidual compute_contribution(
        std::size_t index, const LocalValues& state, LocalJacobian* jacobian = nullptr) const;
    std::pair<std::size_t, std::array<std::size_t, 2>> edge_parent(
        std::size_t region, const Line2BoundaryElement& edge) const;

  private:
    struct ContributionRanges final {
        std::size_t thermal_begin, mechanical_begin, boundary_begin, end;
    };

    struct ContributionLocation final {
        SpatialContributionType type;
        std::size_t local_index;
    };

    struct BoundaryContribution final {
        SpatialContributionType type;
        std::size_t kernel;
        std::array<std::size_t, 4> nodes;
        Line2RzBoundaryGeometry geometry;
    };

    struct ThermalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        Line2RzHeatPointGeometry geometry;
        std::size_t integration_point, primary;
    };

    struct MechanicalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        NodeToLineRzContactGeometry geometry;
        std::size_t secondary, primary, point;
    };

    struct MechanicalPoint final {
        std::size_t contact, secondary;
    };

    void build_volume_geometries();
    ResolvedBoundary resolve_boundary(const UnstructuredQuad4Mesh& source_mesh, const std::string& name) const;
    void build_boundaries(const UnstructuredQuad4Mesh& source_mesh);
    void refresh_controlled_values();
    LocalDofs local_dofs(const std::array<std::size_t, 4>& nodes) const;
    LocalValues contact_state(const std::array<std::size_t, 4>& nodes, const std::vector<double>& state) const;
    void build_contacts(const UnstructuredQuad4Mesh& source_mesh);
    void initialize_contact_search_workspace();
    void update_contact_search_trees(const std::vector<double>& state) const;
    void update_thermal_candidates(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void update_mechanical_candidates(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void update_large_thermal_candidates(const std::vector<double>& state) const;
    void update_large_mechanical_candidates(const std::vector<double>& state) const;
    bool mark_touched_thermal_points(std::size_t first, std::size_t last) const;
    bool mark_touched_mechanical_nodes(std::size_t first, std::size_t last) const;
    std::size_t thermal_point_index(std::size_t contact, std::size_t point) const noexcept;
    std::size_t mechanical_node_index(std::size_t contact, std::size_t node) const noexcept;
    ContributionRanges contribution_ranges() const noexcept;
    ContributionLocation locate_contribution(std::size_t index) const;
    LocalValues contribution_state(std::size_t index, const std::vector<double>& global_state) const;
    std::vector<RegionMesh> _meshes;
    std::vector<std::vector<Quad4RzGeometry>> _region_geometries;
    std::vector<Line2RzBoundaryData> _boundary_data;
    std::vector<std::size_t> _boundary_definition_indices;
    std::vector<BoundaryContribution> _boundary_contributions;
    std::vector<GapHeatProperties> _thermal_properties;
    std::vector<NormalContactProperties> _mechanical_properties;
    std::vector<ThermalContribution> _thermal_contributions;
    std::vector<MechanicalContribution> _mechanical_contributions;
    std::vector<MechanicalPoint> _mechanical_points;
    std::vector<std::size_t> _thermal_point_counts;
    std::vector<std::size_t> _thermal_contact_offsets;
    std::vector<std::size_t> _mechanical_contact_offsets;
    std::vector<std::size_t> _thermal_candidate_offsets, _mechanical_candidate_offsets;
    mutable std::vector<unsigned char> _touched_thermal_points;
    mutable std::vector<double> _thermal_minimum_distance;
    mutable std::vector<std::size_t> _thermal_active_candidates;
    mutable std::vector<std::size_t> _thermal_cached_primary;
    mutable std::vector<unsigned char> _touched_mechanical_nodes;
    mutable std::vector<unsigned char> _projected_mechanical_candidates;
    mutable std::vector<double> _mechanical_minimum_distance;
    mutable std::vector<std::size_t> _mechanical_selected_primary;
    mutable std::vector<std::size_t> _mechanical_cached_primary;
    mutable std::vector<std::size_t> _mechanical_active_candidates;
    mutable std::vector<spatial_detail::ContactSearchTree> _contact_search_trees;
    mutable std::vector<spatial_detail::ContactSearchBox> _contact_search_boxes;
    mutable spatial_detail::ContactSearchQuery _contact_search_query;
    bool _uses_contact_search_tree = false;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<ResolvedBoundary> _primary_boundaries, _secondary_boundaries;
};
} // namespace fuelsim::rz
