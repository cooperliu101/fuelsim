#pragma once
#include "fuelsim/hex8.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "spatial_common.hpp"
#include <cstddef>
#include <vector>

namespace fuelsim::cartesian {
struct ResolvedBoundary final {
    std::size_t region;
    Hex8RegionBoundary boundary;
};

class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);

    const Hex8RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }

    SpatialContributionType contribution_type(std::size_t index) const;
    const Hex8Geometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    void set_load_factor(double load_factor);
    void set_time(double time);

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept {
        return _contact_histories;
    }

    void commit_contact_state(const std::vector<double>& state);
    void restore_contact_state(
        const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories);
    std::vector<CartesianContactNodeSummary> summarize_contact_nodes(
        std::size_t contact_index, const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;

    std::size_t contribution_count() const noexcept { return contribution_ranges().end; }

    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void compute_contribution(std::size_t index, const std::vector<double>& state,
        const std::vector<double>* committed_solution, const Hex8MaterialHistory* committed_material, double time_step,
        std::vector<double>& residual, std::vector<double>* jacobian) const;
    Hex8MaterialHistory transient_update(std::size_t region, std::size_t element, const Hex8LocalValues& state,
        const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material, double time_step) const;
    Hex8LocalValues volume_state(std::size_t index, const std::vector<double>& global_state) const;
    std::array<SymmetricTensor3Values, 8> stress(
        std::size_t region, std::size_t element, const std::vector<double>& state) const;
    double heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const;

  private:
    struct ContributionRanges final {
        std::size_t thermal_begin, mechanical_begin, boundary_begin, end;
    };

    struct BoundaryContribution final {
        SpatialContributionType type;
        std::size_t kernel;
        std::array<std::size_t, 4> nodes;
        Quad4FaceGeometry geometry;
    };

    struct ThermalContribution final {
        std::size_t contact;
        std::array<std::size_t, 8> nodes;
        Quad4ToQuad4HeatGeometry geometry;
        std::size_t integration_point, primary;
    };

    struct MechanicalContribution final {
        std::size_t contact;
        std::array<std::size_t, 8> nodes;
        NodeToQuad4ContactGeometry geometry;
        std::size_t secondary, primary;
    };

    ContributionRanges contribution_ranges() const noexcept;
    ResolvedBoundary resolve_boundary(const UnstructuredHex8Mesh& source_mesh, const std::string& name) const;
    void build_contacts(const UnstructuredHex8Mesh& source_mesh);
    Quad4SurfaceContactLocalDofs contact_dofs(const std::array<std::size_t, 8>& nodes) const;
    Quad4SurfaceContactLocalValues contribution_state(std::size_t index, const std::vector<double>& global_state) const;
    void update_thermal_candidates(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void update_mechanical_candidates(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    bool mark_touched_thermal_points(std::size_t first, std::size_t last) const;
    bool mark_touched_mechanical_nodes(std::size_t first, std::size_t last) const;
    std::size_t thermal_point_index(std::size_t contact, std::size_t point) const noexcept;
    std::size_t mechanical_node_index(std::size_t contact, std::size_t node) const noexcept;
    void refresh_controls();
    std::vector<Hex8RegionMesh> _meshes;
    std::vector<std::vector<Hex8Geometry>> _geometries;
    std::vector<Hex8ThermoelasticData> _kernel_data;
    std::vector<std::size_t> _boundary_definition_indices;
    std::vector<Quad4FaceBoundaryData> _boundary_data;
    std::vector<BoundaryContribution> _boundary_contributions;
    std::vector<GapHeatProperties> _thermal_properties;
    std::vector<NormalContactProperties> _mechanical_properties;
    std::vector<ThermalContribution> _thermal_contributions;
    std::vector<MechanicalContribution> _mechanical_contributions;
    std::vector<std::size_t> _thermal_point_counts, _thermal_contact_offsets, _mechanical_contact_offsets;
    mutable std::vector<unsigned char> _touched_thermal_points, _projected_thermal_candidates,
        _touched_mechanical_nodes, _projected_mechanical_candidates;
    mutable std::vector<double> _thermal_minimum_distance, _mechanical_minimum_distance;
    mutable std::vector<std::size_t> _thermal_selected_primary, _mechanical_selected_primary;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<ResolvedBoundary> _primary_boundaries, _secondary_boundaries;
};
} // namespace fuelsim::cartesian
