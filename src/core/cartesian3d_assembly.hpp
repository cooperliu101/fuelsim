#pragma once
#include "fuelsim/core/cartesian3d_hex20.hpp"
#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/contact.hpp"
#include "fuelsim/core/mesh.hpp"
#include "fuelsim/core/nonlinear_problem.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include "spatial_layout.hpp"
#include <cstddef>
#include <vector>

namespace fuelsim::cartesian {
struct ResolvedBoundary final {
    std::size_t region;
    Hex8RegionBoundary boundary;
};

struct ResolvedHex20Boundary final {
    std::size_t region;
    Hex20RegionBoundary boundary;
};

class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh);

    bool uses_hex20() const noexcept { return _uses_hex20; }

    const Hex8RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }

    const Hex20RegionMesh& hex20_region_mesh(std::size_t index) const { return _hex20_meshes.at(index); }

    SpatialContributionType contribution_type(std::size_t index) const;
    const Hex8Geometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    const Hex20Geometry& hex20_region_element_geometry(std::size_t region_index, std::size_t element_index) const;
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

    std::size_t sparsity_contribution_count() const noexcept;
    std::pair<std::size_t, std::size_t> contribution_partition(
        std::size_t partition, std::size_t partition_count) const;

    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const;
    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const;
    void compute_contribution(std::size_t index, const std::vector<double>& state,
        const std::vector<double>* committed_solution, const CartesianMaterialHistory* committed_material,
        double time_step, std::vector<double>& residual, std::vector<double>* jacobian,
        bool include_thermal_time_term = true) const;
    CartesianMaterialHistory transient_update(std::size_t region, std::size_t element, const Hex8LocalValues& state,
        const Hex8LocalValues& committed_state, const CartesianMaterialHistory& committed_material,
        double time_step) const;
    CartesianMaterialHistory transient_update(std::size_t region, std::size_t element, const Hex20LocalValues& state,
        const Hex20LocalValues& committed_state, const CartesianMaterialHistory& committed_material,
        double time_step) const;
    Hex8LocalValues volume_state(std::size_t index, const std::vector<double>& global_state) const;
    Hex20LocalValues hex20_volume_state(std::size_t index, const std::vector<double>& global_state) const;
    std::array<SymmetricTensor3Values, 8> stress(
        std::size_t region, std::size_t element, const std::vector<double>& state) const;
    std::array<SymmetricTensor3Values, 27> hex20_stress(
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

    struct Hex20BoundaryContribution final {
        SpatialContributionType type;
        std::size_t kernel;
        std::array<std::size_t, 4> temperature_nodes;
        std::array<std::size_t, 8> displacement_nodes;
        Quad8FaceGeometry geometry;
    };

    struct ThermalCandidate final {
        std::size_t contact;
        std::array<std::size_t, 8> nodes;
        Quad4ToQuad4HeatGeometry geometry;
        std::size_t primary;
    };

    struct MechanicalCandidate final {
        std::size_t contact;
        std::array<std::size_t, 8> nodes;
        NodeToQuad4ContactGeometry geometry;
        std::size_t secondary, primary;
    };

    struct Hex20ThermalCandidate final {
        std::size_t contact;
        std::array<std::size_t, 4> secondary_temperature_nodes, primary_temperature_nodes;
        std::array<std::size_t, 8> secondary_displacement_nodes, primary_displacement_nodes;
        Quad8ToQuad8HeatGeometry geometry;
        std::size_t primary;
    };

    struct Hex20MechanicalCandidate final {
        std::size_t contact;
        std::array<std::size_t, 4> secondary_temperature_nodes, primary_temperature_nodes;
        std::array<std::size_t, 8> secondary_displacement_nodes, primary_displacement_nodes;
        NodeToQuad8ContactGeometry node_geometry;
        Quad8ToQuad8MechanicalGeometry surface_geometry;
        bool surface_to_surface;
        std::size_t secondary, primary, secondary_face, secondary_local_point;
    };

    struct PrimaryContactFace final {
        std::array<std::size_t, 4> nodes;
        Quad4FaceCoordinates coordinates;
        CartesianPoint3 parent_centroid;
    };

    struct SecondaryContactFace final {
        std::array<std::size_t, 4> nodes;
        Quad4FaceCoordinates coordinates;
        Quad4FaceGeometry geometry;
        CartesianPoint3 parent_centroid;
    };

    struct SparsityContact final {
        std::array<std::size_t, 8> nodes;
        bool thermal, mechanical;
    };

    struct MechanicalPoint final {
        std::size_t contact, secondary, secondary_face, secondary_local_node;
    };

    struct Hex20PrimaryContactFace final {
        std::array<std::size_t, 4> temperature_nodes;
        std::array<std::size_t, 8> displacement_nodes;
        Quad8FaceCoordinates coordinates;
        CartesianPoint3 parent_centroid;
    };

    struct Hex20SecondaryContactFace final {
        std::array<std::size_t, 4> temperature_nodes;
        std::array<std::size_t, 8> displacement_nodes;
        Quad8FaceCoordinates coordinates;
        Quad8FaceGeometry geometry;
        CartesianPoint3 parent_centroid;
        std::vector<Quad8FaceMechanicalQuadraturePoint> contact_points;
        std::vector<std::size_t> contact_primary_faces;
        std::vector<std::array<double, 8>> contact_primary_shapes, contact_primary_derivatives_xi,
            contact_primary_derivatives_eta;
        std::vector<double> contact_normal_orientations;
    };

    struct Hex20MechanicalPoint final {
        std::size_t contact, secondary, secondary_face, secondary_local_point, reference_primary;
    };

    struct Hex20AveragedConstraint final {
        std::size_t contact, secondary;
        std::vector<std::size_t> nodes, secondary_output_nodes;
        std::vector<double> gap_coefficients, secondary_coefficients;
        CartesianPoint3 normal;
        double reference_gap, area;
    };

    struct Hex20AveragedConstraintValue final {
        double gap, pressure, force;
    };

    ContributionRanges contribution_ranges() const noexcept;
    std::size_t contribution_work(std::size_t index, std::size_t partition_count) const;
    ResolvedBoundary resolve_boundary(const UnstructuredHex8Mesh& source_mesh, const std::string& name) const;
    ResolvedHex20Boundary resolve_boundary(const UnstructuredHex20Mesh& source_mesh, const std::string& name) const;
    void build_contacts(const UnstructuredHex8Mesh& source_mesh);
    void build_hex20_contacts(const UnstructuredHex20Mesh& source_mesh);
    ThermalCandidate thermal_candidate(std::size_t point, std::size_t primary) const;
    MechanicalCandidate mechanical_candidate(std::size_t point, std::size_t primary) const;
    Hex20ThermalCandidate hex20_thermal_candidate(std::size_t point, std::size_t primary) const;
    Hex20MechanicalCandidate hex20_mechanical_candidate(std::size_t point, std::size_t primary) const;
    SparsityContact sparsity_contact(std::size_t index) const;

    struct Hex20SparsityContact final {
        std::array<std::size_t, 4> temperature_nodes, primary_temperature_nodes;
        std::array<std::size_t, 8> displacement_nodes, primary_displacement_nodes;
        bool thermal, mechanical;
    };

    Hex20SparsityContact hex20_sparsity_contact(std::size_t index) const;
    void update_contact_search_trees(const std::vector<double>& state) const;
    Quad4SurfaceContactLocalDofs contact_dofs(const std::array<std::size_t, 8>& nodes) const;
    Quad8SurfaceContactLocalDofs hex20_contact_dofs(const Hex20ThermalCandidate& candidate) const;
    Quad8SurfaceContactLocalDofs hex20_contact_dofs(const Hex20MechanicalCandidate& candidate) const;
    void hex20_averaged_constraint_dofs(
        const Hex20AveragedConstraint& constraint, std::vector<std::size_t>& dofs) const;
    Hex20AveragedConstraintValue hex20_averaged_constraint_value(
        const Hex20AveragedConstraint& constraint, const std::vector<double>& state) const;
    void compute_hex20_averaged_constraint(const Hex20AveragedConstraint& constraint, const std::vector<double>& state,
        std::vector<double>& residual, std::vector<double>* jacobian) const;
    Quad4SurfaceContactLocalValues contribution_state(std::size_t index, const std::vector<double>& global_state) const;
    Quad4SurfaceContactLocalValues contact_state(
        const std::array<std::size_t, 8>& nodes, const std::vector<double>& global_state) const;
    Quad8SurfaceContactLocalValues hex20_contact_state(
        const Hex20ThermalCandidate& candidate, const std::vector<double>& global_state) const;
    Quad8SurfaceContactLocalValues hex20_contact_state(
        const Hex20MechanicalCandidate& candidate, const std::vector<double>& global_state) const;
    void update_thermal_candidates(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void update_mechanical_candidates(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    bool mark_touched_thermal_points(std::size_t first, std::size_t last) const;
    bool mark_touched_mechanical_nodes(std::size_t first, std::size_t last) const;
    std::size_t mechanical_node_index(std::size_t contact, std::size_t node) const noexcept;
    void refresh_controls();
    std::vector<Hex8RegionMesh> _meshes;
    std::vector<Hex20RegionMesh> _hex20_meshes;
    std::vector<std::vector<Hex8Geometry>> _geometries;
    std::vector<std::vector<Hex20Geometry>> _hex20_geometries;
    std::vector<CartesianThermoelasticData> _kernel_data;
    std::vector<std::size_t> _boundary_definition_indices;
    std::vector<Quad4FaceBoundaryData> _boundary_data;
    std::vector<BoundaryContribution> _boundary_contributions;
    std::vector<Hex20BoundaryContribution> _hex20_boundary_contributions;
    std::vector<GapHeatProperties> _thermal_properties;
    std::vector<NormalContactProperties> _mechanical_properties;
    std::vector<std::vector<PrimaryContactFace>> _primary_contact_faces;
    std::vector<std::vector<SecondaryContactFace>> _secondary_contact_faces;
    std::vector<MechanicalPoint> _mechanical_points;
    std::vector<std::vector<Hex20PrimaryContactFace>> _hex20_primary_contact_faces;
    std::vector<std::vector<Hex20SecondaryContactFace>> _hex20_secondary_contact_faces;
    std::vector<Hex20MechanicalPoint> _hex20_mechanical_points;
    std::vector<Hex20AveragedConstraint> _hex20_averaged_constraints;
    std::vector<std::size_t> _thermal_point_counts, _thermal_contact_offsets, _mechanical_contact_offsets,
        _sparsity_contact_offsets;
    mutable std::vector<unsigned char> _touched_thermal_points, _touched_mechanical_nodes;
    mutable std::vector<double> _thermal_minimum_distance, _mechanical_minimum_distance;
    mutable std::vector<std::size_t> _mechanical_selected_primary, _mechanical_cached_primary;
    mutable std::vector<std::size_t> _thermal_active_primary, _thermal_cached_primary, _mechanical_active_primary;
    mutable std::vector<spatial_detail::ContactSearchTree> _contact_search_trees;
    mutable std::vector<spatial_detail::ContactSearchBox> _contact_search_boxes;
    mutable spatial_detail::ContactSearchQuery _contact_search_query;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<ResolvedBoundary> _primary_boundaries, _secondary_boundaries;
    std::vector<ResolvedHex20Boundary> _hex20_primary_boundaries, _hex20_secondary_boundaries;
    bool _uses_hex20 = false;
};
} // namespace fuelsim::cartesian
