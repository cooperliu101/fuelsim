#pragma once
#include "c3d20_types.hpp"
#include "c3d8_types.hpp"
#include "contact_types.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "quad4_face.hpp"
#include "quad8_face.hpp"

#include "core/mesh.hpp"
#include "core/nonlinear_problem.hpp"
#include "core/spatial_definition.hpp"
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

    std::size_t region_material_point_count(std::size_t region) const {
        if (_uses_hex20)
            return this->region(region).hex20_element_formulation == Hex20ElementFormulation::c3d20rt ? 8 : 27;
        return this->region(region).hex8_element_formulation == Hex8ElementFormulation::c3d8rt ? 1 : 8;
    }

    const Hex8RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }

    const Hex20RegionMesh& hex20_region_mesh(std::size_t index) const { return _hex20_meshes.at(index); }

    SpatialContributionType contribution_type(std::size_t index) const;
    const Hex8Geometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    const Hex20Geometry& hex20_region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    void set_load_factor(double load_factor);
    void set_time(double time);
    void set_heat_source_interval(double begin_time, double end_time);
    void set_region_strain_formulation(std::size_t region, StrainFormulation formulation);

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept {
        return _contact_histories;
    }

    double commit_contact_state(const std::vector<double>& state);
    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories);
    std::vector<CartesianContactNodeSummary> summarize_contact_nodes(std::size_t contact_index,
        const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact_index) const;
    InterfaceSummary summarize_interface(std::size_t contact_index, const std::vector<double>& state) const;

    struct FiniteRegionPartitionSummary final {
        std::size_t constraint_count = 0, integration_point_count = 0, active_primary_face_count = 0,
                    cross_face_constraint_count = 0, maximum_owners_per_integration_point = 0;
        bool all_projected = true;
    };

    FiniteRegionPartitionSummary finite_region_partition_summary(std::size_t contact_index,
        const std::vector<double>& state) const;

    std::size_t contribution_count() const noexcept { return contribution_ranges().end; }

    std::size_t sparsity_contribution_count() const noexcept;
    bool jacobian_sparsity_is_state_dependent() const noexcept;
    std::pair<std::size_t, std::size_t> contribution_partition(std::size_t partition,
        std::size_t partition_count) const;

    void validate_state(const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const;
    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const;
    void compute_contribution(std::size_t index,
        const std::vector<double>& state,
        const std::vector<double>* committed_solution,
        const CartesianMaterialHistory* committed_material,
        double time_step,
        std::vector<double>& residual,
        std::vector<double>* jacobian,
        bool include_thermal_time_term = true) const;
    elements::C3d8Result transient_update(std::size_t region,
        std::size_t element,
        const Hex8LocalValues& state,
        const Hex8LocalValues& committed_state,
        const CartesianMaterialHistory& committed_material,
        double time_step) const;
    CartesianMaterialHistory transient_update(std::size_t region,
        std::size_t element,
        const Hex20LocalValues& state,
        const Hex20LocalValues& committed_state,
        const CartesianMaterialHistory& committed_material,
        double time_step) const;
    Hex8LocalValues volume_state(std::size_t index, const std::vector<double>& global_state) const;
    Hex20LocalValues hex20_volume_state(std::size_t index, const std::vector<double>& global_state) const;
    std::array<SymmetricTensor3Values, 8>
    stress(std::size_t region, std::size_t element, const std::vector<double>& state) const;
    std::vector<SymmetricTensor3Values>
    hex20_stress(std::size_t region, std::size_t element, const std::vector<double>& state) const;
    double heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const;
    double mechanical_hourglass_energy(std::size_t region, std::size_t element, const Hex8LocalValues& state) const;

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
        NodeToQuad4ContactGeometry node_geometry;
        Quad4ToQuad4MechanicalGeometry surface_geometry;
        bool surface_to_surface;
        std::size_t secondary, primary, secondary_face, secondary_local_point;
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
        std::array<bool, 4> shared_edges{};
    };

    struct SecondaryContactFace final {
        std::array<std::size_t, 4> nodes;
        Quad4FaceCoordinates coordinates;
        Quad4FaceGeometry geometry;
        CartesianPoint3 parent_centroid;
        std::array<Quad4FaceQuadraturePoint, 4> contact_points;
    };

    struct SparsityContact final {
        std::array<std::size_t, 8> nodes;
        bool thermal, mechanical;
    };

    struct MechanicalPoint final {
        std::size_t contact, secondary, secondary_face, secondary_local_point;
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

    struct Hex20ThermalPoint final {
        std::size_t contact, secondary_face;
        Quad8FaceThermalQuadraturePoint quadrature;
    };

    struct Hex20ThermalPatch final {
        std::size_t contact;
        std::vector<std::size_t> points;
        std::vector<double> fractions;
    };

    struct Hex20MechanicalPoint final {
        std::size_t contact, secondary, secondary_face, secondary_local_point, reference_primary;
    };

    struct AbaqusAveragedConstraint final {
        struct FiniteSlidingSample final {
            struct NormalPoint final {
                double xi, eta, weight;
                std::size_t primary_face;
            };

            std::size_t secondary_face, secondary_local_point;
            double normal_orientation, tangent_orientation;
            std::vector<std::size_t> primary_faces;
            std::vector<std::vector<std::size_t>> primary_transfer_faces;
            std::vector<NormalPoint> normal_points;
        };

        std::size_t contact, secondary, history;
        std::vector<std::size_t> nodes, active_nodes, active_node_indices, secondary_output_nodes;
        std::vector<double> gap_coefficients, secondary_coefficients;
        std::vector<std::array<double, 3>> normal_gap_coefficients, secondary_normal_coefficients;
        std::vector<std::array<double, 3>> tangent_first_coefficients, tangent_second_coefficients,
            traction_first_coefficients, traction_second_coefficients, secondary_tangent_first_coefficients,
            secondary_tangent_second_coefficients;
        std::vector<CartesianPoint3> reference_coordinates;
        std::vector<FiniteSlidingSample> finite_sliding_samples;
        CartesianPoint3 normal, tangent_first, reference_normal, reference_tangent_first;
        double reference_gap, area;
        std::size_t primary_face = 0;
        bool finite_sliding = false, finite_region_normal = false, friction_only = false, projected = true;
        mutable std::vector<double> finite_region_cached_state, finite_region_cached_pressure_derivative;
        mutable double finite_region_cached_gap = 0.0, finite_region_cached_pressure = 0.0,
                       finite_region_cached_force = 0.0;
        mutable bool finite_region_cache_valid = false, finite_region_cached_derivative_valid = false;
    };

    struct AbaqusAveragedConstraintValue final {
        double gap, pressure, force, stick_stiffness, trial_tangential_magnitude, tangential_force,
            friction_dissipation;
        std::array<double, 3> trial_tangential_traction, tangential_traction, tangential_slip, elastic_tangential_slip,
            tangent_first;
        std::array<double, 2> trial_tangential_traction_components, tangential_traction_components;
        bool sliding;
    };

    ContributionRanges contribution_ranges() const noexcept;
    std::size_t contribution_work(std::size_t index, std::size_t partition_count) const;
    ResolvedBoundary resolve_boundary(const UnstructuredHex8Mesh& source_mesh, const std::string& name) const;
    ResolvedHex20Boundary resolve_boundary(const UnstructuredHex20Mesh& source_mesh, const std::string& name) const;
    void build_contacts(const UnstructuredHex8Mesh& source_mesh);
    void build_hex20_contacts(const UnstructuredHex20Mesh& source_mesh);
    ThermalCandidate thermal_candidate(std::size_t point, std::size_t primary) const;
    MechanicalCandidate mechanical_candidate(std::size_t point, std::size_t primary) const;
    NormalContactProperties mechanical_contact_properties(const MechanicalCandidate& candidate) const;
    Hex20ThermalCandidate hex20_thermal_candidate(std::size_t point, std::size_t primary) const;
    void hex20_thermal_patch_dofs(std::size_t patch, std::vector<std::size_t>& dofs, bool all_candidates = false) const;
    std::vector<Quad8HeatPatchSample> hex20_thermal_patch_samples(std::size_t patch,
        const std::vector<std::size_t>& dofs) const;
    void compute_c3d20_finite_constraint_jacobian(const AbaqusAveragedConstraint& constraint,
        const std::vector<double>& state,
        const std::vector<double>& committed_state,
        const ContactPointHistory& history,
        const AbaqusAveragedConstraintValue& value,
        std::vector<double>& jacobian) const;
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
    void averaged_constraint_dofs(const AbaqusAveragedConstraint& constraint, std::vector<std::size_t>& dofs) const;
    std::size_t averaged_sparsity_contribution_count() const noexcept;
    void averaged_sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    AbaqusAveragedConstraintValue averaged_constraint_value(const AbaqusAveragedConstraint& constraint,
        const std::vector<double>& state,
        const std::vector<double>& committed_state,
        const ContactPointHistory& history) const;
    AbaqusAveragedConstraintValue finite_region_normal_value(const AbaqusAveragedConstraint& constraint,
        const std::vector<std::size_t>& local_nodes,
        const std::vector<double>& state,
        std::vector<double>* residual = nullptr,
        std::vector<double>* jacobian = nullptr,
        std::vector<double>* pressure_derivative = nullptr) const;
    double equivalent_normal_pressure(const AbaqusAveragedConstraint& constraint,
        const std::vector<double>& state,
        std::vector<double>* derivative = nullptr,
        const std::vector<double>* friction_area_derivative = nullptr) const;
    void compute_averaged_constraint(const AbaqusAveragedConstraint& constraint,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const;
    void compute_averaged_friction_geometry(const AbaqusAveragedConstraint& constraint,
        const std::vector<double>& state,
        const std::array<double, 2>& traction,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const;
    void compute_averaged_friction_traction_derivatives(const AbaqusAveragedConstraint& constraint,
        const std::vector<double>& state,
        const std::vector<double>& committed_state,
        const ContactPointHistory& history,
        const AbaqusAveragedConstraintValue& value,
        std::array<std::vector<double>, 2>& derivatives) const;
    void refresh_hex20_finite_averaged_constraints(const std::vector<double>& state) const;
    void refresh_finite_averaged_constraints(const std::vector<double>& state) const;
    bool summarize_averaged_contact(std::size_t contact,
        const std::vector<double>& state,
        std::vector<CartesianContactNodeSummary>& summaries) const;
    Quad4SurfaceContactLocalValues contribution_state(std::size_t index, const std::vector<double>& global_state) const;
    Quad4SurfaceContactLocalValues contact_state(const std::array<std::size_t, 8>& nodes,
        const std::vector<double>& global_state) const;
    Quad8SurfaceContactLocalValues hex20_contact_state(const Hex20ThermalCandidate& candidate,
        const std::vector<double>& global_state) const;
    Quad8SurfaceContactLocalValues hex20_contact_state(const Hex20MechanicalCandidate& candidate,
        const std::vector<double>& global_state) const;
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
    std::vector<CartesianRegionData> _kernel_data;
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
    std::vector<Hex20ThermalPoint> _hex20_thermal_points;
    std::vector<Hex20ThermalPatch> _hex20_thermal_patches;
    mutable std::vector<AbaqusAveragedConstraint> _abaqus_averaged_constraints;
    std::vector<std::size_t> _thermal_point_counts, _thermal_contact_offsets, _mechanical_contact_offsets,
        _sparsity_contact_offsets;
    mutable std::vector<unsigned char> _touched_thermal_points, _touched_mechanical_nodes;
    mutable std::vector<double> _thermal_minimum_distance, _mechanical_minimum_distance;
    mutable std::vector<std::size_t> _mechanical_selected_primary, _mechanical_cached_primary;
    mutable std::vector<std::size_t> _thermal_active_primary, _thermal_cached_primary, _mechanical_active_primary;
    mutable std::vector<spatial_detail::ContactSearchTree> _contact_search_trees;
    mutable std::vector<spatial_detail::ContactSearchBox> _contact_search_boxes;
    mutable spatial_detail::ContactSearchQuery _contact_search_query;
    mutable std::vector<double> _fully_validated_contact_state;
    mutable bool _fully_validated_contact_state_current = false;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<ResolvedBoundary> _primary_boundaries, _secondary_boundaries;
    std::vector<ResolvedHex20Boundary> _hex20_primary_boundaries, _hex20_secondary_boundaries;
    bool _uses_hex20 = false;
};
} // namespace fuelsim::cartesian
