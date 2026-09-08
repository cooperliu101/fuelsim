#pragma once
#include "detail/line3_rz_contact.hpp"
#include "fuelsim/core/rz_quad8.hpp"
#include "spatial_layout.hpp"

namespace fuelsim::rz8 {
class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad8Mesh& source);

    const Quad8RegionMesh& region_mesh(std::size_t region) const { return _meshes.at(region); }

    const Quad8RzGeometry& region_element_geometry(std::size_t region, std::size_t element) const {
        return _geometries.at(region).at(element);
    }

    std::size_t contribution_count() const noexcept {
        return volume_contribution_count() + _boundaries.size() + _candidates.size();
    }

    std::size_t sparsity_contribution_count() const noexcept { return contribution_count(); }

    SpatialContributionType contribution_type(std::size_t index) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;

    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
        contribution_dofs(index, dofs);
    }

    void compute_boundary(std::size_t index,
        const std::vector<double>& state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;

    void validate_state(const std::vector<double>& state) const {
        validate_local_state(0, contribution_count(), state);
    }

    void set_load_factor(double value) {
        set_load_factor_value(value);
        refresh_dirichlet_values();
    }

    void set_time(double value) {
        set_time_value(value);
        refresh_dirichlet_values();
    }

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept {
        return _contact_histories;
    }

    void commit_contact_state(const std::vector<double>& state);
    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories);
    bool uses_augmented_contact() const noexcept;
    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>&, std::size_t);
    InterfaceSummary summarize_interface(std::size_t contact, const std::vector<double>& state) const;
    std::vector<ContactNodeSummary> summarize_contact_nodes(std::size_t contact,
        const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact) const;
    // Output-only recovery: positive compression and shear in the primary tangent
    // direction (physical secondary force, opposite the residual traction sign).
    std::vector<std::array<double, 2>> recover_contact_tractions(std::size_t contact,
        const std::vector<ContactNodeSummary>& nodes,
        const std::vector<double>& state) const;

  private:
    struct Boundary final {
        std::size_t definition, region;
        Line3BoundaryElement edge;
    };

    std::vector<Quad8RegionMesh> _meshes;
    std::vector<std::vector<Quad8RzGeometry>> _geometries;
    std::vector<Boundary> _boundaries;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;

    struct ContactBoundary final {
        std::size_t region;
        Quad8RegionBoundary boundary;
    };

    struct Candidate final {
        Line3ContactGeometry geometry;
        std::array<std::size_t, 16> dofs;
        std::size_t contact, point, secondary, primary;
    };

    struct ContactPoint final {
        std::size_t first, last;
    };

    std::vector<Candidate> _candidates;
    std::vector<ContactPoint> _contact_points;
    mutable std::vector<std::size_t> _active_candidates;
    std::vector<ContactBoundary> _primary, _secondary;
    std::vector<GapHeatProperties> _heat;
    std::vector<NormalContactProperties> _mechanical;
    std::vector<double> _committed_contact_solution;
    void build_contacts(const UnstructuredQuad8Mesh& source);
    Line3ContactResult
    candidate_value(std::size_t candidate, const std::vector<double>& state, bool jacobian = false) const;
    void validate_contact_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
};
} // namespace fuelsim::rz8
