#pragma once
#include "core/plane_mesh.hpp"
#include "line3_plane.hpp"
#include "spatial_layout.hpp"

namespace fuelsim::plane {
class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredPlaneQuad8Mesh& mesh);

    std::size_t contribution_count() const noexcept { return contact_offset() + _candidates.size(); }

    std::size_t contact_offset() const noexcept { return volume_contribution_count() + _boundaries.size(); }

    SpatialContributionType contribution_type(std::size_t index) const;
    elements::Cpeg8Result compute_boundary(std::size_t index, const std::vector<double>& local, bool jacobian) const;
    elements::Line3PlaneContactResult
    compute_contact(std::size_t index, const std::vector<double>& local, bool jacobian) const;
    InterfaceSummary summarize_interface(std::size_t contact, const std::vector<double>& state) const;
    std::vector<elements::Line3PlaneContactResult> contact_points(std::size_t contact,
        const std::vector<double>& state) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void validate_state(const std::vector<double>& state) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;

    void set_load_factor(double value) {
        set_load_factor_value(value);
        refresh_constraints();
    }

    void set_time(double value) {
        set_time_value(value);
        refresh_constraints();
    }

    const elements::Cpeg8Geometry& geometry(std::size_t r, std::size_t e) const { return _geometry.at(r).at(e); }

    const std::vector<std::size_t>& source_nodes(std::size_t r) const { return _source_nodes.at(r); }

    std::size_t section_count() const noexcept { return _section_origins.size(); }

    const std::array<double, 2>& section_origin(std::size_t section) const { return _section_origins.at(section); }

    const std::vector<std::size_t>& source_elements(std::size_t r) const { return _source_elements.at(r); }

    const std::vector<bool>& thermal_nodes(std::size_t r) const { return _thermal_nodes.at(r); }

    std::size_t section_dof(std::size_t component, std::size_t control) const;

    std::pair<std::size_t, std::size_t> boundary_identity(std::size_t index) const {
        const auto& boundary = _boundaries.at(index - volume_contribution_count());
        return {boundary.definition, boundary.side};
    }

  private:
    std::vector<std::vector<std::size_t>> _source_nodes, _source_elements;
    std::vector<std::vector<bool>> _thermal_nodes;
    std::vector<std::size_t> _region_sections;
    std::vector<std::array<double, 2>> _section_origins;
    void refresh_constraints();
    std::vector<std::vector<elements::Cpeg8Geometry>> _geometry;
    std::vector<std::array<std::size_t, 23>> _dofs;

    struct Boundary final {
        std::size_t definition, region, element, side;
    };

    std::vector<Boundary> _boundaries;

    struct ContactSide final {
        std::size_t region, element, side;
    };

    struct Candidate final {
        ContactSide secondary, primary;
        std::array<std::size_t, 22> dofs;
        std::size_t contact;
        double coordinate, weight, penalty;
        std::size_t segment;
    };

    std::vector<Candidate> _candidates;
    std::vector<std::pair<std::size_t, std::size_t>> _constraints;
    void build_contacts(const UnstructuredPlaneQuad8Mesh& mesh);
    void validate_contacts(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    elements::Line3PlaneContactResult
    evaluate_candidate(std::size_t candidate, const elements::Line3PlaneValues& state, bool jacobian) const;
};
} // namespace fuelsim::plane
