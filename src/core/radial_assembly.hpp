#pragma once
#include "cax2t_gps.hpp"
#include "ring_gps.hpp"
#include "spatial_layout.hpp"

namespace fuelsim::radial {
class RegionMesh final : public RegionMeshMapping {
  public:
    RegionMesh(const UnstructuredBar2Mesh& source, std::int64_t block);

    const std::vector<RzPoint>& nodes() const noexcept { return _nodes; }

    // Radial connectivity is local; axial connectivity retains source-node IDs.
    const std::vector<Bar2Element>& elements() const noexcept { return _elements; }

  private:
    std::vector<RzPoint> _nodes;
    std::vector<Bar2Element> _elements;
};

struct LocalContribution final {
    std::vector<double> residual, jacobian;
    elements::Cax2tGpsResult volume;
    elements::RingGpsResult contact;
    double convection_heat_rate = 0.0, surface_heat_input_rate = 0.0;
};

class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredBar2Mesh& source);

    const UnstructuredBar2Mesh& source_mesh() const noexcept { return _source; }

    const RegionMesh& region_mesh(std::size_t region) const { return _meshes.at(region); }

    const std::vector<std::size_t>& region_source_element_ids(std::size_t region) const {
        return region_mesh(region).source_element_ids();
    }

    const std::vector<std::size_t>& region_source_node_ids(std::size_t region) const {
        return region_mesh(region).source_node_ids();
    }

    const Cax2tGpsGeometry& region_element_geometry(std::size_t region, std::size_t element) const {
        return _geometries.at(region).at(element);
    }

    std::array<std::size_t, 2> region_axial_dofs(std::size_t region, std::size_t element) const;
    std::size_t axial_dof(std::size_t source_node) const;
    std::size_t radial_node(std::size_t source_node) const;
    std::vector<double> reference_state() const;
    Cax2tGpsLocalValues volume_state(std::size_t index, const std::vector<double>& global) const;

    std::size_t contribution_count() const noexcept {
        return volume_contribution_count() + _contacts.size() + _boundaries.size();
    }

    std::size_t sparsity_contribution_count() const noexcept { return contribution_count(); }

    bool jacobian_sparsity_is_state_dependent() const noexcept { return false; }

    bool contribution_metadata_is_fixed() const noexcept { return true; }

    std::pair<std::size_t, std::size_t> contribution_partition(std::size_t rank, std::size_t ranks) const;
    SpatialContributionType contribution_type(std::size_t index) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;

    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
        contribution_dofs(index, dofs);
    }

    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;

    void validate_state(const std::vector<double>& state) const {
        validate_local_state(0, contribution_count(), state);
    }

    void compute_contribution(std::size_t index,
        const std::vector<double>& local_state,
        const std::vector<double>& committed_global,
        const Cax2tGpsMaterialHistory* history,
        double dt,
        double time,
        bool include_thermal,
        bool linearize,
        LocalContribution& result) const;
    elements::Cax2tGpsResult evaluate_volume(std::size_t region,
        std::size_t element,
        const std::vector<double>& global,
        const std::vector<double>& committed,
        const Cax2tGpsMaterialHistory* history,
        double dt,
        double time,
        bool include_thermal,
        bool linearize = false) const;
    void compute_boundary(std::size_t index,
        const std::vector<double>& local_state,
        std::vector<double>& residual,
        std::vector<double>* jacobian) const;
    void set_load_factor(double value);
    void set_time(double value);

    double time() const noexcept { return _time; }

    void set_heat_source_interval(double begin, double end);

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept {
        return _contact_histories;
    }

    double commit_contact_state(const std::vector<double>& state);
    double prepare_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>>& histories);

    void publish_contact_state(std::vector<double>& solution,
        std::vector<std::vector<ContactPointHistory>>& histories) noexcept {
        _contact_histories.swap(histories);
        _committed_contact_solution.swap(solution);
    }

    void restore_contact_state(const std::vector<double>& state,
        std::vector<std::vector<ContactPointHistory>> histories);

    bool uses_augmented_contact() const noexcept { return false; }

    AugmentedContactUpdate update_augmented_contact_multipliers(const std::vector<double>&, std::size_t) { return {}; }

    InterfaceSummary summarize_interface(std::size_t contact, const std::vector<double>& state) const;
    std::vector<ContactNodeSummary> summarize_contact_nodes(std::size_t contact,
        const std::vector<double>& state) const;
    std::vector<std::size_t> contact_secondary_source_nodes(std::size_t contact) const;
    std::vector<std::pair<std::size_t, std::size_t>> contact_source_elements(std::size_t contact) const;
    elements::RingGpsResult
    evaluate_contact(std::size_t contact, std::size_t pair, const std::vector<double>& state) const;

  private:
    struct Contact final {
        RingGpsGeometry geometry;
        std::array<std::size_t, ring_gps_local_dof_count> dofs;
        std::size_t definition, history_offset, primary_source, secondary_source, primary_element, secondary_element;
        StrainFormulation strain_formulation;
    };

    struct Boundary final {
        std::size_t definition, element, side, axial_source;
    };

    UnstructuredBar2Mesh _source;
    std::vector<RegionMesh> _meshes;
    std::vector<IsotropicThermoelasticMaterial> _materials;
    std::vector<std::vector<Cax2tGpsGeometry>> _geometries;
    std::vector<std::size_t> _source_to_radial, _source_to_axial, _source_element_to_region;
    std::vector<double> _heat_sources;
    std::vector<Contact> _contacts;
    std::vector<Boundary> _boundaries;
    std::vector<GapHeatProperties> _heat;
    std::vector<NormalContactProperties> _normal;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    void refresh_heat_sources();
    void build_boundaries();
    void build_contacts();
    std::vector<ElementSide> boundary_sides(const std::string& name) const;
    elements::RingGpsResult
    evaluate_contact(std::size_t index, const RingGpsLocalValues& state, bool linearize, bool history = true) const;
    elements::RingGpsResult
    evaluate_global_contact(std::size_t index, const std::vector<double>& state, bool history = true) const;
    void compute_surface(std::size_t index,
        const std::vector<double>& state,
        bool linearize,
        LocalContribution& result) const;
};
} // namespace fuelsim::radial
