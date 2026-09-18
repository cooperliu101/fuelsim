#pragma once
#include "spatial_layout.hpp"
#include "thermal_interface.hpp"
#include "thermal_types.hpp"

namespace fuelsim::thermal {
class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& mesh);
    SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad8Mesh& mesh);
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& mesh);
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex20Mesh& mesh);

    void set_time(double value) {
        set_time_value(value);
        refresh_dirichlet_values();
    }

    void set_load_factor(double value) {
        set_load_factor_value(value);
        refresh_dirichlet_values();
    }

    std::size_t contribution_count() const noexcept { return _contributions.size(); }

    SpatialContributionType contribution_type(std::size_t index) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void validate_state(const std::vector<double>& state) const;
    void validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const;
    std::vector<std::size_t> required_state_dofs(std::size_t first, std::size_t last) const;
    elements::ThermalResult evaluate(std::size_t index,
        const std::vector<double>& state,
        const std::vector<double>& previous,
        double step,
        double begin,
        bool jacobian) const;

    const std::vector<std::size_t>& source_nodes(std::size_t region) const { return _source_nodes.at(region); }

    const std::vector<CartesianPoint3>& coordinates() const noexcept { return _coordinates; }

    const std::vector<std::vector<std::size_t>>& connectivity() const noexcept { return _connectivity; }

    const elements::ThermalGeometry& geometry(std::size_t index) const { return _contributions.at(index).geometry; }

    std::size_t source_element(std::size_t index) const { return _source_elements.at(index); }

    InterfaceSummary summarize_interface(std::size_t contact, const std::vector<double>& state) const;

    bool axisymmetric() const noexcept { return _axisymmetric; }

  private:
    struct Contribution final {
        std::vector<std::size_t> dofs;
        elements::ThermalGeometry geometry;
        elements::ThermalInterfacePoint interface;
        std::size_t region = 0, boundary = 0;
        SpatialContributionType type = SpatialContributionType::volume;
    };

    void initialize(const UnstructuredMeshMetadata& mesh, ThermalElement topology);
    bool _axisymmetric;
    std::vector<CartesianPoint3> _coordinates;
    std::vector<std::vector<std::size_t>> _connectivity, _source_nodes;
    std::vector<Contribution> _contributions;
    std::vector<std::size_t> _source_elements;
};
} // namespace fuelsim::thermal
