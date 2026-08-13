#pragma once
#include "fuelsim/dof_map.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>
namespace fuelsim::spatial_detail {
struct ConvectionValues final {
    double coefficient, ambient;
};
std::vector<std::int64_t> resolve_block_ids(const SpatialDefinition& definition,
    const UnstructuredMeshMetadata& source_mesh, bool allow_contacts, bool allow_finite_strain);
void validate_dirichlet_conditions(
    std::vector<DirichletCondition>& conditions, const char* conflict_message, const char* duplicate_message);
double function_value(const SpatialDefinition& definition, double time, const std::string& name);
double controlled_value(const SpatialDefinition& definition, double time, double load_factor, double value,
    bool scale_with_load, const std::string& function);
class SpatialLayout {
  public:
    const SpatialDefinition& definition() const noexcept { return _definition; }
    const DofMap& dof_map() const noexcept { return *_dof_map; }
    std::size_t region_count() const noexcept { return _definition.regions.size(); }
    const RegionDefinition& region(std::size_t index) const { return _definition.regions.at(index); }
    std::size_t region_node_offset(std::size_t index) const;
    std::size_t region_element_count(std::size_t index) const;
    std::size_t region_element_offset(std::size_t index) const;
    std::size_t volume_contribution_count() const noexcept { return _element_offsets.back(); }
    std::pair<std::size_t, std::size_t> element_location(std::size_t index) const;
    double region_heat_source(std::size_t index) const;
    double load_factor() const noexcept { return _load_factor; }
    std::vector<double> initial_state() const;
    std::size_t dof_count() const noexcept { return dof_map().dof_count(); }
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept { return _dirichlet_conditions; }
    std::size_t global_node(std::size_t region, std::size_t local_node) const;

  protected:
    SpatialLayout(SpatialDefinition definition, std::vector<std::int64_t> block_ids, DofLayout layout);
    void initialize_counts(const std::vector<std::size_t>& node_counts, const std::vector<std::size_t>& element_counts);
    void add_dirichlet(std::size_t dof, std::size_t boundary_index);
    void set_load_factor_value(double value);
    void set_time_value(double value);
    void refresh_dirichlet_values();
    ConvectionValues convection_values(const BoundaryConditionDefinition& boundary) const;
    SpatialDefinition _definition;
    std::vector<std::int64_t> _block_ids;
    std::optional<DofMap> _dof_map;
    std::vector<std::size_t> _node_offsets, _element_offsets;
    std::vector<DirichletCondition> _dirichlet_conditions;
    double _load_factor = 0.0, _time = 0.0;

  private:
    struct ControlledDirichlet final {
        std::size_t dof, boundary_index;
    };
    DofLayout _layout;
    std::vector<ControlledDirichlet> _controlled_dirichlet_conditions;
};
} // namespace fuelsim::spatial_detail
