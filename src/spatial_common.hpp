#pragma once
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"
#include <array>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace fuelsim::spatial_detail {
inline constexpr std::size_t contact_search_tree_minimum_items = 64;

enum class DofLayout {
    axisymmetric_rz,
    cartesian_3d,
};

struct ConvectionValues final {
    double coefficient, ambient;
};

struct ContactSearchBox final {
    std::array<double, 3> minimum{}, maximum{};
    std::size_t item = 0;
};

struct ContactSearchQueueEntry final {
    std::size_t node = 0;
    double distance_squared = 0.0;
};

struct ContactSearchQuery final {
    std::array<double, 3> point{};
    std::vector<ContactSearchQueueEntry> queue;
};

class ContactSearchTree final {
  public:
    void build(std::vector<ContactSearchBox> boxes);

    bool can_refit(std::size_t box_count) const noexcept { return box_count == _boxes.size() && !_nodes.empty(); }

    void refit(const std::vector<ContactSearchBox>& boxes);
    void begin_query(const std::array<double, 3>& point, ContactSearchQuery& query) const;
    bool next_candidate(ContactSearchQuery& query, double maximum_distance, std::size_t& item) const;

  private:
    struct Node final {
        std::array<double, 3> minimum{}, maximum{};
        std::size_t left = std::numeric_limits<std::size_t>::max();
        std::size_t right = std::numeric_limits<std::size_t>::max();
        std::size_t item = std::numeric_limits<std::size_t>::max();
    };

    std::size_t build_node(std::vector<std::size_t>& indices, std::size_t begin, std::size_t end);
    double distance_squared(std::size_t node, const std::array<double, 3>& point) const;
    void push_node(ContactSearchQuery& query, std::size_t node) const;

    std::vector<ContactSearchBox> _boxes;
    std::vector<Node> _nodes;
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

    std::size_t node_count() const noexcept { return _node_offsets.back(); }

    std::size_t dof_count() const noexcept { return _field_layout.size() * node_count(); }

    const std::vector<FieldDescriptor>& field_layout() const noexcept { return _field_layout; }

    std::size_t dof(Field field, std::size_t node) const;

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
    std::vector<std::size_t> _node_offsets, _element_offsets;
    std::vector<FieldDescriptor> _field_layout;
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
