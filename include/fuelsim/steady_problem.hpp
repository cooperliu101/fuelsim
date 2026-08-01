#ifndef FUELSIM_STEADY_PROBLEM_HPP
#define FUELSIM_STEADY_PROBLEM_HPP

#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {

struct RegionDefinition final {
    std::string name;
    std::string block;
    ThermoelasticProperties material;
    double volumetric_heat_source;
    double initial_temperature;
};

struct ContactDefinition final {
    std::string name;
    std::string primary;
    std::string secondary;
    bool thermal;
    bool mechanical;
    double gap_conductivity;
    double minimum_gap;
    double penalty;
};

enum class BoundaryConditionType {
    dirichlet,
    pressure,
};

struct BoundaryConditionDefinition final {
    std::string name;
    BoundaryConditionType type;
    std::string boundary;
    Field field;
    double value;
};

struct SteadyProblemDefinition final {
    std::vector<RegionDefinition> regions;
    std::vector<ContactDefinition> contacts;
    std::vector<BoundaryConditionDefinition> boundary_conditions;
};

struct ContactNodeSummary final {
    double z;
    bool projected;
    double gap;
    double pressure;
    double tributary_area;
    double tributary_length;
    double contact_force;
};

struct InterfaceSummary final {
    double minimum_gap;
    double maximum_gap;
    double minimum_contact_gap;
    double maximum_contact_pressure;
    double total_heat_rate;
    double total_contact_force;
    std::size_t projected_contact_nodes;
    std::size_t active_contact_nodes;
    double active_contact_length;
};

class SteadyProblem final : public NonlinearProblem {
  public:
    SteadyProblem(SteadyProblemDefinition definition,
                  const UnstructuredQuad4Mesh& source_mesh);

    const SteadyProblemDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;

    std::size_t region_count() const noexcept;
    std::size_t region_index(const std::string& name) const;
    const RegionDefinition& region(std::size_t region_index) const;
    const RegionMesh& region_mesh(std::size_t region_index) const;
    const Quad4RzThermoelasticKernel&
    region_kernel(std::size_t region_index) const;
    std::size_t region_node_offset(std::size_t region_index) const;
    std::size_t region_element_count(std::size_t region_index) const;
    std::size_t region_element_offset(std::size_t region_index) const;
    std::size_t volume_contribution_count() const noexcept;
    std::pair<std::size_t, std::size_t>
    element_location(std::size_t contribution_index) const;
    const Quad4RzGeometry&
    region_element_geometry(std::size_t region_index,
                            std::size_t element_index) const;

    std::size_t contact_count() const noexcept;
    const ContactDefinition& contact(std::size_t contact_index) const;

    void set_load_factor(double load_factor);
    double load_factor() const noexcept;

    std::vector<double> initial_state() const;
    std::vector<ContactNodeSummary>
    summarize_contact_nodes(std::size_t contact_index,
                            const std::vector<double>& state) const;
    InterfaceSummary
    summarize_interface(std::size_t contact_index,
                        const std::vector<double>& state) const;
    void add_external_residual(std::vector<double>& residual) const;

    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    const std::vector<DirichletCondition>&
    dirichlet_conditions() const noexcept override;
    LocalDofs contribution_dofs(std::size_t contribution_index) const override;
    LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const LocalValues& state) const override;
    LocalSystem linearize_contribution(std::size_t contribution_index,
                                       const LocalValues& state) const override;

  protected:
    void add_state_independent_residual(
        std::vector<double>& residual) const override;

  private:
    struct ResolvedBoundary final {
        std::size_t region;
        RegionBoundary boundary;
    };

    struct PressureLoad final {
        std::size_t region;
        RegionBoundary boundary;
        double pressure;
    };

    SteadyProblem(SteadyProblemDefinition definition,
                  const UnstructuredQuad4Mesh& source_mesh,
                  std::vector<std::int64_t> block_ids,
                  std::vector<RegionMesh> meshes);

    static std::vector<std::int64_t>
    resolve_block_ids(const SteadyProblemDefinition& definition,
                      const UnstructuredQuad4Mesh& source_mesh);
    static std::vector<RegionMesh>
    build_meshes(const SteadyProblemDefinition& definition,
                 const UnstructuredQuad4Mesh& source_mesh);

    ResolvedBoundary resolve_boundary(const UnstructuredQuad4Mesh& source_mesh,
                                      const std::string& name) const;
    std::size_t global_node(std::size_t region_index,
                            std::size_t local_node) const;
    void build_volume_geometries();
    void build_contacts(const UnstructuredQuad4Mesh& source_mesh);
    void build_boundary_conditions(const UnstructuredQuad4Mesh& source_mesh);
    void add_pressure_residual(std::vector<double>& residual) const;

    SteadyProblemDefinition _definition;
    std::vector<std::int64_t> _block_ids;
    std::vector<RegionMesh> _meshes;
    std::vector<std::size_t> _node_offsets;
    std::vector<std::size_t> _element_offsets;
    DofMap _dof_map;
    std::vector<Quad4RzThermoelasticKernel> _region_kernels;
    std::vector<std::vector<Quad4RzGeometry>> _region_geometries;

    std::vector<Line2RzGapHeatKernel> _thermal_kernels;
    std::vector<NodeToLineRzContactKernel> _mechanical_kernels;
    std::vector<std::size_t> _thermal_contact_indices;
    std::vector<std::array<std::size_t, 4>> _thermal_nodes;
    std::vector<Line2RzHeatGeometry> _thermal_geometries;
    std::vector<std::size_t> _mechanical_contact_indices;
    std::vector<std::array<std::size_t, 4>> _mechanical_nodes;
    std::vector<NodeToLineRzContactGeometry> _mechanical_geometries;
    std::vector<std::size_t> _mechanical_secondary_indices;
    std::vector<ResolvedBoundary> _primary_boundaries;
    std::vector<ResolvedBoundary> _secondary_boundaries;

    std::vector<DirichletCondition> _dirichlet_conditions;
    std::vector<PressureLoad> _pressure_loads;
    double _load_factor;
};

} // namespace fuelsim

#endif
