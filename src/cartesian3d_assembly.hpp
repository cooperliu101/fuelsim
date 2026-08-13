#pragma once
#include "fuelsim/dof_map.hpp"
#include "fuelsim/hex8.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"
#include <cstddef>
#include <string>
#include <utility>
#include <variant>
#include <vector>
namespace fuelsim::cartesian {
class SpatialAssembly final {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh,
        std::vector<double> heat_capacities = {});
    const SpatialDefinition& definition() const noexcept { return _definition; }
    const DofMap& dof_map() const noexcept { return _dof_map; }
    std::size_t region_count() const noexcept { return _meshes.size(); }
    std::size_t region_index(const std::string& name) const;
    const RegionDefinition& region(std::size_t index) const { return _definition.regions.at(index); }
    const Hex8RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }
    std::size_t region_node_offset(std::size_t index) const;
    std::size_t region_element_count(std::size_t index) const;
    std::size_t region_element_offset(std::size_t index) const;
    std::size_t volume_contribution_count() const noexcept { return _element_offsets.back(); }
    SpatialContributionType contribution_type(std::size_t index) const;
    std::pair<std::size_t, std::size_t> element_location(std::size_t index) const;
    const Hex8Geometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    double region_heat_source(std::size_t region_index) const;
    void set_load_factor(double load_factor);
    double load_factor() const noexcept { return _load_factor; }
    void set_time(double time);
    std::vector<double> initial_state() const;
    std::size_t dof_count() const noexcept { return _dof_map.dof_count(); }
    std::size_t contribution_count() const noexcept {
        return volume_contribution_count() + _boundary_contributions.size();
    }
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept { return _dirichlet_conditions; }
    void validate_state(const std::vector<double>& state) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void contribution_residual(std::size_t index, const std::vector<double>& state,
        const std::vector<double>* committed_solution, double time_step, std::vector<double>& residual) const;
    void contribution_system(std::size_t index, const std::vector<double>& state,
        const std::vector<double>* committed_solution, double time_step, std::vector<double>& residual,
        std::vector<double>& jacobian) const;
    Hex8LocalValues volume_state(std::size_t index, const std::vector<double>& global_state) const;
    std::array<SymmetricTensor3Values, 8> stress(
        std::size_t region, std::size_t element, const std::vector<double>& state) const;
    double heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const;

  private:
    struct ControlledScalar final {
        double base_value;
        bool scale_with_load;
        std::string function;
    };
    struct ControlledDirichlet final {
        std::size_t dof;
        ControlledScalar control;
    };
    struct ConvectionControl final {
        double coefficient, ambient;
        std::string coefficient_function, ambient_function;
    };
    struct BoundaryContribution final {
        SpatialContributionType type;
        std::size_t kernel;
        std::array<std::size_t, 4> nodes;
        Quad4FaceGeometry geometry;
    };
    double function_value(const std::string& name) const;
    double controlled_value(const ControlledScalar& control) const;
    void refresh_controls();
    std::size_t global_node(std::size_t region, std::size_t local_node) const;
    SpatialDefinition _definition;
    std::vector<std::int64_t> _block_ids;
    std::vector<Hex8RegionMesh> _meshes;
    std::vector<std::size_t> _node_offsets, _element_offsets;
    DofMap _dof_map;
    std::vector<std::vector<Hex8Geometry>> _geometries;
    std::vector<Hex8ThermoelasticKernel> _kernels;
    std::vector<DirichletCondition> _dirichlet_conditions;
    std::vector<ControlledDirichlet> _controlled_dirichlet;
    std::vector<std::variant<ControlledScalar, ConvectionControl>> _boundary_controls;
    std::vector<Quad4FaceBoundaryKernel> _boundary_kernels;
    std::vector<BoundaryContribution> _boundary_contributions;
    double _load_factor, _time;
};
} // namespace fuelsim::cartesian
