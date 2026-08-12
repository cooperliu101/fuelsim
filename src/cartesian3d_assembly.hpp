#ifndef FUELSIM_CARTESIAN3D_ASSEMBLY_HPP
#define FUELSIM_CARTESIAN3D_ASSEMBLY_HPP

#include "fuelsim/hex8.hpp"
#include "fuelsim/hex8_dof_map.hpp"
#include "fuelsim/hex8_thermoelastic.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace cartesian3d {

class SpatialAssembly final {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);

    const SpatialDefinition& definition() const noexcept;
    const Hex8DofMap& dof_map() const noexcept;
    std::size_t region_count() const noexcept;
    std::size_t region_index(const std::string& name) const;
    const RegionDefinition& region(std::size_t index) const;
    const Hex8RegionMesh& region_mesh(std::size_t index) const;
    std::size_t region_node_offset(std::size_t index) const;
    std::size_t region_element_count(std::size_t index) const;
    std::size_t region_element_offset(std::size_t index) const;
    std::size_t volume_contribution_count() const noexcept;
    SpatialContributionType contribution_type(std::size_t contribution_index) const;
    std::pair<std::size_t, std::size_t> element_location(std::size_t contribution_index) const;
    const Hex8Geometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    double region_heat_source(std::size_t region_index) const;

    void set_load_factor(double load_factor);
    double load_factor() const noexcept;
    void set_time(double time);
    std::vector<double> initial_state() const;

    std::size_t dof_count() const noexcept;
    std::size_t contribution_count() const noexcept;
    const std::vector<DirichletCondition>& dirichlet_conditions() const noexcept;
    void validate_state(const std::vector<double>& state) const;
    std::size_t contribution_dof_count(std::size_t contribution_index) const;
    void fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const;
    Quad4FaceLocalResidual boundary_residual(std::size_t contribution_index, const Quad4FaceLocalValues& state) const;
    Quad4FaceLocalSystem boundary_system(std::size_t contribution_index, const Quad4FaceLocalValues& state) const;

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
        double coefficient;
        double ambient;
        std::string coefficient_function;
        std::string ambient_function;
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
    std::vector<std::size_t> _node_offsets;
    std::vector<std::size_t> _element_offsets;
    Hex8DofMap _dof_map;
    std::vector<std::vector<Hex8Geometry>> _geometries;
    std::vector<DirichletCondition> _dirichlet_conditions;
    std::vector<ControlledDirichlet> _controlled_dirichlet;
    std::vector<ControlledScalar> _pressure_controls;
    std::vector<ControlledScalar> _traction_controls;
    std::vector<ConvectionControl> _convection_controls;
    std::vector<Quad4FacePressureKernel> _pressure_kernels;
    std::vector<Quad4FaceTractionKernel> _traction_kernels;
    std::vector<Quad4FaceConvectionKernel> _convection_kernels;
    std::vector<BoundaryContribution> _boundary_contributions;
    double _load_factor;
    double _time;
};

} // namespace cartesian3d
} // namespace fuelsim

#endif
