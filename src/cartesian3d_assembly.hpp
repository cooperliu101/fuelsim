#pragma once
#include "fuelsim/hex8.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "spatial_common.hpp"
#include <cstddef>
#include <vector>
namespace fuelsim::cartesian {
class SpatialAssembly final : public spatial_detail::SpatialLayout {
  public:
    SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh);
    const Hex8RegionMesh& region_mesh(std::size_t index) const { return _meshes.at(index); }
    SpatialContributionType contribution_type(std::size_t index) const;
    const Hex8Geometry& region_element_geometry(std::size_t region_index, std::size_t element_index) const;
    void set_load_factor(double load_factor);
    void set_time(double time);
    std::size_t contribution_count() const noexcept {
        return volume_contribution_count() + _boundary_contributions.size();
    }
    void validate_state(const std::vector<double>& state) const;
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const;
    void compute_contribution(std::size_t index, const std::vector<double>& state,
        const std::vector<double>* committed_solution, const Hex8MaterialHistory* committed_material, double time_step,
        std::vector<double>& residual, std::vector<double>* jacobian) const;
    Hex8MaterialHistory transient_update(std::size_t region, std::size_t element, const Hex8LocalValues& state,
        const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material, double time_step) const;
    Hex8LocalValues volume_state(std::size_t index, const std::vector<double>& global_state) const;
    std::array<SymmetricTensor3Values, 8> stress(
        std::size_t region, std::size_t element, const std::vector<double>& state) const;
    double heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const;

  private:
    struct BoundaryContribution final {
        SpatialContributionType type;
        std::size_t kernel;
        std::array<std::size_t, 4> nodes;
        Quad4FaceGeometry geometry;
    };
    void refresh_controls();
    std::vector<Hex8RegionMesh> _meshes;
    std::vector<std::vector<Hex8Geometry>> _geometries;
    std::vector<Hex8ThermoelasticData> _kernel_data;
    std::vector<std::size_t> _boundary_definition_indices;
    std::vector<Quad4FaceBoundaryData> _boundary_data;
    std::vector<BoundaryContribution> _boundary_contributions;
};
} // namespace fuelsim::cartesian
