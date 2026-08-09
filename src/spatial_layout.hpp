#ifndef FUELSIM_SPATIAL_LAYOUT_HPP
#define FUELSIM_SPATIAL_LAYOUT_HPP

#include "fuelsim/dof_map.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz.hpp"
#include "fuelsim/spatial_definition.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {

class SpatialAssembly;

class SpatialLayout final {
  public:
    struct ResolvedBoundary final {
        std::size_t region;
        RegionBoundary boundary;
    };

    SpatialLayout(SpatialDefinition definition,
                  std::vector<std::int64_t> block_ids,
                  std::vector<RegionMesh> meshes);

    const SpatialDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;
    std::size_t region_count() const noexcept;
    const RegionDefinition& region(std::size_t index) const;
    const RegionMesh& region_mesh(std::size_t index) const;
    std::size_t region_node_offset(std::size_t index) const;
    std::size_t region_element_offset(std::size_t index) const;
    std::size_t volume_contribution_count() const noexcept;
    const Quad4RzGeometry& element_geometry(std::size_t region,
                                            std::size_t element) const;
    ResolvedBoundary resolve_boundary(
        const UnstructuredQuad4Mesh& source_mesh,
        const std::string& name) const;
    std::size_t global_node(std::size_t region,
                            std::size_t local_node) const;
    std::pair<std::size_t, std::array<std::size_t, 2>>
    edge_parent(std::size_t region,
                const Line2BoundaryElement& edge) const;

  private:
    friend class SpatialAssembly;

    void build_volume_geometries();

    SpatialDefinition _definition;
    std::vector<std::int64_t> _block_ids;
    std::vector<RegionMesh> _meshes;
    std::vector<std::size_t> _node_offsets;
    std::vector<std::size_t> _element_offsets;
    DofMap _dof_map;
    std::vector<std::vector<Quad4RzGeometry>> _region_geometries;
};

} // namespace fuelsim

#endif
