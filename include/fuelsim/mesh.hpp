#ifndef FUELSIM_MESH_HPP
#define FUELSIM_MESH_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim {

struct RzPoint final {
    double r;
    double z;
};

struct Quad4Element final {
    std::array<std::size_t, 4> nodes;
};

struct Line2BoundaryElement final {
    std::array<std::size_t, 2> nodes;
};

struct ElementBlockInfo final {
    std::int64_t id;
    std::string name;
};

struct NodeSet final {
    std::int64_t id;
    std::string name;
    std::vector<std::size_t> nodes;
};

struct ElementSide final {
    std::size_t element;
    std::size_t local_side;
};

struct SideSet final {
    std::int64_t id;
    std::string name;
    std::vector<ElementSide> sides;
};

enum class BoundaryId {
    radial_inner,
    radial_outer,
    bottom,
    top,
};

class UnstructuredQuad4Mesh final {
  public:
    UnstructuredQuad4Mesh(std::vector<RzPoint> nodes,
                          std::vector<Quad4Element> elements,
                          std::vector<std::int64_t> element_block_ids,
                          std::vector<ElementBlockInfo> element_blocks,
                          std::vector<NodeSet> node_sets,
                          std::vector<SideSet> side_sets);

    const std::vector<RzPoint>& nodes() const noexcept;
    const std::vector<Quad4Element>& elements() const noexcept;
    const std::vector<std::int64_t>& element_block_ids() const noexcept;
    const std::vector<ElementBlockInfo>& element_blocks() const noexcept;
    const std::vector<NodeSet>& node_sets() const noexcept;
    const std::vector<SideSet>& side_sets() const noexcept;

    const ElementBlockInfo& element_block(const std::string& name) const;
    const NodeSet& node_set(const std::string& name) const;
    const SideSet& side_set(const std::string& name) const;
    std::int64_t side_set_block_id(const std::string& name) const;

  private:
    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;
    std::vector<std::int64_t> _element_block_ids;
    std::vector<ElementBlockInfo> _element_blocks;
    std::vector<NodeSet> _node_sets;
    std::vector<SideSet> _side_sets;
};

enum class RegionBoundaryKind {
    radial_inner,
    radial_outer,
    bottom,
    top,
    general,
};

struct RegionBoundary final {
    RegionBoundaryKind kind;
    std::vector<std::size_t> nodes;
    std::vector<Line2BoundaryElement> elements;
};

class RegionMesh final {
  public:
    static RegionMesh
    from_unstructured_block(const UnstructuredQuad4Mesh& source,
                            const std::string& block_name);
    static RegionMesh
    from_unstructured_block(const UnstructuredQuad4Mesh& source,
                            std::int64_t block_id);

    const std::vector<RzPoint>& nodes() const noexcept;
    const std::vector<Quad4Element>& elements() const noexcept;
    const std::vector<std::size_t>& source_node_ids() const noexcept;
    const std::vector<std::size_t>& source_element_ids() const noexcept;
    std::int64_t block_id() const noexcept;

    RegionBoundary map_side_set(const UnstructuredQuad4Mesh& source,
                                const std::string& side_set_name) const;

  private:
    std::int64_t _block_id = -1;
    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;
    std::vector<std::size_t> _source_node_ids;
    std::vector<std::size_t> _source_element_ids;
    std::vector<std::size_t> _source_node_to_local;
    std::vector<std::size_t> _source_element_to_local;
};

struct RzBoundaryNames final {
    std::string radial_inner;
    std::string radial_outer;
    std::string bottom;
    std::string top;
};

struct StructuredRzBoundary final {
    BoundaryId id;
    std::vector<std::size_t> nodes;
    std::vector<Line2BoundaryElement> elements;
};

class StructuredRzMesh final {
  public:
    static StructuredRzMesh make_annulus(double inner_radius,
                                         double outer_radius, double length,
                                         std::size_t radial_elements,
                                         std::size_t axial_elements);
    static StructuredRzMesh
    from_unstructured_block(const UnstructuredQuad4Mesh& source,
                            const std::string& block_name,
                            const RzBoundaryNames& boundary_names);
    static StructuredRzMesh
    from_unstructured_block(const UnstructuredQuad4Mesh& source,
                            const std::string& block_name);
    static StructuredRzMesh
    from_unstructured_block(const UnstructuredQuad4Mesh& source,
                            std::int64_t block_id,
                            const RzBoundaryNames& boundary_names);

    const std::vector<RzPoint>& nodes() const noexcept;
    const std::vector<Quad4Element>& elements() const noexcept;

    const std::vector<std::size_t>& boundary_nodes(BoundaryId boundary) const;
    const std::vector<Line2BoundaryElement>&
    boundary_elements(BoundaryId boundary) const;

    std::size_t radial_elements() const noexcept;
    std::size_t axial_elements() const noexcept;
    std::size_t node_id(std::size_t radial_index,
                        std::size_t axial_index) const;

    double inner_radius() const noexcept;
    double outer_radius() const noexcept;
    double length() const noexcept;

    StructuredRzBoundary map_side_set(const UnstructuredQuad4Mesh& source,
                                      const std::string& side_set_name,
                                      std::int64_t expected_block_id) const;

  private:
    double _inner_radius = 0.0;
    double _outer_radius = 0.0;
    double _length = 0.0;
    std::size_t _radial_elements = 0;
    std::size_t _axial_elements = 0;

    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;

    std::vector<std::size_t> _radial_inner_nodes;
    std::vector<std::size_t> _radial_outer_nodes;
    std::vector<std::size_t> _bottom_nodes;
    std::vector<std::size_t> _top_nodes;

    std::vector<Line2BoundaryElement> _radial_inner_elements;
    std::vector<Line2BoundaryElement> _radial_outer_elements;
    std::vector<Line2BoundaryElement> _bottom_elements;
    std::vector<Line2BoundaryElement> _top_elements;
};

} // namespace fuelsim

#endif
