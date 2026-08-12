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

struct CartesianPoint3 final {
    double x;
    double y;
    double z;
};

struct Quad4Element final {
    std::array<std::size_t, 4> nodes;
};

struct Line2BoundaryElement final {
    std::array<std::size_t, 2> nodes;
};

struct Hex8Element final {
    std::array<std::size_t, 8> nodes;
};

struct Quad4FaceElement final {
    std::array<std::size_t, 4> nodes;
    std::size_t parent_element;
    std::size_t local_face;
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

class UnstructuredQuad4Mesh final {
  public:
    UnstructuredQuad4Mesh(std::vector<RzPoint> nodes, std::vector<Quad4Element> elements, std::vector<std::int64_t> element_block_ids, std::vector<ElementBlockInfo> element_blocks, std::vector<NodeSet> node_sets, std::vector<SideSet> side_sets);

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

class UnstructuredHex8Mesh final {
  public:
    UnstructuredHex8Mesh(std::vector<CartesianPoint3> nodes, std::vector<Hex8Element> elements, std::vector<std::int64_t> element_block_ids, std::vector<ElementBlockInfo> element_blocks, std::vector<NodeSet> node_sets, std::vector<SideSet> side_sets);

    const std::vector<CartesianPoint3>& nodes() const noexcept;
    const std::vector<Hex8Element>& elements() const noexcept;
    const std::vector<std::int64_t>& element_block_ids() const noexcept;
    const std::vector<ElementBlockInfo>& element_blocks() const noexcept;
    const std::vector<NodeSet>& node_sets() const noexcept;
    const std::vector<SideSet>& side_sets() const noexcept;

    const ElementBlockInfo& element_block(const std::string& name) const;
    const NodeSet& node_set(const std::string& name) const;
    const SideSet& side_set(const std::string& name) const;
    std::int64_t side_set_block_id(const std::string& name) const;

  private:
    std::vector<CartesianPoint3> _nodes;
    std::vector<Hex8Element> _elements;
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
    static RegionMesh from_unstructured_block(const UnstructuredQuad4Mesh& source, const std::string& block_name);
    static RegionMesh from_unstructured_block(const UnstructuredQuad4Mesh& source, std::int64_t block_id);

    const std::vector<RzPoint>& nodes() const noexcept;
    const std::vector<Quad4Element>& elements() const noexcept;
    const std::vector<std::size_t>& source_node_ids() const noexcept;
    const std::vector<std::size_t>& source_element_ids() const noexcept;
    std::int64_t block_id() const noexcept;

    RegionBoundary map_side_set(const UnstructuredQuad4Mesh& source, const std::string& side_set_name) const;

  private:
    std::int64_t _block_id = -1;
    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;
    std::vector<std::size_t> _source_node_ids;
    std::vector<std::size_t> _source_element_ids;
    std::vector<std::size_t> _source_node_to_local;
    std::vector<std::size_t> _source_element_to_local;
};

struct Hex8RegionBoundary final {
    std::vector<std::size_t> nodes;
    std::vector<Quad4FaceElement> faces;
};

class Hex8RegionMesh final {
  public:
    static Hex8RegionMesh from_unstructured_block(const UnstructuredHex8Mesh& source, const std::string& block_name);
    static Hex8RegionMesh from_unstructured_block(const UnstructuredHex8Mesh& source, std::int64_t block_id);

    const std::vector<CartesianPoint3>& nodes() const noexcept;
    const std::vector<Hex8Element>& elements() const noexcept;
    const std::vector<std::size_t>& source_node_ids() const noexcept;
    const std::vector<std::size_t>& source_element_ids() const noexcept;
    std::int64_t block_id() const noexcept;

    Hex8RegionBoundary map_side_set(const UnstructuredHex8Mesh& source, const std::string& side_set_name) const;

  private:
    std::int64_t _block_id = -1;
    std::vector<CartesianPoint3> _nodes;
    std::vector<Hex8Element> _elements;
    std::vector<std::size_t> _source_node_ids;
    std::vector<std::size_t> _source_element_ids;
    std::vector<std::size_t> _source_node_to_local;
    std::vector<std::size_t> _source_element_to_local;
};

} // namespace fuelsim

#endif
