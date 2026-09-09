#pragma once
#include "fuelsim/elements/coordinates.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim {
struct Quad4Element final {
    std::array<std::size_t, 4> nodes;
};

struct Quad8Element final {
    std::array<std::size_t, 8> nodes;
};

struct Line3BoundaryElement final {
    std::array<std::size_t, 3> nodes;
    std::size_t parent_element, local_side;
};

struct Line2BoundaryElement final {
    std::array<std::size_t, 2> nodes;
};

struct Hex8Element final {
    std::array<std::size_t, 8> nodes;
};

struct Hex20Element final {
    std::array<std::size_t, 20> nodes;
};

struct Quad4FaceElement final {
    std::array<std::size_t, 4> nodes;
    std::size_t parent_element, local_face;
};

struct Quad8FaceElement final {
    std::array<std::size_t, 8> nodes;
    std::size_t parent_element, local_face;
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
    std::size_t element, local_side;
};

struct SideSet final {
    std::int64_t id;
    std::string name;
    std::vector<ElementSide> sides;
};

class UnstructuredMeshMetadata {
  public:
    const std::vector<std::int64_t>& element_block_ids() const noexcept { return _element_block_ids; }

    const std::vector<ElementBlockInfo>& element_blocks() const noexcept { return _element_blocks; }

    const std::vector<NodeSet>& node_sets() const noexcept { return _node_sets; }

    const std::vector<SideSet>& side_sets() const noexcept { return _side_sets; }

    const ElementBlockInfo& element_block(const std::string& name) const;
    const SideSet& side_set(const std::string& name) const;
    std::int64_t side_set_block_id(const std::string& name) const;

  protected:
    UnstructuredMeshMetadata(std::size_t node_count,
        std::size_t element_count,
        std::size_t sides_per_element,
        std::vector<std::int64_t> element_block_ids,
        std::vector<ElementBlockInfo> element_blocks,
        std::vector<NodeSet> node_sets,
        std::vector<SideSet> side_sets,
        const std::string& geometry_name);

  private:
    std::vector<std::int64_t> _element_block_ids;
    std::vector<ElementBlockInfo> _element_blocks;
    std::vector<NodeSet> _node_sets;
    std::vector<SideSet> _side_sets;
};

class UnstructuredQuad4Mesh final : public UnstructuredMeshMetadata {
  public:
    UnstructuredQuad4Mesh(std::vector<RzPoint> nodes,
        std::vector<Quad4Element> elements,
        std::vector<std::int64_t> element_block_ids,
        std::vector<ElementBlockInfo> element_blocks,
        std::vector<NodeSet> node_sets,
        std::vector<SideSet> side_sets);

    const std::vector<RzPoint>& nodes() const noexcept { return _nodes; }

    const std::vector<Quad4Element>& elements() const noexcept { return _elements; }

  private:
    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;
};

class UnstructuredQuad8Mesh final : public UnstructuredMeshMetadata {
  public:
    UnstructuredQuad8Mesh(std::vector<RzPoint> nodes,
        std::vector<Quad8Element> elements,
        std::vector<std::int64_t> element_block_ids,
        std::vector<ElementBlockInfo> element_blocks,
        std::vector<NodeSet> node_sets,
        std::vector<SideSet> side_sets);

    const std::vector<RzPoint>& nodes() const noexcept { return _nodes; }

    const std::vector<Quad8Element>& elements() const noexcept { return _elements; }

  private:
    std::vector<RzPoint> _nodes;
    std::vector<Quad8Element> _elements;
};

class UnstructuredHex8Mesh final : public UnstructuredMeshMetadata {
  public:
    UnstructuredHex8Mesh(std::vector<CartesianPoint3> nodes,
        std::vector<Hex8Element> elements,
        std::vector<std::int64_t> element_block_ids,
        std::vector<ElementBlockInfo> element_blocks,
        std::vector<NodeSet> node_sets,
        std::vector<SideSet> side_sets);

    const std::vector<CartesianPoint3>& nodes() const noexcept { return _nodes; }

    const std::vector<Hex8Element>& elements() const noexcept { return _elements; }

  private:
    std::vector<CartesianPoint3> _nodes;
    std::vector<Hex8Element> _elements;
};

class UnstructuredHex20Mesh final : public UnstructuredMeshMetadata {
  public:
    UnstructuredHex20Mesh(std::vector<CartesianPoint3> nodes,
        std::vector<Hex20Element> elements,
        std::vector<std::int64_t> element_block_ids,
        std::vector<ElementBlockInfo> element_blocks,
        std::vector<NodeSet> node_sets,
        std::vector<SideSet> side_sets);

    const std::vector<CartesianPoint3>& nodes() const noexcept { return _nodes; }

    const std::vector<Hex20Element>& elements() const noexcept { return _elements; }

  private:
    std::vector<CartesianPoint3> _nodes;
    std::vector<Hex20Element> _elements;
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

class RegionMeshMapping {
  public:
    const std::vector<std::size_t>& source_node_ids() const noexcept { return _source_node_ids; }

    const std::vector<std::size_t>& source_element_ids() const noexcept { return _source_element_ids; }

  protected:
    RegionMeshMapping(const UnstructuredMeshMetadata& source, std::size_t node_count, std::int64_t block_id);
    void select_nodes(const std::vector<bool>& used_nodes);
    static constexpr std::size_t invalid_index = static_cast<std::size_t>(-1);
    std::int64_t _block_id;
    std::vector<std::size_t> _source_node_ids, _source_element_ids;
    std::vector<std::size_t> _source_node_to_local, _source_element_to_local;
};

class RegionMesh final : public RegionMeshMapping {
  public:
    static RegionMesh from_unstructured_block(const UnstructuredQuad4Mesh& source, std::int64_t block_id);

    const std::vector<RzPoint>& nodes() const noexcept { return _nodes; }

    const std::vector<Quad4Element>& elements() const noexcept { return _elements; }

    RegionBoundary map_side_set(const UnstructuredQuad4Mesh& source, const std::string& side_set_name) const;

  private:
    RegionMesh(const UnstructuredQuad4Mesh& source, std::int64_t block_id);
    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;
};

struct Quad8RegionBoundary final {
    std::vector<std::size_t> temperature_nodes, displacement_nodes;
    std::vector<Line3BoundaryElement> elements;
};

class Quad8RegionMesh final : public RegionMeshMapping {
  public:
    static Quad8RegionMesh from_unstructured_block(const UnstructuredQuad8Mesh& source, std::int64_t block_id);

    const std::vector<RzPoint>& nodes() const noexcept { return _nodes; }

    const std::vector<Quad8Element>& elements() const noexcept { return _elements; }

    const std::vector<bool>& temperature_nodes() const noexcept { return _temperature_nodes; }

    Quad8RegionBoundary map_side_set(const UnstructuredQuad8Mesh& source, const std::string& name) const;

  private:
    Quad8RegionMesh(const UnstructuredQuad8Mesh& source, std::int64_t block_id);
    std::vector<RzPoint> _nodes;
    std::vector<Quad8Element> _elements;
    std::vector<bool> _temperature_nodes;
};

struct Hex8RegionBoundary final {
    std::vector<std::size_t> nodes;
    std::vector<Quad4FaceElement> faces;
};

class Hex8RegionMesh final : public RegionMeshMapping {
  public:
    static Hex8RegionMesh from_unstructured_block(const UnstructuredHex8Mesh& source, std::int64_t block_id);

    const std::vector<CartesianPoint3>& nodes() const noexcept { return _nodes; }

    const std::vector<Hex8Element>& elements() const noexcept { return _elements; }

    Hex8RegionBoundary map_side_set(const UnstructuredHex8Mesh& source, const std::string& side_set_name) const;

  private:
    Hex8RegionMesh(const UnstructuredHex8Mesh& source, std::int64_t block_id);
    std::vector<CartesianPoint3> _nodes;
    std::vector<Hex8Element> _elements;
};

struct Hex20RegionBoundary final {
    std::vector<std::size_t> temperature_nodes, displacement_nodes;
    std::vector<Quad8FaceElement> faces;
};

class Hex20RegionMesh final : public RegionMeshMapping {
  public:
    static Hex20RegionMesh from_unstructured_block(const UnstructuredHex20Mesh& source, std::int64_t block_id);

    const std::vector<CartesianPoint3>& nodes() const noexcept { return _nodes; }

    const std::vector<Hex20Element>& elements() const noexcept { return _elements; }

    const std::vector<bool>& temperature_nodes() const noexcept { return _temperature_nodes; }

    Hex20RegionBoundary map_side_set(const UnstructuredHex20Mesh& source, const std::string& side_set_name) const;

  private:
    Hex20RegionMesh(const UnstructuredHex20Mesh& source, std::int64_t block_id);
    std::vector<CartesianPoint3> _nodes;
    std::vector<Hex20Element> _elements;
    std::vector<bool> _temperature_nodes;
};
} // namespace fuelsim
