#pragma once
#include "core/mesh.hpp"
#include "cpeg8t.hpp"

namespace fuelsim {
struct PlaneQuad8Element final {
    std::array<std::size_t, 8> nodes;
};

class UnstructuredPlaneQuad8Mesh final : public UnstructuredMeshMetadata {
  public:
    UnstructuredPlaneQuad8Mesh(std::vector<std::array<double, 2>> nodes,
        std::vector<PlaneQuad8Element> elements,
        std::vector<std::int64_t> block_ids,
        std::vector<ElementBlockInfo> blocks,
        std::vector<NodeSet> node_sets,
        std::vector<SideSet> side_sets);

    const std::vector<std::array<double, 2>>& nodes() const noexcept { return _nodes; }

    const std::vector<PlaneQuad8Element>& elements() const noexcept { return _elements; }

  private:
    std::vector<std::array<double, 2>> _nodes;
    std::vector<PlaneQuad8Element> _elements;
};
} // namespace fuelsim
