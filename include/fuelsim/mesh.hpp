#ifndef FUELSIM_MESH_HPP
#define FUELSIM_MESH_HPP

#include <array>
#include <cstddef>
#include <cstdint>
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
                          std::vector<std::int64_t> element_block_ids);

    const std::vector<RzPoint>& nodes() const noexcept;
    const std::vector<Quad4Element>& elements() const noexcept;
    const std::vector<std::int64_t>& element_block_ids() const noexcept;

  private:
    std::vector<RzPoint> _nodes;
    std::vector<Quad4Element> _elements;
    std::vector<std::int64_t> _element_block_ids;
};

class StructuredRzMesh final {
  public:
    static StructuredRzMesh make_annulus(double inner_radius,
                                         double outer_radius, double length,
                                         std::size_t radial_elements,
                                         std::size_t axial_elements);

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
