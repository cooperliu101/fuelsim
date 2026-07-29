#ifndef FUELSIM_MESH_HPP
#define FUELSIM_MESH_HPP

#include <array>
#include <cstddef>
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
    double inner_radius_ = 0.0;
    double outer_radius_ = 0.0;
    double length_ = 0.0;
    std::size_t radial_elements_ = 0;
    std::size_t axial_elements_ = 0;

    std::vector<RzPoint> nodes_;
    std::vector<Quad4Element> elements_;

    std::vector<std::size_t> radial_inner_nodes_;
    std::vector<std::size_t> radial_outer_nodes_;
    std::vector<std::size_t> bottom_nodes_;
    std::vector<std::size_t> top_nodes_;

    std::vector<Line2BoundaryElement> radial_inner_elements_;
    std::vector<Line2BoundaryElement> radial_outer_elements_;
    std::vector<Line2BoundaryElement> bottom_elements_;
    std::vector<Line2BoundaryElement> top_elements_;
};

} // namespace fuelsim

#endif
