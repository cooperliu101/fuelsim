#include "fuelsim/mesh.hpp"

#include <limits>
#include <stdexcept>

namespace fuelsim {

StructuredRzMesh StructuredRzMesh::make_annulus(double inner_radius,
                                                double outer_radius,
                                                double length,
                                                std::size_t radial_elements,
                                                std::size_t axial_elements) {
    if (!(inner_radius >= 0.0))
        throw std::invalid_argument(
            "StructuredRzMesh inner_radius must be nonnegative");
    if (!(outer_radius > inner_radius))
        throw std::invalid_argument(
            "StructuredRzMesh outer_radius must exceed inner_radius");
    if (!(length > 0.0))
        throw std::invalid_argument("StructuredRzMesh length must be positive");
    if (radial_elements == 0 || axial_elements == 0)
        throw std::invalid_argument(
            "StructuredRzMesh element counts must be positive");

    const std::size_t radial_nodes = radial_elements + 1;
    const std::size_t axial_nodes = axial_elements + 1;
    if (radial_nodes > std::numeric_limits<std::size_t>::max() / axial_nodes)
        throw std::length_error("StructuredRzMesh node count overflows");

    StructuredRzMesh mesh;
    mesh.inner_radius_ = inner_radius;
    mesh.outer_radius_ = outer_radius;
    mesh.length_ = length;
    mesh.radial_elements_ = radial_elements;
    mesh.axial_elements_ = axial_elements;

    mesh.nodes_.reserve(radial_nodes * axial_nodes);
    for (std::size_t iz = 0; iz < axial_nodes; ++iz) {
        const double z = length * static_cast<double>(iz) /
                         static_cast<double>(axial_elements);
        for (std::size_t ir = 0; ir < radial_nodes; ++ir) {
            const double fraction =
                static_cast<double>(ir) / static_cast<double>(radial_elements);
            const double r =
                inner_radius + (outer_radius - inner_radius) * fraction;
            mesh.nodes_.push_back({r, z});
        }
    }

    mesh.elements_.reserve(radial_elements * axial_elements);
    for (std::size_t iz = 0; iz < axial_elements; ++iz) {
        for (std::size_t ir = 0; ir < radial_elements; ++ir) {
            mesh.elements_.push_back(
                {{mesh.node_id(ir, iz), mesh.node_id(ir + 1, iz),
                  mesh.node_id(ir + 1, iz + 1), mesh.node_id(ir, iz + 1)}});
        }
    }

    mesh.radial_inner_nodes_.reserve(axial_nodes);
    mesh.radial_outer_nodes_.reserve(axial_nodes);
    for (std::size_t iz = 0; iz < axial_nodes; ++iz) {
        mesh.radial_inner_nodes_.push_back(mesh.node_id(0, iz));
        mesh.radial_outer_nodes_.push_back(mesh.node_id(radial_elements, iz));
    }

    mesh.bottom_nodes_.reserve(radial_nodes);
    mesh.top_nodes_.reserve(radial_nodes);
    for (std::size_t ir = 0; ir < radial_nodes; ++ir) {
        mesh.bottom_nodes_.push_back(mesh.node_id(ir, 0));
        mesh.top_nodes_.push_back(mesh.node_id(ir, axial_elements));
    }

    mesh.radial_inner_elements_.reserve(axial_elements);
    mesh.radial_outer_elements_.reserve(axial_elements);
    for (std::size_t iz = 0; iz < axial_elements; ++iz) {
        mesh.radial_inner_elements_.push_back(
            {{{mesh.node_id(0, iz), mesh.node_id(0, iz + 1)}}});
        mesh.radial_outer_elements_.push_back(
            {{{mesh.node_id(radial_elements, iz),
               mesh.node_id(radial_elements, iz + 1)}}});
    }

    mesh.bottom_elements_.reserve(radial_elements);
    mesh.top_elements_.reserve(radial_elements);
    for (std::size_t ir = 0; ir < radial_elements; ++ir) {
        mesh.bottom_elements_.push_back(
            {{{mesh.node_id(ir, 0), mesh.node_id(ir + 1, 0)}}});
        mesh.top_elements_.push_back(
            {{{mesh.node_id(ir, axial_elements),
               mesh.node_id(ir + 1, axial_elements)}}});
    }

    return mesh;
}

const std::vector<RzPoint>& StructuredRzMesh::nodes() const noexcept {
    return nodes_;
}

const std::vector<Quad4Element>& StructuredRzMesh::elements() const noexcept {
    return elements_;
}

const std::vector<std::size_t>&
StructuredRzMesh::boundary_nodes(BoundaryId boundary) const {
    switch (boundary) {
    case BoundaryId::radial_inner:
        return radial_inner_nodes_;
    case BoundaryId::radial_outer:
        return radial_outer_nodes_;
    case BoundaryId::bottom:
        return bottom_nodes_;
    case BoundaryId::top:
        return top_nodes_;
    }
    throw std::invalid_argument("Unknown StructuredRzMesh boundary");
}

const std::vector<Line2BoundaryElement>&
StructuredRzMesh::boundary_elements(BoundaryId boundary) const {
    switch (boundary) {
    case BoundaryId::radial_inner:
        return radial_inner_elements_;
    case BoundaryId::radial_outer:
        return radial_outer_elements_;
    case BoundaryId::bottom:
        return bottom_elements_;
    case BoundaryId::top:
        return top_elements_;
    }
    throw std::invalid_argument("Unknown StructuredRzMesh boundary");
}

std::size_t StructuredRzMesh::radial_elements() const noexcept {
    return radial_elements_;
}

std::size_t StructuredRzMesh::axial_elements() const noexcept {
    return axial_elements_;
}

std::size_t StructuredRzMesh::node_id(std::size_t radial_index,
                                      std::size_t axial_index) const {
    if (radial_index > radial_elements_ || axial_index > axial_elements_)
        throw std::out_of_range("StructuredRzMesh node index is out of range");
    return axial_index * (radial_elements_ + 1) + radial_index;
}

double StructuredRzMesh::inner_radius() const noexcept {
    return inner_radius_;
}

double StructuredRzMesh::outer_radius() const noexcept {
    return outer_radius_;
}

double StructuredRzMesh::length() const noexcept {
    return length_;
}

} // namespace fuelsim
