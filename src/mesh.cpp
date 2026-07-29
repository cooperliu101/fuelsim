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
    mesh._inner_radius = inner_radius;
    mesh._outer_radius = outer_radius;
    mesh._length = length;
    mesh._radial_elements = radial_elements;
    mesh._axial_elements = axial_elements;

    mesh._nodes.reserve(radial_nodes * axial_nodes);
    for (std::size_t iz = 0; iz < axial_nodes; ++iz) {
        const double z = length * static_cast<double>(iz) /
                         static_cast<double>(axial_elements);
        for (std::size_t ir = 0; ir < radial_nodes; ++ir) {
            const double fraction =
                static_cast<double>(ir) / static_cast<double>(radial_elements);
            const double r =
                inner_radius + (outer_radius - inner_radius) * fraction;
            mesh._nodes.push_back({r, z});
        }
    }

    mesh._elements.reserve(radial_elements * axial_elements);
    for (std::size_t iz = 0; iz < axial_elements; ++iz) {
        for (std::size_t ir = 0; ir < radial_elements; ++ir) {
            mesh._elements.push_back(
                {{mesh.node_id(ir, iz), mesh.node_id(ir + 1, iz),
                  mesh.node_id(ir + 1, iz + 1), mesh.node_id(ir, iz + 1)}});
        }
    }

    mesh._radial_inner_nodes.reserve(axial_nodes);
    mesh._radial_outer_nodes.reserve(axial_nodes);
    for (std::size_t iz = 0; iz < axial_nodes; ++iz) {
        mesh._radial_inner_nodes.push_back(mesh.node_id(0, iz));
        mesh._radial_outer_nodes.push_back(mesh.node_id(radial_elements, iz));
    }

    mesh._bottom_nodes.reserve(radial_nodes);
    mesh._top_nodes.reserve(radial_nodes);
    for (std::size_t ir = 0; ir < radial_nodes; ++ir) {
        mesh._bottom_nodes.push_back(mesh.node_id(ir, 0));
        mesh._top_nodes.push_back(mesh.node_id(ir, axial_elements));
    }

    mesh._radial_inner_elements.reserve(axial_elements);
    mesh._radial_outer_elements.reserve(axial_elements);
    for (std::size_t iz = 0; iz < axial_elements; ++iz) {
        mesh._radial_inner_elements.push_back(
            {{{mesh.node_id(0, iz), mesh.node_id(0, iz + 1)}}});
        mesh._radial_outer_elements.push_back(
            {{{mesh.node_id(radial_elements, iz),
               mesh.node_id(radial_elements, iz + 1)}}});
    }

    mesh._bottom_elements.reserve(radial_elements);
    mesh._top_elements.reserve(radial_elements);
    for (std::size_t ir = 0; ir < radial_elements; ++ir) {
        mesh._bottom_elements.push_back(
            {{{mesh.node_id(ir, 0), mesh.node_id(ir + 1, 0)}}});
        mesh._top_elements.push_back(
            {{{mesh.node_id(ir, axial_elements),
               mesh.node_id(ir + 1, axial_elements)}}});
    }

    return mesh;
}

const std::vector<RzPoint>& StructuredRzMesh::nodes() const noexcept {
    return _nodes;
}

const std::vector<Quad4Element>& StructuredRzMesh::elements() const noexcept {
    return _elements;
}

const std::vector<std::size_t>&
StructuredRzMesh::boundary_nodes(BoundaryId boundary) const {
    switch (boundary) {
    case BoundaryId::radial_inner:
        return _radial_inner_nodes;
    case BoundaryId::radial_outer:
        return _radial_outer_nodes;
    case BoundaryId::bottom:
        return _bottom_nodes;
    case BoundaryId::top:
        return _top_nodes;
    }
    throw std::invalid_argument("Unknown StructuredRzMesh boundary");
}

const std::vector<Line2BoundaryElement>&
StructuredRzMesh::boundary_elements(BoundaryId boundary) const {
    switch (boundary) {
    case BoundaryId::radial_inner:
        return _radial_inner_elements;
    case BoundaryId::radial_outer:
        return _radial_outer_elements;
    case BoundaryId::bottom:
        return _bottom_elements;
    case BoundaryId::top:
        return _top_elements;
    }
    throw std::invalid_argument("Unknown StructuredRzMesh boundary");
}

std::size_t StructuredRzMesh::radial_elements() const noexcept {
    return _radial_elements;
}

std::size_t StructuredRzMesh::axial_elements() const noexcept {
    return _axial_elements;
}

std::size_t StructuredRzMesh::node_id(std::size_t radial_index,
                                      std::size_t axial_index) const {
    if (radial_index > _radial_elements || axial_index > _axial_elements)
        throw std::out_of_range("StructuredRzMesh node index is out of range");
    return axial_index * (_radial_elements + 1) + radial_index;
}

double StructuredRzMesh::inner_radius() const noexcept {
    return _inner_radius;
}

double StructuredRzMesh::outer_radius() const noexcept {
    return _outer_radius;
}

double StructuredRzMesh::length() const noexcept {
    return _length;
}

} // namespace fuelsim
