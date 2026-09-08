"""Generate tracked CAX8T material, contact, and sliding validation meshes."""
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


ROOT = Path(__file__).resolve().parent
MESH_ROOT = ROOT.parent / "meshes"


def append_nset(lines, name, nodes):
    lines.append("*Nset, nset=" + name.upper())
    for first in range(0, len(nodes), 16):
        lines.append(", ".join(str(node) for node in nodes[first:first + 16]))


def add_names(file, variable_name, dimension, values):
    output = np.zeros((len(values), 33), dtype="S1")
    for index, value in enumerate(values):
        output[index, :len(value)] = np.frombuffer(value.encode("ascii"), dtype="S1")
    variable = file.createVariable(variable_name, "c", (dimension, "len_name"))
    variable[:] = output


def add_variable(file, name, dimensions, values, kind="i"):
    variable = file.createVariable(name, kind, dimensions)
    variable[:] = values
    return variable


def write_exodus(path, title, coordinates, blocks, boundaries, node_sets):
    with netcdf_file(str(path), "w") as file:
        file.api_version = np.float32(7.22)
        file.version = np.float32(7.22)
        file.floating_point_word_size = np.int32(8)
        file.file_size = np.int32(1)
        file.title = title
        dimensions = [
            ("time_step", None), ("len_name", 33), ("num_dim", 2),
            ("num_nodes", len(coordinates)),
            ("num_elem", sum(len(elements) for _, elements in blocks)),
            ("num_el_blk", len(blocks)), ("num_side_sets", len(boundaries)),
            ("num_node_sets", len(node_sets)),
        ]
        for name, size in dimensions:
            file.createDimension(name, size)

        add_variable(file, "coordx", ("num_nodes",), [point[0] for point in coordinates], "d")
        add_variable(file, "coordy", ("num_nodes",), [point[1] for point in coordinates], "d")
        for prefix, dimension, count in [
                ("eb", "num_el_blk", len(blocks)),
                ("ss", "num_side_sets", len(boundaries)),
                ("ns", "num_node_sets", len(node_sets))]:
            ids = add_variable(file, prefix + "_prop1", (dimension,), list(range(1, count + 1)))
            ids.name = "ID"
            add_variable(file, prefix + "_status", (dimension,), [1] * count)

        add_names(file, "eb_names", "num_el_blk", [name for name, _ in blocks])
        add_names(file, "ss_names", "num_side_sets", [name for name, _, _ in boundaries])
        add_names(file, "ns_names", "num_node_sets", [name for name, _ in node_sets])
        for index, (_, elements) in enumerate(blocks, 1):
            element_dimension = "num_el_in_blk%d" % index
            node_dimension = "num_nod_per_el%d" % index
            file.createDimension(element_dimension, len(elements))
            file.createDimension(node_dimension, len(elements[0]))
            connectivity = add_variable(
                file, "connect%d" % index, (element_dimension, node_dimension), elements)
            connectivity.elem_type = "QUAD8" if len(elements[0]) == 8 else "QUAD4"
        for index, (_, elements, sides) in enumerate(boundaries, 1):
            dimension = "num_side_ss%d" % index
            file.createDimension(dimension, len(sides))
            add_variable(file, "elem_ss%d" % index, (dimension,), elements)
            add_variable(file, "side_ss%d" % index, (dimension,), sides)
        for index, (_, nodes) in enumerate(node_sets, 1):
            dimension = "num_nod_ns%d" % index
            file.createDimension(dimension, len(nodes))
            add_variable(file, "node_ns%d" % index, (dimension,), nodes)


def write_material_mesh():
    coordinates = [
        (0.0, 0.0), (0.001, 0.0), (0.001, 0.001), (0.0, 0.001),
        (0.0005, 0.0), (0.001, 0.0005), (0.0005, 0.001), (0.0, 0.0005),
    ]
    blocks = [("solid", [[1, 2, 3, 4, 5, 6, 7, 8]])]
    boundaries = [("left", [1], [4]), ("bottom", [1], [1]), ("top", [1], [3])]
    node_sets = [("all", list(range(1, 9))), ("left", [1, 4, 8]),
                 ("bottom", [1, 2, 5]), ("top", [3, 4, 7])]
    write_exodus(MESH_ROOT / "b10_cax8t_material.e", "CAX8T material validation", coordinates,
                  blocks, boundaries, node_sets)
    lines = ["** CAX8T material mesh matching b10_cax8t_material.e", "*Node"]
    for label, point in enumerate(coordinates, 1):
        lines.append("%d, %.17g, %.17g" % (label, point[0], point[1]))
    lines += ["*Element, type=CAX8T, elset=SOLID", "1, 1, 2, 3, 4, 5, 6, 7, 8"]
    for name, nodes in node_sets:
        append_nset(lines, name, nodes)
    append_nset(lines, "axis", [1, 4, 8])
    (ROOT / "b10_cax8t_material_mesh.inc").write_text("\n".join(lines) + "\n")


def write_contact_mesh():
    lower = [
        (0.001, 0.0), (0.002, 0.0), (0.002, 0.001), (0.001, 0.001),
        (0.0015, 0.0), (0.002, 0.0005), (0.0015, 0.001), (0.001, 0.0005),
    ]
    upper = [(r, z + 0.001001) for r, z in lower]
    coordinates = lower + upper
    blocks = [("lower", [[1, 2, 3, 4, 5, 6, 7, 8]]),
              ("upper", [[9, 10, 11, 12, 13, 14, 15, 16]])]
    boundaries = [("bottom", [1], [1]), ("lower_contact", [1], [3]),
                  ("upper_contact", [2], [1]), ("top", [2], [3]),
                  ("lower_radial", [1, 1], [2, 4]),
                  ("upper_radial", [2, 2], [2, 4])]
    node_sets = [
        ("all", list(range(1, 17))), ("bottom", [1, 2, 5]),
        ("lower_contact", [3, 4, 7]), ("upper_contact", [9, 10, 13]),
        ("top", [11, 12, 15]), ("lower_radial", list(range(1, 9))),
        ("upper_radial", list(range(9, 17))), ("radial_fixed", list(range(1, 17))),
    ]
    write_exodus(MESH_ROOT / "b10_cax8t_contact.e", "CAX8T contact validation", coordinates,
                  blocks, boundaries, node_sets)
    lines = ["** CAX8T contact mesh matching b10_cax8t_contact.e", "*Node"]
    for label, point in enumerate(coordinates, 1):
        lines.append("%d, %.17g, %.17g" % (label, point[0], point[1]))
    lines += ["*Element, type=CAX8T, elset=LOWER", "1, 1, 2, 3, 4, 5, 6, 7, 8",
              "*Element, type=CAX8T, elset=UPPER", "2, 9, 10, 11, 12, 13, 14, 15, 16"]
    for name, nodes in node_sets:
        append_nset(lines, name, nodes)
    lines += ["*Surface, type=ELEMENT, name=S_LOWER", "LOWER, S3",
              "*Surface, type=ELEMENT, name=S_UPPER", "UPPER, S1"]
    (ROOT / "b10_cax8t_contact_mesh.inc").write_text("\n".join(lines) + "\n")


def enrich_quad4_mesh(corner_coordinates, block_corner_elements):
    coordinates = list(corner_coordinates)
    edge_nodes = {}
    blocks = []
    for name, elements in block_corner_elements:
        enriched = []
        for element in elements:
            midsides = []
            for first, second in zip(element, element[1:] + element[:1]):
                edge = tuple(sorted((first, second)))
                if edge not in edge_nodes:
                    point_a, point_b = coordinates[first - 1], coordinates[second - 1]
                    coordinates.append(((point_a[0] + point_b[0]) / 2,
                                        (point_a[1] + point_b[1]) / 2))
                    edge_nodes[edge] = len(coordinates)
                midsides.append(edge_nodes[edge])
            enriched.append(element + midsides)
        blocks.append((name, enriched))
    return coordinates, blocks


def write_sliding_mesh():
    coordinates = []
    corner_blocks = []
    boundaries = []
    element_label = 1
    for name, radii, heights in [
            ("inner", (0.001, 0.002), (0.00125, 0.00225)),
            ("outer", (0.002001, 0.003001), (0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006))]:
        start = len(coordinates) + 1
        for axial in heights:
            coordinates.extend((radius, axial) for radius in radii)
        elements = [[start + 2 * row, start + 2 * row + 1,
                     start + 2 * row + 3, start + 2 * row + 2]
                    for row in range(len(heights) - 1)]
        labels = list(range(element_label, element_label + len(elements)))
        boundaries += [(name + "_left", labels, [4] * len(labels)),
                       (name + "_right", labels, [2] * len(labels)),
                       (name + "_bottom", [labels[0]], [1]),
                       (name + "_top", [labels[-1]], [3])]
        corner_blocks.append((name, elements))
        element_label += len(elements)
    coordinates, blocks = enrich_quad4_mesh(coordinates, corner_blocks)

    all_elements = {}
    label = 1
    for _, elements in blocks:
        for element in elements:
            all_elements[label] = element
            label += 1
    side_local_nodes = {1: (0, 1, 4), 2: (1, 2, 5), 3: (2, 3, 6), 4: (3, 0, 7)}
    node_sets = []
    for name, elements, sides in boundaries:
        nodes = []
        for element, side in zip(elements, sides):
            for local in side_local_nodes[side]:
                node = all_elements[element][local]
                if node not in nodes:
                    nodes.append(node)
        node_sets.append((name, nodes))
    node_sets.append(("all", list(range(1, len(coordinates) + 1))))
    write_exodus(MESH_ROOT / "b10_cax8t_sliding.e", "CAX8T finite sliding validation",
                  coordinates, blocks, boundaries, node_sets)

    lines = ["** CAX8T sliding mesh matching b10_cax8t_sliding.e", "*Node"]
    for node, point in enumerate(coordinates, 1):
        lines.append("%d, %.17g, %.17g" % (node, point[0], point[1]))
    label = 1
    for name, elements in blocks:
        lines.append("*Element, type=CAX8T, elset=" + name.upper())
        for element in elements:
            lines.append(", ".join(str(value) for value in [label] + element))
            label += 1
    for name, nodes in node_sets:
        append_nset(lines, name, nodes)
    lines += ["*Elset, elset=ALL_ELEMENTS", "INNER, OUTER",
              "*Surface, type=ELEMENT, name=S_INNER", "INNER, S2",
              "*Surface, type=ELEMENT, name=S_OUTER", "OUTER, S4"]
    (ROOT / "b10_cax8t_sliding_mesh.inc").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    write_material_mesh()
    write_contact_mesh()
    write_sliding_mesh()
