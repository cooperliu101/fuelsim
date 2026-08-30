#!/usr/bin/env python3
"""Convert a three-dimensional Exodus HEX8 mesh to an Abaqus mesh include.

The converter preserves Exodus node and element numbers, element-block names,
node-set names, and side-set names.  Exodus HEX8 side numbers are translated
to Abaqus C3D8 face labels.  The generated JSON manifest keeps the exact
reference coordinates and connectivity for deterministic result extraction.
"""

import argparse
import json
import math
import re
from pathlib import Path

try:
    from scipy.io import netcdf_file
except ImportError as error:
    raise RuntimeError("exodus_to_abaqus.py requires scipy.io.netcdf_file") from error


EXODUS_TO_ABAQUS_HEX8_FACE = ("S3", "S4", "S5", "S6", "S1", "S2")


def _dimension_size(database, name):
    if name not in database.dimensions:
        return 0
    return int(database.dimensions[name])


def _decode_names(values):
    result = []
    for row in values:
        raw = b"".join(bytes(value) for value in row).split(b"\0", 1)[0]
        result.append(raw.decode("utf-8").strip())
    return result


def _abaqus_name(name, fallback):
    candidate = re.sub(r"[^A-Za-z0-9_]", "_", name.strip()).strip("_").upper()
    if not candidate:
        candidate = fallback
    if candidate[0].isdigit():
        candidate = "N_" + candidate
    return candidate[:80]


def _unique_names(names, category):
    result = []
    used = set()
    for index, name in enumerate(names, 1):
        candidate = _abaqus_name(name, "%s_%d" % (category, index))
        if candidate in used:
            raise RuntimeError("duplicate Abaqus %s name after sanitization: %s" % (category, candidate))
        used.add(candidate)
        result.append(candidate)
    return result


def _mapped_label(labels, internal_value, entity):
    internal = int(internal_value)
    if internal < 1 or internal > len(labels):
        raise RuntimeError("Exodus %s internal number is out of range: %d" % (entity, internal))
    return labels[internal - 1]


def read_exodus_hex8(path):
    """Read the mesh entities needed by an Abaqus HEX8 input deck."""
    database = netcdf_file(str(path), "r", mmap=False)
    try:
        if _dimension_size(database, "num_dim") != 3:
            raise RuntimeError("Exodus-to-Abaqus conversion requires a three-dimensional mesh")
        node_count = _dimension_size(database, "num_nodes")
        element_count = _dimension_size(database, "num_elem")
        block_count = _dimension_size(database, "num_el_blk")
        if node_count == 0 or element_count == 0 or block_count == 0:
            raise RuntimeError("Exodus mesh must contain nodes, HEX8 elements, and element blocks")

        variables = database.variables
        node_labels = [int(value) for value in variables.get("node_num_map", range(1, node_count + 1))[:]]
        element_labels = [int(value) for value in variables.get("elem_num_map", range(1, element_count + 1))[:]]
        if len(set(node_labels)) != node_count or len(set(element_labels)) != element_count:
            raise RuntimeError("Exodus node and element number maps must be unique")
        if min(node_labels) < 1 or min(element_labels) < 1:
            raise RuntimeError("Abaqus conversion requires positive node and element numbers")
        coordinates = [
            [float(variables[axis][index]) for axis in ("coordx", "coordy", "coordz")]
            for index in range(node_count)
        ]
        if any(not math.isfinite(value) for point in coordinates for value in point):
            raise RuntimeError("Exodus coordinates must be finite")

        block_ids = [int(value) for value in variables["eb_prop1"][:]]
        block_source_names = _decode_names(variables["eb_names"][:])
        block_names = _unique_names(block_source_names, "BLOCK")
        blocks = []
        elements = []
        element_offset = 0
        for block_index in range(block_count):
            variable = variables["connect%d" % (block_index + 1)]
            element_type = getattr(variable, "elem_type", b"").decode("ascii").upper().strip()
            connectivity = variable[:]
            if connectivity.ndim != 2 or connectivity.shape[1] != 8 or not element_type.startswith("HEX"):
                raise RuntimeError("all Exodus element blocks must contain eight-node HEX elements")
            labels = element_labels[element_offset : element_offset + connectivity.shape[0]]
            block_elements = []
            for local_index, internal_nodes in enumerate(connectivity):
                nodes = [_mapped_label(node_labels, internal_node, "connectivity node") for internal_node in internal_nodes]
                record = {"label": labels[local_index], "nodes": nodes, "block": block_names[block_index]}
                elements.append(record)
                block_elements.append(record)
            blocks.append(
                {
                    "id": block_ids[block_index],
                    "source_name": block_source_names[block_index],
                    "name": block_names[block_index],
                    "elements": block_elements,
                }
            )
            element_offset += connectivity.shape[0]
        if element_offset != element_count:
            raise RuntimeError("Exodus block element counts do not match num_elem")

        node_set_count = _dimension_size(database, "num_node_sets")
        node_sets = []
        if node_set_count:
            ids = [int(value) for value in variables["ns_prop1"][:]]
            source_names = _decode_names(variables["ns_names"][:])
            names = _unique_names(source_names, "NODE_SET")
            if "ALL_NODES" in names:
                raise RuntimeError("Exodus node-set name ALL_NODES conflicts with the generated Abaqus set")
            for index in range(node_set_count):
                internal = variables["node_ns%d" % (index + 1)][:]
                node_sets.append(
                    {
                        "id": ids[index],
                        "source_name": source_names[index],
                        "name": names[index],
                        "nodes": [_mapped_label(node_labels, value, "node-set node") for value in internal],
                    }
                )

        side_set_count = _dimension_size(database, "num_side_sets")
        side_sets = []
        if side_set_count:
            ids = [int(value) for value in variables["ss_prop1"][:]]
            source_names = _decode_names(variables["ss_names"][:])
            names = _unique_names(source_names, "SIDE_SET")
            for index in range(side_set_count):
                internal_elements = variables["elem_ss%d" % (index + 1)][:]
                exodus_sides = variables["side_ss%d" % (index + 1)][:]
                faces = []
                for internal_element, exodus_side in zip(internal_elements, exodus_sides):
                    side = int(exodus_side)
                    if side < 1 or side > 6:
                        raise RuntimeError("Exodus HEX8 side number must be in [1, 6]")
                    faces.append(
                        {
                            "element": _mapped_label(element_labels, internal_element, "side-set element"),
                            "exodus_side": side,
                            "abaqus_face": EXODUS_TO_ABAQUS_HEX8_FACE[side - 1],
                        }
                    )
                side_sets.append(
                    {
                        "id": ids[index],
                        "source_name": source_names[index],
                        "name": names[index],
                        "faces": faces,
                    }
                )

        return {
            "source": str(Path(path)),
            "nodes": [
                {"label": label, "coordinates": point} for label, point in zip(node_labels, coordinates)
            ],
            "elements": elements,
            "blocks": blocks,
            "node_sets": node_sets,
            "side_sets": side_sets,
        }
    finally:
        database.close()


def _append_labels(lines, labels):
    for begin in range(0, len(labels), 16):
        lines.append(", ".join(str(label) for label in labels[begin : begin + 16]))


def write_abaqus_hex8_mesh(mesh, output_path, element_type="C3D8T"):
    """Write an Abaqus mesh include from the dictionary returned above."""
    lines = [
        "** Generated by exodus_to_abaqus.py; do not edit by hand.",
        "** Exodus HEX8 sides 1..6 map to Abaqus S3,S4,S5,S6,S1,S2.",
        "*Node",
    ]
    for node in mesh["nodes"]:
        lines.append(
            "%d, %.17g, %.17g, %.17g" % ((node["label"],) + tuple(node["coordinates"]))
        )
    for block in mesh["blocks"]:
        lines.append("*Element, type=%s, elset=%s" % (element_type, block["name"]))
        for element in block["elements"]:
            lines.append("%d, %s" % (element["label"], ", ".join(str(node) for node in element["nodes"])))
    lines.append("*Nset, nset=ALL_NODES")
    _append_labels(lines, [node["label"] for node in mesh["nodes"]])
    for node_set in mesh["node_sets"]:
        lines.append("*Nset, nset=%s" % node_set["name"])
        _append_labels(lines, node_set["nodes"])
    for side_set in mesh["side_sets"]:
        lines.append("*Surface, type=ELEMENT, name=%s" % side_set["name"])
        for face in side_set["faces"]:
            lines.append("%d, %s" % (face["element"], face["abaqus_face"]))
    Path(output_path).write_text("\n".join(lines) + "\n", encoding="ascii")


def convert(exodus_path, abaqus_path, manifest_path=None, element_type="C3D8T"):
    mesh = read_exodus_hex8(exodus_path)
    write_abaqus_hex8_mesh(mesh, abaqus_path, element_type)
    if manifest_path:
        Path(manifest_path).write_text(
            json.dumps(mesh, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8"
        )
    return mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exodus", help="input Exodus HEX8 mesh")
    parser.add_argument("abaqus", help="output Abaqus mesh include")
    parser.add_argument("--manifest", help="optional JSON mesh manifest")
    parser.add_argument("--element-type", default="C3D8T", help="Abaqus HEX8 element type")
    arguments = parser.parse_args()
    mesh = convert(arguments.exodus, arguments.abaqus, arguments.manifest, arguments.element_type)
    print(
        "converted %d nodes, %d elements, %d blocks, %d node sets, and %d side sets"
        % (
            len(mesh["nodes"]),
            len(mesh["elements"]),
            len(mesh["blocks"]),
            len(mesh["node_sets"]),
            len(mesh["side_sets"]),
        )
    )


if __name__ == "__main__":
    main()
