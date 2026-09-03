#!/usr/bin/env python3
"""Generate the B5.51 M5.8 C3D20T finite-sliding contact case."""

import json
from pathlib import Path

from generate_b546 import deck


JOB = "b551_m58_c3d20t_finite_sliding"
MESH_JOB = "b548_m58_c3d20t_integrated"

HEX20_FACE_NODES = (
    (0, 1, 5, 4, 8, 13, 16, 12),
    (1, 2, 6, 5, 9, 14, 17, 13),
    (2, 3, 7, 6, 10, 15, 18, 14),
    (3, 0, 4, 7, 11, 12, 19, 15),
    (0, 3, 2, 1, 11, 10, 9, 8),
    (4, 5, 6, 7, 16, 17, 18, 19),
)


def full_displacement_nodes(mesh, side_set_name):
    elements = {element["label"]: element["nodes"] for element in mesh["elements"]}
    matching = [side_set for side_set in mesh["side_sets"] if side_set["name"] == side_set_name]
    if len(matching) != 1:
        raise RuntimeError("expected one %s side set in the B5.51 mesh manifest" % side_set_name)
    nodes = set()
    for face in matching[0]["faces"]:
        connectivity = elements[face["element"]]
        for local in HEX20_FACE_NODES[face["exodus_side"] - 1]:
            nodes.add(connectivity[local])
    return sorted(nodes)


def displacement_boundary_sets(mesh):
    lines = []
    expected_counts = {"FUEL_BOTTOM": 257, "FUEL_TOP": 257, "CLAD_BOTTOM": 128, "CLAD_TOP": 128}
    for name, expected_count in expected_counts.items():
        labels = full_displacement_nodes(mesh, name)
        if len(labels) != expected_count:
            raise RuntimeError(
                "B5.51 %s displacement boundary has %d nodes, expected %d"
                % (name, len(labels), expected_count)
            )
        lines.append("*Nset, nset=B551_%s_U" % name)
        for begin in range(0, len(labels), 16):
            lines.append(", ".join(str(label) for label in labels[begin : begin + 16]))
    return "\n".join(lines) + "\n"


def main():
    directory = Path(__file__).resolve().parent
    for suffix in ("_mesh.inc", "_mesh.json"):
        if not (directory / (MESH_JOB + suffix)).is_file():
            raise RuntimeError("generate B5.48 mesh artifacts before B5.51")

    with (directory / (MESH_JOB + "_mesh.json")).open() as source:
        mesh = json.load(source)

    contents = deck(MESH_JOB, "B5.51", "C3D20T")
    include = "*Include, input=%s_mesh.inc\n" % MESH_JOB
    if contents.count(include) != 1:
        raise RuntimeError("unexpected B5.51 mesh include")
    contents = contents.replace(include, include + displacement_boundary_sets(mesh))
    for source_name, target_name in (
        ("FUEL_BOTTOM, 1, 2, 0", "B551_FUEL_BOTTOM_U, 1, 2, 0"),
        ("CLAD_BOTTOM, 1, 3, 0", "B551_CLAD_BOTTOM_U, 1, 3, 0"),
        ("FUEL_TOP, 3, 3, 1", "B551_FUEL_TOP_U, 3, 3, 1"),
        ("CLAD_TOP, 3, 3, 1", "B551_CLAD_TOP_U, 3, 3, 1"),
    ):
        if contents.count(source_name) != 1:
            raise RuntimeError("unexpected B5.51 mechanical boundary: %s" % source_name)
        contents = contents.replace(source_name, target_name)
    nominal_heat_source = "FUEL, BF, 2e7"
    if contents.count(nominal_heat_source) != 1:
        raise RuntimeError("unexpected B5.51 fuel heat source")
    contents = contents.replace(nominal_heat_source, "FUEL, BF, 2e8")
    requested = "CEEQ, IVOL, PEEQ, S"
    if contents.count(requested) != 1:
        raise RuntimeError("unexpected B5.51 element output request")
    contents = contents.replace(requested, "CEEQ, COORD, IVOL, PEEQ, S")
    friction = "*Friction\n0.002\n"
    if contents.count(friction) != 1:
        raise RuntimeError("unexpected B5.51 friction definition")
    contents = contents.replace(friction, "")
    node_to_surface = "type=NODE TO SURFACE"
    if contents.count(node_to_surface) != 1:
        raise RuntimeError("unexpected B5.51 contact-pair definition")
    contents = contents.replace(node_to_surface, "type=SURFACE TO SURFACE")

    path = directory / (JOB + ".inp")
    path.write_text(contents, encoding="ascii")
    print("wrote %s using the tracked B5.48 mesh artifacts" % path)


if __name__ == "__main__":
    main()
