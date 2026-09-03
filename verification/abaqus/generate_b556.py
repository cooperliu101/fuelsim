#!/usr/bin/env python3
"""Generate the B5.56 full-size C3D20T finite-sliding friction case."""

import json
from pathlib import Path

from generate_b546 import deck
from generate_b551 import displacement_boundary_sets


JOB = "b556_m58_c3d20t_finite_sliding_friction"
MESH_JOB = "b548_m58_c3d20t_integrated"


def replace_once(contents, source, target, description):
    if contents.count(source) != 1:
        raise RuntimeError("unexpected B5.56 %s" % description)
    return contents.replace(source, target)


def main():
    directory = Path(__file__).resolve().parent
    for suffix in ("_mesh.inc", "_mesh.json"):
        if not (directory / (MESH_JOB + suffix)).is_file():
            raise RuntimeError("generate B5.48 mesh artifacts before B5.56")

    with (directory / (MESH_JOB + "_mesh.json")).open() as source:
        mesh = json.load(source)

    contents = deck(MESH_JOB, "B5.56", "C3D20T")
    include = "*Include, input=%s_mesh.inc\n" % MESH_JOB
    contents = replace_once(
        contents,
        include,
        include + displacement_boundary_sets(mesh).replace("B551_", "B556_"),
        "mesh include",
    )
    for source_name, target_name in (
        ("FUEL_BOTTOM, 1, 2, 0", "B556_FUEL_BOTTOM_U, 1, 2, 0"),
        ("CLAD_BOTTOM, 1, 3, 0", "B556_CLAD_BOTTOM_U, 1, 3, 0"),
        ("FUEL_TOP, 3, 3, 1", "B556_FUEL_TOP_U, 3, 3, 1"),
        ("CLAD_TOP, 3, 3, 1", "B556_CLAD_TOP_U, 3, 3, 1"),
    ):
        contents = replace_once(contents, source_name, target_name, "mechanical boundary")
    contents = replace_once(contents, "FUEL, BF, 2e7", "FUEL, BF, 2e8", "fuel heat source")
    contents = replace_once(
        contents, "CEEQ, IVOL, PEEQ, S", "CEEQ, COORD, IVOL, PEEQ, S", "element output request"
    )
    contents = replace_once(
        contents,
        "*Friction\n0.002\n",
        "*Friction\n0.002\n*Gap Heat Generation\n0, 0.5\n",
        "friction definition",
    )
    contents = replace_once(
        contents, "type=NODE TO SURFACE", "type=SURFACE TO SURFACE", "contact-pair definition"
    )
    contents = replace_once(
        contents,
        "CDISP, CFORCE, CSTRESS, HFL",
        "CDISP, CFORCE, CSTRESS, CTANDIR, HFL",
        "contact output request",
    )

    path = directory / (JOB + ".inp")
    path.write_text(contents, encoding="ascii")
    print("wrote %s using the tracked B5.48 mesh artifacts" % path)


if __name__ == "__main__":
    main()
