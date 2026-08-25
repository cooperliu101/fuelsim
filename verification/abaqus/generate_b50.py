#!/usr/bin/env python3
"""Generate the B5.0 C3D8T transient heat-capacity probe."""

from pathlib import Path


UNIT_COORDINATES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)


def node_label(case, local_node):
    return 8 * case + local_node + 1


def deck():
    lines = [
        "*Heading",
        "** B5.0 Abaqus/Standard C3D8T transient heat-capacity matrix probe.",
        "** SI units: metre, second, kelvin, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for case in range(8):
        for local_node, coordinate in enumerate(UNIT_COORDINATES):
            shifted = (coordinate[0] + 2.0 * case, coordinate[1], coordinate[2])
            lines.append("%d, %.16e, %.16e, %.16e" % ((node_label(case, local_node),) + shifted))
    lines.append("*Element, type=C3D8T, elset=VOLUME")
    for case in range(8):
        connectivity = [node_label(case, local_node) for local_node in range(8)]
        lines.append("%d, %s" % (case + 1, ", ".join(str(label) for label in connectivity)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, 64, 1",
            "*Material, name=THERMAL",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "0.0",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3",
            "*Specific Heat",
            "3.0000000000000000e3",
            "*Solid Section, elset=VOLUME, material=THERMAL",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
            "*Step, name=CAPACITY, nlgeom=NO, inc=4",
            "*Coupled Temperature-Displacement",
            "1.0, 1.0, 1.0, 1.0",
            "*Boundary",
            "ALL_NODES, 1, 3, 0.0",
        ]
    )
    for case in range(8):
        for local_node in range(8):
            temperature = 301.0 if local_node == case else 300.0
            lines.append("%d, 11, 11, %.16e" % (node_label(case, local_node), temperature))
    lines.extend(["*Output, field, frequency=1", "*Node Output", "NT, RFL", "*End Step"])
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b50_hex8_c3d8t_capacity_probe.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
