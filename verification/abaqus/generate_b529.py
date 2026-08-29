#!/usr/bin/env python3
"""Generate regular and warped Abaqus C3D8RT heat-capacity probes."""

import argparse
from pathlib import Path


REGULAR = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)

WARPED = (
    (0.00, 0.00, 0.00),
    (1.20, 0.10, -0.05),
    (1.10, 1.00, 0.10),
    (-0.10, 0.90, 0.00),
    (0.05, -0.05, 1.00),
    (1.15, 0.00, 1.20),
    (1.00, 1.10, 1.10),
    (-0.05, 1.00, 0.90),
)


def node_label(case, local_node):
    return 8 * case + local_node + 1


def deck(coordinates, name):
    lines = [
        "*Heading",
        "** Abaqus/Standard C3D8RT transient heat-capacity matrix probe: %s." % name,
        "** SI units: metre, second, kelvin, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for case in range(8):
        for local_node, coordinate in enumerate(coordinates):
            shifted = (coordinate[0] + 3.0 * case, coordinate[1], coordinate[2])
            lines.append("%d, %.16e, %.16e, %.16e" % ((node_label(case, local_node),) + shifted))
    lines.append("*Element, type=C3D8RT, elset=VOLUME")
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--warped", action="store_true")
    args = parser.parse_args()
    name = "warped" if args.warped else "regular"
    path = Path(__file__).with_name("b529_hex8_c3d8rt_%s_capacity_probe.inp" % name)
    path.write_text(deck(WARPED if args.warped else REGULAR, name), encoding="utf-8")
    print("wrote %s" % path)
