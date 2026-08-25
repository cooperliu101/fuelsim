#!/usr/bin/env python3
"""Generate the controlled B4.4 finite-sliding averaging-operator probe."""

from pathlib import Path


directory = Path(__file__).resolve().parent
source = directory / "b44_hex8_nonmatching_partial_contact_swapped.inp"
target = directory / "b44_hex8_swapped_operator_probe.inp"
prefix = []
for line in source.read_text().splitlines():
    if line.lower().startswith("*step"):
        break
    prefix.append(line)

tiles = [list(range(first, first + 8)) for first in (51, 59, 67, 75, 83, 91)]
contact_nodes = [7, 11, 15, 23, 25, 27, 33, 35, 37]
noncontact_nodes = [node for node in range(1, 51) if node not in contact_nodes]


def labels(name, nodes):
    return ["*Nset, nset=" + name, ", ".join(str(node) for node in nodes)]


lines = prefix
lines[1] = "** B4.4 C3D8 finite-sliding stepped-primary averaging-operator probe."
for index, nodes in enumerate(tiles, 1):
    lines.extend(labels("TILE_%d" % index, nodes))
lines.extend(labels("PATCH_NONCONTACT", noncontact_nodes))

base = -16.0 / 1024.0
perturbation = 1.0e-6
steps = [("BASE", "base", 0, 0.0)]
for index in range(len(tiles)):
    steps.append(("T%d_PLUS" % (index + 1), "tile", index, perturbation))
    steps.append(("T%d_MINUS" % (index + 1), "tile", index, -perturbation))
for index in range(len(contact_nodes)):
    steps.append(("S%d_PLUS" % (index + 1), "secondary", index, perturbation))
    steps.append(("S%d_MINUS" % (index + 1), "secondary", index, -perturbation))

for name, kind, selected, delta in steps:
    lines.extend(
        [
            "*Step, name=%s, nlgeom=YES, inc=1" % name,
            "*Static",
            "1., 1., 1., 1.",
            "*Boundary, op=NEW",
            "PRIMARY_ALL, 2, 3, 0.",
            "PATCH_NONCONTACT, 1, 1, 0.",
            "SECONDARY_ALL, 2, 3, 0.",
        ]
    )
    for index, node in enumerate(contact_nodes):
        value = delta if kind == "secondary" and index == selected else 0.0
        lines.append("%d, 1, 1, %.16g" % (node, value))
    for index in range(len(tiles)):
        value = base + (delta if kind == "tile" and index == selected else 0.0)
        lines.append("TILE_%d, 1, 1, %.16g" % (index + 1, value))
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, RF, U",
            "*Contact Output",
            "CSTRESS, CDISP, CFORCE",
            "*End Step",
        ]
    )

target.write_text("\n".join(lines) + "\n")
print(target)
