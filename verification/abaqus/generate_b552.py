#!/usr/bin/env python3
"""Generate the lightweight B5.52 nonuniform curved HEX8 STS operator case."""

import math
from pathlib import Path


FACETS = 3
PRIMARY_RADIUS = 1.0
SECONDARY_RADIUS = 1.005


def cylindrical(radius, angle, z):
    return (radius * math.cos(angle), radius * math.sin(angle), z)


def append_annular(nodes, node_map, r0, r1, a0, a1):
    logical = (
        (r0, a0, 0.0), (r1, a0, 0.0), (r1, a1, 0.0), (r0, a1, 0.0),
        (r0, a0, 1.0), (r1, a0, 1.0), (r1, a1, 1.0), (r0, a1, 1.0),
    )
    element = []
    for point in logical:
        if point not in node_map:
            node_map[point] = len(nodes) + 1
            nodes.append(cylindrical(*point))
        element.append(node_map[point])
    return element


def inward_displacement(angle, z):
    return 0.014 + 0.003 * math.cos(2.0 * angle) + 0.002 * (z - 0.5)


def deck():
    angles = tuple(index * math.pi / (2.0 * FACETS) for index in range(FACETS + 1))
    nodes = []
    primary_map = {}
    secondary_map = {}
    primary = []
    secondary = []
    for index in range(FACETS):
        primary.append(append_annular(nodes, primary_map, 0.8, PRIMARY_RADIUS, angles[index], angles[index + 1]))
    for index in range(FACETS):
        secondary.append(
            append_annular(nodes, secondary_map, SECONDARY_RADIUS, 1.2, angles[index], angles[index + 1])
        )
    primary_node_count = 4 * (FACETS + 1)
    all_node_count = 2 * primary_node_count
    lines = [
        "*Heading",
        "** B5.52 lightweight nonuniform curved finite-sliding surface-to-surface operator reference.",
        "** Six C3D8T elements and eight secondary constraints; no equilibrium solve is compared.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Node",
    ]
    for label, point in enumerate(nodes, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + point))
    lines.append("*Element, type=C3D8T, elset=PRIMARY")
    for label, element in enumerate(primary, 1):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in element)))
    lines.append("*Element, type=C3D8T, elset=SECONDARY")
    for label, element in enumerate(secondary, FACETS + 1):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in element)))
    lines.extend([
        "*Nset, nset=ALL_NODES, generate", "1, %d, 1" % all_node_count,
        "*Nset, nset=PRIMARY_ALL, generate", "1, %d, 1" % primary_node_count,
        "*Nset, nset=SECONDARY_ALL, generate", "%d, %d, 1" % (primary_node_count + 1, all_node_count),
        "*Surface, type=ELEMENT, name=PRIMARY_CONTACT", "PRIMARY, S4",
        "*Surface, type=ELEMENT, name=SECONDARY_CONTACT", "SECONDARY, S6",
        "*Material, name=ELASTIC", "*Elastic", "1.0000000000000000e9, 0.0",
        "*Conductivity", "1.0000000000000000e1", "*Density", "1.0", "*Specific Heat", "1.0",
        "*Solid Section, elset=PRIMARY, material=ELASTIC", ",",
        "*Solid Section, elset=SECONDARY, material=ELASTIC", ",",
        "*Surface Interaction, name=CONTACT", "*Surface Behavior, pressure-overclosure=LINEAR",
        "1.0000000000000000e7",
        "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.",
        "SECONDARY_CONTACT, PRIMARY_CONTACT",
        "*Initial Conditions, type=TEMPERATURE", "ALL_NODES, 3.0000000000000000e2",
        "*Step, name=LOAD, nlgeom=YES, inc=20", "*Coupled Temperature-Displacement, steady state",
        "1.0, 1.0, 1.0e-8, 1.0", "*Boundary, op=NEW",
    ])
    for node in range(1, primary_node_count + 1):
        lines.append("%d, 1, 3, 0.0" % node)
        lines.append("%d, 11, 11, 3.0000000000000000e2" % node)
    for node in range(primary_node_count + 1, all_node_count + 1):
        point = nodes[node - 1]
        radius = math.hypot(point[0], point[1])
        angle = math.atan2(point[1], point[0])
        inward = inward_displacement(angle, point[2])
        displacement = (-inward * point[0] / radius, -inward * point[1] / radius, 0.0)
        for component in range(3):
            lines.append("%d, %d, %d, %.16e" % (node, component + 1, component + 1, displacement[component]))
        lines.append("%d, 11, 11, 3.0000000000000000e2" % node)
    lines.extend([
        "*Output, field, frequency=1", "*Node Output", "COORD, NT, RF, RFL, U",
        "*Contact Output", "CSTRESS, CDISP, CFORCE, HFL", "*End Step",
    ])
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    path = Path(__file__).resolve().parent / "b552_hex8_sts_cross_face.inp"
    path.write_text(deck(), encoding="utf-8")
    print("wrote %s" % path)
