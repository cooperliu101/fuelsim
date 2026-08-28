#!/usr/bin/env python3
"""Generate the B5.22 faceted-cylinder C3D8T thermal-contact probe."""

import math
from pathlib import Path


def cylindrical(radius, angle, z):
    return (radius * math.cos(angle), radius * math.sin(angle), z)


def append_annular(nodes, node_map, r0, r1, a0, a1):
    logical = (
        (r0, a0, 0.0),
        (r1, a0, 0.0),
        (r1, a1, 0.0),
        (r0, a1, 0.0),
        (r0, a0, 1.0),
        (r1, a0, 1.0),
        (r1, a1, 1.0),
        (r0, a1, 1.0),
    )
    element = []
    for point in logical:
        if point not in node_map:
            node_map[point] = len(nodes) + 1
            nodes.append(cylindrical(*point))
        element.append(node_map[point])
    return element


def deck(facets=3, primary_radius=1.0, secondary_radius=1.005, swap=False):
    angles = tuple(index * math.pi / (2.0 * facets) for index in range(facets + 1))
    nodes = []
    primary_map = {}
    secondary_map = {}
    primary = []
    secondary = []
    for index in range(facets):
        primary.append(append_annular(
            nodes, primary_map, primary_radius - 0.2, primary_radius, angles[index], angles[index + 1]))
    for index in range(facets):
        secondary.append(append_annular(
            nodes, secondary_map, secondary_radius, secondary_radius + 0.195, angles[index], angles[index + 1]))
    primary_node_count = 4 * (facets + 1)
    all_node_count = 2 * primary_node_count
    lines = [
        "*Heading",
        "** B5.22 faceted-cylinder C3D8T pressure-dependent thermal contact.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Node",
    ]
    for label, point in enumerate(nodes, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + point))
    lines.append("*Element, type=C3D8T, elset=PRIMARY")
    for label, element in enumerate(primary, 1):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in element)))
    lines.append("*Element, type=C3D8T, elset=SECONDARY")
    for label, element in enumerate(secondary, facets + 1):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in element)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, %d, 1" % all_node_count,
            "*Nset, nset=PRIMARY_ALL, generate",
            "1, %d, 1" % primary_node_count,
            "*Nset, nset=SECONDARY_ALL, generate",
            "%d, %d, 1" % (primary_node_count + 1, all_node_count),
            "*Surface, type=ELEMENT, name=PRIMARY_CONTACT",
            "PRIMARY, S4",
            "*Surface, type=ELEMENT, name=SECONDARY_CONTACT",
            "SECONDARY, S6",
            "*Material, name=THERMOELASTIC",
            "*Elastic",
            "1.0000000000000000e9, 0.0",
            "*Expansion, zero=3.0000000000000000e2",
            "0.0",
            "*Conductivity",
            "1.0000000000000000e1",
            "*Density",
            "1.0",
            "*Specific Heat",
            "1.0",
            "*Solid Section, elset=PRIMARY, material=THERMOELASTIC",
            ",",
            "*Solid Section, elset=SECONDARY, material=THERMOELASTIC",
            ",",
            "*Surface Interaction, name=CONTACT",
            "*Surface Behavior, pressure-overclosure=LINEAR",
            "1.0000000000000000e7",
            "*Friction",
            "1.0000000000000001e-1",
            "*Gap Conductance, pressure",
            "0.0000000000000000e0, 0.0000000000000000e0",
            "2.0000000000000000e2, 2.0000000000000000e5",
            "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.",
            "PRIMARY_CONTACT, SECONDARY_CONTACT" if swap else "SECONDARY_CONTACT, PRIMARY_CONTACT",
            "*Initial Conditions, type=TEMPERATURE",
            "PRIMARY_ALL, 3.0000000000000000e2",
            "SECONDARY_ALL, 4.0000000000000000e2",
            "*Step, name=LOAD, nlgeom=YES, inc=20",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0e-8, 1.0",
            "*Boundary, op=NEW",
        ]
    )
    for node in range(1, primary_node_count + 1):
        lines.append("%d, 1, 3, 0.0" % node)
        lines.append("%d, 11, 11, 3.0000000000000000e2" % node)
    for node in range(primary_node_count + 1, all_node_count + 1):
        point = nodes[node - 1]
        radius = math.hypot(point[0], point[1])
        displacement = (-0.015 * point[0] / radius, -0.015 * point[1] / radius, 0.0)
        for component in range(3):
            lines.append("%d, %d, %d, %.16e" % (node, component + 1, component + 1, displacement[component]))
        lines.append("%d, 11, 11, 4.0000000000000000e2" % node)
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Contact Output",
            "CSTRESS, CDISP, CFORCE, HFL",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    cases = (
        ("b522_hex8_c3d8t_faceted_thermal_contact", 3, 1.0, 1.005, False),
        ("b522_hex8_c3d8t_faceted_thermal_contact_f6", 6, 1.0, 1.005, False),
        ("b522_hex8_c3d8t_faceted_thermal_contact_f12", 12, 1.0, 1.005, False),
        ("b522_hex8_c3d8t_faceted_thermal_contact_swapped", 3, 1.0, 1.005, True),
        ("b522_hex8_c3d8t_faceted_thermal_contact_tight", 3, 0.5, 0.505, False),
    )
    for name, facets, primary_radius, secondary_radius, swap in cases:
        output = Path(__file__).with_name(name + ".inp")
        output.write_text(deck(facets, primary_radius, secondary_radius, swap), encoding="utf-8")
        print("wrote %s" % output)
