#!/usr/bin/env python3
"""Generate B5.21 C3D8T release, sliding, and recontact thermal-contact probe."""

from pathlib import Path


PATH = (
    ("CLOSE", (-0.01, 0.05, 0.04)),
    ("SLIDE", (-0.01, 0.95, 0.20)),
    ("OPEN", (0.02, 0.95, 0.20)),
    ("OPEN_CROSS", (0.02, 0.05, 0.40)),
    ("RECONTACT", (-0.01, 0.05, 0.40)),
)


def append_cuboid(nodes, node_map, bounds):
    x0, x1, y0, y1, z0, z1 = bounds
    points = (
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    )
    element = []
    for point in points:
        if point not in node_map:
            node_map[point] = len(nodes) + 1
            nodes.append(point)
        element.append(node_map[point])
    return element


def deck():
    nodes = []
    primary_nodes = {}
    secondary_nodes = {}
    elements = (
        append_cuboid(nodes, primary_nodes, (0.0, 1.0, 0.0, 1.0, -1.0, 2.0)),
        append_cuboid(nodes, primary_nodes, (0.0, 1.0, 1.0, 2.0, -1.0, 2.0)),
        append_cuboid(nodes, secondary_nodes, (1.0, 2.0, 0.1, 0.9, 0.1, 0.9)),
    )
    lines = [
        "*Heading",
        "** B5.21 C3D8T pressure-dependent thermal contact through close, slide, open, and recontact.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Node",
    ]
    for label, point in enumerate(nodes, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + point))
    lines.extend(["*Element, type=C3D8T, elset=PRIMARY"])
    for label in (1, 2):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in elements[label - 1])))
    lines.extend(
        [
            "*Element, type=C3D8T, elset=SECONDARY",
            "3, %s" % ", ".join(str(node) for node in elements[2]),
            "*Nset, nset=ALL_NODES, generate",
            "1, %d, 1" % len(nodes),
            "*Nset, nset=PRIMARY_ALL",
            ", ".join(str(node) for node in sorted(primary_nodes.values())),
            "*Nset, nset=SECONDARY_ALL",
            ", ".join(str(node) for node in sorted(secondary_nodes.values())),
            "*Surface, type=ELEMENT, name=PRIMARY_CONTACT",
            "PRIMARY, S4",
            "*Surface, type=ELEMENT, name=SECONDARY_CONTACT",
            "SECONDARY, S6",
            "*Material, name=THERMOELASTIC",
            "*Elastic",
            "1.0000000000000000e9, 2.5000000000000000e-1",
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
            "1.0000000000000000e5",
            "*Friction, slip tolerance=0.025",
            "3.0000000000000000e-1",
            "*Gap Conductance, pressure",
            "0.0000000000000000e0, 0.0000000000000000e0",
            "4.0000000000000000e1, 2.0000000000000000e3",
            "*Contact Pair, interaction=CONTACT, type=SURFACE TO SURFACE, adjust=0.",
            "SECONDARY_CONTACT, PRIMARY_CONTACT",
            "*Initial Conditions, type=TEMPERATURE",
            "PRIMARY_ALL, 3.0000000000000000e2",
            "SECONDARY_ALL, 4.0000000000000000e2",
        ]
    )
    for index, (name, displacement) in enumerate(PATH):
        lines.extend(
            [
                "*Step, name=%s, nlgeom=YES, inc=100" % name,
                "*Coupled Temperature-Displacement, steady state",
                "0.1, 1.0, 1.0e-8, 0.1",
                "*Boundary, op=NEW",
                "PRIMARY_ALL, 1, 3, 0.0",
                "PRIMARY_ALL, 11, 11, 3.0000000000000000e2",
                "SECONDARY_ALL, 1, 1, %.16e" % displacement[0],
                "SECONDARY_ALL, 2, 2, %.16e" % displacement[1],
                "SECONDARY_ALL, 3, 3, %.16e" % displacement[2],
                "SECONDARY_ALL, 11, 11, 4.0000000000000000e2",
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
    output = Path(__file__).with_name("b521_hex8_c3d8t_thermal_contact_path.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
