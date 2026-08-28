#!/usr/bin/env python3
"""Generate the B5.23 multi-element fully coupled C3D8T path."""

from pathlib import Path


TEMPERATURES = (350.0, 400.0, 450.0, 500.0)
PRESSURES = (8.75e4, 1.75e5, 2.625e5, 3.5e5)
TANGENTIAL_Y = (5.0e-3, 1.0e-2, 1.5e-2, 2.0e-2)
TANGENTIAL_Z = (2.5e-3, 5.0e-3, 7.5e-3, 1.0e-2)


def append_block(nodes, elements, x0, x1):
    node_map = {}
    x_values = (x0, 0.5 * (x0 + x1), x1)
    for z in (0.0, 1.0):
        for y in (0.0, 1.0):
            for x in x_values:
                node_map[(x, y, z)] = len(nodes) + 1
                nodes.append((x, y, z))
    for lower_x, upper_x in zip(x_values[:-1], x_values[1:]):
        elements.append(
            (
                node_map[(lower_x, 0.0, 0.0)],
                node_map[(upper_x, 0.0, 0.0)],
                node_map[(upper_x, 1.0, 0.0)],
                node_map[(lower_x, 1.0, 0.0)],
                node_map[(lower_x, 0.0, 1.0)],
                node_map[(upper_x, 0.0, 1.0)],
                node_map[(upper_x, 1.0, 1.0)],
                node_map[(lower_x, 1.0, 1.0)],
            )
        )
    return node_map


def labels_at(nodes, predicate):
    return [index for index, point in enumerate(nodes, 1) if predicate(*point)]


def append_set(lines, name, labels):
    lines.append("*Nset, nset=%s" % name)
    for begin in range(0, len(labels), 16):
        lines.append(", ".join(str(label) for label in labels[begin : begin + 16]))


def step_lines():
    return [
        "*Step, name=PATH, nlgeom=YES, inc=80",
        "*Coupled Temperature-Displacement",
        "2.0000000000000000e-2, 4.0000000000000002e-1, 2.0000000000000000e-2, 2.0000000000000000e-2",
        "*Controls, parameters=FIELD, field=DISPLACEMENT",
        "1.0000000000000000e-10, 1.0000000000000000e-10, , , , 1.0000000000000000e-10",
        "*Boundary",
        "PRIMARY_OUTER, 1, 3, 0.0",
        "PRIMARY_OUTER, 11, 11, 3.0000000000000000e2",
        "*Boundary, amplitude=TEMPERATURE_PATH",
        "SECONDARY_OUTER_NODES, 11, 11, 1.0",
        "*Boundary, amplitude=TANGENTIAL_Y_PATH",
        "SECONDARY_OUTER_NODES, 2, 2, 1.0",
        "*Boundary, amplitude=TANGENTIAL_Z_PATH",
        "SECONDARY_OUTER_NODES, 3, 3, 1.0",
        "*Dsload, follower=YES, amplitude=PRESSURE_PATH",
        "SECONDARY_OUTER, P, 1.0",
        "*Output, field, frequency=1",
        "*Node Output",
        "COORD, NT, RF, RFL, U",
        "*Element Output, directions=YES",
        "CE, CEEQ, COORD, EE, HFL, IVOL, LE, PE, PEEQ, S, TEMP",
        "*Contact Output",
        "CDISP, CFORCE, CSTRESS, HFL",
        "*Output, history, frequency=1",
        "*Energy Output",
        "ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK",
        "*End Step",
    ]


def deck():
    nodes = []
    primary_elements = []
    secondary_elements = []
    append_block(nodes, primary_elements, 0.0, 1.0)
    append_block(nodes, secondary_elements, 1.0, 2.0)
    primary_all = list(range(1, 13))
    secondary_all = list(range(13, 25))
    primary_outer = labels_at(nodes, lambda x, y, z: x == 0.0)
    secondary_outer = labels_at(nodes, lambda x, y, z: x == 2.0)
    secondary_contact = [label for label in secondary_all if nodes[label - 1][0] == 1.0]
    secondary_y0 = [label for label in secondary_all if nodes[label - 1][1] == 0.0]
    secondary_z0 = [label for label in secondary_all if nodes[label - 1][2] == 0.0]

    lines = [
        "*Heading",
        "** B5.23 multi-element transient C3D8T finite-strain thermo-inelastic frictional-contact path.",
        "** SI units: metre, second, kelvin, pascal, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Node",
    ]
    for label, point in enumerate(nodes, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + point))
    lines.append("*Element, type=C3D8T, elset=PRIMARY")
    for label, connectivity in enumerate(primary_elements, 1):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in connectivity)))
    lines.append("*Element, type=C3D8T, elset=SECONDARY")
    for label, connectivity in enumerate(secondary_elements, 3):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in connectivity)))
    lines.extend(
        [
            "*Elset, elset=PRIMARY_CONTACT_ELEMENTS",
            "2",
            "*Elset, elset=SECONDARY_CONTACT_ELEMENTS",
            "3",
            "*Elset, elset=SECONDARY_OUTER_ELEMENTS",
            "4",
        ]
    )
    append_set(lines, "ALL_NODES", list(range(1, 25)))
    append_set(lines, "PRIMARY_ALL", primary_all)
    append_set(lines, "SECONDARY_ALL", secondary_all)
    append_set(lines, "PRIMARY_OUTER", primary_outer)
    append_set(lines, "SECONDARY_OUTER_NODES", secondary_outer)
    append_set(lines, "SECONDARY_CONTACT_NODES", secondary_contact)
    append_set(lines, "SECONDARY_Y0", secondary_y0)
    append_set(lines, "SECONDARY_Z0", secondary_z0)
    lines.extend(
        [
            "*Surface, type=ELEMENT, name=PRIMARY_CONTACT",
            "PRIMARY_CONTACT_ELEMENTS, S4",
            "*Surface, type=ELEMENT, name=SECONDARY_CONTACT",
            "SECONDARY_CONTACT_ELEMENTS, S6",
            "*Surface, type=ELEMENT, name=SECONDARY_OUTER",
            "SECONDARY_OUTER_ELEMENTS, S4",
            "*Material, name=PRIMARY_MATERIAL",
            "*Elastic",
            "1.2000000000000000e8, 2.8000000000000003e-1, 3.0000000000000000e2",
            "1.0000000000000000e8, 2.8000000000000003e-1, 5.0000000000000000e2",
            "*Expansion, zero=3.0000000000000000e2",
            "8.0000000000000007e-6",
            "*Conductivity",
            "1.5000000000000000e1, 3.0000000000000000e2",
            "1.8000000000000000e1, 5.0000000000000000e2",
            "*Density",
            "1.0000000000000000e2",
            "*Specific Heat",
            "1.0000000000000000e0, 3.0000000000000000e2",
            "1.2000000000000000e0, 5.0000000000000000e2",
            "*Material, name=SECONDARY_MATERIAL",
            "*Elastic",
            "1.0000000000000000e8, 3.0000000000000000e-1, 3.0000000000000000e2",
            "9.0000000000000000e7, 3.0000000000000000e-1, 5.0000000000000000e2",
            "*Expansion, zero=3.0000000000000000e2",
            "1.0000000000000001e-5",
            "*Plastic",
            "2.0000000000000000e5, 0.0, 3.0000000000000000e2",
            "1.0200000000000000e7, 1.0, 3.0000000000000000e2",
            "1.8000000000000000e5, 0.0, 5.0000000000000000e2",
            "9.1800000000000000e6, 1.0, 5.0000000000000000e2",
            "*Creep, law=TIME",
            "1.2500000000000001e-20, 3.0, 0.0, 3.0000000000000000e2",
            "1.7500000000000000e-20, 3.0, 0.0, 5.0000000000000000e2",
            "*Conductivity",
            "1.0000000000000000e1, 3.0000000000000000e2",
            "1.2000000000000000e1, 5.0000000000000000e2",
            "*Density",
            "1.0000000000000000e2",
            "*Specific Heat",
            "1.0000000000000000e0, 3.0000000000000000e2",
            "1.2000000000000000e0, 5.0000000000000000e2",
            "*Solid Section, elset=PRIMARY, material=PRIMARY_MATERIAL",
            ",",
            "*Solid Section, elset=SECONDARY, material=SECONDARY_MATERIAL",
            ",",
            "*Surface Interaction, name=COUPLED_CONTACT",
            "*Surface Behavior, pressure-overclosure=LINEAR",
            "1.0000000000000000e9",
            "*Friction, slip tolerance=5.0000000000000001e-3",
            "5.0000000000000003e-2",
            "*Gap Conductance, pressure",
            "5.0000000000000000e1, 0.0",
            "5.5000000000000000e2, 5.0000000000000000e5",
            "*Contact Pair, interaction=COUPLED_CONTACT, type=SURFACE TO SURFACE, adjust=0.",
            "SECONDARY_CONTACT, PRIMARY_CONTACT",
            "*Amplitude, name=TEMPERATURE_PATH, time=TOTAL TIME",
            "0.0, 300.0, 0.1, 350.0, 0.2, 400.0, 0.3, 450.0",
            "0.4, 500.0",
            "*Amplitude, name=PRESSURE_PATH, time=TOTAL TIME",
            "0.0, 0.0, 0.1, 8.75e4, 0.2, 1.75e5, 0.3, 2.625e5",
            "0.4, 3.5e5",
            "*Amplitude, name=TANGENTIAL_Y_PATH, time=TOTAL TIME",
            "0.0, 0.0, 0.1, 5.0e-3, 0.2, 1.0e-2, 0.3, 1.5e-2",
            "0.4, 2.0e-2",
            "*Amplitude, name=TANGENTIAL_Z_PATH, time=TOTAL TIME",
            "0.0, 0.0, 0.1, 2.5e-3, 0.2, 5.0e-3, 0.3, 7.5e-3",
            "0.4, 1.0e-2",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
        ]
    )
    lines.extend(step_lines())
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b523_hex8_c3d8t_integrated_path.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
