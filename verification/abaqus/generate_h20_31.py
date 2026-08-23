from __future__ import print_function

import math
import sys


if len(sys.argv) != 3:
    raise RuntimeError("usage: generate_h20_31.py <h20_30_quadratic.inp> <h20_31_probe.inp>")


def parse_source(path):
    nodes = {}
    primary = None
    secondary = None
    mode = None
    source = open(path, "r")
    for raw_line in source:
        line = raw_line.strip()
        lower = line.lower()
        if lower == "*node":
            mode = "node"
            continue
        if lower.startswith("*element"):
            mode = "primary" if "elset=primary" in lower else "secondary"
            continue
        if line.startswith("*"):
            mode = None
            continue
        if not line:
            continue
        values = [value.strip() for value in line.split(",")]
        if mode == "node":
            nodes[int(values[0])] = tuple(float(value) for value in values[1:4])
        elif mode == "primary" and primary is None:
            primary = [int(value) for value in values]
        elif mode == "secondary" and secondary is None:
            secondary = [int(value) for value in values]
    source.close()
    if primary is None or secondary is None:
        raise RuntimeError("H20.31 could not read source elements")
    return nodes, primary, secondary


def write_set(output, name, labels):
    output.write("*Nset, nset=%s\n" % name)
    for first in range(0, len(labels), 16):
        output.write(", ".join(str(value) for value in labels[first:first + 16]) + "\n")


nodes, primary, secondary = parse_source(sys.argv[1])
primary_connectivity = primary[1:]
secondary_connectivity = secondary[1:]
# Abaqus C3D20 face S6 is (1, 4, 8, 5, 12, 20, 16, 17).
secondary_face_indices = [0, 3, 7, 4, 11, 19, 15, 16]
secondary_face = [secondary_connectivity[index] for index in secondary_face_indices]
needed = sorted(set(primary_connectivity + secondary_connectivity))
normal = {}
for label in secondary_face:
    coordinates = nodes[label]
    radius = math.hypot(coordinates[0], coordinates[1])
    normal[label] = (coordinates[0] / radius, coordinates[1] / radius)

output = open(sys.argv[2], "w")
output.write(
    "*Heading\n"
    "** H20.31 Abaqus curved C3D20 secondary operator probe.\n"
    "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
    "*Node\n"
)
for label in needed:
    coordinates = nodes[label]
    output.write("%d, %.16g, %.16g, %.16g\n" % ((label,) + coordinates))
output.write("*Element, type=C3D20, elset=PRIMARY\n%s\n" % ", ".join(str(value) for value in primary))
output.write("*Element, type=C3D20, elset=SECONDARY\n%s\n" % ", ".join(str(value) for value in secondary))
write_set(output, "PRIMARY_ALL", primary_connectivity)
write_set(output, "SECONDARY_ALL", secondary_connectivity)
for index, label in enumerate(secondary_face, 1):
    write_set(output, "S%d" % index, [label])
output.write(
    "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
    "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
    "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.25\n"
    "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
    "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
    "*Surface Interaction, name=LINEAR_PENALTY\n"
    "*Surface Behavior, pressure-overclosure=LINEAR\n1.e8,\n"
    "*Contact Pair, interaction=LINEAR_PENALTY, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
    "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
)


def write_step(name, perturbed, component, delta):
    output.write("*Step, name=%s, nlgeom=NO, inc=40\n*Static\n0.1, 1., 1.e-8, 0.1\n" % name)
    output.write("*Boundary%s\n" % ("" if name == "BASE" else ", op=MOD"))
    if name == "BASE":
        output.write("PRIMARY_ALL, 1, 3, 0.\n")
    for index, label in enumerate(secondary_face, 1):
        displacement = [-1.0e-4 * normal[label][0], -1.0e-4 * normal[label][1], 0.0]
        if index == perturbed:
            displacement[component] += delta
        for direction in range(3):
            output.write("S%d, %d, %d, %.16g\n" % (index, direction + 1, direction + 1,
                                                    displacement[direction]))
    output.write(
        "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
        "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n"
    )


write_step("BASE", 0, 0, 0.0)
for local_node in range(1, 9):
    for component, name in enumerate(("X", "Y", "Z")):
        write_step("S%d_%s_PLUS" % (local_node, name), local_node, component, 1.0e-6)
        write_step("S%d_%s_MINUS" % (local_node, name), local_node, component, -1.0e-6)
output.close()
