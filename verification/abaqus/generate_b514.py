#!/usr/bin/env python3
"""Generate the B5.14 finite-strain elastic C3D8T kinematics probe."""

import math
from pathlib import Path


NODES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)


def path():
    result = []
    for index in range(1, 6):
        result.append((1.0 + 0.006 * index, 0.0, 0.0))
    for index in range(1, 6):
        result.append((1.03, 0.008 * index, 0.0))
    for index in range(1, 6):
        result.append((1.03, 0.04, math.radians(6.0 * index)))
    for index in range(1, 6):
        result.append((1.03, 0.04, math.radians(30.0 - 6.0 * index)))
    return tuple(result)


def displacement(node, stretch, shear, rotation):
    x, y, z = node
    base_x = stretch * x + shear * y
    base_y = y
    cosine = math.cos(rotation)
    sine = math.sin(rotation)
    current_x = cosine * base_x - sine * base_y
    current_y = sine * base_x + cosine * base_y
    return current_x - x, current_y - y, 0.0 * z


def step_lines(index, deformation):
    stretch, shear, rotation = deformation
    lines = [
        "*Step, name=STEP%d, nlgeom=YES, inc=32" % index,
        "*Coupled Temperature-Displacement",
        "1.0000000000000000e-3, 1.0000000000000000e-3, 1.0000000000000000e-3, 1.0000000000000000e-3",
        "*Controls, parameters=FIELD, field=DISPLACEMENT",
        "1.0000000000000000e-12, 1.0000000000000000e-12, , , , 1.0000000000000000e-12",
        "*Boundary, op=NEW",
        "ALL_NODES, 11, 11, 6.0000000000000000e2",
    ]
    for label, node in enumerate(NODES, 1):
        ux, uy, uz = displacement(node, stretch, shear, rotation)
        lines.extend(
            (
                "N%d, 1, 1, %.16e" % (label, ux),
                "N%d, 2, 2, %.16e" % (label, uy),
                "N%d, 3, 3, %.16e" % (label, uz),
            )
        )
    lines.extend(
        (
            "*Output, field, frequency=1",
            "*Node Output",
            "NT, RFL, RF, U",
            "*Element Output, directions=YES",
            "COORD, E, EE, IVOL, LE, S, TEMP",
            "*Output, history, frequency=1",
            "*Energy Output",
            "ALLIE, ALLSE, ALLWK",
            "*End Step",
        )
    )
    return lines


def deck():
    lines = [
        "*Heading",
        "** B5.14 Abaqus/Standard finite-strain elastic C3D8T kinematics probe.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for label, node in enumerate(NODES, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + node))
    lines.extend(("*Element, type=C3D8T, elset=VOLUME", "1, 1, 2, 3, 4, 5, 6, 7, 8", "*Nset, nset=ALL_NODES, generate", "1, 8, 1"))
    for label in range(1, 9):
        lines.extend(("*Nset, nset=N%d" % label, "%d" % label))
    lines.extend(
        (
            "*Material, name=FINITE_ELASTIC",
            "*Elastic",
            "2.0000000000000000e11, 3.0000000000000000e-1",
            "*Conductivity",
            "1.0",
            "*Density",
            "1.0",
            "*Specific Heat",
            "1.0",
            "*Solid Section, elset=VOLUME, material=FINITE_ELASTIC",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 6.0000000000000000e2",
        )
    )
    for index, deformation in enumerate(path(), 1):
        lines.extend(step_lines(index, deformation))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b514_hex8_c3d8t_finite_elastic.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
