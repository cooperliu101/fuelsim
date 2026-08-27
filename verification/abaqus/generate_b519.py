#!/usr/bin/env python3
"""Generate the B5.19 finite-strain C3D8T selective-volume probe."""

from pathlib import Path

import generate_b49


SCALE = 200.0
TEMPERATURE = (600.0,) * 8
DISPLACEMENT_X = tuple(SCALE * value for value in generate_b49.DISPLACEMENT_X)
DISPLACEMENT_Y = tuple(SCALE * value for value in generate_b49.DISPLACEMENT_Y)
DISPLACEMENT_Z = tuple(SCALE * value for value in generate_b49.DISPLACEMENT_Z)
STATE = TEMPERATURE + DISPLACEMENT_X + DISPLACEMENT_Y + DISPLACEMENT_Z


def deck():
    lines = [
        "*Heading",
        "** B5.19 Abaqus/Standard finite-strain C3D8T selective-volume probe.",
        "** SI units: metre, kelvin, newton, pascal.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for label, coordinate in enumerate(generate_b49.COORDINATES, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + coordinate))
    lines.extend(
        [
            "*Element, type=C3D8T, elset=VOLUME",
            "1, 1, 2, 3, 4, 5, 6, 7, 8",
            "*Nset, nset=ALL_NODES, generate",
            "1, 8, 1",
            "*Material, name=ELASTIC",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Conductivity",
            "1.0",
            "*Density",
            "1.0",
            "*Specific Heat",
            "1.0",
            "*Solid Section, elset=VOLUME, material=ELASTIC",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 6.0000000000000000e2",
            "*Step, name=PROBE, nlgeom=YES, inc=1",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0, 1.0",
            "*Boundary, op=NEW",
        ]
    )
    lines.extend(generate_b49.boundary_lines(STATE))
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Element Output, directions=YES",
            "COORD, E, IVOL, LE, S, TEMP",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b519_hex8_c3d8t_finite_selective.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
