#!/usr/bin/env python3
"""Generate the B5.11 small-strain Norton creep C3D8T path probe."""

from pathlib import Path


STRESS = 3.0e8


def step_lines(index):
    return [
        "*Step, name=STEP%d, nlgeom=NO, inc=64" % index,
        "*Coupled Temperature-Displacement",
        "1.0000000000000001e-1, 1.0000000000000001e-1, 1.0000000000000001e-1, 1.0000000000000001e-1",
        "*Boundary, op=NEW",
        "LEFT, 11, 11, 6.0000000000000000e2",
        "LEFT, 1, 1, 0.0",
        "Y0, 2, 2, 0.0",
        "Z0, 3, 3, 0.0",
        "*Cload, op=NEW, amplitude=CONSTANT",
        "2, 1, %.16e" % (0.25 * STRESS),
        "3, 1, %.16e" % (0.25 * STRESS),
        "6, 1, %.16e" % (0.25 * STRESS),
        "7, 1, %.16e" % (0.25 * STRESS),
        "*Output, field, frequency=1",
        "*Node Output",
        "NT, RFL, RF, U",
        "*Element Output, directions=YES",
        "CE, CEEQ, COORD, E, EE, IVOL, S, TEMP",
        "*Output, history, frequency=1",
        "*Energy Output",
        "ALLCD, ALLIE, ALLWK",
        "*End Step",
    ]


def deck():
    lines = [
        "*Heading",
        "** B5.11 Abaqus/Standard small-strain Norton creep C3D8T constant-stress probe.",
        "** Creep dissipation heat generation is intentionally disabled.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
        "1, 0.0, 0.0, 0.0",
        "2, 1.0, 0.0, 0.0",
        "3, 1.0, 1.0, 0.0",
        "4, 0.0, 1.0, 0.0",
        "5, 0.0, 0.0, 1.0",
        "6, 1.0, 0.0, 1.0",
        "7, 1.0, 1.0, 1.0",
        "8, 0.0, 1.0, 1.0",
        "*Element, type=C3D8T, elset=VOLUME",
        "1, 1, 2, 3, 4, 5, 6, 7, 8",
        "*Nset, nset=ALL_NODES, generate",
        "1, 8, 1",
        "*Nset, nset=LEFT",
        "1, 4, 5, 8",
        "*Nset, nset=RIGHT",
        "2, 3, 6, 7",
        "*Nset, nset=Y0",
        "1, 2, 5, 6",
        "*Nset, nset=Z0, generate",
        "1, 4, 1",
        "*Material, name=NORTON_MATERIAL",
        "*Elastic",
        "2.0000000000000000e11, 3.0000000000000000e-1",
        "*Creep, law=TIME",
        "1.0000000000000001e-28, 3.0000000000000000e0, 0.0",
        "*Conductivity",
        "1.0",
        "*Density",
        "1.0",
        "*Specific Heat",
        "1.0",
        "*Solid Section, elset=VOLUME, material=NORTON_MATERIAL",
        ",",
        "*Initial Conditions, type=TEMPERATURE",
        "ALL_NODES, 6.0000000000000000e2",
        "*Amplitude, name=CONSTANT, time=STEP TIME",
        "0.0, 1.0, 1.0, 1.0",
        "*Step, name=PRELOAD, nlgeom=NO",
        "*Static",
        "1.0, 1.0, 1.0000000000000000e-5, 1.0",
        "*Boundary",
        "LEFT, 1, 1, 0.0",
        "Y0, 2, 2, 0.0",
        "Z0, 3, 3, 0.0",
        "*Cload, amplitude=CONSTANT",
        "2, 1, %.16e" % (0.25 * STRESS),
        "3, 1, %.16e" % (0.25 * STRESS),
        "6, 1, %.16e" % (0.25 * STRESS),
        "7, 1, %.16e" % (0.25 * STRESS),
        "*End Step",
    ]
    for index in range(1, 11):
        lines.extend(step_lines(index))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b511_hex8_c3d8t_small_norton.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
