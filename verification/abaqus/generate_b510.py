#!/usr/bin/env python3
"""Generate the B5.10 small-strain J2 plasticity C3D8T path probe."""

from pathlib import Path


DISPLACEMENTS = (0.0005, 0.002, 0.004, 0.002, 0.0, -0.002, -0.004, -0.001, 0.003, 0.0)


def step_lines(index, displacement):
    return [
        "*Step, name=STEP%d, nlgeom=NO, inc=8" % index,
        "*Coupled Temperature-Displacement",
        "1.0000000000000001e-1, 1.0000000000000001e-1, 1.0000000000000001e-1, 1.0000000000000001e-1",
        "*Boundary, op=NEW",
        "LEFT, 11, 11, 6.0000000000000000e2",
        "LEFT, 1, 1, 0.0",
        "Y0, 2, 2, 0.0",
        "Z0, 3, 3, 0.0",
        "RIGHT, 1, 1, %.16e" % displacement,
        "*Output, field, frequency=1",
        "*Node Output",
        "NT, RFL, RF, U",
        "*Element Output, directions=YES",
        "COORD, E, EE, IVOL, PE, PEEQ, S, TEMP",
        "*Output, history, frequency=1",
        "*Energy Output",
        "ALLIE, ALLPD, ALLWK",
        "*End Step",
    ]


def deck():
    lines = [
        "*Heading",
        "** B5.10 Abaqus/Standard small-strain J2 plasticity C3D8T load/unload/reversal probe.",
        "** Plastic dissipation heat generation is intentionally disabled.",
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
        "*Material, name=J2_MATERIAL",
        "*Elastic",
        "2.0000000000000000e11, 3.0000000000000000e-1",
        "*Plastic",
        "2.0000000000000000e8, 0.0",
        "2.2000000000000000e9, 1.0000000000000000e0",
        "*Conductivity",
        "1.0",
        "*Density",
        "1.0",
        "*Specific Heat",
        "1.0",
        "*Solid Section, elset=VOLUME, material=J2_MATERIAL",
        ",",
        "*Initial Conditions, type=TEMPERATURE",
        "ALL_NODES, 6.0000000000000000e2",
    ]
    for index, displacement in enumerate(DISPLACEMENTS, 1):
        lines.extend(step_lines(index, displacement))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b510_hex8_c3d8t_small_j2.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
