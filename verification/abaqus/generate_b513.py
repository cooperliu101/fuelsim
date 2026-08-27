#!/usr/bin/env python3
"""Generate the B5.13 small-strain noncoaxial coupled plastic-creep C3D8T probe."""

from pathlib import Path


AXIAL = tuple(0.00075 + 0.00045 * i for i in range(1, 6))
AXIAL += (0.003,) * 5
AXIAL += tuple(0.003 - 0.001 * i for i in range(1, 6))
AXIAL += (-0.002,) * 5
SHEAR = (0.0,) * 5
SHEAR += tuple(0.0008 * i for i in range(1, 6))
SHEAR += (0.004,) * 5
SHEAR += tuple(0.004 - 0.0016 * i for i in range(1, 6))


def step_lines(index, axial, shear):
    return [
        "*Step, name=STEP%d, nlgeom=NO, inc=32" % index,
        "*Coupled Temperature-Displacement",
        "1.0000000000000000e-3, 1.0000000000000000e-3, 1.0000000000000000e-3, 1.0000000000000000e-3",
        "*Controls, parameters=FIELD, field=DISPLACEMENT",
        "1.0000000000000000e-12, 1.0000000000000000e-12, , , , 1.0000000000000000e-12",
        "*Boundary, op=NEW",
        "ALL_NODES, 11, 11, 6.0000000000000000e2",
        "LEFT, 1, 2, 0.0",
        "RIGHT, 1, 1, %.16e" % axial,
        "RIGHT, 2, 2, %.16e" % shear,
        "Z0, 3, 3, 0.0",
        "*Output, field, frequency=1",
        "*Node Output",
        "NT, RFL, RF, U",
        "*Element Output, directions=YES",
        "CE, CEEQ, COORD, E, EE, IVOL, PE, PEEQ, S, TEMP",
        "*Output, history, frequency=1",
        "*Energy Output",
        "ALLCD, ALLIE, ALLPD, ALLWK",
        "*End Step",
    ]


def deck():
    lines = [
        "*Heading",
        "** B5.13 Abaqus/Standard small-strain noncoaxial coupled plastic-creep C3D8T probe.",
        "** Inelastic dissipation heat generation is intentionally disabled.",
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
        "*Nset, nset=Z0, generate",
        "1, 4, 1",
        "*Material, name=NONCOAXIAL_MATERIAL",
        "*Elastic",
        "2.0000000000000000e11, 3.0000000000000000e-1",
        "*Plastic",
        "2.0000000000000000e8, 0.0",
        "2.2000000000000000e9, 1.0000000000000000e0",
        "*Creep, law=TIME",
        "1.0000000000000001e-29, 3.0000000000000000e0, 0.0",
        "*Conductivity",
        "1.0",
        "*Density",
        "1.0",
        "*Specific Heat",
        "1.0",
        "*Solid Section, elset=VOLUME, material=NONCOAXIAL_MATERIAL",
        ",",
        "*Initial Conditions, type=TEMPERATURE",
        "ALL_NODES, 6.0000000000000000e2",
    ]
    for index, (axial, shear) in enumerate(zip(AXIAL, SHEAR), 1):
        lines.extend(step_lines(index, axial, shear))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b513_hex8_c3d8t_small_noncoaxial.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
