#!/usr/bin/env python3
"""Generate the B5.8 multi-step temperature-dependent C3D8T full-field probe."""

from pathlib import Path


NODES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (1.0, 1.0, 0.0),
    (2.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (2.0, 0.0, 1.0),
    (0.0, 1.0, 1.0),
    (1.0, 1.0, 1.0),
    (2.0, 1.0, 1.0),
)
ELEMENTS = ((1, 2, 5, 4, 7, 8, 11, 10), (2, 3, 6, 5, 8, 9, 12, 11))


def step_lines(index, right_temperature, body_heat, surface_heat):
    return [
        "*Step, name=STEP%d, nlgeom=NO, inc=4" % index,
        "*Coupled Temperature-Displacement",
        "1.0, 1.0, 1.0, 1.0",
        "*Boundary, op=NEW",
        "LEFT, 11, 11, 3.0000000000000000e2",
        "RIGHT, 11, 11, %.16e" % right_temperature,
        "LEFT, 1, 1, 0.0",
        "Y0, 2, 2, 0.0",
        "Z0, 3, 3, 0.0",
        "*Dflux, op=NEW",
        "VOLUME, BF, %.16e" % body_heat,
        "1, S2, %.16e" % surface_heat,
        "2, S2, %.16e" % surface_heat,
        "*Film, op=NEW",
        "1, F5, 2.8000000000000000e2, 2.0000000000000000e1",
        "2, F5, 2.8000000000000000e2, 2.0000000000000000e1",
        "*Output, field, frequency=1",
        "*Node Output",
        "NT, RFL, RF, U",
        "*Element Output, directions=YES",
        "COORD, E, HFL, IVOL, S, TEMP",
        "*End Step",
    ]


def deck():
    lines = [
        "*Heading",
        "** B5.8 Abaqus/Standard multi-step temperature-dependent C3D8T full-field probe.",
        "** The four one-second increments cover heating, continued heating, holding and cooling.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for label, coordinate in enumerate(NODES, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + coordinate))
    lines.append("*Element, type=C3D8T, elset=VOLUME")
    for label, connectivity in enumerate(ELEMENTS, 1):
        lines.append("%d, %s" % (label, ", ".join(str(node) for node in connectivity)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, 12, 1",
            "*Nset, nset=LEFT",
            "1, 4, 7, 10",
            "*Nset, nset=RIGHT",
            "3, 6, 9, 12",
            "*Nset, nset=Y0",
            "1, 2, 3, 7, 8, 9",
            "*Nset, nset=Z0, generate",
            "1, 6, 1",
            "*Material, name=THERMAL_ELASTIC",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1, 3.0000000000000000e2",
            "1.7000000000000000e11, 2.8000000000000003e-1, 6.0000000000000000e2",
            "*Expansion, zero=3.0000000000000000e2",
            "1.2000000000000000e-5, 3.0000000000000000e2",
            "1.8000000000000000e-5, 6.0000000000000000e2",
            "*Conductivity",
            "4.0000000000000000e0, 3.0000000000000000e2",
            "7.0000000000000000e0, 6.0000000000000000e2",
            "*Density",
            "2.0000000000000000e3, 3.0000000000000000e2",
            "1.7000000000000000e3, 6.0000000000000000e2",
            "*Specific Heat",
            "3.0000000000000000e3, 3.0000000000000000e2",
            "4.2000000000000000e3, 6.0000000000000000e2",
            "*Solid Section, elset=VOLUME, material=THERMAL_ELASTIC",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
        ]
    )
    paths = (
        (360.0, 2.0e6, 2.0e3),
        (420.0, 4.0e6, 5.0e3),
        (420.0, 0.0, 5.0e3),
        (330.0, 0.0, 0.0),
    )
    for index, values in enumerate(paths, 1):
        lines.extend(step_lines(index, *values))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b58_hex8_c3d8t_multistep.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
