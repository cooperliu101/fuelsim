#!/usr/bin/env python3
"""Generate the B5.5 two-element transient C3D8T full-field probe."""

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


def deck():
    lines = [
        "*Heading",
        "** B5.5 Abaqus/Standard two-element transient C3D8T full-field probe.",
        "** SI units: metre, second, kelvin, watt, pascal.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for label, coordinate in enumerate(NODES, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + coordinate))
    lines.extend(["*Element, type=C3D8T, elset=VOLUME"])
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
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "1.2000000000000000e-5",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3",
            "*Specific Heat",
            "3.0000000000000000e3",
            "*Solid Section, elset=VOLUME, material=THERMAL_ELASTIC",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
            "*Step, name=TRANSIENT, nlgeom=NO, inc=4",
            "*Coupled Temperature-Displacement",
            "1.0, 1.0, 1.0, 1.0",
            "*Boundary",
            "LEFT, 11, 11, 3.0000000000000000e2",
            "RIGHT, 11, 11, 4.0000000000000000e2",
            "LEFT, 1, 1, 0.0",
            "Y0, 2, 2, 0.0",
            "Z0, 3, 3, 0.0",
            "*Dflux",
            "VOLUME, BF, 4.0000000000000000e6",
            "1, S2, 5.0000000000000000e3",
            "2, S2, 5.0000000000000000e3",
            "*Film",
            "1, F5, 2.8000000000000000e2, 2.0000000000000000e1",
            "2, F5, 2.8000000000000000e2, 2.0000000000000000e1",
            "*Output, field, frequency=1",
            "*Node Output",
            "NT, RFL, RF, U",
            "*Element Output, directions=YES",
            "COORD, E, HFL, IVOL, S, TEMP",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b55_hex8_c3d8t_transient_full_field.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
