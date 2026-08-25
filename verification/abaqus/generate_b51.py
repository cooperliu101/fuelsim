#!/usr/bin/env python3
"""Generate the B5.1 finite-deformation C3D8T heat-configuration probe."""

from pathlib import Path


COORDINATES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)


def deck():
    lines = [
        "*Heading",
        "** B5.1 Abaqus/Standard C3D8T finite-deformation heat-configuration probe.",
        "** SI units: metre, second, kelvin, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for node, coordinate in enumerate(COORDINATES, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((node,) + coordinate))
    lines.extend(
        [
            "*Element, type=C3D8T, elset=VOLUME",
            "1, 1, 2, 3, 4, 5, 6, 7, 8",
            "*Nset, nset=ALL_NODES, generate",
            "1, 8, 1",
            "*Nset, nset=COLD",
            "1, 4, 5, 8",
            "*Nset, nset=HOT",
            "2, 3, 6, 7",
            "*Material, name=THERMAL",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "0.0",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3",
            "*Specific Heat",
            "3.0000000000000000e3",
            "*Solid Section, elset=VOLUME, material=THERMAL",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
            "*Step, name=DEFORM, nlgeom=YES, inc=10",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0, 1.0",
            "*Boundary",
            "ALL_NODES, 11, 11, 3.0000000000000000e2",
        ]
    )
    for node, coordinate in enumerate(COORDINATES, 1):
        displacement = (0.5 * coordinate[0], 0.25 * coordinate[1], -0.2 * coordinate[2])
        for component, value in enumerate(displacement, 1):
            lines.append("%d, %d, %d, %.16e" % (node, component, component, value))
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RFL, RF, U",
            "*Element Output, directions=YES",
            "COORD, HFL",
            "*End Step",
            "*Step, name=HEAT, nlgeom=YES, inc=10",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0, 1.0",
            "*Boundary",
            "COLD, 11, 11, 3.0000000000000000e2",
            "HOT, 11, 11, 4.0000000000000000e2",
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RFL, RF, U",
            "*Element Output, directions=YES",
            "COORD, HFL",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b51_hex8_c3d8t_finite_heat_probe.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
