#!/usr/bin/env python3
"""Generate the B5.2 fixed-gap C3D8T thermal-contact probe."""

from pathlib import Path


UNIT_COORDINATES = (
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
        "** B5.2 Abaqus/Standard C3D8T fixed-gap thermal-contact probe.",
        "** SI units: metre, second, kelvin, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for block, shift in enumerate((0.0, 1.1)):
        for local_node, coordinate in enumerate(UNIT_COORDINATES):
            node = 8 * block + local_node + 1
            lines.append(
                "%d, %.16e, %.16e, %.16e" % (node, coordinate[0] + shift, coordinate[1], coordinate[2])
            )
    lines.extend(
        [
            "*Element, type=C3D8T, elset=BLOCKS",
            "1, 1, 2, 3, 4, 5, 6, 7, 8",
            "2, 9, 10, 11, 12, 13, 14, 15, 16",
            "*Nset, nset=ALL_NODES, generate",
            "1, 16, 1",
            "*Nset, nset=HOT",
            "1, 4, 5, 8",
            "*Nset, nset=COLD",
            "10, 11, 14, 15",
            "*Surface, type=ELEMENT, name=SECONDARY_CONTACT",
            "1, S4",
            "*Surface, type=ELEMENT, name=PRIMARY_CONTACT",
            "2, S6",
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
            "*Solid Section, elset=BLOCKS, material=THERMAL",
            ",",
            "*Surface Interaction, name=THERMAL_GAP",
            "*Gap Conductance",
            "4.0000000000000000e1, 0.0000000000000000e0",
            "4.0000000000000000e1, 2.0000000000000001e-1",
            "*Contact Pair, interaction=THERMAL_GAP, type=SURFACE TO SURFACE, small sliding, adjust=0.",
            "SECONDARY_CONTACT, PRIMARY_CONTACT",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
            "*Step, name=HEAT, nlgeom=NO, inc=20",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0e-8, 1.0",
            "*Boundary",
            "ALL_NODES, 1, 3, 0.0",
            "HOT, 11, 11, 4.0000000000000000e2",
            "COLD, 11, 11, 3.0000000000000000e2",
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RFL, U",
            "*Element Output, directions=YES",
            "COORD, HFL",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b52_hex8_c3d8t_thermal_contact.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
