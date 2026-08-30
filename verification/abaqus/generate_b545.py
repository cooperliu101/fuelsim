#!/usr/bin/env python3
"""Generate the fixed 100-step B5.45 finite-strain noncoaxial C3D8RT path."""

from pathlib import Path

import generate_b513


def axial(stage):
    one_based = float(stage + 1)
    if stage < 25:
        return 0.03 * one_based / 25.0
    return 0.03


def shear(stage):
    one_based = float(stage + 1)
    if stage < 25:
        return 0.0
    return (one_based - 25.0) / 75.0


if __name__ == "__main__":
    generate_b513.AXIAL = tuple(axial(stage) for stage in range(100))
    generate_b513.SHEAR = tuple(shear(stage) for stage in range(100))
    contents = generate_b513.deck()
    contents = contents.replace(
        "B5.13 Abaqus/Standard small-strain noncoaxial",
        "B5.45 Abaqus/Standard 100-step finite-strain noncoaxial",
    )
    contents = contents.replace("C3D8T probe", "C3D8RT probe")
    contents = contents.replace("type=C3D8T", "type=C3D8RT")
    contents = contents.replace("nlgeom=NO", "nlgeom=YES")
    contents = contents.replace("COORD, E, EE", "COORD, EE, LE")
    output = Path(__file__).with_name("b545_hex8_c3d8rt_finite_noncoaxial_100step.inp")
    output.write_text(contents, encoding="utf-8")
    print("wrote %s" % output)
