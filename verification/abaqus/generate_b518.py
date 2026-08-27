#!/usr/bin/env python3
"""Generate the B5.18 finite-strain noncoaxial coupled plastic-creep C3D8T probe."""

from pathlib import Path

import generate_b513


if __name__ == "__main__":
    output = Path(__file__).with_name("b518_hex8_c3d8t_finite_noncoaxial.inp")
    contents = generate_b513.deck()
    contents = contents.replace(
        "B5.13 Abaqus/Standard small-strain noncoaxial",
        "B5.18 Abaqus/Standard finite-strain noncoaxial",
    ).replace("nlgeom=NO", "nlgeom=YES").replace("COORD, E, EE", "COORD, EE, LE")
    output.write_text(contents, encoding="utf-8")
    print("wrote %s" % output)
