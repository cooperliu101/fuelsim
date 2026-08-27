#!/usr/bin/env python3
"""Generate the B5.15 finite-strain J2 plasticity C3D8T path probe."""

from pathlib import Path

import generate_b510


if __name__ == "__main__":
    output = Path(__file__).with_name("b515_hex8_c3d8t_finite_j2.inp")
    contents = generate_b510.deck()
    contents = contents.replace(
        "B5.10 Abaqus/Standard small-strain J2 plasticity",
        "B5.15 Abaqus/Standard finite-strain J2 plasticity",
    ).replace("nlgeom=NO", "nlgeom=YES").replace("COORD, E, EE", "COORD, EE, LE")
    output.write_text(contents, encoding="utf-8")
    print("wrote %s" % output)
