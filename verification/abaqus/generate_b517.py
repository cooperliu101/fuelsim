#!/usr/bin/env python3
"""Generate the B5.17 finite-strain coupled plastic-creep C3D8T probe."""

from pathlib import Path

import generate_b512


if __name__ == "__main__":
    output = Path(__file__).with_name("b517_hex8_c3d8t_finite_coupled.inp")
    contents = generate_b512.deck()
    contents = contents.replace(
        "B5.12 Abaqus/Standard small-strain active coupled",
        "B5.17 Abaqus/Standard finite-strain active coupled",
    ).replace("nlgeom=NO", "nlgeom=YES").replace("COORD, E, EE", "COORD, EE, LE").replace(
        "*Boundary, op=NEW",
        "*Controls, parameters=FIELD, field=DISPLACEMENT\n"
        "1.0000000000000000e-12, 1.0000000000000000e-12, , , , 1.0000000000000000e-12\n"
        "*Boundary, op=NEW",
    )
    output.write_text(contents, encoding="utf-8")
    print("wrote %s" % output)
