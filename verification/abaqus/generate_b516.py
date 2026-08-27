#!/usr/bin/env python3
"""Generate the B5.16 finite-strain Norton creep C3D8T constant-force probe."""

from pathlib import Path

import generate_b511


if __name__ == "__main__":
    output = Path(__file__).with_name("b516_hex8_c3d8t_finite_norton.inp")
    contents = generate_b511.deck()
    contents = contents.replace(
        "B5.11 Abaqus/Standard small-strain Norton creep",
        "B5.16 Abaqus/Standard finite-strain Norton creep",
    ).replace("nlgeom=NO", "nlgeom=YES").replace("COORD, E, EE", "COORD, EE, LE").replace(
        "*Boundary, op=NEW",
        "*Controls, parameters=FIELD, field=DISPLACEMENT\n"
        "1.0000000000000000e-12, 1.0000000000000000e-12, , , , 1.0000000000000000e-12\n"
        "*Boundary, op=NEW",
    )
    output.write_text(contents, encoding="utf-8")
    print("wrote %s" % output)
