#!/usr/bin/env python3
"""Generate finite-strain C3D8RT plasticity, creep, and coupled material probes."""

from pathlib import Path

import generate_b510
import generate_b511
import generate_b512


def finite_reduced(contents, old_title, new_title):
    return (
        contents.replace(old_title, new_title)
        .replace("type=C3D8T", "type=C3D8RT")
        .replace("nlgeom=NO", "nlgeom=YES")
        .replace("COORD, E, EE", "COORD, EE, LE")
        .replace("C3D8RT C3D8T", "C3D8RT")
        .replace("C3D8RT J2 plasticity-Norton creep C3D8T", "C3D8RT J2 plasticity-Norton creep")
    )


def controlled(contents):
    return contents.replace(
        "*Boundary, op=NEW",
        "*Controls, parameters=FIELD, field=DISPLACEMENT\n"
        "1.0000000000000000e-12, 1.0000000000000000e-12, , , , 1.0000000000000000e-12\n"
        "*Boundary, op=NEW",
    )


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    cases = {
        "b535_hex8_c3d8rt_finite_j2.inp": finite_reduced(
            generate_b510.deck(),
            "B5.10 Abaqus/Standard small-strain J2 plasticity",
            "B5.35 Abaqus/Standard finite-strain J2 plasticity C3D8RT",
        ),
        "b536_hex8_c3d8rt_finite_norton.inp": controlled(
            finite_reduced(
                generate_b511.deck(),
                "B5.11 Abaqus/Standard small-strain Norton creep",
                "B5.36 Abaqus/Standard finite-strain Norton creep C3D8RT",
            )
        ),
        "b537_hex8_c3d8rt_finite_coupled.inp": controlled(
            finite_reduced(
                generate_b512.deck(),
                "B5.12 Abaqus/Standard small-strain active coupled",
                "B5.37 Abaqus/Standard finite-strain active coupled C3D8RT",
            )
        ),
    }
    for name, contents in cases.items():
        output = directory / name
        output.write_text(contents, encoding="utf-8")
        print("wrote %s" % output)
