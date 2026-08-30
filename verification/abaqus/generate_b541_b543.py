#!/usr/bin/env python3
"""Generate small-strain C3D8RT plasticity, creep, and coupled material paths."""

from pathlib import Path

import generate_b510
import generate_b511
import generate_b512


def reduced(contents, old_title, new_title):
    return contents.replace(old_title, new_title).replace("type=C3D8T", "type=C3D8RT")


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    cases = {
        "b541_hex8_c3d8rt_small_j2.inp": reduced(
            generate_b510.deck(),
            "B5.10 Abaqus/Standard small-strain J2 plasticity",
            "B5.41 Abaqus/Standard small-strain J2 plasticity C3D8RT",
        ),
        "b542_hex8_c3d8rt_small_norton.inp": reduced(
            generate_b511.deck(),
            "B5.11 Abaqus/Standard small-strain Norton creep",
            "B5.42 Abaqus/Standard small-strain Norton creep C3D8RT",
        ),
        "b543_hex8_c3d8rt_small_coupled.inp": reduced(
            generate_b512.deck(),
            "B5.12 Abaqus/Standard small-strain active coupled",
            "B5.43 Abaqus/Standard small-strain active coupled C3D8RT",
        ),
    }
    for name, contents in cases.items():
        output = directory / name
        output.write_text(contents, encoding="utf-8")
        print("wrote %s" % output)
