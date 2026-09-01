#!/usr/bin/env python3
"""Generate the B5.48 M5.8 mesh-upgraded Abaqus C3D20T case."""

from pathlib import Path

from generate_b546 import generate


def main():
    generate(
        "b548_m58_c3d20t_integrated",
        "B5.48",
        "C3D20T",
        "m58_integrated_hex20_mesh.e",
        5969,
        1152,
    )
    path = Path(__file__).resolve().parent / "b548_m58_c3d20t_integrated.inp"
    with path.open("r") as source:
        contents = source.read()
    requested = "CEEQ, IVOL, PEEQ, S"
    if contents.count(requested) != 1:
        raise RuntimeError("unexpected B5.48 element output request")
    with path.open("w") as output:
        output.write(contents.replace(requested, "CEEQ, COORD, IVOL, PEEQ, S"))


if __name__ == "__main__":
    main()
