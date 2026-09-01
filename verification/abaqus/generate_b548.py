#!/usr/bin/env python3
"""Generate the B5.48 M5.8 mesh-upgraded Abaqus C3D20T case."""

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


if __name__ == "__main__":
    main()
