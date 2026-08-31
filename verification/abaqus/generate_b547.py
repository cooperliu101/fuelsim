#!/usr/bin/env python3
"""Generate the B5.47 full-size M5.8-equivalent Abaqus C3D8RT case."""

from generate_b546 import generate


def main():
    generate("b547_m58_c3d8rt_integrated", "B5.47", "C3D8RT")


if __name__ == "__main__":
    main()
