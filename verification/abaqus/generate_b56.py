#!/usr/bin/env python3
"""Generate the B5.6 temperature-dependent C3D8T operator probe."""

from pathlib import Path

import generate_b49


def deck():
    result = generate_b49.deck()
    result = result.replace(
        "** B4.9 Abaqus/Standard C3D8T full-integration thermo-mechanical element operator.",
        "** B5.6 Abaqus/Standard temperature-dependent C3D8T element operator.",
    )
    result = result.replace(
        "*Elastic\n2.0000000000000000e11, 2.5000000000000000e-1",
        "*Elastic\n2.0000000000000000e11, 2.5000000000000000e-1, 3.0000000000000000e2\n"
        "1.7000000000000000e11, 2.8000000000000003e-1, 6.0000000000000000e2",
    )
    result = result.replace(
        "*Expansion, zero=3.0000000000000000e2\n1.2000000000000000e-5",
        "*Expansion, zero=3.0000000000000000e2\n"
        "1.2000000000000000e-5, 3.0000000000000000e2\n"
        "1.8000000000000000e-5, 6.0000000000000000e2",
    )
    result = result.replace(
        "*Conductivity\n4.0000000000000000e0",
        "*Conductivity\n4.0000000000000000e0, 3.0000000000000000e2\n"
        "7.0000000000000000e0, 6.0000000000000000e2",
    )
    return result


if __name__ == "__main__":
    output = Path(__file__).with_name("b56_hex8_c3d8t_temperature_operator.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
