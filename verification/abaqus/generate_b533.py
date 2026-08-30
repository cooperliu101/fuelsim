#!/usr/bin/env python3
"""Generate a paired default/near-zero hourglass finite C3D8RT probe."""

from pathlib import Path

from generate_b532 import DISPLACEMENT, GEOMETRIES, TEMPERATURE


RAW_MODE = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]


def pure_hourglass(amplitude=(0.012, -0.007, 0.009)):
    result = []
    for node in range(8):
        result.append(
            (
                amplitude[0] * RAW_MODE[node],
                amplitude[1] * RAW_MODE[node],
                amplitude[2] * RAW_MODE[node],
            )
        )
    return result


def probe_cases():
    cases = []
    states = (
        ("general", DISPLACEMENT),
        ("hourglass_vector", pure_hourglass()),
        ("hourglass_x_small", pure_hourglass((1.0e-4, 0.0, 0.0))),
        ("hourglass_x_large", pure_hourglass((0.012, 0.0, 0.0))),
        ("hourglass_y_large", pure_hourglass((0.0, -0.012, 0.0))),
        ("hourglass_z_large", pure_hourglass((0.0, 0.0, 0.012))),
    )
    for geometry_name in ("regular", "warped"):
        for state_name, displacement in states:
            cases.append((geometry_name + "_" + state_name, GEOMETRIES[geometry_name], displacement))
    return cases


cases = probe_cases()

lines = [
    "*Heading",
    "** B5.33 Abaqus/Standard finite-strain C3D8RT isolated default hourglass probe.",
    "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
    "*Node",
]
labels = []
for case_index, (name, coordinates, displacement) in enumerate(cases):
    for control in ("default", "weak"):
        copy_index = len(labels)
        labels.append((name, control, coordinates, displacement))
        shift = 3.0 * copy_index
        for local, coordinate in enumerate(coordinates):
            node = 8 * copy_index + local + 1
            lines.append("%d, %.16e, %.16e, %.16e" % (node, coordinate[0] + shift, coordinate[1], coordinate[2]))
lines.append("*Element, type=C3D8RT")
for copy_index in range(len(labels)):
    nodes = [8 * copy_index + local + 1 for local in range(8)]
    lines.append("%d, %s" % (copy_index + 1, ", ".join(str(node) for node in nodes)))
default_elements = [str(index + 1) for index, label in enumerate(labels) if label[1] == "default"]
weak_elements = [str(index + 1) for index, label in enumerate(labels) if label[1] == "weak"]
lines.extend(
    [
        "*Elset, elset=DEFAULT_HOURGLASS",
        ", ".join(default_elements),
        "*Elset, elset=WEAK_HOURGLASS",
        ", ".join(weak_elements),
        "*Nset, nset=ALL_NODES, generate",
        "1, %d, 1" % (8 * len(labels)),
        "*Material, name=ELASTIC",
        "*Elastic",
        "2.0000000000000000e11, 2.5000000000000000e-1",
        "*Expansion, zero=3.0000000000000000e2",
        "1.2000000000000000e-5",
        "*Conductivity",
        "4.0",
        "*Density",
        "2.0e3",
        "*Specific Heat",
        "3.0e3",
        "*Solid Section, elset=DEFAULT_HOURGLASS, material=ELASTIC",
        ",",
        "*Solid Section, elset=WEAK_HOURGLASS, material=ELASTIC",
        ",",
        "*Hourglass Stiffness",
        "1.0",
        "*Initial Conditions, type=TEMPERATURE",
        "ALL_NODES, 3.0e2",
        "*Step, name=FINITE_HOURGLASS, nlgeom=YES, inc=1",
        "*Coupled Temperature-Displacement",
        "1.0, 1.0, 1.0e-8, 1.0",
        "*Boundary",
    ]
)
for copy_index, (_, _, _, displacement) in enumerate(labels):
    for local in range(8):
        node = 8 * copy_index + local + 1
        for component in range(3):
            lines.append("%d, %d, %d, %.16e" % (node, component + 1, component + 1, displacement[local][component]))
        lines.append("%d, 11, 11, %.16e" % (node, TEMPERATURE[local]))
lines.extend(
    [
        "*Output, field, frequency=1",
        "*Node Output",
        "RF, U",
        "*Element Output",
        "IVOL, S, TEMP",
        "*End Step",
        "",
    ]
)
output = Path(__file__).with_name("b533_hex8_c3d8rt_finite_hourglass_probe.inp")
if __name__ == "__main__":
    output.write_text("\n".join(lines))
    print(output)
