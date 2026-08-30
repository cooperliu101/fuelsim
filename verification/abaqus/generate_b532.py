#!/usr/bin/env python3
"""Generate independent finite-strain C3D8RT operator probes."""

from pathlib import Path


GEOMETRIES = {
    "regular": [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ],
    "warped": [
        (0.00, 0.00, 0.00),
        (1.20, 0.10, -0.05),
        (1.10, 1.00, 0.10),
        (-0.10, 0.90, 0.00),
        (0.05, -0.05, 1.00),
        (1.15, 0.00, 1.20),
        (1.00, 1.10, 1.10),
        (-0.05, 1.00, 0.90),
    ],
}

TEMPERATURE = [360.0, 410.0, 445.0, 385.0, 470.0, 430.0, 515.0, 455.0]
DISPLACEMENT = [
    (0.000, 0.000, 0.000),
    (0.048, -0.008, 0.012),
    (0.068, 0.036, -0.016),
    (-0.016, 0.044, 0.008),
    (0.012, -0.020, 0.056),
    (0.040, 0.016, 0.040),
    (0.092, 0.064, 0.072),
    (-0.032, 0.028, 0.044),
]


def state(case):
    temperature = list(TEMPERATURE)
    displacement = [list(value) for value in DISPLACEMENT]
    if case:
        column, sign = case
        if column < 8:
            temperature[column] += sign * 1.0e-3
        else:
            local = column - 8
            component = local // 8
            node = local % 8
            displacement[node][component] += sign * 1.0e-7
    return temperature, displacement


def generate(name, coordinates):
    cases = [("BASE", None)]
    for column in range(32):
        cases.append(("D%02d_PLUS" % column, (column, 1.0)))
        cases.append(("D%02d_MINUS" % column, (column, -1.0)))

    lines = [
        "*Heading",
        "** B5.32 Abaqus/Standard finite-strain C3D8RT complete operator probe: %s." % name,
        "** Each perturbation uses a disconnected element so every state follows one increment from the reference state.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for case_index, (_, _) in enumerate(cases):
        shift = 3.0 * case_index
        for local, coordinate in enumerate(coordinates):
            node = 8 * case_index + local + 1
            lines.append("%d, %.16e, %.16e, %.16e" % (node, coordinate[0] + shift, coordinate[1], coordinate[2]))
    lines.append("*Element, type=C3D8RT, elset=VOLUME")
    for case_index, (_, _) in enumerate(cases):
        nodes = [8 * case_index + local + 1 for local in range(8)]
        lines.append("%d, %s" % (case_index + 1, ", ".join(str(node) for node in nodes)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, %d, 1" % (8 * len(cases)),
            "*Material, name=THERMOELASTIC",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "1.2000000000000000e-5",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3",
            "*Specific Heat",
            "3.0000000000000000e3",
            "*Solid Section, elset=VOLUME, material=THERMOELASTIC",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
            "*Step, name=FINITE_OPERATOR, nlgeom=YES, inc=1",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0e-8, 1.0",
            "*Boundary",
        ]
    )
    for case_index, (_, perturbation) in enumerate(cases):
        temperature, displacement = state(perturbation)
        for local in range(8):
            node = 8 * case_index + local + 1
            lines.append("%d, 11, 11, %.16e" % (node, temperature[local]))
            for component in range(3):
                lines.append("%d, %d, %d, %.16e" % (node, component + 1, component + 1, displacement[local][component]))
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Element Output, directions=YES",
            "COORD, HFL, IVOL, LE, S, TEMP",
            "*End Step",
            "",
        ]
    )
    output = Path(__file__).with_name("b532_hex8_c3d8rt_finite_%s_operator_probe.inp" % name)
    output.write_text("\n".join(lines))
    print(output)


if __name__ == "__main__":
    for geometry_name, geometry_coordinates in GEOMETRIES.items():
        generate(geometry_name, geometry_coordinates)
