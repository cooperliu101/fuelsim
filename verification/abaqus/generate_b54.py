#!/usr/bin/env python3
"""Generate the B5.4 finite-deformation C3D8T thermal-load probe."""

from pathlib import Path


COORDINATES = (
    (0.00, 0.00, 0.00),
    (1.20, 0.10, -0.05),
    (1.10, 1.00, 0.10),
    (-0.10, 0.90, 0.00),
    (0.05, -0.10, 1.00),
    (1.30, 0.00, 1.10),
    (1.00, 1.20, 0.90),
    (-0.20, 1.00, 1.20),
)
CASE_COUNT = 11
FILM_TEMPERATURE = (300.0, 300.0, 300.0, 300.0, 360.0, 410.0, 445.0, 385.0)


def node_label(case, local_node):
    return 8 * case + local_node + 1


def displacement(coordinate):
    x, y, z = coordinate
    return (0.5 * x + 0.1 * y, 0.25 * y + 0.05 * z, 0.08 * x - 0.2 * z)


def boundary_lines(temperatures):
    lines = ["*Boundary, op=NEW"]
    for case in range(CASE_COUNT):
        for local_node, coordinate in enumerate(COORDINATES):
            label = node_label(case, local_node)
            for component, value in enumerate(displacement(coordinate), 1):
                lines.append("%d, %d, %d, %.16e" % (label, component, component, value))
            lines.append("%d, 11, 11, %.16e" % (label, temperatures[case][local_node]))
    return lines


def steady_step(name, temperatures, load_lines):
    return [
        "*Step, name=%s, nlgeom=YES, inc=20" % name,
        "*Coupled Temperature-Displacement, steady state",
        "1.0, 1.0, 1.0e-8, 1.0",
        *boundary_lines(temperatures),
        *load_lines,
        "*Output, field, frequency=1",
        "*Node Output",
        "COORD, NT, RFL, U",
        "*Element Output",
        "IVOL, TEMP",
        "*End Step",
    ]


def deck():
    lines = [
        "*Heading",
        "** B5.4 Abaqus/Standard C3D8T finite-deformation thermal-load probe.",
        "** SI units: metre, second, kelvin, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for case in range(CASE_COUNT):
        for local_node, coordinate in enumerate(COORDINATES):
            shifted = (coordinate[0] + 3.0 * case, coordinate[1], coordinate[2])
            lines.append("%d, %.16e, %.16e, %.16e" % ((node_label(case, local_node),) + shifted))
    lines.append("*Element, type=C3D8T, elset=VOLUME")
    for case in range(CASE_COUNT):
        connectivity = [node_label(case, local_node) for local_node in range(8)]
        lines.append("%d, %s" % (case + 1, ", ".join(str(label) for label in connectivity)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, %d, 1" % (8 * CASE_COUNT),
            "*Material, name=THERMAL",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "0.0",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3",
            "*Specific Heat",
            "3.0000000000000000e3",
            "*Solid Section, elset=VOLUME, material=THERMAL",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
        ]
    )
    uniform = (300.0,) * 8
    uniform_cases = (uniform,) * CASE_COUNT
    lines.extend(steady_step("DEFORM", uniform_cases, []))
    capacity_temperatures = []
    for case in range(CASE_COUNT):
        if case < 8:
            capacity_temperatures.append(tuple(301.0 if node == case else 300.0 for node in range(8)))
        else:
            capacity_temperatures.append(uniform)
    lines.extend(
        [
            "*Step, name=CAPACITY, nlgeom=YES, inc=4",
            "*Coupled Temperature-Displacement",
            "1.0, 1.0, 1.0, 1.0",
            *boundary_lines(tuple(capacity_temperatures)),
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RFL, U",
            "*Element Output",
            "IVOL, TEMP",
            "*End Step",
        ]
    )
    lines.extend(steady_step("BODY", uniform_cases, ["*Dflux, op=NEW", "9, BF, 8.0000000000000000e1"]))
    lines.extend(steady_step("SURFACE", uniform_cases, ["*Dflux, op=NEW", "10, S2, 4.0000000000000000e1"]))
    film_cases = list(uniform_cases)
    film_cases[10] = FILM_TEMPERATURE
    lines.extend(steady_step("FILM_BASE", tuple(film_cases), ["*Dflux, op=NEW", "11, BF, 0.0"]))
    lines.extend(
        steady_step(
            "FILM_ACTIVE",
            tuple(film_cases),
            [
                "*Dflux, op=NEW",
                "11, BF, 0.0",
                "*Film, op=NEW",
                "11, F2, 2.5000000000000000e2, 1.0000000000000000e1",
            ],
        )
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b54_hex8_c3d8t_finite_thermal_load_probe.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
