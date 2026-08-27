#!/usr/bin/env python3
"""Generate the B5.3 C3D8T reference-configuration thermal-load probe."""

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
PLANAR_TRAPEZOID = (
    (0.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (1.5, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (2.0, 0.0, 1.0),
    (1.5, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)
WARPED_VERTICAL = tuple(
    (x, y, z - 1.0) for x, y, z in COORDINATES[4:]
) + COORDINATES[4:]


def single_node_warp(height):
    top = ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0 + height), (0.0, 1.0, 1.0))
    return tuple((x, y, z - 1.0) for x, y, z in top) + top


GEOMETRIES = (
    COORDINATES,
    COORDINATES,
    COORDINATES,
    PLANAR_TRAPEZOID,
    WARPED_VERTICAL,
    single_node_warp(0.1),
    single_node_warp(0.5),
    single_node_warp(1.0),
)
FILM_TEMPERATURE = (300.0, 300.0, 300.0, 300.0, 360.0, 410.0, 445.0, 385.0)


def node_label(case, local_node):
    return 8 * case + local_node + 1


def boundary_lines(temperature):
    lines = []
    for case in range(len(GEOMETRIES)):
        for local_node in range(8):
            label = node_label(case, local_node)
            lines.append("%d, 1, 3, 0.0" % label)
            lines.append("%d, 11, 11, %.16e" % (label, temperature[local_node]))
    return lines


def step_lines(name, temperature, load_lines):
    return [
        "*Step, name=%s, nlgeom=NO, inc=20" % name,
        "*Coupled Temperature-Displacement, steady state",
        "1.0, 1.0, 1.0e-8, 1.0",
        "*Boundary, op=NEW",
        *boundary_lines(temperature),
        *load_lines,
        "*Output, field, frequency=1",
        "*Node Output",
        "NT, RFL",
        "*Element Output",
        "IVOL, TEMP",
        "*End Step",
    ]


def deck():
    lines = [
        "*Heading",
        "** B5.3 Abaqus/Standard C3D8T corner-integrated thermal-load probe.",
        "** SI units: metre, kelvin, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for case, geometry in enumerate(GEOMETRIES):
        for local_node, coordinate in enumerate(geometry):
            shifted = (coordinate[0] + 3.0 * case, coordinate[1], coordinate[2])
            lines.append("%d, %.16e, %.16e, %.16e" % ((node_label(case, local_node),) + shifted))
    lines.append("*Element, type=C3D8T, elset=VOLUME")
    for case in range(len(GEOMETRIES)):
        connectivity = [node_label(case, local_node) for local_node in range(8)]
        lines.append("%d, %s" % (case + 1, ", ".join(str(label) for label in connectivity)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, %d, 1" % (8 * len(GEOMETRIES)),
            "*Material, name=THERMAL",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
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
    lines.extend(step_lines("BODY", uniform, ["*Dflux, op=NEW", "1, BF, 8.0000000000000000e1"]))
    lines.extend(step_lines("SURFACE", uniform, ["*Dflux, op=NEW", "2, S2, 4.0000000000000000e1"]))
    lines.extend(step_lines("FILM_BASE", FILM_TEMPERATURE, ["*Dflux, op=NEW", "3, BF, 0.0"]))
    lines.extend(
        step_lines(
            "FILM_ACTIVE",
            FILM_TEMPERATURE,
            [
                "*Dflux, op=NEW",
                "3, BF, 0.0",
                "*Film, op=NEW",
                "3, F2, 2.5000000000000000e2, 1.0000000000000000e1",
            ],
        )
    )
    lines.extend(
        step_lines(
            "SURFACE_GEOMETRY",
            uniform,
            [
                "*Dflux, op=NEW",
                "2, S2, 4.0000000000000000e1",
                "4, S2, 4.0000000000000000e1",
                "5, S2, 4.0000000000000000e1",
                "6, S2, 4.0000000000000000e1",
                "7, S2, 4.0000000000000000e1",
                "8, S2, 4.0000000000000000e1",
            ],
        )
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b53_hex8_c3d8t_thermal_load_probe.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
