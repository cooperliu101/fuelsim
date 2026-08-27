#!/usr/bin/env python3
"""Generate the B5.9 multi-material mesh-refinement C3D8T probe."""

import math
from pathlib import Path


LENGTH = 2.0
HALF_THICKNESS = 0.1
WIDTH = 0.25
STEP_TIME = 1.0e7


CASES = (
    ("coarse", 0, 0, 2, 1, 1, False),
    ("refined", 1000, 1000, 4, 2, 2, False),
    ("distorted", 2000, 2000, 4, 2, 2, True),
)


def wrapped_labels(labels):
    return [", ".join(str(value) for value in labels[index : index + 12]) for index in range(0, len(labels), 12)]


def coordinates(nx, subdivisions_per_layer, nz, distorted):
    result = []
    ny = 2 * subdivisions_per_layer
    for iz in range(nz + 1):
        z = WIDTH * iz / nz
        for iy in range(ny + 1):
            y = -HALF_THICKNESS + 2.0 * HALF_THICKNESS * iy / ny
            for ix in range(nx + 1):
                x = LENGTH * ix / nx
                point = [x, y, z]
                if distorted:
                    sx = math.sin(math.pi * x / LENGTH)
                    sy = math.sin(math.pi * (y + HALF_THICKNESS) / (2.0 * HALF_THICKNESS))
                    sz = math.sin(math.pi * z / WIDTH)
                    point[0] += 0.040 * sx * sy * (2.0 * z / WIDTH - 1.0)
                    point[1] += 0.006 * sx * sy * sz
                    point[2] += 0.010 * sx * sy * sz * (2.0 * y / (2.0 * HALF_THICKNESS))
                result.append(tuple(point))
    return result


def connectivity(nx, subdivisions_per_layer, nz):
    ny = 2 * subdivisions_per_layer

    def node(ix, iy, iz):
        return iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix

    result = []
    for lower in (True, False):
        first_y = 0 if lower else subdivisions_per_layer
        last_y = subdivisions_per_layer if lower else ny
        for iz in range(nz):
            for iy in range(first_y, last_y):
                for ix in range(nx):
                    result.append(
                        (
                            node(ix, iy, iz),
                            node(ix + 1, iy, iz),
                            node(ix + 1, iy + 1, iz),
                            node(ix, iy + 1, iz),
                            node(ix, iy, iz + 1),
                            node(ix + 1, iy, iz + 1),
                            node(ix + 1, iy + 1, iz + 1),
                            node(ix, iy + 1, iz + 1),
                        )
                    )
    return result


def case_data(case):
    name, node_base, element_base, nx, subdivisions_per_layer, nz, distorted = case
    points = coordinates(nx, subdivisions_per_layer, nz, distorted)
    elements = connectivity(nx, subdivisions_per_layer, nz)
    ny = 2 * subdivisions_per_layer

    def node(ix, iy, iz):
        return node_base + 1 + iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix

    bottom = [node(ix, 0, iz) for iz in range(nz + 1) for ix in range(nx + 1)]
    top = [node(ix, ny, iz) for iz in range(nz + 1) for ix in range(nx + 1)]
    left_lower = [
        node(0, iy, iz) for iz in range(nz + 1) for iy in range(subdivisions_per_layer + 1)
    ]
    lower_count = nx * subdivisions_per_layer * nz
    return {
        "name": name,
        "node_base": node_base,
        "element_base": element_base,
        "points": points,
        "elements": elements,
        "bottom": bottom,
        "top": top,
        "left_lower": left_lower,
        "lower_count": lower_count,
    }


def step_lines(index, bottom_temperature):
    return [
        "*Step, name=STEP%d, nlgeom=NO, inc=4" % index,
        "*Coupled Temperature-Displacement",
        "%.16e, %.16e, %.16e, %.16e" % ((STEP_TIME,) * 4),
        "*Boundary, op=NEW",
        "BOTTOM, 11, 11, %.16e" % bottom_temperature,
        "TOP, 11, 11, 3.0000000000000000e2",
        "LEFT_LOWER, 1, 3, 0.0",
        "*Output, field, frequency=1",
        "*Node Output",
        "NT, RFL, RF, U",
        "*Element Output, directions=YES",
        "COORD, E, HFL, IVOL, S, TEMP",
        "*End Step",
    ]


def deck():
    data = [case_data(case) for case in CASES]
    lines = [
        "*Heading",
        "** B5.9 Abaqus/Standard multi-material C3D8T mesh-refinement and distortion probe.",
        "** Two bonded layers share nodes; three meshes undergo heating, holding and cooling.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for case in data:
        for local, point in enumerate(case["points"], 1):
            lines.append("%d, %.16e, %.16e, %.16e" % ((case["node_base"] + local,) + point))

    lines.append("*Element, type=C3D8T, elset=LOWER")
    for case in data:
        for local, element in enumerate(case["elements"][: case["lower_count"]], 1):
            labels = tuple(case["node_base"] + 1 + node for node in element)
            lines.append("%d, %s" % (case["element_base"] + local, ", ".join(str(node) for node in labels)))
    lines.append("*Element, type=C3D8T, elset=UPPER")
    for case in data:
        for offset, element in enumerate(case["elements"][case["lower_count"] :], case["lower_count"] + 1):
            labels = tuple(case["node_base"] + 1 + node for node in element)
            lines.append("%d, %s" % (case["element_base"] + offset, ", ".join(str(node) for node in labels)))

    for set_name, key in (("BOTTOM", "bottom"), ("TOP", "top"), ("LEFT_LOWER", "left_lower")):
        lines.append("*Nset, nset=%s" % set_name)
        labels = []
        for case in data:
            labels.extend(case[key])
        lines.extend(wrapped_labels(labels))
    all_nodes = []
    for case in data:
        all_nodes.extend(case["node_base"] + local for local in range(1, len(case["points"]) + 1))
    lines.append("*Nset, nset=ALL_NODES")
    lines.extend(wrapped_labels(all_nodes))

    lines.extend(
        [
            "*Material, name=LOWER_MATERIAL",
            "*Elastic",
            "1.2000000000000000e11, 2.8000000000000003e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "1.8000000000000000e-5",
            "*Conductivity",
            "8.0000000000000000e0",
            "*Density",
            "1.0000000000000000e4",
            "*Specific Heat",
            "3.0000000000000000e2",
            "*Material, name=UPPER_MATERIAL",
            "*Elastic",
            "2.0000000000000000e11, 3.0000000000000000e-1",
            "*Expansion, zero=3.0000000000000000e2",
            "7.0000000000000001e-6",
            "*Conductivity",
            "2.4000000000000000e1",
            "*Density",
            "6.5000000000000000e3",
            "*Specific Heat",
            "5.0000000000000000e2",
            "*Solid Section, elset=LOWER, material=LOWER_MATERIAL",
            ",",
            "*Solid Section, elset=UPPER, material=UPPER_MATERIAL",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.0000000000000000e2",
        ]
    )
    for index, temperature in enumerate((480.0, 600.0, 600.0, 330.0), 1):
        lines.extend(step_lines(index, temperature))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b59_hex8_c3d8t_multimaterial.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
