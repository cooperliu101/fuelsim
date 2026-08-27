#!/usr/bin/env python3
"""Generate the B5.7 temperature-dependent C3D8T lumped-capacity probe."""

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
OLD_TEMPERATURES = (310.0, 330.0, 350.0, 370.0, 390.0, 410.0, 430.0, 450.0)
BASE_TEMPERATURES = (360.0, 410.0, 445.0, 385.0, 470.0, 430.0, 515.0, 455.0)
PERTURBATION = 1.0e-3


def state_name(state):
    if state == 0:
        return "BASE"
    node = (state - 1) // 2
    sign = "PLUS" if (state - 1) % 2 == 0 else "MINUS"
    return "D%02d_%s" % (node, sign)


def state_temperatures(state):
    result = list(BASE_TEMPERATURES)
    if state != 0:
        node = (state - 1) // 2
        result[node] += PERTURBATION if (state - 1) % 2 == 0 else -PERTURBATION
    return tuple(result)


def element_index(state, transient):
    return 2 * state + (1 if transient else 0)


def node_label(element, local_node):
    return 8 * element + local_node + 1


def boundary_lines(transient_step):
    lines = []
    for state in range(17):
        final = state_temperatures(state)
        for transient_copy in (False, True):
            element = element_index(state, transient_copy)
            temperatures = final if (transient_step or not transient_copy) else OLD_TEMPERATURES
            for local_node, temperature in enumerate(temperatures):
                label = node_label(element, local_node)
                lines.append("%d, 1, 3, 0.0" % label)
                lines.append("%d, 11, 11, %.16e" % (label, temperature))
    return lines


def deck():
    element_count = 34
    lines = [
        "*Heading",
        "** B5.7 Abaqus/Standard temperature-dependent C3D8T lumped-capacity probe.",
        "** Each state has a steady copy and a transient copy on the same distorted geometry.",
        "** Their final reaction-flux difference isolates the transient capacity residual.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for element in range(element_count):
        shift_x = 3.0 * (element % 9)
        shift_y = 3.0 * (element // 9)
        for local_node, coordinate in enumerate(COORDINATES):
            shifted = (coordinate[0] + shift_x, coordinate[1] + shift_y, coordinate[2])
            lines.append("%d, %.16e, %.16e, %.16e" % ((node_label(element, local_node),) + shifted))
    lines.append("*Element, type=C3D8T, elset=VOLUME")
    for element in range(element_count):
        connectivity = [node_label(element, local_node) for local_node in range(8)]
        lines.append("%d, %s" % (element + 1, ", ".join(str(label) for label in connectivity)))
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, %d, 1" % (8 * element_count),
            "*Material, name=THERMAL",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3, 3.0000000000000000e2",
            "1.7000000000000000e3, 6.0000000000000000e2",
            "*Specific Heat",
            "3.0000000000000000e3, 3.0000000000000000e2",
            "4.2000000000000000e3, 6.0000000000000000e2",
            "*Solid Section, elset=VOLUME, material=THERMAL",
            ",",
            "*Initial Conditions, type=TEMPERATURE",
        ]
    )
    for element in range(element_count):
        for local_node, temperature in enumerate(OLD_TEMPERATURES):
            lines.append("%d, %.16e" % (node_label(element, local_node), temperature))
    lines.extend(
        [
            "*Step, name=STEADY, nlgeom=NO, inc=20",
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0e-8, 1.0",
            "*Boundary, op=NEW",
            *boundary_lines(False),
            "*Output, field, frequency=1",
            "*Node Output",
            "NT, RFL",
            "*End Step",
            "*Step, name=TRANSIENT, nlgeom=NO, inc=4",
            "*Coupled Temperature-Displacement",
            "1.0, 1.0, 1.0, 1.0",
            "*Boundary, op=NEW",
            *boundary_lines(True),
            "*Output, field, frequency=1",
            "*Node Output",
            "NT, RFL",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b57_hex8_c3d8t_temperature_capacity.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
