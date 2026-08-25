#!/usr/bin/env python3
"""Generate the B4.9 C3D8T thermo-mechanical element-operator probe."""

from pathlib import Path


COORDINATES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)

# The non-affine displacement and nonuniform temperature fields make all eight
# integration-point strains, stresses, temperatures, and heat fluxes observable.
TEMPERATURE = (360.0, 410.0, 445.0, 385.0, 470.0, 430.0, 515.0, 455.0)
DISPLACEMENT_X = (0.0, 1.2e-4, 1.7e-4, -0.4e-4, 0.3e-4, 1.0e-4, 2.3e-4, -0.8e-4)
DISPLACEMENT_Y = (0.0, -0.2e-4, 0.9e-4, 1.1e-4, -0.5e-4, 0.4e-4, 1.6e-4, 0.7e-4)
DISPLACEMENT_Z = (0.0, 0.3e-4, -0.4e-4, 0.2e-4, 1.4e-4, 1.0e-4, 1.8e-4, 1.1e-4)
BASE_STATE = TEMPERATURE + DISPLACEMENT_X + DISPLACEMENT_Y + DISPLACEMENT_Z
PERTURBATION = (1.0e-3,) * 8 + (1.0e-7,) * 24


def boundary_lines(state):
    lines = []
    for local_node in range(8):
        lines.append("%d, 11, 11, %.16e" % (local_node + 1, state[local_node]))
        for component in range(3):
            offset = 8 * (component + 1)
            lines.append(
                "%d, %d, %d, %.16e"
                % (local_node + 1, component + 1, component + 1, state[offset + local_node])
            )
    return lines


def step_lines(name, state, include_output=False):
    lines = [
        "*Step, name=%s, nlgeom=NO, inc=20" % name,
        "*Coupled Temperature-Displacement, steady state",
        "1.0, 1.0, 1.0e-8, 1.0",
        "*Boundary, op=NEW",
        *boundary_lines(state),
    ]
    if include_output:
        lines.extend(
            [
                "*Output, field, frequency=1",
                "*Node Output",
                "COORD, NT, RF, RFL, U",
                "*Element Output, directions=YES",
                "COORD, E, HFL, S, TEMP",
            ]
        )
    lines.append("*End Step")
    return lines


def deck():
    lines = [
        "*Heading",
        "** B4.9 Abaqus/Standard C3D8T full-integration thermo-mechanical element operator.",
        "** SI units: metre, kelvin, newton, pascal, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for label, coordinate in enumerate(COORDINATES, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + coordinate))
    lines.extend(
        [
            "*Element, type=C3D8T, elset=VOLUME",
            "1, 1, 2, 3, 4, 5, 6, 7, 8",
            "*Nset, nset=ALL_NODES, generate",
            "1, 8, 1",
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
        ]
    )
    lines.extend(step_lines("BASE", BASE_STATE, include_output=True))
    for column, step in enumerate(PERTURBATION):
        for suffix, sign in (("PLUS", 1.0), ("MINUS", -1.0)):
            state = list(BASE_STATE)
            state[column] += sign * step
            lines.extend(step_lines("D%02d_%s" % (column, suffix), state))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    output = Path(__file__).with_name("b49_hex8_c3d8t_operator_probe.inp")
    output.write_text(deck(), encoding="utf-8")
    print("wrote %s" % output)
