#!/usr/bin/env python3
"""Generate B5.20 C3D8T gap-conductance derivative probes."""

from pathlib import Path


UNIT_COORDINATES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)


def probe_steps(base_displacement):
    return (
        ("BASE", base_displacement, 400.0, 300.0),
        ("GAP_PLUS", base_displacement - 1.0e-5, 400.0, 300.0),
        ("GAP_MINUS", base_displacement + 1.0e-5, 400.0, 300.0),
        ("SECONDARY_PLUS", base_displacement, 401.0, 300.0),
        ("SECONDARY_MINUS", base_displacement, 399.0, 300.0),
        ("PRIMARY_PLUS", base_displacement, 400.0, 301.0),
        ("PRIMARY_MINUS", base_displacement, 400.0, 299.0),
    )


def append_step(lines, name, displacement, secondary_temperature, primary_temperature):
    lines.extend(
        [
            "*Step, name=%s, nlgeom=YES, inc=40" % name,
            "*Coupled Temperature-Displacement, steady state",
            "1.0, 1.0, 1.0e-8, 1.0",
            "*Boundary, op=NEW",
            "SECONDARY, 1, 1, %.16e" % displacement,
            "SECONDARY, 2, 3, 0.0",
            "PRIMARY, 1, 3, 0.0",
            "SECONDARY, 11, 11, %.16e" % secondary_temperature,
            "PRIMARY, 11, 11, %.16e" % primary_temperature,
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*End Step",
        ]
    )


def deck(pressure_dependent):
    dependency = "pressure" if pressure_dependent else "clearance"
    lines = [
        "*Heading",
        "** B5.20 Abaqus/Standard C3D8T %s- and temperature-dependent gap conductance." % dependency,
        "** SI units: metre, second, kelvin, watt, pascal.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for block, shift in enumerate((0.0, 1.1)):
        for local_node, coordinate in enumerate(UNIT_COORDINATES):
            node = 8 * block + local_node + 1
            lines.append(
                "%d, %.16e, %.16e, %.16e" % (node, coordinate[0] + shift, coordinate[1], coordinate[2])
            )
    lines.extend(
        [
            "*Element, type=C3D8T, elset=BLOCKS",
            "1, 1, 2, 3, 4, 5, 6, 7, 8",
            "2, 9, 10, 11, 12, 13, 14, 15, 16",
            "*Nset, nset=ALL_NODES, generate",
            "1, 16, 1",
            "*Nset, nset=SECONDARY, generate",
            "1, 8, 1",
            "*Nset, nset=PRIMARY, generate",
            "9, 16, 1",
            "*Surface, type=ELEMENT, name=SECONDARY_CONTACT",
            "1, S4",
            "*Surface, type=ELEMENT, name=PRIMARY_CONTACT",
            "2, S6",
            "*Material, name=THERMAL",
            "*Elastic",
            "2.0000000000000000e11, 2.5000000000000000e-1",
            "*Conductivity",
            "4.0000000000000000e0",
            "*Density",
            "2.0000000000000000e3",
            "*Specific Heat",
            "3.0000000000000000e3",
            "*Solid Section, elset=BLOCKS, material=THERMAL",
            ",",
            "*Surface Interaction, name=THERMAL_GAP",
        ]
    )
    if pressure_dependent:
        lines.extend(
            [
                "*Surface Behavior, pressure-overclosure=LINEAR",
                "1.0000000000000000e5",
                "*Gap Conductance, pressure",
                "5.9000000000000000e1, 5.0000000000000000e2, 3.4000000000000000e2",
                "7.9000000000000000e1, 1.5000000000000000e3, 3.4000000000000000e2",
                "6.1000000000000000e1, 5.0000000000000000e2, 3.6000000000000000e2",
                "8.1000000000000000e1, 1.5000000000000000e3, 3.6000000000000000e2",
            ]
        )
    else:
        lines.extend(
            [
                "*Gap Conductance",
                "9.9000000000000000e1, 0.0000000000000000e0, 3.4000000000000000e2",
                "3.9000000000000000e1, 5.9999999999999998e-2, 3.4000000000000000e2",
                "1.0100000000000000e2, 0.0000000000000000e0, 3.6000000000000000e2",
                "4.1000000000000000e1, 5.9999999999999998e-2, 3.6000000000000000e2",
            ]
        )
    lines.extend(
        [
            "*Contact Pair, interaction=THERMAL_GAP, type=SURFACE TO SURFACE, small sliding, adjust=0.",
            "SECONDARY_CONTACT, PRIMARY_CONTACT",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 3.5000000000000000e2",
        ]
    )
    for step in probe_steps(0.11 if pressure_dependent else 0.05):
        append_step(lines, *step)
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    directory = Path(__file__).parent
    for suffix, pressure_dependent in (("clearance", False), ("pressure", True)):
        output = directory / ("b520_hex8_c3d8t_gap_conductance_%s.inp" % suffix)
        output.write_text(deck(pressure_dependent), encoding="utf-8")
        print("wrote %s" % output)
