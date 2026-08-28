#!/usr/bin/env python3
"""Generate all B5.27 engineering fuel-cladding Abaqus full-field cases."""

import math
from pathlib import Path


CASES = {
    "b527_coarse": ((1, 1, 3, 2), 1000.0, 1.0e14),
    "b527_medium": ((2, 1, 4, 3), 1000.0, 1.0e14),
    "b527_fine": ((2, 2, 6, 4), 1000.0, 1.0e14),
    "b527_half_step": ((2, 1, 4, 3), 500.0, 1.0e14),
    "b527_low_penalty": ((2, 1, 4, 3), 1000.0, 5.0e13),
    "b527_high_penalty": ((2, 1, 4, 3), 1000.0, 2.0e14),
}


def append_set(lines, keyword, name, labels):
    lines.append("*%s, %s=%s" % (keyword, keyword.lower(), name))
    for begin in range(0, len(labels), 16):
        lines.append(", ".join(str(value) for value in labels[begin : begin + 16]))


def append_region(nodes, elements, inner_radius, outer_radius, radial_elements, angular_elements, axial_elements):
    node_map = {}
    for axial in range(axial_elements + 1):
        for angular in range(angular_elements + 1):
            for radial in range(radial_elements + 1):
                radius = inner_radius + (outer_radius - inner_radius) * float(radial) / radial_elements
                angle = 0.5 * math.pi * float(angular) / angular_elements
                label = len(nodes) + 1
                node_map[(radial, angular, axial)] = label
                nodes.append(
                    (
                        radius * math.cos(angle),
                        radius * math.sin(angle),
                        4.0e-2 * float(axial) / axial_elements,
                    )
                )
    for axial in range(axial_elements):
        for angular in range(angular_elements):
            for radial in range(radial_elements):
                elements.append(
                    (
                        node_map[(radial, angular, axial)],
                        node_map[(radial + 1, angular, axial)],
                        node_map[(radial + 1, angular + 1, axial)],
                        node_map[(radial, angular + 1, axial)],
                        node_map[(radial, angular, axial + 1)],
                        node_map[(radial + 1, angular, axial + 1)],
                        node_map[(radial + 1, angular + 1, axial + 1)],
                        node_map[(radial, angular + 1, axial + 1)],
                    )
                )
    return node_map


def deck(name, divisions, step, penalty):
    fuel_radial, clad_radial, angular_elements, axial_elements = divisions
    nodes = []
    fuel_elements = []
    clad_elements = []
    fuel_nodes = append_region(
        nodes, fuel_elements, 1.0e-3, 4.0e-3, fuel_radial, angular_elements, axial_elements
    )
    clad_nodes = append_region(
        nodes, clad_elements, 4.005e-3, 4.7e-3, clad_radial, angular_elements, axial_elements
    )
    fuel_labels = list(range(1, len(fuel_elements) + 1))
    clad_labels = list(range(len(fuel_elements) + 1, len(fuel_elements) + len(clad_elements) + 1))

    def element_label(labels, radial_count, radial, angular, axial):
        return labels[(axial * angular_elements + angular) * radial_count + radial]

    fuel_contact_elements = [
        element_label(fuel_labels, fuel_radial, fuel_radial - 1, angular, axial)
        for axial in range(axial_elements)
        for angular in range(angular_elements)
    ]
    clad_contact_elements = [
        element_label(clad_labels, clad_radial, 0, angular, axial)
        for axial in range(axial_elements)
        for angular in range(angular_elements)
    ]
    fuel_inner_nodes = [
        fuel_nodes[(0, angular, axial)]
        for axial in range(axial_elements + 1)
        for angular in range(angular_elements + 1)
    ]
    clad_outer_nodes = [
        clad_nodes[(clad_radial, angular, axial)]
        for axial in range(axial_elements + 1)
        for angular in range(angular_elements + 1)
    ]
    secondary_contact_nodes = [
        clad_nodes[(0, angular, axial)]
        for axial in range(axial_elements + 1)
        for angular in range(angular_elements + 1)
    ]

    def face_nodes(node_map, radial_count, face):
        if face == "theta0":
            return [node_map[(radial, 0, axial)] for axial in range(axial_elements + 1) for radial in range(radial_count + 1)]
        if face == "theta90":
            return [node_map[(radial, angular_elements, axial)] for axial in range(axial_elements + 1) for radial in range(radial_count + 1)]
        return [node_map[(radial, angular, 0)] for angular in range(angular_elements + 1) for radial in range(radial_count + 1)]

    lines = [
        "*Heading",
        "** %s engineering quarter-cylinder fuel-cladding C3D8T full-field reference." % name,
        "** fuel_radial=%d clad_radial=%d angular=%d axial=%d step=%.16e penalty=%.16e"
        % (fuel_radial, clad_radial, angular_elements, axial_elements, step, penalty),
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES",
        "*Node",
    ]
    for label, point in enumerate(nodes, 1):
        lines.append("%d, %.16e, %.16e, %.16e" % ((label,) + point))
    lines.append("*Element, type=C3D8T, elset=FUEL")
    for label, connectivity in zip(fuel_labels, fuel_elements):
        lines.append("%d, %s" % (label, ", ".join(str(value) for value in connectivity)))
    lines.append("*Element, type=C3D8T, elset=CLAD")
    for label, connectivity in zip(clad_labels, clad_elements):
        lines.append("%d, %s" % (label, ", ".join(str(value) for value in connectivity)))
    append_set(lines, "Elset", "FUEL_CONTACT_ELEMENTS", fuel_contact_elements)
    append_set(lines, "Elset", "CLAD_CONTACT_ELEMENTS", clad_contact_elements)
    append_set(lines, "Nset", "ALL_NODES", list(range(1, len(nodes) + 1)))
    append_set(lines, "Nset", "FUEL_INNER_NODES", fuel_inner_nodes)
    append_set(lines, "Nset", "CLAD_OUTER_NODES", clad_outer_nodes)
    append_set(lines, "Nset", "SECONDARY_CONTACT_NODES", secondary_contact_nodes)
    append_set(lines, "Nset", "FUEL_SYMMETRY_Y_NODES", face_nodes(fuel_nodes, fuel_radial, "theta0"))
    append_set(lines, "Nset", "FUEL_SYMMETRY_X_NODES", face_nodes(fuel_nodes, fuel_radial, "theta90"))
    append_set(lines, "Nset", "FUEL_BOTTOM_NODES", face_nodes(fuel_nodes, fuel_radial, "bottom"))
    append_set(lines, "Nset", "CLAD_SYMMETRY_Y_NODES", face_nodes(clad_nodes, clad_radial, "theta0"))
    append_set(lines, "Nset", "CLAD_SYMMETRY_X_NODES", face_nodes(clad_nodes, clad_radial, "theta90"))
    append_set(lines, "Nset", "CLAD_BOTTOM_NODES", face_nodes(clad_nodes, clad_radial, "bottom"))
    lines.extend(
        [
            "*Surface, type=ELEMENT, name=FUEL_CONTACT",
            "FUEL_CONTACT_ELEMENTS, S4",
            "*Surface, type=ELEMENT, name=CLAD_CONTACT",
            "CLAD_CONTACT_ELEMENTS, S6",
            "*Material, name=FUEL_MATERIAL",
            "*Elastic",
            "2.0000000000000000e11, 3.0000000000000000e-1",
            "*Expansion, zero=6.0000000000000000e2",
            "1.0000000000000001e-5",
            "*Conductivity",
            "3.0000000000000000e0",
            "*Density",
            "1.0000000000000000e4",
            "*Specific Heat",
            "3.0000000000000000e2",
            "*Material, name=CLAD_MATERIAL",
            "*Elastic",
            "1.0000000000000000e11, 3.2000000000000001e-1",
            "*Expansion, zero=6.0000000000000000e2",
            "5.0000000000000004e-6",
            "*Conductivity",
            "1.5000000000000000e1",
            "*Density",
            "6.5000000000000000e3",
            "*Specific Heat",
            "3.3000000000000000e2",
            "*Solid Section, elset=FUEL, material=FUEL_MATERIAL",
            ",",
            "*Solid Section, elset=CLAD, material=CLAD_MATERIAL",
            ",",
            "*Surface Interaction, name=FUEL_CLAD_CONTACT",
            "*Surface Behavior, pressure-overclosure=LINEAR",
            "%.16e" % penalty,
            "*Friction, slip tolerance=5.0000000000000001e-3",
            "1.0000000000000001e-1",
            "*Gap Conductance, pressure",
            "1.0000000000000000e2, 0.0",
            "1.1000000000000000e3, 1.0e8",
            "*Contact Pair, interaction=FUEL_CLAD_CONTACT, type=SURFACE TO SURFACE, adjust=0.",
            "CLAD_CONTACT, FUEL_CONTACT",
            "*Amplitude, name=FUEL_TEMPERATURE, time=TOTAL TIME",
            "0.0, 600.0, 10000.0, 1200.0",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 6.0000000000000000e2",
            "*Step, name=PATH, nlgeom=YES, inc=80",
            "*Coupled Temperature-Displacement",
            "%.16e, 1.0000000000000000e4, %.16e, %.16e" % (step, step, step),
            "*Controls, parameters=FIELD, field=DISPLACEMENT",
            "1.0000000000000000e-10, 1.0000000000000000e-10, , , , 1.0000000000000000e-10",
            "*Boundary",
            "FUEL_INNER_NODES, 11, 11, 6.0000000000000000e2",
            "CLAD_OUTER_NODES, 11, 11, 6.0000000000000000e2",
            "FUEL_SYMMETRY_Y_NODES, 2, 2, 0.0",
            "FUEL_SYMMETRY_X_NODES, 1, 1, 0.0",
            "FUEL_BOTTOM_NODES, 3, 3, 0.0",
            "CLAD_SYMMETRY_Y_NODES, 2, 2, 0.0",
            "CLAD_SYMMETRY_X_NODES, 1, 1, 0.0",
            "CLAD_BOTTOM_NODES, 3, 3, 0.0",
            "*Boundary, amplitude=FUEL_TEMPERATURE",
            "FUEL_INNER_NODES, 11, 11, 1.0",
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Element Output, directions=YES",
            "CE, CEEQ, COORD, EE, HFL, IVOL, LE, PE, PEEQ, S, TEMP",
            "*Contact Output",
            "CDISP, CFORCE, CSTRESS, CTANDIR, CSTATUS, HFL",
            "*Output, history, frequency=1",
            "*Energy Output",
            "ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    manifest = ["case\tfuel_radial\tclad_radial\tangular\taxial\ttime_step_s\tpenalty_pa_m\texpected_frames"]
    for case_name in sorted(CASES):
        divisions, step, penalty = CASES[case_name]
        path = directory / (case_name + ".inp")
        path.write_text(deck(case_name, divisions, step, penalty), encoding="utf-8")
        manifest.append(
            "%s\t%d\t%d\t%d\t%d\t%.16g\t%.16g\t%d"
            % ((case_name,) + divisions + (step, penalty, int(round(10000.0 / step))))
        )
        print("wrote %s" % path)
    manifest_path = directory / "b527_cases.tsv"
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print("wrote %s" % manifest_path)
