#!/usr/bin/env python3
"""Generate the 376-degree-of-freedom B5.49 Abaqus C3D20T isolation case."""

from pathlib import Path

from exodus_to_abaqus import convert


JOB = "b549_small_c3d20t"


def append_labels(lines, name, labels):
    lines.append("*Nset, nset=%s" % name)
    for begin in range(0, len(labels), 16):
        lines.append(", ".join(str(value) for value in labels[begin : begin + 16]))


def amplitude(lines, name, final_value):
    lines.extend(
        [
            "*Amplitude, name=%s, time=TOTAL TIME" % name,
            "0, 0, 0.4, %.16g" % final_value,
        ]
    )


def deck(mesh):
    corner_labels = sorted(
        set(node for element in mesh["elements"] for node in element["nodes"][:8])
    )
    primary_labels = set(
        node
        for block in mesh["blocks"]
        if block["name"] == "PRIMARY"
        for element in block["elements"]
        for node in element["nodes"]
    )
    secondary_labels = set(
        node
        for block in mesh["blocks"]
        if block["name"] == "SECONDARY"
        for element in block["elements"]
        for node in element["nodes"]
    )
    coordinates = dict((node["label"], node["coordinates"]) for node in mesh["nodes"])
    primary_outer = sorted(
        label for label in primary_labels if abs(coordinates[label][0]) < 1.0e-14
    )
    secondary_outer = sorted(
        label for label in secondary_labels if abs(coordinates[label][0] - 2.0) < 1.0e-14
    )
    primary_contact = sorted(
        label for label in primary_labels if abs(coordinates[label][0] - 1.0) < 1.0e-14
    )
    secondary_contact = sorted(
        label for label in secondary_labels if abs(coordinates[label][0] - 1.0) < 1.0e-14
    )
    primary_temperature = sorted(set(primary_outer).intersection(corner_labels))
    secondary_temperature = sorted(set(secondary_outer).intersection(corner_labels))
    primary_contact_temperature = sorted(set(primary_contact).intersection(corner_labels))
    secondary_contact_temperature = sorted(set(secondary_contact).intersection(corner_labels))
    if (
        len(mesh["nodes"]) != 112
        or len(mesh["elements"]) != 8
        or len(corner_labels) != 40
        or len(primary_outer) != 8
        or len(secondary_outer) != 8
        or len(primary_temperature) != 4
        or len(secondary_temperature) != 4
        or len(primary_contact) != 8
        or len(secondary_contact) != 8
        or len(primary_contact_temperature) != 4
        or len(secondary_contact_temperature) != 4
    ):
        raise RuntimeError("unexpected B5.49 mixed-order mesh topology")

    lines = [
        "*Heading",
        "** B5.49 small C3D20T finite-strain thermo-inelastic volume-element isolation case.",
        "** Eight elements, 112 displacement nodes, 40 temperature nodes, and 376 total degrees of freedom.",
        "** SI units: metre, second, kelvin, pascal, watt.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Include, input=%s_mesh.inc" % JOB,
    ]
    append_labels(lines, "ALL_TEMPERATURE_NODES", corner_labels)
    append_labels(lines, "PRIMARY_OUTER_NODES", primary_outer)
    append_labels(lines, "PRIMARY_OUTER_TEMPERATURE", primary_temperature)
    append_labels(lines, "SECONDARY_OUTER_NODES", secondary_outer)
    append_labels(lines, "SECONDARY_OUTER_TEMPERATURE", secondary_temperature)
    append_labels(lines, "PRIMARY_CONTACT_NODES", primary_contact)
    append_labels(lines, "PRIMARY_CONTACT_TEMPERATURE", primary_contact_temperature)
    append_labels(lines, "SECONDARY_CONTACT_NODES", secondary_contact)
    append_labels(lines, "SECONDARY_CONTACT_TEMPERATURE", secondary_contact_temperature)
    lines.extend(
        [
            "*Material, name=PRIMARY_MATERIAL",
            "*Elastic",
            "1.2e8, 0.28",
            "*Expansion, zero=300",
            "8e-6",
            "*Conductivity",
            "15",
            "*Density",
            "100",
            "*Specific Heat",
            "1",
            "*Material, name=SECONDARY_MATERIAL",
            "*Elastic",
            "1e8, 0.3",
            "*Expansion, zero=300",
            "1e-5",
            "*Plastic",
            "1e5, 0",
            "1.01e7, 1",
            "*Creep, law=TIME",
            "1.25e-20, 3, 0",
            "*Conductivity",
            "10",
            "*Density",
            "100",
            "*Specific Heat",
            "1",
            "*Solid Section, elset=PRIMARY, material=PRIMARY_MATERIAL",
            ",",
            "*Solid Section, elset=SECONDARY, material=SECONDARY_MATERIAL",
            ",",
        ]
    )
    amplitude(lines, "SECONDARY_TEMPERATURE", 400.0)
    lines[-1] = "0, 300, 0.4, 400"
    amplitude(lines, "OUTER_PRESSURE", 1.75e5)
    amplitude(lines, "TANGENTIAL_Y", 0.0)
    amplitude(lines, "TANGENTIAL_Z", 0.0)
    lines.extend(
        [
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_TEMPERATURE_NODES, 300",
            "*Step, name=PATH, nlgeom=YES, inc=40",
            "*Coupled Temperature-Displacement",
            "0.02, 0.4, 0.02, 0.02",
            "*Controls, parameters=FIELD, field=DISPLACEMENT",
            "1e-10, 1e-10, , , , 1e-10",
            "*Boundary",
            "PRIMARY_OUTER_NODES, 1, 3, 0",
            "PRIMARY_OUTER_TEMPERATURE, 11, 11, 300",
            "SECONDARY_CONTACT_NODES, 1, 1, 0",
            "SECONDARY_CONTACT_TEMPERATURE, 11, 11, 300",
            "*Boundary, amplitude=SECONDARY_TEMPERATURE",
            "PRIMARY_CONTACT_TEMPERATURE, 11, 11, 1",
            "SECONDARY_OUTER_TEMPERATURE, 11, 11, 1",
            "*Boundary, amplitude=TANGENTIAL_Y",
            "SECONDARY_OUTER_NODES, 2, 2, 1",
            "*Boundary, amplitude=TANGENTIAL_Z",
            "SECONDARY_OUTER_NODES, 3, 3, 1",
            "*Dsload, follower=YES, amplitude=OUTER_PRESSURE",
            "SECONDARY_OUTER, P, 1",
            "*Output, field, frequency=20",
            "*Node Output",
            "COORD, NT, RF, RFL, U",
            "*Element Output, directions=YES",
            "CEEQ, COORD, IVOL, PEEQ, S, TEMP",
            "*Output, history, frequency=20",
            "*Energy Output",
            "ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK",
            "*End Step",
        ]
    )
    return "\n".join(lines) + "\n"


def main():
    directory = Path(__file__).resolve().parent
    mesh_path = directory / (JOB + "_mesh.e")
    include_path = directory / (JOB + "_mesh.inc")
    manifest_path = directory / (JOB + "_mesh.json")
    mesh = convert(mesh_path, include_path, manifest_path, "C3D20T")
    input_path = directory / (JOB + ".inp")
    input_path.write_text(deck(mesh), encoding="ascii")
    print("wrote %s, %s, and %s" % (input_path, include_path, manifest_path))


if __name__ == "__main__":
    main()
