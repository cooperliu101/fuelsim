#!/usr/bin/env python3
"""Generate fixed-clearance C3D20 pressure-recovery probes for H20.42."""

from pathlib import Path
import csv
import random
import sys


CONTACT_NODES = [
    71, 74, 75, 78, 82, 83, 86, 90, 92, 94, 97, 99, 102, 103, 106,
    107, 110, 114, 116, 118, 121, 122, 125, 126, 129, 133, 135, 137, 140,
]


def read_mesh(path):
    coordinates = {}
    primary = []
    secondary = []
    mode = None
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if line.startswith("*"):
            lower = line.lower()
            if lower == "*node":
                mode = "node"
            elif lower == "*element, type=c3d20, elset=primary":
                mode = "primary"
            elif lower == "*element, type=c3d20, elset=secondary":
                mode = "secondary"
            else:
                mode = None
            continue
        if not line:
            continue
        values = [value.strip() for value in line.split(",")]
        if mode == "node":
            coordinates[int(values[0])] = tuple(float(value) for value in values[1:4])
        elif mode == "primary":
            primary.append([int(value) for value in values])
        elif mode == "secondary":
            secondary.append([int(value) for value in values])
    return coordinates, primary, secondary


def h2041_clearance(path):
    result = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            label = int(row["id"]) + 1
            if label in CONTACT_NODES:
                result[label] = float(row["gap"])
    if len(result) != len(CONTACT_NODES):
        raise RuntimeError("H20.41 reference does not contain every secondary contact node")
    return result


def build_states(reference_path):
    base = -2.0e-4
    delta = 1.0e-5
    states = [("BASE", dict((label, base) for label in CONTACT_NODES))]
    for index, label in enumerate(CONTACT_NODES, 1):
        plus = dict((node, base) for node in CONTACT_NODES)
        minus = dict(plus)
        plus[label] = base - delta
        minus[label] = base + delta
        states.append(("N%02d_PLUS" % index, plus))
        states.append(("N%02d_MINUS" % index, minus))
    states.append(("H2041", h2041_clearance(reference_path)))

    rng = random.Random(2047)
    for sample in range(12):
        values = {}
        active_pressure = {}
        for label in CONTACT_NODES:
            if rng.random() < 0.55:
                pressure = (0.15 + 1.85 * rng.random()) * 1.0e7
                values[label] = -pressure / 1.0e11
                active_pressure[label] = pressure
            else:
                values[label] = (0.2 + 1.8 * rng.random()) * 1.0e-4
                active_pressure[label] = 0.0
        states.append(("RANDOM_%02d" % (sample + 1), values))
        if sample < 4:
            alternate = {}
            for label in CONTACT_NODES:
                alternate[label] = (-active_pressure[label] / 1.0e11 if active_pressure[label] > 0.0
                                    else (2.2 + 1.8 * rng.random()) * 1.0e-4)
            states.append(("RANDOM_%02d_OPEN2" % (sample + 1), alternate))

    for sample in range(12):
        values = {}
        for label in CONTACT_NODES:
            pressure = (0.05 + 2.95 * rng.random()) * 1.0e7
            values[label] = -pressure / 1.0e11
        states.append(("ACTIVE_%02d" % (sample + 1), values))
    for sample in range(40):
        values = {}
        for label in CONTACT_NODES:
            pressure = (0.001 + 3.999 * rng.random()) * 1.0e7
            values[label] = -pressure / 1.0e11
        states.append(("EXTREME_%02d" % (sample + 1), values))
    return states


def main():
    if len(sys.argv) != 2:
        raise RuntimeError("usage: generate_h20_42.py <probe-directory>")
    source_directory = Path(__file__).resolve().parent
    probe_directory = Path(sys.argv[1]).resolve()
    probe_directory.mkdir(parents=True, exist_ok=True)
    source = source_directory / "h20_41_hex20_finite_sliding_partial_contact.inp"
    reference = source_directory / "h20_41_hex20_finite_sliding_partial_contact_contact.csv"
    target = probe_directory / "h20_42_hex20_pressure_recovery_probe.inp"
    manifest_path = probe_directory / "h20_42_hex20_pressure_recovery_manifest.csv"
    coordinates, primary_elements, secondary_elements = read_mesh(source)
    states = build_states(reference)

    pieces = [
        "*Heading\n",
        "** H20.42: fixed-clearance C3D20 contact-pressure recovery probe.\n",
        "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n",
        "*Node\n",
    ]
    all_nodes = []
    all_primary_elements = []
    all_secondary_elements = []
    for copy_index, (name, values) in enumerate(states):
        node_offset = 1000 * copy_index
        element_offset = 100 * copy_index
        y_offset = 2.0 * copy_index
        for label in sorted(coordinates):
            x, y, z = coordinates[label]
            output_label = node_offset + label
            all_nodes.append(output_label)
            pieces.append("%d, %.17g, %.17g, %.17g\n" % (output_label, x, y + y_offset, z))
        for element in primary_elements:
            all_primary_elements.append([element_offset + element[0]] +
                                        [node_offset + node for node in element[1:]])
        for element in secondary_elements:
            all_secondary_elements.append([element_offset + element[0]] +
                                          [node_offset + node for node in element[1:]])

    pieces.append("*Element, type=C3D20, elset=PRIMARY_ALL\n")
    for element in all_primary_elements:
        pieces.append(", ".join(str(value) for value in element) + "\n")
    pieces.append("*Element, type=C3D20, elset=SECONDARY_ALL\n")
    for element in all_secondary_elements:
        pieces.append(", ".join(str(value) for value in element) + "\n")

    for copy_index, (name, values) in enumerate(states):
        element_offset = 100 * copy_index
        node_offset = 1000 * copy_index
        pieces.extend([
            "*Elset, elset=P%03d\n" % copy_index,
            ", ".join(str(element_offset + element[0]) for element in primary_elements) + "\n",
            "*Elset, elset=S%03d\n" % copy_index,
            ", ".join(str(element_offset + element[0]) for element in secondary_elements) + "\n",
            "*Surface, type=ELEMENT, name=PC%03d\n" % copy_index,
            "P%03d, S4\n" % copy_index,
            "*Surface, type=ELEMENT, name=SC%03d\n" % copy_index,
            "S%03d, S6\n" % copy_index,
            "*Clearance, slave=SC%03d, master=PC%03d, tabular\n" % (copy_index, copy_index),
        ])
        for label in CONTACT_NODES:
            pieces.append("%d, %.17g\n" % (node_offset + label, values[label]))

    pieces.extend([
        "*Material, name=PRIMARY_MATERIAL\n",
        "*Elastic\n",
        "2000000000, 0.25\n",
        "*Material, name=SECONDARY_MATERIAL\n",
        "*Elastic\n",
        "650000000, 0.25\n",
        "*Solid Section, elset=PRIMARY_ALL, material=PRIMARY_MATERIAL\n",
        ",\n",
        "*Solid Section, elset=SECONDARY_ALL, material=SECONDARY_MATERIAL\n",
        ",\n",
        "*Surface Interaction, name=PENALTY_CONTACT\n",
        "*Surface Behavior, pressure-overclosure=LINEAR\n",
        "100000000000,\n",
    ])
    for copy_index in range(len(states)):
        pieces.extend([
            "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.\n",
            "SC%03d, PC%03d\n" % (copy_index, copy_index),
        ])
    pieces.append("*Nset, nset=ALL_FIXED\n")
    for start in range(0, len(all_nodes), 16):
        pieces.append(", ".join(str(value) for value in all_nodes[start:start + 16]) + "\n")
    pieces.extend([
        "*Step, name=LOAD, nlgeom=NO, inc=20\n",
        "*Static\n",
        "1., 1., 1.e-10, 1.\n",
        "*Boundary\n",
        "ALL_FIXED, 1, 3, 0.\n",
        "*Output, field, frequency=1\n",
        "*Contact Output\n",
        "CSTRESS, CDISP\n",
        "*End Step\n",
    ])
    target.write_text("".join(pieces), encoding="ascii")
    with manifest_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["copy", "state"])
        for copy_index, (name, values) in enumerate(states):
            writer.writerow([copy_index, name])
    print("wrote %s with %d independent states" % (target, len(states)))


if __name__ == "__main__":
    main()
