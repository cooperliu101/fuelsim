#!/usr/bin/env python3
"""Compare the H20.42 Abaqus pressure output with the reconstructed recovery rule."""

from pathlib import Path
import csv
import math
import sys


CONTACT_NODES = [
    71, 74, 75, 78, 82, 83, 86, 90, 92, 94, 97, 99, 102, 103, 106,
    107, 110, 114, 116, 118, 121, 122, 125, 126, 129, 133, 135, 137, 140,
]
PROJECTOR = [
    [17, -3, 1, -3, 7, -1, -1, 7],
    [-3, 17, -3, 1, 7, 7, -1, -1],
    [1, -3, 17, -3, -1, 7, 7, -1],
    [-3, 1, -3, 17, -1, -1, 7, 7],
    [7, 7, -1, -1, 7, 3, -1, 3],
    [-1, 7, 7, -1, 3, 7, 3, -1],
    [-1, -1, 7, 7, -1, 3, 7, 3],
    [7, -1, -1, 7, 3, -1, 3, 7],
]


def secondary_faces(path):
    elements = []
    reading = False
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if line.startswith("*"):
            reading = line.lower() == "*element, type=c3d20, elset=secondary"
            continue
        if reading and line:
            elements.append([int(value.strip()) for value in line.split(",")])
    face_indices = [0, 3, 7, 4, 11, 19, 15, 16]
    return [[element[index + 1] for index in face_indices] for element in elements]


def read_states(path):
    states = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            states.setdefault(row["state"], {})[int(row["label"])] = {
                "gap": float(row["copen"]),
                "pressure": float(row["cpress"]),
            }
    for name, rows in states.items():
        if sorted(rows) != sorted(CONTACT_NODES):
            raise RuntimeError("state %s does not contain every H20.42 contact node" % name)
    return states


def raw_recovery(rows, faces):
    constraint = dict((node, max(-1.0e11 * value["gap"], 0.0)) for node, value in rows.items())
    recovered = dict((node, 0.0) for node in rows)
    count = dict((node, 0) for node in rows)
    for face in faces:
        for output, node in enumerate(face):
            recovered[node] += sum(PROJECTOR[output][input_node] * constraint[face[input_node]]
                                   for input_node in range(8)) / 24.0
            count[node] += 1
    return constraint, dict((node, recovered[node] / count[node]) for node in rows)


def recovered_pressure(rows, faces):
    constraint, raw = raw_recovery(rows, faces)
    minimum = min(constraint.values())
    maximum = max(constraint.values())
    return dict((node, min(max(value, minimum), maximum)) for node, value in raw.items())


def metrics(rows, predicted):
    differences = [predicted[node] - rows[node]["pressure"] for node in CONTACT_NODES]
    references = [rows[node]["pressure"] for node in CONTACT_NODES]
    relative_l2 = math.sqrt(sum(value * value for value in differences) /
                            sum(value * value for value in references))
    relative_peak = max(abs(value) for value in differences) / max(abs(value) for value in references)
    pointwise = max(abs(differences[index] / references[index]) for index in range(len(references))
                    if references[index] != 0.0)
    zero_difference = max([abs(differences[index]) for index in range(len(references))
                           if references[index] == 0.0] or [0.0])
    return relative_l2, relative_peak, pointwise, zero_difference


def main():
    if len(sys.argv) not in (2, 3):
        raise RuntimeError("usage: compare_h20_42.py <result.csv> [metrics.tsv]")
    directory = Path(__file__).resolve().parent
    faces = secondary_faces(directory / "h20_41_hex20_finite_sliding_partial_contact.inp")
    states = read_states(Path(sys.argv[1]))
    if len(states) != 128:
        raise RuntimeError("H20.42 has %d states, expected 128" % len(states))

    baseline = states["BASE"]
    baseline_pressure = dict((node, baseline[node]["pressure"]) for node in CONTACT_NODES)
    maximum_operator_difference = 0.0
    for column, input_node in enumerate(CONTACT_NODES, 1):
        expected_constraint = dict((node, 1.0 if node == input_node else 0.0) for node in CONTACT_NODES)
        expected_rows = dict((node, {"gap": -expected_constraint[node] / 1.0e11, "pressure": 0.0})
                             for node in CONTACT_NODES)
        expected = raw_recovery(expected_rows, faces)[1]
        actual = states["N%02d_PLUS" % column]
        for output_node in CONTACT_NODES:
            coefficient = (actual[output_node]["pressure"] - baseline_pressure[output_node]) / 1.0e6
            maximum_operator_difference = max(maximum_operator_difference,
                                              abs(coefficient - expected[output_node]))

    state_metrics = dict((name, metrics(rows, recovered_pressure(rows, faces))) for name, rows in states.items())
    unbounded_metrics = dict((name, metrics(rows, raw_recovery(rows, faces)[1])) for name, rows in states.items())
    maxima = [max(values[index] for values in state_metrics.values()) for index in range(4)]
    unbounded_maxima = [max(values[index] for values in unbounded_metrics.values()) for index in range(4)]
    h2041 = state_metrics["H2041"]
    partial_names = [name for name in states if name.startswith("RANDOM_")]
    partial_maxima = [max(state_metrics[name][index] for name in partial_names) for index in range(4)]
    active_names = [name for name in states if name.startswith("ACTIVE_") or name.startswith("EXTREME_")]
    active_maxima = [max(state_metrics[name][index] for name in active_names) for index in range(4)]
    open_gap_difference = 0.0
    for index in range(1, 5):
        first = states["RANDOM_%02d" % index]
        second = states["RANDOM_%02d_OPEN2" % index]
        open_gap_difference = max(open_gap_difference,
                                  max(abs(first[node]["pressure"] - second[node]["pressure"])
                                      for node in CONTACT_NODES))

    values = [
        ("state_count", float(len(states)), "count"),
        ("assembled_projector_maximum_absolute_coefficient_difference", maximum_operator_difference, "dimensionless"),
        ("all_state_maximum_relative_l2", 100.0 * maxima[0], "percent"),
        ("all_state_maximum_relative_absolute_peak", 100.0 * maxima[1], "percent"),
        ("all_state_maximum_pointwise_relative", 100.0 * maxima[2], "percent"),
        ("all_state_maximum_zero_reference_absolute_difference", maxima[3], "Pa"),
        ("unbounded_all_state_maximum_relative_l2", 100.0 * unbounded_maxima[0], "percent"),
        ("unbounded_all_state_maximum_relative_absolute_peak", 100.0 * unbounded_maxima[1], "percent"),
        ("unbounded_all_state_maximum_pointwise_relative", 100.0 * unbounded_maxima[2], "percent"),
        ("unbounded_all_state_maximum_zero_reference_absolute_difference", unbounded_maxima[3], "Pa"),
        ("h2041_relative_l2", 100.0 * h2041[0], "percent"),
        ("h2041_relative_absolute_peak", 100.0 * h2041[1], "percent"),
        ("h2041_maximum_pointwise_relative", 100.0 * h2041[2], "percent"),
        ("partial_state_maximum_relative_l2", 100.0 * partial_maxima[0], "percent"),
        ("active_state_maximum_relative_l2", 100.0 * active_maxima[0], "percent"),
        ("alternate_open_gap_maximum_pressure_difference", open_gap_difference, "Pa"),
    ]
    for name, value, unit in values:
        print("%s=%.16e %s" % (name, value, unit))
    if len(sys.argv) == 3:
        with Path(sys.argv[2]).open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["quantity", "value", "unit"])
            for row in values:
                writer.writerow([row[0], "%.16e" % row[1], row[2]])


if __name__ == "__main__":
    main()
