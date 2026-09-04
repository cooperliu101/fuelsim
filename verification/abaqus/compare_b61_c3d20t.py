#!/usr/bin/env python3
"""Compare final C3D20T plasticity-creep bending fields at nodes and material points."""

import argparse
import csv
import json
import math
from pathlib import Path


COLUMNS = (
    "relative_l2_percent",
    "relative_absolute_peak_percent",
    "maximum_pointwise_relative_percent",
    "maximum_absolute_difference",
    "zero_reference_count",
    "maximum_zero_reference_absolute_difference",
)


def read_rows(path):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source))


def metric(name, actual, reference):
    if len(actual) != len(reference):
        raise RuntimeError("count differs for " + name)
    differences = [abs(a - b) for a, b in zip(actual, reference)]
    reference_squared = sum(value * value for value in reference)
    maximum_reference = max((abs(value) for value in reference), default=0.0)
    pointwise = [difference / abs(value) for difference, value in zip(differences, reference) if value != 0.0]
    zero = [difference for difference, value in zip(differences, reference) if value == 0.0]
    result = {
        "relative_l2_percent": 100.0 * math.sqrt(sum(value * value for value in differences) / reference_squared)
        if reference_squared
        else float("nan"),
        "relative_absolute_peak_percent": 100.0 * max(differences, default=0.0) / maximum_reference
        if maximum_reference
        else float("nan"),
        "maximum_pointwise_relative_percent": 100.0 * max(pointwise, default=0.0),
        "maximum_absolute_difference": max(differences, default=0.0),
        "zero_reference_count": len(zero),
        "maximum_zero_reference_absolute_difference": max(zero, default=0.0),
    }
    print("%s: " % name + " ".join("%s=%.9g" % (column, result[column]) for column in COLUMNS))
    return result


def vector_metric(name, actual, reference):
    difference = [
        math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
        for left, right in zip(actual, reference)
    ]
    reference_norm = [math.sqrt(sum(value * value for value in vector)) for vector in reference]
    reference_squared = sum(value * value for value in reference_norm)
    maximum_reference = max(reference_norm, default=0.0)
    pointwise = [delta / value for delta, value in zip(difference, reference_norm) if value != 0.0]
    zero = [delta for delta, value in zip(difference, reference_norm) if value == 0.0]
    result = {
        "relative_l2_percent": 100.0 * math.sqrt(sum(value * value for value in difference) / reference_squared)
        if reference_squared else float("nan"),
        "relative_absolute_peak_percent": 100.0 * max(difference, default=0.0) / maximum_reference
        if maximum_reference else float("nan"),
        "maximum_pointwise_relative_percent": 100.0 * max(pointwise, default=0.0),
        "maximum_absolute_difference": max(difference, default=0.0),
        "zero_reference_count": len(zero),
        "maximum_zero_reference_absolute_difference": max(zero, default=0.0),
    }
    print("%s: " % name + " ".join("%s=%.9g" % (column, result[column]) for column in COLUMNS))
    return result


parser = argparse.ArgumentParser()
parser.add_argument("mesh_manifest")
parser.add_argument("fuelsim_nodal")
parser.add_argument("fuelsim_material")
parser.add_argument("abaqus_nodal")
parser.add_argument("abaqus_integration")
parser.add_argument("comparison_output")
args = parser.parse_args()

with Path(args.mesh_manifest).open() as source:
    mesh = json.load(source)
corner_labels = {node for element in mesh["elements"] for node in element["nodes"][:8]}
minimum_x = min(node["coordinates"][0] for node in mesh["nodes"])
left_labels = {
    node["label"] for node in mesh["nodes"] if abs(node["coordinates"][0] - minimum_x) < 1.0e-14
}
fuelsim = {int(row["id"]): row for row in read_rows(args.fuelsim_nodal)}
abaqus = {int(row["id"]): row for row in read_rows(args.abaqus_nodal)}
if set(fuelsim) != set(abaqus):
    raise RuntimeError("Fuelsim and Abaqus node labels differ")

metrics = []
labels = sorted(corner_labels)
metrics.append(
    (
        "temperature",
        metric(
            "temperature",
            [float(fuelsim[label]["temperature"]) for label in labels],
            [float(abaqus[label]["temperature"]) for label in labels],
        ),
    )
)
free_labels = sorted(set(fuelsim) - left_labels)
fuelsim_vectors = [
    tuple(float(fuelsim[label]["displacement_" + component]) for component in "xyz")
    for label in free_labels
]
abaqus_vectors = [
    tuple(float(abaqus[label]["displacement_" + component]) for component in "xyz")
    for label in free_labels
]
metrics.append(
    (
        "free_node_displacement_vector",
        vector_metric("free_node_displacement_vector", fuelsim_vectors, abaqus_vectors),
    )
)
for index, component in enumerate("xyz"):
    metric(
        "diagnostic_displacement_" + component,
        [value[index] for value in fuelsim_vectors],
        [value[index] for value in abaqus_vectors],
    )

fuelsim_points = {
    (int(row["element"]), int(row["integration_point"])): row for row in read_rows(args.fuelsim_material)
}
abaqus_points = {
    (int(row["element"]), int(row["integration_point"])): row for row in read_rows(args.abaqus_integration)
}
if set(fuelsim_points) != set(abaqus_points) or len(fuelsim_points) != 27 * len(mesh["elements"]):
    raise RuntimeError("Fuelsim and Abaqus material-point labels differ")
keys = sorted(fuelsim_points)
for field in ("vonmises_stress", "peeq", "ceeq"):
    metrics.append(
        (
            field,
            metric(
                field,
                [float(fuelsim_points[key][field]) for key in keys],
                [float(abaqus_points[key][field]) for key in keys],
            ),
        )
    )

with Path(args.comparison_output).open("w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("field",) + COLUMNS)
    for name, result in metrics:
        writer.writerow((name,) + tuple(result[column] for column in COLUMNS))

failed = [
    (name, column, result[column])
    for name, result in metrics
    for column in COLUMNS[:3]
    if result[column] >= 0.5
]
if failed:
    raise SystemExit(
        "B6.1 C3D20T comparison exceeds 0.5 percent: "
        + ", ".join("%s %s %.9g" % item for item in failed)
    )
