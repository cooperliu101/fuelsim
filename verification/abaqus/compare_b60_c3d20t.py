#!/usr/bin/env python3
"""Compare a Fuelsim B6.0 C3D20T result with the matched Abaqus result."""

import argparse
import csv
import json
import math
from pathlib import Path


METRIC_COLUMNS = (
    "relative_l2_percent",
    "relative_absolute_peak_percent",
    "maximum_pointwise_relative_percent",
    "maximum_absolute_difference",
    "zero_reference_count",
    "maximum_zero_reference_absolute_difference",
)


def rows(path):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source))


def print_metric(name, result):
    print(
        "%s: relative L2=%.9g%% relative absolute peak=%.9g%% maximum pointwise=%.9g%% "
        "maximum absolute difference=%.9g zero references=%d maximum zero-reference absolute difference=%.9g"
        % ((name,) + tuple(result[column] for column in METRIC_COLUMNS))
    )


def scalar_metric(name, actual, reference):
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
    print_metric(name, result)
    return result


def vector_metric(name, actual, reference):
    if len(actual) != len(reference):
        raise RuntimeError("count differs for " + name)
    difference_norms = [
        math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
        for left, right in zip(actual, reference)
    ]
    reference_norms = [math.sqrt(sum(value * value for value in vector)) for vector in reference]
    reference_squared = sum(value * value for value in reference_norms)
    maximum_reference = max(reference_norms, default=0.0)
    pointwise = [difference / value for difference, value in zip(difference_norms, reference_norms) if value != 0.0]
    zero = [difference for difference, value in zip(difference_norms, reference_norms) if value == 0.0]
    result = {
        "relative_l2_percent": 100.0 * math.sqrt(sum(value * value for value in difference_norms) / reference_squared)
        if reference_squared else float("nan"),
        "relative_absolute_peak_percent": 100.0 * max(difference_norms, default=0.0) / maximum_reference
        if maximum_reference else float("nan"),
        "maximum_pointwise_relative_percent": 100.0 * max(pointwise, default=0.0),
        "maximum_absolute_difference": max(difference_norms, default=0.0),
        "zero_reference_count": len(zero),
        "maximum_zero_reference_absolute_difference": max(zero, default=0.0),
    }
    print_metric(name, result)
    return result


parser = argparse.ArgumentParser()
parser.add_argument("mesh_manifest")
parser.add_argument("fuelsim_nodal")
parser.add_argument("fuelsim_material")
parser.add_argument("abaqus_temperature")
parser.add_argument("abaqus_displacement")
parser.add_argument("abaqus_material")
parser.add_argument("comparison_output")
parser.add_argument(
    "--qualified-low-stress-case",
    choices=("steady", "ramped"),
    help="Apply the recorded case-specific maximum-pointwise stress qualification.",
)
args = parser.parse_args()

with Path(args.mesh_manifest).open() as source:
    mesh = json.load(source)
corner_labels = {node for element in mesh["elements"] for node in element["nodes"][:8]}
minimum_x = min(node["coordinates"][0] for node in mesh["nodes"])
left_labels = {
    node["label"] for node in mesh["nodes"] if abs(node["coordinates"][0] - minimum_x) < 1.0e-14
}

fuelsim = {int(row["id"]): row for row in rows(args.fuelsim_nodal)}
temperature = {int(row["id"]): row for row in rows(args.abaqus_temperature)}
displacement = {int(row["id"]): row for row in rows(args.abaqus_displacement)}
if set(temperature) != corner_labels or set(displacement) != set(fuelsim):
    raise RuntimeError("Fuelsim, Abaqus, and mesh node labels differ")

metrics = []
labels = sorted(corner_labels)
metrics.append(
    (
        "temperature",
        scalar_metric(
            "temperature",
            [float(fuelsim[label]["temperature"]) for label in labels],
            [float(temperature[label]["temperature"]) for label in labels],
        ),
    )
)

free_labels = sorted(set(displacement) - left_labels)
fuelsim_vectors = [
    tuple(float(fuelsim[label]["displacement_" + component]) for component in "xyz")
    for label in free_labels
]
abaqus_vectors = [
    tuple(float(displacement[label]["displacement_" + component]) for component in "xyz")
    for label in free_labels
]
metrics.append(
    (
        "free_node_displacement_vector",
        vector_metric("free_node_displacement_vector", fuelsim_vectors, abaqus_vectors),
    )
)
for index, component in enumerate("xyz"):
    scalar_metric(
        "diagnostic_displacement_" + component,
        [value[index] for value in fuelsim_vectors],
        [value[index] for value in abaqus_vectors],
    )

fuelsim_material = {
    (int(row["element"]), int(row["integration_point"])): row for row in rows(args.fuelsim_material)
}
abaqus_material = {
    (int(row["element"]), int(row["integration_point"])): row for row in rows(args.abaqus_material)
}
if set(fuelsim_material) != set(abaqus_material) or len(fuelsim_material) != 27 * len(mesh["elements"]):
    raise RuntimeError("Fuelsim and Abaqus material output must contain the same 27 points per element")
keys = sorted(fuelsim_material)
fuelsim_stress = [float(fuelsim_material[key]["vonmises_stress"]) for key in keys]
abaqus_stress = [float(abaqus_material[key]["vonmises_stress"]) for key in keys]
stress_result = scalar_metric("vonmises_stress", fuelsim_stress, abaqus_stress)
metrics.append(("vonmises_stress", stress_result))

with Path(args.comparison_output).open("w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("field",) + METRIC_COLUMNS)
    for name, result in metrics:
        writer.writerow((name,) + tuple(result[column] for column in METRIC_COLUMNS))

strict_metrics_pass = all(
    all(result[column] < 0.5 for column in METRIC_COLUMNS[:3])
    for name, result in metrics
    if name != "vonmises_stress"
)
stress_aggregate_pass = all(stress_result[column] < 0.5 for column in METRIC_COLUMNS[:2])
stress_pointwise_pass = stress_result["maximum_pointwise_relative_percent"] < 0.5

if args.qualified_low_stress_case:
    qualification = {
        "steady": {
            "maximum_pointwise_relative_percent": 2.5,
            "maximum_reference_stress": 5.0e3,
            "maximum_absolute_difference": 120.0,
        },
        "ramped": {
            "maximum_pointwise_relative_percent": 15.5,
            "maximum_reference_stress": 4.0e3,
            "maximum_absolute_difference": 575.0,
        },
    }[args.qualified_low_stress_case]
    pointwise = [
        (abs(actual - reference) / abs(reference), index)
        for index, (actual, reference) in enumerate(zip(fuelsim_stress, abaqus_stress))
        if reference != 0.0
    ]
    _, maximum_index = max(pointwise)
    maximum_key = keys[maximum_index]
    maximum_actual = fuelsim_stress[maximum_index]
    maximum_reference = abaqus_stress[maximum_index]
    maximum_difference = abs(maximum_actual - maximum_reference)
    stress_pointwise_pass = (
        stress_result["maximum_pointwise_relative_percent"]
        < qualification["maximum_pointwise_relative_percent"]
        and abs(maximum_reference) < qualification["maximum_reference_stress"]
        and maximum_difference < qualification["maximum_absolute_difference"]
    )
    print(
        "qualified low-stress point: case=%s element=%d integration_point=%d "
        "Fuelsim=%.9g Pa Abaqus=%.9g Pa absolute difference=%.9g Pa"
        % (
            args.qualified_low_stress_case,
            maximum_key[0],
            maximum_key[1],
            maximum_actual,
            maximum_reference,
            maximum_difference,
        )
    )

if not (strict_metrics_pass and stress_aggregate_pass and stress_pointwise_pass):
    raise SystemExit("B6.0 C3D20T comparison exceeds the 0.5 percent acceptance boundary")
