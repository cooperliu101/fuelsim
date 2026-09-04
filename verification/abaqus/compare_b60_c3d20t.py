#!/usr/bin/env python3
"""Compare a Fuelsim B6.0 C3D20T result with an Abaqus result."""

import argparse
import csv
import math
from pathlib import Path


def rows(path):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source))


def metric(name, actual, reference):
    if len(actual) != len(reference):
        raise RuntimeError("count differs for " + name)
    difference_squared = sum((a - b) ** 2 for a, b in zip(actual, reference))
    reference_squared = sum(b * b for b in reference)
    maximum_difference = max((abs(a - b) for a, b in zip(actual, reference)), default=0.0)
    maximum_reference = max((abs(b) for b in reference), default=0.0)
    pointwise = [abs(a - b) / abs(b) for a, b in zip(actual, reference) if b != 0.0]
    zero = [abs(a - b) for a, b in zip(actual, reference) if b == 0.0]
    result = {
        "relative_l2_percent": 100.0 * math.sqrt(difference_squared / reference_squared)
        if reference_squared else float("nan"),
        "relative_absolute_peak_percent": 100.0 * maximum_difference / maximum_reference
        if maximum_reference else float("nan"),
        "maximum_pointwise_relative_percent": 100.0 * max(pointwise, default=0.0),
        "maximum_absolute_difference": maximum_difference,
        "zero_reference_count": len(zero),
        "maximum_zero_reference_absolute_difference": max(zero, default=0.0),
    }
    print(
        "%s: relative L2=%.9g%% relative absolute peak=%.9g%% maximum pointwise=%.9g%% "
        "maximum absolute difference=%.9g zero references=%d maximum zero-reference absolute difference=%.9g"
        % (
            name,
            result["relative_l2_percent"],
            result["relative_absolute_peak_percent"],
            result["maximum_pointwise_relative_percent"],
            result["maximum_absolute_difference"],
            result["zero_reference_count"],
            result["maximum_zero_reference_absolute_difference"],
        )
    )
    return result


parser = argparse.ArgumentParser()
parser.add_argument("fuelsim_nodal")
parser.add_argument("abaqus_temperature")
parser.add_argument("abaqus_displacement")
parser.add_argument("abaqus_material")
parser.add_argument("comparison_output")
args = parser.parse_args()

fuelsim = {int(row["id"]): row for row in rows(args.fuelsim_nodal)}
temperature = {int(row["id"]): row for row in rows(args.abaqus_temperature)}
displacement = {int(row["id"]): row for row in rows(args.abaqus_displacement)}
if set(displacement) != set(fuelsim):
    raise RuntimeError("Fuelsim and Abaqus displacement node labels differ")

metrics = []
metrics.append(
    ("temperature", metric(
        "temperature",
        [float(fuelsim[label]["temperature"]) for label in sorted(temperature)],
        [float(temperature[label]["temperature"]) for label in sorted(temperature)],
    ))
)
for component in ("x", "y", "z"):
    name = "displacement_" + component
    metrics.append(
        (name, metric(
            name,
            [float(fuelsim[label]["displacement_" + component]) for label in sorted(displacement)],
            [float(displacement[label]["displacement_" + component]) for label in sorted(displacement)],
        ))
    )

material = rows(args.abaqus_material)
if len(material) != 27 * len({(int(row["element"]),) for row in material}):
    raise RuntimeError("Abaqus material output does not contain 27 points per element")

with Path(args.comparison_output).open("w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("field", "relative_l2_percent", "relative_absolute_peak_percent", "maximum_pointwise_relative_percent", "maximum_absolute_difference", "zero_reference_count", "maximum_zero_reference_absolute_difference"))
    for name, result in metrics:
        writer.writerow((name,) + tuple(result[key] for key in ("relative_l2_percent", "relative_absolute_peak_percent", "maximum_pointwise_relative_percent", "maximum_absolute_difference", "zero_reference_count", "maximum_zero_reference_absolute_difference")))

passed = all(
    result["relative_l2_percent"] < 0.5
    and result["relative_absolute_peak_percent"] < 0.5
    and result["maximum_pointwise_relative_percent"] < 0.5
    for _, result in metrics
)
if not passed:
    raise SystemExit("B6.0 C3D20T comparison exceeds the 0.5 percent acceptance boundary")
