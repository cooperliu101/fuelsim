#!/usr/bin/env python3
"""Compare the B6.0 C3D8RT fuel-plate bending nodal fields and timings."""

import argparse
import csv
import math


def read_rows(path):
    with open(path, newline="") as source:
        rows = list(csv.DictReader(source))
    by_id = {int(row["id"]): row for row in rows}
    if len(by_id) != len(rows):
        raise RuntimeError("duplicate node id in " + path)
    return by_id


def metrics(name, actual, reference):
    if len(actual) != len(reference):
        raise RuntimeError("node count differs for " + name)
    difference_squared = sum((a - b) ** 2 for a, b in zip(actual, reference))
    reference_squared = sum(b * b for b in reference)
    maximum_difference = max((abs(a - b) for a, b in zip(actual, reference)), default=0.0)
    maximum_reference = max((abs(b) for b in reference), default=0.0)
    nonzero = [(abs(a - b) / abs(b), a, b) for a, b in zip(actual, reference) if b != 0.0]
    maximum_pointwise = max(nonzero, default=(0.0, 0.0, 0.0))
    zero = [(abs(a - b)) for a, b in zip(actual, reference) if b == 0.0]
    relative_l2 = math.sqrt(difference_squared / reference_squared) if reference_squared else float("nan")
    relative_peak = maximum_difference / maximum_reference if maximum_reference else float("nan")
    print(
        "%s: relative L2=%.9g%% relative absolute peak=%.9g%% maximum pointwise=%.9g%% "
        "maximum absolute difference=%.9g zero references=%d maximum zero-reference absolute difference=%.9g"
        % (
            name,
            100.0 * relative_l2,
            100.0 * relative_peak,
            100.0 * maximum_pointwise[0],
            maximum_difference,
            len(zero),
            max(zero, default=0.0),
        )
    )


def vector_metrics(name, actual, reference):
    if len(actual) != len(reference):
        raise RuntimeError("node count differs for " + name)
    differences = [math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3))) for a, b in zip(actual, reference)]
    references = [math.sqrt(sum(value * value for value in b)) for b in reference]
    difference_squared = sum(value * value for value in differences)
    reference_squared = sum(value * value for value in references)
    maximum_difference = max(differences, default=0.0)
    maximum_reference = max(references, default=0.0)
    nonzero = [difference / reference for difference, reference in zip(differences, references) if reference != 0.0]
    zero = [difference for difference, reference in zip(differences, references) if reference == 0.0]
    print(
        "%s: relative L2=%.9g%% relative absolute peak=%.9g%% maximum pointwise=%.9g%% "
        "maximum absolute difference=%.9g zero references=%d maximum zero-reference absolute difference=%.9g"
        % (
            name,
            100.0 * math.sqrt(difference_squared / reference_squared),
            100.0 * maximum_difference / maximum_reference,
            100.0 * max(nonzero, default=0.0),
            maximum_difference,
            len(zero),
            max(zero, default=0.0),
        )
    )


def timing(path):
    values = {}
    if not path:
        return values
    with open(path) as source:
        for line in source:
            if "=" in line:
                key, value = line.rstrip().split("=", 1)
                try:
                    values[key] = float(value)
                except ValueError:
                    pass
    return values


parser = argparse.ArgumentParser()
parser.add_argument("fuelsim_csv")
parser.add_argument("abaqus_csv")
parser.add_argument("--fuelsim-timing")
parser.add_argument("--abaqus-timing")
parser.add_argument("--fuelsim-external-seconds", type=float,
                    help="external wall-clock time measured around the Fuelsim process")
args = parser.parse_args()

fuelsim = read_rows(args.fuelsim_csv)
abaqus = read_rows(args.abaqus_csv)
if set(fuelsim) != set(abaqus):
    raise RuntimeError("Fuelsim and Abaqus node labels differ")

for field in ("temperature", "displacement_x", "displacement_y", "displacement_z"):
    actual = [float(fuelsim[node][field]) for node in sorted(fuelsim)]
    reference = [float(abaqus[node][field]) for node in sorted(abaqus)]
    metrics(field, actual, reference)

minimum_x = min(float(fuelsim[node]["x"]) for node in fuelsim)
free_nodes = [node for node in sorted(fuelsim) if abs(float(fuelsim[node]["x"]) - minimum_x) >= 1.0e-12]
vector_metrics(
    "free-node displacement vector",
    [tuple(float(fuelsim[node]["displacement_" + component]) for component in ("x", "y", "z"))
     for node in free_nodes],
    [tuple(float(abaqus[node]["displacement_" + component]) for component in ("x", "y", "z"))
     for node in free_nodes],
)

right_x = max(float(fuelsim[node]["x"]) for node in fuelsim)
right = [node for node in sorted(fuelsim) if abs(float(fuelsim[node]["x"]) - right_x) < 1.0e-12]
fuelsim_tip = max(float(fuelsim[node]["displacement_z"]) for node in right) - min(
    float(fuelsim[node]["displacement_z"]) for node in right
)
abaqus_tip = max(float(abaqus[node]["displacement_z"]) for node in right) - min(
    float(abaqus[node]["displacement_z"]) for node in right
)
print("free-right-edge bending range (x=%.12g m): Fuelsim=%.12g m Abaqus=%.12g m absolute difference=%.12g m" %
      (right_x, fuelsim_tip, abaqus_tip, abs(fuelsim_tip - abaqus_tip)))

fuelsim_time = timing(args.fuelsim_timing)
abaqus_time = timing(args.abaqus_timing)
fuelsim_seconds = args.fuelsim_external_seconds
if fuelsim_seconds is None:
    fuelsim_seconds = fuelsim_time.get("external_wall_seconds", fuelsim_time.get("total_seconds"))
if fuelsim_seconds is not None and "abaqus_external_wall_seconds" in abaqus_time:
    ratio = fuelsim_seconds / abaqus_time["abaqus_external_wall_seconds"]
    print("external wall time: Fuelsim=%.9g s Abaqus=%.9g s Fuelsim/Abaqus=%.9g" %
          (fuelsim_seconds, abaqus_time["abaqus_external_wall_seconds"], ratio))
