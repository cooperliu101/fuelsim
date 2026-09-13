#!/usr/bin/env python3
"""Identify CAX4T Norton parameter temperatures from native material-point data.

This module reads results only. It does not generate or modify Abaqus inputs.
Small-strain inputs explicitly request forward Euler creep; finite-strain inputs
use NLGEOM=YES and the default backward Euler creep integration. The prescribed
temperature and displacement are held fixed during both relaxation steps.

The primary metric compares each native CEEQ increment with dt*A(T)*q**n(T),
using the appropriate stress time level. The competing integration time level
is reported separately; it is not accepted as a substitute for the requested
one. Parameter inversion provides an additional diagnostic. Stress is in Pa;
the exponent probes intentionally keep the tabulated numerical A constant.

Official integration and law definitions:
https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-ratedepcreep.htm
https://docs.software.vt.edu/abaqusv2025/English/SIMACAEKEYRefMap/simakey-r-coupledtemperature-displacement.htm
"""

import argparse
import csv
import json
import math
from pathlib import Path


CASES = (
    "creep_coefficient",
    "creep_coefficient_finite",
    "creep_exponent",
    "creep_exponent_finite",
)
RELAXATION_STEPS = ("RELAX_FAST", "RELAX_SLOW")
NODE_TEMPERATURES = {1: 600.0, 2: 650.0, 3: 700.0, 4: 750.0}
# Native Abaqus CAX4T point numbering: (-,-), (+,-), (-,+), (+,+).
POINT_NODES = {1: 1, 2: 2, 3: 4, 4: 3}
GAUSS_SIGNS = {1: (-1, -1), 2: (1, -1), 3: (-1, 1), 4: (1, 1)}
RELATIVE_TOLERANCE = 0.001


def _read(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("No result rows in {}".format(path))
    return rows


def _number(row, field):
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError("Non-finite {} in {}".format(field, row))
    return value


def _frame(row):
    return row["step"], int(row["frame"])


def _q(row):
    rr, zz, hoop, rz = (_number(row, "s_" + name) for name in ("rr", "zz", "hoop", "rz"))
    return math.sqrt(0.5 * ((rr - zz) ** 2 + (zz - hoop) ** 2 + (hoop - rr) ** 2) + 3.0 * rz * rz)


def _parameters(case, temperature):
    if case.startswith("creep_coefficient"):
        return 1.0e-17 + (temperature - 600.0) * (8.0e-17 / 150.0), 2.0
    if case.startswith("creep_exponent"):
        return 1.0e-19, 2.2 + (temperature - 600.0) * 0.001
    raise ValueError("Unknown creep probe {}".format(case))


def _temperatures(point, nodes):
    xi_sign, eta_sign = GAUSS_SIGNS[point]
    xi, eta = xi_sign / math.sqrt(3.0), eta_sign / math.sqrt(3.0)
    shape = (
        (1.0 - xi) * (1.0 - eta) / 4.0,
        (1.0 + xi) * (1.0 - eta) / 4.0,
        (1.0 + xi) * (1.0 + eta) / 4.0,
        (1.0 - xi) * (1.0 + eta) / 4.0,
    )
    values = [_number(nodes[node], "temperature") for node in range(1, 5)]
    return {
        "corner": values[POINT_NODES[point] - 1],
        "gauss_interpolated": sum(weight * value for weight, value in zip(shape, values)),
        "element_average": sum(values) / 4.0,
    }


def _samples(case, point_rows, node_rows):
    points = {}
    nodes = {}
    for row in point_rows:
        if row["step"] not in RELAXATION_STEPS:
            continue
        frame = _frame(row)
        point = int(row["point"])
        if int(row["element"]) != 1 or point not in POINT_NODES:
            raise ValueError("Unexpected material-point identity in {}".format(case))
        target = points.setdefault(frame, {})
        if point in target:
            raise ValueError("Duplicate material point in {} {}".format(case, frame))
        target[point] = row
    for row in node_rows:
        if row["step"] not in RELAXATION_STEPS:
            continue
        frame = _frame(row)
        node = int(row["node"])
        target = nodes.setdefault(frame, {})
        if node in target or node not in NODE_TEMPERATURES:
            raise ValueError("Duplicate or unexpected node in {} {}".format(case, frame))
        target[node] = row
        if abs(_number(row, "temperature") - NODE_TEMPERATURES[node]) > 1.0e-8:
            raise ValueError("The prescribed temperature was not held fixed in {}".format(case))
        expected_uz = 0.0 if node in (1, 2) else 1.0e-4
        if abs(_number(row, "ur")) > 1.0e-12 or abs(_number(row, "uz") - expected_uz) > 1.0e-12:
            raise ValueError("The prescribed displacement was not held fixed in {}".format(case))
    samples = []
    implicit = case.endswith("_finite")
    for step in RELAXATION_STEPS:
        frames = sorted(frame for name, frame in points if name == step)
        if len(frames) < 3 or frames != list(range(len(frames))):
            raise ValueError("Need frame zero and every relaxation increment in {} {}".format(case, step))
        for frame in frames:
            key = (step, frame)
            if set(points[key]) != set(POINT_NODES) or set(nodes.get(key, {})) != set(NODE_TEMPERATURES):
                raise ValueError("Incomplete four-node/four-point coverage in {} {}".format(case, key))
        for frame in frames[1:]:
            for point in range(1, 5):
                old = points[(step, frame - 1)][point]
                new = points[(step, frame)][point]
                dt = _number(new, "time") - _number(old, "time")
                increment = _number(new, "ceeq") - _number(old, "ceeq")
                if dt <= 0.0 or increment <= 0.0:
                    raise ValueError("Non-positive time/creep increment in {} {} point {}".format(case, frame, point))
                stress_row = new if implicit else old
                stress = _q(stress_row)
                other_stress = _q(old if implicit else new)
                if stress <= 1.0 or other_stress <= 1.0:
                    raise ValueError("Stress too small for Norton parameter inversion in {}".format(case))
                samples.append({
                    "step": step,
                    "frame": frame,
                    "point": point,
                    "dt": dt,
                    "increment": increment,
                    "stress": stress,
                    "other_stress": other_stress,
                    "temperatures": _temperatures(point, nodes[(step, frame)]),
                })
    increments = {round(sample["dt"], 12) for sample in samples}
    if len(increments) < 2:
        raise ValueError("Two different increment sizes were not retained in {}".format(case))
    return samples


def analyze_case(case, point_rows, node_rows):
    """Return per-rule metrics from already extracted native rows."""
    if case not in CASES:
        raise ValueError("Unknown creep probe {}".format(case))
    samples = _samples(case, point_rows, node_rows)
    results = []
    for rule in ("corner", "gauss_interpolated", "element_average", "corner_wrong_time_level"):
        relative_errors, absolute_errors = [], []
        parameter_relative_errors, parameter_absolute_errors = [], []
        point_errors = {point: [] for point in POINT_NODES}
        for sample in samples:
            temperature_rule = "corner" if rule == "corner_wrong_time_level" else rule
            temperature = sample["temperatures"][temperature_rule]
            coefficient, exponent = _parameters(case, temperature)
            stress = sample["other_stress"] if rule == "corner_wrong_time_level" else sample["stress"]
            expected = sample["dt"] * coefficient * stress ** exponent
            observed = sample["increment"]
            absolute_error = abs(expected - observed)
            relative_error = absolute_error / abs(observed)
            absolute_errors.append(absolute_error)
            relative_errors.append(relative_error)
            point_errors[sample["point"]].append(relative_error)
            if case.startswith("creep_coefficient"):
                inferred = observed / (sample["dt"] * stress ** exponent)
                expected_parameter = coefficient
            else:
                inferred = math.log(observed / (sample["dt"] * coefficient)) / math.log(stress)
                expected_parameter = exponent
            parameter_absolute_errors.append(abs(inferred - expected_parameter))
            parameter_relative_errors.append(abs(inferred - expected_parameter) / abs(expected_parameter))
        maximum_relative_error = max(relative_errors)
        expected_to_match = rule == "corner"
        rule_matches = maximum_relative_error < RELATIVE_TOLERANCE
        results.append({
            "case": case,
            "rule": rule,
            "quantity": "equivalent_creep_increment",
            "integration": "backward_euler" if case.endswith("_finite") else "forward_euler",
            "samples": len(samples),
            "max_relative_error": maximum_relative_error,
            "max_absolute_error": max(absolute_errors),
            "parameter": "equivalent_creep_increment",
            "inferred_parameter": "A" if case.startswith("creep_coefficient") else "n",
            "max_parameter_relative_error": max(parameter_relative_errors),
            "max_parameter_absolute_error": max(parameter_absolute_errors),
            "max_relative_error_by_point": {str(point): max(errors) for point, errors in point_errors.items()},
            "relative_tolerance": RELATIVE_TOLERANCE,
            "expected_to_match": expected_to_match,
            "rule_matches": rule_matches,
            "pass": rule_matches,
            "accepted": rule_matches == expected_to_match,
        })
    return results


def analyze(directory):
    """Read <case>_points.csv and <case>_nodes.csv for all four probes."""
    directory = Path(directory)
    results = []
    for case in CASES:
        results.extend(analyze_case(case, _read(directory / (case + "_points.csv")),
                                    _read(directory / (case + "_nodes.csv"))))
    for result in results:
        if not result["accepted"]:
            raise AssertionError("{}: ambiguous or rejected creep temperature/time rule {} (error {})".format(
                result["case"], result["rule"], result["max_relative_error"]))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    results = analyze(args.directory)
    print(json.dumps(results, indent=2, sort_keys=True))
    return 0 if all(result["accepted"] for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
