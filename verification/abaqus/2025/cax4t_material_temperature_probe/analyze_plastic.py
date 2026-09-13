#!/usr/bin/env python3
"""Identify CAX4T plastic-parameter temperatures from native S and PEEQ.

No material update is replayed: active Mises plasticity directly gives
q = sigma0(T) + H(T) * PEEQ.  Separate jobs vary only sigma0 or only H.
The callable analyze(directory) returns one metric row per job and hypothesis.
"""

import argparse
import csv
import math
from pathlib import Path
import sys


CASES = (
    "plastic_yield", "plastic_hardening",
    "plastic_yield_finite", "plastic_hardening_finite",
)
STEPS = ("LOAD_FIXED_T", "LOAD_SWAPPED_T")
CORNER = {1: 1, 2: 2, 3: 4, 4: 3}
GAUSS = 1.0 / math.sqrt(3.0)
NATURAL = {
    1: (-GAUSS, -GAUSS), 2: (GAUSS, -GAUSS),
    3: (-GAUSS, GAUSS), 4: (GAUSS, GAUSS),
}
INITIAL_T = {1: 600.0, 2: 650.0, 3: 700.0, 4: 750.0}
FINAL_T = {1: 750.0, 2: 700.0, 3: 650.0, 4: 600.0}
RULES = ("corner", "gauss_interpolated", "element_average")
TOLERANCE = 1.0e-3
# Abaqus stores frameValue as binary32 even when material fields use binary64.
# This tolerance applies only to output-time metadata, not material comparisons.
FRAME_TIME_TOLERANCE = 2.0**-23


def _read(path):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise AssertionError("empty native output: " + str(path))
    return rows


def _key(row):
    return row["step"], int(row["frame"])


def _number(row, field):
    value = float(row[field])
    if not math.isfinite(value):
        raise AssertionError("nonfinite native field " + field)
    return value


def _mises(row):
    rr, zz, hoop, rz = (_number(row, "s_" + component)
                        for component in ("rr", "zz", "hoop", "rz"))
    return math.sqrt(0.5 * ((rr - zz)**2 + (zz - hoop)**2
                           + (hoop - rr)**2) + 3.0 * rz**2)


def _weights(point):
    xi, eta = NATURAL[point]
    return (0.25 * (1 - xi) * (1 - eta),
            0.25 * (1 + xi) * (1 - eta),
            0.25 * (1 + xi) * (1 + eta),
            0.25 * (1 - xi) * (1 + eta))


def _temperature_candidates(point, nodes):
    values = {node: _number(nodes[node], "temperature") for node in range(1, 5)}
    return {
        "corner": values[CORNER[point]],
        "gauss_interpolated": sum(weight * values[node]
                                  for node, weight in enumerate(_weights(point), 1)),
        "element_average": sum(values.values()) / 4.0,
    }


def _check_nodes(key, nodes):
    if set(nodes) != {1, 2, 3, 4}:
        raise AssertionError("incomplete nodal coverage at " + str(key))
    step, frame = key
    progress = frame / 5.0 if step == STEPS[1] else 0.0
    axial = 0.0004 * (frame / 5.0 + STEPS.index(step))
    for node, row in nodes.items():
        expected_t = INITIAL_T[node] + progress * (FINAL_T[node] - INITIAL_T[node])
        if abs(_number(row, "temperature") - expected_t) > 1.0e-6:
            raise AssertionError("prescribed temperature mismatch at " + str((key, node)))
        if abs(_number(row, "ur")) > 1.0e-12:
            raise AssertionError("unexpected radial motion")
        expected_u = axial if node in (3, 4) else 0.0
        if abs(_number(row, "uz") - expected_u) > 1.0e-10:
            raise AssertionError("prescribed axial displacement mismatch")


def _check_point_location(point, row, nodes):
    """Check the IP numbering independently against native coordinates."""
    r, z = _number(row, "r"), _number(row, "z")
    distances = {}
    for node, values in nodes.items():
        current_r = _number(values, "r") + _number(values, "ur")
        current_z = _number(values, "z") + _number(values, "uz")
        distances[node] = (r - current_r)**2 + (z - current_z)**2
    nearest = min(distances, key=distances.get)
    if nearest != CORNER[point]:
        raise AssertionError("native integration-point coordinate contradicts corner mapping")
    if sum(distance == distances[nearest] for distance in distances.values()) != 1:
        raise AssertionError("integration point has no unique nearest corner")


def analyze_case(directory, case):
    """Return hypothesis metrics; reject incomplete or inactive native samples."""
    if case not in CASES:
        raise ValueError("unknown plastic probe: " + case)
    directory = Path(directory)
    point_rows = _read(directory / (case + "_points.csv"))
    node_rows = _read(directory / (case + "_nodes.csv"))
    expected_keys = {(step, frame) for step in STEPS for frame in range(1, 6)}
    point_frames = {}
    node_frames = {}
    for rows, frames, label in ((point_rows, point_frames, "point"),
                                (node_rows, node_frames, "node")):
        for row in rows:
            if int(row["frame"]) == 0:
                continue
            key = _key(row)
            if key not in expected_keys:
                raise AssertionError("unexpected output frame: " + str(key))
            number = int(row[label])
            frame = frames.setdefault(key, {})
            if number in frame:
                raise AssertionError("duplicate native entry: " + str((key, number)))
            frame[number] = row
        if set(frames) != expected_keys:
            raise AssertionError("missing output frames in " + case)

    hardening = "hardening" in case
    errors = {rule: [] for rule in RULES}
    flow_errors = {rule: [] for rule in RULES}
    previous_peeq = {point: 0.0 for point in range(1, 5)}
    min_increment = math.inf
    for key in sorted(expected_keys, key=lambda value: (STEPS.index(value[0]), value[1])):
        points, nodes = point_frames[key], node_frames[key]
        if set(points) != {1, 2, 3, 4}:
            raise AssertionError("incomplete integration-point coverage at " + str(key))
        _check_nodes(key, nodes)
        for point, row in sorted(points.items()):
            if int(row["element"]) != 1:
                raise AssertionError("unexpected element label")
            _check_point_location(point, row, nodes)
            expected_time = STEPS.index(key[0]) + key[1] / 5.0
            if abs(_number(row, "time") - expected_time) > FRAME_TIME_TOLERANCE:
                raise AssertionError("unexpected physical output time")
            peeq = _number(row, "peeq")
            increment = peeq - previous_peeq[point]
            if increment <= 1.0e-8:
                raise AssertionError("plasticity inactive at " + str((case, key, point)))
            if not 0.0 < peeq < 0.2:
                raise AssertionError("PEEQ outside the specified plastic table")
            min_increment = min(min_increment, increment)
            previous_peeq[point] = peeq
            q = _mises(row)
            inferred = (q - 2.0e6) / peeq if hardening else q
            for rule, temperature in _temperature_candidates(point, nodes).items():
                parameter = (2.0e7 + 4.0e5 * (temperature - 600.0) if hardening
                             else 4.0e6 - 2.0e4 * (temperature - 600.0))
                expected_q = 2.0e6 + parameter * peeq if hardening else parameter
                errors[rule].append((abs(inferred - parameter), parameter))
                flow_errors[rule].append(abs(q - expected_q) / abs(expected_q))

    result = []
    for rule in RULES:
        pairs = errors[rule]
        absolute = max(error for error, _ in pairs)
        relative = max(error / abs(reference) for error, reference in pairs)
        result.append({
            "case": case, "parameter": "hardening_modulus" if hardening else "yield_stress",
            "rule": rule, "samples": len(pairs), "max_relative_error": relative,
            "max_absolute_error": absolute, "pass": relative < TOLERANCE,
            "max_flow_stress_relative_error": max(flow_errors[rule]),
            "minimum_peeq_increment": min_increment,
        })
    if not result[0]["pass"]:
        raise AssertionError("paired-corner hypothesis failed: " + str(result[0]))
    if any(row["pass"] for row in result[1:]):
        raise AssertionError("competing temperature hypothesis was not excluded: " + str(result))
    return result


def analyze(directory):
    """Analyze all four complete small/finite strain native plastic probes."""
    return [row for case in CASES for row in analyze_case(directory, case)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=Path(__file__).parent)
    parser.add_argument("--case", choices=CASES)
    args = parser.parse_args()
    rows = analyze_case(args.directory, args.case) if args.case else analyze(args.directory)
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
