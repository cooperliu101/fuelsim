#!/usr/bin/env python3
"""Compare isolated native reactions against independent geometric weights."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from analyze_thermal_weights import NATIVE_POINTS, SIGNS, metrics, read_csv, write_tsv


HERE = Path(__file__).resolve().parent
PREFIX = HERE / "c3d8t_nonaffine_capacity_probe"


def read_geometry():
    coordinates, elements = {}, {}
    section = ""
    for line in PREFIX.with_suffix(".inp").read_text().splitlines():
        if line.startswith("**"):
            continue
        if line.startswith("*"):
            section = line.split(",")[0].lower()
        elif section in ("*node", "*element"):
            values = [float(value) for value in line.split(",")]
            if section == "*node":
                coordinates[int(values[0])] = values[1:]
            else:
                elements[int(values[0])] = [int(value) for value in values[1:]]
    assert len(coordinates) == 128 and len(elements) == 16
    return coordinates, elements


def main():
    coordinates, elements = read_geometry()
    nodal_path = Path(str(PREFIX) + "_nodal.csv")
    points_path = Path(str(PREFIX) + "_integration.csv")
    nodes = {int(row["node"]): row for row in read_csv(nodal_path)}
    points = {(int(row["element"]), int(row["integration_point"])): row for row in read_csv(points_path)}
    assert len(nodes) == len(points) == 128
    rows, volume_rows = [], []
    maximum_inactive_reaction = 0.0
    for element, labels in elements.items():
        kind = "regular" if element <= 8 else "warped"
        reference = np.array([coordinates[node] for node in labels])
        current = reference + np.array([[nodes[node]["u%d_m" % direction] for direction in range(1, 4)]
                                        for node in labels])
        ref = [metrics(reference, sign / np.sqrt(3.0)) for sign in SIGNS]
        cur = [metrics(current, sign / np.sqrt(3.0)) for sign in SIGNS]
        reference_volume = sum(metric[2] for metric in ref)
        current_volume = sum(metric[2] for metric in cur)
        ratio = current_volume / reference_volume
        average_f = sum(current.T @ metric[1] * metric[2] for metric in ref) / reference_volume
        determinant_average_f = np.linalg.det(average_f)
        volume_rows.append(dict(geometry=kind, element=element, reference_volume_m3=reference_volume,
                                current_volume_m3=current_volume, current_reference_ratio=ratio,
                                determinant_average_f=determinant_average_f,
                                native_ivol_sum_m3=sum(points[element, point]["ivol_m3"] for point in range(1, 9))))
        active = (element - 1) % 8
        for local_node, label in enumerate(labels):
            expected_temperature = 301.0 if local_node == active else 300.0
            assert abs(nodes[label]["temperature_k"] - expected_temperature) < 1.0e-10
            if local_node != active:
                maximum_inactive_reaction = max(maximum_inactive_reaction, abs(nodes[label]["reaction_heat_flux_w"]))
                continue
            capacity_factor = 2000.0 * 3000.0
            ref_corner = metrics(reference, SIGNS[local_node])[2]
            weights = {"reference_gauss": ref[local_node][2], "reference_corner": ref_corner,
                       "current_gauss": cur[local_node][2],
                       "current_corner": metrics(current, SIGNS[local_node])[2],
                       "current_row_sum": sum(metric[0][local_node] * metric[2] for metric in cur),
                       "reference_gauss_volume_ratio": ref[local_node][2] * ratio,
                       "reference_corner_volume_ratio": ref_corner * ratio,
                       "reference_gauss_determinant_average_f": ref[local_node][2] * determinant_average_f,
                       "native_ivol": points[element, NATIVE_POINTS[local_node]]["ivol_m3"]}
            for rule, weight in weights.items():
                prediction = weight * capacity_factor
                native = nodes[label]["reaction_heat_flux_w"]
                rows.append(dict(geometry=kind, element=element, local_node=local_node + 1, node=label, rule=rule,
                                 weight_m3=weight, prediction_w=prediction, native_rfl_w=native,
                                 relative_error=abs(prediction - native) / abs(native), absolute_error_w=abs(prediction - native)))
    summary = []
    for kind in ["regular", "warped"]:
        for rule in sorted({row["rule"] for row in rows}):
            selected = [row for row in rows if row["geometry"] == kind and row["rule"] == rule]
            summary.append(dict(geometry=kind, rule=rule, samples=len(selected),
                                maximum_relative_error=max(row["relative_error"] for row in selected),
                                maximum_absolute_error_w=max(row["absolute_error_w"] for row in selected)))
    write_tsv("nonaffine_capacity_candidates.tsv", rows)
    write_tsv("nonaffine_capacity_summary.tsv", summary)
    write_tsv("nonaffine_capacity_volumes.tsv", volume_rows)
    provenance = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in [PREFIX.with_suffix(".inp"), nodal_path, points_path, Path(__file__)]}
    provenance["maximum_inactive_reaction_w"] = maximum_inactive_reaction
    (HERE / "nonaffine_capacity_analysis_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    for row in summary:
        print(row)
    print("maximum_inactive_reaction_w", maximum_inactive_reaction)


if __name__ == "__main__":
    main()
