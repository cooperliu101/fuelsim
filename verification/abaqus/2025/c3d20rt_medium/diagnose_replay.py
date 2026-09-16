#!/usr/bin/env python3
"""Summarize prescribed-state diagnostics; this is not an acceptance checker."""

import argparse
import csv
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from compare import columns, hex20_shapes, read_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path, help="Output prefix from fuelsim_c3d20rt_state_replay")
    parser.add_argument("mesh", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference-directory", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.reference_directory
    nodes = read_reference(root / "reference_nodes.csv.gz").reshape(20, 5969)
    contact = read_reference(root / "reference_contact.csv.gz").reshape(20, 416)
    points = read_reference(root / "reference_points.csv.gz").reshape(20, 1152, 8)
    thermal = read_reference(str(args.prefix) + "_thermal.csv")
    stresses = read_reference(str(args.prefix) + "_stress.csv").reshape(20, 1152, 8)
    residual = read_reference(str(args.prefix) + "_residual.csv").reshape(20, 1617)
    times = nodes[:, 0]["time"]
    if len(thermal) != 20 or not np.allclose(thermal["time"], times, atol=5e-8, rtol=0):
        raise RuntimeError("Incomplete prescribed-state history")
    with Dataset(args.mesh) as database:
        connectivity = np.vstack([database["connect1"][:], database["connect2"][:]]) - 1
        coordinates = np.column_stack([database["coord" + axis][:] for axis in "xyz"])
    gauss = 1 / np.sqrt(3)
    shapes = np.asarray([hex20_shapes(xi, eta, zeta) for zeta in (-gauss, gauss)
                         for eta in (-gauss, gauss) for xi in (-gauss, gauss)])
    rows = []
    components = ("xx", "yy", "zz", "xy", "yz", "xz")
    for frame, time in enumerate(times):
        current = coordinates + columns(nodes[frame], ["u" + axis for axis in "xyz"])
        expected = np.einsum("qn,enc->eqc", shapes, current[connectivity])
        native = np.stack([points[frame][axis] for axis in "xyz"], axis=2)
        distances = np.linalg.norm(expected[:, :, None, :] - native[:, None, :, :], axis=3)
        permutation = np.argmin(distances, axis=2)
        if not np.all(np.sort(permutation, axis=1) == np.arange(8)) or np.max(np.min(distances, axis=2)) > 1e-8:
            raise RuntimeError("Native material-point matching failed")
        reference = np.take_along_axis(points[frame], permutation, axis=1)
        if not np.all(stresses[frame]["element"] == np.arange(1, 1153)[:, None]) or not np.all(
                stresses[frame]["q"] == np.arange(8)[None, :]):
            raise RuntimeError("Replay material-point coverage differs")
        actual = np.stack([stresses[frame][axis] for axis in components], axis=2)
        reference = np.stack([reference["s_" + axis] for axis in components], axis=2)
        delta = np.linalg.norm(actual - reference, axis=2)
        scale = np.linalg.norm(reference, axis=2)
        if np.any(scale == 0) or not np.all(np.isfinite(delta)):
            raise RuntimeError("This diagnostic expects finite nonzero stress tensors")
        native_heat = abs(np.sum(contact[frame]["heat_rate"]))
        interior = (residual[frame]["node"] <= 801) & ~np.isin(residual[frame]["node"], contact[frame]["node"])
        row = {"time": time, "maximum_native_temperature_K": float(np.max(nodes[frame]["temperature"])),
               "native_interface_heat_W": native_heat,
               "replayed_interface_heat_W": thermal[frame]["interface_heat_rate"],
               "interface_heat_relative_difference_percent":
                   100 * (thermal[frame]["interface_heat_rate"] / native_heat - 1),
               "maximum_fuel_interior_thermal_residual_W": float(np.max(np.abs(residual[frame]["total"][interior])))}
        for name, subset in (("fuel", slice(0, 640)), ("cladding", slice(640, 1152))):
            row[name + "_stress_relative_l2_percent"] = 100 * np.linalg.norm(delta[subset]) / np.linalg.norm(scale[subset])
            row[name + "_stress_maximum_pointwise_percent"] = 100 * np.max(delta[subset] / scale[subset])
            row[name + "_stress_maximum_absolute_difference_Pa"] = float(np.max(delta[subset]))
        rows.append(row)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("Diagnostic only: summarized all 20 frames; no equilibrium solve or acceptance decision")


if __name__ == "__main__":
    main()
