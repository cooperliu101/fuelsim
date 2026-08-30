#!/usr/bin/env python3
"""Compare isolated Abaqus finite hourglass forces with candidate updates."""

import csv
import sys

import numpy as np

from analyze_b532_candidates import geometry
from generate_b533 import probe_cases


def read(path):
    result = {}
    with open(path, newline="") as source:
        for row in csv.DictReader(source):
            values = result.setdefault(row["case"], np.zeros((8, 3)))
            values[int(row["node"]) - 1] = [float(row["rf_x_n"]), float(row["rf_y_n"]), float(row["rf_z_n"])]
    return result


def finite_total_stiffness(coordinates, displacement):
    volume, gradient, _, gamma = geometry(coordinates)
    deformation = np.eye(3) + displacement.T @ gradient
    metrics = np.sqrt(2.0) * volume * np.sum(gradient * gradient, axis=0) / 6.0
    amplitude = gamma.T @ displacement
    transported = amplitude @ deformation
    weighted = transported * metrics
    modal_force = weighted @ deformation.T
    result = gamma @ modal_force
    for node in range(8):
        gradient_force = weighted @ gradient[node]
        for component in range(3):
            result[node, component] += np.sum(amplitude[:, component] * gradient_force)
    return 0.005 * (2.0e11 / 2.5) * result


def main(path):
    reference = read(path)
    for case_name, case_coordinates, displacement in probe_cases():
        geometry_name = case_name.split("_", 1)[0]
        state_name = case_name[len(geometry_name) + 1 :]
        coordinates = np.asarray(case_coordinates)
        actual = reference["%s_%s_default" % (geometry_name, state_name)] - reference[
            "%s_%s_weak" % (geometry_name, state_name)
        ]
        predicted = finite_total_stiffness(coordinates, np.asarray(displacement))
        error = np.linalg.norm(predicted - actual) / np.linalg.norm(actual)
        scale = np.sum(predicted * actual) / np.sum(predicted * predicted)
        print(
            "%s/%s relative=%.9e projected_scale=%.9e actual_norm=%.9e"
            % (geometry_name, state_name, error, scale, np.linalg.norm(actual))
        )


if __name__ == "__main__":
    main(sys.argv[1])
